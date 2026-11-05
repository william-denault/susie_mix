# ------------------------------------------------------------
# Extract the lead SNP and TSS distance for every credible set
# ------------------------------------------------------------

extract_snp_position <- function(snp_names) {

  snp_names <- as.character(snp_names)

  # Expected GTEx/PLINK format:
  # chr14_92014703_A_C_b38_A
  snp_position <- suppressWarnings(
    as.numeric(
      sub(
        "^chr[^_]+_([0-9]+)_.*$",
        "\\1",
        snp_names
      )
    )
  )

  if (anyNA(snp_position)) {
    stop(
      "Could not extract genomic position from SNP name(s): ",
      paste(
        snp_names[is.na(snp_position)],
        collapse = ", "
      )
    )
  }

  snp_position
}


get_cs_lead_tss_distance <- function(
    fit,
    predictor_map,
    tss) {

  cs_list <- fit$sets$cs

  # Required output when SuSiE reports no credible sets.
  if (
    is.null(cs_list) ||
    length(cs_list) == 0L
  ) {
    return(
      setNames(
        numeric(0),
        character(0)
      )
    )
  }

  cs_index <- fit$sets$cs_index

  # Identify the lead predictor from the single-effect
  # component responsible for each credible set.
  lead_predictor_index <- vapply(
    seq_along(cs_list),
    function(i) {

      members <- as.integer(
        cs_list[[i]]
      )

      members <- members[
        members >= 1L &
          members <= ncol(fit$alpha)
      ]

      if (length(members) == 0L) {
        return(NA_integer_)
      }

      if (
        !is.null(cs_index) &&
        length(cs_index) >= i &&
        !is.na(cs_index[i]) &&
        cs_index[i] >= 1L &&
        cs_index[i] <= nrow(fit$alpha)
      ) {

        component_probability <- fit$alpha[
          cs_index[i],
          members
        ]

        members[
          which.max(component_probability)
        ]

      } else {

        # Fallback for older SuSiE objects without cs_index.
        members[
          which.max(fit$pip[members])
        ]
      }
    },
    integer(1L)
  )

  if (anyNA(lead_predictor_index)) {
    stop(
      "A lead predictor could not be identified for every ",
      "credible set."
    )
  }

  map_row <- match(
    lead_predictor_index,
    predictor_map$predictor_index
  )

  if (anyNA(map_row)) {
    stop(
      "A lead predictor was not found in the predictor map."
    )
  }

  # For the mixed model, this maps the coding-specific
  # predictor back to its biological SNP.
  lead_snp <- predictor_map$snp[
    map_row
  ]

  lead_snp_position <- extract_snp_position(
    lead_snp
  )

  distance_to_tss <-  (
    lead_snp_position - as.numeric(tss)
  )

  # Names contain the lead biological SNP.
  names(distance_to_tss) <- lead_snp

  distance_to_tss
}


recode_matrix_by_freq <- function(X) {

  storage.mode(X) <- "integer"

  counts <- rbind(
    colSums(X == 0L),
    colSums(X == 1L),
    colSums(X == 2L)
  )

  new_code <- apply(
    counts,
    2,
    function(cnt) {
      rank(
        -cnt,
        ties.method = "first"
      ) - 1L
    }
  )

  out <- X

  for (v in 0:2) {

    idx <- which(
      X == v,
      arr.ind = TRUE
    )

    if (nrow(idx) > 0L) {
      out[idx] <- new_code[
        v + 1L,
        idx[, "col"]
      ]
    }
  }

  dimnames(out) <- dimnames(X)

  out
}


qc_filter_geno <- function(
    X,
    hwe_thresh = 1e-8,
    maf_min = 0.05) {

  X <- as.matrix(X)
  storage.mode(X) <- "integer"

  # Orient all SNPs to the minor allele.
  af <- colMeans(X) / 2
  flip <- is.finite(af) & af > 0.5

  if (any(flip)) {
    X[, flip] <- 2L - X[, flip, drop = FALSE]
  }

  # Filter by minor-allele frequency.
  maf <- colMeans(X) / 2
  keep_maf <- is.finite(maf) & maf > maf_min

  X <- X[
    ,
    keep_maf,
    drop = FALSE
  ]

  if (ncol(X) == 0L) {
    stop(
      "No SNPs remained after MAF filtering for ",
      target_gene
    )
  }

  n <- nrow(X)
  maf <- colMeans(X) / 2

  # HWE chi-square test.
  count0 <- colSums(X == 0L)
  count1 <- colSums(X == 1L)
  count2 <- colSums(X == 2L)

  exp0 <- n * (1 - maf)^2
  exp1 <- n * 2 * maf * (1 - maf)
  exp2 <- n * maf^2

  valid_hwe <- (
    exp0 > 0 &
      exp1 > 0 &
      exp2 > 0
  )

  chisq <- rep(
    NA_real_,
    ncol(X)
  )

  chisq[valid_hwe] <- (
    (count0[valid_hwe] - exp0[valid_hwe])^2 /
      exp0[valid_hwe] +
      (count1[valid_hwe] - exp1[valid_hwe])^2 /
      exp1[valid_hwe] +
      (count2[valid_hwe] - exp2[valid_hwe])^2 /
      exp2[valid_hwe]
  )

  hwe_p <- pchisq(
    chisq,
    df = 1,
    lower.tail = FALSE
  )

  # Remove SNPs with HWE P < hwe_thresh.
  fail_hwe <- (
    !is.na(hwe_p) &
      hwe_p < hwe_thresh
  )

  removed_hwe_snps <- colnames(X)[fail_hwe]

  if (any(fail_hwe)) {
    cat(
      sum(fail_hwe),
      "SNP(s) removed because HWE P <",
      format(hwe_thresh, scientific = TRUE),
      ":\n"
    )

    print(removed_hwe_snps)

    X <- X[
      ,
      !fail_hwe,
      drop = FALSE
    ]

    maf <- maf[!fail_hwe]
    hwe_p <- hwe_p[!fail_hwe]
  }

  if (ncol(X) == 0L) {
    stop(
      "No SNPs remained after HWE filtering for ",
      target_gene
    )
  }

  list(
    X = X,
    maf = maf,
    hwe_p = hwe_p,
    flipped = which(flip),
    removed_hwe_snps = removed_hwe_snps,
    n_removed_hwe = length(removed_hwe_snps)
  )
}


recode_snp_matrix <- function(
    X,
    warn = TRUE) {

  X <- as.matrix(X)

  if (
    warn &&
    !all(X %in% c(0L, 1L, 2L))
  ) {
    warning(
      paste0(
        "recode_snp_matrix: entries are not all in {0,1,2}; ",
        "expecting additive 0/1/2 dosages."
      ),
      call. = FALSE
    )
  }

  storage.mode(X) <- "integer"

  additive <- X
  dominant <- (X >= 1L) * 1L
  recessive <- (X == 2L) * 1L

  dimnames(additive) <- dimnames(X)
  dimnames(dominant) <- dimnames(X)
  dimnames(recessive) <- dimnames(X)

  list(
    additive = additive,
    recessive = recessive,
    dominant = dominant
  )
}

