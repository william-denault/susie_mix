# Resolve helper paths from this source file, on both the cluster and locally.
.sim_sources <- vapply(sys.frames(), function(f) if (is.null(f$ofile)) "" else f$ofile, "")
.sim_source <- normalizePath(tail(.sim_sources[nzchar(.sim_sources)], 1L),
                             winslash = "/", mustWork = TRUE)
source(file.path(dirname(.sim_source), "simulation_design.R"), local = TRUE)
source(file.path(dirname(.sim_source), "../scan_tissue_attempt/workhorse_utils.R"), local = TRUE)
rm(.sim_source, .sim_sources)

sim_mix <- function(
    pve = 0.6,
    n = 500,
    L_add = 0,
    L_rec = 1,
    L_dom = 0,
    L = 10,
    seed = 1,
    all_additive = FALSE,
    min_maf = 0.05,
    hwe_thresh = 1e-8,
    min_n_rec = 5,
    temp_dir = "/project2/mstephens/wdenault/susie_mix/temp_plink/",
    L_prec = 0,
    L_pdom = 0,
    slide_min_obs = 5,
    return_data = FALSE
) {

  for (package in c("data.table", "matrixStats", "susieR", "susieSlide")) {
    if (!requireNamespace(package, quietly = TRUE)) stop("Required package is missing: ", package)
  }

  set.seed(seed)
  counts <- c(L_add, L_rec, L_dom, L_prec, L_pdom)
  sim_scenario_name(counts) # Validate counts, K <= 5 and at most three effect types.
  K <- sum(counts)
  stopifnot(pve > 0, pve < 1, K > 0,
            all(c(L_add, L_rec, L_dom) >= 0),
            all(c(L_add, L_rec, L_dom) == as.integer(c(L_add, L_rec, L_dom))),
            n >= 3, n == as.integer(n),
            slide_min_obs >= 0, slide_min_obs == as.integer(slide_min_obs))

  # ------------------------------------------------------------
  # Read one randomly selected genotype file.
  # ------------------------------------------------------------

  lf <- list.files(temp_dir, pattern = "\\.raw$", full.names = TRUE)
  if (!length(lf)) stop("No .raw files found in temp_dir.")
  raw_file <- lf[sample.int(length(lf), size = 1)]

  geno_all <- data.table::fread(raw_file, data.table = FALSE)
  rownames(geno_all) <- geno_all$IID
  geno_all <- as.matrix(geno_all[, -(1:6), drop = FALSE])
  storage.mode(geno_all) <- "double"

  # Remove SNPs with missing genotypes or no variation.
  keep <- colSums(is.na(geno_all)) == 0
  geno_all <- geno_all[, keep, drop = FALSE]
  keep <- matrixStats::colSds(geno_all) > 0
  geno_all <- geno_all[, keep, drop = FALSE]

  # Orient to the minor allele and filter MAF/HWE using all donors.
  geno_all <- qc_filter_geno(
    X = geno_all,
    hwe_thresh = hwe_thresh,
    maf_min = min_maf
  )$X

  # Actually use the requested sample size.
  if (n > nrow(geno_all)) stop("Not enough donors for the requested n.")
  rows <- sample.int(nrow(geno_all), size = n)
  geno_all <- geno_all[rows, , drop = FALSE]

  # ------------------------------------------------------------
  # Construct and filter each coding separately.
  # ------------------------------------------------------------

  geno_mix_parts <- recode_snp_matrix(geno_all)

  for (coding in names(geno_mix_parts)) {
    X <- geno_mix_parts[[coding]]
    keep <- colSums(X) >= min_n_rec & matrixStats::colSds(X) > 0
    geno_mix_parts[[coding]] <- X[, keep, drop = FALSE]
  }

  geno_all <- geno_mix_parts$additive

  # Causal SNPs must be available under the requested codings.
  # This also allows a paired all-additive control on exactly the same SNPs.
  eligible <- colnames(geno_mix_parts$additive)
  if (L_rec > 0) {
    eligible <- intersect(eligible, colnames(geno_mix_parts$recessive))
  }
  if (L_dom > 0) {
    eligible <- intersect(eligible, colnames(geno_mix_parts$dominant))
  }
  # Partial effects require all three genotype classes to identify their shape.
  # Keep the original eligibility rule unchanged for the existing scenarios.
  if (L_prec + L_pdom > 0) {
    supported <- Reduce(`&`, lapply(0:2, function(g) colSums(geno_all == g) >= min_n_rec))
    eligible <- intersect(eligible, colnames(geno_all)[supported])
  }
  if (length(eligible) < K) stop("Not enough eligible causal SNPs in this locus.")

  for (coding in names(geno_mix_parts)) {
    if (ncol(geno_mix_parts[[coding]]) > 0) {
      colnames(geno_mix_parts[[coding]]) <- paste0(
        colnames(geno_mix_parts[[coding]]), "__", coding
      )
    }
  }

  geno_mix_all <- cbind(
    geno_mix_parts$additive,
    geno_mix_parts$recessive,
    geno_mix_parts$dominant
  )
  storage.mode(geno_all) <- "double"
  storage.mode(geno_mix_all) <- "double"

  # The coding blocks can have different numbers of columns after filtering.
  mix_to_add <- match(
    sub("__(additive|recessive|dominant)$", "", colnames(geno_mix_all)),
    colnames(geno_all)
  )
  stopifnot(!anyNA(mix_to_add))

  # ------------------------------------------------------------
  # Select distinct causal SNPs and map their column numbers.
  # ------------------------------------------------------------

  causal_snps <- eligible[sample.int(length(eligible), size = K)]

  causal_coding <- c(
    rep("additive", L_add),
    rep("recessive", L_rec),
    rep("dominant", L_dom),
    rep("partial_recessive", L_prec),
    rep("partial_dominant", L_pdom)
  )

  # Change only the generating coding for the paired additive control.
  if (all_additive) causal_coding[] <- "additive"

  true_pos <- match(causal_snps, colnames(geno_all))
  true_pos_mix <- match(
    paste0(causal_snps, "__", causal_coding),
    colnames(geno_mix_all)
  )
  stopifnot(!anyNA(true_pos))
  causal_delta <- unname(sim_effect_delta[causal_coding])

  # ------------------------------------------------------------
  # Simulate phenotype with equal individual contribution variances.
  # ------------------------------------------------------------

  # Preserve the exact old endpoint predictors. Partial codings are generating
  # predictors only: SuSiE-mix still fits its original three coding blocks.
  X_causal <- geno_all[, true_pos, drop = FALSE]
  for (k in seq_len(K)) {
    if (!is.na(true_pos_mix[k])) {
      X_causal[, k] <- geno_mix_all[, true_pos_mix[k]]
    } else {
      X_causal[, k] <- X_causal[, k] + causal_delta[k] * (X_causal[, k] == 1)
    }
  }
  Z <- scale(X_causal)

  beta <- sample(c(-1, 1), size = K, replace = TRUE)
  g <- drop(Z %*% beta)
  if (!is.finite(var(g)) || var(g) <= 0) stop("The genetic signal has zero variance.")

  # Account for LD: var(g) = pve; Gaussian error variance = 1 - pve.
  beta <- beta * sqrt(pve / var(g))
  g <- drop(Z %*% beta)
  noise <- rnorm(n, mean = 0, sd = sqrt(1 - pve))
  y <- g + noise

  # Let the initialization experiment reuse this exact generator. Normal
  # simulation jobs keep the existing fits and compact return value below.
  if (return_data) return(list(
    X = geno_all, y = y, true_pos = true_pos, causal_snps = causal_snps,
    causal_coding = causal_coding, raw_file = raw_file,
    beta_standardized = beta, genetic_variance = var(g),
    phenotype_variance = var(y), seed = seed
  ))

  # ------------------------------------------------------------
  # Fit all three methods to the same phenotype; namespace calls avoid masking.
  # ------------------------------------------------------------

  susie_res <- susieR::susie(
    X = geno_all, y = y, L = L,
    standardize = TRUE, estimate_prior_method = "optim",
    coverage = 0.95, min_abs_corr = 0.5, max_iter = 1000
  )

  susie_res_mix <- susieR::susie(
    X = geno_mix_all, y = y, L = L,
    standardize = TRUE, estimate_prior_method = "optim",
    coverage = 0.95, min_abs_corr = 0.5, max_iter = 1000
  )

  susie_res_slide <- susieSlide::susie(
    X = geno_all, y = y, L = L, min_obs = slide_min_obs,
    standardize = TRUE, estimate_prior_method = "optim",
    coverage = 0.95, min_abs_corr = 0.5, max_iter = 1000
  )

  # ------------------------------------------------------------
  # Calculate one mixed PIP per biological SNP.
  # ------------------------------------------------------------

  # Combine codings within each single effect, then calculate the PIPs.
  # Adding the final coding-specific PIPs would give a different quantity.
  alpha_snp <- matrix(
    0,
    nrow = nrow(susie_res_mix$alpha),
    ncol = ncol(geno_all)
  )

  for (j in seq_along(mix_to_add)) {
    k <- mix_to_add[j]
    alpha_snp[, k] <- alpha_snp[, k] + susie_res_mix$alpha[, j]
  }

  # Keep the fitted prior variances so inactive effects are excluded as usual.
  fit_snp <- susie_res_mix
  fit_snp$alpha <- alpha_snp
  pip_mix_snp <- susieR::susie_get_pip(fit_snp)

  # ------------------------------------------------------------
  # Map mixed credible sets back to biological SNPs in geno_all.
  # ------------------------------------------------------------

  cs_add <- susie_res$sets$cs
  cs_mix <- lapply(susie_res_mix$sets$cs, function(idx) {
    unique(mix_to_add[idx])
  })

  # A false CS contains none of the generating causal SNPs.
  methods <- c("SuSiE", "SuSiE-mix", "SuSiE-slide")
  fits <- list(susie_res, susie_res_mix, susie_res_slide)
  sets <- list(cs_add, cs_mix, susie_res_slide$sets$cs)
  metrics <- do.call(rbind, lapply(seq_along(methods), function(m) {
    hit <- vapply(sets[[m]], function(cs) any(cs %in% true_pos), logical(1))
    data.frame(method = methods[m], converged = fits[[m]]$converged,
               n_cs = length(hit), false_cs = sum(!hit),
               cs_coverage = if (length(hit)) mean(hit) else NA_real_,
               causal_recall = mean(true_pos %in% unlist(sets[[m]])))
  }))
  metrics$cs_fdp <- metrics$false_cs / pmax(1, metrics$n_cs)

  # Keep every PIP for calibration, ROC and power-FDR plots.
  # Omit long PIP names, phenotype vectors and full fitted models from the save.
  return(list(
    settings = list(
      schema_version = sim_schema_version,
      pve = pve, n = n, L = L,
      L_add = L_add, L_rec = L_rec, L_dom = L_dom,
      L_prec = L_prec, L_pdom = L_pdom,
      delta_prec = -0.5, delta_pdom = 0.5, slide_min_obs = slide_min_obs,
      all_additive = all_additive,
      min_maf = min_maf, hwe_thresh = hwe_thresh, min_n_rec = min_n_rec
    ),
    seed = seed,
    metrics = metrics,
    raw_file = raw_file,
    causal_snps = causal_snps,
    causal_coding = causal_coding,
    causal_delta = causal_delta,
    true_pos = true_pos,
    true_pos_mix = true_pos_mix,
    mix_to_add = mix_to_add,
    mix_coding = sub(".*__", "", colnames(geno_mix_all)),
    susie_cs = susie_res$sets,
    susie_pip = unname(susie_res$pip),
    susie_mix_cs = susie_res_mix$sets,
    susie_mix_pip = unname(susie_res_mix$pip),
    susie_mix_pip_snp = unname(pip_mix_snp),
    cs_mix_as_additive_indices = cs_mix,
    susie_slide_cs = susie_res_slide$sets,
    susie_slide_pip = unname(susie_res_slide$pip),
    susie_slide_delta_cs = susie_res_slide$delta_cs,
    susie_slide_delta_causal = susie_res_slide$delta[, true_pos, drop = FALSE],
    susie_slide_delta_forced = unname(susie_res_slide$delta_forced),
    genetic_variance = var(g),
    beta_standardized = beta
  ))
}
