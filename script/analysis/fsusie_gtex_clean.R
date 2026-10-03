# GTEx coverage-profile fine-mapping with fSuSiE and flashier + mvSuSiE.
# On the cluster:
#   source("/path/to/fsusie_gtex_clean.R")
#   out <- run_fsusie_gtex()
# Or run: Rscript fsusie_gtex_clean.R
#
# Edit the settings below. Coverage and genotypes must use the same genome
# assembly (the existing GTEx v8 analysis uses GRCh38). Coverage columns must
# represent consecutive genomic bases; omitted bases are not assumed zero.
# This preserves the supplied analysis's library-size normalization and model
# settings. It does not add ancestry, batch, or other covariate adjustment.

fsusie_gtex_config <- list(
  ensembl_id = "ENSG00000112081",
  tissue = "Brain_Cortex",
  window_bp = 200000,
  n_bins = 1024L,
  min_maf = 0.05,
  seed = 1L,
  counts_dir = "/project2/mstephens/cfbuenabadn/gtex-stm/code/coverage/counts_filtered",
  genotype_prefix = paste0(
    "/project2/mstephens/gtex/",
    "GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv"
  ),
  plink_exec = "/project2/mstephens/gtex/plink2",
  library_size_file = "/project2/mstephens/wdenault/GTEX_analysis_Fsusie/sum_count.RData",
  # If sum_count has several columns, explicitly name the library-total column.
  library_size_column = NULL,
  temp_dir = "/project2/mstephens/fsusie_gtex/temp",
  output_dir = "/project2/mstephens/wdenault/ERC_case_study",
  run_ebmf = TRUE,
  max_ebmf_factors = 20L
)

# Return individuals x variants, counting the ALT allele in this VCF.
fsusie_decode_gt <- function(gt) {
  if (!is.matrix(gt)) stop("GT must be a variants-by-samples matrix.")
  calls <- chartr("|", "/", as.vector(gt))
  missing <- is.na(calls) | calls %in% c(".", "./.", "0/.", "1/.", "./0", "./1")
  valid <- calls %in% c("0/0", "0/1", "1/0", "1/1")
  if (any(!missing & !valid)) {
    stop("Unexpected GT calls: this script requires diploid, biallelic genotypes.")
  }
  lookup <- c("0/0" = 0, "0/1" = 1, "1/0" = 1, "1/1" = 2)
  dosage <- unname(lookup[calls])
  t(matrix(dosage, nrow = nrow(gt), ncol = ncol(gt), dimnames = dimnames(gt)))
}

# Equal-width genomic bins, with each observed base assigned exactly once.
fsusie_bin_counts <- function(counts, bp, requested_bins) {
  if (!is.matrix(counts) || !is.numeric(counts) || ncol(counts) != length(bp) ||
      length(bp) < 2L || any(!is.finite(bp)) || any(diff(bp) != 1)) {
    stop("Coverage must be a numeric matrix on a sorted, consecutive genomic grid.")
  }
  if (any(!is.finite(counts)) || any(counts < 0)) {
    stop("Coverage contains missing, nonfinite, or negative counts.")
  }
  if (length(requested_bins) != 1L || !is.finite(requested_bins) || requested_bins < 2) {
    stop("n_bins must be at least 2.")
  }
  n_bins <- as.integer(2^floor(log2(min(requested_bins, length(bp)))))
  bin_id <- floor((seq_along(bp) - 1) * n_bins / length(bp)) + 1L
  Y <- t(rowsum(t(counts), group = bin_id, reorder = FALSE))
  colnames(Y) <- paste0("bin_", seq_len(n_bins))
  # Boundaries are half-open genomic intervals; positions are bin centers.
  edges <- seq(min(bp), max(bp) + 1, length.out = n_bins + 1L)
  stopifnot(isTRUE(all.equal(unname(rowSums(Y)), unname(rowSums(counts)))))
  list(
    Y = Y,
    start = head(edges, -1L),
    end = tail(edges, -1L),
    pos = (head(edges, -1L) + tail(edges, -1L)) / 2,
    n_positions = tabulate(bin_id, nbins = n_bins)
  )
}

# Match full RNA sample IDs, never donor IDs, to the library-size table.
fsusie_match_library_sizes <- function(sum_count, sample_ids, column = NULL) {
  if (is.numeric(sum_count) && is.null(dim(sum_count))) {
    ids <- names(sum_count)
    values <- as.numeric(sum_count)
  } else if (is.matrix(sum_count) || is.data.frame(sum_count)) {
    ids <- rownames(sum_count)
    if (is.null(column)) {
      if (ncol(sum_count) != 1L) {
        stop("sum_count has several columns; set library_size_column explicitly.")
      }
      column <- 1L
    }
    values <- sum_count[, column, drop = TRUE]
    if (!is.numeric(values)) stop("The library-size column must be numeric.")
  } else {
    stop("sum_count must be a named numeric vector or a table with sample row names.")
  }
  if (is.null(ids) || anyNA(ids) || any(!nzchar(ids))) {
    stop("sum_count must have RNA sample IDs as names/row names.")
  }
  # The supplied sum_count table uses dots where coverage IDs use hyphens.
  keys <- chartr(".", "-", ids)
  wanted <- chartr(".", "-", sample_ids)
  if (anyDuplicated(keys) || anyDuplicated(wanted)) {
    stop("Duplicate RNA sample IDs in the library-size table or selected coverage.")
  }
  idx <- match(wanted, keys)
  if (anyNA(idx)) {
    stop("No library size for: ", paste(sample_ids[is.na(idx)], collapse = ", "))
  }
  sizes <- as.numeric(values[idx])
  if (any(!is.finite(sizes)) || any(sizes <= 0)) {
    stop("Matched library sizes must be finite and strictly positive.")
  }
  setNames(sizes, sample_ids)
}

run_fsusie_gtex <- function(config = fsusie_gtex_config) {
  packages <- c("data.table", "vcfR", "fsusieR")
  if (isTRUE(config$run_ebmf)) packages <- c(packages, "flashier", "mvsusieR", "ebnm")
  available <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  if (any(!available)) stop("Missing R packages: ", paste(packages[!available], collapse = ", "))
  if (length(config$window_bp) != 1L || !is.finite(config$window_bp) ||
      config$window_bp < 0 || config$window_bp != floor(config$window_bp)) {
    stop("window_bp must be a nonnegative integer.")
  }
  if (length(config$min_maf) != 1L || !is.finite(config$min_maf) ||
      config$min_maf < 0 || config$min_maf > 0.5) stop("min_maf must be in [0, 0.5].")
  if (config$max_ebmf_factors < 1) stop("max_ebmf_factors must be positive.")
  set.seed(config$seed)

  required <- c(config$plink_exec, config$library_size_file,
                paste0(config$genotype_prefix, c(".bed", ".bim", ".fam")))
  missing <- required[file.access(required, 4) != 0L]
  if (length(missing)) stop("Missing or unreadable files:\n", paste(missing, collapse = "\n"))
  if (file.access(config$plink_exec, 1) != 0L) stop("PLINK is not executable: ", config$plink_exec)
  for (directory in c(config$temp_dir, config$output_dir)) {
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(directory) || file.access(directory, 2) != 0L) {
      stop("Cannot write to directory: ", directory)
    }
  }

  # 1. Read coverage and select the requested tissue.
  files <- list.files(config$counts_dir, full.names = TRUE)
  hit <- grepl(paste0(config$ensembl_id, "([^0-9]|$)"), basename(files))
  files <- files[hit & !dir.exists(files)]
  if (length(files) != 1L) {
    stop("Expected exactly one coverage file for ", config$ensembl_id,
         "; found ", length(files), ". Check counts_dir.")
  }
  message("Reading coverage: ", files)
  tt <- data.table::fread(files, sep = ",", header = TRUE, check.names = FALSE)
  if (!"Sample_ID" %in% names(tt)) stop("Coverage must contain a Sample_ID column.")
  coverage_cols <- setdiff(names(tt), "Sample_ID")
  if (!length(coverage_cols) || any(!grepl("^chr[^:]+:[0-9]+$", coverage_cols))) {
    stop("Every coverage column must have a name such as chr6:36594312.")
  }
  chromosome <- unique(sub(":.*$", "", coverage_cols))
  if (length(chromosome) != 1L || !chromosome %in% paste0("chr", 1:22)) {
    stop("This script expects one autosomal locus (diploid genotype coding).")
  }
  bp <- as.numeric(sub("^.*:", "", coverage_cols))
  ord <- order(bp)
  bp <- bp[ord]
  coverage_cols <- coverage_cols[ord]
  if (length(bp) < 2L || any(diff(bp) != 1)) {
    stop("Coverage positions must be consecutive and unique; omitted bases cannot be assumed zero.")
  }
  sample_names <- as.character(tt$Sample_ID)
  if (anyNA(sample_names) || any(!grepl("^GTEX-[^-]+-.+[.].+$", sample_names))) {
    stop("Unexpected Sample_ID format; expected a GTEx RNA sample ID followed by .tissue.")
  }
  sample_info <- data.frame(
    ind = sub("^(GTEX-[^-]+)-.*$", "\\1", sample_names),
    tissue = sub("^.*[.]", "", sample_names),
    sample_id = sub("[.].*$", "", sample_names),
    full_name = sample_names,
    stringsAsFactors = FALSE
  )
  rows <- which(sample_info$tissue == config$tissue)
  if (!length(rows)) stop("No coverage samples for tissue ", config$tissue)
  sample_info <- sample_info[rows, , drop = FALSE]
  if (anyDuplicated(sample_info$ind)) stop("More than one RNA sample per donor in the selected tissue.")
  counts <- as.matrix(tt[rows, coverage_cols, with = FALSE])
  if (!is.numeric(counts)) stop("Coverage columns must be numeric.")
  storage.mode(counts) <- "double"
  rm(tt)

  # 2. Extract the locus from the binary files used by the existing GTEx analysis.
  lower <- max(1, min(bp) - config$window_bp)
  upper <- max(bp) + config$window_bp
  prefix <- tempfile(paste0(config$ensembl_id, "_"), tmpdir = config$temp_dir)
  vcf_file <- paste0(prefix, ".vcf")
  message("Extracting ", chromosome, ":", lower, "-", upper)
  status <- system2(config$plink_exec, args = c(
    "--bfile", shQuote(config$genotype_prefix),
    "--chr", sub("^chr", "", chromosome),
    "--from-bp", format(lower, scientific = FALSE, trim = TRUE),
    "--to-bp", format(upper, scientific = FALSE, trim = TRUE),
    "--snps-only", "--max-alleles", "2", "--rm-dup", "exclude-all",
    "--threads", "2", "--memory", "8000",
    "--export", "vcf-4.2", "id-paste=iid", "--out", shQuote(prefix)
  ))
  if (status != 0L || !file.exists(vcf_file)) stop("PLINK extraction failed. See ", prefix, ".log")
  vcf <- vcfR::read.vcfR(vcf_file, verbose = FALSE)
  info_SNP <- as.data.frame(vcfR::getFIX(vcf), stringsAsFactors = FALSE)
  gt <- vcfR::extract.gt(vcf, element = "GT", as.numeric = FALSE, IDtoRowNames = FALSE)
  if (!nrow(info_SNP) || nrow(gt) != nrow(info_SNP)) stop("Empty or inconsistent VCF.")
  X <- fsusie_decode_gt(gt)
  rm(vcf, gt)
  info_SNP$POS <- as.numeric(info_SNP$POS)
  info_SNP$variant_key <- with(info_SNP, paste(CHROM, POS, REF, ALT, sep = ":"))
  if (anyDuplicated(info_SNP$variant_key)) stop("Duplicate physical variants in the extracted VCF.")
  colnames(X) <- info_SNP$variant_key
  rownames(info_SNP) <- info_SNP$variant_key
  # REF/ALT in a VCF reconstructed from BED may be provisional. X counts the
  # exported ALT; verify against the genome reference before interpreting signs.
  info_SNP$counted_allele <- info_SNP$ALT

  # 3. Align donors explicitly; keep only the selected tissue's matched samples.
  genotype_ids <- rownames(X)
  if (is.null(genotype_ids) || anyNA(genotype_ids) || anyDuplicated(genotype_ids)) {
    stop("Genotype sample IDs must be present and unique.")
  }
  matched <- sample_info$ind %in% genotype_ids
  message(sum(matched), "/", nrow(sample_info), " tissue samples have genotypes.")
  sample_info <- sample_info[matched, , drop = FALSE]
  counts <- counts[matched, , drop = FALSE]
  if (nrow(sample_info) < 3L) stop("Fewer than three matched donors.")
  X <- X[match(sample_info$ind, genotype_ids), , drop = FALSE]
  rownames(counts) <- sample_info$ind
  rownames(sample_info) <- sample_info$ind
  stopifnot(identical(rownames(X), rownames(counts)))

  # Complete-case SNP filtering is applied within the actual analysis samples.
  complete <- colSums(is.na(X)) == 0L
  af <- colMeans(X, na.rm = TRUE) / 2
  maf <- pmin(af, 1 - af)
  keep <- complete & is.finite(maf) & maf >= config$min_maf & maf > 0
  message(sum(keep), "/", ncol(X), " SNPs pass missingness and tissue-specific MAF filters.")
  if (!any(keep)) stop("No SNPs remain after filtering.")
  X <- X[, keep, drop = FALSE]
  info_SNP$AF <- af
  info_SNP$MAF <- maf
  info_SNP <- info_SNP[keep, , drop = FALSE]
  stopifnot(identical(colnames(X), info_SNP$variant_key))

  # 4. Bin the coverage, then match and apply the library-size normalization.
  bins <- fsusie_bin_counts(counts, bp, config$n_bins)
  message("Using ", nrow(X), " donors, ", ncol(X), " SNPs, and ", ncol(bins$Y), " bins.")
  sizes_env <- new.env(parent = emptyenv())
  load(config$library_size_file, envir = sizes_env)
  if (!exists("sum_count", envir = sizes_env, inherits = FALSE)) {
    stop("library_size_file does not contain an object named sum_count.")
  }
  library_size <- fsusie_match_library_sizes(
    sizes_env$sum_count, sample_info$sample_id, config$library_size_column
  )
  size_factor <- setNames(as.numeric(library_size / mean(library_size)), sample_info$ind)
  Y_cor <- sweep(bins$Y, 1L, size_factor, "/")
  Y_log <- log1p(Y_cor)
  if (any(!is.finite(Y_log)) || all(Y_log == 0)) stop("Normalized profiles are invalid or all zero.")
  stopifnot(identical(rownames(X), rownames(Y_log)), ncol(Y_log) == length(bins$pos))

  # 5. Keep the fSuSiE settings from the supplied script.
  message("Fitting fSuSiE.")
  res <- fsusieR::susiF(
    Y = Y_log, X = X, L_start = 5, L = 20, pos = bins$pos,
    cal_obj = TRUE, nullweight = 1, verbose = FALSE, maxit = 40,
    max_SNP_EM = 100, cor_small = TRUE, post_processing = "HMM"
  )
  output_file <- file.path(config$output_dir, paste0(config$ensembl_id, "_", config$tissue, ".RData"))
  out <- list(
    res = res, Y = Y_cor, Y_log = Y_log, Y_counts = bins$Y, X = X,
    library_size = library_size, size_factor_global = size_factor,
    info_ind_gen = sample_info, info_SNP = info_SNP,
    locus = range(bp), genotype_window = c(lower, upper),
    pos = bins$pos, start_bin = bins$start, end_bin = bins$end,
    n_positions_per_bin = bins$n_positions, chr = sub("^chr", "", chromosome),
    facto = NULL, mv_res = NULL, Y_EBMF = NULL,
    comparison_status = if (isTRUE(config$run_ebmf)) "pending" else "not_requested",
    config = config, coverage_file = files, extracted_vcf = vcf_file,
    plink_log = paste0(prefix, ".log"),
    normalization = "Coverage bin sums divided by matched library size / mean library size; log1p for fitting.",
    session_info = utils::sessionInfo()
  )
  # Save the main fit first, so a comparison-model failure does not lose it.
  save(out, file = output_file)
  message("Saved fSuSiE result: ", output_file)

  # 6. Optional comparison: flashier loadings followed by mvSuSiE.
  if (isTRUE(config$run_ebmf)) {
    message("Fitting flashier and mvSuSiE comparison.")
    comparison_error <- NULL
    out <- tryCatch({
      out$facto <- flashier::flash(Y_log, ebnm_fn = ebnm::ebnm_point_exponential)
      loadings <- out$facto$L_pm
      if (is.null(loadings) || !is.matrix(loadings) || ncol(loadings) == 0L) {
        out$comparison_status <- "no_factors"
      } else {
        # Constant factors cannot be phenotypes for the association comparison.
        variances <- apply(loadings, 2L, var)
        usable <- which(is.finite(variances) & variances > 0)
        usable <- head(usable, config$max_ebmf_factors)
        if (!length(usable)) {
          out$comparison_status <- "no_variable_factors"
        } else {
          out$Y_EBMF <- loadings[, usable, drop = FALSE]
          rownames(out$Y_EBMF) <- rownames(X)
          out$ebmf_factor_indices <- usable
          prior <- mvsusieR::create_mixture_prior(R = ncol(out$Y_EBMF))
          out$mv_res <- mvsusieR::mvsusie(Y = out$Y_EBMF, X = X, L = 10, prior_variance = prior)
          out$comparison_status <- "complete"
        }
      }
      out
    }, error = function(e) {
      comparison_error <<- conditionMessage(e)
      out$comparison_status <- "failed"
      out$comparison_error <- comparison_error
      out
    })
    out$session_info <- utils::sessionInfo()
    save(out, file = output_file)
    if (!is.null(comparison_error)) {
      warning("Comparison failed, but fSuSiE was saved: ", comparison_error, call. = FALSE)
    } else if (out$comparison_status != "complete") {
      warning("mvSuSiE skipped: ", out$comparison_status, call. = FALSE)
    }
  }
  message("Results saved to ", output_file)
  invisible(out)
}

if (sys.nframe() == 0L) out <- run_fsusie_gtex()
