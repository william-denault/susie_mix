# Genotype simulation and evaluation helpers. No analysis runs when sourced.

simulation_scenarios <- function(pve = c(0.05, 0.10, 0.20)) {
  architectures <- list(
    additive = "additive", recessive = "recessive", dominant = "dominant",
    additive_recessive = c("additive", "recessive"),
    additive_dominant = c("additive", "dominant"),
    recessive_dominant = c("recessive", "dominant"),
    additive_recessive_dominant = c("additive", "recessive", "dominant")
  )
  out <- do.call(rbind, lapply(names(architectures), function(a) {
    expand.grid(architecture = a, k = seq.int(length(architectures[[a]]), 5L),
                pve = pve, stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out$scenario_id <- seq_len(nrow(out))
  out[, c("scenario_id", "architecture", "k", "pve")]
}

simulation_seed <- function(seed, scenario_id, replication) {
  # Independent deterministic streams for each trial, regardless of shard/order.
  as.integer((as.double(seed) + 1000003 * scenario_id +
                9176 * replication) %% 2147483646 + 1)
}

simulation_helpers <- function(root) {
  env <- new.env(parent = globalenv())
  # get_gene_annotations uses fread without a namespace.
  env$fread <- function(...) data.table::fread(...)
  sys.source(file.path(root, "script/scan_tissue_attempt/workhorse_utils.R"), env)
  sys.source(file.path(root, "script/scan_tissue_attempt/get_gene_annotations.R"), env)
  env
}

simulation_validate <- function(config) {
  positive_integer <- function(x) is.numeric(x) && length(x) == 1L &&
    is.finite(x) && x >= 1 && x == as.integer(x)
  for (nm in c("n", "replications", "seed", "L", "max_iter", "n_purity",
               "max_locus_attempts", "max_causal_attempts")) {
    if (!positive_integer(config[[nm]])) stop(nm, " must be a positive integer.")
  }
  if (config$n != 200L) stop("This benchmark requires N = 200.")
  if (config$L != 10L) stop("Both methods must use L = 10.")
  if (!config$pilot && config$replications < 200L)
    stop("Production requires at least 200 replications per scenario; use --pilot for tests.")
  if (!config$pilot && !is.null(config$raw_file))
    stop("A fixed --raw-file is allowed only for a pilot, never production.")
  if (!identical(sort(as.numeric(config$pve)), c(0.05, 0.1, 0.2)))
    stop("Use all three PVE levels: 0.05, 0.10, 0.20.")
  if (!identical(config$effect_distribution, "equal_standardized_random_sign"))
    stop("Unknown effect distribution.")
  for (nm in c("thresholds", "pip_threshold", "coverage", "min_abs_corr")) {
    x <- config[[nm]]
    if (!is.numeric(x) || !length(x) || any(!is.finite(x) | x < 0 | x > 1))
      stop(nm, " must be between 0 and 1.")
  }
  b <- config$calibration_breaks
  if (length(b) < 2L || any(!is.finite(b)) || any(diff(b) <= 0) ||
      b[1] != 0 || tail(b, 1) != 1) stop("Invalid calibration breaks.")
  invisible(config)
}

simulation_read_raw <- function(path) {
  raw <- data.table::fread(path, header = TRUE, check.names = FALSE,
                          data.table = FALSE, showProgress = FALSE)
  if (ncol(raw) <= 6L || !all(c("FID", "IID") %in% names(raw)[1:6]))
    stop("Invalid PLINK additive .raw file: ", path)
  if (anyDuplicated(raw$IID)) stop("Genotype donor IDs must be unique.")
  X <- as.matrix(raw[, -(1:6), drop = FALSE])
  storage.mode(X) <- "double"
  rownames(X) <- as.character(raw$IID)
  if (anyDuplicated(colnames(X))) stop("Duplicate genotype predictor names.")
  if (any(!is.na(X) & !(X %in% 0:2))) stop("Expected hard-call genotypes 0/1/2.")
  X
}

simulation_global_qc <- function(X, config, helpers, gene) {
  counts <- c(n_donors = nrow(X), n_extracted = ncol(X))
  X <- X[, colSums(is.na(X)) == 0L, drop = FALSE]
  counts <- c(counts, n_complete = ncol(X))
  X <- X[, matrixStats::colSds(X) > 0, drop = FALSE]
  counts <- c(counts, n_variable = ncol(X))
  if (!ncol(X)) return(NULL)
  # The existing QC helper uses target_gene in its empty-result messages.
  helpers$target_gene <- gene
  qc <- tryCatch(helpers$qc_filter_geno(X, hwe_thresh = config$hwe_thresh,
                                      maf_min = config$min_maf),
                 error = function(e) {
                   if (grepl("No SNPs remained after", conditionMessage(e))) return(NULL)
                   stop(e)
                 })
  if (is.null(qc) || !ncol(qc$X)) return(NULL)
  list(X = qc$X, counts = c(counts, n_global_qc = ncol(qc$X)))
}

simulation_design <- function(X, config, helpers, ids = NULL) {
  if (nrow(X) < config$n) return(NULL)
  if (is.null(ids)) ids <- sample(rownames(X), config$n, replace = FALSE)
  if (length(ids) != config$n || anyDuplicated(ids) || any(!ids %in% rownames(X)))
    stop("Invalid analysis donor sample.")
  # Orientation/MAF/HWE have already been established using ALL donors.
  X <- X[ids, , drop = FALSE]
  parts <- helpers$recode_snp_matrix(X)
  snps <- colnames(X)
  coding <- rep(c("additive", "recessive", "dominant"), each = length(snps))
  mix <- do.call(cbind, parts[c("additive", "recessive", "dominant")])
  mix_snps <- rep(snps, 3L)
  colnames(mix) <- paste(mix_snps, coding, sep = "__")
  keep_add <- colSums(X) >= config$min_n_rec & matrixStats::colSds(X) > 0
  keep_mix <- colSums(mix) >= config$min_n_rec & matrixStats::colSds(mix) > 0
  add <- X[, keep_add, drop = FALSE]
  mix <- mix[, keep_mix, drop = FALSE]
  storage.mode(add) <- "double"
  storage.mode(mix) <- "double"
  if (!ncol(add) || !ncol(mix)) return(NULL)
  map <- data.frame(predictor = colnames(mix), snp = mix_snps[keep_mix],
                    coding = coding[keep_mix], stringsAsFactors = FALSE)
  # Under the workhorse filter, every retained nonadditive SNP also retains additive.
  stopifnot(setequal(colnames(add), unique(map$snp)))
  list(add = add, mix = mix, map = map, ids = ids,
       counts = c(n_add = ncol(add), n_mix = ncol(mix),
                  n_recessive = sum(map$coding == "recessive"),
                  n_dominant = sum(map$coding == "dominant")))
}

simulation_causals <- function(design, scenario, config) {
  types <- strsplit(scenario$architecture, "_", fixed = TRUE)[[1]]
  # Uniformly sample labelled coding assignments conditional on every requested
  # class being represented. Each causal SNP has exactly one generating coding.
  repeat {
    codes <- sample(types, scenario$k, replace = TRUE)
    if (all(types %in% codes)) break
  }
  available <- lapply(types, function(type) which(design$map$coding == type))
  names(available) <- types
  if (any(vapply(types, function(type) length(available[[type]]) < sum(codes == type),
                 logical(1)))) return(NULL)
  # Rejection sampling avoids order-dependent selection of overlapping SNP pools.
  for (attempt in seq_len(config$max_causal_attempts)) {
    idx <- integer(scenario$k)
    for (type in types) {
      pool <- available[[type]]
      idx[codes == type] <- pool[sample.int(length(pool), sum(codes == type))]
    }
    if (!anyDuplicated(design$map$snp[idx])) return(idx)
  }
  NULL
}

simulation_phenotype <- function(design, causal_index, pve) {
  W <- design$mix[, causal_index, drop = FALSE]
  signs <- sample(c(-1, 1), ncol(W), replace = TRUE)
  beta <- signs / matrixStats::colSds(W)
  g <- drop(scale(W, center = TRUE, scale = FALSE) %*% beta)
  if (!is.finite(var(g)) || var(g) <= 1e-12) return(NULL)
  multiplier <- sqrt(pve / var(g))
  beta <- beta * multiplier
  g <- g * multiplier
  noise <- rnorm(length(g), sd = sqrt(1 - pve))
  y <- g + noise
  truth <- design$map[causal_index, , drop = FALSE]
  truth$beta <- beta
  truth$sample_effect_allele_frequency <- colMeans(design$add[, truth$snp, drop = FALSE]) / 2
  truth$sample_maf <- pmin(truth$sample_effect_allele_frequency,
                          1 - truth$sample_effect_allele_frequency)
  list(y = y, genetic_value = g, truth = truth,
       variance = c(target_pve = pve, genetic_variance = var(g),
         noise_variance_target = 1 - pve, noise_variance_observed = var(noise),
         phenotype_variance = var(y), genetic_noise_covariance = cov(g, noise),
         realized_variance_ratio = var(g) / var(y),
         realized_r_squared = cor(g, y)^2))
}

simulation_catalogue <- function(config, helpers) {
  needed <- c(config$gtf_file, config$gene_list,
              paste0(config$genotype_prefix, c(".bed", ".bim", ".fam")))
  if (any(!file.exists(needed)))
    stop("Missing GTEx inputs: ", paste(needed[!file.exists(needed)], collapse = ", "))
  if (!file.exists(config$plink_exec) && !nzchar(Sys.which(config$plink_exec)))
    stop("PLINK2 executable not found: ", config$plink_exec)
  genes <- helpers$get_gene_annotations(config$gtf_file)
  genes <- genes[genes$gene_name %in% readLines(config$gene_list), , drop = FALSE]
  genes$chr <- sub("^chr", "", as.character(genes$chromosome))
  # Diploid autosomes only; sex chromosomes need a different generating model.
  genes <- genes[genes$chr %in% as.character(1:22), , drop = FALSE]
  genes <- genes[!duplicated(genes$gene_name), , drop = FALSE]
  if (!nrow(genes)) stop("No eligible autosomal protein-coding genes.")
  genes$tss <- ifelse(genes$strand == "+", genes$start, genes$end)
  genes
}

simulation_extract <- function(gene, config) {
  task_dir <- tempfile("plink_locus_")
  dir.create(task_dir)
  on.exit(unlink(task_dir, recursive = TRUE), add = TRUE)
  prefix <- file.path(task_dir, "genotype")
  args <- c("--bfile", shQuote(config$genotype_prefix), "--chr", gene$chr,
    "--from-bp", max(0, gene$tss - config$cis_window),
    "--to-bp", gene$tss + config$cis_window,
    "--snps-only", "--max-alleles", "2", "--rm-dup", "exclude-all",
    "--threads", config$plink_threads, "--memory", config$plink_memory_mb,
    "--maf", config$min_maf_plink, "--recode", "A", "--out", shQuote(prefix))
  log <- suppressWarnings(system2(config$plink_exec, as.character(args),
                                  stdout = TRUE, stderr = TRUE))
  status <- attr(log, "status")
  if ((!is.null(status) && status != 0L) || !file.exists(paste0(prefix, ".raw"))) {
    if (any(grepl("No variants remaining|No variants in", log, ignore.case = TRUE)))
      return(NULL)
    stop("PLINK extraction failed for ", gene$gene_name, ": ",
         paste(tail(log, 15), collapse = "\n"))
  }
  simulation_read_raw(paste0(prefix, ".raw"))
}

simulation_dataset <- function(scenario, config, helpers, catalogue, fixed_qc = NULL) {
  rejected <- list()
  for (attempt in seq_len(config$max_locus_attempts)) {
    if (!is.null(fixed_qc)) {
      gene <- data.frame(gene_name = "fixed_raw_pilot", chr = NA_character_, tss = NA_real_)
      qc <- fixed_qc
    } else {
      gene <- catalogue[sample.int(nrow(catalogue), 1L), , drop = FALSE]
      raw <- simulation_extract(gene, config)
      qc <- if (is.null(raw)) NULL else
        simulation_global_qc(raw, config, helpers, gene$gene_name)
    }
    reason <- "no_snps_after_global_qc"
    if (!is.null(qc)) {
      design <- simulation_design(qc$X, config, helpers)
      reason <- "insufficient_donors_or_retained_predictors"
      if (!is.null(design)) {
        causal <- simulation_causals(design, scenario, config)
        reason <- "insufficient_distinct_causal_snps_for_requested_codings"
        if (!is.null(causal)) {
          pheno <- simulation_phenotype(design, causal, scenario$pve)
          reason <- "zero_genetic_variance"
          if (!is.null(pheno)) return(list(
            gene = gene, design = design, phenotype = pheno,
            qc = c(qc$counts, design$counts), rejected = rejected))
        }
      }
    }
    rejected[[attempt]] <- data.frame(gene = gene$gene_name, reason = reason)
  }
  stop(structure(list(message = paste0("No eligible locus/sample after ",
    config$max_locus_attempts, " attempts. Last reason: ", reason),
    call = NULL, rejected_loci = rejected),
    class = c("simulation_locus_error", "error", "condition")))
}

simulation_snp_pip <- function(fit, map, snps, prior_tol = 1e-9) {
  active <- if (is.numeric(fit$V)) which(fit$V > prior_tol) else seq_len(nrow(fit$alpha))
  if (!length(active)) return(setNames(numeric(length(snps)), snps))
  alpha <- fit$alpha[active, , drop = FALSE]
  # Sum mutually exclusive codings WITHIN each single effect, then combine effects.
  # Taking max/sum of final predictor PIPs is not the SNP inclusion probability.
  if (!is.null(fit$slot_weights)) alpha <- alpha * fit$slot_weights[active]
  by_snp <- t(rowsum(t(alpha), group = factor(map$snp, levels = snps), reorder = TRUE))
  pip <- -expm1(colSums(log1p(-pmin(pmax(by_snp, 0), 1))))
  setNames(as.numeric(pip[match(snps, colnames(by_snp))]), snps)
}

simulation_purity <- function(X, indices, max_snps) {
  indices <- unique(indices)
  n <- length(indices)
  # Same random subsampling principle as susie_get_cs. Report when approximate.
  if (n > max_snps) indices <- indices[sample.int(n, max_snps)]
  purity <- if (length(indices) <= 1L) 1 else {
    corr <- abs(cor(X[, indices, drop = FALSE]))
    min(corr[upper.tri(corr)])
  }
  c(purity = purity, n_evaluated = length(indices), approximate = n > max_snps)
}

simulation_cs <- function(fit, X, map, additive, truth, config) {
  cs <- fit$sets$cs
  if (!length(cs)) return(data.frame())
  out <- lapply(seq_along(cs), function(i) {
    idx <- cs[[i]]
    snps <- unique(map$snp[idx])
    add_idx <- match(snps, colnames(additive))
    purity <- simulation_purity(additive, add_idx, config$n_purity)
    reported_purity <- fit$sets$purity$min.abs.corr[i]
    if (length(reported_purity) != 1L)
      reported_purity <- simulation_purity(X, idx, config$n_purity)["purity"]
    data.frame(cs_id = names(cs)[i],
      component = fit$sets$cs_index[i], size_snp = length(snps),
      size_predictor = length(idx), covered_snp = any(snps %in% truth$snp),
      covered_coding = any(paste(map$snp[idx], map$coding[idx]) %in%
                           paste(truth$snp, truth$coding)),
      purity_predictor = unname(reported_purity), purity_snp = unname(purity["purity"]),
      purity_snp_n = unname(purity["n_evaluated"]),
      purity_snp_approximate = as.logical(purity["approximate"]),
      purity_predictor_approximate = length(idx) > config$n_purity,
      snps = paste(snps, collapse = ";"),
      predictors = paste(map$predictor[idx], collapse = ";"), stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

simulation_calibration <- function(pip, truth, breaks) {
  bin <- cut(pip, breaks, include.lowest = TRUE, right = TRUE, labels = FALSE)
  do.call(rbind, lapply(seq_len(length(breaks) - 1L), function(i) {
    take <- which(bin == i)
    data.frame(bin = i, lower = breaks[i], upper = breaks[i + 1L],
      n = length(take), sum_pip = sum(pip[take]), n_causal = sum(truth[take]),
      sum_squared_error = sum((pip[take] - truth[take])^2))
  }))
}

simulation_evaluate <- function(fit, X, map, design, pheno, config) {
  snps <- colnames(design$add)
  pip <- simulation_snp_pip(fit, map, snps, config$prior_tol)
  if (any(!is.finite(pip)) || any(pip < 0 | pip > 1)) stop("Invalid SNP PIPs.")
  causal <- snps %in% pheno$truth$snp
  cs <- simulation_cs(fit, X, map, design$add, pheno$truth, config)
  cs_snps <- unique(map$snp[unlist(fit$sets$cs, use.names = FALSE)])
  called <- pip >= config$pip_threshold
  tp <- sum(causal & called)
  fp <- sum(!causal & called)
  n_cs <- nrow(cs)
  safe_mean <- function(x) if (length(x)) mean(x) else NA_real_
  metrics <- data.frame(n_snps = length(snps), n_causal = sum(causal),
    n_cs = n_cs, tp = tp, fp = fp, discoveries = tp + fp,
    power = tp / sum(causal), fdp = fp / max(1, tp + fp),
    cs_power = sum(pheno$truth$snp %in% cs_snps) / sum(causal),
    cs_coverage = safe_mean(cs$covered_snp),
    cs_coding_coverage = safe_mean(cs$covered_coding),
    cs_size_snp = safe_mean(cs$size_snp), cs_size_predictor = safe_mean(cs$size_predictor),
    cs_purity_snp = safe_mean(cs$purity_snp),
    cs_purity_predictor = safe_mean(cs$purity_predictor),
    brier_score = mean((pip - causal)^2))
  curve <- do.call(rbind, lapply(config$thresholds, function(t) {
    called <- pip >= t
    tp <- sum(called & causal)
    fp <- sum(called & !causal)
    data.frame(threshold = t, tp = tp, fp = fp, discoveries = tp + fp,
      power = tp / sum(causal), fdp = fp / max(1, tp + fp))
  }))
  list(metrics = metrics, cs = cs, curve = curve,
    calibration = simulation_calibration(pip, causal, config$calibration_breaks),
    snp_pip = data.frame(snp = snps, pip = unname(pip), causal = causal),
    predictor_pip = data.frame(map, pip = as.numeric(susieR::susie_get_pip(
      fit, prune_by_cs = FALSE, prior_tol = config$prior_tol))))
}

simulation_fit <- function(dataset, config) {
  design <- dataset$design
  add_map <- data.frame(predictor = colnames(design$add), snp = colnames(design$add),
                         coding = "additive", stringsAsFactors = FALSE)
  out <- list()
  for (method in c("susie", "susie_mix")) {
    X <- if (method == "susie") design$add else design$mix
    map <- if (method == "susie") add_map else design$map
    warnings <- character()
    started <- proc.time()[["elapsed"]]
    result <- tryCatch(withCallingHandlers({
      fit <- susieR::susie(X, dataset$phenotype$y, L = config$L,
        standardize = config$standardize, estimate_prior_method = config$estimate_prior_method,
        min_abs_corr = config$min_abs_corr, coverage = config$coverage,
        max_iter = config$max_iter, tol = config$tol, prior_tol = config$prior_tol,
        n_purity = config$n_purity, verbose = FALSE)
      evaluated <- simulation_evaluate(fit, X, map, design, dataset$phenotype, config)
      evaluated$converged <- isTRUE(fit$converged)
      evaluated$niter <- fit$niter
      evaluated$fit <- if (config$save_fits) fit else NULL
      evaluated
    }, warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) list(error = conditionMessage(e), converged = FALSE))
    result$seconds <- proc.time()[["elapsed"]] - started
    result$warnings <- unique(warnings)
    out[[method]] <- result
  }
  out
}

simulation_save <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile("checkpoint_", tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(object, tmp)
  # On Windows rename cannot replace an existing file. Existing successes are
  # never rewritten; failed trials are rewritten only with --retry-failed.
  if (file.exists(path) && !file.remove(path)) stop("Cannot replace ", path)
  if (!file.rename(tmp, path)) stop("Cannot finalize checkpoint: ", path)
  invisible(path)
}
