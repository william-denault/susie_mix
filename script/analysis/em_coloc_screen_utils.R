# First-stage eQTL screening only: no GWAS selection, fitting, or colocalization.
em_screen_lse <- function(x) {
  m <- max(x)
  if (m == -Inf) return(-Inf)
  m + log(sum(exp(x - m)))
}

em_screen_map <- function(fit, map) {
  needed <- c("predictor_index", "predictor_name", "snp", "coding")
  if (!is.matrix(fit$alpha) || !is.data.frame(map) || !all(needed %in% names(map)) ||
      nrow(map) != ncol(fit$alpha) ||
      !identical(as.integer(map$predictor_index), seq_len(nrow(map))) ||
      anyNA(map[needed]) || any(!nzchar(map$snp)) || any(!nzchar(map$predictor_name)) ||
      anyDuplicated(map$predictor_name) || anyDuplicated(map[c("snp", "coding")]) ||
      !all(map$coding %in% c("additive", "recessive", "dominant"))) {
    stop("Invalid or misaligned predictor map.")
  }
  if (!is.null(colnames(fit$alpha)) && !identical(colnames(fit$alpha), map$predictor_name))
    stop("alpha column names disagree with predictor map.")
  if (!isTRUE(fit$converged) || any(!is.finite(fit$alpha)) || any(fit$alpha < 0) ||
      any(abs(rowSums(fit$alpha) - 1) > 1e-6)) stop("Fit is not converged or has invalid alpha.")
  invisible(TRUE)
}

em_screen_cs <- function(fit, map) {
  em_screen_map(fit, map)
  cs <- get_cs(fit)
  components <- fit$sets$cs_index
  if (length(cs) && (length(components) != length(cs) || anyNA(components) ||
      any(components != as.integer(components)) || anyDuplicated(components) ||
      any(components < 1 | components > nrow(fit$alpha)))) stop("Invalid CS component indices.")
  coords <- get_predictor_coordinates(map)
  lapply(seq_along(cs), function(i) {
    members <- cs[[i]]
    if (!length(members) || anyNA(members) || any(members != as.integer(members)) ||
        any(members < 1 | members > nrow(map))) stop("Invalid CS members.")
    lead <- get_cs_lead_predictor(fit, cs, i)
    component <- as.integer(components[i])
    list(cs = i, component = component, lead = lead, snp = map$snp[lead],
         coding = map$coding[lead], chromosome = sub("^chr", "", coords$chromosome[lead]),
         position = coords$position[lead], snps = unique(map$snp[members]),
         alpha = fit$alpha[component, lead],
         nonadditive_mass = sum(fit$alpha[component, map$coding != "additive"]),
         cs_mass = sum(fit$alpha[component, unique(members)]))
  })
}

# GTEx v8 IDs retain REF/ALT and the PLINK counted allele. Never join on position alone.
em_screen_variants <- function(snps) {
  parts <- strsplit(snps, "_", fixed = TRUE)
  valid <- vapply(parts, function(x) length(x) %in% c(5L, 6L) &&
    grepl("^(chr)?([0-9]+|X|Y|M|MT)$", x[1]) && grepl("^[0-9]+$", x[2]) &&
    grepl("^[ACGT]+$", x[3]) && grepl("^[ACGT]+$", x[4]) && x[5] == "b38" &&
    (length(x) == 5L || x[6] %in% x[3:4]), logical(1))
  if (!all(valid)) stop("Cannot identify GRCh38 REF/ALT from SNP IDs: ",
                        paste(head(snps[!valid], 3), collapse = ", "))
  at <- function(i) vapply(parts, `[`, character(1), i)
  out <- data.frame(source_snp = snps, chromosome = sub("^chr", "", at(1)),
                    position = as.numeric(at(2)), ref = at(3), alt = at(4),
                    counted_allele = vapply(parts, function(x) if (length(x) == 6L) x[6] else NA_character_, ""),
                    build = "GRCh38", stringsAsFactors = FALSE)
  out$snp <- paste(out$chromosome, format(out$position, scientific = FALSE, trim = TRUE), out$ref, out$alt, sep = ":")
  if (anyDuplicated(out$snp) || any(!is.finite(out$position)) || any(out$position <= 0))
    stop("Duplicate or invalid biological variant IDs.")
  out
}

em_screen_bfs <- function(fit, map, weights = NULL) {
  em_screen_map(fit, map)
  bf <- fit$lbf_variable
  if (!is.matrix(bf) || !identical(dim(bf), dim(fit$alpha)) || anyNA(bf) || any(bf == Inf))
    stop("Missing or invalid lbf_variable; BFs cannot be reconstructed from PIPs.")
  if (!is.null(colnames(bf)) && !identical(colnames(bf), map$predictor_name))
    stop("lbf_variable column names disagree with predictor map.")
  if (is.null(weights)) weights <- fit$pi
  if (is.null(weights)) stop("Saved predictor priors are absent; refusing to guess BF mixture weights.")
  if (length(weights) != nrow(map) || any(!is.finite(weights)) || any(weights < 0) || sum(weights) <= 0)
    stop("Invalid predictor priors.")
  if (!is.null(names(weights)) && !identical(names(weights), map$predictor_name))
    stop("Predictor prior names disagree with predictor map.")
  weights <- unname(weights / sum(weights))
  if (!is.null(fit$pi) && !isTRUE(all.equal(unname(fit$pi), weights, tolerance = 1e-10)))
    stop("Saved weights disagree with fit$pi.")
  # This also detects stale BFs, wrong predictor order, and prior double-counting.
  joint <- sweep(bf, 2, log(weights), "+")
  norm <- apply(joint, 1, em_screen_lse)
  if (any(!is.finite(norm)) || max(abs(exp(joint - norm) - fit$alpha)) > 1e-6)
    stop("Saved alpha is inconsistent with lbf_variable and predictor priors.")
  snps <- unique(map$snp)
  groups <- split(seq_len(nrow(map)), factor(map$snp, levels = snps))
  q <- vapply(groups, function(ii) sum(weights[ii]), numeric(1))
  supported <- q > 0
  groups <- groups[supported]
  q <- q[supported]
  variants <- em_screen_variants(snps[supported])
  lbf <- matrix(NA_real_, nrow(bf), length(groups))
  alpha <- lbf
  for (j in seq_along(groups)) {
    ii <- groups[[j]]
    # Marginal likelihood conditional on biological SNP J=j, NOT its posterior.
    lbf[, j] <- apply(sweep(bf[, ii, drop = FALSE], 2, log(weights[ii] / q[j]), "+"),
                      1, em_screen_lse)
    alpha[, j] <- rowSums(fit$alpha[, ii, drop = FALSE])
  }
  dimnames(lbf) <- dimnames(alpha) <- list(paste0("L", seq_len(nrow(bf))), variants$snp)
  names(q) <- variants$snp
  components <- as.integer(fit$sets$cs_index)
  list(lbf = lbf[components, , drop = FALSE], all_component_lbf = lbf,
       alpha = alpha[components, , drop = FALSE], snp_prior = q, variants = variants,
       components = components, zero_prior_snps = snps[!supported])
}

em_screen_bundle <- function(x, em_cs, vanilla_cs, metadata) {
  em <- em_screen_bfs(x$weighted_fit_mix, x$mix_predictor_map, x$weighted_mix_prior_weights)
  vanilla <- em_screen_bfs(x$susie_add, x$add_predictor_map)
  common <- intersect(colnames(em$lbf), colnames(vanilla$lbf))
  if (!length(common)) stop("No common biological SNPs between EM and vanilla.")
  # Keep native matrices as well: intersecting SNPs later with GWAS loses more mass.
  paired <- list(snp = common, em_lbf = em$lbf[, common, drop = FALSE],
                 vanilla_lbf = vanilla$lbf[, common, drop = FALSE],
                 em_retained_mass = rowSums(em$alpha[, common, drop = FALSE]),
                 vanilla_retained_mass = rowSums(vanilla$alpha[, common, drop = FALSE]))
  list(schema_version = 1L, metadata = metadata, em = em, vanilla = vanilla,
       em_cs = em_cs, vanilla_cs = vanilla_cs, paired = paired,
       notes = c("Natural-log BFs; rows are CS components; columns are biological SNPs.",
                 "SNP priors are separate from BFs. Declare coloc priors explicitly.",
                 "Harmonize GWAS build/alleles and intersect all THREE datasets before paired coloc.",
                 "all_component_lbf includes components without a CS for diagnostics only."))
}

em_screen_write_bfs <- function(bundle, path) {
  con <- gzfile(path, "wt")
  on.exit(close(con))
  writeLines("method\tcomponent\tsnp\tlog_bf\tsnp_prior\tcomponent_alpha", con)
  for (method in c("em", "vanilla")) {
    b <- bundle[[method]]
    for (i in seq_len(nrow(b$lbf))) {
      write.table(data.frame(method, component = b$components[i], snp = colnames(b$lbf),
                             log_bf = b$lbf[i, ], snp_prior = unname(b$snp_prior),
                             component_alpha = b$alpha[i, ]),
                  con, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
  }
}

# Lazy, one-gene cache. Reuse the exact GTEx sample/predictor reconstruction used by plots.
em_screen_ld_reader <- function(project_dir, datadir, output_dir) {
  source(file.path(project_dir, "script/analysis/fit_mix_expression_plot_utils.R"), local = TRUE)
  shared <- NULL
  cached_gene <- NULL
  gene_data <- NULL
  temp_dir <- file.path(output_dir, "temp_plink")
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
  function(gene, tissue, x) {
    if (is.null(shared)) shared <<- load_fit_mix_expression_inputs(
      subject_pheno_file = file.path(datadir, "GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt.gz"),
      sample_attr_file = file.path(datadir, "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt.gz"),
      expr_file = file.path(datadir, "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz"),
      gtf_file = file.path(datadir, "Homo_sapiens.GRCh38.103.chr.reformatted.collapse_only.gene.gtf.gz"),
      gene_annotation_function = file.path(project_dir, "script/scan_tissue_attempt/get_gene_annotations.R"))
    if (!identical(gene, cached_gene)) {
      gene_data <<- prepare_fit_mix_gene_data(gene_name = gene, shared_inputs = shared,
        results_dir = file.path(project_dir, "results"), temp_dir = temp_dir,
        plink_exec = file.path(datadir, "plink2"),
        geno_file = file.path(datadir, "GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv"),
        cis_window = 5e5, min_maf_plink = 0, min_maf = .05, hwe_thresh = 1e-8)
      cached_gene <<- gene
    }
    x$susie_mix <- x$weighted_fit_mix
    prepared <- prepare_fit_mix_tissue_data(gene_data, gene, tissue, x,
                                            min_samples = 0L, min_n_rec = 5)
    prepared$geno_for_counts
  }
}

em_screen_ld <- function(genotypes, em, vanilla, n_ind) {
  snps <- unique(c(vapply(em, `[[`, "", "snp"), vapply(vanilla, `[[`, "", "snp")))
  if (!is.matrix(genotypes) || is.null(colnames(genotypes)) ||
      anyDuplicated(colnames(genotypes)) || !all(snps %in% colnames(genotypes)) ||
      nrow(genotypes) < 3 || length(n_ind) != 1L || nrow(genotypes) != n_ind)
    stop("LD genotypes do not match saved leads/sample count.")
  g <- genotypes[, snps, drop = FALSE]
  if (any(!is.finite(g)) || any(!g %in% 0:2) || any(apply(g, 2, sd) == 0))
    stop("LD requires complete, polymorphic 0/1/2 genotypes from fitted samples.")
  g
}

em_screen_bind <- function(rows, columns) {
  if (!length(rows)) return(as.data.frame(setNames(rep(list(character()), length(columns)), columns)))
  columns <- unique(c(columns, unlist(lapply(rows, names))))
  do.call(rbind, lapply(rows, function(row) {
    for (nm in setdiff(columns, names(row))) row[[nm]] <- NA
    row[columns]
  }))
}

em_screen_coloc <- function(project_dir, iteration = "latest", datadir = "/project2/mstephens/gtex",
                            association_threshold = 5e-8, shift_bp = 10000,
                            compute_ld = TRUE, genotype_reader = NULL, export_tsv = TRUE,
                            output_dir = NULL) {
  if (length(association_threshold) != 1L || !is.finite(association_threshold) ||
      association_threshold <= 0 || association_threshold > 1 || length(shift_bp) != 1L ||
      !is.finite(shift_bp) || shift_bp < 0) stop("Invalid screening thresholds.")
  context <- em_analysis_context(project_dir, iteration)
  source(file.path(context$project_dir, "script/analysis/generate_summary_results.R"), local = TRUE)
  join <- em_join_gene
  environment(join) <- environment()
  cs_info <- em_screen_cs
  environment(cs_info) <- environment()
  # Each run has its own manifest and exports; no stale bundles from earlier filters.
  if (is.null(output_dir)) {
    parent <- file.path(context$iteration_dir, "coloc_screen")
    dir.create(parent, recursive = TRUE, showWarnings = FALSE)
    output_dir <- tempfile(paste0("screen_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"), tmpdir = parent)
  }
  if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE)))
    stop("output_dir must be new or empty to prevent stale exports.")
  dir.create(file.path(output_dir, "bf"), recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  metadata <- list(iteration = context$iteration, iteration_dir = context$iteration_dir,
    baseline_dir = context$baseline_dir, priors = context$priors,
    association_threshold = association_threshold, shift_bp = shift_bp, compute_ld = compute_ld,
    p_value_definition = "Saved minimum marginal additive cis-SNP P value (min_pv); gene-tissue screen",
    lead_definition = "Highest component alpha predictor within each saved credible set",
    ld_definition = "Squared Pearson correlation of raw additive dosages in reconstructed tissue samples",
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    session_info = sessionInfo())
  saveRDS(metadata, file.path(output_dir, "screen_metadata.rds"))
  write.csv(context$priors, file.path(output_dir, "em_priors_used.csv"), row.names = FALSE)
  reader_error <- NULL
  if (compute_ld && is.null(genotype_reader)) {
    genotype_reader <- tryCatch(em_screen_ld_reader(context$project_dir, datadir, output_dir), error = identity)
    if (inherits(genotype_reader, "error")) reader_error <- conditionMessage(genotype_reader)
  }
  screen <- pairs <- audit <- exports <- list()
  record_audit <- function(gene, tissue, stage, message) {
    audit[[length(audit) + 1L]] <<- data.frame(gene, tissue, stage, message)
  }
  counts <- c(n_gene_files = length(context$result_files), n_tissues = 0L,
              n_p_pass = 0L, n_no_em_cs = 0L)
  message("Screening saved EM fits; summaries are not required. Output: ", output_dir)
  for (file_index in seq_along(context$result_files)) {
    file <- context$result_files[file_index]
    gene <- sub("\\.rds$", "", basename(file))
    joined <- tryCatch(join(file, context), error = identity)
    if (inherits(joined, "error") || !is.null(joined$error)) {
      record_audit(gene, NA_character_, "join", if (inherits(joined, "error")) conditionMessage(joined) else joined$error)
      next
    }
    for (tissue in names(attr(joined, "tissue_errors")))
      record_audit(gene, tissue, "saved_tissue_error", as.character(attr(joined, "tissue_errors")[[tissue]]))
    for (tissue in names(joined)) {
      x <- joined[[tissue]]
      counts["n_tissues"] <- counts["n_tissues"] + 1L
      p <- x$min_pv
      if (length(p) != 1L || !is.finite(p) || p < 0 || p > 1) {
        record_audit(gene, tissue, "p_value", "Missing or invalid saved min_pv")
        next
      }
      if (p >= association_threshold) next
      counts["n_p_pass"] <- counts["n_p_pass"] + 1L
      compared <- tryCatch({
        if (!all(x$add_predictor_map$coding == "additive")) stop("Vanilla predictor map is not additive.")
        list(em = cs_info(x$weighted_fit_mix, x$mix_predictor_map),
             vanilla = cs_info(x$susie_add, x$add_predictor_map))
      }, error = identity)
      if (inherits(compared, "error")) {
        record_audit(gene, tissue, "cs", conditionMessage(compared))
        next
      }
      em <- compared$em
      vanilla <- compared$vanilla
      if (!length(em)) {
        counts["n_no_em_cs"] <- counts["n_no_em_cs"] + 1L
        next
      }
      local_screen <- local_pairs <- list()
      for (e in em) {
        pair_rows <- lapply(vanilla, function(v) {
          valid <- !anyNA(c(e$chromosome, v$chromosome)) &&
            e$chromosome == v$chromosome && all(is.finite(c(e$position, v$position))) &&
            min(e$position, v$position) > 0
          shared <- length(intersect(e$snps, v$snps))
          data.frame(gene, tissue, em_cs = e$cs, em_component = e$component,
            vanilla_cs = v$cs, vanilla_component = v$component,
            em_lead_snp = e$snp, vanilla_lead_snp = v$snp,
            distance_bp = if (valid) abs(e$position - v$position) else NA_real_,
            n_shared_snps = shared, jaccard = shared / length(union(e$snps, v$snps)),
            em_cs_overlap_fraction = shared / length(e$snps),
            vanilla_cs_overlap_fraction = shared / length(v$snps),
            lead_r2 = NA_real_, coded_lead_r2 = NA_real_, ld_status = "not_requested_for_unshifted")
        })
        pr <- if (length(pair_rows)) do.call(rbind, pair_rows) else NULL
        distance_ok <- length(vanilla) > 0L && all(is.finite(pr$distance_bp))
        shifted <- distance_ok && all(pr$distance_bp > shift_bp)
        classification <- if (!length(vanilla)) "no_vanilla_cs" else if (!distance_ok)
          "invalid_coordinates" else if (shifted) "shifted_from_all_vanilla_leads" else "not_shifted"
        local_screen[[length(local_screen) + 1L]] <- data.frame(gene, tissue,
          em_iteration = context$iteration, min_pv = p, em_cs = e$cs, em_component = e$component,
          em_lead_snp = e$snp, em_lead_coding = e$coding, em_lead_alpha = e$alpha,
          em_nonadditive_mass = e$nonadditive_mass, em_cs_mass = e$cs_mass,
          em_cs_snps = length(e$snps), n_vanilla_cs = length(vanilla),
          min_distance_bp = if (distance_ok) min(pr$distance_bp) else NA_real_,
          max_cs_jaccard = if (length(vanilla)) max(pr$jaccard) else NA_real_,
          max_shared_snps = if (length(vanilla)) max(pr$n_shared_snps) else NA_integer_,
          shifted = shifted, classification = classification, max_lead_r2 = NA_real_,
          n_ind = if (length(x$n_ind) == 1L) x$n_ind else NA_integer_,
          mean_read = if (length(x$mean_read) == 1L) x$mean_read else NA_real_,
          lead_n_0 = NA_integer_, lead_n_1 = NA_integer_, lead_n_2 = NA_integer_,
          ld_status = if (!length(vanilla)) "no_vanilla_cs" else "not_requested_for_unshifted",
          bf_status = "not_exported", bundle = NA_character_)
        local_pairs <- c(local_pairs, pair_rows)
      }
      s <- do.call(rbind, local_screen)
      pr <- if (length(local_pairs)) do.call(rbind, local_pairs) else NULL
      export_case <- any(s$shifted) || !length(vanilla)
      if (export_case) {
        id <- sprintf("case_%06d", length(exports) + 1L)
        relative <- paste0("bf/", id, ".rds")
        bundle <- tryCatch(em_screen_bundle(x, em, vanilla,
          c(metadata[c("iteration", "association_threshold", "shift_bp")], list(gene = gene, tissue = tissue, min_pv = p))), error = identity)
        if (inherits(bundle, "error")) {
          s$bf_status <- "error"
          record_audit(gene, tissue, "bf_export", conditionMessage(bundle))
        } else {
          saveRDS(bundle, file.path(output_dir, relative))
          if (export_tsv) em_screen_write_bfs(bundle, file.path(output_dir, "bf", paste0(id, ".tsv.gz")))
          s$bf_status <- if (length(vanilla)) "exported" else "exported_no_vanilla_cs"
          s$bundle <- relative
        }
        exports[[length(exports) + 1L]] <- data.frame(gene, tissue, bundle = s$bundle[1],
          bf_status = s$bf_status[1], n_em_cs = length(em), n_vanilla_cs = length(vanilla),
          n_common_snps = if (inherits(bundle, "error")) NA_integer_ else length(bundle$paired$snp),
          min_em_retained_mass = if (inherits(bundle, "error")) NA_real_ else min(bundle$paired$em_retained_mass),
          min_vanilla_retained_mass = if (inherits(bundle, "error") || !length(vanilla)) NA_real_ else
            min(bundle$paired$vanilla_retained_mass))
        if (length(vanilla)) {
          ld <- if (!compute_ld) NULL else tryCatch({
            if (!is.null(reader_error)) stop(reader_error)
            em_screen_ld(genotype_reader(gene, tissue, x), em, vanilla, x$n_ind)
          }, error = identity)
          status <- if (!compute_ld) "disabled" else if (inherits(ld, "error")) "error" else "computed"
          s$ld_status <- pr$ld_status <- status
          if (inherits(ld, "error")) record_audit(gene, tissue, "ld", conditionMessage(ld))
          if (identical(status, "computed")) {
            for (i in seq_len(nrow(pr))) {
              eg <- ld[, pr$em_lead_snp[i]]
              vg <- ld[, pr$vanilla_lead_snp[i]]
              coding <- em[[pr$em_cs[i]]]$coding
              coded <- switch(coding, additive = eg, dominant = as.numeric(eg > 0), recessive = as.numeric(eg == 2))
              pr$lead_r2[i] <- cor(eg, vg)^2
              pr$coded_lead_r2[i] <- if (sd(coded) > 0) cor(coded, vg)^2 else NA_real_
            }
            s$max_lead_r2 <- vapply(s$em_cs, function(i) max(pr$lead_r2[pr$em_cs == i]), numeric(1))
            for (k in 0:2) s[[paste0("lead_n_", k)]] <- vapply(s$em_lead_snp,
              function(snp) sum(ld[, snp] == k), integer(1))
          }
        }
      }
      screen[[length(screen) + 1L]] <- s
      if (!is.null(pr)) pairs[[length(pairs) + 1L]] <- pr
    }
    if (file_index %% 100L == 0L) message("Screened ", file_index, "/", length(context$result_files), " gene files")
  }
  s <- em_screen_bind(screen, c("gene", "tissue", "em_cs", "shifted", "classification", "bf_status", "ld_status"))
  pr <- em_screen_bind(pairs, c("gene", "tissue", "em_cs", "vanilla_cs", "distance_bp", "lead_r2", "jaccard"))
  a <- em_screen_bind(audit, c("gene", "tissue", "stage", "message"))
  ex <- em_screen_bind(exports, c("gene", "tissue", "bundle", "bf_status", "n_common_snps"))
  outputs <- list(all_em_cs = s, candidates = s[which(s$shifted %in% TRUE), , drop = FALSE],
    no_vanilla_cs = s[s$classification == "no_vanilla_cs", , drop = FALSE],
    cs_pair_comparisons = pr, audit = a, bf_manifest = ex)
  for (nm in names(outputs)) write.csv(outputs[[nm]], file.path(output_dir, paste0(nm, ".csv")), row.names = FALSE)
  stats <- c(as.list(counts), list(n_em_cs = nrow(s), n_shifted_cs = nrow(outputs$candidates),
    n_no_vanilla_cs = nrow(outputs$no_vanilla_cs), n_exported_gene_tissues = sum(grepl("^exported", ex$bf_status)),
    n_audit_records = nrow(a)))
  write.csv(as.data.frame(stats), file.path(output_dir, "screen_counts.csv"), row.names = FALSE)
  saveRDS(c(metadata, stats, list(completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE))),
          file.path(output_dir, "screen_metadata.rds"))
  print(as.data.frame(stats))
  message("Screen complete: ", output_dir)
  if (nrow(a)) warning(nrow(a), " audit records; inspect audit.csv before interpreting the screen.", call. = FALSE)
  invisible(c(outputs, list(output_dir = output_dir, counts = stats)))
}
