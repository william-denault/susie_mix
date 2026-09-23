# Post-processing of workhorse.R's fit_slide / fit_slide_perm; no refitting.
# Jobs and other source wrappers add call frames, so frame 1 need not belong
# to source(). Find this helper's innermost source frame instead.
.slide_source_dir <- (function() {
  for (frame in rev(sys.frames())) {
    path <- get0("ofile", envir = frame, inherits = FALSE)
    if (is.character(path) && length(path) == 1L && !is.na(path) &&
        basename(path) == "slide_analysis_utils.R") {
      return(dirname(normalizePath(path, winslash = "/", mustWork = TRUE)))
    }
  }
  stop("Cannot locate slide_analysis_utils.R; load it with source() using its full path.")
})()
source(file.path(.slide_source_dir, "generate_summary_results.R"), local = TRUE)
source(file.path(.slide_source_dir, "tss_disagreement_utils.R"), local = TRUE)
source(file.path(.slide_source_dir, "one_cs_lead_distance_plot_utils.R"), local = TRUE)
rm(.slide_source_dir)

slide_coding <- function(delta, tolerance = 1e-6) {
  if (length(tolerance) != 1L || !is.finite(tolerance) || tolerance <= 0 || tolerance >= .5)
    stop("delta_tolerance must be between 0 and 0.5.")
  if (any(is.finite(delta) & abs(delta) > 1 + tolerance)) stop("Slide delta outside [-1, 1].")
  out <- rep("unresolved", length(delta))
  out[is.finite(delta) & delta > tolerance] <- "partial_dominant"
  out[is.finite(delta) & delta < -tolerance] <- "partial_recessive"
  out[is.finite(delta) & abs(delta) <= tolerance] <- "additive"
  out[is.finite(delta) & abs(delta - 1) <= tolerance] <- "dominant"
  out[is.finite(delta) & abs(delta + 1) <= tolerance] <- "recessive"
  out
}

slide_validate_fit <- function(fit, map, slide = FALSE) {
  if (!is.list(fit) || !is.matrix(fit$alpha)) stop("Missing fit/alpha matrix.")
  validate_predictor_map(fit, map, "add_predictor_map")
  if (anyNA(map$snp) || any(!nzchar(map$snp)) || anyDuplicated(map$snp))
    stop("Biological SNP identifiers must be unique and nonempty.")
  if (!is.null(colnames(fit$alpha)) && !identical(colnames(fit$alpha), as.character(map$predictor_name)))
    stop("Fit predictor names/order differ from add_predictor_map.")
  if (length(fit$pip) != nrow(map) || any(!is.finite(fit$pip)) || any(fit$pip < 0 | fit$pip > 1))
    stop("Invalid SNP PIP vector.")
  get_cs_snp_sets(fit, map) # Check all saved CS member indices.
  cs <- get_cs(fit)
  if (length(cs) && (length(fit$sets$cs_index) != length(cs) ||
      anyNA(fit$sets$cs_index) || any(fit$sets$cs_index < 1 | fit$sets$cs_index > nrow(fit$alpha) |
                                     fit$sets$cs_index != floor(fit$sets$cs_index))))
    stop("Missing or invalid CS component indices.")
  if (any(lengths(cs) == 0L)) stop("Empty reported credible set.")
  if (slide && (!is.matrix(fit$delta) || !identical(dim(fit$delta), dim(fit$alpha)) ||
                any(!is.finite(fit$delta)))) stop("Slide delta must be a finite L-by-p matrix aligned with alpha.")
  if (slide) slide_coding(as.vector(fit$delta))
  if (slide && !is.null(fit$delta_forced) && length(fit$delta_forced) != nrow(map))
    stop("delta_forced is not aligned with the SNP map.")
  invisible(TRUE)
}

slide_leads <- function(fit, map, slide = FALSE, delta_tolerance = 1e-6) {
  cs <- get_cs(fit)
  index <- vapply(seq_along(cs), function(i) get_cs_lead_predictor(fit, cs, i), integer(1))
  component <- as.integer(fit$sets$cs_index)
  delta <- if (slide) fit$delta[cbind(component, index)] else rep(0, length(index))
  data.frame(cs_number = seq_along(cs), component = component, index = index,
             lead_snp = as.character(map$snp[index]), lead_pip = fit$pip[index],
             lead_delta = delta, lead_coding = slide_coding(delta, delta_tolerance),
             lead_delta_forced = if (slide && !is.null(fit$delta_forced))
               as.logical(fit$delta_forced[index]) else rep(NA, length(index)))
}

slide_bind_rows <- function(rows, empty = data.frame()) {
  rows <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, rows)
  if (!length(rows)) return(empty)
  columns <- unique(unlist(lapply(rows, names)))
  do.call(rbind, lapply(rows, function(x) {
    for (nm in setdiff(columns, names(x))) x[[nm]] <- NA
    x[columns]
  }))
}

slide_gene_annotations <- function(project_dir, gtf_file) {
  if (is.null(gtf_file) || !nzchar(gtf_file) || !file.exists(gtf_file)) return(NULL)
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Reading the GTF requires data.table.")
  fread <- data.table::fread
  source(file.path(project_dir, "script/scan_tissue_attempt/get_gene_annotations.R"), local = TRUE)
  annot <- get_gene_annotations(gtf_file)
  annot[!duplicated(annot$gene_name), , drop = FALSE] # Same first-row rule as the workhorse.
}

slide_tss_metadata <- function(x, gene, annotations = NULL) {
  info <- get_tss_metadata(x, x$add_predictor_map)
  if (is.finite(info$tss_position)) return(info)
  if (!is.null(annotations) && gene %in% annotations$gene_name) {
    a <- annotations[match(gene, annotations$gene_name), ]
    return(list(tss_position = if (a$strand == "-") a$end else a$start,
                chromosome = a$chromosome, strand = a$strand, tss_source = "workhorse_GTF"))
  }
  # Current workhorse_utils.R saves named SIGNED genomic offsets (pos - TSS).
  # Use the SNP names on these vectors, not CS order or newly computed leads.
  values <- unlist(lapply(c("susie_add", "susie_mix", "weighted_fit_mix"), function(key) {
    v <- x[[paste0(key, "_lead_snp_tss_distance")]]
    if (!is.numeric(v) || !is.null(dim(v)) || is.null(names(v))) return(numeric())
    pos <- parse_variant_position(names(v))
    (pos - v)[is.finite(pos) & is.finite(v)]
  }), use.names = FALSE)
  if (length(values) && all(values > 0) && diff(range(values)) <= 1) {
    info$tss_position <- values[1]
    info$tss_source <- "recovered_from_workhorse_signed_offsets"
  } else if (length(values)) info$tss_source <- "inconsistent_saved_offsets"
  info
}

slide_summarize_tissue <- function(x, gene, tissue, file, annotations = NULL, delta_tolerance = 1e-6) {
  map <- x$add_predictor_map
  slide_validate_fit(x[["susie_add"]], map)
  slide_validate_fit(x[["fit_slide"]], map, TRUE)
  add <- x$susie_add
  fit <- x$fit_slide
  leads <- slide_leads(fit, map, TRUE, delta_tolerance)
  add_sets <- get_cs_snp_sets(add, map)
  slide_sets <- get_cs_snp_sets(fit, map)
  add_snps <- unique(unlist(add_sets, use.names = FALSE))
  slide_snps <- unique(unlist(slide_sets, use.names = FALSE))
  shared <- length(intersect(add_snps, slide_snps))
  comparison <- compare_mixed_fits(add, fit, map)
  row <- data.frame(gene = gene, tissue = tissue, result_file = basename(file),
    n_ind = as.numeric(value_or_na(x$n_ind)), n_SNP = nrow(map),
    min_pv = as.numeric(value_or_na(x$min_pv)), mean_count = as.numeric(value_or_na(x$mean_read)),
    median_count = as.numeric(value_or_na(x$median_read)), ncs_susie = length(add_sets), ncs_slide = length(slide_sets),
    converged_susie = as.logical(value_or_na(add$converged)), converged_slide = as.logical(value_or_na(fit$converged)),
    elbo_susie = max_elbo(add), elbo_slide = max_elbo(fit), dif_elbo_slide_vs_add = max_elbo(fit) - max_elbo(add),
    log_lik_add = get_log_lik_metric(add), log_lik_slide = get_log_lik_metric(fit),
    n_add_cs_snps = length(add_snps), n_slide_cs_snps = length(slide_snps), overlap_snp = shared,
    jaccard_snp = if (length(union(add_snps, slide_snps))) shared / length(union(add_snps, slide_snps)) else NA_real_,
    pct_add_cs_snps_retained = safe_percentage(shared, length(add_snps)),
    same_cs_count = length(add_sets) == length(slide_sets),
    same_cs_snp_sets = comparison$mix_weighted_same_cs_snp_sets,
    same_lead_snp_set = comparison$mix_weighted_same_lead_snp_set,
    delta_tolerance = delta_tolerance)
  levels <- c("additive", "dominant", "recessive", "partial_dominant", "partial_recessive")
  for (coding in levels) row[[paste0("n_", coding, "_slide")]] <- sum(leads$lead_coding == coding)
  row$n_forced_additive_slide <- sum(leads$lead_delta_forced %in% TRUE)
  perm_errors <- character()
  for (key in c("susie_add_perm", "fit_slide_perm")) {
    is_slide <- key == "fit_slide_perm"
    perm <- x[[key]]
    error <- tryCatch({ slide_validate_fit(perm, map, is_slide); NULL }, error = conditionMessage)
    label <- if (is_slide) "slide" else "susie"
    row[[paste0("has_", key)]] <- is.null(error)
    row[[paste0("perm_cs_", label)]] <- if (is.null(error)) length(get_cs(perm)) else NA_integer_
    row[[paste0("elbo_", label, "_perm")]] <- if (is.null(error)) max_elbo(perm) else NA_real_
    if (is_slide) for (coding in levels) row[[paste0("n_", coding, "_slide_perm")]] <- if (is.null(error))
      sum(slide_leads(perm, map, TRUE, delta_tolerance)$lead_coding == coding) else NA_integer_
    if (!is.null(error)) perm_errors <- c(perm_errors, paste(key, error, sep = ": "))
  }
  row$dif_elbo_slide_vs_add_perm <- row$elbo_slide_perm - row$elbo_susie_perm
  info <- slide_tss_metadata(x, gene, annotations)
  x$tss <- info$tss_position
  x$strand <- info$strand
  x$chromosome <- info$chromosome
  # Supply only the two models being compared. Slide codings vary by component.
  inputs <- list(list(model = "SuSiE", model_key = "susie_add", fit = add, predictor_map = map,
                      saved_tss_summary = x[["susie_add_lead_snp_tss_distance"]]))
  cs_rows <- list(summarize_credible_sets(x, gene, tissue, file, row, model_inputs = inputs))
  for (i in seq_len(nrow(leads))) {
    component_map <- map
    component_map$coding <- slide_coding(fit$delta[leads$component[i], ], delta_tolerance)
    one_fit <- fit
    # Select one CS without changing its component row in alpha/delta.
    one_fit$sets$cs <- get_cs(fit)[i]
    one_fit$sets$cs_index <- fit$sets$cs_index[i]
    one_fit$sets$coverage <- fit$sets$coverage[i]
    if (!is.null(fit$sets$purity)) one_fit$sets$purity <- fit$sets$purity[i, , drop = FALSE]
    s <- summarize_credible_sets(x, gene, tissue, file, row, model_inputs = list(list(
      model = "SuSiE-slide", model_key = "fit_slide", fit = one_fit,
      predictor_map = component_map, saved_tss_summary = NULL)))
    s$cs_number <- i
    s$cs_id <- paste(gene, tissue, "fit_slide", i, sep = "|")
    s$lead_delta <- leads$lead_delta[i]
    s$lead_delta_forced <- leads$lead_delta_forced[i]
    cs_rows[[length(cs_rows) + 1L]] <- s
  }
  cs <- slide_bind_rows(cs_rows)
  if (nrow(cs)) {
    if (!"lead_delta" %in% names(cs)) cs$lead_delta <- NA_real_
    if (!"lead_delta_forced" %in% names(cs)) cs$lead_delta_forced <- NA
    cs$lead_delta[cs$model_key == "susie_add"] <- 0
    cs$tss_source <- info$tss_source
  }
  row$tss_position <- info$tss_position
  row$tss_source <- info$tss_source
  list(row = row, cs = cs, permutation_errors = perm_errors)
}

generate_summary_results_slide <- function(project_dir,
    results_dir = file.path(project_dir, "results"), output_dir = file.path(project_dir, "results_slide/summary"),
    gtf_file = "/project2/mstephens/gtex/Homo_sapiens.GRCh38.103.chr.reformatted.collapse_only.gene.gtf.gz",
    delta_tolerance = 1e-6) {
  slide_coding(0, delta_tolerance)
  files <- list.files(results_dir, pattern = "\\.rds$", full.names = TRUE)
  if (!length(files)) stop("No gene RDS files found in ", results_dir)
  annotations <- slide_gene_annotations(project_dir, gtf_file)
  rows <- cs_rows <- errors <- list()
  record <- function(file, gene, tissue, stage, message) {
    errors[[length(errors) + 1L]] <<- data.frame(result_file = basename(file), gene = gene,
      tissue = tissue, stage = stage, error = message)
  }
  for (i in seq_along(files)) {
    file <- files[i]
    gene <- tools::file_path_sans_ext(basename(file))
    out <- tryCatch(readRDS(file), error = identity)
    if (inherits(out, "error") || !is.list(out) || !is.null(out[["error"]])) {
      record(file, gene, NA_character_, "gene", if (inherits(out, "error")) conditionMessage(out) else
        if (is.list(out)) as.character(out$error) else "Not a tissue list")
      next
    }
    saved_errors <- attr(out, "tissue_errors", exact = TRUE)
    for (tissue in names(saved_errors)) record(file, gene, tissue, "workhorse", saved_errors[[tissue]])
    if (!length(out) || is.null(names(out)) || anyNA(names(out)) || any(!nzchar(names(out))) || anyDuplicated(names(out))) {
      record(file, gene, NA_character_, "gene", "No tissues or invalid tissue names")
      next
    }
    for (tissue in names(out)) {
      result <- tryCatch(slide_summarize_tissue(out[[tissue]], gene, tissue, file,
                                                annotations, delta_tolerance), error = identity)
      if (inherits(result, "error")) {
        record(file, gene, tissue, "summary", conditionMessage(result))
        next
      }
      rows[[length(rows) + 1L]] <- result$row
      cs_rows[[length(cs_rows) + 1L]] <- result$cs
      for (err in result$permutation_errors) record(file, gene, tissue, "permutation", err)
    }
    if (i %% 100L == 0L || i == length(files)) message("Processed ", i, "/", length(files), " gene files")
  }
  res_summary <- slide_bind_rows(rows)
  res_cs_summary <- slide_bind_rows(cs_rows, data.frame(gene = character(), tissue = character(),
    model = character(), model_key = character(), cs_number = integer(), component = integer(),
    lead_snp = character(), lead_coding = character(), lead_delta = numeric(), lead_delta_forced = logical(),
    lead_pip = numeric(), cs_size_snps = integer(), distance_to_tss_kb = numeric()))
  res_errors <- slide_bind_rows(errors, data.frame(result_file = character(), gene = character(),
    tissue = character(), stage = character(), error = character()))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  save(res_summary, file = file.path(output_dir, "res_summary.RData"))
  save(res_cs_summary, file = file.path(output_dir, "res_cs_summary.RData"))
  write.csv(res_summary, file.path(output_dir, "res_summary.csv"), row.names = FALSE)
  write.csv(res_cs_summary, file.path(output_dir, "res_cs_summary.csv"), row.names = FALSE)
  write.csv(res_errors, file.path(output_dir, "res_errors.csv"), row.names = FALSE)
  saveRDS(list(schema_version = 1L, results_dir = normalizePath(results_dir, winslash = "/"),
    generated_at = Sys.time(), n_files = length(files), delta_tolerance = delta_tolerance,
    gtf_file = gtf_file, models = c("susie_add", "fit_slide")), file.path(output_dir, "slide_summary_metadata.rds"))
  message("Slide summary: ", nrow(res_summary), " compared tissues, ", nrow(res_cs_summary),
          " CSs, ", nrow(res_errors), " audit records. Saved to ", output_dir)
  if (!nrow(res_summary)) stop("No valid SuSiE-slide comparisons; see res_errors.csv.")
  invisible(list(res_summary = res_summary, res_cs_summary = res_cs_summary, res_errors = res_errors))
}

slide_load_summary <- function(summary_dir) {
  e <- new.env(parent = emptyenv())
  load(file.path(summary_dir, "res_summary.RData"), envir = e)
  load(file.path(summary_dir, "res_cs_summary.RData"), envir = e)
  if (!is.data.frame(e$res_summary) || !nrow(e$res_summary) || !"ncs_slide" %in% names(e$res_summary))
    stop("Run generate_summary_results_slide.R first.")
  list(res = e$res_summary, cs = e$res_cs_summary)
}

slide_primary <- function(res, association_threshold, minimum_mean_reads) {
  res[which(is.finite(res$min_pv) & res$min_pv < association_threshold &
              is.finite(res$mean_count) & res$mean_count >= minimum_mean_reads), , drop = FALSE]
}

slide_region_key <- function(gene, tissue) paste0(nchar(gene), ":", gene, nchar(tissue), ":", tissue)
