# Opt-in tissue jobs for an already established, frozen gene/tissue cohort.
# All units use the existing controller and worker slot accounting.
sps_config_file <- function(ctx) file.path(ctx$project, "results_slide_prior_em/tissue_split_genes.rds")
sps_dir <- function(ctx, gene) file.path(ctx$dir, "tissue_results", gene)
sps_genes <- function(ctx) {
  f <- sps_config_file(ctx)
  if (!file.exists(f)) return(character())
  x <- readRDS(f)
  if (!identical(x$schema, 1L) || !is.character(x$genes) || anyNA(x$genes) ||
      anyDuplicated(x$genes) || any(!x$genes %in% ctx$manifest$gene) ||
      any(x$genes %in% ctx$scope$excluded)) stop("Invalid tissue-splitting selection.")
  if (ctx$iteration == 0L) stop("Tissue splitting requires a validated preceding iteration.")
  x$genes
}
sps_enable <- function(ctx, genes) {
  if (ctx$iteration == 0L) stop("Establish the initial gene/tissue cohort before enabling splitting.")
  genes <- unique(as.character(genes))
  if (!length(genes) || anyNA(genes) || any(!genes %in% ctx$manifest$gene)) stop("Unknown gene selection.")
  if (any(genes %in% ctx$scope$excluded))
    stop("Selected gene is excluded by the frozen chromosome scope. Tissue splitting cannot repair PLINK inputs or expand this run's cohort.")
  for (g in genes) {
    previous <- spp_previous_cohort(ctx, g)
    if (!is.null(previous$error) || !length(previous)) stop("No previously fitted tissues for ", g)
  }
  selected <- sort(unique(c(sps_genes(ctx), genes)))
  em_atomic_write(list(schema = 1L, genes = selected,
    updated_at = format(Sys.time(), tz = "UTC", usetz = TRUE)), sps_config_file(ctx))
  message("Tissue splitting enabled for: ", paste(genes, collapse = ", "),
          ". Existing jobs finish normally; the controller adopts splitting when the gene is idle. No new controller is needed.")
  invisible(selected)
}
sps_key <- function(ctx, gene) {
  sources <- file.path(ctx$project, "script/scan_tissue_attempt",
    c("workhorse_slide_prior_em.R", "workhorse_utils.R", "slide_prior_tissue_parallel.R",
      "slide_prior_em_utils.R", "get_gene_annotations.R"))
  hashes <- tools::md5sum(sources); names(hashes) <- basename(sources)
  if (anyNA(hashes)) stop("Missing tissue-worker source.")
  list(schema = 1L, key = spp_key(ctx), gene = gene, source_md5 = hashes,
       previous_stamp = spp_stamp(spp_previous_file(ctx, gene)))
}
sps_manifest <- function(ctx, gene, key = sps_key(ctx, gene)) {
  f <- file.path(sps_dir(ctx, gene), "manifest.rds")
  if (!file.exists(f)) return(NULL)
  m <- readRDS(f)
  if (!identical(m$key, key) || !is.character(m$tissues) || !length(m$tissues) ||
      anyNA(m$tissues) || anyDuplicated(m$tissues) || any(!m$tissues %in% ctx$priors$tissue) ||
      !identical(m$units, sprintf("tissue_%03d", seq_along(m$tissues))) ||
      length(m$input_md5) != length(m$tissues)) stop("Stale or invalid tissue manifest: ", gene)
  paths <- file.path(sps_dir(ctx, gene), "inputs", paste0(m$units, ".rds"))
  if (!identical(m$input_stamps, spp_stamp(paths))) stop("Prepared tissue inputs changed: ", gene)
  m
}
sps_done <- function(ctx, gene, unit, m) {
  f <- file.path(sps_dir(ctx, gene), "completed", paste0(unit, ".rds"))
  if (!file.exists(f)) return(FALSE)
  x <- readRDS(f); i <- match(unit, m$units)
  path <- file.path(sps_dir(ctx, gene), "results", paste0(unit, ".rds"))
  if (is.na(i) || !identical(x$key, m$key) || !identical(x$tissue, m$tissues[i]) ||
      !identical(x$input_md5, m$input_md5[i]) || !identical(x$result_stamp, spp_stamp(path)) ||
      !is.character(x$result_md5) || !grepl("^[a-f0-9]{32}$", x$result_md5))
    stop("Completed tissue result changed: ", gene, "/", unit)
  TRUE
}
sps_attempt_path <- function(ctx, gene, unit)
  file.path(sps_dir(ctx, gene), "attempts", paste0(unit, ".rds"))
sps_attempt <- function(ctx, gene, unit, key) {
  f <- sps_attempt_path(ctx, gene, unit)
  if (!file.exists(f)) return(0L)
  a <- readRDS(f)
  if (!identical(a$key, key)) stop("Tissue attempt inputs changed: ", gene, "/", unit)
  sre_int(a$attempt, "tissue attempt")
}
sps_plan <- function(ctx, genes, live, capacity) {
  work <- data.frame(task = integer(), gene = character(), action = character(), unit = character())
  for (g in genes) {
    active <- live$units$unit[live$units$gene == g]
    if (any(!grepl("^tissue_", active))) next
    key <- sps_key(ctx, g); m <- sps_manifest(ctx, g, key)
    if (is.null(m)) {
      if (length(active)) stop("Active tissue worker without prepared inputs: ", g)
      units <- "prepare"
    } else {
      done <- vapply(m$units, function(u) sps_done(ctx, g, u, m), logical(1))
      units <- setdiff(m$units[!done], active)
      if (all(done) && !length(active)) units <- "assemble"
    }
    units <- head(units, capacity - nrow(work))
    for (u in units) {
      if (sps_attempt(ctx, g, u, key) >= 3L)
        stop("Tissue-unit retry budget exhausted: ", g, "/", u,
             ". Completed tissues are retained; no prior update is allowed.")
      work <- rbind(work, data.frame(task = nrow(work) + 1L, gene = g, action = "fit", unit = u))
    }
    if (nrow(work) >= capacity) break
  }
  work
}
sps_prepare <- function(ctx, gene, key, run_gene = run_slide_prior_gene) {
  if (!is.null(sps_manifest(ctx, gene, key))) return(invisible(TRUE))
  if (sre_ready(spp_inspect(ctx, gene)$status)) return(invisible(TRUE))
  previous <- sre_previous(ctx, gene)
  if (is.null(previous) || !is.null(previous$error) || !length(previous)) stop("No fitted tissue cohort.")
  tissues <- names(previous); units <- sprintf("tissue_%03d", seq_along(tissues))
  seen <- character(); root <- sps_dir(ctx, gene)
  sink <- function(tissue, input) {
    i <- match(tissue, tissues)
    if (is.na(i) || tissue %in% seen) stop("Unexpected or duplicate prepared tissue.")
    if (!identical(input$metadata$sample_ids, previous[[tissue]]$sample_ids) ||
        !identical(colnames(input$args$X), colnames(previous[[tissue]]$fit_slide_prior$alpha)) ||
        !isTRUE(all.equal(unname(input$args$delta_prior), unname(spe_prior(ctx$priors, tissue)))))
      stop("Prepared tissue samples, SNPs or prior changed.")
    em_atomic_write(list(key = key, tissue = tissue, input = input),
                    file.path(root, "inputs", paste0(units[i], ".rds")))
    seen <<- c(seen, tissue)
  }
  settings <- ctx$settings; settings$max_iter <- spe_fit_budget(settings$max_iter)
  args <- c(list(target_gene = gene, tissue_priors = ctx$priors, project_dir = ctx$project,
    previous_result = previous, prepare_sink = sink,
    temp_dir = file.path(ctx$dir, "temp_plink", "tissue_prepare", gene)),
    settings[setdiff(names(settings), c("schema", "package_version", "cohort_source"))])
  out <- do.call(run_gene, args)
  if (!identical(seen, tissues) || !identical(names(out), tissues) ||
      !identical(attr(out, "tissue_errors", exact = TRUE), attr(previous, "tissue_errors", exact = TRUE)))
    stop("Frozen tissue cohort changed during preparation.")
  if (!identical(key, sps_key(ctx, gene))) stop("Inputs changed during tissue preparation.")
  paths <- file.path(root, "inputs", paste0(units, ".rds"))
  m <- list(key = key, tissues = tissues, units = units,
    tissue_errors = attr(previous, "tissue_errors", exact = TRUE), fit_max_iter = settings$max_iter,
    input_md5 = unname(tools::md5sum(paths)), input_stamps = spp_stamp(paths))
  em_atomic_write(data.frame(unit = units, tissue = tissues), file.path(root, "tissues.csv"), csv = TRUE)
  em_atomic_write(m, file.path(root, "manifest.rds"))
  invisible(TRUE)
}
sps_fit_tissue <- function(ctx, gene, unit, m, fit_tissue = spe_fit_prepared_slide_tissue) {
  if (sps_done(ctx, gene, unit, m)) return(invisible(TRUE))
  i <- match(unit, m$units)
  if (is.na(i)) stop("Unknown tissue unit.")
  path <- file.path(sps_dir(ctx, gene), "inputs", paste0(unit, ".rds"))
  if (!identical(unname(tools::md5sum(path)), m$input_md5[i])) stop("Prepared tissue checksum changed.")
  prepared <- readRDS(path)
  if (!identical(prepared$key, m$key) || !identical(prepared$tissue, m$tissues[i])) stop("Wrong tissue input.")
  # Fitting is deterministic. This seed also makes CS purity subsampling
  # reproducible across retries of this tissue (it is not part of EM counts).
  assign(".Random.seed", prepared$input$rng_state, envir = .GlobalEnv)
  result <- fit_tissue(prepared$input, ctx$project)
  out <- setNames(list(result), m$tissues[i]); attr(out, "em_iteration") <- ctx$iteration
  sre_validate(out, ctx, require_converged = TRUE)
  if (!identical(result$sample_ids, prepared$input$metadata$sample_ids) ||
      !identical(colnames(result$fit_slide_prior$alpha), colnames(prepared$input$args$X)))
    stop("Fitted tissue samples or SNPs changed.")
  if (!identical(m$key, sps_key(ctx, gene))) stop("Inputs changed during tissue fitting.")
  path <- file.path(sps_dir(ctx, gene), "results", paste0(unit, ".rds"))
  em_atomic_write(out, path)
  em_atomic_write(list(key = m$key, tissue = m$tissues[i], input_md5 = m$input_md5[i],
    result_md5 = unname(tools::md5sum(path)), result_stamp = spp_stamp(path)),
    file.path(sps_dir(ctx, gene), "completed", paste0(unit, ".rds")))
  invisible(TRUE)
}
sps_assemble <- function(ctx, gene, m) {
  if (!all(vapply(m$units, function(u) sps_done(ctx, gene, u, m), logical(1))))
    stop("Unfinished tissues; gene cannot be assembled.")
  out <- setNames(vector("list", length(m$tissues)), m$tissues)
  for (i in seq_along(m$units)) {
    path <- file.path(sps_dir(ctx, gene), "results", paste0(m$units[i], ".rds"))
    receipt <- readRDS(file.path(sps_dir(ctx, gene), "completed", paste0(m$units[i], ".rds")))
    if (!identical(unname(tools::md5sum(path)), receipt$result_md5)) stop("Tissue result checksum changed.")
    part <- readRDS(path)
    if (!identical(names(part), m$tissues[i])) stop("Wrong tissue result in assembly.")
    sre_validate(part, ctx, require_converged = TRUE)
    out[[i]] <- part[[1]]
  }
  attr(out, "tissue_errors") <- m$tissue_errors
  attr(out, "em_iteration") <- ctx$iteration
  attr(out, "fit_max_iter") <- m$fit_max_iter
  sre_validate(out, ctx, spp_previous_cohort(ctx, gene))
  if (!identical(m$key, sps_key(ctx, gene))) stop("Inputs changed during assembly.")
  em_atomic_write(out, sre_result(ctx, gene))
  row <- spp_inspect(ctx, gene, force = TRUE)
  if (!sre_ready(row$status)) stop("Assembled gene failed validation: ", row$detail)
  invisible(TRUE)
}
sps_worker <- function(ctx, gene, unit, run_gene = run_slide_prior_gene,
                       fit_tissue = spe_fit_prepared_slide_tissue) {
  if (ctx$iteration == 0L || !gene %in% sps_genes(ctx)) stop("Gene was not enabled for tissue splitting.")
  key <- sps_key(ctx, gene); m <- sps_manifest(ctx, gene, key)
  if (unit == "prepare" && !is.null(m)) return(invisible(TRUE))
  if (grepl("^tissue_", unit) && !is.null(m) && sps_done(ctx, gene, unit, m)) return(invisible(TRUE))
  if (unit == "assemble" && sre_ready(spp_inspect(ctx, gene)$status)) return(invisible(TRUE))
  attempt <- sps_attempt(ctx, gene, unit, key) + 1L
  if (attempt > 3L) stop("Tissue-unit retry budget exhausted: ", gene, "/", unit)
  receipt <- list(key = key, unit = unit, attempt = attempt, slurm_job_id = Sys.getenv("SLURM_JOB_ID"),
    started_at = format(Sys.time(), tz = "UTC", usetz = TRUE), error = "Interrupted before receipt")
  path <- sps_attempt_path(ctx, gene, unit); em_atomic_write(receipt, path)
  error <- tryCatch({
    if (unit == "prepare") sps_prepare(ctx, gene, key, run_gene)
    else {
      if (is.null(m)) stop("No prepared tissue manifest.")
      if (unit == "assemble") sps_assemble(ctx, gene, m)
      else if (grepl("^tissue_[0-9]{3}$", unit)) sps_fit_tissue(ctx, gene, unit, m, fit_tissue)
      else stop("Unknown tissue unit.")
    }
    NULL
  }, error = conditionMessage)
  receipt$error <- error; receipt$finished_at <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  em_atomic_write(receipt, path)
  if (!is.null(error)) stop(error)
  invisible(TRUE)
}
sps_report <- function(ctx) {
  files <- list.files(file.path(ctx$root, "attempts"), "\\.rds$", full.names = TRUE)
  rows <- lapply(files, function(f) {
    a <- readRDS(f)
    start <- as.POSIXct(a$started_at, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    finished <- if (is.null(a$finished_at)) as.POSIXct(NA) else
      as.POSIXct(a$finished_at, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    data.frame(gene = a$gene, attempt = a$attempt,
      runtime_hours = round(as.numeric(difftime(finished, start, units = "hours")), 2),
      state = if (is.null(a$error)) "validated" else a$error)
  })
  x <- if (length(rows)) do.call(rbind, rows) else data.frame()
  if (nrow(x)) print(head(x[order(x$runtime_hours, decreasing = TRUE, na.last = TRUE), ], 30L), row.names = FALSE)
  message("Unfinished receipts do not prove a timeout; inspect the matching Slurm log. Enabled genes: ",
          paste(if (ctx$iteration > 0L) sps_genes(ctx) else character(), collapse = ", "))
  invisible(x)
}
