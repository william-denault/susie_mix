# Parallel execution of the existing tissue-wide slider-prior EM.
# Large posteriors are read by workers only. The coordinator uses atomic,
# per-gene summaries and a compact index. Completed outputs are immutable under
# the gene lock; size/mtime changes invalidate summaries, whose creator also
# records the result checksum. Frozen priors/settings/scope are checked always.
spp_root <- function(ctx) file.path(ctx$root, "parallel")
spp_summary_file <- function(ctx, gene) file.path(spp_root(ctx), "summaries", paste0(gene, ".rds"))
spp_stamp <- function(files) {
  x <- file.info(files)
  data.frame(size = unname(x$size), mtime = unname(as.numeric(x$mtime)), row.names = NULL)
}
spp_review_hash <- function(ctx) {
  f <- file.path(ctx$root, "reviewed_exclusions.csv")
  if (file.exists(f)) unname(tools::md5sum(f)) else NULL
}
spp_key <- function(ctx) list(schema = 1L, frozen = ctx$frozen, scope = ctx$scope$md5,
                             review = spp_review_hash(ctx))
spp_previous_file <- function(ctx, gene) {
  if (ctx$iteration == 0L) return(character())
  file.path(ctx$project, "results_slide_prior_em", sprintf("iteration_%03d", ctx$iteration - 1L),
            "results", paste0(gene, ".rds"))
}
spp_cohort <- function(out) {
  if (!is.null(out$error)) return(list(error = out$error))
  x <- setNames(rep(list(list()), length(out)), names(out))
  attr(x, "tissue_errors") <- attr(out, "tissue_errors", exact = TRUE)
  x
}
spp_init <- function(ctx) {
  if (!file.exists(file.path(ctx$root, "state.rds"))) {
    # A new iteration has no fits. Never call the old serial audit on existing
    # results: those must be adopted by the parallel workers below.
    plan <- data.frame(task = ceiling(seq_len(nrow(ctx$manifest))/10),
                       gene = ctx$manifest$gene, original_chunk = ctx$manifest$chunk)
    em_atomic_write(plan, file.path(ctx$root, "manifest.csv"), csv = TRUE)
    state <- list(schema = 1L, frozen = ctx$frozen, genes_per_task = 10L, max_attempts = 3L,
                  scope_md5 = ctx$scope$md5,
                  manifest_md5 = unname(tools::md5sum(file.path(ctx$root, "manifest.csv"))))
    em_atomic_write(state, file.path(ctx$root, "state.rds"))
  }
  sre_load(ctx)
  dir.create(spp_root(ctx), recursive = TRUE, showWarnings = FALSE)
  marker <- file.path(ctx$project, "results_slide_prior_em/parallel_adopted.rds")
  if (!file.exists(marker))
    em_atomic_write(list(schema = 1L, adopted_at = format(Sys.time(), tz = "UTC", usetz = TRUE)), marker)
  invisible(ctx)
}
spp_valid_summary <- function(x, ctx, gene, stamp, previous_stamp, key = spp_key(ctx)) {
  if (!is.list(x) || !identical(x$key, key) || !identical(x$gene, gene) ||
      !identical(x$result_stamp, stamp) || !identical(x$previous_stamp, previous_stamp) ||
      !is.data.frame(x$row) || nrow(x$row) != 1L || !identical(x$row$gene, gene) ||
      !x$row$status %in% c("valid", "qc_excluded", "reviewed_exclusion", "scope_excluded", "retry"))
    return(FALSE)
  if (sre_ready(x$row$status) && x$row$status != "scope_excluded") {
    if (!is.list(x$summary) || !is.matrix(x$summary$counts) ||
        !identical(dimnames(x$summary$counts), list(ctx$priors$tissue, spe_columns())) ||
        any(!is.finite(x$summary$counts)) || any(x$summary$counts < 0) ||
        !is.data.frame(x$summary$audit) || !nrow(x$summary$audit) ||
        any(x$summary$audit$gene != gene) ||
        sum(x$summary$audit$issue == "fit") != x$row$n_fits ||
        !grepl("^[a-f0-9]{32}$", x$row$result_md5)) return(FALSE)
    for (t in ctx$priors$tissue) {
      a <- x$summary$audit
      if (abs(sum(x$summary$counts[t, ]) - sum(a$free_mass[a$tissue == t & a$issue == "fit"])) > 1e-7)
        return(FALSE)
    }
  }
  TRUE
}
spp_previous_cohort <- function(ctx, gene) {
  if (ctx$iteration == 0L) return(NULL)
  d <- dirname(dirname(spp_previous_file(ctx, gene)))
  prev_ctx <- sre_context(ctx$project, d)
  x <- sre_read(spp_summary_file(prev_ctx, gene))
  if (!spp_valid_summary(x, prev_ctx, gene, spp_stamp(sre_result(prev_ctx, gene)),
                        spp_stamp(spp_previous_file(prev_ctx, gene))) || !sre_ready(x$row$status))
    stop("Previous iteration summary is missing or stale: ", gene)
  x$cohort
}
spp_inspect <- function(ctx, gene, force = FALSE) {
  stamp <- spp_stamp(sre_result(ctx, gene)); prev_stamp <- spp_stamp(spp_previous_file(ctx, gene))
  key <- spp_key(ctx)
  cached <- sre_read(spp_summary_file(ctx, gene))
  if (!force && spp_valid_summary(cached, ctx, gene, stamp, prev_stamp, key)) return(cached$row)
  summary <- cohort <- NULL
  sink <- function(out, row) {
    summary <<- spe_summarize_gene(out, gene, ctx$priors)
    cohort <<- spp_cohort(out)
  }
  row <- sre_inspect(ctx, gene, summary_sink = sink,
                     previous = spp_previous_cohort(ctx, gene))
  if (!identical(stamp, spp_stamp(sre_result(ctx, gene))) ||
      !identical(prev_stamp, spp_stamp(spp_previous_file(ctx, gene))))
    stop("Result changed during validation: ", gene)
  x <- list(key = key, gene = gene, result_stamp = stamp, previous_stamp = prev_stamp,
            row = row, summary = summary, cohort = cohort,
            validated_at = format(Sys.time(), tz = "UTC", usetz = TRUE))
  em_atomic_write(x, spp_summary_file(ctx, gene))
  row
}

# Load only changed small summaries. Never checksum/decompress full gene fits
# here. A worker's atomic rename changes the stamp and forces recollection.
spp_index <- function(ctx, active_genes = character()) {
  genes <- ctx$manifest$gene; key <- spp_key(ctx)
  files <- vapply(genes, function(g) spp_summary_file(ctx, g), character(1))
  summary_stamps <- spp_stamp(files)
  result_stamps <- spp_stamp(file.path(ctx$dir, "results", paste0(genes, ".rds")))
  previous_stamps <- if (ctx$iteration > 0L)
    spp_stamp(vapply(genes, function(g) spp_previous_file(ctx, g), character(1))) else NULL
  f <- file.path(spp_root(ctx), "index.rds")
  index <- sre_read(f)
  if (is.null(index) || !identical(index$key, key) || !identical(names(index$entries), genes))
    index <- list(key = key, entries = setNames(vector("list", length(genes)), genes),
                  stamps = data.frame(size = rep(NA_real_, length(genes)), mtime = NA_real_))
  for (i in seq_along(genes)) {
    g <- genes[i]
    if (g %in% active_genes) next
    rs <- result_stamps[i, , drop = FALSE]; rownames(rs) <- NULL
    ps <- if (is.null(previous_stamps)) spp_stamp(character()) else previous_stamps[i, , drop = FALSE]
    rownames(ps) <- NULL
    if (g %in% ctx$scope$excluded) {
      index$entries[[i]] <- list(key = key, gene = g, result_stamp = rs, previous_stamp = ps,
        row = data.frame(gene = g, status = "scope_excluded", detail = "Frozen chromosome exclusion",
                         result_md5 = "", n_fits = 0L))
      next
    }
    changed <- !identical(unname(as.numeric(index$stamps[i, ])), unname(as.numeric(summary_stamps[i, ])))
    x <- if (changed || is.null(index$entries[[i]])) sre_read(files[i]) else index$entries[[i]]
    valid <- spp_valid_summary(x, ctx, g, rs, ps, key)
    if (!valid && !changed) {
      x <- sre_read(files[i]); valid <- spp_valid_summary(x, ctx, g, rs, ps, key)
    }
    index$entries[i] <- list(if (valid) x else NULL)
  }
  index$stamps <- summary_stamps
  em_atomic_write(index, f)
  index
}
spp_batch <- function(ctx, batch) {
  batch <- normalizePath(batch, winslash = "/", mustWork = TRUE)
  if (!identical(dirname(batch), paste0(spp_root(ctx), "/batches"))) stop("Foreign parallel batch.")
  settings <- readRDS(file.path(batch, "settings.rds"))
  f <- file.path(batch, "manifest.csv")
  if (!identical(settings$key, spp_key(ctx)) || !identical(settings$manifest_md5, unname(tools::md5sum(f))))
    stop("Parallel batch inputs changed.")
  work <- read.csv(f, stringsAsFactors = FALSE)
  if (!"unit" %in% names(work)) work$unit <- "gene"
  if (!nrow(work) || anyNA(work) ||
      anyDuplicated(paste(work$gene, work$unit, sep = ":")) || any(!work$gene %in% ctx$manifest$gene) ||
      any(!work$action %in% c("validate", "fit")) || length(unique(work$action)) != 1L ||
      any(!grepl("^(gene|prepare|assemble|tissue_[0-9]{3})$", work$unit)) ||
      any(work$action == "validate" & work$unit != "gene") ||
      !identical(sort(unique(as.integer(work$task))), seq_len(max(work$task))) || max(work$task) > 298L ||
      any(table(work$task) > if (work$action[1] == "fit") 1L else 25L)) stop("Invalid parallel assignment.")
  spp_check_units(work)
  work
}
spp_check_units <- function(work) {
  for (g in unique(work$gene)) {
    u <- work$unit[work$gene == g]
    if (anyDuplicated(u) || (length(u) > 1L && any(!grepl("^tissue_", u))))
      stop("Overlapping active gene assignments: ", g)
  }
  invisible(TRUE)
}
spp_queue <- function(ctx, queue_file) {
  lines <- readLines(queue_file, warn = FALSE); lines <- lines[nzchar(lines)]
  q <- data.frame(id = character(), array = character(), name = character())
  if (length(lines)) {
    parts <- strsplit(lines, "|", fixed = TRUE)
    if (any(lengths(parts) != 3L)) stop("Invalid queue snapshot.")
    q <- as.data.frame(do.call(rbind, parts), stringsAsFactors = FALSE); names(q) <- c("id", "array", "name")
  }
  if (anyDuplicated(q$id)) stop("Duplicate queued job.")
  f <- list.files(file.path(spp_root(ctx), "batches"), "^array_job_id\\.txt$", recursive = TRUE, full.names = TRUE)
  ids <- vapply(f, function(p) {
    id <- readLines(p, warn = FALSE)
    if (length(id) != 1L || !grepl("^[0-9]+$", id)) stop("Invalid array registration.")
    id
  }, character(1))
  if (anyDuplicated(ids)) stop("Duplicate array registration.")
  workers <- q$name %in% c("slide_parallel_validate", "slide_parallel_fit")
  if (any(!q$array[workers] %in% ids)) stop("Unregistered parallel worker; inspect submissions.")
  active <- data.frame(gene = character(), unit = character()); tasks <- 0L
  for (i in seq_along(f)) {
    batch <- dirname(f[i]); rows <- which(q$array == ids[i])
    # Finished batches are inspected once, then a tiny marker avoids repeated
    # reads of thousands of historical manifests during five-minute polling.
    if (!length(rows) && file.exists(file.path(batch, "closed.rds"))) next
    w <- spp_batch(ctx, batch)
    expected <- paste0(ids[i], "_", unique(w$task))
    if (length(rows)) {
      if (any(!q$id[rows] %in% expected) || any(!workers[rows])) stop("Unexpected registered worker.")
      ts <- unique(w$task)[match(q$id[rows], expected)]
      active <- rbind(active, w[w$task %in% ts, c("gene", "unit"), drop = FALSE])
      tasks <- tasks + length(rows)
    } else {
      receipts <- file.path(batch, paste0("started_", unique(w$task), ".rds"))
      if (any(!file.exists(receipts))) stop("A completed worker never entered R: ", batch,
                                           ". Inspect its Slurm log; no automatic resubmission.")
      em_atomic_write(TRUE, file.path(batch, "closed.rds"))
    }
  }
  spp_check_units(active)
  list(genes = unique(active$gene), units = active, tasks = tasks)
}
spp_validation_attempts <- function(ctx, gene) {
  f <- file.path(spp_root(ctx), "validation_attempts", paste0(gene, ".rds"))
  x <- sre_read(f)
  if (is.null(x)) return(0L)
  if (!identical(x$key, spp_key(ctx)) || !identical(x$stamp, spp_stamp(sre_result(ctx, gene)))) return(0L)
  sre_int(x$attempt, "validation attempt")
}
spp_task <- function(ctx, batch, task) {
  w <- spp_batch(ctx, batch); w <- w[w$task == sre_int(task, "task", 1L, 298L), , drop = FALSE]
  if (!nrow(w)) stop("Unknown parallel task.")
  em_atomic_write(list(time = format(Sys.time(), tz = "UTC", usetz = TRUE),
                       job = Sys.getenv("SLURM_JOB_ID")), file.path(batch, paste0("started_", task, ".rds")))
  w
}
spp_worker_gene <- function(ctx, gene, action, run_gene = NULL) {
  if (action == "validate") {
    # Validation interruptions do not consume the scientific fitting retry budget.
    n <- spp_validation_attempts(ctx, gene) + 1L
    limit <- sre_int(Sys.getenv("SUSIE_SLIDE_VALIDATION_MAX_ATTEMPTS", "3"), "validation retry limit", 1L, 20L)
    if (n > limit) stop("Validation retry budget exhausted: ", gene)
    em_atomic_write(list(key = spp_key(ctx), stamp = spp_stamp(sre_result(ctx, gene)), attempt = n),
                    file.path(spp_root(ctx), "validation_attempts", paste0(gene, ".rds")))
    spp_inspect(ctx, gene, force = TRUE)
  } else if (action == "fit") {
    if (is.null(run_gene)) run_gene <- run_slide_prior_gene
    sre_run_gene(ctx, gene, run_gene = run_gene, inspect = spp_inspect)
  } else stop("Unknown worker action.")
}
spp_complete <- function(ctx, index) {
  if (any(vapply(index$entries, is.null, logical(1))) ||
      any(!vapply(index$entries, function(x) sre_ready(x$row$status), logical(1))))
    stop("Incomplete gene summaries; no prior update.")
  entries <- index$entries[!names(index$entries) %in% ctx$scope$excluded]
  entries <- entries[order(names(entries))] # Same addition order as the serial file reader.
  result_stamps <- do.call(rbind, lapply(index$entries, `[[`, "result_stamp")); rownames(result_stamps) <- NULL
  previous_stamps <- if (ctx$iteration > 0L) do.call(rbind, lapply(index$entries, `[[`, "previous_stamp")) else NULL
  if (!is.null(previous_stamps)) rownames(previous_stamps) <- NULL
  # Unexpected files are not silently included or used to redefine the cohort.
  sre_pool_files(ctx)
  pooled <- spe_pool_summaries(lapply(entries, `[[`, "summary"), ctx$priors)
  if (ctx$iteration > 0L) {
    d <- file.path(ctx$project, "results_slide_prior_em", sprintf("iteration_%03d", ctx$iteration - 1L))
    prev <- readRDS(file.path(d, "recovery/parallel/pooled.rds"))
    cohort <- function(a) {
      a <- a[a$issue == "fit", c("gene", "tissue"), drop = FALSE]
      a <- a[order(a$gene, a$tissue), , drop = FALSE]; rownames(a) <- NULL; a
    }
    if (!identical(cohort(prev$audit), cohort(pooled$audit))) stop("Frozen gene/tissue cohort changed.")
  }
  audit <- do.call(rbind, lapply(index$entries, `[[`, "row")); rownames(audit) <- NULL
  em_atomic_write(audit, file.path(ctx$root, "audit.csv"), csv = TRUE)
  em_atomic_write(pooled, file.path(spp_root(ctx), "pooled.rds"))
  em_atomic_write(pooled$counts, file.path(ctx$dir, "slider_counts.csv"), csv = TRUE)
  em_atomic_write(pooled$priors, file.path(ctx$dir, "proposed_priors.csv"), csv = TRUE)
  if (Sys.getenv("SUSIE_SLIDE_RECORD_ELBO", "0") == "1") spe_record_objective(ctx$project, ctx$dir, pooled)
  # Refuse to certify a different generation of files than the summaries just
  # pooled. This also catches an external change during coordinator aggregation.
  current_results <- spp_stamp(file.path(ctx$dir, "results", paste0(ctx$manifest$gene, ".rds")))
  current_summaries <- spp_stamp(vapply(ctx$manifest$gene, function(g) spp_summary_file(ctx, g), character(1)))
  current_previous <- if (ctx$iteration > 0L)
    spp_stamp(vapply(ctx$manifest$gene, function(g) spp_previous_file(ctx, g), character(1))) else NULL
  if (!identical(result_stamps, current_results) || !identical(index$stamps, current_summaries) ||
      !identical(previous_stamps, current_previous)) stop("Inputs changed while pooling; resume validation.")
  cert <- list(key = spp_key(ctx), genes = ctx$manifest$gene,
    result_stamps = current_results, summary_stamps = current_summaries, previous_stamps = current_previous,
    pooled_md5 = unname(tools::md5sum(file.path(spp_root(ctx), "pooled.rds"))),
    completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE))
  em_atomic_write(cert, file.path(spp_root(ctx), "COMPLETE.rds"))
  message("Validated iteration ", ctx$iteration, "; collected counts for ", sum(audit$n_fits), " fits.")
  invisible(pooled)
}
spp_require_complete <- function(project, iteration_dir) {
  ctx <- sre_context(project, iteration_dir)
  c <- sre_read(file.path(spp_root(ctx), "COMPLETE.rds"))
  previous_stamps <- if (ctx$iteration > 0L)
    spp_stamp(vapply(ctx$manifest$gene, function(g) spp_previous_file(ctx, g), character(1))) else NULL
  if (is.null(c) || !identical(c$key, spp_key(ctx)) || !identical(c$genes, ctx$manifest$gene) ||
      !identical(c$previous_stamps, previous_stamps) ||
      !identical(c$result_stamps, spp_stamp(file.path(ctx$dir, "results", paste0(ctx$manifest$gene, ".rds")))) ||
      !identical(c$summary_stamps, spp_stamp(vapply(ctx$manifest$gene, function(g) spp_summary_file(ctx, g), character(1)))) ||
      !identical(c$pooled_md5, unname(tools::md5sum(file.path(spp_root(ctx), "pooled.rds")))))
    stop("Parallel completion is missing or changed; resume parallel validation.")
  readRDS(file.path(spp_root(ctx), "pooled.rds"))
}
spp_plan <- function(ctx, capacity, queue_file) {
  capacity <- sre_int(capacity, "capacity", 0L, 298L)
  live <- spp_queue(ctx, queue_file)
  index <- spp_index(ctx, live$genes)
  genes <- ctx$manifest$gene
  status <- vapply(index$entries, function(x) if (is.null(x)) "unvalidated" else x$row$status, character(1))
  message("Iteration ", ctx$iteration, ": ", sum(sre_ready(status)), "/", length(genes),
          " resolved; ", live$tasks, " active tasks; ", capacity, " submission slots.")
  if (!live$tasks && all(sre_ready(status))) {
    spp_complete(ctx, index)
    return(list(status = "complete", iteration = ctx$iteration))
  }
  split_genes <- if (exists("sps_genes", mode = "function")) sps_genes(ctx) else character()
  if (capacity > 0L && length(split_genes)) {
    work <- sps_plan(ctx, split_genes[!sre_ready(status[match(split_genes, genes)])], live, capacity)
    if (nrow(work)) return(spp_write_batch(ctx, work))
  }
  pending <- which(!sre_ready(status) & !genes %in% live$genes)
  pending <- pending[!genes[pending] %in% split_genes]
  if (!length(pending) || capacity == 0L) return(list(status = "waiting", iteration = ctx$iteration))
  validate <- pending[status[pending] == "unvalidated" & file.exists(sre_result(ctx, genes[pending]))]
  action <- if (length(validate)) "validate" else "fit"
  choose <- if (action == "validate") validate else pending
  # Bound reads/retry checks to the number of genes this submission can own.
  choose <- head(choose, capacity * if (action == "validate") 25L else 1L)
  if (action == "validate") {
    limits <- sre_int(Sys.getenv("SUSIE_SLIDE_VALIDATION_MAX_ATTEMPTS", "3"), "validation retry limit", 1L, 20L)
    blocked <- genes[choose][vapply(genes[choose], function(g) spp_validation_attempts(ctx, g) >= limits, logical(1))]
  } else {
    counts <- vapply(genes[choose], function(g) {
      f <- sre_attempt_file(ctx, g); a <- sre_read(f)
      if (is.null(a)) { if (file.exists(f)) stop("Unreadable fitting receipt: ", g); return(0L) }
      if (!identical(a$frozen, ctx$frozen)) stop("Fitting receipt inputs changed.")
      sre_int(a$attempt, "attempt")
    }, integer(1))
    blocked <- genes[choose][counts >= sre_attempt_limits(ctx, sre_load(ctx)$state, genes[choose])]
    details <- vapply(index$entries[choose], function(x) if (is.null(x)) "" else x$row$detail, character(1))
    resource <- grepl("cannot allocate|out.of.memory|oom|system call failed", details, ignore.case = TRUE)
    if (sum(resource & counts >= 1L) >= 2L || any(resource & counts >= 2L))
      stop("Repeated resource failures; inspect workers before resuming.")
  }
  if (length(blocked)) {
    em_atomic_write(data.frame(gene = blocked, action = action), file.path(spp_root(ctx), "blocked.csv"), csv = TRUE)
    stop("Retry budget exhausted; inspect parallel/blocked.csv.")
  }
  work <- data.frame(task = if (action == "validate") ceiling(seq_along(choose)/25) else seq_along(choose),
                     gene = genes[choose], action = action)
  spp_write_batch(ctx, work)
}
spp_write_batch <- function(ctx, work) {
  dir.create(file.path(spp_root(ctx), "batches"), showWarnings = FALSE)
  batch <- tempfile("batch_", tmpdir = file.path(spp_root(ctx), "batches")); dir.create(batch)
  em_atomic_write(work, file.path(batch, "manifest.csv"), csv = TRUE)
  em_atomic_write(list(key = spp_key(ctx), manifest_md5 = unname(tools::md5sum(file.path(batch, "manifest.csv")))),
                  file.path(batch, "settings.rds"))
  list(status = "batch", iteration = ctx$iteration, batch = batch, tasks = max(work$task), action = work$action[1])
}
