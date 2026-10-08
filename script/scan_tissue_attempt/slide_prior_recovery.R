# Execution and recovery only: spe_pool(), spe_settings() and the fitter are unchanged.
sre_int <- function(x, name, low = 1L, high = .Machine$integer.max) {
  y <- suppressWarnings(as.numeric(x))
  if (length(y) != 1L || !is.finite(y) || y != floor(y) || y < low || y > high)
    stop(name, " must be an integer from ", low, " to ", high, ".")
  as.integer(y)
}

sre_same_priors <- function(a, b) {
  identical(names(a), names(b)) && nrow(a) == nrow(b) && all(vapply(names(a), function(k) {
    # CSV inference changes an all-NA column from logical to numeric once later
    # history rows contain numbers; that is not a change in the frozen values.
    if (all(is.na(a[[k]])) && all(is.na(b[[k]]))) return(TRUE)
    isTRUE(all.equal(a[[k]], b[[k]], check.attributes = FALSE))
  }, logical(1)))
}

sre_scope <- function(project, manifest) {
  file <- file.path(project, "results_slide_prior_em/chromosome_scope.rds")
  if (!file.exists(file)) return(list(md5 = NULL, excluded = character(), chromosomes = NULL))
  policy <- readRDS(file)
  map <- policy$genes
  if (!identical(policy$schema, 1L) || !identical(policy$exclude, c("X", "Y", "MT")) ||
      !is.data.frame(map) || !identical(map$gene, manifest$gene) ||
      !is.character(map$chromosome) || anyNA(map$chromosome) || any(!nzchar(map$chromosome)))
    stop("Invalid frozen chromosome scope.")
  list(md5 = unname(tools::md5sum(file)), excluded = map$gene[map$chromosome %in% policy$exclude],
       chromosomes = setNames(map$chromosome, map$gene))
}

sre_context <- function(project, iteration_dir, check_package = FALSE, allow_scope_adoption = FALSE) {
  project <- normalizePath(project, winslash = "/", mustWork = TRUE)
  iteration_dir <- normalizePath(iteration_dir, winslash = "/", mustWork = TRUE)
  if (!identical(dirname(iteration_dir), paste0(project, "/results_slide_prior_em")) ||
      !grepl("^iteration_[0-9]{3,}$", basename(iteration_dir))) stop("Invalid iteration path.")
  priors <- read.csv(file.path(iteration_dir, "priors.csv"), stringsAsFactors = FALSE)
  spe_validate_priors(priors)
  iteration <- unique(priors$iteration)
  if (length(iteration) != 1L || iteration < 0 || iteration != floor(iteration) ||
      !identical(basename(iteration_dir), sprintf("iteration_%03d", iteration)) ||
      !all(priors$source_iteration == iteration - 1L)) stop("Invalid iteration metadata.")
  history <- read.csv(file.path(project, "results_slide_prior_em/prior_history.csv"))
  if (!sre_same_priors(priors, history[history$iteration == iteration, , drop = FALSE]))
    stop("Frozen priors differ from history.")
  settings <- readRDS(file.path(iteration_dir, "settings.rds"))
  if (!identical(spe_settings(), settings[names(spe_settings())]))
    stop("Fitting settings changed since preparation.")
  if (check_package && !identical(spe_check_packages(), settings$package_version))
    stop("Package version changed since preparation.")
  manifest <- read.csv(file.path(iteration_dir, "manifest.csv"), stringsAsFactors = FALSE)
  if (!nrow(manifest) || !all(c("gene", "chunk") %in% names(manifest)) ||
      anyNA(manifest$gene) || anyDuplicated(manifest$gene) ||
      any(!grepl("^[A-Za-z0-9_.-]+$", manifest$gene)) ||
      any(manifest$gene %in% c(".", ".."))) stop("Invalid gene manifest.")
  frozen <- tools::md5sum(file.path(iteration_dir, c("priors.csv", "settings.rds", "manifest.csv", "grid.csv")))
  names(frozen) <- basename(names(frozen))
  if (anyNA(frozen)) stop("Missing frozen iteration file.")
  scope <- sre_scope(project, manifest)
  root <- file.path(iteration_dir, "recovery")
  if (file.exists(file.path(root, "state.rds"))) {
    state <- readRDS(file.path(root, "state.rds"))
    if (!identical(state$frozen, frozen)) stop("Frozen iteration files changed after recovery audit.")
    if (!identical(state$scope_md5, scope$md5) &&
        !(allow_scope_adoption && iteration == 0L && is.null(state$scope_md5)))
      stop("Chromosome scope changed after recovery audit.")
  }
  list(project = project, dir = iteration_dir, root = root, priors = priors,
       settings = settings, manifest = manifest, iteration = iteration, frozen = frozen, scope = scope)
}

sre_read <- function(file) {
  if (!file.exists(file)) return(NULL)
  suppressWarnings(tryCatch(readRDS(file), error = function(e) NULL))
}
sre_result <- function(ctx, gene) file.path(ctx$dir, "results", paste0(gene, ".rds"))
sre_attempt_file <- function(ctx, gene) file.path(ctx$root, "attempts", paste0(gene, ".rds"))
sre_ready <- function(status) status %in% c("valid", "qc_excluded", "reviewed_exclusion", "scope_excluded")

# Explicit scope selection before the first learned-prior update. Keep the
# original manifests/results; the policy governs scheduling AND pooling.
sre_exclude_sex_mt <- function(ctx, annotations = NULL,
    gtf_file = "/project2/mstephens/gtex/Homo_sapiens.GRCh38.103.chr.reformatted.collapse_only.gene.gtf.gz") {
  history <- read.csv(file.path(ctx$project, "results_slide_prior_em/prior_history.csv"))
  if (ctx$iteration != 0L || max(history$iteration) != 0L)
    stop("Select chromosome scope during initialization, before learned-prior updates.")
  loaded <- sre_load(ctx, allow_scope_adoption = TRUE)
  if (is.null(ctx$scope$md5)) {
    if (is.null(annotations)) {
      helper <- new.env(parent = globalenv())
      helper$fread <- data.table::fread
      sys.source(file.path(ctx$project, "script/scan_tissue_attempt/get_gene_annotations.R"), helper)
      annotations <- helper$get_gene_annotations(gtf_file)
      source_info <- list(path = gtf_file, md5 = unname(tools::md5sum(gtf_file)))
    } else source_info <- list(path = "supplied annotation table", md5 = NULL)
    if (!all(c("gene_name", "chromosome") %in% names(annotations))) stop("Invalid gene annotation table.")
    # match() chooses the first annotation, exactly as the existing worker does.
    chromosome <- as.character(annotations$chromosome[match(ctx$manifest$gene, annotations$gene_name)])
    chromosome <- toupper(sub("^chr", "", chromosome, ignore.case = TRUE))
    aliases <- c(M = "MT", `23` = "X", `24` = "Y", `26` = "MT")
    use <- !is.na(chromosome) & chromosome %in% names(aliases)
    chromosome[use] <- unname(aliases[chromosome[use]])
    if (anyNA(chromosome) || any(!nzchar(chromosome)))
      stop("Cannot assign chromosomes to all manifest genes using the worker's annotation.")
    policy <- list(schema = 1L, exclude = c("X", "Y", "MT"),
                   reason = "User requested exclusion of X/Y/MT genes",
                   genes = data.frame(gene = ctx$manifest$gene, chromosome = chromosome),
                   annotation = source_info, created_at = format(Sys.time(), tz = "UTC", usetz = TRUE))
    em_atomic_write(policy, file.path(ctx$project, "results_slide_prior_em/chromosome_scope.rds"))
    ctx$scope <- sre_scope(ctx$project, ctx$manifest)
  }
  # This also safely finishes an interrupted policy activation. Existing
  # autosomal attempt receipts retain their budgets and fitting settings.
  loaded$state$scope_md5 <- ctx$scope$md5
  em_atomic_write(loaded$state, file.path(ctx$root, "state.rds"))
  audit_file <- file.path(ctx$root, "audit.csv")
  backup <- file.path(ctx$root, "audit_before_chromosome_scope.csv")
  if (!file.exists(backup)) em_atomic_write(read.csv(audit_file), backup, csv = TRUE)
  audit <- sre_refresh(ctx)
  excluded <- data.frame(gene = ctx$scope$excluded,
                          chromosome = unname(ctx$scope$chromosomes[ctx$scope$excluded]),
                          reason = rep("User requested exclusion of X/Y/MT genes", length(ctx$scope$excluded)))
  em_atomic_write(excluded, file.path(ctx$project, "results_slide_prior_em/chromosome_exclusions.csv"), csv = TRUE)
  message("Chromosome scope recorded: ", nrow(excluded), " genes excluded; ",
          sum(!sre_ready(audit$status)), " genes still need recovery.")
  invisible(ctx)
}

sre_pool_files <- function(ctx) {
  files <- sort(list.files(file.path(ctx$dir, "results"), "\\.rds$", full.names = TRUE))
  genes <- sub("\\.rds$", "", basename(files))
  required <- setdiff(ctx$manifest$gene, ctx$scope$excluded)
  if (length(setdiff(required, genes)) || length(setdiff(genes, ctx$manifest$gene)))
    stop("Result files do not match the in-scope gene cohort.")
  files[!genes %in% ctx$scope$excluded]
}

# Exclusions must be explicitly reviewed against the exact saved error record.
# No error-message heuristic is allowed to redefine the initialization cohort.
sre_exclusions <- function(ctx) {
  file <- file.path(ctx$root, "reviewed_exclusions.csv")
  if (!file.exists(file)) return(data.frame(gene = character(), result_md5 = character(), reason = character()))
  x <- read.csv(file, stringsAsFactors = FALSE, colClasses = "character")
  if (!all(c("gene", "result_md5", "reason") %in% names(x)) || anyNA(x) ||
      anyDuplicated(x$gene) || any(!x$gene %in% ctx$manifest$gene) ||
      any(!grepl("^[a-f0-9]{32}$", x$result_md5)) || any(!nzchar(trimws(x$reason))))
    stop("Invalid reviewed_exclusions.csv; supply gene, result_md5 and a review reason.")
  x
}

sre_validate <- function(out, ctx, previous = NULL, require_converged = TRUE) {
  if (!is.list(out) || !identical(as.integer(attr(out, "em_iteration")), as.integer(ctx$iteration)))
    stop("Missing or wrong em_iteration on saved result.")
  if (!is.null(previous)) {
    if (!is.null(previous$error)) {
      if (!identical(out$error, previous$error)) stop("Frozen gene exclusion changed.")
    } else if (!is.null(out$error) || !identical(names(out), names(previous)) ||
               !identical(attr(out, "tissue_errors", exact = TRUE),
                          attr(previous, "tissue_errors", exact = TRUE))) {
      stop("Previously successful gene/tissue cohort changed.")
    }
  }
  if (!is.null(out$error)) {
    if (!is.character(out$error) || length(out$error) != 1L || is.na(out$error) || !nzchar(out$error))
      stop("Invalid gene-error record.")
    return(invisible(TRUE))
  }
  if (length(out) && (is.null(names(out)) || anyNA(names(out)) || anyDuplicated(names(out)) ||
                      any(!names(out) %in% ctx$priors$tissue))) stop("Invalid tissue names.")
  errors <- attr(out, "tissue_errors", exact = TRUE)
  if (length(errors) && (is.null(names(errors)) || anyDuplicated(names(errors)) ||
                         any(!names(errors) %in% ctx$priors$tissue) ||
                         any(names(errors) %in% names(out)) ||
                         !all(vapply(errors, function(x) is.character(x) && length(x) == 1L &&
                                       !is.na(x) && nzchar(x), logical(1))))) stop("Invalid tissue-error metadata.")
  for (tissue in names(out)) {
    fit <- out[[tissue]]$fit_slide_prior
    spe_fit_counts(fit, spe_prior(ctx$priors, tissue), require_converged)
    if (any(!is.finite(fit$alpha_delta)) || any(fit$alpha_delta < 0) ||
        any(abs(apply(fit$alpha_delta, c(1, 2), sum) - fit$alpha) > 1e-7))
      stop("Invalid joint SNP/slider posterior for warm starts.")
    if (!length(fit$elbo) || !is.finite(tail(fit$elbo, 1))) stop("Missing finite ELBO.")
    if (nrow(fit$alpha) != min(ctx$settings$L, ncol(fit$alpha)))
      stop("Saved L differs from frozen settings and predictor count.")
  }
  invisible(TRUE)
}

sre_previous <- function(ctx, gene) {
  if (ctx$iteration == 0L) return(NULL)
  readRDS(file.path(ctx$project, "results_slide_prior_em", sprintf("iteration_%03d", ctx$iteration - 1L),
                    "results", paste0(gene, ".rds")))
}

sre_inspect <- function(ctx, gene, exclusions = sre_exclusions(ctx)) {
  file <- sre_result(ctx, gene)
  row <- data.frame(gene = gene, status = "retry", detail = "Missing result", result_md5 = "",
                    n_fits = 0L, stringsAsFactors = FALSE)
  if (gene %in% ctx$scope$excluded) {
    row$status <- "scope_excluded"
    row$detail <- paste("User-requested exclusion: chromosome", ctx$scope$chromosomes[[gene]])
    return(row)
  }
  if (!file.exists(file)) return(row)
  row$result_md5 <- unname(tools::md5sum(file))
  tryCatch({
    out <- readRDS(file)
    previous <- sre_previous(ctx, gene)
    sre_validate(out, ctx, previous)
    errors <- c(if (!is.null(out$error)) as.character(out$error),
                unlist(attr(out, "tissue_errors", exact = TRUE), use.names = TRUE))
    row$n_fits <- if (is.null(out$error)) length(out) else 0L
    if (length(errors)) {
      approval <- exclusions[exclusions$gene == gene & exclusions$result_md5 == row$result_md5, , drop = FALSE]
      # In later iterations the prior iteration's validated cohort is frozen.
      inherited <- ctx$iteration > 0L && !is.null(previous) &&
        identical(out$error, previous$error) &&
        identical(attr(out, "tissue_errors", exact = TRUE), attr(previous, "tissue_errors", exact = TRUE))
      row$status <- if (nrow(approval) || inherited) "reviewed_exclusion" else "retry"
      row$detail <- paste(errors, collapse = " | ")
      if (nrow(approval)) row$detail <- paste(row$detail, "Review:", approval$reason)
    } else {
      row$status <- if (length(out)) "valid" else "qc_excluded"
      row$detail <- if (length(out)) "Validated converged fits" else "No eligible tissues under existing input QC"
    }
    row
  }, error = function(e) { row$detail <- conditionMessage(e); row })
}

sre_audit <- function(ctx, genes_per_task = 10L, max_attempts = 3L) {
  genes_per_task <- sre_int(genes_per_task, "genes_per_task", 1L, 15L)
  max_attempts <- sre_int(max_attempts, "max_attempts", 1L, 10L)
  state_file <- file.path(ctx$root, "state.rds")
  if (file.exists(state_file)) {
    sre_load(ctx)
    state <- readRDS(state_file)
    if (state$genes_per_task != genes_per_task || state$max_attempts != max_attempts)
      stop("Existing recovery plan settings differ; do not replace its manifest.")
  } else {
    state <- list(schema = 1L, frozen = ctx$frozen, genes_per_task = genes_per_task,
                  max_attempts = max_attempts, scope_md5 = ctx$scope$md5,
                  created_at = format(Sys.time(), tz = "UTC", usetz = TRUE))
    plan <- data.frame(task = ceiling(seq_len(nrow(ctx$manifest)) / genes_per_task),
                       gene = ctx$manifest$gene, original_chunk = ctx$manifest$chunk)
    # State is committed last, so a failed audit can be restarted without trusting partial output.
    em_atomic_write(plan, file.path(ctx$root, "manifest.csv"), csv = TRUE)
  }
  exclusions <- sre_exclusions(ctx)
  rows <- vector("list", nrow(ctx$manifest))
  for (i in seq_len(nrow(ctx$manifest))) {
    rows[[i]] <- sre_inspect(ctx, ctx$manifest$gene[i], exclusions)
    if (i %% 100L == 0L || i == length(rows)) {
      message("Audited ", i, "/", length(rows), " genes")
      gc(verbose = FALSE)
    }
  }
  audit <- do.call(rbind, rows)
  em_atomic_write(audit, file.path(ctx$root, "audit.csv"), csv = TRUE)
  state$manifest_md5 <- unname(tools::md5sum(file.path(ctx$root, "manifest.csv")))
  em_atomic_write(state, state_file)
  message(paste(names(table(audit$status)), table(audit$status), collapse = "; "))
  invisible(audit)
}

sre_load <- function(ctx, allow_scope_adoption = FALSE) {
  state <- readRDS(file.path(ctx$root, "state.rds"))
  if (!identical(state$frozen, ctx$frozen) ||
      !identical(state$manifest_md5, unname(tools::md5sum(file.path(ctx$root, "manifest.csv")))))
    stop("Recovery manifest or frozen files changed.")
  if (!identical(state$scope_md5, ctx$scope$md5) &&
      !(allow_scope_adoption && ctx$iteration == 0L && is.null(state$scope_md5)))
    stop("Chromosome scope changed after recovery audit.")
  plan <- read.csv(file.path(ctx$root, "manifest.csv"), stringsAsFactors = FALSE)
  if (!identical(plan$gene, ctx$manifest$gene) || any(table(plan$task) > state$genes_per_task))
    stop("Recovery plan does not cover the original genes exactly once.")
  list(state = state, plan = plan)
}

# Receipts accelerate intermediate scheduling. Final completion always rereads every fit.
sre_refresh <- function(ctx) {
  audit <- read.csv(file.path(ctx$root, "audit.csv"), stringsAsFactors = FALSE)
  if (!identical(audit$gene, ctx$manifest$gene)) stop("Audit gene universe changed.")
  if (!"attempts" %in% names(audit)) audit$attempts <- 0L
  exclusions <- sre_exclusions(ctx)
  for (i in seq_len(nrow(audit))) {
    if (audit$gene[i] %in% ctx$scope$excluded) {
      row <- sre_inspect(ctx, audit$gene[i], exclusions)
      audit[i, names(row)] <- row
      next
    }
    a <- sre_read(sre_attempt_file(ctx, audit$gene[i]))
    if (is.null(a)) {
      if (file.exists(sre_attempt_file(ctx, audit$gene[i]))) stop("Unreadable attempt receipt: ", audit$gene[i])
      next
    }
    if (!identical(a$frozen, ctx$frozen)) stop("Attempt receipt belongs to different settings.")
    sre_int(a$attempt, "attempt receipt count")
    # Completed receipts are immutable under this scheduler. Recheck all fit
    # files during final verification, rather than rereading old fits each wave.
    if (audit$attempts[i] == a$attempt && audit$status[i] %in% c("valid", "qc_excluded")) next
    audit$attempts[i] <- a$attempt
    if (!is.null(a$result) && a$result$status %in% c("valid", "qc_excluded") &&
        identical(a$result$result_md5,
                                       unname(tools::md5sum(sre_result(ctx, audit$gene[i]))))) {
      audit[i, names(a$result)] <- a$result
    } else {
      inspected <- sre_inspect(ctx, audit$gene[i], exclusions)
      audit[i, names(inspected)] <- inspected
      if (!sre_ready(audit$status[i]) && !is.null(a$error)) audit$detail[i] <- a$error
    }
  }
  em_atomic_write(audit, file.path(ctx$root, "audit.csv"), csv = TRUE)
  audit
}

sre_plan_batch <- function(ctx, capacity, pilot_genes = NULL, automatic = FALSE,
                           genes_per_task = 4L) {
  capacity <- sre_int(capacity, "capacity", 1L, 298L)
  genes_per_task <- sre_int(genes_per_task, "genes_per_task", 1L, 15L)
  loaded <- sre_load(ctx)
  audit <- sre_refresh(ctx)
  pending <- which(!sre_ready(audit$status))
  if (!length(pending)) return(list(status = "verify"))
  if (automatic && file.exists(file.path(ctx$root, "last_batch.txt"))) {
    last <- readLines(file.path(ctx$root, "last_batch.txt"), warn = FALSE)
    if (length(last) != 1L || !identical(dirname(last), paste0(ctx$root, "/batches")))
      stop("Invalid last batch record.")
    saved <- readRDS(file.path(last, "settings.rds"))
    if (!identical(saved$frozen, ctx$frozen) || !identical(saved$scope_md5, ctx$scope$md5) ||
        !identical(saved$manifest_md5, unname(tools::md5sum(file.path(last, "manifest.csv")))))
      stop("Last batch manifest changed.")
    prior_work <- read.csv(file.path(last, "manifest.csv"), stringsAsFactors = FALSE)
    ix <- match(prior_work$gene, audit$gene)
    if (!any(audit$attempts[ix] > prior_work$attempts_before | sre_ready(audit$status[ix]))) {
      message("No worker entered fitting in the last batch; inspect runtime/Slurm logs before manual resume.")
      return(list(status = "blocked"))
    }
  }
  blocked <- pending[audit$attempts[pending] >= loaded$state$max_attempts]
  if (length(blocked)) {
    em_atomic_write(audit[blocked, ], file.path(ctx$root, "blocked.csv"), csv = TRUE)
    return(list(status = "blocked"))
  }
  # A repeated memory/system failure pauses the run before it creates more exclusion records.
  systemic <- pending[audit$attempts[pending] >= 1L &
                        grepl("cannot allocate|out.of.memory|oom|system call failed", audit$detail[pending], ignore.case = TRUE)]
  if (automatic && (length(systemic) >= 2L || any(audit$attempts[systemic] >= 2L))) {
    em_atomic_write(audit[systemic, ], file.path(ctx$root, "blocked.csv"), csv = TRUE)
    return(list(status = "blocked"))
  }
  if (!is.null(pilot_genes)) {
    if (!length(pilot_genes) || anyDuplicated(pilot_genes) || any(!pilot_genes %in% audit$gene))
      stop("Pilot gene list must contain unique genes from the original manifest.")
    pending <- pending[audit$gene[pending] %in% pilot_genes]
    if (!length(pending)) return(list(status = "pilot_empty"))
  }
  # Subdivide only the new batch. The audited plan, previous batch assignments
  # and per-gene receipts remain immutable when the execution cap changes.
  work <- loaded$plan[pending, , drop = FALSE]
  work$attempts_before <- audit$attempts[pending]
  part <- ave(seq_len(nrow(work)), work$task,
              FUN = function(i) ceiling(seq_along(i) / genes_per_task))
  key <- paste(work$task, part, sep = "/")
  groups <- match(key, unique(key))
  attempts <- tapply(work$attempts_before, groups, min)
  selected <- as.integer(names(attempts)[order(attempts, as.integer(names(attempts)))])
  selected <- head(selected, capacity)
  work$array_task <- match(groups, selected)
  work <- work[!is.na(work$array_task), , drop = FALSE]
  # Local array IDs always stay below 299, independent of MaxArraySize for the full pass.
  batch_root <- file.path(ctx$root, "batches")
  dir.create(batch_root, recursive = TRUE, showWarnings = FALSE)
  batch <- tempfile("batch_", tmpdir = batch_root)
  dir.create(batch)
  em_atomic_write(work, file.path(batch, "manifest.csv"), csv = TRUE)
  em_atomic_write(list(frozen = ctx$frozen, scope_md5 = ctx$scope$md5,
                       genes_per_task = genes_per_task,
                       manifest_md5 = unname(tools::md5sum(file.path(batch, "manifest.csv")))),
                  file.path(batch, "settings.rds"))
  list(status = "batch", path = normalizePath(batch, winslash = "/"), tasks = length(selected))
}

sre_task_genes <- function(ctx, batch, task) {
  sre_load(ctx)
  batch <- normalizePath(batch, winslash = "/", mustWork = TRUE)
  if (!identical(dirname(batch), paste0(ctx$root, "/batches"))) stop("Invalid batch path.")
  saved <- readRDS(file.path(batch, "settings.rds"))
  if (!identical(saved$frozen, ctx$frozen) || !identical(saved$scope_md5, ctx$scope$md5) ||
      !identical(saved$manifest_md5, unname(tools::md5sum(file.path(batch, "manifest.csv")))))
    stop("Batch manifest changed.")
  work <- read.csv(file.path(batch, "manifest.csv"), stringsAsFactors = FALSE)
  # Older batches predate the configurable cap and remain readable as saved.
  cap <- if (is.null(saved$genes_per_task)) 15L else
    sre_int(saved$genes_per_task, "batch genes_per_task", 1L, 15L)
  genes <- work$gene[work$array_task == sre_int(task, "array task", 1L, 298L)]
  if (!length(genes) || any(table(work$array_task) > cap) || anyDuplicated(work$gene) ||
      any(!work$gene %in% ctx$manifest$gene)) stop("Invalid batch gene assignment.")
  genes
}

sre_run_gene <- function(ctx, gene, run_gene = run_slide_prior_gene) {
  loaded <- sre_load(ctx)
  if (!gene %in% ctx$manifest$gene) stop("Gene is not in the frozen manifest.")
  before <- sre_inspect(ctx, gene)
  if (sre_ready(before$status)) return(invisible(before))
  gc(verbose = FALSE)
  a <- sre_read(sre_attempt_file(ctx, gene))
  if (is.null(a) && file.exists(sre_attempt_file(ctx, gene))) stop("Unreadable gene attempt receipt.")
  if (!is.null(a)) {
    if (!identical(a$frozen, ctx$frozen)) stop("Attempt receipt belongs to different settings.")
    sre_int(a$attempt, "attempt receipt count")
  }
  attempt <- if (is.null(a)) 1L else a$attempt + 1L
  if (attempt > loaded$state$max_attempts) stop("Gene retry budget exhausted: ", gene)
  record <- list(gene = gene, attempt = attempt, frozen = ctx$frozen,
                 started_at = format(Sys.time(), tz = "UTC", usetz = TRUE), error = "Interrupted before receipt")
  em_atomic_write(record, sre_attempt_file(ctx, gene))
  warnings <- character()
  tryCatch(withCallingHandlers({
    previous <- sre_previous(ctx, gene)
    settings <- ctx$settings
    settings$max_iter <- spe_fit_budget(settings$max_iter)
    if (!is.null(previous$error) || (!is.null(previous) && !length(previous))) {
      out <- previous
    } else {
      args <- c(list(target_gene = gene, tissue_priors = ctx$priors, project_dir = ctx$project,
                     previous_result = previous,
                     temp_dir = file.path(ctx$dir, "temp_plink", "recovery", gene)),
                settings[setdiff(names(settings), c("schema", "package_version", "cohort_source"))])
      out <- do.call(run_gene, args)
    }
    attr(out, "em_iteration") <- ctx$iteration
    attr(out, "fit_max_iter") <- settings$max_iter
    # Keep inspectable nonconverged fits, but never certify them as complete.
    sre_validate(out, ctx, previous, require_converged = FALSE)
    em_atomic_write(out, sre_result(ctx, gene))
    # Release fitted and warm-start objects before reading the checkpoint back.
    # Otherwise verification can hold two copies of each full gene posterior.
    rm(out, previous)
    if (exists("args", inherits = FALSE)) rm(args)
    gc(verbose = FALSE)
    record$result <- sre_inspect(ctx, gene)
    record$error <- if (sre_ready(record$result$status)) NULL else record$result$detail
  }, warning = function(w) { warnings <<- c(warnings, conditionMessage(w)) }), error = function(e) {
    record$error <<- conditionMessage(e)
  })
  record$warnings <- warnings
  if (!is.null(record$error) && length(warnings)) record$error <- paste(record$error, paste(warnings, collapse = " | "), sep = " | ")
  if (!is.null(record$error) && ctx$iteration == 0L && !file.exists(sre_result(ctx, gene))) {
    # An error is saved for review, never counted as completed automatically.
    out <- structure(list(gene = gene, error = record$error), em_iteration = ctx$iteration)
    em_atomic_write(out, sre_result(ctx, gene))
    record$result <- sre_inspect(ctx, gene)
  }
  record$finished_at <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  em_atomic_write(record, sre_attempt_file(ctx, gene))
  if (!is.null(record$error)) stop(gene, ": ", record$error)
  invisible(record$result)
}

sre_verify <- function(ctx) {
  loaded <- sre_load(ctx)
  audit <- sre_audit(ctx, loaded$state$genes_per_task, loaded$state$max_attempts)
  if (any(!sre_ready(audit$status))) stop("Iteration has unresolved genes; no completion certificate written.")
  if (!sum(audit$n_fits)) stop("No successful tissue fits; cannot complete iteration.")
  em_atomic_write(list(frozen = ctx$frozen, scope_md5 = ctx$scope$md5,
                       audit_md5 = unname(tools::md5sum(file.path(ctx$root, "audit.csv"))),
                       completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE)), file.path(ctx$root, "COMPLETE.rds"))
  message("Recovery complete: all genes validated or explicitly audited; no new EM update submitted.")
}

sre_require_complete <- function(project, iteration_dir) {
  ctx <- sre_context(project, iteration_dir)
  certificate <- sre_read(file.path(ctx$root, "COMPLETE.rds"))
  if (is.null(certificate) || !identical(certificate$frozen, ctx$frozen) ||
      !identical(certificate$scope_md5, ctx$scope$md5) ||
      !identical(certificate$audit_md5, unname(tools::md5sum(file.path(ctx$root, "audit.csv")))))
    stop("Recovery is unfinished or its audit changed; finish with recover_slide_prior_em before advancing.")
  # Revalidate all results immediately before allowing the existing scientific M-step.
  sre_verify(ctx)
  invisible(TRUE)
}
