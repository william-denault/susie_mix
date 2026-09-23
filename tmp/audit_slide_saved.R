manifest <- read.csv('script/sim/jobs_slide/manifest.csv', stringsAsFactors = FALSE)
present <- file.exists(manifest$output_file)
cat('CHECKPOINT INVENTORY:', sum(present), '/', nrow(manifest), '\n')
cat('MISSING JOB IDS:', paste(manifest$job_id[!present], collapse = ','), '\n')
selected_jobs <- suppressWarnings(as.integer(commandArgs(trailingOnly = TRUE)))
if (length(selected_jobs)) {
  present <- present & manifest$job_id %in% selected_jobs
  cat('SELECTED CHECKPOINT AUDIT:', paste(manifest$job_id[present], collapse = ','), '\n')
}
methods <- c('SuSiE', 'SuSiE-mix', 'SuSiE-slide')
totals <- setNames(rep(0, 9), c('files', 'records', 'success', 'errors', 'complete', 'error_files', 'clean_files', 'invalid', 'nonconverged_replicates'))
nonconv <- setNames(rep(0, 3), methods)
error_counts <- integer()
rows <- list()
raw_by_seed <- list()
for (i in which(present)) {
  env <- new.env(parent = emptyenv())
  problem <- tryCatch({ load(manifest$output_file[i], envir = env); NULL }, error = conditionMessage)
  if (!is.null(problem)) {cat('UNREADABLE:', i, problem, '\n'); next}
  results <- env$results
  errors <- vapply(results, function(x) !is.null(x$error), logical(1))
  ok <- results[!errors]
  valid <- vapply(ok, function(x) identical(x$metrics$method, methods) &&
    all(c('susie_slide_pip', 'susie_slide_cs', 'susie_slide_delta_causal', 'susie_slide_delta_forced') %in% names(x)), logical(1))
  badconv <- vapply(ok, function(x) any(!x$metrics$converged | is.na(x$metrics$converged)), logical(1))
  for (x in ok) {
    nonconv <- nonconv + as.integer(!x$metrics$converged | is.na(x$metrics$converged))
    seed <- as.character(x$seed)
    raw_by_seed[[seed]] <- union(raw_by_seed[[seed]], x$raw_file)
  }
  for (x in results[errors]) {
    msg <- x$error
    if (!msg %in% names(error_counts)) error_counts[msg] <- 0L
    error_counts[msg] <- error_counts[msg] + 1L
  }
  nr <- length(results); ne <- sum(errors)
  totals <- totals + c(1, nr, sum(!errors), ne, nr == 400, ne > 0, nr == 400 && ne == 0, sum(!valid), sum(badconv))
  rows[[length(rows) + 1L]] <- data.frame(job_id = i, records = nr, success = sum(!errors), errors = ne, nonconverged = sum(badconv))
  if (length(rows) == 1L) {cat('PACKAGES:\n'); print(env$checkpoint_settings$package_versions); cat('INPUT:', env$checkpoint_settings$arguments$temp_dir, '\n')}
  if (length(rows) %% 50L == 0L) {cat('PROGRESS', length(rows), '\n'); flush.console()}
  rm(env, results, ok)
}
cat('FINAL TOTALS:\n'); print(totals)
cat('NONCONVERGENCE BY METHOD:\n'); print(nonconv)
cat('ERROR TYPES:\n'); print(sort(error_counts, decreasing = TRUE))
audit <- do.call(rbind, rows)
if (length(selected_jobs)) {cat('SELECTED CHECKPOINT DETAILS:\n'); print(audit)}
cat('BY BATCH:\n'); print(aggregate(audit[-1], list(batch = (audit$job_id - 1L) %/% 300L + 1L), sum))
cat('RANGE OF UNIQUE GENOTYPE FILES PER SEED ACROSS CONDITIONS:\n'); print(range(lengths(raw_by_seed)))
cat('INCOMPLETE CHECKPOINTS:\n'); print(audit[audit$records < 400, ])
