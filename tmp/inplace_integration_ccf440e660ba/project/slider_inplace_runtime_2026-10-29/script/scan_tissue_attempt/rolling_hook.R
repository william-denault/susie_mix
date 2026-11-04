# Called only by the controller while its exclusive .prepare.lock is held.
# This barrier also runs on a resumed prepared iteration, before any planning.
sr_barrier <- function(ctx, queue_file) {
  receipt <- file.path(ctx$project,"results_slide_prior_em/rolling_cleanup_inplace",
                       sprintf("before_%03d",ctx$iteration),"COMPLETE.json")
  if(file.exists(receipt)) return(invisible(NULL))
  if(ctx$iteration < 6L) stop("Invalid expanded iteration")
  args <- c(shQuote(file.path(ctx$project,"slider_inplace_runtime_2026-10-29/retention.py")),"step",shQuote(ctx$project),
            as.character(ctx$iteration),shQuote(queue_file))
  output <- system2(Sys.getenv("PYTHON","python3"),args,stdout=TRUE,stderr=TRUE)
  if(length(output)) message(paste(output,collapse="\n"))
  status <- attr(output,"status"); if(is.null(status)) status <- 0L
  if(status != 0L) stop("Rolling cleanup is unfinished. No workers submitted; resume latest after inspecting the log.")
  if(!file.exists(receipt)) stop("Missing cleanup receipt")
  invisible(NULL)
}
