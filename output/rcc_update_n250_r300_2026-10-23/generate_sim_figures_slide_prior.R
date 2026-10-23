# Compatibility entry point for the new three-method benchmark.
local({
  files <- vapply(sys.frames(), function(f) if (is.null(f$ofile)) "" else f$ofile, "")
  entry <- if (any(nzchar(files))) tail(files[nzchar(files)], 1L) else
    sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
  source(file.path(dirname(normalizePath(entry, winslash = "/", mustWork = TRUE)),
                   "generate_sim_figures.R"), local = TRUE)
})
