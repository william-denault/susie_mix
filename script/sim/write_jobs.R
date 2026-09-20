# Generate the job manifest; this does not fit models or submit cluster jobs.
.job_sources <- vapply(sys.frames(), function(f) if (is.null(f$ofile)) "" else f$ofile, "")
.job_source <- if (any(nzchar(.job_sources))) tail(.job_sources[nzchar(.job_sources)], 1L) else
  sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
.job_source <- normalizePath(.job_source, winslash = "/", mustWork = TRUE)
source(file.path(dirname(.job_source), "simulation_design.R"), local = TRUE)
sim_project_dir <- dirname(dirname(dirname(.job_source)))
rm(.job_source, .job_sources)

write_simulation_jobs <- function(project_dir = sim_project_dir,
                                  pve_values = c(.05, .10, .20, .30, .40),
                                  n = 500L, L = 10L, chunks_per_cell = 1L,
                                  reps_per_chunk = 400L, seed_base = 1000000L,
                                  array_batch_size = 400L) {
  stopifnot(length(pve_values) > 0L, !anyDuplicated(pve_values),
            all(pve_values > 0 & pve_values < 1), n >= 3, L >= 1,
            chunks_per_cell >= 1, reps_per_chunk >= 1,
            length(array_batch_size) == 1L, is.finite(array_batch_size),
            array_batch_size >= 1L, array_batch_size <= 400L,
            array_batch_size == floor(array_batch_size),
            all(c(n, L, chunks_per_cell, reps_per_chunk, seed_base) ==
                  floor(c(n, L, chunks_per_cell, reps_per_chunk, seed_base))))
  conditions <- sim_conditions()
  jobs <- list()
  for (i in seq_len(nrow(conditions))) {
    for (pve in pve_values) {
      for (chunk in seq_len(chunks_per_cell)) {
        condition <- conditions[i, , drop = FALSE]
        jobs[[length(jobs) + 1L]] <- data.frame(
          job_id = length(jobs) + 1L, schema_version = sim_schema_version,
          condition, n = n, L = L, pve = pve, seed_base = seed_base,
          reps_per_chunk = reps_per_chunk, chunk = chunk,
          min_maf = .05, hwe_thresh = 1e-8, min_n_rec = 5L, slide_min_obs = 5L,
          delta_prec = -.5, delta_pdom = .5,
          output_file = paste0("simulation results/slide_v1/chunks/",
            sim_checkpoint_name(condition, n, L, pve, seed_base, reps_per_chunk, chunk)))
      }
    }
  }
  manifest <- do.call(rbind, jobs)
  stopifnot(!anyDuplicated(manifest$output_file))
  job_dir <- file.path(project_dir, "script/sim/jobs_slide")
  dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(conditions, file.path(job_dir, "conditions.csv"), row.names = FALSE)
  write.csv(manifest, file.path(job_dir, "manifest.csv"), row.names = FALSE)
  offsets <- seq.int(0L, nrow(manifest) - 1L, by = array_batch_size)
  sizes <- pmin(array_batch_size, nrow(manifest) - offsets)
  batches <- data.frame(batch = seq_along(offsets), array_start = 0L,
                        array_end = sizes - 1L, offset = offsets,
                        first_job = offsets + 1L, last_job = offsets + sizes)
  write.csv(batches, file.path(job_dir, "submission_batches.csv"), row.names = FALSE)
  cat("Generated", nrow(conditions), "configurations and", nrow(manifest),
      "jobs (", sum(manifest$reps_per_chunk), "replicates).\n")
  cat("Submit ONE batch at a time; wait for it to finish before submitting the next:\n")
  for (i in seq_len(nrow(batches))) {
    cat("  sbatch --array=0-", batches$array_end[i], " job/run_simulation ",
        batches$offset[i], "  # manifest jobs ", batches$first_job[i], "-",
        batches$last_job[i], "\n", sep = "")
  }
  invisible(manifest)
}

if (sys.nframe() == 0L) write_simulation_jobs()
