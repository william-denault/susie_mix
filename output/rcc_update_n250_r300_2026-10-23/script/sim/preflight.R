# Read-only RCC validation: packages, design, genotype sources and job manifest.
local({
  files <- vapply(sys.frames(), function(f) if (is.null(f$ofile)) "" else f$ofile, "")
  entry <- if (any(nzchar(files))) tail(files[nzchar(files)], 1L) else
    sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
  root <- dirname(dirname(dirname(normalizePath(entry, winslash = "/", mustWork = TRUE))))
  source(file.path(root, "script/sim/run_job.R"), local = TRUE)
  manifest <- read.csv(file.path(root, "script/sim/jobs_slide/manifest.csv"))
  expected <- sim_conditions()
  stopifnot(nrow(manifest) == 900L, all(manifest$n == 250L),
    all(manifest$reps_per_chunk == 300L), all(manifest$L == 10L),
    all(manifest$schema_version == sim_schema_version),
    identical(sort(unique(manifest$pve)), c(.1, .2, .3, .4)),
    setequal(unique(manifest$K), 1:5),
    nrow(unique(manifest[c("name", sim_count_columns)])) == nrow(expected),
    !anyDuplicated(manifest$output_file),
    all(startsWith(manifest$output_file, paste0(sim_results_dir, "/chunks/"))))
  for (package in c("susieR", "susieRSlidePrior", "data.table", "matrixStats")) {
    if (!requireNamespace(package, quietly = TRUE)) stop("Install missing package: ", package)
    cat(package, as.character(utils::packageVersion(package)), "\n")
  }
  stopifnot(all(c("min_obs", "delta_grid", "delta_prior") %in%
                  names(formals(susieRSlidePrior::susie))))
  source <- prepare_simulation_genotypes(root)
  if (!identical(source$mode, "gtex") || nzchar(source$settings$target_gene))
    stop("This benchmark requires random GTEx genes: unset SUSIE_MIX_SIM_INPUT and SUSIE_MIX_SIM_GENE.")
  stopifnot(all(file.exists(file.path(root, "script/sim/jobs_slide",
                                    paste0("sim_job_", manifest$job_id, ".R")))))
  cat("PASS: 225 allocations, 900 conditions, 270000 datasets, 810000 fits.\n")
  cat("Genotype source:", nrow(source$loci), "GTEx gene regions; no extraction or fitting performed.\n")
  cat("Output:", file.path(root, sim_results_dir), "\n")
})
