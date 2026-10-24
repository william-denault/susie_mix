# Read-only diagnosis of checkpoint compatibility. No fitting, writing or submission.
# From the project root: Rscript --vanilla script/sim/audit_checkpoint_resume.R 17 72

checkpoint_field_differences <- function(saved, current, field = "settings") {
  if (identical(saved, current)) return(character())
  describe <- function(x) {
    if (is.null(x)) return("NULL")
    if (is.atomic(x)) {
      shown <- paste(capture.output(dput(head(x, 5L))), collapse = " ")
      return(paste0(typeof(x), "[", length(x), "]: ", shown))
    }
    paste0(typeof(x), "[", length(x), "]")
  }
  if (is.list(saved) && is.list(current) &&
      identical(names(saved), names(current)) && length(saved) == length(current)) {
    out <- character()
    if (!identical(attributes(saved), attributes(current)))
      out <- c(out, checkpoint_field_differences(attributes(saved), attributes(current),
                                                paste0(field, " [attributes]")))
    for (i in seq_along(saved)) {
      label <- if (is.null(names(saved))) paste0("[[", i, "]]") else paste0("$", names(saved)[i])
      out <- c(out, checkpoint_field_differences(saved[[i]], current[[i]], paste0(field, label)))
    }
    return(out)
  }
  out <- paste0(field, "\n  saved:   ", describe(saved), "\n  current: ", describe(current))
  if (is.atomic(saved) && is.atomic(current) && length(saved) == length(current) &&
      typeof(saved) == typeof(current)) {
    changed <- which(!vapply(seq_along(saved), function(i) identical(saved[i], current[i]), logical(1)))
    if (length(changed)) {
      selected <- head(changed, 8L)
      out <- c(out, paste0("  changed positions (first 8 of ", length(changed), "): ",
        paste(selected, collapse = ", ")),
        paste0("  saved at those positions: ", describe(saved[selected])),
        paste0("  current at those positions: ", describe(current[selected])))
    }
    if (!identical(attributes(saved), attributes(current)))
      out <- c(out, paste0("  attribute differences: ", paste(all.equal(attributes(saved),
        attributes(current), tolerance = 0), collapse = "; ")))
  }
  out
}

audit_checkpoint_resume <- function(project_dir, job_ids) {
  project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)
  runtime <- new.env(parent = globalenv())
  source(file.path(project_dir, "script/sim/run_job.R"), local = runtime)
  packages <- c("susieR", "susieRSlidePrior", "data.table", "matrixStats")
  for (p in packages) if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p)
  versions <- setNames(vapply(packages, function(p) as.character(utils::packageVersion(p)), ""), packages)
  genotype_source <- runtime$prepare_simulation_genotypes(project_dir)
  source_paths <- c("script/sim/sim_workhorse.R", "script/sim/simulation_design.R",
    "script/sim/run_job.R", "script/sim/additive_init_genotypes.R",
    "script/scan_tissue_attempt/workhorse_utils.R", "script/scan_tissue_attempt/get_gene_annotations.R")
  argument_names <- c("pve", "n", "L", runtime$sim_count_columns, "all_additive",
                       "min_maf", "hwe_thresh", "min_n_rec", "slide_min_obs")
  for (job_id in job_ids) {
    cat("\nCONDITION ", job_id, "\n", sep = "")
    # Read the same embedded settings that the numbered Slurm driver uses.
    expressions <- as.list(parse(file.path(project_dir, "script/sim/jobs_slide",
                                           paste0("sim_job_", job_id, ".R"))))
    assignments <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
                            identical(x[[2]], as.name("job")), expressions)
    if (length(assignments) != 1L) stop("Cannot identify embedded settings for job ", job_id)
    row <- eval(assignments[[1]][[3]], envir = runtime)
    stopifnot(is.data.frame(row), nrow(row) == 1L, row$job_id == job_id)
    args <- as.list(row[argument_names])
    if (identical(genotype_source$mode, "raw")) args$temp_dir <- genotype_source$directory
    current <- list(design = row, arguments = args, package_versions = versions,
      methods = runtime$sim_methods, delta_grid = runtime$sim_delta_grid,
      delta_prior = runtime$sim_delta_prior, genotype_inputs = genotype_source$settings,
      R_version = R.version.string, RNG_kind = RNGkind(),
      source_md5 = unname(tools::md5sum(file.path(project_dir, source_paths))))
    checkpoint <- file.path(project_dir, row$output_file)
    if (!file.exists(checkpoint)) {
      cat("Checkpoint not present: ", checkpoint, "\n", sep = "")
      next
    }
    saved <- new.env(parent = emptyenv())
    load(checkpoint, envir = saved)
    cat("Saved replicate records: ", length(saved$results), "\n", sep = "")
    differences <- checkpoint_field_differences(
      runtime$sim_comparable_checkpoint(saved$checkpoint_settings),
      runtime$sim_comparable_checkpoint(current))
    if (!length(differences)) cat("MATCH: checkpoint metadata matches the current runtime.\n")
    else {
      cat("MISMATCH:\n", paste(differences, collapse = "\n"), "\n", sep = "")
      if (any(grepl("genotype_inputs\\$files", differences))) {
        cat("Input signature row mapping:\n")
        print(current$genotype_inputs$files, row.names = TRUE)
      }
      if (any(grepl("source_md5", differences))) {
        cat("Source hash position mapping:\n")
        print(data.frame(position = seq_along(source_paths), path = source_paths), row.names = FALSE)
      }
    }
  }
  cat("\nRead-only audit finished. No checkpoints or scripts were changed; no jobs submitted.\n")
  invisible(NULL)
}

if (sys.nframe() == 0L) {
  entry <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
  root <- dirname(dirname(dirname(normalizePath(entry, winslash = "/", mustWork = TRUE))))
  ids <- commandArgs(trailingOnly = TRUE)
  if (!length(ids)) ids <- c("17", "72")
  if (any(!grepl("^[1-9][0-9]*$", ids)) || any(as.numeric(ids) > 900))
    stop("Supply condition numbers between 1 and 900.")
  audit_checkpoint_resume(root, as.integer(ids))
}
