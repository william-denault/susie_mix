# Integration fixture from the count/history tests; no real scan files touched.
source("script/scan_tissue_attempt/tests/test_slide_prior_parallel.R")
source("script/scan_tissue_attempt/slide_prior_run_control.R")
original_target <- Sys.getenv("SUSIE_SLIDE_TARGET_ITERATION", unset = NA_character_)
Sys.setenv(SUSIE_SLIDE_TARGET_ITERATION = "20")
stopifnot(spt_target(root, 1L) == 20L)
frozen_files <- c(file.path(root, "results_slide_prior_em/prior_history.csv"),
                  list.files(launch$iteration_dir, full.names = TRUE, pattern = "\\.(csv|rds)$"),
                  list.files(file.path(launch$iteration_dir, "results"), full.names = TRUE))
before <- tools::md5sum(frozen_files)
spt_set_target(root, 1L, "40")
stopifnot(spt_target(root, 1L) == 40L, identical(before, tools::md5sum(frozen_files)),
          Sys.getenv("SUSIE_SLIDE_TARGET_ITERATION") == "20")
spt_set_target(root, 1L, "40")
fails(spt_set_target(root, 1L, "39"), "only extends")
fails(spt_set_target(root, 1L, "40.5"), "integer")
saved <- readRDS(spt_file(root)); saveRDS(list(schema = 7L), spt_file(root))
fails(spt_target(root, 1L), "Invalid")
saveRDS(saved, spt_file(root))
result <- spt_audit(root)
stopifnot(result$iteration == 1L)
# Catch corruption independently from the production M-step itself.
counts_file <- file.path(launch$iteration_dir, "pooled_counts.csv")
counts <- read.csv(counts_file); bad_counts <- counts
bad_counts$w_01[1] <- bad_counts$w_01[1] + 5
write.csv(bad_counts, counts_file, row.names = FALSE)
fails(spt_audit(root), "M-step count/weight mismatch")
write.csv(counts, counts_file, row.names = FALSE)
# Both old inherited environment and existing saved data are retained by target changes.
if (is.na(original_target)) Sys.unsetenv("SUSIE_SLIDE_TARGET_ITERATION") else
  Sys.setenv(SUSIE_SLIDE_TARGET_ITERATION = original_target)
cat("PASS: 20 -> 40 override, inherited environment, idempotence, invalid/decreased targets, unchanged fits/history, metadata audit and corrupt pooled-count rejection.\n")
