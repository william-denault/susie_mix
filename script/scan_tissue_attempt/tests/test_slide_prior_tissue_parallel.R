# Real-fitter integration plus interrupted tissue work and scheduler ownership.
source("script/scan_tissue_attempt/tests/test_slide_prior_em_worker.R")
repo <- normalizePath(".", winslash = "/")
for (f in c("slide_prior_parallel.R", "slide_prior_tissue_parallel.R"))
  source(file.path(repo, "script/scan_tissue_attempt", f))
source(file.path(repo, "script/scan_tissue_attempt/workhorse_slide_prior_em.R"))
invisible(file.copy(list.files(file.path(repo, "script/scan_tissue_attempt"), "\\.R$", full.names = TRUE),
                    scripts, overwrite = TRUE))
fails <- function(expr, pattern) {
  e <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(e, "error"), grepl(pattern, conditionMessage(e)))
}
queue <- file.path(project, "queue.txt"); writeLines(character(), queue)
for (context in list(ctx, ctx1)) {
  spp_init(context)
  spp_inspect(context, "GENE", force = TRUE)
  stopifnot(spp_plan(context, 298L, queue)$status == "complete")
}
launch <- spe_prepare_iteration(project, "new", completion_reader = spp_require_complete)
ctx2 <- sre_context(project, launch$iteration_dir, check_package = TRUE); spp_init(ctx2)
previous <- sre_previous(ctx2, "GENE")
call_args <- c(list(target_gene = "GENE", tissue_priors = ctx2$priors, previous_result = previous),
               ctx2$settings[setdiff(names(ctx2$settings), c("schema", "package_version", "cohort_source"))])
expected <- do.call(run_real, call_args)
# Compare the refactored whole-gene path against the actually deployed worker.
old_env <- new.env(parent = globalenv()); old_env$system <- worker_env$system
connection <- unz(file.path(repo, "output/slide_prior_shared_rcc_2026-10-23.zip"),
                  "script/scan_tissue_attempt/workhorse_slide_prior_em.R")
source(connection, local = old_env); close(connection)
old_args <- call_args
paths <- args[c("datadir", "gene_annot_fun", "subject_pheno_file", "sample_attr_file", "expr_file")]
old_args[names(paths)] <- paths; old_args$project_dir <- repo; old_args$temp_dir <- file.path(root, "old_plink")
expected_old <- do.call(old_env$run_slide_prior_gene, old_args)
stopifnot(isTRUE(all.equal(expected, expected_old, tolerance = 1e-12)))

sps_enable(ctx2, "GENE")
bad <- ctx2; bad$scope$excluded <- "GENE"
fails(sps_enable(bad, "GENE"), "excluded by the frozen chromosome")
fails(sps_enable(ctx, "GENE"), "initial gene/tissue cohort")
register <- function(plan, id) {
  writeLines(id, file.path(plan$batch, "array_job_id.txt"))
  for (task in seq_len(plan$tasks)) spp_task(ctx2, plan$batch, task)
}
# An already submitted whole-gene task retains exclusive ownership.
whole <- spp_write_batch(ctx2, data.frame(task = 1L, gene = "GENE", action = "fit"))
register(whole, "100")
writeLines("100_1|100|slide_parallel_fit", queue)
stopifnot(spp_plan(ctx2, 2L, queue)$status == "waiting")
writeLines(character(), queue)
prepare <- spp_plan(ctx2, 2L, queue); register(prepare, "101")
stopifnot(spp_batch(ctx2, prepare$batch)$unit == "prepare")
sps_worker(ctx2, "GENE", "prepare", run_gene = run_real)
m <- sps_manifest(ctx2, "GENE")
stopifnot(identical(m$tissues, names(previous)), length(m$units) == 2L,
          !file.exists(sre_result(ctx2, "GENE")))

first <- spp_plan(ctx2, 1L, queue); register(first, "102")
stopifnot(first$tasks == 1L, spp_batch(ctx2, first$batch)$unit == "tissue_001")
writeLines("102_1|102|slide_parallel_fit", queue)
second <- spp_plan(ctx2, 1L, queue); register(second, "103")
stopifnot(spp_batch(ctx2, second$batch)$unit == "tissue_002")
writeLines(c("102_1|102|slide_parallel_fit", "103_1|103|slide_parallel_fit"), queue)
stopifnot(spp_queue(ctx2, queue)$tasks == 2L, spp_plan(ctx2, 298L, queue)$status == "waiting")
fails(spp_check_units(data.frame(gene = c("GENE", "GENE"), unit = c("gene", "tissue_001"))), "Overlapping")
sps_worker(ctx2, "GENE", "tissue_001")
good_paths <- file.path(sps_dir(ctx2, "GENE"), c("results/tissue_001.rds", "completed/tissue_001.rds"))
hashes <- tools::md5sum(good_paths)
fails(sps_worker(ctx2, "GENE", "tissue_002", fit_tissue = function(...) stop("simulated interruption")), "simulated interruption")
fails(sps_assemble(ctx2, "GENE", m), "Unfinished tissues")
stopifnot(!file.exists(sre_result(ctx2, "GENE")), !file.exists(file.path(spp_root(ctx2), "COMPLETE.rds")))
writeLines(character(), queue)
retry <- spp_plan(ctx2, 298L, queue); register(retry, "104")
stopifnot(retry$tasks == 1L, spp_batch(ctx2, retry$batch)$unit == "tissue_002")
not_converged <- function(...) {
  result <- spe_fit_prepared_slide_tissue(...)
  result$fit_slide_prior$converged <- FALSE
  result
}
fails(sps_worker(ctx2, "GENE", "tissue_002", fit_tissue = not_converged), "Unconverged")
stopifnot(!sps_done(ctx2, "GENE", "tissue_002", m))
receipt_path <- sps_attempt_path(ctx2, "GENE", "tissue_002")
receipt <- readRDS(receipt_path); exhausted <- receipt; exhausted$attempt <- 3L
saveRDS(exhausted, receipt_path)
fails(spp_plan(ctx2, 298L, queue), "retry budget exhausted")
saveRDS(receipt, receipt_path)
sps_worker(ctx2, "GENE", "tissue_002")
sps_worker(ctx2, "GENE", "tissue_001") # Idempotent, without another fitting attempt.
stopifnot(identical(hashes, tools::md5sum(good_paths)),
          readRDS(sps_attempt_path(ctx2, "GENE", "tissue_001"))$attempt == 1L,
          readRDS(sps_attempt_path(ctx2, "GENE", "tissue_002"))$attempt == 3L)
assembly <- spp_plan(ctx2, 298L, queue); register(assembly, "105")
stopifnot(spp_batch(ctx2, assembly$batch)$unit == "assemble")
sps_worker(ctx2, "GENE", "assemble")
assembled <- readRDS(sre_result(ctx2, "GENE"))
for (tissue in names(expected))
  stopifnot(isTRUE(all.equal(assembled[[tissue]], expected[[tissue]], tolerance = 1e-12)))
attr(expected, "em_iteration") <- ctx2$iteration
serial <- spe_summarize_gene(expected, "GENE", ctx2$priors)
split <- spe_summarize_gene(assembled, "GENE", ctx2$priors)
stopifnot(isTRUE(all.equal(serial, split, tolerance = 1e-12)),
          spp_plan(ctx2, 298L, queue)$status == "complete")
f <- file.path(sps_dir(ctx2, "GENE"), "inputs/tissue_001.rds")
writeLines("damaged", f)
fails(sps_manifest(ctx2, "GENE"), "inputs changed")
cat("PASS: deployed-vs-refactored worker, real whole-gene-vs-tissue fits and EM counts, frozen cohort, concurrent tissues, whole-gene ownership, partial retry, idempotence, assembly barrier and changed-input detection.\n")
