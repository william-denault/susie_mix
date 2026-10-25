# Base-R integration fixtures. Never reads or writes the real GTEx results.
repo <- normalizePath(".", winslash = "/")
dir.create(file.path(repo, "tmp"), showWarnings = FALSE)
for (f in c("em_utils.R", "slide_prior_em_utils.R", "slide_prior_recovery.R",
            "prepare_slide_prior_em_iteration.R", "slide_prior_parallel.R"))
  source(file.path(repo, "script/scan_tissue_attempt", f))
fails <- function(expr, pattern) {
  e <- tryCatch({force(expr); NULL}, error = identity)
  stopifnot(inherits(e, "error"), grepl(pattern, conditionMessage(e)))
}
root <- tempfile("parallel_em_", tmpdir = file.path(repo, "tmp")); dir.create(root)
scripts <- file.path(root, "script/scan_tissue_attempt"); dir.create(scripts, recursive = TRUE)
file.copy(file.path(repo, "script/scan_tissue_attempt/slide_prior_recovery.R"), scripts)
d <- file.path(root, "results_slide_prior_em/iteration_000"); dir.create(file.path(d, "results"), recursive = TRUE)
genes <- sprintf("G%03d", 1:55)
manifest <- data.frame(chunk = 1L, gene = genes)
w <- matrix(1/17, 2, 17, dimnames = list(NULL, spe_columns()))
p <- data.frame(iteration = 0L, source_iteration = -1L, tissue = c("Brain", "Liver"), w,
                n_fits = 0L, n_active = 0L, free_mass = 0, prior_max_change = 0,
                m_step_q_gain = 0, source_elbo = NA_real_,
                update_method = "slider_prior_active_counts_v1", created_at = "fixture")
write.csv(p, file.path(d, "priors.csv"), row.names = FALSE)
write.csv(p, file.path(root, "results_slide_prior_em/prior_history.csv"), row.names = FALSE)
write.csv(manifest, file.path(d, "manifest.csv"), row.names = FALSE)
write.csv(data.frame(delta = spe_grid()), file.path(d, "grid.csv"), row.names = FALSE)
saveRDS(c(spe_settings(), list(package_version = "fixture")), file.path(d, "settings.rds"))
scope <- list(schema = 1L, exclude = c("X", "Y", "MT"),
              genes = data.frame(gene = genes, chromosome = c(rep("1", 54), "X")))
saveRDS(scope, file.path(root, "results_slide_prior_em/chromosome_scope.rds"))
ctx <- sre_context(root, d)
fit <- function(prior = rep(1/17, 17)) {
  # Forced SNP, null candidate and V=0 component must not add slider counts.
  a <- rbind(c(.2,.5,.3), c(.8,.1,.1), c(.5,.3,.2))
  q <- matrix(0, 3, 17); q[1,1] <- 1; q[2,17] <- 1; q[3,9] <- 1
  joint <- array(0, c(3L,3L,17L))
  for (l in 1:3) { joint[l,1,] <- a[l,1]*q[l,]; joint[l,2:3,9] <- a[l,2:3] }
  counts <- a[,1]*q; counts[3,] <- 0
  structure(list(alpha=a, alpha_delta=joint, delta_prior_counts=counts,
    delta_forced=c(FALSE,TRUE,TRUE), input_p=2L, delta_grid=spe_grid(), delta_prior=prior,
    V=c(1,1,0), converged=TRUE, elbo=-100, sets=NULL), class=c("susie_slide","susie"))
}
out <- function(prior=rep(1/17,17), iter=0L) structure(list(Brain=list(fit_slide_prior=fit(prior))), em_iteration=iter)
for (g in genes[1:53]) saveRDS(out(), sre_result(ctx,g))
saveRDS(structure(list(), em_iteration=0L), sre_result(ctx,"G054"))
# Existing successful/failing files on excluded chromosomes remain untouched.
saveRDS(structure(list(error="old sex-chromosome error"), em_iteration=0L), sre_result(ctx,"G055"))
files <- list.files(file.path(d,"results"), full.names=TRUE)
hashes <- tools::md5sum(c(files,file.path(d,names(ctx$frozen))))
spp_init(ctx)
state_hash <- tools::md5sum(file.path(ctx$root,c("state.rds","manifest.csv")))
queue <- file.path(root,"queue.txt"); writeLines(character(), queue)
plan <- spp_plan(ctx,2L,queue)
stopifnot(plan$action=="validate", plan$tasks==2L)
work <- spp_batch(ctx,plan$batch)
stopifnot(nrow(work)==50L, max(table(work$task))==25L, !"G055" %in% work$gene)
writeLines("100",file.path(plan$batch,"array_job_id.txt"))
writeLines("100_1|100|slide_parallel_validate",queue)
# Finish task 2 while task 1 still owns its entire group.
for(g in spp_task(ctx,plan$batch,2)$gene) spp_worker_gene(ctx,g,"validate")
next_plan <- spp_plan(ctx,1L,queue)
next_work <- spp_batch(ctx,next_plan$batch)
stopifnot(!any(next_work$gene %in% work$gene), nrow(next_work)==4L)
for(g in spp_task(ctx,next_plan$batch,1)$gene) spp_worker_gene(ctx,g,"validate")
# Simulate an interrupted validation task: first three receipts survive.
first <- spp_task(ctx,plan$batch,1)$gene
for(g in first[1:3]) spp_worker_gene(ctx,g,"validate")
saved <- tools::md5sum(vapply(first[1:3],function(g)spp_summary_file(ctx,g),character(1)))
writeLines(character(),queue)
resumed <- spp_plan(ctx,298L,queue)
stopifnot(resumed$action=="validate")
rw <- spp_batch(ctx,resumed$batch)
stopifnot(setequal(rw$gene,first[-(1:3)]))
for(g in spp_task(ctx,resumed$batch,1)$gene) spp_worker_gene(ctx,g,"validate")
stopifnot(identical(saved,tools::md5sum(names(saved))),
          identical(hashes,tools::md5sum(names(hashes))),
          identical(state_hash,tools::md5sum(names(state_hash))))
expected <- spe_pool(sre_pool_files(ctx),ctx$priors)
# Prove coordinator completion and the M-step do not load full posteriors.
real_read <- readRDS
readRDS <- function(file,...) {
  if(is.character(file) && grepl("/results/",gsub("\\\\","/",file),fixed=TRUE)) stop("Coordinator read a full posterior")
  real_read(file,...)
}
done <- spp_plan(ctx,298L,queue)
pooled <- spp_require_complete(root,d)
stopifnot(done$status=="complete", isTRUE(all.equal(pooled,expected)),
          abs(pooled$priors$w_01[1]-.2)<1e-12, abs(pooled$priors$w_17[1]-.8)<1e-12,
          all(pooled$priors[2,spe_columns()]==ctx$priors[2,spe_columns()]),
          pooled$priors$n_fits[1]==53L, abs(pooled$priors$free_mass[1]-53)<1e-10,
          !file.exists(file.path(root,"results_slide_prior_em/objective_history.csv")))
launch <- spe_prepare_iteration(root,"new",package_version="fixture",completion_reader=spp_require_complete)
rm(readRDS)
ctx1 <- sre_context(root,launch$iteration_dir); spp_init(ctx1)
stopifnot(ctx1$iteration==1L, isTRUE(all.equal(spe_prior(ctx1$priors,"Brain"),as.numeric(expected$priors[1,spe_columns()]))))
fits <- spp_plan(ctx1,2L,queue); fw <- spp_batch(ctx1,fits$batch)
stopifnot(fits$action=="fit",nrow(fw)==2L,all(table(fw$task)==1L))
calls <- 0L
runner <- function(target_gene,tissue_priors,previous_result,...) {
  calls <<- calls+1L
  stopifnot(!is.null(previous_result$Brain$fit_slide_prior))
  out(spe_prior(tissue_priors,"Brain"),1L)
}
spp_worker_gene(ctx1,fw$gene[1],"fit",runner)
spp_worker_gene(ctx1,fw$gene[1],"fit",runner)
stopifnot(calls==1L,readRDS(sre_attempt_file(ctx1,fw$gene[1]))$attempt==1L)
# A cohort change, nonconvergence, corrupt output or stale result cannot certify.
bad <- out(spe_prior(ctx1$priors,"Brain"),1L); bad$Liver <- bad$Brain
saveRDS(bad,sre_result(ctx1,"G003"))
stopifnot(spp_inspect(ctx1,"G003",force=TRUE)$status=="retry")
bad <- out(); bad$Brain$fit_slide_prior$converged <- FALSE
saveRDS(bad,sre_result(ctx,"G001"))
fails(spp_require_complete(root,d),"missing or changed")
ix <- spp_index(ctx); stopifnot(is.null(ix$entries$G001))
stopifnot(spp_inspect(ctx,"G001",force=TRUE)$status=="retry")
writeLines("corrupt",sre_result(ctx,"G002"))
stopifnot(spp_inspect(ctx,"G002",force=TRUE)$status=="retry")
# Successful validation never spends or resets fitting attempts.
attempt <- list(gene="G004",attempt=3L,frozen=ctx$frozen,error="Interrupted before receipt")
em_atomic_write(attempt,sre_attempt_file(ctx,"G004"))
spp_worker_gene(ctx,"G004","validate")
stopifnot(identical(readRDS(sre_attempt_file(ctx,"G004")),attempt))
# Unknown workers, duplicate ownership and launcher failures fail closed.
writeLines("999_1|999|slide_parallel_fit",queue)
fails(spp_queue(ctx,queue),"Unregistered")
writeLines(character(),queue)
orphan <- spp_plan(ctx1,1L,queue)
writeLines("777",file.path(orphan$batch,"array_job_id.txt"))
fails(spp_queue(ctx1,queue),"never entered R")
cat("PASS: parallel adoption/resume, identical EM counts/weights, no coordinator posterior reads, one-gene fitting, frozen cohort, checksums, exclusions, warm starts and failure guards.\n")
