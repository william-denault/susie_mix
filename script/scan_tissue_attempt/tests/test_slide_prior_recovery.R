# Base-R tests: real metadata/checkpoint handling, no cluster or GTEx data needed.
repo <- normalizePath(".", winslash = "/")
for (f in c("em_utils.R", "slide_prior_em_utils.R", "slide_prior_recovery.R", "prepare_slide_prior_em_iteration.R"))
  source(file.path(repo, "script/scan_tissue_attempt", f))
fails <- function(expr, pattern) {
  e <- tryCatch({force(expr); NULL}, error = identity)
  stopifnot(inherits(e, "error"), grepl(pattern, conditionMessage(e)))
}
dir.create(file.path(repo, "tmp"), showWarnings = FALSE)
root <- tempfile("slide_recovery_", tmpdir = file.path(repo, "tmp")); dir.create(root)
scripts <- file.path(root, "script/scan_tissue_attempt"); dir.create(scripts, recursive = TRUE)
file.copy(file.path(repo, "script/scan_tissue_attempt/slide_prior_recovery.R"), scripts)
d <- file.path(root, "results_slide_prior_em/iteration_000")
dir.create(file.path(d, "results"), recursive = TRUE)
dir.create(file.path(d, "completed"))
genes <- sprintf("G%03d", 1:23)
manifest <- data.frame(chunk = rep(1:2, c(12, 11)), gene = genes)
w <- matrix(1/17, 2, 17, dimnames = list(NULL, spe_columns()))
priors <- data.frame(iteration = 0L, source_iteration = -1L, tissue = c("Brain", "Liver"), w,
                     n_fits = 0L, n_active = 0L, free_mass = 0, prior_max_change = 0,
                     m_step_q_gain = 0, source_elbo = NA_real_,
                     update_method = "slider_prior_active_counts_v1", created_at = "fixture")
write.csv(priors, file.path(d, "priors.csv"), row.names = FALSE)
write.csv(priors, file.path(root, "results_slide_prior_em/prior_history.csv"), row.names = FALSE)
write.csv(manifest, file.path(d, "manifest.csv"), row.names = FALSE)
write.csv(data.frame(delta = spe_grid()), file.path(d, "grid.csv"), row.names = FALSE)
saveRDS(c(spe_settings(), list(package_version = "fixture")), file.path(d, "settings.rds"))
saveRDS(list(chunk = 1L), file.path(d, "completed/chunk_001.done"))
ctx <- sre_context(root, d)
frozen_files <- c(file.path(d, names(ctx$frozen)), file.path(d, "completed/chunk_001.done"))
hashes <- tools::md5sum(frozen_files)
make_fit <- function(prior = rep(1/17, 17), converged = TRUE) {
  L <- 2L
  a <- matrix(rep(c(.8, .2), each = L), L, 2)
  joint <- array(0, c(L, 2L, 17L))
  for (l in 1:L) {joint[l, 1, ] <- .8 * prior; joint[l, 2, 9] <- .2}
  counts <- matrix(0, L, 17); counts[1, ] <- .8 * prior
  structure(list(alpha = a, alpha_delta = joint, delta_prior_counts = counts,
                 delta_forced = c(FALSE, TRUE), input_p = 2L, delta_grid = spe_grid(),
                 delta_prior = prior, V = c(1, 0), converged = converged, elbo = -100),
            class = c("susie_slide", "susie"))
}
make_out <- function(converged = TRUE) {
  x <- list(Brain = list(fit_slide_prior = make_fit(converged = converged)))
  attr(x, "em_iteration") <- 0L
  x
}
save_gene <- function(g, x) saveRDS(x, sre_result(ctx, g))
save_gene("G001", make_out())
error_out <- structure(list(gene = "G002", error = "PLINK failed: Cannot allocate memory"), em_iteration = 0L)
save_gene("G002", error_out)
save_gene("G003", structure(list(), em_iteration = 0L))
partial <- make_out(); attr(partial, "tissue_errors") <- list(Liver = "Allocation failure")
save_gene("G004", partial)
save_gene("G005", structure(list(error = "Unknown gene"), em_iteration = 0L))
save_gene("G006", make_out(FALSE))
bad <- make_out(); bad$Brain$fit_slide_prior$delta_prior <- c(1, rep(0, 16)); save_gene("G007", bad)
writeLines("corrupt", sre_result(ctx, "G008"))
good_hash <- tools::md5sum(sre_result(ctx, "G001"))
audit <- sre_audit(ctx)
stopifnot(sum(sre_ready(audit$status)) == 2L, audit$status[2] == "retry",
          audit$status[4] == "retry", audit$status[6] == "retry")
plan <- sre_load(ctx)$plan
stopifnot(identical(plan$gene, genes), identical(as.integer(table(plan$task)), c(10L, 10L, 3L)))
full_plan <- data.frame(task = ceiling(seq_len(18468)/10), gene = seq_len(18468))
stopifnot(length(unique(full_plan$task)) == 1847L, max(table(full_plan$task)) == 10L)
batch <- sre_plan_batch(ctx, 1L)
selected <- sre_task_genes(ctx, batch$path, 1L)
stopifnot(length(selected) == 8L, "G002" %in% selected, !"G001" %in% selected)
fails(sre_task_genes(ctx, batch$path, 2L), "Invalid batch")
fails(sre_audit(ctx, 16L), "genes_per_task")
fails(sre_audit(ctx, 15L), "settings differ")
fails(sre_verify(ctx), "unresolved")
fails(sre_require_complete(root, d), "unfinished")
fails(spe_prepare_iteration(root, "resume", package_version = "fixture"), "recovery-managed")

# Targeted pilots respect the frozen assignment, and launch failures cannot loop.
pilot <- sre_plan_batch(ctx, 5L, pilot_genes = c("G001", "G011", "G023"))
stopifnot(pilot$tasks == 2L,
          identical(sre_task_genes(ctx, pilot$path, 1L), "G011"),
          identical(sre_task_genes(ctx, pilot$path, 2L), "G023"),
          sre_plan_batch(ctx, 5L, pilot_genes = "G001")$status == "pilot_empty")
fails(sre_plan_batch(ctx, 5L, pilot_genes = "UNKNOWN"), "Pilot gene list")
writeLines(pilot$path, file.path(ctx$root, "last_batch.txt"))
stopifnot(sre_plan_batch(ctx, 5L, automatic = TRUE)$status == "blocked",
          sre_plan_batch(ctx, 5L)$status == "batch")
unlink(file.path(ctx$root, "last_batch.txt"))

# A reviewed error cannot conceal a malformed successful tissue posterior.
invalid <- partial; invalid$Brain$fit_slide_prior$alpha_delta[1, 1, 1] <- -1
save_gene("G004", invalid)
write.csv(data.frame(gene = "G004", result_md5 = unname(tools::md5sum(sre_result(ctx, "G004"))),
                     reason = "Fixture review"), file.path(ctx$root, "reviewed_exclusions.csv"), row.names = FALSE)
stopifnot(sre_inspect(ctx, "G004")$status == "retry")
save_gene("G004", partial)
unlink(file.path(ctx$root, "reviewed_exclusions.csv"))

calls <- character()
fake <- function(target_gene, ...) {calls <<- c(calls, target_gene); make_out()}
sre_run_gene(ctx, "G001", fake)
stopifnot(!length(calls), identical(good_hash, tools::md5sum(sre_result(ctx, "G001"))))
fails(sre_run_gene(ctx, "G002", function(...) {warning("system call failed: Cannot allocate memory"); stop("PLINK failed")}), "PLINK failed")
a <- readRDS(sre_attempt_file(ctx, "G002"))
stopifnot(a$attempt == 1L, grepl("Cannot allocate", a$error), identical(readRDS(sre_result(ctx, "G002")), error_out))
sre_run_gene(ctx, "G002", fake)
sre_run_gene(ctx, "G002", fake)
stopifnot(readRDS(sre_attempt_file(ctx, "G002"))$attempt == 2L, sum(calls == "G002") == 1L)
# Missing output gets an inspectable error, but never a completion status.
for (g in c("G020", "G021"))
  fails(sre_run_gene(ctx, g, function(...) stop("Cannot allocate memory")), "Cannot allocate")
stopifnot(file.exists(sre_result(ctx, "G020")), sre_inspect(ctx, "G020")$status == "retry",
          sre_plan_batch(ctx, 2L, automatic = TRUE)$status == "blocked",
          sre_plan_batch(ctx, 2L)$status == "batch")
for (g in c("G020", "G021")) sre_run_gene(ctx, g, fake)
# Corrupt or mismatched receipts must be rejected, including cached successes.
invisible(sre_refresh(ctx))
receipt <- readRDS(sre_attempt_file(ctx, "G002"))
writeLines("corrupt", sre_attempt_file(ctx, "G002"))
fails(sre_refresh(ctx), "Unreadable attempt")
wrong <- receipt; wrong$frozen <- "different"
saveRDS(wrong, sre_attempt_file(ctx, "G002"))
fails(sre_refresh(ctx), "different settings")
saveRDS(receipt, sre_attempt_file(ctx, "G002"))
fails(sre_run_gene(ctx, "G006", function(...) make_out(FALSE)), "Unconverged")
sre_run_gene(ctx, "G006", fake)
# A killed process consumes an attempt; an output saved before its receipt is recovered.
em_atomic_write(list(attempt = 1L, frozen = ctx$frozen, error = "Interrupted"), sre_attempt_file(ctx, "G009"))
save_gene("G009", make_out())
audit <- sre_refresh(ctx)
stopifnot(audit$status[audit$gene == "G009"] == "valid")
# Explicit approval is bound to the exact error output; it can resolve an exhausted retry.
for (i in 1:3) fails(sre_run_gene(ctx, "G005", function(...) stop("Unknown gene")), "Unknown gene")
stopifnot(sre_plan_batch(ctx, 298)$status == "blocked")
write.csv(data.frame(gene = "G005", result_md5 = unname(tools::md5sum(sre_result(ctx, "G005"))),
                     reason = "Fixture: reviewed missing annotation"),
          file.path(ctx$root, "reviewed_exclusions.csv"), row.names = FALSE)
audit <- sre_refresh(ctx)
stopifnot(audit$status[audit$gene == "G005"] == "reviewed_exclusion")
for (g in genes) if (!sre_ready(sre_inspect(ctx, g)$status)) sre_run_gene(ctx, g, fake)
stopifnot(identical(hashes, tools::md5sum(frozen_files)))
sre_verify(ctx)
# Scope is based on annotation, not failed-job text: it excludes missing genes,
# successful fits and exhausted failures alike, while retaining autosomal errors.
annotations <- data.frame(gene_name = genes, chromosome = rep("chr1", length(genes)))
annotations$chromosome[10:16] <- c("chrX", "chrY", "chrMT", "chrM", "23", "24", "26")
annotations <- rbind(annotations, data.frame(gene_name = "G001", chromosome = "X"))
fails(sre_exclude_sex_mt(ctx, annotations[1:3, ]), "Cannot assign chromosomes")
scope_file <- file.path(root, "results_slide_prior_em/chromosome_scope.rds")
stopifnot(!file.exists(scope_file))
unlink(sre_result(ctx, "G010"))
excluded_hash <- tools::md5sum(sre_result(ctx, "G011"))
save_gene("G012", structure(list(error = "PLINK failed on MT"), em_iteration = 0L))
save_gene("G017", structure(list(error = "Autosomal failure"), em_iteration = 0L))
writeLines("corrupt receipt on an excluded chromosome", sre_attempt_file(ctx, "G013"))
auto_attempt <- readRDS(sre_attempt_file(ctx, "G002"))$attempt
ctx <- sre_exclude_sex_mt(ctx, annotations)
scope_hash <- tools::md5sum(scope_file)
stopifnot(identical(ctx$scope$excluded, genes[10:16]),
          sre_inspect(ctx, "G001")$status == "valid",
          sre_inspect(ctx, "G010")$status == "scope_excluded",
          sre_inspect(ctx, "G012")$status == "scope_excluded",
          sre_inspect(ctx, "G017")$status == "retry",
          readRDS(sre_attempt_file(ctx, "G002"))$attempt == auto_attempt)
fails(sre_task_genes(ctx, batch$path, 1L), "Batch manifest changed")
fails(sre_require_complete(root, d), "unfinished")
scope_batch <- sre_plan_batch(ctx, 298L)
scope_work <- read.csv(file.path(scope_batch$path, "manifest.csv"))
stopifnot(!any(scope_work$gene %in% genes[10:16]))
for (g in genes[10:16]) sre_run_gene(ctx, g, function(...) stop("Excluded gene must not fit"))
stopifnot(!file.exists(sre_result(ctx, "G010")),
          identical(excluded_hash, tools::md5sum(sre_result(ctx, "G011"))))
ctx <- sre_exclude_sex_mt(ctx, annotations)
stopifnot(identical(scope_hash, tools::md5sum(scope_file)),
          identical(hashes, tools::md5sum(frozen_files)))
sre_run_gene(ctx, "G017", fake)
# Resume an activation interrupted after saving the policy but before its state.
active_state <- readRDS(file.path(ctx$root, "state.rds"))
interrupted_state <- active_state; interrupted_state$scope_md5 <- NULL
saveRDS(interrupted_state, file.path(ctx$root, "state.rds"))
fails(sre_context(root, d), "Chromosome scope changed")
ctx <- sre_exclude_sex_mt(sre_context(root, d, allow_scope_adoption = TRUE))
scope_backup <- file.path(root, "scope_backup.rds")
invisible(file.copy(scope_file, scope_backup))
unlink(scope_file)
fails(sre_context(root, d), "Chromosome scope changed")
invisible(file.copy(scope_backup, scope_file))
tampered <- readRDS(scope_file); tampered$genes$chromosome[1] <- "X"
saveRDS(tampered, scope_file)
fails(sre_context(root, d), "Chromosome scope changed")
invisible(file.copy(scope_backup, scope_file, overwrite = TRUE))
sre_verify(ctx)
sre_require_complete(root, d)
stopifnot(file.exists(file.path(ctx$root, "COMPLETE.rds")), all(sre_ready(sre_refresh(ctx)$status)))
# Refresh changes audit columns; recertification is required before a new update.
sre_verify(ctx)
next_iteration <- spe_prepare_iteration(root, "new", package_version = "fixture")
stopifnot(basename(next_iteration$iteration_dir) == "iteration_001")
updated <- read.csv(file.path(next_iteration$iteration_dir, "priors.csv"))
stopifnot(sum(updated$n_fits) == 14L,
          !any(sub("\\.rds$", "", basename(sre_pool_files(ctx))) %in% genes[10:16]))
ctx1 <- sre_context(root, next_iteration$iteration_dir)
sre_audit(ctx1)
fails(sre_exclude_sex_mt(ctx1, annotations), "before learned-prior updates")
sre_run_gene(ctx1, "G010", function(...) stop("Excluded genes must stay excluded in later iterations"))
stopifnot(sre_inspect(ctx1, "G010")$status == "scope_excluded")
sre_run_gene(ctx1, "G005", function(...) stop("Must inherit reviewed cohort exclusion"))
stopifnot(sre_inspect(ctx1, "G005")$status == "reviewed_exclusion")
fails(sre_run_gene(ctx1, "G001", function(...) list()), "cohort changed")
# The certificate cannot hide corruption of a previously validated fit.
writeLines("now corrupt", sre_result(ctx, "G001"))
fails(sre_require_complete(root, d), "unresolved")
# Neither re-audit nor a worker is allowed to silently bless manifest tampering.
write.csv(data.frame(task = 1, gene = "G001"), file.path(ctx$root, "manifest.csv"), row.names = FALSE)
fails(sre_audit(ctx), "manifest.*changed")
cat("PASS: recovery audit, batching, retries, frozen chromosome scope, excluded missing/successful/failed genes, pooling, warm-start cohort and completion gates.\n")
stopifnot(identical(dirname(normalizePath(root, winslash = "/")), paste0(repo, "/tmp")))
unlink(root, recursive = TRUE)
