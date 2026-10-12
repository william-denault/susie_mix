# Rscript --vanilla script/analysis/tests/test_em_gwas_coloc.R
# Optional local testing library: Sys.setenv(COLOC_TEST_LIBRARY="tmp/coloc-test-library")
if (nzchar(Sys.getenv("COLOC_TEST_LIBRARY"))) .libPaths(c(Sys.getenv("COLOC_TEST_LIBRARY"), .libPaths()))
source("script/analysis/em_gwas_coloc_utils.R")
cg_need(c("coloc", "data.table", "jsonlite", "digest", "susieR"))
expect_error <- function(expr, pattern) {
  err <- tryCatch(force(expr), error = identity)
  stopifnot(inherits(err, "error"), grepl(pattern, conditionMessage(err)))
}
p <- 120L
variants <- data.frame(chromosome = "1", position = seq_len(p)*1000, ref = "A", alt = "G", build = "GRCh38")
key <- cg_key(variants$chromosome, variants$position, variants$ref, variants$alt)
variants$snp <- key
make_qtl <- function(bf, components = seq_len(nrow(bf))) {
  colnames(bf) <- key
  list(lbf = bf, alpha = exp(bf - apply(bf, 1, cg_lse)), variants = variants, components = components)
}
em_bf <- matrix(-2, 2, p); em_bf[1, 10] <- 25; em_bf[2, 70] <- 23
vanilla_bf <- em_bf; vanilla_bf[1, 10] <- -2; vanilla_bf[1, 40] <- 25
bundle <- list(schema_version = 1L, metadata = list(iteration = 18L),
               em = make_qtl(em_bf, c(7, 2)), vanilla = make_qtl(vanilla_bf, c(3, 8)))
gw <- make_qtl(em_bf[1, , drop = FALSE], "GWAS_1")
gw$variants$key <- key; gw$analysis <- "synthetic"
priors <- data.frame(p1 = 1e-4, p2 = 1e-4, p12 = c(1e-6, 5e-6, 1e-5))
res <- cg_pair(bundle, gw, priors)
stopifnot(nrow(res$pairs) == 12L, all(res$pairs$status == "ok"),
          all(abs(rowSums(res$pairs[, paste0("H",0:4), with=FALSE]) - 1) < 1e-10))
delta <- cg_changes(res$pairs)
stopifnot(all(delta$classification == "gain"), all(delta$em_component == 7L),
          all(delta$vanilla_top_hypothesis == "H3"), all(delta$em_top_hypothesis == "H4"))
swapped <- bundle; swapped$em <- bundle$vanilla; swapped$vanilla <- bundle$em
loss <- cg_changes(cg_pair(swapped, gw, priors)$pairs)
stopifnot(all(loss$classification == "loss"), max(abs(delta$delta_H4 + loss$delta_H4)) < 1e-12)

# Enumerate the five hypotheses independently of coloc, including absolute BF scale.
manual <- function(a, b, p1, p2, p12) {
  a <- exp(a); b <- exp(b)
  joint <- outer(a, b); diag(joint) <- 0
  x <- c(1, p1*sum(a), p2*sum(b), p1*p2*sum(joint), p12*sum(a*b))
  x/sum(x)
}
for (k in seq_len(nrow(priors))) {
  row <- res$pairs[method == "em" & qtl_component == 7 & p12 == priors$p12[k]]
  expected <- manual(em_bf[1,], gw$lbf[1,], priors$p1[k], priors$p2[k], priors$p12[k])
  stopifnot(max(abs(unlist(row[, paste0("H",0:4), with=FALSE]) - expected)) < 1e-10)
}
weak <- bundle; weak$em <- make_qtl(matrix(-10, 1, p)); weak$vanilla <- weak$em
weak_gw <- gw; weak_gw$lbf[] <- -10; weak_gw$alpha[] <- 1/p
weak_pairs <- cg_pair(weak, weak_gw, priors)$pairs
stopifnot(all(weak_pairs$H0 > .99))

# Reordering GWAS alleles/columns cannot change a BF-based answer.
reordered <- gw
reordered$lbf <- gw$lbf[, p:1, drop=FALSE]; reordered$alpha <- gw$alpha[, p:1, drop=FALSE]
stopifnot(isTRUE(all.equal(cg_pair(bundle, reordered, priors)$pairs$H4, res$pairs$H4)))
# Losing the causal SNP must fail the PRE-intersection posterior coverage check.
missing <- gw; missing$lbf <- gw$lbf[, -10, drop=FALSE]; missing$alpha <- gw$alpha[, -10, drop=FALSE]
bad <- cg_pair(bundle, missing, priors, min_snps = 100)$pairs
stopifnot(any(bad$status == "insufficient_posterior_coverage"),
          all(cg_changes(bad)$classification == "unassessed"))
no_vanilla <- bundle; no_vanilla$vanilla$lbf <- vanilla_bf[FALSE,,drop=FALSE]
colnames(no_vanilla$vanilla$lbf) <- key
no_vanilla$vanilla$alpha <- no_vanilla$vanilla$lbf; no_vanilla$vanilla$components <- integer()
stopifnot(all(cg_changes(cg_pair(no_vanilla, gw, priors)$pairs)$classification == "unassessed"))

g <- data.frame(chromosome = "1", position = variants$position, effect_allele = "G", other_allele = "A",
                beta = .1, se = .05, p = .05, n = 1000, eaf = .3, info = 1, log_bf = NA_real_)
g$key <- key
abf <- cg_abf(g, list(columns = list(), trait_type = "cc"))
r <- .2^2/(.2^2+.05^2)
stopifnot(max(abs(abf$lbf - .5*(log(1-r)+r*4))) < 1e-12)
# Quantitative beta units: sdY estimation must scale with beta and SE.
quant <- cg_abf(g, list(columns=list(),trait_type="quant"))
scaled <- g; scaled$beta <- 10*g$beta; scaled$se <- 10*g$se
stopifnot(max(abs(cg_abf(scaled,list(columns=list(),trait_type="quant"))$lbf - quant$lbf)) < 1e-10)
na_g <- g; na_g$eaf <- NA_real_
expect_error(cg_abf(na_g, list(columns=list(),trait_type="quant")), "effect prior SD")
dupe <- rbind(g, g[1,]); dupe$other_allele[2] <- "C" # G/C palindrome
clean <- cg_clean_gwas(dupe)
stopifnot(!key[1] %in% clean$data$key, !variants$position[2] %in% clean$data$position)

# Actual SuSiE RSS with independent, simulated reference genotypes, two signals.
set.seed(312)
X <- matrix(rbinom(1000*p, 2, .3), 1000, p)
y <- 1.2*X[,10] + .9*X[,70] + rnorm(1000)
zs <- vapply(seq_len(p), function(j) coef(summary(lm(y ~ X[,j])))[2,3], 0)
g$beta <- zs*.05
fake_ld <- function(g, dir) list(gwas=g, R=cor(X), n_reference=nrow(X))
dir <- tempfile("rss-test-"); dir.create(dir)
fitted <- cg_susie(g, list(primary_eligible=TRUE), fake_ld, dir, L=3, max_s=1)
stopifnot(nrow(fitted$lbf) >= 2, file.exists(file.path(dir,"gwas_fit.rds")),
          all(is.finite(fitted$lbf)), all(abs(rowSums(fitted$alpha)-1)<1e-8))
bad_n <- g; bad_n$n[1] <- 100
expect_error(cg_susie(bad_n,list(primary_eligible=TRUE),fake_ld,dir), "varies")

# Complete on-disk run, checksum validation, all traits/cases/modes audited.
job <- tempfile("coloc-job-"); dir.create(job); dir.create(file.path(job,"gwas")); dir.create(file.path(job,"bundles"))
saveRDS(list(iteration=18L), file.path(job,"job_metadata.rds"))
bundle_path <- file.path(job,"bundles/test.rds"); saveRDS(bundle,bundle_path)
write.csv(data.frame(case_id="case_1",gene="GENE",tissue="TISSUE",bundle="bundles/test.rds",md5=unname(tools::md5sum(bundle_path))),file.path(job,"cases.csv"),row.names=FALSE)
write.csv(data.frame(case_id="case_1",chromosome="1",start=1,end=200000),file.path(job,"regions.csv"),row.names=FALSE)
source <- list(trait_id="test",trait_type="cc",ancestry="EUR",primary_eligible=TRUE,columns=list(beta="beta"))
jsonlite::write_json(list(sources=list(source,within(source,trait_id<-"absent"))),file.path(job,"sources.json"),auto_unbox=TRUE)
gf <- file.path(job,"gwas/test.tsv.gz");data.table::fwrite(g,gf,sep="\t")
jsonlite::write_json(list(source=source,output_sha256=digest::digest(file=gf,algo="sha256")),file.path(job,"gwas/test.json"),auto_unbox=TRUE)
output <- run_em_gwas_coloc(job, modes=c("abf","susie"))
audit <- read.csv(file.path(output,"audit.csv"))
stopifnot(nrow(audit)==4L, any(audit$status=="pending_ld"),sum(audit$status=="missing_or_invalid_gwas")==2L,
          any(audit$status=="completed"),file.exists(file.path(output,"paired_changes.csv")))

# Select a completed screening run; ignore a newer unfinished run. Freeze inputs.
project <- tempfile("prepare-test-")
dir.create(file.path(project,"script/analysis"),recursive=TRUE)
file.copy("script/analysis/coloc_gwas_sources.json",file.path(project,"script/analysis"))
iteration_dir <- file.path(project,"results_em/iteration_018")
screen <- file.path(iteration_dir,"coloc_screen")
for (name in c("older","newer","unfinished")) {
  path <- file.path(screen,name);dir.create(file.path(path,"bf"),recursive=TRUE)
  saveRDS(bundle,file.path(path,"bf/test.rds"))
  write.csv(data.frame(gene="GENE",tissue="TISSUE",bundle="bf/test.rds",bf_status="exported"),file.path(path,"bf_manifest.csv"),row.names=FALSE)
  saveRDS(list(iteration=18L,completed_at=if(name=="unfinished") NULL else
    if(name=="newer") "2026-10-11 20:00:00 UTC" else "2026-10-11 19:00:00 UTC"),file.path(path,"screen_metadata.rds"))
}
# The iteration selector itself has separate integration tests in test_em_coloc_screen.R.
em_analysis_context <- function(project_dir,iteration) list(iteration=18L,iteration_dir=iteration_dir)
prepared <- prepare_em_gwas_coloc(project)
stopifnot(grepl("newer$",readRDS(file.path(prepared,"job_metadata.rds"))$screen_dir),
          file.exists(file.path(prepared,"bundles/case_000001.rds")),
          read.csv(file.path(prepared,"regions.csv"))$start==1)
cat("PASS: exact H0-H4 enumeration, gains and losses, nulls, allele ordering,\n",
    "coverage loss, no-CS comparator, effect units, duplicate QC, real SuSiE RSS and disk workflow.\n",sep="")
