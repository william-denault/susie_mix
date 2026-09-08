# Rscript --vanilla script/simulation/tests/test_simulation.R (from project root)
source("script/simulation/config.R")
source("script/simulation/simulation_utils.R")
config <- simulation_config(getwd())
config$pilot <- TRUE
helpers <- simulation_helpers(getwd())
expect_error <- function(expr) stopifnot(inherits(tryCatch({force(expr); NULL},
                                                  error = identity), "error"))

# All requested architectures/counts/PVEs, minimum replication count, N and L.
scenarios <- simulation_scenarios()
stopifnot(nrow(scenarios) == 90L, nrow(scenarios) * config$replications == 18000L,
          config$L == 10L, config$n == 200L,
          all(table(scenarios$architecture) == c(15, 12, 12, 9, 15, 15, 12)[
            match(names(table(scenarios$architecture)), c("additive", "additive_dominant",
              "additive_recessive", "additive_recessive_dominant", "dominant", "recessive",
              "recessive_dominant"))]))
bad <- config
bad$pilot <- FALSE
bad$replications <- 199L
expect_error(simulation_validate(bad))
bad <- config
bad$L <- 5L
expect_error(simulation_validate(bad))
bad <- config
bad$pilot <- FALSE
bad$raw_file <- "any.raw"
expect_error(simulation_validate(bad))

# Guard against drifting away from the actual workhorse settings.
env <- new.env()
program <- parse("script/scan_tissue_attempt/workhorse.R")
for (expr in program) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name("run_susie_gene"))) eval(expr, env)
}
for (nm in c("cis_window", "min_maf_plink", "min_maf", "hwe_thresh", "min_n_rec",
             "L", "standardize", "estimate_prior_method", "min_abs_corr"))
  stopifnot(isTRUE(all.equal(config[[nm]], eval(formals(env$run_susie_gene)[[nm]]))))

for (pkg in c("data.table", "matrixStats", "susieR"))
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Install ", pkg, " to run tests.")
set.seed(112)
X <- sapply(seq_len(40), function(i) rbinom(400, 2, runif(1, 0.2, 0.45)))
colnames(X) <- paste0("snp", seq_len(ncol(X)))
rownames(X) <- paste0("donor", seq_len(nrow(X)))
X <- cbind(X, flip = rbinom(400, 2, 0.8), missing = c(NA, rep(1, 399)),
           invariant = 0, low_maf = c(rep(1, 10), rep(0, 390)),
           hwe_failure = rep(c(0, 2), 200))
qc <- simulation_global_qc(X, config, helpers, "fixture")
stopifnot(!any(c("missing", "invariant", "low_maf", "hwe_failure") %in% colnames(qc$X)),
          all(colMeans(qc$X) <= 1), all(colMeans(qc$X) > 0.1),
          identical(unname(qc$X[, "flip"]), as.integer(2 - X[, "flip"])))
design <- simulation_design(qc$X, config, helpers, rownames(qc$X)[1:200])
expected <- helpers$recode_snp_matrix(qc$X[1:200, ])
expected_mix <- do.call(cbind, expected[c("additive", "recessive", "dominant")])
storage.mode(expected_mix) <- "double"
keep <- colSums(expected_mix) >= config$min_n_rec & matrixStats::colSds(expected_mix) > 0
stopifnot(identical(unname(design$mix), unname(expected_mix[, keep, drop = FALSE])),
          nrow(design$add) == 200, setequal(colnames(design$add), unique(design$map$snp)))

# A rare recessive coding must disappear in the analysis sample while its additive
# and dominant codings survive; no extra allele reorientation occurs in this stage.
Z <- matrix(rep(c(rep(0, 120), rep(1, 76), rep(2, 4)), 2), ncol = 1,
             dimnames = list(paste0("z", 1:400), "rare_rec"))
z <- simulation_design(Z, config, helpers, rownames(Z)[1:200])
stopifnot(identical(sort(z$map$coding), c("additive", "dominant")))

for (i in seq_len(nrow(scenarios))) {
  s <- scenarios[i, ]
  set.seed(simulation_seed(config$seed, s$scenario_id, 1))
  idx <- simulation_causals(design, s, config)
  stopifnot(length(idx) == s$k, !anyDuplicated(design$map$snp[idx]),
            setequal(design$map$coding[idx], strsplit(s$architecture, "_")[[1]]))
  ph <- simulation_phenotype(design, idx, s$pve)
  stopifnot(abs(var(ph$genetic_value) - s$pve) < 1e-12,
            abs(ph$variance["noise_variance_target"] - (1 - s$pve)) < 1e-12,
            max(abs(drop(scale(design$mix[, idx, drop = FALSE], scale = FALSE) %*%
              ph$truth$beta) - ph$genetic_value)) < 1e-12)
}
seeds <- outer(1:90, 1:200, function(s, r) simulation_seed(config$seed, s, r))
stopifnot(!anyDuplicated(as.vector(seeds)))

# Analytical grouped PIP: A has two mutually exclusive codings per component.
map <- data.frame(predictor = c("A_a", "A_r", "B_a"),
                   snp = c("A", "A", "B"), coding = c("additive", "recessive", "additive"))
fit <- structure(list(alpha = rbind(c(.3, .2, .5), c(.1, .1, .8), c(.5, .4, .1)),
                       V = c(1, 1, 0)), class = "susie")
pip <- simulation_snp_pip(fit, map, c("B", "A"))
stopifnot(isTRUE(all.equal(unname(pip), c(.9, .6))), identical(names(pip), c("B", "A")))
fit$V[] <- 0
stopifnot(all(simulation_snp_pip(fit, map, c("A", "B")) == 0))

# Boundary PIPs are included once in calibration; no-CS coverage is undefined.
cal <- simulation_calibration(c(0, 0.05, 1), c(FALSE, TRUE, TRUE), config$calibration_breaks)
stopifnot(sum(cal$n) == 3, sum(cal$n_causal) == 2, cal$n[1] == 2, tail(cal$n, 1) == 1)
small <- list(add = design$add[, 1:3])
map <- data.frame(predictor = colnames(small$add), snp = colnames(small$add), coding = "additive")
fake <- structure(list(alpha = matrix(c(.8, .1, .1), 1), V = 1, sets = list(cs = NULL)), class = "susie")
ph <- list(truth = map[1, , drop = FALSE])
ev <- simulation_evaluate(fake, small$add, map, small, ph, config)
stopifnot(ev$metrics$n_cs == 0, is.na(ev$metrics$cs_coverage),
          ev$metrics$cs_power == 0, ev$metrics$power == 0,
          ev$metrics$fdp == 0, ev$curve$tp[1] == 1, ev$curve$fp[1] == 2,
          tail(ev$curve$discoveries, 1) == 0)

# Real paired model fit: no truth-dependent L and the grouped additive PIPs agree
# with the package's predictor PIPs exactly.
config$save_fits <- TRUE
set.seed(921)
s <- subset(scenarios, architecture == "additive_recessive_dominant" & k == 3 & pve == .2)
idx <- simulation_causals(design, s, config)
ph <- simulation_phenotype(design, idx, s$pve)
out <- simulation_fit(list(design = design, phenotype = ph), config)
for (method in names(out)) {
  if (!is.null(out[[method]]$error)) stop(out[[method]]$error)
  stopifnot(out[[method]]$converged, nrow(out[[method]]$fit$alpha) == 10L,
            nrow(out[[method]]$snp_pip) == ncol(design$add))
}
stopifnot(max(abs(out$susie$snp_pip$pip - out$susie$predictor_pip$pip)) < 1e-12)
cat("Simulation tests passed: QC, 90 scenarios, causal truth, PVE, grouped PIPs, metrics, L=10 fits.\n")
