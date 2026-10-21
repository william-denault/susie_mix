source("script/analysis/summarize_slide_prior_delta.R")
set.seed(42)
make_fit <- function(grid = c(-1, -.5, 0, .5, 1)) {
  alpha <- matrix(c(.1, .2, .3, .25, .15,
                    .2, .2, .2, .2, .2,
                    .05, .1, .5, .2, .15), 3, 5, byrow = TRUE)
  joint <- array(0, c(3, 5, length(grid)))
  forced <- c(FALSE, TRUE, FALSE, FALSE, TRUE)
  for (ell in 1:3) for (j in 1:5) {
    w <- runif(length(grid))
    if (forced[j]) w <- as.numeric(grid == 0)
    joint[ell, j, ] <- alpha[ell, j] * w / sum(w)
  }
  list(alpha = alpha, alpha_delta = joint, delta_grid = grid,
       sets = list(cs = list(L3 = c(1L, 3L), L1 = 2L), cs_index = c(3L, 1L)),
       delta_forced = forced, min_obs = 5,
       genotype_counts = matrix(1:15, 5, 3, dimnames = list(NULL, c("0", "1", "2"))))
}
reference <- function(file) {
  out <- readRDS(file)
  rows <- list()
  summarize_delta <- function(weights, grid, prefix) {
    prob <- weights / sum(weights)
    ans <- c(delta_mean = sum(prob * grid),
             p_recessive_like = sum(prob[grid < 0]),
             p_recessive = sum(prob[grid == -1]),
             p_additive = sum(prob[grid == 0]),
             p_dominant_like = sum(prob[grid > 0]))
    setNames(ans, paste0(prefix, names(ans)))
  }
  if (length(out) <= 2L) return(data.frame())
  for (o in seq_along(out)) {
    fit <- out[[o]]$fit_slide_prior
    if (!length(fit$sets$cs)) next
    tissue <- if (is.null(names(out))) as.character(o) else names(out)[o]
    for (s in seq_along(fit$sets$cs)) {
      cs <- fit$sets$cs[[s]]
      ell <- fit$sets$cs_index[s]
      alpha_cs <- fit$alpha[ell, cs]
      lead <- cs[which.max(alpha_cs)]
      cs_weights <- apply(fit$alpha_delta[ell, cs, , drop = FALSE], 3, sum)
      rows[[length(rows) + 1L]] <- data.frame(
        file = basename(file), tissue = tissue, cs = names(fit$sets$cs)[s],
        component = ell, n_snps = length(cs), cs_coverage = sum(alpha_cs),
        lead_index = lead, lead_alpha = unname(fit$alpha[ell, lead]),
        lead_delta_forced = unname(fit$delta_forced[lead]),
        lead_n0 = unname(fit$genotype_counts[lead, "0"]),
        lead_n1 = unname(fit$genotype_counts[lead, "1"]),
        lead_n2 = unname(fit$genotype_counts[lead, "2"]), min_obs = fit$min_obs,
        cs_forced_mass = sum(alpha_cs * fit$delta_forced[cs]) / sum(alpha_cs),
        as.list(summarize_delta(cs_weights, fit$delta_grid, "cs_")),
        as.list(summarize_delta(fit$alpha_delta[ell, lead, ], fit$delta_grid, "lead_")),
        row.names = NULL)
    }
  }
  if (length(rows)) do.call(rbind, rows) else data.frame()
}
fixture <- tempfile("slide_delta_", tmpdir = "tmp")
dir.create(fixture)
f <- make_fit()
single <- make_fit(0)
valid <- list(brain = list(fit_slide_prior = f), liver = list(fit_slide_prior = f), empty = list())
saveRDS(valid, file.path(fixture, "valid.rds"))
saveRDS(list(list(fit_slide_prior = single), list(), list()), file.path(fixture, "single.rds"))
saveRDS(valid[1:2], file.path(fixture, "filtered.rds"))
saveRDS(list(list(), list(), list()), file.path(fixture, "no_cs.rds"))
writeLines("bad RDS", file.path(fixture, "broken.rds"))
bad <- valid
bad$liver$fit_slide_prior$sets$cs_index <- c(99L, 1L)
saveRDS(bad, file.path(fixture, "bad_index.rds"))
files <- list.files(fixture, full.names = TRUE, pattern = "rds$")
serial <- suppressWarnings(summarize_slide_prior_delta(fixture, batch_size = 2,
                          batch_dir = file.path(fixture, "batches"), verbose = FALSE))
psock <- suppressWarnings(summarize_slide_prior_delta(fixture, workers = 2,
                          batch_size = 3, verbose = FALSE))
stopifnot(identical(serial, psock), nrow(serial$results) == 6L,
          nrow(serial$errors) == 2L,
          serial$counts["filtered"] == 1L, serial$counts["ok"] == 3L)
expected <- .slide_delta_bind(lapply(files[!basename(files) %in% c("broken.rds", "bad_index.rds")], reference))
stopifnot(isTRUE(all.equal(serial$results, expected, tolerance = 1e-14)))
stopifnot(all(serial$results$lead_index == rep(c(3L, 2L), 3L)),
          all(serial$results$lead_p_additive[serial$results$lead_delta_forced] == 1))
parts <- lapply(list.files(file.path(fixture, "batches"), full.names = TRUE), readRDS)
stopifnot(identical(.slide_delta_bind(lapply(parts, `[[`, "results")), serial$results))
unfiltered <- summarize_slide_prior_delta(fixture, files = file.path(fixture, "filtered.rds"),
                                         min_entries = 0, verbose = FALSE)
stopifnot(nrow(unfiltered$results) == 4L)
cat("PASS: all columns agree with original; singleton CS/grid; reordered components;\n",
    "forced SNPs; serial/PSOCK equality; filtering; file errors; saved batches.\n")
