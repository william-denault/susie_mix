# Hand-calculated synthetic posteriors; exercises the complete plotting script.
fixture <- tempfile("posterior_cdf_", tmpdir = "tmp")
dir.create(fixture)
root <- file.path(fixture, "results_slide_prior_em")
for (i in c("iteration_009", "iteration_024", "iteration_025"))
  dir.create(file.path(root, i, "results"), recursive = TRUE)
alpha <- rbind(c(.2, .3, .5), c(.1, .6, .3))
g <- c(1, -1, 0)  # Deliberately unsorted.
joint <- array(0, c(2, 3, 3))
for (ell in 1:2) {
  joint[ell, 1, ] <- alpha[ell, 1] * c(.2, .5, .3)
  joint[ell, 2, ] <- alpha[ell, 2] * c(.6, .1, .3)
  joint[ell, 3, ] <- alpha[ell, 3] * c(0, 0, 1)
}
joint[2, 1, ] <- .1 * c(.1, .7, .2)
fit <- list(alpha = alpha, alpha_delta = joint, delta_grid = g,
  delta_forced = c(FALSE, FALSE, TRUE), converged = TRUE,
  sets = list(cs = list(L2 = c(1L, 2L), forced = 3L, singleton = 1L),
              cs_index = c(2L, 1L, 1L)))
out <- list(Brain = list(fit_slide_prior = fit), Liver = list(fit_slide_prior = fit),
            Empty = list())
saveRDS(out, file.path(root, "iteration_024/results/gene1.rds"))
saveRDS(out[1], file.path(root, "iteration_024/results/filtered.rds"))
writeLines("should not be read", file.path(root, "iteration_009/results/old.rds"))
old_env <- Sys.getenv("SUSIE_MIX_PROJECT_DIR", unset = NA_character_)
Sys.setenv(SUSIE_MIX_PROJECT_DIR = normalizePath(fixture, winslash = "/"))
e <- new.env()
source("plot_posterior_cdf_last_iteration.R", local = e)
if (is.na(old_env)) Sys.unsetenv("SUSIE_MIX_PROJECT_DIR") else
  Sys.setenv(SUSIE_MIX_PROJECT_DIR = old_env)
x <- e$posterior_cdf
close_to <- function(a, b) isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = 1e-12))
stopifnot(x$iteration == "iteration_024", nrow(x$signals) == 4,
          x$filtered_files == 1, x$forced_leads_excluded == 2,
          identical(x$grid, c(-1, 0, 1)),
          identical(x$signals$component, c(2L, 1L, 2L, 1L)),
          close_to(x$lead_cdf[1, ], c(.1, .4, 1)),
          close_to(x$cs_pmf[1, ], c(.13, .20, .37) / .7),
          close_to(x$cs_cdf[1, ], c(.13, .33, .7) / .7),
          close_to(x$lead_cdf[2, ], c(.5, .8, 1)),
          close_to(x$cs_cdf[2, ], x$lead_cdf[2, ]),
          close_to(x$tissue_cdf$lead_cdf[x$tissue_cdf$pooled], c(.3, .6, 1)),
          close_to(x$tissue_cdf$cs_cdf[x$tissue_cdf$pooled],
                   (c(.13, .33, .7) / .7 + c(.5, .8, 1)) / 2),
          all(x$lead_cdf >= 0 & x$lead_cdf <= 1),
          all(x$cs_cdf >= 0 & x$cs_cdf <= 1))
y <- e$read_slider_posteriors(root, "iteration_024", FALSE, min_entries = 0)
stopifnot(nrow(y$signals) == 9,
          close_to(y$lead_cdf[y$signals$lead_delta_forced, ][1, ], c(0, 1, 1)))
files <- list.files(e$output_dir, full.names = TRUE)
stopifnot(length(files) == 7L, all(file.info(files)$size > 0))
stopifnot(identical(readRDS(file.path(e$output_dir, "posterior_cdf.rds")), x))
bad <- out
bad$Brain$fit_slide_prior$alpha_delta <- NULL
dir.create(file.path(root, "iteration_030/results"), recursive = TRUE)
saveRDS(bad, file.path(root, "iteration_030/results/missing.rds"))
err <- tryCatch(e$read_slider_posteriors(root), error = conditionMessage)
stopifnot(is.character(err), grepl("alpha_delta", err), grepl("missing.rds", err))
bad <- out
bad$Brain$fit_slide_prior$alpha_delta[2, 1, 1] <- 1
saveRDS(bad, file.path(root, "iteration_030/results/missing.rds"))
err <- tryCatch(e$read_slider_posteriors(root), error = conditionMessage)
stopifnot(is.character(err), grepl("does not agree with alpha", err))
cat("PASS: analytic CDFs, CS conditioning, component mapping, singleton CS, grid sorting,\n",
    "forced filtering, equal CS weighting, latest nonempty iteration, missing/invalid posterior,\n",
    "saved RDS/CSV and both PNG/PDF figures.\n")
cat("PREVIEW:", normalizePath(file.path(e$output_dir, "lead_posterior_cdf.png"), winslash = "/"), "\n")
