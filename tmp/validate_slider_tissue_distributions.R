source("script/analysis/slider_tissue_distributions.R")
make_data <- function(n = 300, shape_difference = FALSE, missing = FALSE) {
  z <- rnorm(n)
  ans <- do.call(rbind, lapply(seq_len(3), function(t) {
    scale <- if (shape_difference && t == 2) 0.9 else 0.4
    data.frame(file = paste0("g", seq_len(n)), tissue = LETTERS[t],
               lead_delta_mean = scale * tanh(z + 0.4 * rnorm(n)),
               lead_delta_forced = FALSE)
  }))
  if (missing) ans <- ans[runif(nrow(ans)) > 0.15, ]
  ans
}
set.seed(10)
d <- make_data(300, TRUE, TRUE)
a <- slider_tissue_distributions(d, B = 499, n_grid = 201, verbose = FALSE)
stopifnot(a$omnibus$p_value < 0.02,
          all(a$pairwise$W1 >= 0),
          all(diag(a$coassignment[["2"]]) == 1),
          isTRUE(all.equal(a$coassignment[["2"]], t(a$coassignment[["2"]]))))
for (tissue in rownames(a$cdf)) {
  x <- d$lead_delta_mean[d$tissue == tissue]
  stopifnot(max(abs(a$cdf[tissue, ] - ecdf(x)(a$grid))) < 1e-12)
}
# Endpoints and exact knots, including unequal numbers of CS per gene.
edge <- data.frame(file = rep(c("a", "b", "c", "d", "e"), 2),
                   tissue = rep(c("A", "B"), each = 5),
                   lead_delta_mean = rep(c(-1, -0.5, 0, 0.5, 1), 2),
                   lead_delta_forced = FALSE)
e <- suppressWarnings(slider_tissue_distributions(edge, B = 99,
                                                  n_grid = 21, verbose = FALSE))
stopifnot(e$omnibus$p_value == 1, all(e$pairwise$W1 == 0),
          all(abs(e$cdf[1, ] - ecdf(edge$lead_delta_mean[1:5])(e$grid)) < 1e-12))
dup <- rbind(edge, edge[c(1, 1, 2), ])
ed <- suppressWarnings(slider_tissue_distributions(dup, B = 99,
                                                   n_grid = 21, verbose = FALSE))
stopifnot(max(abs(e$cdf - ed$cdf)) < 1e-12)
es <- suppressWarnings(slider_tissue_distributions(dup, B = 99,
                    n_grid = 21, weighting = "signal", verbose = FALSE))
stopifnot(es$pairwise$W1 > 0)
forced <- rbind(edge, transform(edge[1, ], file = "extra", lead_delta_forced = TRUE))
ef <- suppressWarnings(slider_tissue_distributions(forced, B = 99,
                                                    n_grid = 21, verbose = FALSE))
stopifnot(ef$settings$n_filtered == 1, identical(ef$cdf, e$cdf))
# Available versus common gene cohorts.
common <- slider_tissue_distributions(d, B = 99, cohort = "common", verbose = FALSE)
counts <- table(unique(d[c("file", "tissue")])$file)
stopifnot(all(common$sample_sizes$n_genes == sum(counts == 3)))
# Exact W1 for a selected pair versus the numerical grid calculation.
x <- d$lead_delta_mean[d$tissue == "A"]
y <- d$lead_delta_mean[d$tissue == "B"]
knots <- sort(unique(c(x, y)))
exact <- sum(abs(ecdf(x)(head(knots, -1)) - ecdf(y)(head(knots, -1))) * diff(knots))
approx <- a$pairwise$W1[a$pairwise$tissue1 == "A" & a$pairwise$tissue2 == "B"]
stopifnot(abs(exact - approx) < 2 * diff(a$grid)[1])
# A modest calibration diagnostic, not a proof: correlated tissues and missing
# records under equal marginals. The null is not independence of tissues.
null_p <- vapply(seq_len(60), function(i) {
  set.seed(1000 + i)
  d0 <- make_data(200, FALSE, TRUE)
  slider_tissue_distributions(d0, B = 199, n_grid = 101, ks = integer(0),
                              seed = 3000 + i, verbose = FALSE)$omnibus$p_value
}, numeric(1))
cat("Edge cases, weights, cohort handling, CDFs, and W1 checks passed.\n")
cat("Same-mean shape alternative p:", a$omnibus$p_value, "\n")
cat("Null rejections at 5%:", sum(null_p <= 0.05), "/", length(null_p), "\n")
