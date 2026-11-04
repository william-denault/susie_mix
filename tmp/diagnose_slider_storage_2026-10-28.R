.libPaths(c("tmp/slider_prior_impl/library", "C:/Users/willi/AppData/Local/R/win-library/4.6", .libPaths()))
cat("Package version:", as.character(packageVersion("susieRSlidePrior")), "\n")
fn <- deparse(body(susieRSlidePrior::susie))
cat(paste(fn[grepl("match.call|fit\\$call", fn)], collapse = "\n"), "\n")
set.seed(821)
X <- matrix(rbinom(200 * 500, 2, 0.3), 200, 500)
colnames(X) <- paste0("snp", seq_len(ncol(X)))
y <- 0.5 * X[, 1] + rnorm(nrow(X))
args <- list(X = X, y = y, L = 3L, delta_grid = seq(-1, 1, length.out = 17L),
             delta_prior = rep(1/17, 17), min_obs = 5L, standardize = FALSE,
             estimate_prior_method = "EM", max_iter = 1000L, tol = 1e-5,
             coverage = NULL, verbose = FALSE)
compressed_bytes <- function(x) length(memCompress(serialize(x, NULL), "gzip"))
depth <- function(x) {
  n <- 0L
  while (is.list(x) && is.call(x$call) && "model_init" %in% names(as.list(x$call))) {
    x <- x$call$model_init
    n <- n + 1L
    if (n > 30L) stop("Unexpected recursion")
  }
  n
}
previous <- NULL
rows <- list()
for (iteration in 0:4) {
  a <- args
  if (!is.null(previous)) a$model_init <- previous
  set.seed(100)
  fit <- do.call(susieRSlidePrior::susie, a)
  compact <- fit
  compact$call <- NULL
  rows[[iteration+1L]] <- data.frame(iteration, converged = fit$converged,
    fit_bytes = as.numeric(object.size(fit)),
    call_bytes = as.numeric(object.size(fit$call)),
    compressed_bytes = compressed_bytes(fit),
    without_call_compressed_bytes = compressed_bytes(compact),
    embedded_previous_fits = depth(fit),
    call_contains_genotypes = is.matrix(fit$call$X))
  previous <- fit
}
print(do.call(rbind, rows), row.names = FALSE)
cat("Largest current-fit fields (excluding call):\n")
print(head(sort(vapply(compact, function(x) as.numeric(object.size(x)), numeric(1)), decreasing = TRUE), 8))
a <- args; a$model_init <- previous
set.seed(100); from_full <- do.call(susieRSlidePrior::susie, a)
a$model_init$call <- NULL
set.seed(100); from_compact <- do.call(susieRSlidePrior::susie, a)
from_full$call <- from_compact$call <- NULL
cat("Next fit identical with/without warm-start call:", identical(from_full, from_compact), "\n")
stopifnot(identical(from_full, from_compact))
