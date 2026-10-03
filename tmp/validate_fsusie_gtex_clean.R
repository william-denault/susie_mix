# Small synthetic checks; does not access GTEx data or run the model packages.
source("script/analysis/fsusie_gtex_clean.R")
must_fail <- function(expr) {
  result <- tryCatch({ force(expr); FALSE }, error = function(e) TRUE)
  stopifnot(result)
}

# Include phased, unphased, missing, and partially missing calls. Preserve
# orientation and dimensions even with only one retained SNP or one sample.
gt <- matrix(c("0/0", "0|1", "1/0", "1|1", NA, "0/."), nrow = 2,
             dimnames = list(c("v1", "v2"), c("d3", "d1", "d2")))
X <- fsusie_decode_gt(gt)
expected <- matrix(c(0, 1, NA, 1, 2, NA), nrow = 3,
                   dimnames = list(c("d3", "d1", "d2"), c("v1", "v2")))
stopifnot(identical(X, expected))
stopifnot(identical(dim(fsusie_decode_gt(gt[1, , drop = FALSE])), c(3L, 1L)))
stopifnot(identical(dim(fsusie_decode_gt(gt[, 1, drop = FALSE])), c(1L, 2L)))
must_fail(fsusie_decode_gt(matrix("0/2", 1, 1)))

# The boundary bases must not be double-counted. Cover an indivisible length
# and a short locus where 1024 requested bins must fall back to a power of 2.
counts <- rbind(d1 = 1:11, d2 = rep(1, 11))
b <- fsusie_bin_counts(counts, 100:110, 4)
stopifnot(identical(unname(b$Y), rbind(c(6, 15, 24, 21), c(3, 3, 3, 2))))
stopifnot(identical(b$n_positions, c(3L, 3L, 3L, 2L)))
stopifnot(all(rowSums(b$Y) == rowSums(counts)), length(b$pos) == 4L)
stopifnot(all(diff(b$pos) == 11 / 4), b$start[1] == 100, tail(b$end, 1) == 111)
stopifnot(ncol(fsusie_bin_counts(counts, 100:110, 1024)$Y) == 8L)
stopifnot(ncol(fsusie_bin_counts(counts[, 1:8], 100:107, 1024)$Y) == 8L)
must_fail(fsusie_bin_counts(counts, c(100:109, 111), 4))
bad_counts <- counts
bad_counts[1, 1] <- NA
must_fail(fsusie_bin_counts(bad_counts, 100:110, 4))

# Library sizes follow RNA samples, even with multiple tissues for one donor,
# shuffled rows, and the dot-separated IDs used by the original script.
ids <- c("GTEX-A-001-SM-1", "GTEX-A-002-SM-2", "GTEX-B-001-SM-3")
sizes <- setNames(c(30, 10, 20), chartr("-", ".", ids[c(3, 1, 2)]))
stopifnot(identical(unname(fsusie_match_library_sizes(sizes, ids)), c(10, 20, 30)))
table <- matrix(sizes, ncol = 1, dimnames = list(names(sizes), "total"))
stopifnot(identical(unname(fsusie_match_library_sizes(table, ids)), c(10, 20, 30)))
table2 <- data.frame(unrelated = rep(99, 3), total = sizes, row.names = names(sizes))
must_fail(fsusie_match_library_sizes(table2, ids))
stopifnot(identical(unname(fsusie_match_library_sizes(table2, ids, "total")), c(10, 20, 30)))
must_fail(fsusie_match_library_sizes(sizes, "GTEX-X-001-SM-1"))
must_fail(fsusie_match_library_sizes(c(sizes, sizes[1]), ids))
sizes[1] <- 0
must_fail(fsusie_match_library_sizes(sizes, ids))

cat("PASS: genotype decoding, bin boundaries/count conservation, and RNA library-size matching.\n")
