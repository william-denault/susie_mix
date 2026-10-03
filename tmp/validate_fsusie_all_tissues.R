# Exercise the actual edited preprocessing blocks without cluster files/packages.
lines <- readLines("script/analysis/fsusie_gtex_all_tissues.R")
block <- function(start, end) {
  first <- which(grepl(start, lines, fixed = TRUE))
  last <- which(grepl(end, lines, fixed = TRUE))
  stopifnot(length(first) == 1L, length(last) == 1L, first < last)
  parse(text = lines[first:(last - 1L)])
}
missing_filter <- block("keep_snp <- colSums(is.na(X)) == 0", "  ind_names=  tt[,1]")
alignment <- block("  # Match all tissue samples", "  res <- susiF")

fixture <- function(vector_sizes = FALSE, single_snp = FALSE) {
  e <- new.env(parent = globalenv())
  e$X <- data.frame(
    variant1 = c(2, 0, 1), missing = c(NA, 1, 0),
    fixed_ref = c(0, 0, 0), fixed_alt = c(2, 2, 2),
    variant2 = if (single_snp) c(0, 0, 0) else c(0, 2, 1),
    row.names = c("GTEX-B", "GTEX-A", "GTEX-C")
  )
  e$info_SNP <- data.frame(ID = names(e$X))
  e$tt <- data.frame(Sample_ID = c(
    "GTEX-B-001-SM-1.Brain_Cortex", "GTEX-A-001-SM-2.Brain_Cortex",
    "GTEX-B-002-SM-3.Lung", "GTEX-D-001-SM-4.Lung", "GTEX-C-001-SM-5.Skin"
  ))
  e$binned_data <- cbind(1:5, (1:5) * 10)
  sizes <- c(300, 200, 100, 50, 80)
  ids <- c("GTEX.C.001.SM.5", "GTEX.B.002.SM.3", "GTEX.A.001.SM.2",
           "GTEX.B.001.SM.1", "GTEX.D.001.SM.4")
  e$sum_count <- if (vector_sizes) setNames(sizes, ids) else
    matrix(sizes, ncol = 1, dimnames = list(ids, "total"))
  # Supply the fixture in place of the remote RData load.
  e$load <- function(...) invisible("sum_count")
  e
}

for (vector_sizes in c(FALSE, TRUE)) {
  for (single_snp in c(FALSE, TRUE)) {
    e <- fixture(vector_sizes, single_snp)
    eval(missing_filter, e)
    eval(alignment, e)
    stopifnot(nrow(e$X) == 4L, nrow(e$Y) == 4L)
    stopifnot(identical(unname(e$X[, 1]), c(2, 0, 2, 1)))
    stopifnot(identical(rownames(e$X), e$tt$Sample_ID[c(1, 2, 3, 5)]))
    stopifnot(identical(rownames(e$X), rownames(e$Y)))
    stopifnot(identical(e$sub_ind$tissue, c("Brain_Cortex", "Brain_Cortex", "Lung", "Skin")))
    stopifnot(identical(e$library_size, c(50, 100, 200, 300)))
    stopifnot(e$af[1] == 0.5)  # Three unique donors, not four tissue rows.
    stopifnot(ncol(e$X) == if (single_snp) 1L else 2L)
    stopifnot(identical(colnames(e$X), e$info_SNP$ID))
    expected <- sweep(e$binned_data[c(1, 2, 3, 5), , drop = FALSE], 1,
                      c(50, 100, 200, 300) / mean(c(50, 100, 200, 300)), "/")
    stopifnot(isTRUE(all.equal(unname(e$Y_cor), unname(expected))))
  }
}

e <- fixture()
e$sum_count <- e$sum_count[-1, , drop = FALSE]
eval(missing_filter, e)
error <- tryCatch({ eval(alignment, e); NULL }, error = conditionMessage)
stopifnot(grepl("Missing library sizes for", error, fixed = TRUE))
cat("PASS: all-tissue matching, repeated donors, unique-donor MAF, SNP metadata, and library-size alignment.\n")
