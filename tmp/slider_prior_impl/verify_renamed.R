lib <- "C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/tmp/slider_prior_impl/library"
root <- "C:/Document/Serieux/Travail/Package/git/susieR"
legacy <- loadNamespace("susieSlide",lib.loc="C:/Users/willi/AppData/Local/R/win-library/4.6")
library(susieRSlidePrior,lib.loc=lib)
stopifnot(identical(unname(getNamespaceName(legacy)),"susieSlide"),
          identical(environmentName(environment(susieRSlidePrior::susie)),"susieRSlidePrior"),
          all(c("susieRSlidePrior","susieSlide") %in% names(getLoadedDLLs())),
          length(utils::help("susieRSlidePrior",package="susieRSlidePrior"))==1L)
cat("Renamed namespace, native library, help and coexistence with susieSlide: PASS\n")
res <- testthat::test_dir(file.path(root,"tests/testthat"),
  filter="^slider",package="susieRSlidePrior",load_package="installed",
  reporter="summary",stop_on_failure=TRUE)
s <- as.data.frame(res)
cat("\nInstalled renamed package test summary:\n")
print(colSums(s[,intersect(c("nb","failed","error","warning","skipped"),names(s)),drop=FALSE]))
