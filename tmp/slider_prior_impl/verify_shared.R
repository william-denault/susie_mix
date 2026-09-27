res <- testthat::test_local(
  path="C:/Document/Serieux/Travail/Package/git/susieR",
  filter="susie_get_functions|single_effect_regression|individual_data_methods",
  reporter="summary", stop_on_failure=TRUE)
s <- as.data.frame(res)
cat("\nShared-engine test summary:\n")
print(colSums(s[,intersect(c("nb","failed","error","warning","skipped"),names(s)),drop=FALSE]))
