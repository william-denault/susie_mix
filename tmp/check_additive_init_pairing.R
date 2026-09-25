.libPaths(c('C:/Users/willi/AppData/Local/R/win-library/4.6', .libPaths()))
try(source('script/sim/tests/test_additive_slide_init.R'))
print(all.equal(unname(fitted$fits$SuSiE$pip), reference$susie_pip))
print(formals(susieR::susie)$tol)
e <- new.env()
for (nm in names(formals(sim_mix))) assign(nm, eval(formals(sim_mix)[[nm]]), e)
e$pve <- .05; e$n <- 500; e$L_add <- 2; e$L_rec <- 0
e$seed <- 1000001; e$temp_dir <- genotypes
for (expr in as.list(body(sim_mix))[-1]) {
  if (is.call(expr) && identical(expr[[1]], as.name('<-')) &&
      identical(expr[[2]], as.name('susie_res'))) break
  eval(expr, e)
}
print(all.equal(e$geno_all, data$X))
print(all.equal(e$y, data$y))
