.libPaths(c(
  "C:/Users/willi/AppData/Local/R/win-library/4.6",
  "C:/Users/willi/AppData/Local/R/win-library/4.5", .libPaths()
))
library(susieR)
source("script/scan_tissue_attempt/workhorse_utils.R")

# Load each function with local utilities in place of its cluster source line.
old_env <- new.env()
new_env <- new.env()
eval(parse("C:/Users/willi/.codex/attachments/1b16c0d3-2a9a-49b4-8c69-a6a5e51da126/pasted-text.txt")[-1], old_env)
eval(parse("script/sim/sim_workhorse.R")[-1], new_env)

fixture_dir <- tempfile("compact_pip_genotypes_")
dir.create(fixture_dir)
set.seed(4321)
X <- sapply(seq(.06, .45, length.out = 30), function(maf) rbinom(600, 2, maf))
colnames(X) <- paste0("chr1_", 100001:100030, "_A_G_b38_A")
raw <- data.frame(FID = 1:600, IID = paste0("donor", 1:600),
                  PAT = 0, MAT = 0, SEX = 0, PHENOTYPE = -9, X)
write.table(raw, file.path(fixture_dir, "fixture.raw"),
            row.names = FALSE, quote = FALSE, sep = "\t")

fits <- list()
new_env$susie <- function(X, y, ...) {
  fit <- susieR::susie(X, y, ...)
  fits[[length(fits) + 1L]] <<- fit
  fit
}
args <- list(pve = .3, n = 500, L_add = 1, L_rec = 1, L_dom = 1,
             seed = 101, temp_dir = fixture_dir)
old <- do.call(old_env$sim_mix, args)
new <- do.call(new_env$sim_mix, args)

stopifnot(
  identical(new$metrics, old$metrics),
  identical(new$causal_snps, old$causal_snps),
  identical(new$causal_coding, old$causal_coding),
  identical(new$true_pos, old$true_pos),
  identical(new$true_pos_mix, old$true_pos_mix),
  identical(new$beta_standardized, old$beta_standardized),
  identical(new$susie_pip, unname(old$susie_pip)),
  identical(new$susie_mix_pip, unname(old$susie_mix_pip)),
  identical(new$cs_mix_as_additive_indices, old$cs_mix_as_additive_indices),
  identical(new$mix_to_add[new$true_pos_mix], new$true_pos),
  length(unique(new$causal_snps)) == 3,
  length(new$susie_mix_pip_snp) == length(new$susie_pip),
  all(new$susie_mix_pip_snp >= 0 & new$susie_mix_pip_snp <= 1),
  is.null(names(new$susie_pip)), is.null(names(new$susie_mix_pip)),
  is.null(names(new$susie_mix_pip_snp)),
  !any(c("y", "g", "noise", "alpha") %in% names(new)),
  new$settings$pve == .3, new$seed == 101
)
mix_names <- names(old$susie_mix_pip)
coding_counts <- table(sub(".*__", "", mix_names))
stopifnot(coding_counts["recessive"] < coding_counts["additive"])
cat("Unequal coding counts:", coding_counts, "\n")

# Check a hand-calculated union probability, including one inactive effect.
fake_fit <- structure(list(
  alpha = rbind(c(.2, .3, .4, .1), c(.1, .5, .25, .15), rep(.25, 4)),
  V = c(1, 1, 0), null_index = 0
), class = "susie")
fit_snp <- fake_fit
fit_snp$alpha <- cbind(rowSums(fake_fit$alpha[, c(1, 3, 4), drop = FALSE]),
                      fake_fit$alpha[, 2])
stopifnot(isTRUE(all.equal(susie_get_pip(fit_snp), c(.85, .65))))
fit_snp$V[] <- 0
stopifnot(identical(susie_get_pip(fit_snp), c(0, 0)))

# A pure additive simulation remains valid when no recessive column survives.
set.seed(4322)
X <- matrix(rbinom(200 * 8, 1, .3), 200, 8)
colnames(X) <- paste0("chr1_", 200001:200008, "_A_G_b38_A")
raw <- data.frame(FID = 1:200, IID = paste0("donor", 1:200),
                  PAT = 0, MAT = 0, SEX = 0, PHENOTYPE = -9, X)
write.table(raw, file.path(fixture_dir, "fixture.raw"),
            row.names = FALSE, quote = FALSE, sep = "\t")
empty_rec <- new_env$sim_mix(pve = .3, n = 200, L_add = 1, L_rec = 0,
                            L_dom = 0, seed = 102, hwe_thresh = 0,
                            temp_dir = fixture_dir)
stopifnot(length(empty_rec$susie_mix_pip) == 2 * length(empty_rec$susie_pip),
          identical(empty_rec$causal_coding, "additive"),
          all(is.finite(empty_rec$susie_mix_pip_snp)))

# Check that compact results survive the same save/load used by jobs.
checkpoint <- file.path(fixture_dir, "results.RData")
results <- list(new, empty_rec)
save(results, file = checkpoint)
loaded <- new.env()
load(checkpoint, envir = loaded)
stopifnot(identical(results, loaded$results))
cat("All compact-output, unchanged-fit, grouped-PIP and empty-coding checks passed.\n")
