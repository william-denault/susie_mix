# Finite slider prior: microbenchmark results

Run on 27 September 2026 using the `susie_slide_prior` branch. These results use a clean installed release build of the new implementation, compiled with `-O2` via `R CMD INSTALL --preclean`. Preliminary development-build timings are excluded.

The comparison uses `susieR::susie` from susieR 0.16.6, the original `susieSlide::susie` from susieSlide 0.2.0, and the finite-prior implementation from susieSlide 0.3.0. The original installed slider entry point was checked against branch `susie_slide`. Separate R processes load the two versions of susieSlide; workers run sequentially.

## Median full-fit time

Five repetitions per method and scenario, after one full warm-up fit. Values in parentheses are the interquartile range. All times are seconds.

| Samples | SNPs | Gaussian variances | Standard SuSiE | Original slider | 17-point prior |
|---:|---:|:---|---:|---:|---:|
| 500 | 1,000 | Fixed | 0.177 (0.176–0.181) | 1.185 (1.027–1.494) | 1.539 (1.506–1.621) |
| 500 | 1,000 | Estimated | 0.276 (0.275–0.282) | 5.257 (5.175–5.350) | 10.433 (10.407–10.774) |
| 1,000 | 3,000 | Fixed | 1.288 (1.281–1.288) | 7.125 (7.018–7.447) | 6.211 (6.096–7.781) |

The finite-prior/original-slider median ratios are 1.30, 1.98 and 0.87, respectively. Relative to standard SuSiE they are 8.68, 37.85 and 4.82. The larger fixed-variance comparison has overlapping interquartile ranges, so its lower finite-prior median should not be interpreted as an established speed advantage.

Variance optimization is expensive in the current implementation: each objective evaluation recomputes the finite mixture. The small estimated-variance fit takes roughly twice as long as the original slider even though it converges in half as many sweeps. This is consistent with additional mixture work during optimization; these timings are not a function-level profile. They do not support a general twofold-runtime claim relative to standard SuSiE.

## Iterations and returned-fit storage

Every warm-up fit converged. The timed calls use the same deterministic fitting inputs and settings.

| Samples × SNPs | Gaussian variances | SuSiE sweeps | Original-slider sweeps | Finite-prior sweeps |
|:---|:---|---:|---:|---:|
| 500 × 1,000 | Fixed | 8 | 8 | 7 |
| 500 × 1,000 | Estimated | 4 | 8 | 4 |
| 1,000 × 3,000 | Fixed | 11 | 8 | 8 |

These are complete fits with different statistical models, not a fixed number of identical update operations. In the estimated-variance case, standard SuSiE and the finite-prior fit retain three positive-variance effects; the original slider retains ten.

Returned-fit sizes at 500 × 1,000 are 0.71, 5.84 and 10.08 MiB for standard SuSiE, original slider and finite prior. At 1,000 × 3,000 they are 2.10, 27.40 and 40.34 MiB. These are `object.size` measurements of the returned objects, not peak memory usage.

## Matched settings and scope

- `L=10`, `tol=1e-6`, `max_iter=300`, intercept and additive-genotype standardization enabled.
- Complete numeric genotypes in `{0,1,2}`, simulated as independent Binomial(2, 0.3) draws. Three causal SNPs have coefficient effects `(0.9,-0.7,0.8)` and slider values `(-0.5,0,0.5)`; noise SD is 0.8. Each size uses one fixed simulated dataset shared by all methods.
- Initial coefficient-prior variance 0.5 and residual variance 0.64. The estimated setting updates both, using `estimate_prior_method="optim"`.
- The 17 slider values are evenly spaced from −1 to 1. Their probabilities remain fixed at `1/17`, including in the estimated-Gaussian-variance scenario. This benchmark does not learn the slider prior by EM.
- Slider defaults: `min_obs=5`, `chunk_size=1000`, `cache_heterozygotes=FALSE`. No SNP in these datasets is forced to additive coding.
- `coverage=NULL` excludes credible-set construction and purity calculations for every method.
- Timing includes input preparation, IBSS and output construction. It excludes simulation, R startup, package loading, file I/O and the explicit pre-measurement garbage collection (`setup=invisible(gc())`). Garbage collection during a fit remains included.
- `microbenchmark` 1.5.0; Windows 11, R 4.6.1; AMD64 Family 26 Model 112 Stepping 0. OpenMP/MKL/OpenBLAS thread environment variables were set to one. Exact session information is in the three `*_session.txt` files.

Five repetitions characterize execution variability on this machine; they do not measure variability across loci, allele frequencies, linkage-disequilibrium patterns or numbers of signals.

## Files and reproduction

- `benchmark_summary.csv`: medians, quartiles, range, iterations and runtime ratios.
- `benchmark_raw.csv`: all 45 timing measurements, in milliseconds.
- `benchmark_diagnostics.csv`: convergence, iterations, variances, active effects and returned-fit sizes.
- `*_microbenchmark.rds`: original microbenchmark objects.
- `benchmark_results.rds`: combined tables and case definitions.
- `benchmark_settings.txt` and `*_session.txt`: settings and environment.

The reusable script is saved in `C:/Document/Serieux/Travail/Package/git/susieR/inst/examples/benchmark_slider_prior.R` on branch `susie_slide_prior`. To reproduce this run in PowerShell:

```powershell
$env:SLIDER_BENCH_PRIOR_LIB = 'C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/tmp/slider_prior_impl/library'
$env:SLIDER_BENCH_LEGACY_LIB = 'C:/Users/willi/AppData/Local/R/win-library/4.6'
$env:SLIDER_BENCH_LEGACY_SOURCE = 'C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/tmp/slider_prior_impl/legacy_slider_source.R'
$env:SLIDER_BENCH_OUT = 'C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/output/benchmarks/slider_prior_2026-09-27/release'
$env:SLIDER_BENCH_PILOT = 'false'
$env:SLIDER_BENCH_REPS = '5'
& 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' --vanilla 'C:/Document/Serieux/Travail/Package/git/susieR/inst/examples/benchmark_slider_prior.R'
```

The new package and microbenchmark were installed in the isolated `SLIDER_BENCH_PRIOR_LIB` directory; the user's regular susieSlide installation remains version 0.2.0. For another machine, change library/output paths and install the required packages. `SLIDER_BENCH_LEGACY_SOURCE` is optional. Appending `--summarize` rebuilds the combined result tables from existing worker CSVs without rerunning fits.

Source provenance: branch base commit `e959d272c4abb426e620e3db4a375d2c8a5c42ab`, with the finite-prior changes uncommitted. SHA-256 hashes at measurement time:

```text
src/slider_prior.cpp
28699D93859EC2EC673811CACFE0D4DDAA8B6DE71E00C99F7550F854AF5A5CD3
R/slider_prior.R
893885476E8454120A78E86187DD9CD9BB8F12927F4FC18B28FA2F6D79A2A6A8
R/slider_engine.R
6522B75847EF4C794FED799B8079809F42280AD613FEFFAFF7CDA4CDA6D5F570
inst/examples/benchmark_slider_prior.R
93238006A9CB41E22D2D3F38DC65525B9C0F0475546C71519F24C92C9A482D2E
```
