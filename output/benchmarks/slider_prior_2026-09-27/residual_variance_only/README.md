# Cost of updating only the residual variance

This follow-up separates the residual variance sigma2 from the Gaussian effect-prior variances V_l. Every V_l is fixed at 0.5 in both settings. The only changed option is `estimate_residual_variance`, which updates sigma2 after each non-final IBSS sweep. The slider prior probabilities remain fixed at 1/17. No numerical prior-variance optimization runs.

All methods use the same 500 samples, 1,000 hard-call SNPs and simulated phenotype as the earlier benchmark, with L=10 and tol=1e-6. Five microbenchmark repetitions per scenario follow one complete warm-up fit. Methods run in separate sequential R processes; data simulation, package loading, I/O and the explicit pre-fit garbage collection are excluded. The source script is `tmp/slider_prior_impl/benchmark_residual_variance.R` relative to the susie_mix workspace.

| Implementation | Fixed sigma2: median seconds | Updated sigma2: median seconds | Sweeps, fixed / updated |
|:---|---:|---:|:---|
| susieR 0.16.6 | 0.179 | 0.171 | 8 / 7 |
| Original susieSlide 0.2.0 | 0.813 | 0.970 | 8 / 10 |
| susieRSlidePrior 0.3.0 | 1.058 | 1.200 | 7 / 8 |

Every warm-up fit converged. For the finite-prior model, full-fit median time divided by the number of sweeps is 0.151 seconds with fixed sigma2 and 0.150 seconds with updated sigma2. These ratios include amortized setup and output construction. The measured total-time increase is consistent with one additional sweep, rather than a costly residual-variance update.

The earlier 10.43-second result also used `estimate_prior_variance=TRUE, estimate_prior_method="optim"`. That setting numerically optimizes each V_l inside each component update, evaluating the 17-point mixture repeatedly. It therefore does not isolate the cost of updating sigma2 after a sweep.

The existing `estimate_prior_method="EM"` option provides a posterior-second-moment update for V_l without the numerical optimization loop. In the current slider implementation it runs after each component's posterior update, followed by one posterior refresh at the new variance. This follow-up did not benchmark that option or change the fitting defaults.

The 30 raw timings, summary statistics, convergence diagnostics, settings, microbenchmark objects and R session information are saved alongside this report. These are measurements on one synthetic dataset and machine; the fits can require different numbers of sweeps.
