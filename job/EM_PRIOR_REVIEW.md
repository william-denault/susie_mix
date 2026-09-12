# Empirical Bayes estimation of SuSiE coding weights

Reviewed 2026-09-12 against Zhao et al., *Adjusting for genetic confounders
in transcriptome-wide association studies improves discovery of risk genes
of complex traits*, Nature Genetics 56, 336-347 (2024),
[doi:10.1038/s41588-023-01648-9](https://www.nature.com/articles/s41588-023-01648-9).
The supplied article's PDF page 13 contains Methods equations (3)-(6).
The supplied supplement's PDF pages 2-5 (printed pages 1-4) contain the EM
derivation and the single-effect approximation, equations (2)-(31).

## Finding

The earlier update, summed coding PIPs divided by total PIP, was a useful
descriptive allocation statistic but was not the EM M-step for the fitted
multi-effect SuSiE model. It implemented the originally requested formula;
it must not be presented as a marginal-likelihood-maximizing EM algorithm.

The corrected implementation retains the requested three tissue-specific
coding mixture weights. It estimates them by variational empirical Bayes,
using SuSiE component responsibilities in the M-step. It does not silently
replace the model with cTWAS's independent Bernoulli inclusion model.

## What the paper estimates

cTWAS assigns a Bernoulli indicator to each predictor. For a predictor in
group k, the probability of a nonzero effect is rho_k. Its complete-data
prior term is

    sum_j [gamma_j log(rho_k) + (1-gamma_j) log(1-rho_k)].

Taking expectations and maximizing gives

    rho_k(new) = sum_{j in group k} PIP_j / number_of_predictors_in_group_k.

These probabilities are per predictor; they do not have to sum to one
across groups. The paper also estimates group-specific effect variances
using posterior second moments (Methods equation 6; supplement equation 16).
The supplement uses a single-effect approximation for posterior calculations,
including a null configuration, and describes the approximation involved in
using SuSiE SER (supplement equations 27 and 30). The real-data analysis uses
L=1 for prior estimation and L=5 for subsequent fine-mapping.

Thus the paper's use of PIPs is not an error. Its latent model and parameter
definition differ from ours. In particular, the paper's alpha_j notation for
a Bernoulli PIP must not be confused with SuSiE's L-by-p `alpha` matrix.

## Our model and its M-step

For each gene g in tissue t, SuSiE writes its effect vector as a sum of L_g
single effects. Component l selects one predictor J_gl. Let c(j) denote its
coding, n_gc the number of retained predictors of coding c, and A_g the set
of coding classes available after QC. We retain the existing prior mapping:

    pi_t,add + pi_t,rec + pi_t,dom = 1
    w_gj = pi_t,c(j) / [n_g,c(j) * sum_{d in A_g} pi_t,d].

The quantity w_gj is the categorical selection probability for a predictor
within one SuSiE component. The pi parameters are coding-class masses, not
per-predictor Bernoulli inclusion probabilities and not phenotype-variance
fractions. Unequal coding block sizes are accounted for by n_gc.

The E-step supplies

    alpha_glj = q(J_gl = j),
    C_gc = sum_l sum_{j: c(j)=c} alpha_glj.

Holding the approximate posterior and other parameters fixed, the part of
the expected complete log likelihood involving the coding weights is

    Q_t(pi) = sum_g sum_{c in A_g} C_gc log(pi_t,c)
              - sum_g L_g log(sum_{c in A_g} pi_t,c) + constant.

The within-class log(n_gc) terms are constant with respect to pi. If every
fit contains all three classes, the second term is zero on the simplex,
and the M-step is

    pi_t,c(new) = sum_g C_gc / sum_g L_g.

Equivalently, normalize the three summed alpha masses. If some classes are
absent, simply normalizing global counts is generally wrong. The corrected
code groups sufficient statistics by the seven nonempty availability sets
and optimizes the conditional-choice objective in log weights. For example,
posterior A:R counts of 8:2 in A/R-only fits and R:D counts of 2:8 in R/D-only
fits imply relative weights 4:1:4, not the naive pooled counts 8:4:8.
Disconnected availability groups have unidentified relative total mass;
the optimizer preserves their previous mass. Every update checks that Q
has not decreased within numerical tolerance.

## Why marginal PIPs differ

SuSiE's marginal PIP combines component assignments through a union:

    PIP_j = 1 - product_l (1 - alpha_lj),

with SuSiE's usual exclusion of components whose prior variance is near
zero. It is not the expected number of component assignments to predictor j.
For two rows both equal to (0.8, 0.1, 0.1), the M-step is (0.8, 0.1, 0.1).
The PIP vector is (0.96, 0.19, 0.19), whose normalized shares are approximately
(0.716, 0.142, 0.142). Those shares do not maximize the categorical Q.

All alpha rows are included, without credible-set, lead-PIP, or variance
filtering. Under the chosen latent-variable representation, a component
with V_l=0 still has a categorical assignment; its posterior is normally
its prior, so it provides no information favoring a coding. Including many
such components can make EM move slowly. Filtering them out would require
a separately derived collapsed update, not an arbitrary threshold.

## Likelihood target and practical limits

The target is the sum of gene-level log marginal likelihoods within a tissue,
with gene fits treated as independent in the working model. Genes can have
correlated expression and overlapping cis-regions, so this is not a fitted
joint expression model. The same SNP associated with different genes is a
different predictor-outcome occurrence; the learned weights do not describe
the prevalence of distinct biological SNPs in each coding class.

For multiple effects, SuSiE supplies a variational posterior. Alternating
IBSS with the corrected M-step therefore targets an evidence lower bound
(ELBO), not an exactly evaluated genome-wide marginal likelihood. Neither
this procedure nor ordinary EM guarantees a global maximum. Tight inner
convergence and small changes across several outer iterations should be
assessed; a fixed requested iteration count is not a convergence guarantee.

The history includes `mstep_q_gain`, `max_abs_prior_change`, and
`source_elbo_sum`. The latter is present only when every included fit supplies
a finite final ELBO; `n_source_elbo` records the count. Compare summed ELBOs
only when the same fits and model settings contribute. It is a variational
objective diagnostic, not an exact marginal-likelihood calculation.

SuSiE's `estimate_prior_method="EM"` controls estimation of its component
effect variances. It does not learn coding inclusion weights by itself.
The new outer update learns those weights; component variances and residual
variances remain gene/tissue-specific nuisance parameters estimated by SuSiE.
We do not add cTWAS-style shared coding-specific slab variances in this change.
The existing standardization and QC choices are retained; the learned coding
weights remain conditional on those effect-size prior and data-model choices.

## Implementation and existing results

- `em_utils.R`: validate full component posteriors, accumulate alpha counts,
  maximize the coding Q, retain PIP sums as diagnostics, and construct safe
  posterior initializations.
- `prepare_em_iteration.R`: append corrected priors and diagnostics; save
  `component_counts.rds`; preserve all old frozen iteration files and values.
- `workhorse_em.R`: initialize from the preceding posterior and nuisance
  variances, use the new coding weights, and verify the weights were retained.
  Initialization omits the old `pi` field because some SuSiE versions would
  otherwise overwrite the new `prior_weights`. All L rows are preserved.
- `run_em_chunk.R`: read the source fit, refuse completion on unconverged
  fits, and propagate invalid initialization or weight errors as job failures.

When a new iteration is prepared, old history rows receive the method label
`legacy_pip_share`; new rows receive `susie_component_alpha_v1`. Old numeric
values are preserved, and unavailable new diagnostics are left empty.
Completed legacy fits are valid starting values conditional on their old
priors. They need not be deleted or rerun just to initialize the corrected
algorithm. If the original unweighted scan is the only source, its fit is
an initialization rather than a previous shared-coding-prior EM iteration.

The locally synchronized iteration 2 contains completion records and logs,
but not the full gene RDS files. Consequently no new tissue weights were
computed in this review. The next cluster preparation will read those full
fits and calculate the corrected update before submitting fine-mapping.

## Validation

- Synthetic integration checks cover pooling, missing/corrupt source files,
  predictor alignment, legacy-history migration, immutable snapshots, and
  resuming unconverged fits.
- Analytic examples distinguish alpha counts from PIP shares and verify the
  optimizer for missing and disconnected coding classes.
- An independent Gaussian model with two effects enumerates every assignment
  and integrates effect sizes analytically. Its exact log marginal likelihood
  increased from -29.216109 to -27.995768 across 50 EM updates, with no decreases.
- Four small fits using the local SuSiE 0.14.21 numerical source exercised ten
  outer updates. Total ELBO increased from -935.103109 to -933.276036, with
  no decreases; each fit retained its new weights and all four component rows.
  Only the `colSds` dependency used base R in this source-only smoke test;
  SuSiE's fitting and objective functions were not substituted.
- Mock Slurm checks exercise the launcher without submitting cluster jobs.

These checks establish the update mathematics and exercise the implementation;
they are not a new GTEx analysis or a simulation study of biological calibration.
