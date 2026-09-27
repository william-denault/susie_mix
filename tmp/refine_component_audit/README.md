# Extra-component and refinement audit — 2026-09-27

## Finding

The concern is real and is not automatically handled by the current coding-prior update. Every component with V > 0 contributes one total unit of alpha mass, irrespective of credible-set membership, purity, duplication, or the magnitude of V. This is the stated component-model M-step; it is not an estimator that counts only distinct, well-supported biological signals.

Changing min_abs_corr alone with refinement disabled changes reported credible sets, not the EM inputs. Enabling refine = TRUE changes the fitting search and can therefore change alpha, V, and the subsequent coding-prior update. Purity controls which sets seed that refinement search.

## Reproduction

Used installed susieR 0.16.6 and the user's phenotype/genotype construction. Removed only unrelated references to absent objects (res_susie_slide and X). Fixed the random seed for post-fit purity sampling. The phenotype contains three effects (SNPs 12, 92, 128); SNP 165 is listed but not added to y.

The supplied unrefined call produced five positive variances and four credible sets:

| Component | V | Reported CS at purity 0 | Additive mass | Dominant mass | Recessive mass |
|---|---:|---|---:|---:|---:|
| 1 | 4.036463 | yes | 1 | 0 | 0 |
| 2 | 0.607551 | yes | 0 | 0 | 1 |
| 3 | 3.995628 | yes | 0 | 1 | 0 |
| 4 | 1.428966 | no: duplicate of component 2 | 0 | 0 | 1 |
| 5 | 0.003470 | yes: 2,063 predictors, purity 0 | 0.424341 | 0.447699 | 0.127961 |

Thus it is component 4, not component 5, that is absent from the original CS list. Components 2 and 4 select the same five recessive predictors. Component 5's max alpha is only 0.02223 and its component log BF is 0.76323.

Post-hoc purity >= 0.5 leaves components 1, 2, and 3. That does not remove rows 4 or 5 from the stored alpha/V or from the current EM update.

## Matched refinement result

Fit susieR::susie(X, y, L=10, min_abs_corr=.5, refine=TRUE), retaining default convergence settings.

| Diagnostic | Unrefined | Refined |
|---|---:|---:|
| Final ELBO | -462.3865 | -452.3653 |
| Positive V components | 5 | 4 |
| CS passing purity .5 | 3 | 3 |
| One-step additive weight | 0.284868 | 0.356121 |
| One-step dominant weight | 0.289540 | 0.361961 |
| One-step recessive weight | 0.425592 | 0.281919 |

These are one-step normalized alpha sums for this one example, not converged genome-wide EM estimates.

Refinement removes the duplicated recessive component. The retained fourth component is diffuse (2,058 predictors), has V=0.00347217, and coding mass 0.424482/0.447843/0.127675. It fails the .5 purity filter but still contributes to the current EM M-step.

Four internal refinement candidates hit the 100-iteration limit. The final retained fit reports converged=TRUE; this does not establish a global optimum. Additional probes using max_iter=1000 and tol=1e-5 were interrupted after the matched default refinement completed; no result from those incomplete probes is used here.

## Interpretation

This establishes sensitivity of coding updates to the component decomposition in the supplied example. It does not establish how frequently this happens in GTEx, nor prove that GTEx recessive weights would increase after refinement. Here refinement actually reduces the recessive one-step weight because the initial split signal counted recessive twice.

It would be premature to interpret the current low recessive weights as biological prevalence without auditing component contributions and refinement sensitivity. A suitable GTEx audit needs saved full fits (alpha, V, component evidence, predictor maps, and CS component indices); the uploaded aggregate summaries do not contain enough detail.

Compare fixed-prior fits before/after refinement on a representative panel, including weak/null regions. Record coding mass for distinct high-purity sets, duplicates, diffuse components, and small positive V. If fitting changes materially affect aggregate coding weights, repeat the outer EM with the revised E-step rather than just rebuilding CS tables. A CS-only update can be shown as a sensitivity analysis, but changes the estimator and may select differently across coding types; it is not an automatic EM correction.

## Separate issues in the example

The genotype matrix is centered/scaled and some columns have four distinct values. Ranking genotype values by frequency does not preserve allele dosage when the heterozygote is most common. In addition, AD[AD==3] <- 2 is applied after DOM/REC construction, leaving 1,466 entries where REC != (AD==2). These issues do not invalidate the demonstrated component-count mechanism, but the example should be rebuilt from valid genotype dosage for biological coding validation.

## Files

- audit_default_refine.R: completed reproduction.
- default_refine_output.txt: complete output including warnings and component tables.
- user_fit.rds: initial fit plus X and y.
- refined_default.rds: final refined fit.
- audit_example.R / audit_output.txt / tight_unrefined_output.txt: exploratory tighter-convergence runs, interrupted as noted above.

Production analysis code and existing results were not modified.

