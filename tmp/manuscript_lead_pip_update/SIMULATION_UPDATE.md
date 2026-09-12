# Coding-discovery manuscript update - 11 September 2026

Updated the simulation Results, Methods, and Discussion using the completed
unweighted simulation archive. Added Table 1 (high-PIP global leads), Figure 5
(exact lead accuracy versus PIP, pooled over causal counts within PVE), and
Figures S16-S17 (additive-only lead frequencies and pooled coding PIP shares).
The earlier localization, coverage, ROC, and power-FDR results remain intact.

The updated text reports exact SNP AND coding matches, distinguishes coding
from localization errors, and separates unique-lead denominators from all-run
additive-control rates. It includes the PVE=10%, K=5 dominant-lead failure case,
the discrepancy between average PIP and observed exact accuracy, and the
nonzero coding PIP-mass baseline under additive truth. These results support
follow-up of confident non-additive candidates; they do not independently
validate GTEx signals, all credible-set labels, biological mixture proportions,
or weighted/EM-reweighted fits. No new simulations or fine-mapping fits ran.

Numerical sources: `simulation results/lead_pip_figures/lead_pip_analysis.rds`
and `simulation results/coding_figures/mixture_recovery.csv` in the analysis
project. `tmp/build_manuscript_lead_evidence.R` reproduces compact figures and
tables. The new lead figures report point estimates; the PIP-mass figure uses
the existing 500-draw shared-seed bootstrap intervals. High-PIP tables retain
both the pooled grid and individual PVE/K/called-coding conditions. New audit
files are prefixed `lead_pip_`, preserving the earlier localization audit.

The simulation review PDF is a reading copy of the updated simulation text,
table, and captioned vector figures, not a compilation of the full manuscript.
No TeX compiler is available locally. Source environments, labels, figure
paths, numerical counts, source preservation, and PDF rendering were checked.
The original template, introduction, bibliography, and unfinished GTEx counts
are preserved. This update supersedes the older statement that coding accuracy
or PIP reliability had not been evaluated.


---

Historical update record:

# Simulation manuscript update - 10 September 2026

This update supersedes the simulation description in the 8 September METHODS_AUDIT.md. The older n=200, 90-scenario benchmark was not the benchmark used for these results. The current completed benchmark is script/sim/sim_workhorse.R with 55 causal allocations, four PVE values, n=500, fitted L=10 and 400 scheduled replicates per allocation/PVE condition.

Updated susie_mix_manuscript.tex in place: simulation Results, Discussion, simulation Methods, four main PDF figures and fifteen supplementary PDF figures with numbered captions and references. Removed the template example figure/table and the template-only supplement to avoid collisions with real figure numbering. Repaired inherited undefined formatting macros (hetero/textcolor) and redirected the old encoding-definition Figure 1 pointer to Methods. The original science_template.tex remains intact. GTEx numerical placeholders were not filled, and unrelated introduction/acknowledgment prose remains a draft.

Evidence and computations:

- Analysis checkout: 8a76904e4dbc1285aea3915992cac42a50eb9889. Manuscript base: f0d21e8599e0e80da85cae7e5641dc064e8a77a1.
- Source methods: script/sim/sim_workhorse.R, script/sim/write_jobs.R, script/sim/jobs/conditions.csv and workhorse_utils.R in the simulation project.
- Results: simulation_tables/metric_summary.csv; status: file_audit.csv and configuration_counts.csv. There are 87,988 unique successes, twelve zero-genetic-variance errors and nine duplicate checkpoint entries; no included nonconverged fits. Total converged fits = 175,976. Conditions sharing a seed are dependent, not independent datasets for uncertainty estimation.
- Curve comparisons near 10% empirical FDR are calculated from roc_counts_all_L.rds. Starting at the no-discovery end, use the last threshold before the first empirical FDR > 0.10. The selected threshold, achieved empirical FDR, TP, FP and power are in power_fdr_first_crossing.csv. This is a descriptive truth-based operating point, not a fitted rule with FDR control.
- Original PDF files were copied byte-for-byte. simulation_tables/figure_manifest.json records the source hashes and manuscript figure numbers. The main figure numbering is coverage pure/mixed (1/2), then pooled power-FDR pure/mixed (3/4). Supplement S1/S2: CS recovery; S3/S4: purity; S5/S6: pooled ROC; S7-S11: pure ROC K=1..5; S12-S15: mixed ROC K=2..5.

Corrections from the obsolete design include: available .raw-file sampling rather than fresh gene extraction; n=500; the common retained-coding intersection for selecting causal SNPs; exhaustive count allocations; no locus redraws after an error; standardize=TRUE; optim prior-variance estimation; minimum purity 0.5; pooled rather than mean per-replicate FDP; and shared seeds 1,000,001 through 1,000,400. The simulation fitting controls differ from the GTEx scan and this is stated explicitly.

PVE uses independent Gaussian noise of variance 1-h and a genetic vector rescaled to sample variance h. It is not forced to equal realized var(g)/var(y). Equal variances of individual standardized contributions are not independent h/K shares in the presence of LD. Mixed SNP PIPs sum assignments within components before combining active components; mixed CSs collapse encodings within sets.

Current figure status and limits:

- The existing full-data figures have no uncertainty intervals and no CS-size comparison. The new plotting script supports seed-block bootstrap intervals and unique-SNP CS sizes, but its updated full run has not been executed. No subset-test figures have been included. The manuscript describes only the figures actually available.
- No PIP reliability/calibration diagram was generated; ROC and power-FDR measure discrimination. No formal FDR guarantee, coding-identification accuracy, or independent validation of human causal signals is claimed.
- Higher coverage is not claimed to establish tighter credible sets. The additive-case CS-recovery cost and remaining undercoverage are reported.
- The current simulation samples a fixed genotype-file collection and restricts recessive eligibility. Cluster software versions and an immutable file inventory are not recoverable from the saved compact objects. The outcomes, sampled donor rows, and full fitted models are not saved.
- The existing introduction still includes broad novelty assertions, unfinished prose and unresolved bibliography entries, and the abstract still contains a second-dataset placeholder. These require a separate literature/cohort revision; no evidence for them was invented.

The simulation_sections_review.pdf is a standalone reading copy of the newly written sections and their captioned vector figures, not a LaTeX compilation of the entire unfinished manuscript. It is generated from simulation_sections.json, the same section/caption content used in the LaTeX update. No TeX compiler was available locally. The LaTeX source is structurally checked, including figure references and paths. No simulations, fits, or full checkpoint plotting run were performed for this manuscript update.
