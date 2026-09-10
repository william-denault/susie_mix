> Simulation description superseded on 10 September 2026; see [SIMULATION_UPDATE.md](SIMULATION_UPDATE.md) for the completed benchmark, results, and figure inventory. The historical audit below describes the earlier design.

Methods drafting audit — 8 September 2026

The revised manuscript is [susie_mix_manuscript.tex](C:/Document/Serieux/Travail/Package/git/Science-SuSiE-mix/susie_mix_manuscript.tex:287). It contains a replacement GTEx Methods section and the implemented simulation design. Numerical results remain placeholders; fixed design parameters are specified. The original science_template.tex is preserved. Its accompanying readme requires modified template files to be renamed.

The author confirmed that the current workhorse is the intended analysis: broad SMTS tissue grouping, two inverse-normal rank transformations, and sex-only expression adjustment. This establishes the intended method; it does not independently establish which code version generated every saved result file.

The analysis checkout was at commit e130e7871f22ebeaf10350103030e395a4c620a6; the original manuscript checkout was at 245722dd6e6331613e83f685f4a1cc58c06affa6. There were no local changes to the inspected analysis scripts. The analysis code and result files were not edited.

| Description | Implementation inspected |
| --- | --- |
| GTEx inputs, sample selection, normalization, filtering, fits, priors, permutations | [workhorse.R](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/script/scan_tissue_attempt/workhorse.R:9) |
| Minor-allele orientation, HWE exclusion, genotype encodings, lead positions | [workhorse_utils.R](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/script/scan_tissue_attempt/workhorse_utils.R:187) |
| Gene annotation and protein-coding gene selection | [get_gene_annotations.R](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/script/scan_tissue_attempt/get_gene_annotations.R:3), [get_protein_coding_gene_names.R](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/script/scan_tissue_attempt/get_protein_coding_gene_names.R:1) |
| Production gene dispatch | [run_chunk_001.R](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/script/scan_job/run_chunk_001.R:1), [write_job.R](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/script/scan_tissue_attempt/write_job.R:1) |
| CS membership, lead selection, coding counts, overlaps | [generate_summary_results.R](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/script/analysis/generate_summary_results.R:613) |
| Descriptive thresholds and analysis subsets | [descriptive results.R](<C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/script/analysis/descriptive results.R:32>), [descriptive_results_weighted.R](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/script/analysis/descriptive_results_weighted.R:29) |
| Simulation constants and scenario grid | [config.R](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/script/simulation/config.R:2), [simulation_utils.R](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/script/simulation/simulation_utils.R:3) |
| Causal selection, phenotype generation, fitting and metrics | [simulation_utils.R](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/script/simulation/simulation_utils.R:126) |
| Completion rules and uncertainty estimates | [summarize_simulation.R](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/script/simulation/summarize_simulation.R:4) |

Material corrections to the original Methods:

- HWE P < 1e-8 variants are removed by the current QC function. The older frequency-recoding helper still exists, but the current workhorse does not call it.
- The expression phenotype is library-size normalized, log1p transformed, rank-inverse-normal transformed, residualized on sex, and rank-inverse-normal transformed again.
- SMTS defines broad tissue groups. One sample per donor is retained by the first match; detailed SMTSD subtypes are not separate analysis units.
- Both descriptive scripts use minimum additive marginal P < 1e-8 and mean raw count >= 100. The corresponding threshold in the revised Results paragraph was corrected from 5e-8.
- The predictor column-sum filter means allele copies for additive predictors, carriers for dominant predictors, and homozygotes for recessive predictors. It is not a requirement that every genotype class contain five donors.
- Ordinary mixed SuSiE uses uniform predictor priors. The additional weighted fit uses total coding-class masses 0.80 additive, 0.15 dominant, and 0.05 recessive, divided by retained class counts.
- There is one paired additive/unweighted-mixed phenotype permutation per gene–tissue pair. There is no weighted permutation fit.
- Lead coding is determined by the largest component-specific assignment probability within the CS, with a predictor-PIP fallback for older objects.
- Stored TSS displacement is genomic position minus TSS, rather than absolute or strand-oriented distance.
- GTEx descriptive summaries do not explicitly exclude nonconverged fits. Simulation comparisons do require both methods to converge.

The implemented simulation has 90 scenarios, n = 200 and L = 10 throughout, three PVE levels, and 200 configured replications per scenario. It compares additive SuSiE and ordinary unweighted mixed SuSiE. Recessive causal selection is conditional on surviving the homozygote-count filter. Equal standardized effect magnitudes with random signs are rescaled jointly to account for LD, and independent Gaussian errors are generated without orthogonalization.

The production directory inspected locally contains a plan, but no completed production fit archive or completion report. The pilot status explicitly reports an incomplete pilot. Accordingly, the manuscript describes the implemented benchmark without claiming that all production simulations have run. Completion may exist on the cluster; it was not established from this checkout.

Items still needed to finalize the Methods are the actual GTEx run's R, susieR and PLINK2 versions; its effective default fitting controls; run/code identifiers linking the reported tables to that run; data-access and ethics information; and the completed production simulation manifest and completion report. Unspecified GTEx defaults should be recovered from that environment rather than copied from the current website. The [susieR documentation](https://stephenslab.github.io/susieR/reference/susie.html) was checked for the meaning of the fitting arguments.

Scientific assessment:

Nature Genetics is the strongest match among the proposed journals. Its [scope](https://www.nature.com/ng/aims) covers genetics and functional genomics, and it has published a directly relevant example of improving fine-mapping by addressing model misspecification: [Improving fine-mapping by modeling infinitesimal effects](https://www.nature.com/articles/s41588-023-01597-3). This is a precedent for the kind of contribution, not a prediction of editorial acceptance.

The strongest possible contribution is evidence that a widespread modeling assumption can confidently prioritize the wrong variants, and that accounting for heterogeneous genotype effects restores calibration and changes validated biological interpretation. A simple correction can support a major paper if that consequence is established. My present assessment is that Nature Genetics is a plausible ambitious target; Science and Nature remain stretch targets requiring broader and independently validated consequences. Nature explicitly considers importance and interdisciplinary interest in its [editorial criteria](https://www.nature.com/nature/for-authors/editorial-criteria-and-processes).

The current manuscript and available local outputs do not yet establish that the mixed model identifies the correct regulatory SNP in humans. Model disagreement, a different lead SNP, or a higher variational objective is not ground truth. The saved [agreement table](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/descriptive_results/agreement_statistics.csv) reports SNP overlap for every single-CS-versus-single-CS pair in the primary subset. These tables therefore do not currently support a claim of completely disjoint single-CS results for that analysis. They can still support changes in set membership, lead variants, or signal counts, once tied to a verified run.

The key work needed to substantiate the strongest claims is:

- Complete the paired simulation benchmark and establish coverage, PIP calibration, false discovery rates, and power. Distinguish the reported empirical fraction of CSs containing any causal SNP from a guarantee about each underlying signal.
- Show that the key GTEx conclusions survive detailed-tissue analysis and adequate covariate control. The present sex-only adjustment and broad tissue grouping leave plausible alternative explanations involving ancestry, technical effects, subtype composition, and donor/sample selection.
- Assess the role of expression scale. Two nonlinear rank transformations can change genotype–phenotype relationships, and the second transformation need not preserve orthogonality to sex. The current simulated phenotype does not undergo this preprocessing.
- Establish robustness to predictor scaling, coding priors, optimization and CS purity. With expanded predictors, G = dominant + recessive, and distinct SuSiE components may select different encodings of the same SNP. Thus CS count is not automatically a count of distinct causal SNPs, and lead coding is not a uniquely identified mechanism.
- Validate high-confidence discordant variants in an independent cohort and/or with credible variant-level functional evidence. Demonstrate a consequential downstream change, such as corrected colocalization or regulatory interpretation. A larger mixed-model objective alone cannot establish this.
- Position the novelty precisely. [Peltola et al.](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0029115) already considered Bayesian variable selection across genetic effect models, and [Palmer et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC10345642/) performed dominance fine-mapping with SuSiE. The manuscript's blanket statement that a framework does not exist should be revised in a later introduction edit. This was a targeted literature check, not an exhaustive novelty review.

Validation of this editing task:

The scenario generator and configuration validation ran successfully in base R and confirmed the stated grid and production design size. The new Methods were checked against the functions rather than inferred from README prose alone. Result numbers were not added. A static check of new LaTeX structure and cited Methods keys was performed. No scientific analysis was rerun.

No LaTeX compiler was found, so no new PDF is claimed. The inherited manuscript remains a draft outside the replaced Methods: it still contains template example figures/tables, unfinished prose, undefined commands (including hetero and textcolor), and unresolved bibliography keys including wenstephens2014, sharepro2024, peltola_bayesian_2012 and sabourin_fine-mapping_2015. These are pre-existing issues, not missing Methods content. The new Methods cite only GTEx and Wang et al. keys already present in library.bib.
