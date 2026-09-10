# Recessive one-CS examples: lead SNP distances and colocalization assessment

Reviewed 2026-09-09. Scope: the 141 plotted gene–tissue examples in `plot/one_cs_recessive`, comparing ordinary additive SuSiE with the unweighted SuSiE-mix fit.

**84 examples have different lead SNPs; 57 have the same lead. Eight changed leads are more than 100 kb apart.** All 84 changed pairs are on the same chromosome. Coordinates below are the b38 coordinates encoded in the saved SNP identifiers; no liftover was performed. Distance is the absolute difference in those base-pair positions, not the separation of columns on the stacked predictor plot.

Source: [plot summary](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/fit_mix_one_cs_recessive_plot_summary.csv). Credible-set overlap comes from the existing [biological-SNP overlap report](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/descriptive_results/one_cs_overlap_and_lead_detail.csv), joined on gene, tissue and both exact lead identifiers. A missing match is reported as unavailable, never as zero overlap.

**Interpretation of the largest changes**

- ZDHHC2 in blood has the largest lead shift, 438,420 bp. Its lead PIPs are 0.285 and 0.267. The mixed CS contains six unique biological SNPs, including all five additive-CS SNPs. Its eight plotted mixed predictors include multiple codings. This is a changed ranking within strongly overlapping candidate sets; it does not establish an independent new signal. The plotted descriptive fit metric is -0.032, providing no positive fit advantage by that metric.
- GUSB in thyroid shifts 366,690 bp, but the mixed lead PIP is only 0.061. HMG20A in blood vessel shifts 275,824 bp, with mixed lead PIP 0.076. Their distances are striking, but localization remains uncertain.
- NEK3 in lung shifts 181,050 bp. The plot reports 54 additive-CS SNPs and 19 mixed-CS predictor columns; the mixed lead PIP is 0.224 and the mixed lead has 65 genotype-2 individuals. The expression pattern merits follow-up, but overlap and LD need examining before claiming a different causal explanation.
- ZNF354B in skin shifts 48,377 bp and is the most distant changed-lead example with a mixed lead PIP at least 0.5. The respective lead PIPs are 0.574 and 0.719, with four CS predictors in each fit and 69 genotype-2 individuals at the mixed lead. This is a useful candidate for checking posterior relocation and colocalization. Its CS overlap is unavailable in the existing overlap report; disjointness has not been established.
- PM20D1 in stomach shifts 85,615 bp, but its mixed lead PIP is 0.084, and all 15 additive-CS SNPs remain among the 24 mixed-CS biological SNPs. The Alzheimer-related biology makes it relevant to investigate, but these results alone do not demonstrate improved localization.

The PIP columns below are the saved PIPs for each model's own lead predictor. For SuSiE-mix these are coding-specific values, not PIPs marginalized over coding. A change between these two values is not the change in PIP for a single identical SNP. CS sizes in the main table are predictor counts; the overlap column uses unique biological SNPs.

## All changed leads, ordered by distance

| Rank | Gene — tissue / plot | Additive lead (GRCh38) | Mixed recessive lead (GRCh38) | Distance (bp) | Lead PIP: additive → mixed | CS predictors: additive → mixed | Shared / additive / mixed unique CS SNPs |
|---:|---|---|---|---:|---|---|---|
| 1 | [ZDHHC2 — Blood](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/ZDHHC2_Blood_add_vs_mix.png) | chr8:17,597,135 | chr8:17,158,715 | 438,420 | 0.285 → 0.267 | 5 → 8 | 5 / 5 / 6 |
| 2 | [GUSB — Thyroid](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/GUSB_Thyroid_add_vs_mix.png) | chr7:65,597,145 | chr7:65,963,835 | 366,690 | 0.230 → 0.061 | 54 → 88 | Unavailable |
| 3 | [HMG20A — Blood Vessel](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/HMG20A_Blood_Vessel_add_vs_mix.png) | chr15:77,247,572 | chr15:77,523,396 | 275,824 | 0.193 → 0.076 | 58 → 72 | Unavailable |
| 4 | [NEK3 — Lung](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/NEK3_Lung_add_vs_mix.png) | chr13:51,978,966 | chr13:52,160,016 | 181,050 | 0.309 → 0.224 | 54 → 19 | Unavailable |
| 5 | [AGAP4 — Stomach](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/AGAP4_Stomach_add_vs_mix.png) | chr10:46,027,444 | chr10:45,848,036 | 179,408 | 0.056 → 0.030 | 44 → 84 | 44 / 44 / 47 |
| 6 | [ANAPC4 — Pancreas](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/ANAPC4_Pancreas_add_vs_mix.png) | chr4:25,227,787 | chr4:25,393,278 | 165,491 | 0.044 → 0.044 | 71 → 116 | 68 / 71 / 76 |
| 7 | [FTSJ3 — Muscle](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/FTSJ3_Muscle_add_vs_mix.png) | chr17:63,877,558 | chr17:63,721,364 | 156,194 | 0.015 → 0.037 | 184 → 241 | 177 / 184 / 211 |
| 8 | [LIMD1 — Pancreas](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/LIMD1_Pancreas_add_vs_mix.png) | chr3:45,768,337 | chr3:45,639,089 | 129,248 | 0.133 → 0.158 | 20 → 26 | Unavailable |
| 9 | [P2RX2 — Spleen](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/P2RX2_Spleen_add_vs_mix.png) | chr12:132,704,876 | chr12:132,610,104 | 94,772 | 0.030 → 0.037 | 76 → 141 | Unavailable |
| 10 | [PM20D1 — Stomach](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/PM20D1_Stomach_add_vs_mix.png) | chr1:205,745,244 | chr1:205,830,859 | 85,615 | 0.123 → 0.084 | 15 → 33 | 15 / 15 / 24 |
| 11 | [ENOX1 — Heart](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/ENOX1_Heart_add_vs_mix.png) | chr13:43,917,093 | chr13:43,839,387 | 77,706 | 0.076 → 0.086 | 42 → 87 | Unavailable |
| 12 | [TMEM150C — Thyroid](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/TMEM150C_Thyroid_add_vs_mix.png) | chr4:82,655,676 | chr4:82,585,679 | 69,997 | 0.089 → 0.064 | 37 → 56 | Unavailable |
| 13 | [RBM6 — Blood Vessel](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/RBM6_Blood_Vessel_add_vs_mix.png) | chr3:49,892,782 | chr3:49,959,631 | 66,849 | 0.058 → 0.072 | 68 → 169 | 62 / 68 / 123 |
| 14 | [LGALS9C — Spleen](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/LGALS9C_Spleen_add_vs_mix.png) | chr17:18,572,316 | chr17:18,638,462 | 66,146 | 0.183 → 0.112 | 12 → 16 | 7 / 12 / 16 |
| 15 | [ASNSD1 — Prostate](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/ASNSD1_Prostate_add_vs_mix.png) | chr2:189,886,437 | chr2:189,829,449 | 56,988 | 0.201 → 0.141 | 65 → 82 | Unavailable |
| 16 | [THBS4 — Lung](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/THBS4_Lung_add_vs_mix.png) | chr5:80,034,736 | chr5:79,978,423 | 56,313 | 0.121 → 0.255 | 66 → 123 | 44 / 66 / 96 |
| 17 | [ZNF567 — Thyroid](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/ZNF567_Thyroid_add_vs_mix.png) | chr19:36,699,208 | chr19:36,643,241 | 55,967 | 0.021 → 0.029 | 92 → 163 | Unavailable |
| 18 | [GGCX — Spleen](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/GGCX_Spleen_add_vs_mix.png) | chr2:85,534,156 | chr2:85,582,866 | 48,710 | 0.041 → 0.106 | 44 → 47 | 42 / 44 / 47 |
| 19 | [ZNF354B — Skin](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/ZNF354B_Skin_add_vs_mix.png) | chr5:179,125,714 | chr5:179,077,337 | 48,377 | 0.574 → 0.719 | 4 → 4 | Unavailable |
| 20 | [ANKRD34B — Blood](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/ANKRD34B_Blood_add_vs_mix.png) | chr5:80,609,592 | chr5:80,654,678 | 45,086 | 0.022 → 0.050 | 211 → 212 | 188 / 211 / 206 |
| 21 | [IRF6 — Nerve](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/IRF6_Nerve_add_vs_mix.png) | chr1:209,805,432 | chr1:209,846,356 | 40,924 | 0.159 → 0.177 | 18 → 24 | 18 / 18 / 20 |
| 22 | [C19orf71 — Ovary](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/C19orf71_Ovary_add_vs_mix.png) | chr19:3,548,059 | chr19:3,511,948 | 36,111 | 0.186 → 0.156 | 16 → 27 | Unavailable |
| 23 | [SLC2A9 — Thyroid](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/SLC2A9_Thyroid_add_vs_mix.png) | chr4:9,994,227 | chr4:10,029,444 | 35,217 | 0.045 → 0.023 | 127 → 230 | Unavailable |
| 24 | [CENPQ — Lung](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/CENPQ_Lung_add_vs_mix.png) | chr6:49,405,083 | chr6:49,438,635 | 33,552 | 0.038 → 0.021 | 41 → 74 | 41 / 41 / 44 |
| 25 | [LY6G5C — Spleen](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/LY6G5C_Spleen_add_vs_mix.png) | chr6:31,694,452 | chr6:31,661,319 | 33,133 | 0.092 → 0.057 | 15 → 24 | 15 / 15 / 16 |
| 26 | [TPSB2 — Adipose Tissue](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/TPSB2_Adipose_Tissue_add_vs_mix.png) | chr16:1,230,423 | chr16:1,255,825 | 25,402 | 0.777 → 0.179 | 8 → 8 | 4 / 8 / 8 |
| 27 | [DOK7 — Spleen](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/DOK7_Spleen_add_vs_mix.png) | chr4:3,508,871 | chr4:3,484,033 | 24,838 | 0.245 → 0.192 | 16 → 5 | 4 / 16 / 5 |
| 28 | [LGALS3 — Pituitary](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/LGALS3_Pituitary_add_vs_mix.png) | chr14:55,120,579 | chr14:55,144,222 | 23,643 | 0.088 → 0.084 | 93 → 140 | 91 / 93 / 92 |
| 29 | [NUDT19 — Muscle](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/NUDT19_Muscle_add_vs_mix.png) | chr19:32,850,210 | chr19:32,826,637 | 23,573 | 0.088 → 0.151 | 37 → 41 | Unavailable |
| 30 | [HLA-DQB2 — Adrenal Gland](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/HLA-DQB2_Adrenal_Gland_add_vs_mix.png) | chr6:32,636,436 | chr6:32,659,473 | 23,037 | 0.236 → 0.226 | 135 → 254 | 134 / 135 / 147 |
| 31 | [CYP2A7 — Vagina](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/CYP2A7_Vagina_add_vs_mix.png) | chr19:40,897,420 | chr19:40,876,275 | 21,145 | 0.070 → 0.030 | 38 → 68 | 37 / 38 / 47 |
| 32 | [CFHR1 — Spleen](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/CFHR1_Spleen_add_vs_mix.png) | chr1:196,856,834 | chr1:196,875,463 | 18,629 | 0.210 → 0.579 | 6 → 10 | 6 / 6 / 7 |
| 33 | [MRPS10 — Breast](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/MRPS10_Breast_add_vs_mix.png) | chr6:42,196,186 | chr6:42,211,649 | 15,463 | 0.086 → 0.026 | 54 → 95 | 53 / 54 / 64 |
| 34 | [NOMO3 — Brain](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/NOMO3_Brain_add_vs_mix.png) | chr16:16,217,459 | chr16:16,232,314 | 14,855 | 0.700 → 0.614 | 3 → 2 | Unavailable |
| 35 | [TPSB2 — Esophagus](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/TPSB2_Esophagus_add_vs_mix.png) | chr16:1,215,720 | chr16:1,230,423 | 14,703 | 0.204 → 0.375 | 9 → 15 | 9 / 9 / 10 |
| 36 | [ST7L — Thyroid](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/ST7L_Thyroid_add_vs_mix.png) | chr1:112,534,134 | chr1:112,520,503 | 13,631 | 0.079 → 0.062 | 32 → 58 | 31 / 32 / 31 |
| 37 | [NPTX1 — Skin](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/NPTX1_Skin_add_vs_mix.png) | chr17:80,681,067 | chr17:80,694,135 | 13,068 | 0.066 → 0.046 | 20 → 35 | 20 / 20 / 20 |
| 38 | [LRRC8B — Nerve](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/LRRC8B_Nerve_add_vs_mix.png) | chr1:89,566,286 | chr1:89,578,606 | 12,320 | 0.052 → 0.049 | 34 → 63 | 34 / 34 / 40 |
| 39 | [DGCR6 — Muscle](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/DGCR6_Muscle_add_vs_mix.png) | chr22:18,971,161 | chr22:18,982,890 | 11,729 | 0.346 → 0.133 | 15 → 19 | Unavailable |
| 40 | [CFHR1 — Nerve](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/CFHR1_Nerve_add_vs_mix.png) | chr1:196,868,521 | chr1:196,856,834 | 11,687 | 0.492 → 0.311 | 6 → 4 | 4 / 6 / 4 |
| 41 | [HLA-DQB2 — Prostate](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/HLA-DQB2_Prostate_add_vs_mix.png) | chr6:32,643,948 | chr6:32,654,396 | 10,448 | 0.037 → 0.126 | 137 → 179 | 133 / 137 / 133 |
| 42 | [PAQR5 — Adipose Tissue](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/PAQR5_Adipose_Tissue_add_vs_mix.png) | chr15:69,304,337 | chr15:69,314,302 | 9,965 | 0.212 → 0.179 | 39 → 46 | 26 / 39 / 27 |
| 43 | [PSD4 — Pancreas](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/PSD4_Pancreas_add_vs_mix.png) | chr2:113,183,104 | chr2:113,174,854 | 8,250 | 0.211 → 0.153 | 8 → 14 | 8 / 8 / 8 |
| 44 | [CRACR2B — Pancreas](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/CRACR2B_Pancreas_add_vs_mix.png) | chr11:817,786 | chr11:825,777 | 7,991 | 0.294 → 0.232 | 5 → 9 | Unavailable |
| 45 | [POLI — Nerve](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/POLI_Nerve_add_vs_mix.png) | chr18:54,261,622 | chr18:54,269,477 | 7,855 | 0.144 → 0.320 | 17 → 26 | Unavailable |
| 46 | [LSM7 — Pituitary](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/LSM7_Pituitary_add_vs_mix.png) | chr19:2,329,604 | chr19:2,322,322 | 7,282 | 0.204 → 0.119 | 11 → 21 | 11 / 11 / 14 |
| 47 | [LIPG — Muscle](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/LIPG_Muscle_add_vs_mix.png) | chr18:49,654,294 | chr18:49,647,298 | 6,996 | 0.078 → 0.056 | 23 → 17 | Unavailable |
| 48 | [REEP1 — Pancreas](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/REEP1_Pancreas_add_vs_mix.png) | chr2:86,342,281 | chr2:86,349,102 | 6,821 | 0.152 → 0.172 | 9 → 11 | 9 / 9 / 9 |
| 49 | [F2R — Adipose Tissue](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/F2R_Adipose_Tissue_add_vs_mix.png) | chr5:76,744,116 | chr5:76,750,006 | 5,890 | 0.143 → 0.114 | 10 → 18 | 10 / 10 / 18 |
| 50 | [L3HYPDH — Skin](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/L3HYPDH_Skin_add_vs_mix.png) | chr14:59,478,818 | chr14:59,473,009 | 5,809 | 0.095 → 0.102 | 53 → 104 | 52 / 53 / 60 |
| 51 | [TRPM4 — Spleen](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/TRPM4_Spleen_add_vs_mix.png) | chr19:49,148,533 | chr19:49,142,749 | 5,784 | 0.099 → 0.079 | 18 → 28 | Unavailable |
| 52 | [TMEM50B — Small Intestine](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/TMEM50B_Small_Intestine_add_vs_mix.png) | chr21:33,424,868 | chr21:33,430,639 | 5,771 | 0.178 → 0.015 | 63 → 73 | 56 / 63 / 68 |
| 53 | [SNX31 — Esophagus](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/SNX31_Esophagus_add_vs_mix.png) | chr8:100,661,688 | chr8:100,667,183 | 5,495 | 0.168 → 0.142 | 11 → 15 | Unavailable |
| 54 | [SLC29A3 — Testis](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/SLC29A3_Testis_add_vs_mix.png) | chr10:71,384,612 | chr10:71,379,576 | 5,036 | 0.119 → 0.093 | 13 → 20 | 13 / 13 / 14 |
| 55 | [MYRF — Liver](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/MYRF_Liver_add_vs_mix.png) | chr11:61,758,258 | chr11:61,753,846 | 4,412 | 0.867 → 0.394 | 2 → 4 | 2 / 2 / 4 |
| 56 | [WASHC1 — Pancreas](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/WASHC1_Pancreas_add_vs_mix.png) | chr9:47,303 | chr9:42,917 | 4,386 | 0.200 → 0.370 | 14 → 11 | 9 / 14 / 9 |
| 57 | [CDC16 — Adipose Tissue](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/CDC16_Adipose_Tissue_add_vs_mix.png) | chr13:114,236,830 | chr13:114,232,601 | 4,229 | 0.076 → 0.081 | 40 → 42 | 31 / 40 / 36 |
| 58 | [L3HYPDH — Prostate](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/L3HYPDH_Prostate_add_vs_mix.png) | chr14:59,428,238 | chr14:59,432,317 | 4,079 | 0.084 → 0.073 | 47 → 72 | Unavailable |
| 59 | [RASGRP3 — Pancreas](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/RASGRP3_Pancreas_add_vs_mix.png) | chr2:33,470,357 | chr2:33,474,404 | 4,047 | 0.390 → 0.702 | 3 → 5 | Unavailable |
| 60 | [RASA3 — Testis](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/RASA3_Testis_add_vs_mix.png) | chr13:114,046,441 | chr13:114,042,599 | 3,842 | 0.131 → 0.121 | 97 → 123 | 97 / 97 / 109 |
| 61 | [HLA-DQB2 — Small Intestine](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/HLA-DQB2_Small_Intestine_add_vs_mix.png) | chr6:32,659,351 | chr6:32,655,698 | 3,653 | 0.018 → 0.025 | 171 → 263 | 156 / 171 / 179 |
| 62 | [CDRT4 — Adipose Tissue](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/CDRT4_Adipose_Tissue_add_vs_mix.png) | chr17:15,540,428 | chr17:15,543,515 | 3,087 | 0.294 → 0.431 | 11 → 16 | 11 / 11 / 11 |
| 63 | [TSPEAR — Heart](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/TSPEAR_Heart_add_vs_mix.png) | chr21:44,704,546 | chr21:44,707,376 | 2,830 | 0.214 → 0.119 | 12 → 27 | Unavailable |
| 64 | [HPD — Blood](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/HPD_Blood_add_vs_mix.png) | chr12:121,856,099 | chr12:121,853,539 | 2,560 | 0.239 → 0.223 | 8 → 9 | 8 / 8 / 9 |
| 65 | [GRAMD4 — Thyroid](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/GRAMD4_Thyroid_add_vs_mix.png) | chr22:46,587,475 | chr22:46,589,890 | 2,415 | 0.315 → 0.693 | 24 → 31 | 19 / 24 / 27 |
| 66 | [RPS9 — Lung](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/RPS9_Lung_add_vs_mix.png) | chr19:54,203,694 | chr19:54,201,720 | 1,974 | 0.039 → 0.028 | 36 → 36 | 34 / 36 / 36 |
| 67 | [ETV7 — Blood](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/ETV7_Blood_add_vs_mix.png) | chr6:36,386,735 | chr6:36,388,545 | 1,810 | 0.549 → 0.304 | 23 → 17 | 15 / 23 / 15 |
| 68 | [GPA33 — Lung](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/GPA33_Lung_add_vs_mix.png) | chr1:167,090,596 | chr1:167,092,160 | 1,564 | 0.170 → 0.070 | 12 → 20 | 12 / 12 / 12 |
| 69 | [METTL7B — Blood](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/METTL7B_Blood_add_vs_mix.png) | chr12:55,676,884 | chr12:55,675,354 | 1,530 | 0.464 → 0.338 | 5 → 4 | Unavailable |
| 70 | [DNAH11 — Skin](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/DNAH11_Skin_add_vs_mix.png) | chr7:21,545,202 | chr7:21,543,901 | 1,301 | 0.095 → 0.081 | 18 → 14 | 14 / 18 / 14 |
| 71 | [SLC6A1 — Thyroid](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/SLC6A1_Thyroid_add_vs_mix.png) | chr3:10,993,660 | chr3:10,994,894 | 1,234 | 0.194 → 0.198 | 11 → 19 | 11 / 11 / 11 |
| 72 | [PSMG1 — Blood Vessel](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/PSMG1_Blood_Vessel_add_vs_mix.png) | chr21:39,307,042 | chr21:39,308,268 | 1,226 | 0.082 → 0.133 | 57 → 96 | 54 / 57 / 58 |
| 73 | [DDX59 — Heart](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/DDX59_Heart_add_vs_mix.png) | chr1:200,623,167 | chr1:200,624,340 | 1,173 | 0.145 → 0.334 | 25 → 18 | Unavailable |
| 74 | [CDKL1 — Thyroid](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/CDKL1_Thyroid_add_vs_mix.png) | chr14:50,397,441 | chr14:50,396,275 | 1,166 | 0.453 → 0.546 | 5 → 5 | 4 / 5 / 5 |
| 75 | [HTR2B — Heart](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/HTR2B_Heart_add_vs_mix.png) | chr2:231,189,862 | chr2:231,191,021 | 1,159 | 0.179 → 0.117 | 50 → 64 | Unavailable |
| 76 | [MRI1 — Muscle](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/MRI1_Muscle_add_vs_mix.png) | chr19:13,759,586 | chr19:13,760,740 | 1,154 | 0.229 → 0.202 | 5 → 10 | 5 / 5 / 6 |
| 77 | [RAET1G — Adrenal Gland](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/RAET1G_Adrenal_Gland_add_vs_mix.png) | chr6:149,888,060 | chr6:149,886,999 | 1,061 | 0.119 → 0.053 | 36 → 45 | 34 / 36 / 45 |
| 78 | [CYP4F11 — Small Intestine](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/CYP4F11_Small_Intestine_add_vs_mix.png) | chr19:15,934,939 | chr19:15,935,840 | 901 | 0.662 → 0.424 | 27 → 7 | 7 / 27 / 7 |
| 79 | [LRPAP1 — Blood Vessel](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/LRPAP1_Blood_Vessel_add_vs_mix.png) | chr4:3,491,488 | chr4:3,492,369 | 881 | 0.402 → 0.208 | 19 → 37 | 18 / 19 / 27 |
| 80 | [NUCB1 — Esophagus](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/NUCB1_Esophagus_add_vs_mix.png) | chr19:48,921,097 | chr19:48,921,668 | 571 | 0.201 → 0.164 | 6 → 11 | 6 / 6 / 6 |
| 81 | [COL18A1 — Blood](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/COL18A1_Blood_add_vs_mix.png) | chr21:45,479,536 | chr21:45,479,205 | 331 | 0.061 → 0.045 | 34 → 37 | 22 / 34 / 22 |
| 82 | [TPD52 — Adrenal Gland](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/TPD52_Adrenal_Gland_add_vs_mix.png) | chr8:80,099,793 | chr8:80,099,657 | 136 | 0.238 → 0.140 | 24 → 33 | 18 / 24 / 19 |
| 83 | [TMEM132B — Esophagus](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/TMEM132B_Esophagus_add_vs_mix.png) | chr12:125,314,635 | chr12:125,314,587 | 48 | 0.056 → 0.074 | 19 → 21 | 19 / 19 / 21 |
| 84 | [RAB29 — Muscle](C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/plot/one_cs_recessive/RAB29_Muscle_add_vs_mix.png) | chr1:205,784,684 | chr1:205,784,680 | 4 | 0.168 → 0.108 | 12 → 22 | Unavailable |

## How these complement the one-CS versus two-CS examples

The two folders address different effects of allowing non-additive coding: finding an additional conditional expression signal, changing the location or weighting of a single signal, and changing the preferred coding while retaining the same variant.

In this selected recessive folder, both models report one CS for every example. The existing biological-SNP overlap report covers 105 of the 141 examples: 26 have identical SNP sets, 79 have partial overlap, and none of these 105 have disjoint sets. Overlap is unavailable for the other 36 examples. These are selected examples, so these proportions are not genome-wide prevalence estimates.

Several examples are useful even without a long-distance lead shift:

| Example | Observation | What it can support |
|---|---|---|
| CYP4F11 — Esophagus | Same lead; PIP 0.500 → 0.901; 27 additive-CS SNPs → 2 mixed-CS biological SNPs (3 predictor columns); genotype counts 118/204/90 | A substantial concentration of localization; audit raw expression and potential floor effects before interpreting the genotype pattern. |
| ADH1C — Adipose | Same lead; PIP 0.140 → 0.424; 82 → 16 biological CS SNPs | A smaller candidate set despite unchanged CS count. |
| DNAJC15 — Blood | Same lead; PIP 0.394 → 0.288; 38 → 6 biological CS SNPs | Changed uncertainty across the locus; the leading coding-specific PIP need not increase when the CS becomes smaller. |
| METTL7B — Blood | Different leads, 1,530 bp apart; PIP 0.464 → 0.338; 39 genotype-2 individuals at mixed lead | A visibly recessive expression pattern worth follow-up, with a relatively local change in variant ranking. |
| RGS11 — Blood Vessel | Same SNP; singleton CS in both; PIP 0.994 → 1.000 | A coding interpretation example with little scope to change the prioritized physical variant. |
| VAT1L — Blood Vessel | Same SNP; singleton CS in both; PIP 0.988 → 0.995 | Another coding interpretation comparison. |

The fit statistic printed on these plots is an ELBO-plus-KL summary, not a calibrated likelihood-ratio test. Nine examples have a negative value. A recessive lead does not, by itself, establish that the mixed model fits substantially better. The paired genotype panels also reuse the same expression observations; they are not independent replications.

## Biological follow-up and colocalization

**A reproducible change in disease–expression colocalization could strengthen the manuscript**, because it would demonstrate a consequence for disease-gene or tissue prioritization. It still needs evidence that the change is better calibrated or independently supported, rather than simply different. Multiple causal signals should be compared at the signal level; this is the motivation for [SuSiE-based colocalization](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1009440).

Priorities supported by these examples and existing biology:

- **CFHR1 in nerve and spleen:** relevant to complement-mediated disease. Lead shifts are 11,687 and 18,629 bp. However, the known CFHR3/CFHR1 deletion is central to interpretation: test whether the apparent recessive SNP association tags a structural variant. A primary fine-mapping study implicates that deletion in [IgA nephropathy protection](https://pmc.ncbi.nlm.nih.gov/articles/PMC5042673/).
- **PM20D1 in stomach:** published work links the PM20D1 regulatory locus to [Alzheimer's disease biology](https://pubmed.ncbi.nlm.nih.gov/29736028/). Its 85.6-kb shift warrants checking, but the low PIPs, retained additive CS and stomach tissue do not establish a new Alzheimer mechanism.
- **LGALS3 in pituitary:** experimental work connects galectin-3 to [TREM2 and inflammatory responses in Alzheimer's disease](https://pubmed.ncbi.nlm.nih.gov/31006066/). In these fits, the lead shift is 23,643 bp but the CS sets are almost identical: 91 shared SNPs among a union of 94. This is currently a weak example of changed localization.
- **LIPG in muscle:** human loss-of-function variants affect [HDL cholesterol](https://pubmed.ncbi.nlm.nih.gov/19287092/), making lipid GWAS a biologically grounded follow-up. The current lead PIPs are diffuse; muscle expression should not automatically be treated as the disease-relevant mechanism.
- **UGT2B17:** the genotype-expression pattern needs checking against its known whole-gene deletion. Human liver data show [deletion-dependent expression and activity](https://pmc.ncbi.nlm.nih.gov/articles/PMC2993461/). A SNP tagging a deletion is biologically meaningful, but it is a different explanation from a newly established recessive regulatory SNP.

For each tested gene–tissue–GWAS combination, report the complete H0–H4 posterior distribution for both approaches, with common variant coverage and stated priors. Distinguish unresolved evidence becoming supportive of a shared variant from a reversal between strong distinct-signal evidence (H3) and strong shared-signal evidence (H4). Include both increases and decreases in colocalization support.

Colocalization must compare the same biological SNP across traits. Stacked additive/recessive/dominant columns cannot be passed as if they were distinct physical SNPs. Coding uncertainty needs appropriate marginalization with the model's priors and signal-specific evidence; simply stripping suffixes or summing plotted overall PIPs is not a Bayes-factor construction. The official [coloc.susie documentation](https://chr1swallace.github.io/coloc/reference/coloc.susie.html) describes signal-level comparisons, and [coloc.bf_bf](https://chr1swallace.github.io/coloc/reference/coloc.bf_bf.html) requires appropriately aligned SNP evidence.

Distance does not measure LD. Before interpreting a lead shift as a different causal explanation, inspect LD between the two SNPs, their posterior probabilities in both fits, credible-set membership and purity, and robustness to influential samples and genotype/structural-variant QC. Colocalization with a shared physical variant would still not prove expression mediates disease, or that the disease itself has recessive inheritance.

No GWAS colocalization, LD calculation or new model fitting was performed in this review.

