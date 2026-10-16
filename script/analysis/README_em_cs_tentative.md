# Tentative EM four-panel comparison

Source `plot_cs_comparison_em_tentative.R` in RStudio after editing its settings,
or run from the project root on the cluster:

```sh
Rscript --vanilla script/analysis/plot_cs_comparison_em_tentative.R /project2/mstephens/wdenault/susie_mix 18
```

This is a separate figure: it does not modify `plot_cs_agreement_em.R`, its
settings, or its outputs. Only base R is needed. It compares original additive
SuSiE with `weighted_fit_mix` from iteration 18, which is pinned by default.
No model fits, simulation runs or expression reconstruction are performed.

Inputs are the completed iteration and provenance-checked region and CS summaries
from `generate_summary_results_em.R`. The CS summary must contain `cs_snps` and
`cs_size_snps`. The existing loader checks completion markers, result-file
availability, iteration identity and prior provenance. The old figure's CSV audits
alone cannot reconstruct the CS membership required for the new panels.

The cohort uses P < 1e-8, mean reads >= 100, and at least one CS in either model.
Each region is a gene/tissue pair. This is the descriptive-analysis cohort, not
the more permissive selected one-CS example cohort. Counts include regions with
no CS in one model. Both-zero regions are excluded and audited.

- **A:** exact CS-count tables for the two tissues with lowest agreement and the
  two with highest agreement, requiring at least 100 regions per tissue. Ranking
  uses equal CS-count proportions; ties use larger sample size, then tissue name.
  A tissue cannot be featured twice. Fewer qualifying tissues produce fewer tables.
  All tissue rankings and descriptive Wilson intervals are exported.
- **B:** ECDFs for changed paired leads, unmatched EM leads to their nearest SuSiE
  lead, and unmatched SuSiE leads to their nearest EM lead. Identical biological
  lead IDs are paired first, ignoring coding. Remaining leads are paired one to
  one to minimize total physical distance. The larger model's leftover CSs are
  unmatched. Exact matches are counted separately; different alleles at the same
  coordinate retain their zero distance. Extras in model-only regions have no
  distance to the other model and are counted separately, not assigned zero.
- **C:** physical distance versus CS Jaccard overlap among changed paired leads,
  coloured by EM lead coding. Jaccard = intersection / union of biological SNP
  sets; 0 means disjoint and 1 means identical. Missing/invalid membership is
  audited and omitted. This panel does **not** estimate LD. Pairing by distance is
  descriptive and should not be interpreted as validated signal correspondence.
- **D:** two small views of resolution and support. Left: EM/SuSiE unique-SNP CS
  size ratios for all matched pairs, grouped by EM lead coding; below 1 means a
  smaller EM set. Boxes show median/IQR, whiskers extend to 1.5 IQR, and plus signs
  mark the full range. Right: unmatched EM CS size versus coding-specific lead PIP,
  including model-only regions. A diffuse CS can have a small lead PIP; lead PIP
  is not the posterior probability that the entire signal exists. Counts at PIP
  >= 0.5 and >= 0.9 are exported without imposing a display filter.

CS members are deduplicated biological SNP IDs, collapsing genetic codings.
Membership must agree with the declared unique-SNP size and contain the lead.
Invalid CS membership does not remove a region from the count comparison.
Unusable coordinates/incomplete CS summaries are audited by the shared matching
helper and omitted from matching-dependent panels.

Neither smaller CSs nor additional CSs establish improved accuracy without truth
or validation. Matched pairs and extra CSs are observation units for B-D; these
are not equally weighted per-region distributions. Multiple tissues can share
genes and donors. Repetition across tissues is not independent replication.

Outputs are saved separately under:

```text
results_em/iteration_018/descriptive_results/cs_comparison_tentative/
  em_cs_comparison_tentative_4panel.png
  em_cs_comparison_tentative_4panel.pdf
  tentative_plot_data.rds
  selected_tissues.csv
  tissue_ranking.csv
  tissue_cs_count_tables.csv
  matched_cs_pairs.csv
  unpaired_cs.csv
  resolution_summary.csv
  extra_em_support_summary.csv
  lead_distance_summary.csv
  selection_summary.csv
  region_audit.csv
  cs_membership_audit.csv
  analysis_settings.csv
```

To redraw locally without copying all original fits, copy `tentative_plot_data.rds`
from the cluster and set `plot_data_file` to its path in the entry script. Exported
cohort settings are retained in this mode. Synthetic validation exports retain
their synthetic label when replotted.

Sync the new entry and helper together with existing `em_analysis_utils.R`,
`em_cs_agreement_utils.R`, and `script/scan_tissue_attempt/em_utils.R`.

Validation (no real data):

```sh
Rscript --vanilla script/analysis/tests/test_em_cs_tentative.R
```

Set `EM_TENTATIVE_TEST_OUTPUT` to retain the labelled synthetic previews in a
chosen directory. Synthetic output is never written to the production results
folder by the test.
