from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / 'plot/one_cs_dominant'
OUT = ROOT / 'tmp/one_cs_dominant_review'
OUT.mkdir(exist_ok=True)
d = pd.read_csv(SRC / 'fit_mix_one_cs_dominant_plot_summary.csv')
assert not d.duplicated(['gene', 'tissue']).any()
assert d.status.eq('plotted').all()
assert d.n_cs.eq(1).all() and d.lead_coding.eq('dominant').all()
for prefix, col in [('additive', 'additive_lead_snp'), ('mixed', 'lead_snp')]:
    parts = d[col].str.extract(r'^(chr[^_]+)_(\d+)_[^_]+_[^_]+_(b\d+)(?:_|$)')
    assert parts.notna().all().all(), col
    d[prefix + '_chromosome'] = parts[0]
    d[prefix + '_position_bp'] = parts[1].astype('int64')
    d[prefix + '_build'] = parts[2]
assert d.additive_chromosome.eq(d.mixed_chromosome).all()
assert d.additive_build.eq('b38').all() and d.mixed_build.eq('b38').all()
d['distance_bp'] = (d.mixed_position_bp - d.additive_position_bp).abs()
d['distance_kb'] = d.distance_bp / 1000
d['changed_lead'] = d.additive_lead_snp.ne(d.lead_snp)
assert d.changed_lead.eq(~d.same_lead_snp).all()
assert np.allclose(d.likelihood_statistic, 2 * (d.log_lik_mix - d.log_lik_add), atol=1e-8)
d['local_plot'] = [str(SRC / Path(x).name) for x in d.output_file]
assert all(Path(x).exists() for x in d.local_plot)
d['minimum_genotype_group'] = d[['n_genotype_0', 'n_genotype_1', 'n_genotype_2']].min(axis=1)
d['sample_size'] = d[['n_genotype_0', 'n_genotype_1', 'n_genotype_2']].sum(axis=1)
assert d.sample_size.eq(d[['additive_n_genotype_0', 'additive_n_genotype_1', 'additive_n_genotype_2']].sum(axis=1)).all()
d['cs_predictor_reduction_fold'] = d.additive_cs_size / d.cs_size
o = pd.read_csv(ROOT / 'descriptive_results/one_cs_overlap_and_lead_detail.csv')
keys = ['gene', 'tissue', 'additive_lead_snp', 'mixed_lead_snp']
assert not o.duplicated(keys).any()
keep = keys + ['additive_cs_snps', 'mixed_cs_snps', 'additive_cs_size_snps', 'mixed_cs_size_snps', 'shared_snp_count', 'union_snp_count', 'jaccard_similarity', 'overlap_class']
d = d.merge(o[keep].rename(columns={'mixed_lead_snp': 'lead_snp'}), on=['gene', 'tissue', 'additive_lead_snp', 'lead_snp'], how='left', validate='one_to_one')
d['overlap_available'] = d.shared_snp_count.notna()
for r in d[d.overlap_available].itertuples():
    a, b = set(r.additive_cs_snps.split(';')), set(r.mixed_cs_snps.split(';'))
    assert len(a) == r.additive_cs_size_snps == r.additive_cs_size
    assert len(b) == r.mixed_cs_size_snps <= r.cs_size
    assert len(a & b) == r.shared_snp_count
    assert len(a | b) == r.union_snp_count
    assert abs(len(a & b) / len(a | b) - r.jaccard_similarity) < 1e-12
d['cs_biological_snp_reduction_fold'] = d.additive_cs_size_snps / d.mixed_cs_size_snps
d['mixed_lead_in_additive_cs'] = pd.Series(pd.NA, index=d.index, dtype='boolean')
d['additive_lead_in_mixed_cs'] = pd.Series(pd.NA, index=d.index, dtype='boolean')
for i, r in d[d.overlap_available].iterrows():
    d.loc[i, 'mixed_lead_in_additive_cs'] = r.lead_snp in r.additive_cs_snps.split(';')
    d.loc[i, 'additive_lead_in_mixed_cs'] = r.additive_lead_snp in r.mixed_cs_snps.split(';')
d = d.sort_values(['distance_bp', 'gene', 'tissue'], ascending=[False, True, True])
d['distance_rank_all_changed'] = pd.Series(range(1, d.changed_lead.sum() + 1), index=d.index[d.changed_lead], dtype='Int64')
one = d.changed_lead & d.additive_n_cs.eq(1)
d['distance_rank_both_one_cs'] = pd.Series(range(1, one.sum()+1), index=d.index[one], dtype='Int64')
d.to_csv(OUT / 'all_candidates.csv', index=False)
stats = {
    'total': len(d), 'changed': int(d.changed_lead.sum()), 'same': int((~d.changed_lead).sum()),
    'additive_cs_counts': d.additive_n_cs.value_counts().to_dict(),
    'changed_both_one_cs': int(one.sum()),
    'changed_additive_two_cs': int((d.changed_lead & d.additive_n_cs.eq(2)).sum()),
    'distance_gt_100kb': int((d.distance_bp > 100000).sum()),
    'distance_gt_500kb': int((d.distance_bp > 500000).sum()),
    'distance_gt_1mb': int((d.distance_bp > 1000000).sum()),
    'overlap_counts': d.overlap_class.fillna('Unavailable').value_counts().to_dict(),
    'negative_fit_metric': int((d.likelihood_statistic < 0).sum()),
    'changed_distance_quantiles_bp': d.loc[d.changed_lead, 'distance_bp'].quantile([.25,.5,.75,.9]).to_dict(),
    'changed_same_position': int((d.changed_lead & d.distance_bp.eq(0)).sum()),
}
(OUT / 'statistics.json').write_text(json.dumps(stats, indent=2))
print(json.dumps(stats, indent=2))
cols = ['gene','tissue','distance_bp','additive_n_cs','additive_lead_pip','lead_pip','additive_cs_size','cs_size','likelihood_statistic','n_genotype_0','n_genotype_1','n_genotype_2','shared_snp_count','additive_cs_size_snps','mixed_cs_size_snps']
groups = {
    'TOP DISTANCE': d.head(20),
    'CHANGED PIP >= .5': d[d.changed_lead & d.lead_pip.ge(.5)].head(20),
    'FIT METRIC, PIP >= .5 AND ALL GROUPS >=20': d[d.lead_pip.ge(.5) & d.minimum_genotype_group.ge(20)].nlargest(15,'likelihood_statistic'),
    'BIOLOGICAL CS SHRINKAGE, PIP >= .3 AND ALL GROUPS >=20': d[d.lead_pip.ge(.3) & d.minimum_genotype_group.ge(20)].nlargest(15,'cs_biological_snp_reduction_fold'),
    'CHANGED BOTH PIPS >=.5 AND ALL GROUPS >=20': d[d.changed_lead & d.lead_pip.ge(.5) & d.additive_lead_pip.ge(.5) & d.minimum_genotype_group.ge(20)].head(15),
}
for label, frame in groups.items():
    print('\n' + label)
    print(frame[cols].to_string(index=False, float_format=lambda x: f'{x:.3f}'))

# Curated after inspecting the original plots. Reasons distinguish physical
# lead shifts, posterior concentration, and coding interpretation.
notes = [
    ('TSLP', 'Esophagus', 'Priority: long-distance lead change',
     'Largest changed-lead distance among examples with mixed lead PIP >= 0.5. The mixed-lead genotype-1 and genotype-2 means are very similar and both exceed genotype 0. All groups are well represented. CS predictor count increases from 20 to 29, so this is a lead-prioritization example, not demonstrated CS contraction. Biological-SNP overlap is unavailable.'),
    ('CMTM6', 'Thyroid', 'Priority: long-distance lead change',
     'A 113.9-kb shift, mixed lead PIP 0.614, and a clear visual plateau across genotypes 1 and 2. Genotype counts are 107/196/87. The plotted CS decreases slightly from 15 to 13 predictors, but unique mixed-CS SNP count and overlap are unavailable.'),
    ('TPSD1', 'Small Intestine', 'Priority: changed lead and concentrated posterior',
     'The lead shifts 1,992 bp and the mixed fit concentrates on one dominant predictor with PIP 0.999. Four additive-CS biological SNPs become one. The selected mixed SNP was already in the additive CS. All groups have at least 26 individuals. Genotype 0 is strongly separated from carriers, although genotype 2 remains somewhat higher than genotype 1.'),
    ('APOBEC3D', 'Nerve', 'Priority: CS contraction at the same lead',
     'The lead stays at chr22:39,024,088. Its PIP changes from 0.258 to 0.922 and the CS contracts from 15 to 2 unique biological SNPs (3 mixed predictor columns). Genotype-1 and genotype-2 means are very similar. This is a strong visual example of reduced localization uncertainty without a lead change.'),
    ('KIAA0930', 'Muscle', 'Priority: coding interpretation at a fixed SNP',
     'Both fits have the same singleton CS and lead PIP approximately 1. The expression shift from genotype 0 to 1 is much larger than the additional shift from 1 to 2. The fit metric is +26.26 and genotype counts are 231/212/34. Useful for coding interpretation, with no change in the prioritized physical SNP.'),
    ('STAP2', 'Nerve', 'Priority: dominant expression pattern and fit metric',
     'Same lead, CS contraction from 5 to 3 biological SNPs, and fit metric +30.52. Genotypes 1 and 2 have similar elevated expression relative to genotype 0. Mixed lead PIP is moderate at 0.573. This is the largest fit metric among examples with mixed PIP >= 0.5 and at least 20 individuals in every genotype group.'),
    ('EIF3C', 'Blood Vessel', 'Priority: maximum physical shift',
     'The largest shift is 480,493 bp, from chr16:28,481,454 to chr16:28,961,947. Lead PIPs are 0.079 and 0.302, with 56 versus 28 CS predictors and fit metric +7.78. Carrier means are similar at the mixed lead. This is worth showing for the physical shift, but localization is less concentrated than in TSLP or CMTM6 and biological-SNP overlap is unavailable.'),
    ('CPXM1', 'Spleen', 'Additional: two confident but different leads',
     'The leads are 4,025 bp apart with own-lead PIPs 0.971 and 0.854. The mixed-lead plot has similar genotype-1 and genotype-2 means and counts 48/54/44. The mixed CS has two SNPs and retains the single additive-CS SNP. This is a useful local change in prioritization, not a disjoint candidate set. The additive lead has only 9 genotype-2 individuals.'),
    ('NOL6', 'Esophagus', 'Additional: same-lead posterior concentration',
     'Same lead, PIP 0.631 to 0.992, and four CS SNPs become one. Fit metric +12.36. Genotype 2 has 22 individuals and its mean is lower than the genotype-1 mean, so the plot is approximately compatible with a carrier effect rather than a perfect plateau.'),
    ('SPINK5', 'Heart', 'Additional: changed lead and smaller candidate set',
     'Lead shift 33,930 bp, CS contraction from 12 to 5 biological SNPs with only 2 shared SNPs, and fit metric +21.78. Genotype counts are 110/154/67. Both lead PIPs remain modest (0.169 and 0.313), and genotype 2 remains higher than genotype 1. Useful as a change in the candidate set with residual uncertainty.'),
    ('CRYM', 'Muscle', 'Additional: strong fit metric with a smaller genotype-2 group',
     'Same lead and the same two biological CS SNPs. PIP 0.746 to 0.882 and fit metric +55.44, the second largest in this folder. The genotype-1 and genotype-2 means are similar, but genotype 2 has 18 individuals and there are repeated values near the lower expression tail. Inspect original counts and influential observations before choosing it as the main expression-pattern example.'),
]
recommendations = []
for gene, tissue, purpose, note in notes:
    row = d[d.gene.eq(gene) & d.tissue.eq(tissue)].iloc[0].to_dict()
    row.update(example_priority=len(recommendations)+1, example_purpose=purpose, example_notes=note)
    recommendations.append(row)
recommended = pd.DataFrame(recommendations)
csv_cols = ['distance_rank_all_changed', 'distance_rank_both_one_cs', 'gene', 'tissue',
            'changed_lead', 'additive_n_cs', 'n_cs', 'additive_lead_snp', 'lead_snp',
            'additive_chromosome', 'additive_position_bp', 'mixed_position_bp', 'distance_bp', 'distance_kb',
            'additive_lead_pip', 'lead_pip', 'additive_cs_size', 'cs_size',
            'overlap_available', 'additive_cs_size_snps', 'mixed_cs_size_snps',
            'shared_snp_count', 'union_snp_count', 'jaccard_similarity', 'overlap_class',
            'mixed_lead_in_additive_cs', 'additive_lead_in_mixed_cs',
            'likelihood_statistic', 'additive_n_genotype_0', 'additive_n_genotype_1', 'additive_n_genotype_2',
            'n_genotype_0', 'n_genotype_1', 'n_genotype_2', 'local_plot']
d.loc[d.changed_lead, csv_cols].to_csv(SRC / 'lead_snp_distance_ranking.csv', index=False)
recommended[['example_priority', 'example_purpose', 'example_notes']+csv_cols].to_csv(SRC / 'recommended_examples.csv', index=False)

def link(r):
    return f'[{r.gene} — {r.tissue}](<{Path(r.local_plot).as_posix()}>)'

def coord(r, prefix):
    return f'{r[prefix + "_chromosome"]}:{int(r[prefix + "_position_bp"]):,}'

def overlap(r):
    if not r.overlap_available:
        return 'Unavailable'
    return f'{int(r.shared_snp_count)} / {int(r.additive_cs_size_snps)} / {int(r.mixed_cs_size_snps)}'

def distance_table(frame):
    lines = ['| Rank | Gene and tissue / plot | Additive lead (GRCh38) | Mixed dominant lead (GRCh38) | Distance (bp) | Own-lead PIPs, additive → mixed | CS predictors, additive → mixed | Additive CS count | Shared / additive / mixed CS SNPs |',
             '|---:|---|---|---|---:|---|---|---:|---|']
    for _, r in frame.iterrows():
        lines.append(f'| {int(r.distance_rank_all_changed)} | {link(r)} | {coord(r,"additive")} | {coord(r,"mixed")} | {r.distance_bp:,} | {r.additive_lead_pip:.3f} → {r.lead_pip:.3f} | {int(r.additive_cs_size)} → {int(r.cs_size)} | {r.additive_n_cs} | {overlap(r)} |')
    return '\n'.join(lines)

source_link = f'[saved plot summary](<{(SRC / "fit_mix_one_cs_dominant_plot_summary.csv").as_posix()}>)'
overlap_link = f'[biological-SNP overlap report](<{(ROOT / "descriptive_results/one_cs_overlap_and_lead_detail.csv").as_posix()}>)'
text = [
    '# Dominant one-CS examples: lead-SNP distance ranking and recommended plots',
    'Reviewed 2026-09-09. Scope: the 1,340 plotted gene–tissue examples in `plot/one_cs_dominant`, comparing additive SuSiE with the saved unweighted SuSiE-mix fit.',
    '**Main findings**',
    '**466 examples have different lead SNPs; 874 have the same lead. The maximum separation is 480,493 bp for EIF3C in blood vessel. Forty-six changed leads are more than 100 kb apart.** All compared pairs are on the same chromosome and use the b38 coordinates encoded in the identifiers. The median separation among changed leads is 15,145 bp.',
    'For a long-distance example with appreciable lead support and well represented genotype groups, start with **TSLP in esophagus** and **CMTM6 in thyroid**. For posterior concentration, start with **TPSD1 in small intestine** and **APOBEC3D in nerve**. For a coding interpretation example at a fixed SNP, use **KIAA0930 in muscle** or **STAP2 in nerve**.',
    '**Scope and measurement**',
    f'Source: {source_link}. Each row is a gene–tissue pair. Distance is `abs(mixed lead base-pair position - additive lead base-pair position)`. The x-axis in the original PIP plots indexes stacked predictors; its apparent separation is not the physical distance. No liftover was performed.',
    'All 1,340 mixed fits have one CS with a dominant lead. Additive SuSiE has one CS in 1,307 examples and two CSs in 33. Of the changed-lead cases, 450 have one CS in both models and 16 have two additive CSs. For a two-CS additive fit, the saved comparison uses the highest-PIP CS lead; it is not a matched signal-level comparison. The first 16 entries in the overall distance ranking all have one CS in both models. Within the strict one-CS-versus-one-CS subset, 44 lead shifts exceed 100 kb. Both rank columns are saved in the CSV.',
    'The selection script applies minimum association p-value < 5e-8, mean count >= 100, one mixed CS, and a dominant lead in that CS. This review uses the saved plotted cohort and does not re-estimate those selection quantities. Counts here describe this selected folder, not genome-wide prevalence.',
    'PIPs refer to each model\'s own saved lead predictor. The mixed PIP is coding-specific, not marginalized over all codings of a biological SNP. When leads differ, the two PIPs are for different SNPs. Plotted CS sizes count predictors; multiple codings can represent the same biological SNP.',
    f'CS overlap was joined from the {overlap_link} on gene, tissue, and both exact lead-SNP identifiers. Available sets and counts were cross-checked. Missing matches remain unavailable. Among all 1,340 examples, 1,045 have matched overlap: 295 identical sets and 750 partially overlapping sets, with no disjoint sets among these matches. Overlap is unavailable for 295 examples. Among the 466 changed-lead cases, 47 have identical sets, 296 partial overlap, and 123 unavailable overlap.',
    'The plotted fit statistic is `2 * [(final ELBO + sum(KL))_mix - (final ELBO + sum(KL))_add]`, as documented in the local plotting utility. It is a descriptive fit metric without a calibrated likelihood-ratio-test p-value. There are 116 negative values in this folder. Positive values should not be described as statistical significance.',
    '**The ten largest physical shifts**',
    distance_table(d.head(10)),
    'EIF3C is the strongest candidate to inspect among the top five distances because it combines a positive fit metric (+7.78), increased own-lead PIP (0.079 to 0.302), and 38 genotype-2 individuals at the mixed lead. Its biological-SNP CS overlap is unavailable, so a disjoint candidate set has not been established.',
    'The next-largest shifts are less compelling localization examples. LRRC37A2 in ovary has mixed lead PIP 0.009, 900 biological SNPs in the mixed CS, and only four genotype-2 individuals. AS3MT in esophagus has PIP 0.024 and 142 shared SNPs across CSs of 164 and 153 biological SNPs. CLP1 in testis has PIP 0.067 and only five genotype-2 individuals. NDUFA6 in blood vessel keeps exactly the same two biological CS SNPs despite a 296,986-bp lead shift; its mixed lead PIP is 0.425, its fit metric is +1.15, and genotype 2 has seven individuals.',
    '**Recommended examples beyond distance alone**',
    'The shortlist weighs the physical lead shift, PIP concentration, biological-SNP CS contraction where available, genotype-group sample sizes, the descriptive fit metric, and the original expression plots. A screening preference of mixed lead PIP >= 0.5 and at least 20 individuals per genotype group helped identify clear examples; these are practical review cutoffs, not significance or power thresholds. EIF3C and the additional candidates retain their explicit limitations.',
    '| Example / plot | Distance (bp) | Own-lead PIPs, additive → mixed | CS predictors, additive → mixed | Shared / additive / mixed CS SNPs | Fit metric | Mixed-lead genotype counts 0 / 1 / 2 |',
    '|---|---:|---|---|---|---:|---|',
]
for _, r in recommended.iterrows():
    text.append(f'| {link(r)} | {r.distance_bp:,} | {r.additive_lead_pip:.3f} → {r.lead_pip:.3f} | {int(r.additive_cs_size)} → {int(r.cs_size)} | {overlap(r)} | {r.likelihood_statistic:+.2f} | {r.n_genotype_0} / {r.n_genotype_1} / {r.n_genotype_2} |')
text.append('')
for _, r in recommended.iterrows():
    text.append(f'- **{r.gene} — {r.tissue}: {r.example_purpose}.** {r.example_notes}')
text += [
    '',
    'Only three changed-lead cases exceed 100 kb and have mixed lead PIP >= 0.5: TSLP in esophagus (143,395 bp), MAN2C1 in uterus (121,542 bp), and CMTM6 in thyroid (113,862 bp). MAN2C1 has only five genotype-2 individuals, so it is a useful secondary candidate with less information for distinguishing the genotype-1 and genotype-2 expression levels.',
    'TPSD1 also appears in esophagus and nerve with higher mixed lead PIPs (0.732 and 0.970), positive fit metrics (+24.40 and +12.52), and CS contraction (4 to 2 and 3 to 1 biological SNPs). These two tissues share the lead chr16:1,244,024, while small intestine selects chr16:1,245,817. The plots show some residual increase from genotype 1 to 2. This is a useful cross-tissue follow-up, not evidence of the exact same causal SNP or independent replication.',
    'MDGA1 in blood has PIP 0.998, a singleton mixed CS and fit metric +24.64, but genotype 2 remains visibly higher than genotype 1. PSG4 in skin has the largest fit metric (+56.64), yet mixed lead PIP is only 0.079. These illustrate why a fit-metric ranking alone does not identify the clearest dominant plateau or the best localized SNP.',
    '**What would strengthen these examples**',
    'Measure LD between each lead pair and inspect both SNPs in both posterior distributions, preferably with coding-marginalized biological-SNP evidence as well as coding-specific PIPs. Compare credible-set membership and purity. For EIF3C, TSLP and CMTM6, recover the missing biological-SNP CS overlap before interpreting their long shifts as changes to the candidate set. Distance alone does not establish low LD or a different causal signal.',
    'At each proposed dominant lead, quantify the two contrasts between genotype-group means (1 minus 0, and 2 minus 1), with uncertainty and sample counts. A dominant-shaped plot should have a sizeable carrier contrast and comparatively little additional change from genotype 1 to 2. The summaries do not contain individual expression values, so this review makes visual assessments rather than estimating those contrasts or declaring equivalence of genotype means.',
    'The two expression panels reuse the same normalized expression observations. Their differences arise from genotype grouping. Original counts, genotype quality and influential observations should be checked in follow-up. No new model fitting, LD calculation, formal genotype-contrast test or GWAS colocalization was performed.',
    '**Files and verification**',
    f'- [All 466 changed leads with both rank columns](<{(SRC / "lead_snp_distance_ranking.csv").as_posix()}>).',
    f'- [Annotated shortlist](<{(SRC / "recommended_examples.csv").as_posix()}>).',
    f'- [Distance and lead-PIP overview](<{(SRC / "lead_snp_distance_overview.png").as_posix()}>).',
    'Checked all 1,340 rows for unique gene–tissue keys, plotted status, one dominant mixed CS, parseable and matched b38 chromosome coordinates, agreement with the saved same-lead flag, equal sample totals across panels, correct fit-statistic arithmetic, and existing local plot files. Matched overlap records were checked against their listed SNP sets and CS sizes. Twenty original plots were visually reviewed: the top five distance cases, all eleven shortlisted cases, MAN2C1 in uterus, MDGA1 in blood, PSG4 in skin, and TPSD1 in esophagus and nerve (with EIF3C counted once across these groups).',
    '**All changed leads, ordered by physical distance**',
    distance_table(d[d.changed_lead]),
]
(SRC / 'lead_snp_distance_ranking.md').write_text('\n\n'.join(text).replace('|\n\n|', '|\n|') + '\n', encoding='utf-8')

# Use the installed R graphics runtime for the scientific overview; the
# bundled analysis Python does not include matplotlib.
import subprocess
subprocess.run([r'C:\Program Files\R\R-4.6.1\bin\Rscript.exe',
                str(ROOT / 'tmp/plot_dominant_distance_overview.R'), str(ROOT)], check=True)
print('\nSaved report, changed-lead ranking, annotated examples, and overview figure in', SRC)
