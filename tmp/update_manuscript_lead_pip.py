"""Stage a targeted manuscript update from the current source and checked results."""
from pathlib import Path
import copy
import csv
import hashlib
import json
import re
import shutil

ROOT = Path(__file__).resolve().parents[1]
TARGET = Path('C:/Document/Serieux/Travail/Package/git/Science-SuSiE-mix')
STAGE = ROOT / 'tmp/manuscript_lead_pip_update'
STAGE.mkdir(exist_ok=True)
raw = (TARGET / 'susie_mix_manuscript.tex').read_bytes()
tex = raw.decode('utf-8-sig').replace('\r\n', '\n')
data = json.loads((TARGET / 'simulation_sections.json').read_text(encoding='utf-8-sig'))
old_data = copy.deepcopy(data)
for folder in ('simulation_figures', 'simulation_tables'):
    for src in (TARGET / folder).iterdir():
        if src.is_file(): shutil.copyfile(src, STAGE / folder / src.name)

def paragraph(s): return {'kind': 'paragraph', 'text': s.strip()}
def heading(s): return {'kind': 'heading', 'text': s}
def replace_text(blocks, old, new):
    found = 0
    for b in blocks:
        if old in b.get('text', ''):
            b['text'] = b['text'].replace(old, new)
            found += 1
    assert found == 1, (old, found)

results = data['blocks']['results']
replace_text(results, 'Evaluation used biological SNPs, combining\nthe retained encodings of each SNP in the mixed-model output (Methods).',
    'SNP-localization analyses combined the retained encodings of each SNP\nin the mixed-model output. We additionally evaluated exact identification\nof the generating SNP--coding pair using the global lead PIP (Methods).')
replace_text(results, 'The\nfigures report pooled point estimates across these shared replication\nblocks.',
    'The figures report pooled point estimates; uncertainty intervals for\nthe coding PIP-mass shares resample these shared replication blocks\n(Methods).')

replace_text(results, 'ten single-effect components. SNP-localization', 'ten single-effect components.\nSNP-localization')
replace_text(results, 'results are dependent. The figures', 'results are dependent.\nThe figures')
rows = list(csv.DictReader((STAGE / 'simulation_tables/lead_pip_high_confidence.csv').open()))
by_code = {r['called_coding']: r for r in rows}
table_rows = []
for cl in ('additive', 'recessive', 'dominant'):
    r = by_code[cl]
    table_rows.append([cl.capitalize(), f"{int(r['n_leads']):,}", f"{int(r['n_exact']):,}",
                       f"{100*float(r['mean_pip']):.1f}", f"{100*float(r['exact_accuracy']):.1f}"])
table = {'kind': 'table', 'label': 'tab:sim-lead-accuracy', 'number': '1',
    'title': 'Exact discovery among high-PIP global leads.',
    'text': r'''One positive, unambiguous global lead was selected per dataset;
only leads with coding-specific PIP $\geq 0.9$ are included. An exact
discovery matches both a generating causal SNP and its generating coding.
Counts pool all causal allocations, PVE values, and causal counts with
one observation per selected dataset. Mean PIP and exact-match frequency
are percentages. These rates measure precision of the selected leads,
not recovery of all causal effects.''',
    'columns': ['Called coding', 'Leads', 'Exact matches', 'Mean PIP (%)', 'Exact (%)'], 'rows': table_rows}
results.extend([
heading('High-PIP leads identify generating SNP--coding pairs'),
paragraph(r'''To assess whether non-additive labels identify the generating
coding, we selected the single largest coding-specific PIP across each
simulated locus and evaluated the exact SNP--coding pair. Among
unambiguous leads with PIP $\geq 0.9$, 6,521 of 6,805 dominant leads
(95.8\%) and 11,361 of 11,570 recessive leads (98.2\%) matched a generating
pair; the corresponding rate for additive leads was 5,458 of 6,063
(90.0\%; Table~\ref{tab:sim-lead-accuracy}). Errors at a causal SNP with
the wrong coding accounted for 92 dominant leads, six recessive leads,
and 100 additive leads; the remaining errors selected noncausal SNPs.'''),
table,
paragraph(r'''Exact-match frequency generally increased with lead PIP
(Fig.~\ref{fig:sim-lead-accuracy}). However, numerical PIPs were not
uniformly calibrated on this simulation grid. Mean PIPs among the
high-PIP additive, dominant, and recessive leads were 97.7\%, 97.8\%,
and 98.9\%, respectively, compared with exact-match frequencies of
90.0\%, 95.8\%, and 98.2\%. Performance also depended on signal strength
and causal count: at PVE of 10\% with five causal effects, only 12 of
26 high-PIP dominant leads were exact matches; all 14 errors selected
noncausal SNPs. Condition-specific counts are supplied with the
simulation tables. Thus, the pooled accuracy does not imply uniformly
reliable inference in weak, multicausal loci.'''),
heading('Additive controls distinguish confident calls from diffuse coding support'),
paragraph(r'''The 7,996 additive-only datasets provided a control for false
non-additive lead calls. Only 27 datasets (0.34\%) produced an unambiguous
dominant global lead with PIP $\geq 0.9$, and three (0.038\%) produced
a recessive global lead at this threshold. Without a PIP threshold,
dominant coding accounted for 11.9--27.9\% of unambiguous leads across
causal counts at PVE of 10\%, but 0.6--5.8\% at PVE of 40\%
(Fig.~\ref{fig:sim-additive-lead}). These per-dataset control rates
support the specificity of confident non-additive lead calls within
the evaluated design, while showing that the highest-ranked coding
alone is less informative when posterior support is weak.'''),
paragraph(r'''Confident lead calls and pooled posterior mass nevertheless
gave different summaries. Even under wholly additive generating models,
the sum of dominant PIPs divided by the sum of PIPs over all codings
ranged from 18.0\% to 35.0\% across PVE and causal-count conditions;
the corresponding recessive shares ranged from 6.6\% to 24.9\%
(Fig.~\ref{fig:sim-additive-pip-mass}). These nonzero shares show that
pooled coding PIP mass cannot directly be interpreted as the fraction
of generating effects following each coding. The lead results provide
support for investigating high-confidence non-additive candidates in
the main analysis, while the mass-allocation controls limit conclusions
about their genome-wide prevalence.''')])

discussion = data['blocks']['discussion']
discussion.insert(2, paragraph(r'''The coding-specific evaluation adds evidence
that high-PIP dominant and recessive leads can identify the generating
SNP and coding. High-PIP non-additive lead calls were infrequent
under additive truth, providing a relevant control for the non-additive
candidates in the GTEx analysis. This evidence concerns the unique
global lead of an unweighted simulated fit; the GTEx summaries label
individual credible sets by their lead predictor, so these numerical
accuracy rates cannot be transferred directly to every reported set.
Moreover, appreciable dominant and recessive PIP mass under additive
truth means that coding-mass proportions are not estimates of biological
effect prevalence without further validation. The present benchmark
does not assess the weighted or iteratively reweighted fits.'''))
replace_text(discussion,
    'The point-estimate\ncomparisons also require uncertainty assessment that preserves shared\nseeds across methods and conditions. ROC and power--FDR curves assess\ndiscrimination; they do not establish that the numerical PIPs equal\nempirical causal probabilities. Finally, the simulated outcome was not',
    'The point-estimate comparisons of coverage, recovery, and lead\naccuracy also require uncertainty assessment that preserves shared\nseeds across methods and conditions. ROC and power--FDR curves assess\ndiscrimination; the additional lead-PIP reliability analysis evaluates\nexact SNP--coding matches and reveals remaining overconfidence,\nparticularly for additive leads and some weak-signal settings. These\ndescriptive checks do not establish formal false discovery control.\nFinally, the simulated outcome was not')

replace_text(discussion, 'this issue. The point-estimate comparisons', 'this issue.\nThe point-estimate comparisons')
methods = data['blocks']['methods']
replace_text(methods,
    'The figures included here display pooled point estimates without\nMonte Carlo uncertainty intervals. The shared seeds define replication\nblocks across conditions and must be retained together when assessing\nuncertainty or comparing methods. The current figure set does not\ninclude a credible-set-size comparison or PIP reliability diagram.',
    'Coverage, recovery, ROC, power--FDR, and lead-PIP figures display\npooled point estimates without Monte Carlo uncertainty intervals.\nThe coding PIP-mass figure additionally reports seed-block bootstrap\nintervals, as described below. Shared seeds define replication blocks\nacross conditions and must be retained together in uncertainty analyses.\nThe current figure set does not include a credible-set-size comparison.')
methods.extend([
heading('Exact discovery and reliability of the global lead PIP'),
paragraph(r'''We used the saved coding-specific PIPs from the unweighted
SuSiE-mix fit to select the largest PIP across all retained SNP--coding
pairs in each dataset. Selection did not use generating truth, credible
sets, or biological-SNP PIPs. A positive lead was unambiguous if no
other predictor was within an absolute PIP difference of $10^{-12}$
of the maximum. Multiple such predictors were recorded as tied leads,
including ties across distinct SNPs or between true effects; they were
not resolved using column order or truth. Datasets with all PIPs zero
were recorded separately. There were 15,547 tied datasets and 71 with
zero support among the 87,988 usable datasets; only 20 tied datasets
had maximum PIP $\geq 0.9$.'''),
paragraph(r'''An exact discovery required the selected predictor to match
both the SNP and coding of any generating effect. Errors were divided
into a wrong coding at a causal SNP and selection of a noncausal SNP.
For multiple causal effects, this evaluates the strongest reported hit;
it does not measure recovery of every effect. For positive, unambiguous
leads, we grouped maximum PIPs into ten bins, $[0,0.1),\ldots,[0.9,1]$,
and compared mean PIP with the fraction of exact matches by called
coding and PVE. Figure~\ref{fig:sim-lead-accuracy} pools all causal
allocations and counts within PVE, with one observation per selected
dataset; tables also retain PVE and causal count separately. Empty
bins have no estimated accuracy. Table~\ref{tab:sim-lead-accuracy}
restricts these leads to PIP $\geq 0.9$ and pools the full grid.'''),
paragraph(r'''In additive-only datasets, we calculated the frequency of
each called coding among positive, unambiguous leads. The dominant
and recessive high-PIP control rates instead used all 7,996 usable
additive-only datasets as the denominator, counting whether a dataset
produced an unambiguous non-additive global lead with PIP $\geq 0.9$.
These are rates of false non-additive lead calls per simulated dataset,
not false discovery proportions among non-additive calls in real data.
The checked, deduplicated checkpoint collection was the same as for
the SNP-localization analysis. Unequal retained coding blocks were
reconstructed from saved predictor-to-SNP maps and checked against
the generating indices. PIP excursions outside $[0,1]$ no larger than
$10^{-10}$ were clipped for numerical roundoff; this affected 24
datasets, including roundoff in biological-SNP PIPs. No dataset was
excluded because of an invalid coding map or PIP.'''),
heading('Coding PIP-mass allocation under additive truth'),
paragraph(r'''For each additive-only PVE and causal-count condition, we
pooled the raw coding-specific PIPs over all retained predictors and
usable datasets. For coding $c$, the mass share was'''),
{'kind': 'equation', 'text': r'''\widehat{\pi}_c =
\frac{\sum_{r}\sum_{j:(j,c)\in\mathcal{P}_r}\mathrm{PIP}_{rjc}}
     {\sum_{r}\sum_{(j,c')\in\mathcal{P}_r}\mathrm{PIP}_{rjc'}},''',
 'plain': 'pi_c = (sum of PIPs for coding c over datasets and predictors)\n       / (sum of PIPs over all codings, datasets and predictors)'},
paragraph(r'''where $\mathcal{P}_r$ is the set of retained SNP--coding
pairs in dataset $r$. No PIP threshold or credible-set filter was
applied, and we pooled masses before taking the ratio instead of
averaging dataset-specific ratios. The generating proportions were
one for additive coding and zero for the other codings. Percentile
95\% intervals used 500 bootstrap resamples of the 400 replication
seeds, with shared resampling weights across conditions (bootstrap
seed 20260911). The lead-PIP and coding-mass analyses summarize
different quantities: a selected predictor's exact-match frequency
and the distribution of total posterior mass, respectively. Neither
analysis refitted SuSiE or tested the convergence or calibration of
the weighted EM procedure.''')])

new_figures = [
{'stem': 'lead_pip_accuracy_pooled_K', 'label': 'fig:sim-lead-accuracy', 'number': '5', 'supplementary': False,
 'title': 'Exact SNP--coding discovery by global lead PIP.', 'caption': r'''Panels show PVE; colors indicate the called coding.
Each dataset contributes its single highest-PIP predictor when positive
and unambiguous. Points compare mean coding-specific PIP with exact-match
frequency in ten PIP bins, pooling all causal allocations and counts
within PVE. Point size increases with bin count; the dashed line is
equality. Correctness requires both the causal SNP and its generating
coding. Ties and all-zero datasets are excluded from these conditional
rates and counted separately in Methods. Curves are descriptive point
estimates; condition-specific counts and weak-signal failures are reported
in the simulation tables and Results.'''.strip()},
{'stem': 'additive_only_lead_coding_compact', 'label': 'fig:sim-additive-lead', 'number': 'S16', 'supplementary': True,
 'title': 'Lead coding under wholly additive generating effects.', 'caption': r'''Panels show PVE and the horizontal axis gives the number of
causal SNPs. Each colored curve is the number of positive, unambiguous
global leads using that coding divided by all such leads, without a
PIP threshold. All SNPs compete for the lead. Among 7,996 usable
additive-only datasets, 6,561 had unambiguous positive leads; tied
datasets are excluded from these shares and reported in the tables.
These ranking frequencies do not represent coding PIP-mass shares.
No uncertainty intervals are shown.'''.strip()},
{'stem': 'additive_only_pip_allocation_compact', 'label': 'fig:sim-additive-pip-mass', 'number': 'S17', 'supplementary': True,
 'title': 'Coding PIP-mass allocation under wholly additive generating effects.', 'caption': r'''Panels show PVE and the horizontal axis gives the number of
causal SNPs. Each share pools all PIPs of one coding across usable
datasets and divides by the total PIP mass over all three codings.
Every retained predictor contributes, without a credible-set or PIP
filter. The true generating proportions are one for additive coding
and zero for recessive and dominant coding. Bars show 95\% percentile
intervals from 500 bootstrap resamples of shared replication seeds.
Nonzero non-additive mass under additive truth limits interpreting these
shares as proportions of biological effects.'''.strip()}]
data['figures'] = ([f for f in data['figures'] if not f['supplementary']] + [new_figures[0]] +
                   [f for f in data['figures'] if f['supplementary']] + new_figures[1:])
data['refs'].update({f['label']: f['number'] for f in new_figures})
data['refs'][table['label']] = table['number']

def latex_table(b):
    return ('\\begin{table}[t]\n\\centering\n\\small\n'
        + '\\caption{\\textbf{' + b['title'] + '} ' + b['text'] + '}\n'
        + '\\label{' + b['label'] + '}\n\\begin{tabular}{lrrrr}\n\\hline\n'
        + ' & '.join(s.replace('%', r'\%') for s in b['columns']) + r' \\' + '\n\\hline\n'
        + '\n'.join(' & '.join(row) + r' \\' for row in b['rows'])
        + '\n\\hline\n\\end{tabular}\n\\end{table}')
def latex_blocks(d, part):
    output = []
    for b in d['blocks'][part]:
        if b['kind'] == 'heading': output.append('\\' + ('subsection' if part == 'methods' else 'paragraph') + '{' + b['text'] + '}')
        elif b['kind'] == 'equation': output.append('\\begin{equation}\n' + b['text'] + '\n\\end{equation}')
        elif b['kind'] == 'table': output.append(latex_table(b))
        else: output.append(b['text'])
    return '\n\n'.join(output) + '\n'
for part in ('results', 'discussion', 'methods'):
    old = latex_blocks(old_data, part)
    assert tex.count(old) == 1, f'Current {part} differ from recorded section text; inspect before editing.'
    tex = tex.replace(old, latex_blocks(data, part), 1)
def latex_figure(f):
    return ('\\begin{figure}[p]\n\\centering\n'
        + '\\includegraphics[width=\\textwidth]{simulation_figures/' + f['stem'] + '.pdf}\n'
        + '\\caption{\\textbf{' + f['title'] + '}\n' + f['caption'] + '}\n'
        + '\\label{' + f['label'] + '}\n\\end{figure}\n\\clearpage\n')
marker = '%%%%%%%%%%%%%%%% REFERENCES'
assert tex.count(marker) == 1
tex = tex.replace(marker, latex_figure(new_figures[0]) + '\n' + marker, 1)
tex = tex.replace(r'\end{document}', '\n'.join(latex_figure(f) for f in new_figures[1:]) + '\n' + r'\end{document}')
tex = tex.replace('Figures S1--S15: simulation recovery, purity, and ROC curves.',
    'Figures S1--S17: simulation recovery, purity, ROC curves, and additive-only coding controls.')
tex = tex.replace('They report point estimates;\nshared replication seeds should be preserved in uncertainty analyses.',
    'Most report point estimates; the coding PIP-mass figure additionally\nshows bootstrap intervals that preserve shared replication seeds.')
tex = tex.replace('% Simulation methods and results updated from completed checkpoints on 10 September 2026.',
    '% Simulation methods and results updated with exact coding discovery on 11 September 2026.')
(STAGE / 'susie_mix_manuscript.tex').write_text(tex, encoding='utf-8')
(STAGE / 'simulation_sections.json').write_text(json.dumps(data, indent=2), encoding='utf-8')

manifest = json.loads((STAGE / 'simulation_tables/figure_manifest.json').read_text())
manifest.extend({'file': f"simulation_figures/{f['stem']}.pdf", 'figure': f['number'],
    'sha256': hashlib.sha256((STAGE / 'simulation_figures' / (f['stem'] + '.pdf')).read_bytes()).hexdigest()} for f in new_figures)
(STAGE / 'simulation_tables/figure_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
note = '''# Coding-discovery manuscript update - 11 September 2026

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
'''
old_note = (TARGET / 'SIMULATION_UPDATE.md').read_text(encoding='utf-8-sig')
(STAGE / 'SIMULATION_UPDATE.md').write_text(note + '\n\n---\n\nHistorical update record:\n\n' + old_note, encoding='utf-8')

# Keep a checkable record for a guarded copy into the separate manuscript repo.
tracked = ['susie_mix_manuscript.tex', 'simulation_sections.json', 'simulation_sections_review.pdf',
           'SIMULATION_UPDATE.md', 'simulation_tables/figure_manifest.json']
hashes = {p: hashlib.sha256((TARGET / p).read_bytes()).hexdigest() for p in tracked}
(STAGE / 'original_hashes.json').write_text(json.dumps(hashes, indent=2))
(STAGE / 'source_before_sha256.txt').write_text(hashlib.sha256(raw).hexdigest())

# Verify all content outside authorized simulation edits against the original.
reconstructed = tex
for f in new_figures:
    reconstructed = reconstructed.replace(latex_figure(f) + '\n', '', 1)
for part in ('results', 'discussion', 'methods'):
    reconstructed = reconstructed.replace(latex_blocks(data, part), latex_blocks(old_data, part), 1)
# Prefix and GTEx sections are checked directly (the added figure spacing is immaterial).
before = raw.decode('utf-8-sig').replace('\r\n', '\n')
assert tex[:tex.index(r'\section{Results}')] == before[:before.index(r'\section{Results}')]
for start, end in [(r'\paragraph{Descriptive comparison of fine-mapping results.}', r'\section{Discussion}'),
                   (r'\subsection{GTEx data and analysis units}', r'\subsection{Simulation design and genotype sampling}'),
                   ('%%%%%%%%%%%%%%%% REFERENCES', '%%%%%%%%%%%%%%%% SUPPLEMENT LIST')]:
    assert tex[tex.index(start):tex.index(end)] == before[before.index(start):before.index(end)]
print('Staged manuscript sections, numerical table, three figure additions, and provenance records.')
