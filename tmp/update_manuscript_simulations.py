"""Stage a code-verified simulation update; never write to the manuscript checkout.

The final, reviewed files are copied separately after inspecting the stage.
Existing simulation PDFs are copied byte-for-byte and remain vector graphics.
"""
from pathlib import Path
import csv
import hashlib
import html
import json
import re
import shutil

ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = Path('C:/Document/Serieux/Travail/Package/git/Science-SuSiE-mix')
STAGE = ROOT / 'tmp/manuscript_simulation_update'
FIGURES = ROOT / 'simulation results/figures'
STAGE.mkdir(parents=True, exist_ok=True)
(STAGE / 'simulation_figures').mkdir(exist_ok=True)
(STAGE / 'simulation_tables').mkdir(exist_ok=True)

blocks = {'results': [], 'methods': [], 'discussion': []}
def heading(part, text):
    blocks[part].append({'kind': 'heading', 'text': text})
def para(part, text):
    blocks[part].append({'kind': 'paragraph', 'text': text.strip()})
def equation(tex, plain):
    blocks['methods'].append({'kind': 'equation', 'text': tex.strip(), 'plain': plain})

heading('results', 'Simulation benchmark using observed genotype structure')
para('results', r'''
We compared additive SuSiE with unweighted SuSiE-mix using simulated
quantitative phenotypes and observed GTEx genotypes. The benchmark used
500 donors per dataset, one to five distinct causal SNPs, and target
proportions of variance explained (PVE) of 10\%, 20\%, 30\%, and 40\%.
We enumerated all 55 allocations of these causal SNPs to additive,
dominant, and recessive effects, giving 220 allocation--PVE conditions
with 400 scheduled replications each. The completed archive contained
87,988 usable datasets and 175,976 converged fits; 12 scheduled datasets
had zero genetic variance and were excluded without replacement.
Both methods were fitted to the same phenotype with an upper bound of
ten single-effect components. Evaluation used biological SNPs, combining
the retained encodings of each SNP in the mixed-model output (Methods).
''')
para('results', r'''
The seven architecture classes comprise the three single-coding models,
the three possible pairs of codings, and all three codings together.
Within each class, results at a given causal count pool every compatible
allocation. Coverage figures retain the causal count on the horizontal
axis; pooled PIP curves combine all causal counts within each
architecture and PVE. The same 400 replication seeds were reused across
conditions, so the condition-specific results are dependent. The
figures report pooled point estimates across these shared replication
blocks.
''')

heading('results', 'Non-additive effects reduce coverage and causal-SNP recovery')
para('results', r'''
Allowing alternative genotype encodings improved credible-set coverage
under non-additive generating models, with the largest gains under
recessive effects (Figs.~\ref{fig:sim-coverage-pure} and
\ref{fig:sim-coverage-mixed}). For example, with three recessive causal
SNPs and PVE of 20\%, additive SuSiE achieved 83.0\% empirical coverage,
whereas SuSiE-mix achieved 96.5\%. The fraction of causal SNPs recovered
in the union of reported credible sets increased from 21.2\% to 79.2\%
(Fig.~\ref{fig:sim-power-pure}). At the same causal count and PVE under
dominant effects, coverage increased from 92.2\% to 95.4\%, and
credible-set recovery increased from 69.5\% to 78.8\%.
''')
para('results', r'''
Gains also occurred when different causal SNPs followed different
encodings. With five causal SNPs spanning all three encodings and PVE
of 40\%, coverage increased from 87.3\% to 93.7\%, and credible-set
recovery increased from 61.1\% to 85.9\%
(Figs.~\ref{fig:sim-coverage-mixed} and \ref{fig:sim-power-mixed}).
These results show that the effect of additive misspecification extends
beyond loss of sensitivity: reported credible sets can also omit every
generating causal SNP more often than the nominal 5\% rate.
''')

heading('results', 'PIP-based detection at comparable empirical false discovery proportions')
para('results', r'''
PIP-based comparisons supported improved discrimination under
non-additive effects (Figs.~\ref{fig:sim-fdr-pure} and
\ref{fig:sim-fdr-mixed}). Pooling all causal counts at PVE of 40\%,
recessive-model detection power near 10\% empirical FDR increased from
6.1\% with additive SuSiE to 31.2\% with SuSiE-mix. Corresponding power
increased from 22.0\% to 30.5\% under dominant effects and from 11.5\%
to 24.8\% for architectures containing all three effect types.
Under additive-only effects, the corresponding powers were similar,
26.9\% and 26.6\%. These descriptive operating points use the last
sampled PIP threshold before the first crossing of 10\% empirical FDR
as the threshold decreases (Methods).
''')
para('results', r'''
The ROC curves likewise showed the clearest separation under recessive
effects, while the additive-only curves were close
(Figs.~\ref{fig:sim-roc-pure} and \ref{fig:sim-roc-mixed}).
ROC curves stratified by causal count are provided in
Figs.~\ref{fig:sim-roc-pure-L1}--\ref{fig:sim-roc-mixed-L5}.
The absolute power at low empirical FDR remained limited, even where
ROC power was high: the false positive rate divides false selections
by all noncausal SNPs, whereas empirical FDR divides them by selected
SNPs. The power--FDR comparison therefore addresses a different and
more selective operating regime. Here empirical FDR is the pooled
false discovery proportion across simulated datasets; these curves
do not establish formal FDR control or PIP calibration.
''')

heading('results', 'Additive-case costs and remaining undercoverage')
para('results', r'''
The expanded model did not improve every performance measure. Under
additive-only effects with three causal SNPs and PVE of 20\%,
credible-set recovery decreased from 77.6\% to 74.7\%, despite similar
coverage (93.1\% versus 93.5\%). Across all additive-only conditions,
pooled credible-set recovery was 65.7\% for SuSiE and 63.3\% for
SuSiE-mix. Furthermore, SuSiE-mix did not uniformly attain nominal
coverage. With five causal SNPs spanning all three encodings and PVE
of 10\%, coverage improved from 81.9\% to 87.6\% but remained below
95\%, and credible-set recovery was only 7.5\% and 8.7\%, respectively.
''')
para('results', r'''
Purity was also architecture dependent. For the three-recessive-SNP,
20\%-PVE example, mean minimum absolute within-set correlation increased
from 0.820 to 0.911; under the corresponding additive-only scenario it
decreased from 0.893 to 0.861 (Figs.~\ref{fig:sim-purity-pure} and
\ref{fig:sim-purity-mixed}). Mixed-model purity is computed on fitted
coding-specific predictors, so it is not a common additive-LD measure
of resolution. These results support a benefit from modeling
non-additivity while showing that broader predictor models retain
power and calibration trade-offs.
''')

heading('methods', 'Simulation design and genotype sampling')
para('methods', r'''
The completed benchmark was generated by \texttt{sim\_workhorse.R} and
the R jobs written by \texttt{write\_jobs.R}. We fixed the sample size
at $n=500$ and the fitted upper bound at $L_{\mathrm{fit}}=10$.
For every true causal count $K=1,\ldots,5$, we enumerated all
nonnegative integer triples $(K_A,K_R,K_D)$ summing to $K$, where zero
excludes that generating coding. There are 55 such allocations. Each
was crossed with PVE $h\in\{0.10,0.20,0.30,0.40\}$, yielding 220
conditions with 400 scheduled replications per condition. Single-coding
architectures permit $K=1,\ldots,5$, two-coding architectures
$K=2,\ldots,5$, and the three-coding architecture $K=3,\ldots,5$.
The true count $K$ is denoted by ``L'' in the PIP figure titles; it is
distinct from the fixed fitted bound of ten components.
''')
para('methods', r'''
For each replication, one pre-extracted GTEx PLINK genotype file was
sampled uniformly from the available \texttt{.raw} files. Thus, the
sampling frame was the available genotype-file collection rather than
a new draw from a genome-wide gene list. SNPs with any missing
genotype or no variation were removed. Across all donors in that file,
genotypes were oriented to the minor allele, SNPs were required to
have minor-allele frequency strictly greater than 0.05, and SNPs with
a Pearson Hardy--Weinberg chi-square $P$ value below $10^{-8}$ were
excluded. We then sampled 500 distinct donors without replacement,
retaining the full-donor allele orientation and local genotype
correlations. No additional MAF/HWE filtering or LD pruning was applied
to this donor subset.
''')
para('methods', r'''
For minor-allele dosage $G_{ij}\in\{0,1,2\}$, the candidate predictors
were $G_{ij}$ (additive), $\mathbf{1}\{G_{ij}\geq1\}$ (dominant), and
$\mathbf{1}\{G_{ij}=2\}$ (recessive). Each coding was filtered separately
to require positive sample standard deviation and a column sum of at
least five. This requires five allele copies for additive predictors,
five carriers for dominant predictors, or five minor-allele homozygotes
for recessive predictors; it does not require five donors in every
genotype class. Consequently, the retained coding blocks can have
different numbers of columns.
''')
para('methods', r'''
The causal-SNP pool was the intersection of the retained additive SNPs
and the retained SNPs of every non-additive coding requested by the
allocation. We sampled $K$ distinct biological SNPs uniformly without
replacement from this common pool and assigned the first $K_A$ additive,
the next $K_R$ recessive, and the final $K_D$ dominant effects. Each
causal SNP had exactly one generating coding, and there were no
interactions. In conditions containing recessive effects, this common
eligibility rule also restricts the additive and dominant causal SNPs
to those retaining a recessive predictor. An ineligible draw or a
zero-variance genetic signal produced a saved error without drawing a
replacement locus. No tissue-specific expression phenotype, covariate
adjustment, or inverse-normal transformation was used in the simulation.
''')

heading('methods', 'Effect sizes, LD and phenotype generation')
para('methods', r'''
Let $W_{ik}$ be the generating encoding of causal SNP $k$ in donor $i$,
and let $\overline{W}_k$ and $s_k$ be its sample mean and standard
deviation. We standardized each causal predictor and independently
sampled signs $a_k\in\{-1,1\}$ with equal probabilities:
''')
equation(r'''
 z_{ik}=\frac{W_{ik}-\overline{W}_k}{s_k},\qquad
 g_i^{(0)}=\sum_{k=1}^{K}a_k z_{ik}.
''', 'z_ik = (W_ik - mean(W_k)) / s_k;     g_i^(0) = sum_k a_k z_ik.')
para('methods', r'''
For target PVE $h$, we applied one common rescaling to all effects:
''')
equation(r'''
 b=\sqrt{\frac{h}{\widehat{\mathrm{Var}}(g^{(0)})}},\qquad
 g_i=b g_i^{(0)},\qquad y_i=g_i+\epsilon_i,\qquad
 \epsilon_i\stackrel{\mathrm{iid}}{\sim}N(0,1-h).
''', 'b = sqrt[h / sample_var(g^(0))];   g_i = b g_i^(0);   y_i = g_i + error_i,\nwith independent error_i ~ Normal(0, 1 - h).')
para('methods', r'''
Sample variance uses denominator $n-1$. This gives
$\widehat{\mathrm{Var}}(g)=h$ while accounting for covariance among
causal predictors. Each individual contribution $b a_k z_k$ has the
same variance $b^2$; in LD this is not an allocation of $h/K$ to each
SNP, because the variance of the sum also includes covariance terms.
The target PVE is defined as sample genetic variance divided by sample
genetic variance plus the generating noise variance. Realized
$\widehat{\mathrm{Var}}(g)/\widehat{\mathrm{Var}}(y)$ can differ from
$h$ because the noise is an independent Gaussian draw. Errors were
neither orthogonalized against the genotypes nor rescaled after drawing.
''')

heading('methods', 'Model fitting and SNP-level evaluation')
para('methods', r'''
We fitted \texttt{susieR::susie} to the retained additive matrix and to
the concatenated additive, recessive, and dominant matrices. Both fits
used the same outcome, $L_{\mathrm{fit}}=10$,
\texttt{standardize = TRUE}, \texttt{estimate\_prior\_method = "optim"},
\texttt{coverage = 0.95}, \texttt{min\_abs\_corr = 0.5}, and
\texttt{max\_iter = 1000}. The true causal count, generating codings,
and residual variance were not supplied. Other controls used the
installed package defaults. The benchmark compares ordinary additive
and unweighted mixed fits, not the weighted or EM-reweighted coding
models. Uniform prior inclusion weights are over retained predictor
columns; in the mixed fit, total prior weight per biological SNP
therefore depends on the number of encodings surviving filtering.
The standardization, prior-variance optimizer, and purity threshold
differ from the GTEx analysis settings described above.
''')
para('methods', r'''
For SuSiE-mix, coding-specific posterior assignment probabilities were
summed within each single-effect component before calculating a PIP
for biological SNP $j$:
''')
equation(r'''
 \mathrm{PIP}_j=1-\prod_{\ell\in\mathcal{A}}
 \left(1-\sum_{c\in C_j}\alpha_{\ell jc}\right).
''', 'PIP_j = 1 - product_l [1 - sum_(c in C_j) alpha_(l,j,c)], over active components l.')
para('methods', r'''
Here $C_j$ is the set of retained encodings of SNP $j$. The fitted
prior variances were retained when calling \texttt{susie\_get\_pip},
so its usual exclusion of inactive components was preserved. Final
coding-specific PIPs were not summed. Mixed credible sets were mapped
to unique biological SNP indices within each set, without reconstructing
new credible sets from SNP-level PIPs. A reported set was considered
covered if it contained at least one generating causal SNP, irrespective
of the member's encoding. This evaluates SNP localization, not correct
identification of a generating coding or a one-to-one correspondence
between components and causal SNPs.
''')
para('methods', r'''
Empirical coverage was the total number of covered reported sets divided
by the total number of reported sets. Credible-set recovery was the
number of generating causal SNPs appearing in the union of reported
sets divided by the total number of generating causal SNPs; each SNP
was counted once within a dataset. Runs with no reported set contributed
zero recovery and no sets to the coverage or purity denominator. Purity
was the mean of the saved minimum absolute within-set correlations,
using the fitted predictor encodings for each method. Singleton purity
was one. We did not recompute mixed-set purity on additive genotypes.
Coverage and purity therefore weight reported sets, whereas recovery
weights causal SNPs.
''')

heading('methods', 'PIP curves, aggregation and computational record')
para('methods', r'''
At each common PIP threshold $t$, a SNP was selected when its PIP was
at least $t$. We summed true positive (TP) and false positive (FP)
counts across the included datasets and calculated power as TP divided
by all causal SNPs, false positive rate as FP divided by all noncausal
SNPs, and pooled empirical FDR as FP divided by TP plus FP. The last
quantity is a pooled false discovery proportion, not the unweighted
mean of per-dataset false discovery proportions. The no-discovery
endpoint was defined as zero power and zero empirical FDR.
''')
para('methods', r'''
Counts were evaluated on a common dense grid with steps of 0.001
between 0.01 and 0.99, logarithmically spaced tails down to $10^{-12}$
and up to $1-10^{-12}$, and endpoints at zero, one, and infinity.
Tied PIPs entered together. Curves join threshold-specific points in
decreasing threshold order; empirical-FDR curves were not sorted by
FDR or replaced by an optimized envelope. The numerical comparisons
near 10\% empirical FDR use the grid point immediately preceding the
first exceedance of 0.10 in this order. Thresholds and achieved rates
are supplied with the summary tables. These truth-based operating
points are descriptive evaluations, not a proposed threshold-selection
rule for data with unknown causal variants.
''')
para('methods', r'''
For plots at fixed $K$, all compatible allocations within an architecture
were pooled. For ``All L pooled'' plots, TP and FP counts were also
summed across all available $K$ before rates were recalculated. This is
equivalent to stacking the SNP observations; it does not average
per-$K$ rates or give each $K$ equal weight. Larger loci contribute
more noncausal SNPs, and higher causal counts contribute more causal
SNPs. In mixed architectures, causal counts admitting more allocations
also contribute more datasets. Available-case weights reflect the
small number of failed draws.
''')
para('methods', r'''
Replication $r$ used seed $1{,}000{,}000+r$, for $r=1,\ldots,400$,
with the same seeds reused across allocations and PVE values. The
archive comprised 221 checkpoints: 220 production-condition files and
an older overlapping checkpoint. Nine duplicate condition--seed records
were removed. Of 88,000 scheduled condition--replication records,
12 had zero genetic variance and 87,988 were usable; no included fit
was flagged as nonconverged. Errors were omitted from both methods
rather than counted as no-discovery outcomes. Saved objects contain
settings, seeds, the genotype-file path, causal SNPs and codings,
standardized effects, convergence status, credible sets, and PIPs;
they do not retain the simulated outcome, sampled donor indices, or
full fitted models.
''')
para('methods', r'''
The figures included here display pooled point estimates without
Monte Carlo uncertainty intervals. The shared seeds define replication
blocks across conditions and must be retained together when assessing
uncertainty or comparing methods. The current figure set does not
include a credible-set-size comparison or PIP reliability diagram.
The analysis scripts, plotting code, checkpoint audit, configuration
counts, and numerical summaries document the reported results. The
cluster's R and package versions and the immutable genotype-file
inventory were not recorded in these compact objects and remain to be
added to the final reproducibility record.
''')

heading('discussion', 'Implications and limits of the simulation evidence')
para('discussion', r'''
These simulations show that an additive restriction can impair both
causal-SNP recovery and credible-set coverage when causal effects are
non-additive. The improvement is most pronounced for recessive effects,
and remains apparent when comparison is based on SNP selection at
similar empirical false discovery proportions. This supports genotype
encoding as a consequential modeling choice, rather than treating a
change in the number of credible sets as sufficient evidence of improved
inference. At the same time, the expanded model incurs a modest
credible-set recovery cost in additive settings and does not uniformly
restore nominal coverage in weak, multicausal settings.
''')
para('discussion', r'''
The benchmark evaluates a defined family of alternatives: Gaussian
phenotypes, a fixed sample size, at most five causal SNPs, and effects
drawn from the same additive, dominant, and recessive encodings offered
to SuSiE-mix. It therefore demonstrates the consequences of omitting
these encodings and the benefit of including the correct family, rather
than robustness to arbitrary nonlinear effects. Recessive effects were
selected only when at least five minor-allele homozygotes were available;
performance for rarer recessive effects is not established. Effect-size
heterogeneity, interactions, untyped causal variants, and additional
sample sizes also remain outside the present design.
''')
para('discussion', r'''
Several distinctions matter for interpretation. Increased credible-set
recovery can reflect larger sets, so comparative set sizes are needed
before attributing every gain to improved resolution. Purity computed
under different codings does not resolve this issue. The point-estimate
comparisons also require uncertainty assessment that preserves shared
seeds across methods and conditions. ROC and power--FDR curves assess
discrimination; they do not establish that the numerical PIPs equal
empirical causal probabilities. Finally, the simulated outcome was not
subjected to the GTEx expression transformations or covariate adjustment,
and its fitting controls differed from the GTEx scan. The simulations
make the real-data disagreements worth investigating, but do not by
themselves validate an additional eQTL signal, a changed colocalization,
or a particular regulatory mechanism.
''')

# Four main figures, then supplementary recovery/purity/ROC figures.
figures = []
def fig(stem, label, title, caption, supplementary=False):
    number = ('S' + str(sum(f['supplementary'] for f in figures) + 1)) if supplementary else str(sum(not f['supplementary'] for f in figures) + 1)
    figures.append(dict(stem=stem, label=label, title=title, caption=caption.strip(),
                        supplementary=supplementary, number=number))

fig('coverage_pure', 'fig:sim-coverage-pure', 'Credible-set coverage under single-coding architectures.', r'''
Rows show additive-only, dominant-only, and recessive-only generating
effects; columns show PVE of 10\%, 20\%, 30\%, and 40\%. The horizontal
axis gives the true number of distinct causal SNPs, $K=1,\ldots,5$.
Blue denotes additive SuSiE and pink denotes SuSiE-mix. Each point is
the number of reported credible sets containing at least one causal
biological SNP divided by all reported sets, pooled over 399--400
usable replications per condition, each with 500 donors. Mixed sets
are mapped to biological SNPs irrespective of coding. Runs without a
set contribute no sets to this denominator. The dashed line marks
nominal 95\% coverage; the shared coverage axis starts at 0.75 and
ends at 1. Points are horizontally offset for visibility. No
uncertainty intervals are shown; both methods use ten fitted components.
''')
fig('coverage_mixed', 'fig:sim-coverage-mixed', 'Credible-set coverage under mixed-coding architectures.', r'''
Rows show additive plus dominant, additive plus recessive, recessive
plus dominant, and all three generating encodings. Columns show PVE.
Within each panel, the horizontal axis gives $K=2,\ldots,5$ for pairs
of encodings or $K=3,\ldots,5$ for all three. Every named encoding is
represented by at least one causal SNP. Points pool all compatible
allocations at fixed $K$ (399--400 usable replications per allocation,
500 donors per dataset). Blue denotes SuSiE and pink SuSiE-mix.
Coverage is the fraction of reported sets containing any causal
biological SNP, as in Fig.~\ref{fig:sim-coverage-pure}; the dashed
line marks 95\%. The y-axis is restricted to 0.75--1. No uncertainty
intervals are shown, and conditions sharing seeds are dependent.
''')
fig('power_fdr_pure_all_L', 'fig:sim-fdr-pure', 'PIP-based power versus empirical FDR under single-coding architectures.', r'''
Rows show additive-only, dominant-only, and recessive-only effects;
columns show PVE. Blue denotes SuSiE and pink SuSiE-mix. Within each
panel, all usable simulations with $K=1,\ldots,5$ are pooled, giving
one curve per method. At each PIP threshold, TP and FP counts are
summed before calculating power as TP divided by all causal SNPs and
empirical FDR as FP/(TP+FP). Thus empirical FDR is a pooled false
discovery proportion, not the mean of per-dataset FDPs. Curves follow
decreasing PIP thresholds, with the no-discovery endpoint set to (0,0).
The displayed empirical-FDR range is 0--0.25. Codings are combined
within components to obtain one mixed PIP per biological SNP.
Each dataset has 500 donors; curves have no uncertainty bands and
do not establish formal FDR control.
''')
fig('power_fdr_mixed_all_L', 'fig:sim-fdr-mixed', 'PIP-based power versus empirical FDR under mixed-coding architectures.', r'''
Rows show additive plus dominant, additive plus recessive, recessive
plus dominant, and all three encodings; columns show PVE. Blue denotes
SuSiE and pink SuSiE-mix. Counts are pooled across every compatible
allocation and all available causal counts ($K=2,\ldots,5$ for pairs;
$K=3,\ldots,5$ for all three), producing one curve per method in each
panel. Higher causal counts and counts admitting more allocations
contribute more causal SNP observations. Axes, SNP-level PIPs,
threshold ordering and empirical-FDR definition are as in
Fig.~\ref{fig:sim-fdr-pure}. Each dataset has 500 donors, the displayed
empirical-FDR range is 0--0.25, and no uncertainty bands are shown.
''')
for group in ('pure', 'mixed'):
    rows = 'additive-only, dominant-only, and recessive-only effects' if group == 'pure' else 'additive plus dominant, additive plus recessive, recessive plus dominant, and all three encodings'
    counts = r'$K=1,\ldots,5$' if group == 'pure' else r'$K=2,\ldots,5$ for pairs of encodings and $K=3,\ldots,5$ for all three'
    fig(f'power_{group}', f'fig:sim-power-{group}', f'Causal-SNP recovery by credible sets under {"single" if group=="pure" else "mixed"}-coding architectures.', rf'''
Rows show {rows}; columns show PVE. The horizontal axis gives {counts}.
Blue denotes SuSiE and pink SuSiE-mix. Recovery is the number of
generating causal biological SNPs found in the union of reported
credible sets divided by all generating causal SNPs. Each SNP counts
once per dataset, regardless of how many sets or codings represent it.
Runs with no credible set have zero recovery. Compatible allocations
are pooled at fixed $K$, with 399--400 usable replications per allocation
and 500 donors per dataset. Points are pooled estimates without
uncertainty intervals. This set-based measure differs from selection
using a PIP threshold and can increase with set size.
''', True)
for group in ('pure', 'mixed'):
    rows = 'additive-only, dominant-only, and recessive-only effects' if group == 'pure' else 'additive plus dominant, additive plus recessive, recessive plus dominant, and all three encodings'
    fig(f'purity_{group}', f'fig:sim-purity-{group}', f'Credible-set purity under {"single" if group=="pure" else "mixed"}-coding architectures.', rf'''
Rows show {rows}; columns show PVE and the horizontal axis gives the
true causal count. Blue denotes SuSiE and pink SuSiE-mix. Points are
means of the saved minimum absolute within-set correlations across
reported credible sets, pooling compatible allocations at fixed count.
Singleton purity is one; no-set runs contribute no sets. Mixed purity
uses fitted coding-specific predictors, whereas additive SuSiE uses
additive genotypes, so the measures do not give both methods a common
additive-LD scale. Both fits use a reporting threshold of 0.5. Each
dataset has 500 donors and no uncertainty intervals are shown.
''', True)
for group in ('pure', 'mixed'):
    rows = 'additive-only, dominant-only, and recessive-only effects' if group == 'pure' else 'additive plus dominant, additive plus recessive, recessive plus dominant, and all three encodings'
    fig(f'roc_{group}_all_L', f'fig:sim-roc-{group}', f'Pooled ROC curves under {"single" if group=="pure" else "mixed"}-coding architectures.', rf'''
Rows show {rows}; columns show PVE. Blue denotes SuSiE and pink
SuSiE-mix. Each panel pools all available causal counts and compatible
allocations into one curve per method by summing TP and FP at each
common PIP threshold. The vertical axis is TP divided by all causal
SNPs; the horizontal axis is FP divided by all noncausal SNPs. The
dashed diagonal indicates equal true and false positive rates.
Only FPR from 0 to 0.25 is displayed. Each dataset has 500 donors;
PIPs and truth are evaluated at the biological-SNP level. Larger loci
contribute more noncausal SNPs. Curves follow threshold order and
have no uncertainty bands. FPR is distinct from empirical FDR.
''', True)
for group, ks in [('pure', range(1,6)), ('mixed', range(2,6))]:
    rows = 'additive-only, dominant-only, and recessive-only effects' if group == 'pure' else 'additive plus dominant, additive plus recessive, recessive plus dominant, and all three encodings'
    for k in ks:
        note = ' The all-three-encodings row is not applicable for two causal SNPs.' if group == 'mixed' and k == 2 else ''
        fig(f'roc_{group}_L{k}', f'fig:sim-roc-{group}-L{k}', f'ROC curves for {k} causal SNP{"s" if k > 1 else ""} under {"single" if group=="pure" else "mixed"}-coding architectures.', rf'''
The true causal count is fixed at $K={k}$, while both fitted models
allow ten components. Rows show {rows}; columns show PVE. Blue
denotes SuSiE and pink SuSiE-mix. Compatible allocations are pooled
within each panel, using 399--400 usable replications per allocation
and 500 donors per dataset. ROC rates, the diagonal reference, and
the displayed FPR range of 0--0.25 are as in
Fig.~\ref{{fig:sim-roc-{group}}}. No uncertainty bands are shown.{note}
''', True)

refs = {f['label']: f['number'] for f in figures}
def latex_blocks(part):
    output = []
    for b in blocks[part]:
        if b['kind'] == 'heading':
            command = 'subsection' if part == 'methods' else 'paragraph'
            output.append('\\' + command + '{' + b['text'] + '}')
        elif b['kind'] == 'equation':
            output.append('\\begin{equation}\n' + b['text'] + '\n\\end{equation}')
        else:
            output.append(b['text'])
    return '\n\n'.join(output) + '\n'

def latex_figure(f):
    return (f'\\begin{{figure}}[p]\n\\centering\n'
            f'\\includegraphics[width=\\textwidth]{{simulation_figures/{f["stem"]}.pdf}}\n'
            f'\\caption{{\\textbf{{{f["title"]}}}\n{f["caption"]}}}\n'
            f'\\label{{{f["label"]}}}\n\\end{{figure}}\n\\clearpage\n')

main = (MANUSCRIPT / 'susie_mix_manuscript.tex').read_text(encoding='utf-8-sig')
# Repair inherited formatting commands and the old example-figure pointer.
main = main.replace(r'\hetero', 'heterogeneous')
main = re.sub(r'\\textcolor\{palette(?:Navy|Red|Orange)\}\{(\\textbf\{[^}]*\})\}', r'\1', main)
main = main.replace('see Fig 1  for description', 'see Methods for definitions')
assert main.count(r'\section{Results}') == 1
main = main.replace(r'\section{Results}', r'\section{Results}' + '\n\n' + latex_blocks('results'), 1)
start = main.index(r'\section{Conclusion}')
end = main.index(r'\section{Methods}', start)
main = main[:start] + '\\section{Discussion}\n\n' + latex_blocks('discussion') + '\n\\newpage\n\n' + main[end:]
start = main.index(r'\subsection{Simulation design and genotype sampling}')
end = main.index('%%%%%%%%%%%%%%%% REFERENCES', start)
main = main[:start] + latex_blocks('methods') + '\n% Simulation figures use the saved full-data PDF outputs.\n\\clearpage\n' + '\n'.join(latex_figure(f) for f in figures if not f['supplementary']) + '\n' + main[end:]
main = main.replace('% Production simulation completion has not been verified.', '% Simulation methods and results updated from completed checkpoints on 10 September 2026.\n% The GTEx numerical results remain placeholders.')
old = 'The statistical consequences of not modeling this heterogeneity are severe: preliminary simulations'
start = main.index(old)
end = main.index('The \\textbf{second}', start)
main = main[:start] + r'''The simulations reported below show that restricting genotype encoding
to additive effects can reduce credible-set coverage and causal-SNP
recovery under non-additive architectures. Including dominant and
recessive predictors improves these outcomes in the evaluated settings,
although the expanded model does not uniformly achieve nominal
coverage (Figs.~\ref{fig:sim-coverage-pure} and
\ref{fig:sim-coverage-mixed}).

''' + main[end:]
# Correct the broken simulation claim, without inventing the unfinished second cohort.
main = main.replace('and prevent credible  from being concentrated on noncausal variants.',
                    'and improve empirical credible-set coverage under non-additive simulation architectures.')
# Replace the template-only supplementary material with actual simulation figures.
start = main.index('%%%%%%%%%%%%%%%% SUPPLEMENT LIST')
main = main[:start] + r'''%%%%%%%%%%%%%%%% SUPPLEMENT LIST %%%%%%%%%%%%%%%
\subsection*{Supplementary materials}
Figures S1--S15: simulation recovery, purity, and ROC curves.\\
Simulation summary tables and checkpoint audit (\texttt{simulation\_tables/}).

\clearpage
\renewcommand{\thefigure}{S\arabic{figure}}
\renewcommand{\thetable}{S\arabic{table}}
\renewcommand{\theequation}{S\arabic{equation}}
\renewcommand{\thepage}{S\arabic{page}}
\setcounter{figure}{0}
\setcounter{table}{0}
\setcounter{equation}{0}
\setcounter{page}{1}
\begin{center}
\section*{Supplementary Materials for\\ \scititle}
William R.P.~Denault, Peter~Carbonetto, Gao~Wang, Matthew~Stephens\\
Corresponding author: williade@uio.no
\end{center}
\subsection*{Simulation figures}
The following figures accompany the simulation methods and results in
the main text. ``L'' in ROC titles denotes the generating causal count
$K$, not the fitted upper bound of ten. The figures use the completed
checkpoint archive described in Methods. They report point estimates;
shared replication seeds should be preserved in uncertainty analyses.
\clearpage
''' + '\n'.join(latex_figure(f) for f in figures if f['supplementary']) + '\n\\end{document}\n'
(STAGE / 'susie_mix_manuscript.tex').write_text(main, encoding='utf-8')

# Snapshot source text for review and provenance before modifying the destination.
(STAGE / 'simulation_sections.json').write_text(json.dumps({'blocks': blocks, 'figures': figures, 'refs': refs}, indent=2), encoding='utf-8')
for f in figures:
    shutil.copyfile(FIGURES / (f['stem'] + '.pdf'), STAGE / 'simulation_figures' / (f['stem'] + '.pdf'))
for name in ('metric_summary.csv', 'file_audit.csv', 'configuration_counts.csv'):
    shutil.copyfile(FIGURES / name, STAGE / 'simulation_tables' / name)
shutil.copyfile(ROOT / 'tmp/power_fdr_first_crossing.csv', STAGE / 'simulation_tables/power_fdr_first_crossing.csv')
manifest = [{'file': 'simulation_figures/' + f['stem'] + '.pdf',
             'figure': f['number'], 'sha256': hashlib.sha256((FIGURES / (f['stem'] + '.pdf')).read_bytes()).hexdigest()} for f in figures]
(STAGE / 'simulation_tables/figure_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
(STAGE / 'source_before_sha256.txt').write_text(hashlib.sha256((MANUSCRIPT / 'susie_mix_manuscript.tex').read_bytes()).hexdigest(), encoding='ascii')

audit = '''# Simulation manuscript update - 10 September 2026

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
'''
(STAGE / 'SIMULATION_UPDATE.md').write_text(audit, encoding='utf-8')
old_audit = (MANUSCRIPT / 'METHODS_AUDIT.md').read_text(encoding='utf-8-sig')
(STAGE / 'METHODS_AUDIT.md').write_text('> Simulation description superseded on 10 September 2026; see [SIMULATION_UPDATE.md](SIMULATION_UPDATE.md) for the completed benchmark, results, and figure inventory. The historical audit below describes the earlier design.\n\n' + old_audit, encoding='utf-8')
print(f'Staged manuscript, {len(figures)} vector figures and evidence tables in {STAGE}')
