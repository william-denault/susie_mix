from pathlib import Path
import hashlib
import json
import re
import difflib

root = Path('C:/Document/Serieux/Travail/Package/git/Science-SuSiE-mix')
stage = Path('C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/tmp/manuscript_em_methods_update')
stage.mkdir(parents=True, exist_ok=True)
names = ['susie_mix_manuscript.tex', 'library.bib']
original = {name: (root / name).read_bytes() for name in names}
texts = {name: data.decode('utf-8').replace('\r\n', '\n') for name, data in original.items()}

main = r'''\subsection{Empirical Bayes estimation of coding priors}
We implemented a tissue-specific empirical Bayes extension to estimate
the additive, recessive, and dominant prior masses from all available
gene-level mixed-coding fits. Following the principle of learning grouped
priors across loci used in cTWAS \cite{zhao_adjusting_2024}, we adapted the
update to SuSiE's categorical single-effect assignments. The three masses
sum to one and are distributed uniformly among retained predictors within
each coding class. The algorithm alternates SuSiE fitting with an update
based on the summed component assignment probabilities (\texttt{alpha}).
When all three classes are represented in every fit, the update normalizes
these sums by the total number of components; otherwise it accounts for
the classes available in each fit. This variational EM procedure targets
a lower bound on the tissue-level working log marginal likelihood.
Supplementary Methods provide the model, derivation, initialization, and
convergence diagnostics.

'''

supplement = r'''\subsection*{Supplementary Methods: empirical Bayes coding priors}
\label{supp:em-coding-priors}

\subsubsection*{Model and empirical Bayes target}
For gene $g$ in tissue $t$, let $X_{gt}$ contain the retained SNP--coding
predictors and let $y_{gt}$ be the processed expression phenotype. SuSiE
represents the regression coefficient vector as the sum of $L_{gt}$
single-effect vectors. Component $\ell$ selects a predictor
$J_{gt\ell}$ and assigns it a zero-mean Gaussian effect with variance $V_{gt\ell}$;
the remaining entries of that component are zero. The regression includes
an intercept and residual variance $\sigma_{gt}^{2}$. Each component has
the same predictor selection probabilities within a gene--tissue fit.
Different components may select the same SNP or different encodings of
that SNP.

Let $c(j)\in\{A,R,D\}$ denote the coding of predictor $j$,
$m_{gtc}$ the number of retained predictors in class $c$, and
$\mathcal{A}_{gt}=\{c:m_{gtc}>0\}$ the available classes. The tissue-specific
parameters satisfy $\pi_{tc}\geq 0$ and $\sum_{c\in\{A,R,D\}}\pi_{tc}=1$.
The predictor selection probability supplied to SuSiE is
\begin{equation}
 w_{gtj}=\Pr(J_{gt\ell}=j\mid\pi_t)
 =\frac{\pi_{t,c(j)}}
 {m_{gt,c(j)}\sum_{d\in\mathcal{A}_{gt}}\pi_{td}}.
 \label{eq:em-predictor-prior}
\end{equation}
Thus, when all classes are available, $\pi_{tc}$ is the total prior mass
for coding $c$ within a single component. Division by $m_{gtc}$ accounts
for unequal numbers of retained columns. If a class is absent, the
denominator conditions selection on the available classes. These
parameters describe categorical prior masses, rather than per-predictor
Bernoulli inclusion probabilities or fractions of phenotypic variance.

We use a working model that treats gene-level likelihoods within each
tissue as independent. Writing $\theta_{gt}$ for the gene-specific
nuisance parameters, the empirical Bayes target is
\begin{equation}
 \mathcal{L}_t(\pi_t,\theta_t)
 =\sum_g\log p(y_{gt}\mid X_{gt},\pi_t,\theta_{gt}).
 \label{eq:em-marginal-target}
\end{equation}
SuSiE supplies a variational approximation $q_{gt}$ to the posterior of
the component assignments and effects, collectively denoted $z_{gt}$
\cite{wang_simple_2020}. The corresponding evidence lower bound (ELBO) is
\begin{equation}
 \mathcal{F}_t(q,\pi_t,\theta_t)
 =\sum_g E_{q_{gt}}\left[
 \log\frac{p(y_{gt},z_{gt}\mid X_{gt},\pi_t,\theta_{gt})}
 {q_{gt}(z_{gt})}\right]
 \leq\mathcal{L}_t(\pi_t,\theta_t).
 \label{eq:em-elbo}
\end{equation}
Alternating SuSiE's iterative Bayesian stepwise selection with the
coding-prior update is therefore a variational EM procedure. It does not
evaluate the exact multi-effect marginal likelihood or guarantee a
global maximum. Correlated expression and overlapping cis-regions also
mean that Eq.~(\ref{eq:em-marginal-target}) is a working objective, rather
than the likelihood of a fitted joint model of all expression traits.

\subsubsection*{Component responsibilities and the M-step}
At outer iteration $r$, the SuSiE fit provides
$\alpha_{gt\ell j}^{(r)}=q_{gt}^{(r)}(J_{gt\ell}=j)$, with each row summing
to one. The expected number of component assignments to coding $c$ is
\begin{equation}
 C_{gtc}^{(r)}
 =\sum_{\ell=1}^{L_{gt}}\sum_{j:c(j)=c}\alpha_{gt\ell j}^{(r)},
 \qquad \sum_{c\in\mathcal{A}_{gt}}C_{gtc}^{(r)}=L_{gt}.
 \label{eq:em-counts}
\end{equation}
Holding these responsibilities and the nuisance parameters fixed, the
part of the expected complete log likelihood that depends on $\pi_t$ is
\begin{equation}
 \begin{array}{rl}
 Q_t(\pi_t\mid q^{(r)})={}&
 \displaystyle\sum_g\sum_{c\in\mathcal{A}_{gt}}
 C_{gtc}^{(r)}\log\pi_{tc}\\[6pt]
 &\displaystyle{}-\sum_g L_{gt}
 \log\left(\sum_{d\in\mathcal{A}_{gt}}\pi_{td}\right)
 +\mathrm{constant}.
 \end{array}
 \label{eq:em-q}
\end{equation}
The terms involving $\log m_{gtc}$ are constant with respect to $\pi_t$.
If every fit contains all three coding classes, the second term in
Eq.~(\ref{eq:em-q}) vanishes on the simplex. Maximizing the first term
subject to $\sum_c\pi_{tc}=1$ gives
\begin{equation}
 \pi_{tc}^{(r+1)}
 =\frac{\sum_g C_{gtc}^{(r)}}{\sum_g L_{gt}}.
 \label{eq:em-closed-update}
\end{equation}
Thus, each component contributes one unit of posterior assignment mass;
the update pools across genes within a tissue. It does not first
normalize each gene's marginal PIPs or select only genes with a reported
credible set.

If some classes are absent, the normalization term in
Eq.~(\ref{eq:em-q}) must be retained. We aggregate the sufficient statistics
over the seven nonempty subsets of $\{A,R,D\}$ and maximize this objective
in log weights using BFGS with an analytic gradient. The objective is
concave in log weights, up to an arbitrary common additive constant.
The implementation uses relative optimization tolerance $10^{-13}$,
at most 2,000 optimization steps, and a maximum absolute gradient
tolerance of $10^{-6}$ for the normalized objective. If groups of classes
are disconnected in their
co-occurrence across fits, their relative total masses are not
identified; those totals retain their previous values. Each M-step is
checked for a nondecrease in $Q_t$ within numerical tolerance. No
pseudocount is added to the posterior assignment counts.

\subsubsection*{Why the update uses component probabilities}
A marginal predictor PIP is a probability of selection by at least one
component, rather than an expected count of component assignments. For
the components included in the PIP calculation, SuSiE combines
responsibilities as
\begin{equation}
 \mathrm{PIP}_{gtj}
 =1-\prod_{\ell\in\mathcal{I}_{gt}}
       (1-\alpha_{gt\ell j}),
 \label{eq:em-pip}
\end{equation}
where $\mathcal{I}_{gt}$ excludes components with negligible prior
variance under the package's PIP convention. Consequently, summed PIPs
by coding, divided by total PIP, do not generally yield the M-step in
Eq.~(\ref{eq:em-closed-update}). For example, with one predictor per
coding and two component rows each equal to $(0.8,0.1,0.1)$ in $(A,R,D)$
order, the M-step gives $(0.8,0.1,0.1)$. The corresponding PIPs are
$(0.96,0.19,0.19)$, whose normalized shares are approximately
$(0.716,0.142,0.142)$.

The M-step includes every row of \texttt{alpha}, without filtering by
credible-set membership, lead PIP, or component variance. In this latent
representation, even a component with $V_{gt\ell}=0$ retains a categorical
assignment. Its posterior assignment equals its prior when the component
has no effect, providing no information favoring a particular coding.
Such components can slow the outer updates, but discarding them would
change the stated EM derivation.

\subsubsection*{Relation to cTWAS}
cTWAS also pools information across loci to learn grouped priors by
empirical Bayes \cite{zhao_adjusting_2024}. Its independent Bernoulli
inclusion model gives a group-specific update equal to the sum of
predictor PIPs divided by the number of predictors in that group
(Methods, Eq.~5 of Zhao et al.). Those per-predictor inclusion
probabilities need not sum to one across groups. Our three coding masses
instead parameterize the categorical selection within each SuSiE
component, leading to Eqs.~(\ref{eq:em-q}) and
(\ref{eq:em-closed-update}). cTWAS additionally estimates shared
group-specific effect variances; this extension estimates coding
selection masses while retaining gene--tissue-specific component and
residual variances. In particular, SuSiE's
\texttt{estimate\_prior\_method = "EM"} setting estimates component
effect variances internally; the outer update described here is what
estimates the coding masses.

\subsubsection*{Initialization, computation, and interpretation}
The original unweighted genome-wide mixed-coding scan supplies the
initial component posteriors. Its uniform weights over retained columns
need not correspond to a common tissue-level coding-mass vector, so this
scan is treated as an initialization. Before the first weighted scan,
the pooled responsibilities determine the first coding-prior update.
Subsequent iterations alternate updates from the preceding weighted
fits and refitting with the new tissue-specific weights in
Eq.~(\ref{eq:em-predictor-prior}). Each refit starts from the preceding
component posterior moments and nuisance variances, while replacing the
old predictor weights with the newly estimated values.

The EM worker fits only the weighted mixed-coding model, using $L=10$,
\texttt{standardize = FALSE}, \texttt{min\_abs\_corr = 0}, and internal
EM estimation of effect variances. It allows 1,000 inner SuSiE iterations
with ELBO tolerance $10^{-5}$. Predictor identities and component
probabilities are validated before pooling, and unconverged fits block
advancement to the next outer iteration. For every tissue and iteration,
the pipeline saves the three priors, assignment counts, number of
contributing fits and components, M-step objective gain, largest absolute
change in a coding mass, and summed source-fit ELBO when available for
every contributing fit. Gene-level fits are retained separately for
each iteration. Outer iterations are requested in batches; a requested
iteration count alone does not establish convergence. Assessment of
stabilization uses the coding-mass changes and ELBO trajectory, with
ELBO comparisons restricted to the same contributing fits and model
settings.

The estimated masses summarize support for coding classes under the
chosen predictor construction, filtering, and effect-size priors. A SNP
appearing in several gene-level fits contributes separately to each
predictor--outcome occurrence; the weights are not estimates of the
fraction of distinct causal SNPs following each mode of inheritance.
The simulation results reported here evaluate the unweighted model and
do not establish calibration of this empirical Bayes extension.

\clearpage
'''

bib = r'''@article{zhao_adjusting_2024,
  title = {Adjusting for genetic confounders in transcriptome-wide association studies improves discovery of risk genes of complex traits},
  author = {Zhao, Siming and Crouse, Wesley and Qian, Sheng and Luo, Kaixuan and Stephens, Matthew and He, Xin},
  journal = {Nature Genetics},
  year = {2024},
  volume = {56},
  pages = {336--347},
  doi = {10.1038/s41588-023-01648-9},
  url = {https://www.nature.com/articles/s41588-023-01648-9}
}
'''

tex = texts[names[0]]
main_anchor = r'\subsection{Marginal association, permutation controls, and GTEx summaries}'
supp_anchor = r'\subsection*{Simulation figures}'
list_anchor = r'\subsection*{Supplementary materials}'
for anchor in [main_anchor, supp_anchor, list_anchor]:
    assert tex.count(anchor) == 1, anchor
assert 'zhao_adjusting_2024' not in texts[names[1]]
assert '10.1038/s41588-023-01648-9' not in texts[names[1]]
updated = {
    names[0]: tex.replace(main_anchor, main + main_anchor).replace(
        supp_anchor, supplement + supp_anchor).replace(
        list_anchor, list_anchor + '\n' + r'Supplementary Methods: empirical Bayes estimation of coding priors.\\'),
    names[1]: texts[names[1]].rstrip('\n') + '\n\n' + bib,
}

for name in names:
    newline = '\r\n' if b'\r\n' in original[name] else '\n'
    (stage / name).write_bytes(updated[name].replace('\n', newline).encode('utf-8'))
    (stage / (name + '.diff')).write_text(''.join(difflib.unified_diff(
        texts[name].splitlines(keepends=True), updated[name].splitlines(keepends=True),
        fromfile='a/' + name, tofile='b/' + name)), encoding='utf-8')
(stage / 'original_hashes.json').write_text(json.dumps({
    name: hashlib.sha256(data).hexdigest() for name, data in original.items()
}, indent=2), encoding='utf-8')
print(f'Staged {len(names)} manuscript files in {stage}')
print('Main Methods words:', len(main.split()))
print('Supplement source words:', len(supplement.split()))
