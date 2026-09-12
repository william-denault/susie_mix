from pathlib import Path
import difflib
import hashlib
import json
import re

root = Path('C:/Document/Serieux/Travail/Package/git/Science-SuSiE-mix')
stage = Path('C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/tmp/manuscript_em_active_update')
stage.mkdir(parents=True, exist_ok=True)
name = 'susie_mix_manuscript.tex'
raw = (root / name).read_bytes()
before = raw.decode('utf-8').replace('\r\n', '\n')
after = before

def replace(old, new):
    global after
    assert after.count(old) == 1, old
    after = after.replace(old, new)

replace(r'''based on the summed component assignment probabilities (\texttt{alpha}).
When all three classes are represented in every fit, the update normalizes
these sums by the total number of components; otherwise it accounts for
the classes available in each fit. This variational EM procedure targets''',
r'''based on the summed assignment probabilities (\texttt{alpha}) of components
with positive fitted variance. Assignments of exactly zero-variance
components are integrated out. When all three classes are represented in
every fit, the update normalizes these sums by the number of active
components; otherwise it accounts for the classes available in each fit.
This collapsed variational EM procedure targets''')
replace(r'''to one. The expected number of component assignments to coding $c$ is''',
r'''to one. Holding the fitted variances fixed, define the active components
as $\mathcal{H}_{gt}^{(r)}=\{\ell:V_{gt\ell}^{(r)}>0\}$ and let
$K_{gt}^{(r)}=|\mathcal{H}_{gt}^{(r)}|$. The expected number of active
component assignments to coding $c$ is''')
replace(r''' =\sum_{\ell=1}^{L_{gt}}\sum_{j:c(j)=c}\alpha_{gt\ell j}^{(r)},
 \qquad \sum_{c\in\mathcal{A}_{gt}}C_{gtc}^{(r)}=L_{gt}.''',
r''' =\sum_{\ell\in\mathcal{H}_{gt}^{(r)}}\sum_{j:c(j)=c}\alpha_{gt\ell j}^{(r)},
 \qquad \sum_{c\in\mathcal{A}_{gt}}C_{gtc}^{(r)}=K_{gt}^{(r)}.''')
replace(r'''part of the expected complete log likelihood that depends on $\pi_t$ is''',
r'''part of the collapsed expected complete log likelihood depending on $\pi_t$ is''')
replace(r''' &\displaystyle{}-\sum_g L_{gt}''', r''' &\displaystyle{}-\sum_g K_{gt}^{(r)}''')
replace(r''' =\frac{\sum_g C_{gtc}^{(r)}}{\sum_g L_{gt}}.''',
r''' =\frac{\sum_g C_{gtc}^{(r)}}{\sum_g K_{gt}^{(r)}}.''')
replace(r'''Thus, each component contributes one unit of posterior assignment mass;
the update pools across genes within a tissue. It does not first''',
r'''Thus, each active component contributes one unit of posterior assignment
mass; the update pools across genes within a tissue. If no active component
is available, the preceding prior is retained; a tissue without an existing
prior receives a uniform initialization rather than a data-driven estimate.
The update does not first''')
replace(r'''coding and two component rows each equal to $(0.8,0.1,0.1)$ in $(A,R,D)$''',
r'''coding and two active component rows each equal to $(0.8,0.1,0.1)$ in $(A,R,D)$''')
replace(r'''The M-step includes every row of \texttt{alpha}, without filtering by
credible-set membership, lead PIP, or component variance. In this latent
representation, even a component with $V_{gt\ell}=0$ retains a categorical
assignment. Its posterior assignment equals its prior when the component
has no effect, providing no information favoring a particular coding.
Such components can slow the outer updates, but discarding them would
change the stated EM derivation.''',
r'''The collapsed M-step integrates out assignments of components whose
fitted variance is exactly zero. For such a component, the effect is zero
regardless of the selected predictor; summing over its assignment gives
$\sum_j w_{gtj}=1$. This removes its coding-prior contribution without
changing the observed-data likelihood at the fixed fitted variances.
Equivalently, the full variational bound is recovered by setting inactive
assignment distributions to their priors. The remaining responsibilities
give Eqs.~(\ref{eq:em-counts})--(\ref{eq:em-closed-update}).

Every positive-variance row contributes, including components without a
reported credible set and those with variance below the package's numerical
PIP-reporting tolerance. The exact-zero boundary is not replaced by a
positive variance threshold. There is no additional screen on association
$P$ values, lead PIP, CS membership, or purity. Positive fitted variance
defines an active model component; it does not establish that the component
represents a true biological signal.''')
replace(r'''component posterior moments and nuisance variances, while replacing the
old predictor weights with the newly estimated values.''',
r'''component posterior moments and nuisance variances, while replacing the
old predictor weights with the newly estimated values. All $L$ component
rows remain in saved fits and initializations; only the coding-prior
sufficient statistics exclude zero-variance rows. The active set is
determined anew from the fitted variances at each outer update.''')
replace(r'''the pipeline saves the three priors, assignment counts, number of
contributing fits and components, M-step objective gain, largest absolute
change in a coding mass, and summed source-fit ELBO when available for
every contributing fit.''',
r'''the pipeline saves the three priors, active assignment counts, total and
active fit counts, total and active component counts, zero-variance counts,
collapsed M-step objective gain, and largest absolute change in a coding
mass. The summed source-fit ELBO includes all valid fits, including fits
with no active components, and is recorded when available for every fit.''')
replace('every fit. Gene-level fits are retained separately for\neach iteration.',
        'every fit. Gene-level fits are retained separately for each iteration.')

def block(text, start, end):
    return text[text.index(start):text.index(end)]

main_start = r'\subsection{Empirical Bayes estimation of coding priors}'
main_end = r'\subsection{Marginal association, permutation controls, and GTEx summaries}'
supp_start = r'\subsection*{Supplementary Methods: empirical Bayes coding priors}'
supp_end = r'\subsection*{Simulation figures}'
restored = after
for start, end in [(main_start, main_end), (supp_start, supp_end)]:
    restored = restored.replace(block(after, start, end), block(before, start, end))
assert restored == before, 'Unexpected edit outside the EM methods.'

text = re.sub(r'(?<!\\)%[^\n]*', '', after)
depth = 0
for brace in re.findall(r'(?<!\\)[{}]', text):
    depth += 1 if brace == '{' else -1
    assert depth >= 0
assert depth == 0
stack = []
for command, env in re.findall(r'\\(begin|end)\{([^}]+)\}', text):
    if command == 'begin':
        stack.append(env)
    else:
        assert stack.pop() == env
assert not stack
for command in ['label', 'ref', 'cite', 'includegraphics']:
    pattern = rf'\\{command}(?:\[[^]]*\])?\{{([^}}]+)\}}'
    if command == 'ref':
        assert set(re.findall(pattern, after)) <= set(re.findall(r'\\label\{([^}]+)\}', after))
    else:
        assert re.findall(pattern, before) == re.findall(pattern, after)
assert len(re.findall(r'(?<!\\)\$', text)) % 2 == 0
newline = '\r\n' if b'\r\n' in raw else '\n'
(stage / name).write_bytes(after.replace('\n', newline).encode('utf-8'))
(stage / (name + '.diff')).write_text(''.join(difflib.unified_diff(
    before.splitlines(keepends=True), after.splitlines(keepends=True),
    fromfile='a/' + name, tofile='b/' + name)), encoding='utf-8')
(stage / 'original_hashes.json').write_text(json.dumps({name: hashlib.sha256(raw).hexdigest()}, indent=2), encoding='utf-8')
report = {'checks': 'PASS', 'files': [name], 'scope': 'EM main and supplementary methods only',
          'latex_compilation': 'Not run: no local TeX compiler.',
          'main_line': after[:after.index(main_start)].count('\n') + 1,
          'supplement_line': after[:after.index(supp_start)].count('\n') + 1}
(stage / 'validation.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
print(json.dumps(report, indent=2))

# Reuse the reviewed guarded-copy workflow with a fresh stage and one target.
apply_source = Path('tmp/apply_manuscript_em_methods.ps1').read_text(encoding='utf-8')
apply_source = apply_source.replace('tmp/manuscript_em_methods_update', 'tmp/manuscript_em_active_update')
apply_source = apply_source.replace("@('susie_mix_manuscript.tex', 'library.bib')", "@('susie_mix_manuscript.tex')")
Path('tmp/apply_manuscript_em_active.ps1').write_text(apply_source, encoding='utf-8')
