from pathlib import Path
from collections import Counter
import re
import json
import hashlib

root = Path('C:/Document/Serieux/Travail/Package/git/Science-SuSiE-mix')
stage = Path('C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/tmp/manuscript_em_methods_update')
before = (root / 'susie_mix_manuscript.tex').read_text(encoding='utf-8')
after = (stage / 'susie_mix_manuscript.tex').read_text(encoding='utf-8')
bib_before = (root / 'library.bib').read_text(encoding='utf-8')
bib_after = (stage / 'library.bib').read_text(encoding='utf-8')

def uncomment(s):
    return re.sub(r'(?<!\\)%[^\n]*', '', s)

def braces(s):
    stack = []
    for match in re.finditer(r'(?<!\\)[{}]', uncomment(s)):
        if match.group() == '{':
            stack.append(match.start())
        else:
            assert stack, f'Unmatched closing brace at {match.start()}'
            stack.pop()
    assert not stack, f'Unclosed braces: {stack}'

def environments(s):
    stack = []
    for match in re.finditer(r'\\(begin|end)\{([^}]+)\}', uncomment(s)):
        action, env = match.groups()
        if action == 'begin':
            stack.append(env)
        else:
            assert stack and stack.pop() == env, f'Environment mismatch: {match.group()}'
    assert not stack, stack

def citations(s):
    return {key.strip() for group in re.findall(r'\\cite(?:\[[^]]*\])?\{([^}]+)\}', uncomment(s))
            for key in group.split(',')}

def bibkeys(s):
    return re.findall(r'^@\w+\s*\{\s*([^,\s]+)', s, flags=re.MULTILINE)

def labels(s):
    return re.findall(r'\\label\{([^}]+)\}', uncomment(s))

def references(s):
    return set(re.findall(r'\\(?:ref|pageref|eqref)\{([^}]+)\}', uncomment(s)))

def figures(s):
    return re.findall(r'\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}', uncomment(s))

main_title = r'\subsection{Empirical Bayes estimation of coding priors}'
main_next = r'\subsection{Marginal association, permutation controls, and GTEx summaries}'
supp_title = r'\subsection*{Supplementary Methods: empirical Bayes coding priors}'
supp_next = r'\subsection*{Simulation figures}'
main_new = after[after.index(main_title):after.index(main_next)]
supp_new = after[after.index(supp_title):after.index(supp_next)]
new_list = 'Supplementary Methods: empirical Bayes estimation of coding priors.\\\\\n'
restored = after.replace(main_new, '').replace(supp_new, '').replace(new_list, '')
assert restored == before, 'Existing manuscript text was modified outside the three insertion sites.'
assert bib_after.startswith(bib_before.rstrip('\n') + '\n\n')
new_entry = bib_after[len(bib_before.rstrip('\n')):].strip()
assert len(bibkeys(new_entry)) == 1 and bibkeys(new_entry)[0] == 'zhao_adjusting_2024'
braces(new_entry)
braces(after)
environments(after)
assert not [x for x, n in Counter(labels(after)).items() if n > 1]
assert not (references(after) - set(labels(after))) - (references(before) - set(labels(before)))
assert not (citations(after) - set(bibkeys(bib_after))) - (citations(before) - set(bibkeys(bib_before)))
assert citations(after) - citations(before) == {'zhao_adjusting_2024'}
assert len(bibkeys(bib_after)) == len(bibkeys(bib_before)) + 1
assert figures(before) == figures(after)
assert all((root / f).exists() for f in figures(after))
assert len(re.findall(r'(?<!\\)\$', uncomment(main_new + supp_new))) % 2 == 0
assert (main_new + supp_new).count(r'\begin{equation}') == 7
assert (main_new + supp_new).count('zhao_adjusting_2024') == 2
hashes = json.loads((stage / 'original_hashes.json').read_text())
assert all(hashlib.sha256((root / name).read_bytes()).hexdigest() == expected
           for name, expected in hashes.items())
report = {
    'checks': 'PASS',
    'files_changed': ['susie_mix_manuscript.tex', 'library.bib'],
    'main_methods_line': after[:after.index(main_title)].count('\n') + 1,
    'supplement_methods_line': after[:after.index(supp_title)].count('\n') + 1,
    'bibliography_entry_line': bib_after[:bib_after.index('@article{zhao_adjusting_2024,')].count('\n') + 1,
    'supplement_equations_added': 7,
    'figure_paths_unchanged_and_present': len(figures(after)),
    'new_citation_key': 'zhao_adjusting_2024',
    'preexisting_missing_citation_keys': sorted(citations(before) - set(bibkeys(bib_before))),
    'preexisting_missing_reference_labels': sorted(references(before) - set(labels(before))),
    'latex_compilation': 'Not run: no TeX compiler found locally.',
}
(stage / 'validation.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
print(json.dumps(report, indent=2))
