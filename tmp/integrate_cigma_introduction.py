"""Prepare the targeted CIGMA introduction/bibliography update for review."""
from pathlib import Path
import hashlib
import json
import re

root = Path(__file__).resolve().parents[1]
target = Path('C:/Document/Serieux/Travail/Package/git/Science-SuSiE-mix')
stage = root / 'tmp/cigma_introduction_update'
stage.mkdir(parents=True, exist_ok=True)
source = (target / 'susie_mix_manuscript.tex').read_bytes()
bib = (target / 'library.bib').read_bytes()
text = source.decode('utf-8-sig').replace('\r\n', '\n')
original_results = text.split(r'\section{Results}', 1)[1]
key = 'chen_cell_type_specific_2026'
assert key not in text and key.encode() not in bib

anchor = 'molecular mechanisms underlying complex traits.\n'
assert text.count(anchor) == 1
context = r'''

Recent work emphasizes the importance of cellular context in this
translation. Chen et al. developed cell-type-informed genetic
mixed-model analysis (CIGMA) to partition gene-expression variation
into contributions from cell-type-shared and cell-type-specific
genetic effects while accounting for cell-to-cell variability
\cite{chen_cell_type_specific_2026}. Using population-scale single-cell
RNA sequencing, they found that genes with cell-type-specific genetic
regulation were enriched for autoimmune-trait heritability and showed
relationships with evolutionary constraint and enhancer complexity,
with replication of broad patterns in a second cohort. Their findings
support the view that bulk-tissue eQTL discovery can underrepresent
regulation relevant to complex traits.
'''
text = text.replace(anchor, anchor + context, 1)

start = text.index('One of the overlooked  assumption')
end = text.index('This misspecification need not result only in reduced power.', start)
bridge = r'''Genotype encoding is a complementary modeling choice. Cell-type
specificity concerns how an effect changes across cellular contexts;
genotype encoding concerns the relationship between allele dosage
and phenotype within a context. A variant can therefore have different
effect sizes across cell types while remaining additive within each.
CIGMA models this context dependence through shared and cell-type-specific
coefficients multiplying genotype dosage and estimates aggregate genetic
variance, rather than assigning posterior probabilities to individual
causal SNPs \cite{chen_cell_type_specific_2026}. It does not address
whether restricting genotype encoding to additive effects distorts
variant-level fine-mapping inference.

Under additive encoding, the modeled effect scales linearly with
allele count. Dominant encoding allows one copy to confer the same
effect as two, whereas recessive encoding requires two copies.
Dominant and recessive effects have been investigated in molecular
and complex traits \cite{cui_dominance_2023,palmer_analysis_2023}.
When such effects generate a phenotype, a fine-mapping model using
only additive predictors may be misspecified, with consequences for
both causal-SNP detection and credible-set calibration.

'''
text = text[:start] + bridge + text[end:]

# Replace unfinished and contradictory closing draft material with a focused
# positioning statement. Preserve the existing prior-work citations above it.
start = text.index('\n However, it remains unclear')
end = text.index(r'\section{Results}', start)
closing = r'''
The practical consequences for credible-set calibration when multiple
causal variants act through different genotype encodings remain an
important question, as does the extent to which these modeling choices
change fine-mapping conclusions in human eQTL data. Here, we compare
additive SuSiE with a mixed-coding analysis that includes additive,
dominant, and recessive predictors. We use simulations based on observed
GTEx genotypes to evaluate causal-SNP recovery, credible-set coverage,
and PIP-based discrimination, and examine how genotype encoding changes
fine-mapping results in GTEx gene-expression data.

This study complements the evidence from CIGMA that cellular context
matters for interpreting genetic regulation. We focus on the form of
the genotype--phenotype relationship and its consequences for localizing
causal variants. The simulations show that allowing alternative
encodings can improve recovery and coverage under non-additive
architectures, while retaining power trade-offs and undercoverage in
some conditions (Figs.~\ref{fig:sim-coverage-pure} and
\ref{fig:sim-coverage-mixed}). Differences between the GTEx fits motivate
further biological investigation; they do not alone establish which
variant or regulatory mechanism is causal.

'''
text = text[:start] + closing + text[end:]

entry = r'''
@article{chen_cell_type_specific_2026,
  title = {Cell-type-specific {eQTLs} underlie the genetic architecture of complex traits},
  author = {Chen, Minhui and Wang, Xinpei and Krockenberger, Lena and Tyebally, Rika and Berg, Jeremy J. and Pott, Sebastian and Flint, Jonathan and Powell, Joseph E. and Balliu, Brunilda and Liu, Xuanyao and Dahl, Andy},
  journal = {Nature},
  year = {2026},
  month = aug,
  doi = {10.1038/s41586-026-10577-6},
  url = {https://www.nature.com/articles/s41586-026-10577-6},
  note = {Published online 26 August 2026}
}
'''
new_bib = bib + (b'' if bib.endswith(b'\n') else b'\n') + entry.encode('utf-8')
assert text.split(r'\section{Results}', 1)[1] == original_results
assert text.count(r'\cite{' + key + '}') == 2
assert new_bib.count(('@article{' + key + ',').encode()) == 1
assert 'the statistical framework to address them does\nnot exist' not in text
assert 'GPT Double CHECK' not in text

newline = '\r\n' if b'\r\n' in source else '\n'
new_source = text.replace('\n', newline).encode('utf-8')
if source.startswith(b'\xef\xbb\xbf'):
    new_source = b'\xef\xbb\xbf' + new_source
assert new_source.split(b'\\section{Results}', 1)[1] == source.split(b'\\section{Results}', 1)[1]
for name, data in [('susie_mix_manuscript.tex', new_source), ('library.bib', new_bib)]:
    (stage / name).write_bytes(data)
(stage / 'original_hashes.json').write_text(json.dumps({
    'susie_mix_manuscript.tex': hashlib.sha256(source).hexdigest(),
    'library.bib': hashlib.sha256(bib).hexdigest()
}, indent=2), encoding='utf-8')

clean = re.sub(r'(?<!\\)%[^\n]*', '', text)
depth = 0
for match in re.finditer(r'(?<!\\)[{}]', clean):
    depth += 1 if match[0] == '{' else -1
    assert depth >= 0
assert depth == 0
labels = set(re.findall(r'\\label\{([^}]+)\}', clean))
assert set(re.findall(r'\\ref\{([^}]+)\}', clean)) <= labels
print('Prepared introduction and reference; results onward unchanged; citation key and LaTeX structure checked.')
