"""Read the public GTEx v8 gene model and count annotated protein-coding genes."""
import collections
import csv
import json
import re
import time
import urllib.request
from pathlib import Path

URL = 'https://storage.googleapis.com/gtex_analysis_v8/reference/gencode.v26.GRCh38.genes.gtf'
counts = collections.Counter()
ids = collections.defaultdict(set)
symbols = collections.defaultdict(set)
started = time.monotonic()
with urllib.request.urlopen(URL, timeout=60) as stream:
    for raw in stream:
        if raw.startswith(b'#'):
            continue
        parts = raw.decode('utf-8').rstrip('\n').split('\t')
        if len(parts) != 9 or parts[2] != 'gene':
            continue
        attributes = dict(re.findall(r'(\w+) "([^"]*)"', parts[8]))
        if attributes.get('gene_type', attributes.get('gene_biotype')) != 'protein_coding':
            continue
        chrom = parts[0].removeprefix('chr')
        counts[chrom] += 1
        ids[chrom].add(attributes['gene_id'])
        symbols[chrom].add(attributes.get('gene_name', ''))
auto_ids = set().union(*(ids[str(c)] for c in range(1, 23)))
auto_symbols = set().union(*(symbols[str(c)] for c in range(1, 23)))
manifest = {r['gene'] for r in csv.DictReader(open('results_slide_prior_em/iteration_000/manifest.csv'))}
excluded = {r['gene'] for r in csv.DictReader(open('results_slide_prior_em/chromosome_exclusions.csv'))}
expression_symbols = set(Path('data/genes_all.txt').read_text().splitlines())
result = {
    'source': URL,
    'definition': 'unique gene_id on gene-level protein_coding rows; autosomes are chr1 through chr22',
    'records_by_chromosome': dict(counts),
    'unique_ids_by_chromosome': {c: len(v) for c, v in ids.items()},
    'all_chromosome_gene_ids': len(set().union(*ids.values())),
    'autosomal_gene_ids': len(auto_ids),
    'autosomal_gene_symbols': len(auto_symbols),
    'autosomal_symbols_in_local_expression_list': len(auto_symbols & expression_symbols),
    'autosomal_symbols_in_current_autosomal_manifest': len(auto_symbols & (manifest - excluded)),
    'autosomal_symbols_missing_from_current_manifest': len(auto_symbols - (manifest - excluded)),
    'elapsed_seconds': round(time.monotonic() - started, 1),
}
Path('tmp/gtex_v8_protein_coding_counts.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
