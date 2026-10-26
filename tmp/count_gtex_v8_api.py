import collections
import concurrent.futures
import csv
import json
import string
import urllib.request
from pathlib import Path

BASE = 'https://gtexportal.org/api/v2/reference/geneSearch'
def query(prefix):
    url = BASE + '?geneId=' + prefix + '&gencodeVersion=v26&genomeBuild=GRCh38%2Fhg38&itemsPerPage=100000'
    with urllib.request.urlopen(url, timeout=45) as response:
        result = json.load(response)
    assert len(result['data']) == result['paging_info']['totalNumberOfItems'], result['paging_info']
    return prefix, result['data']

genes = {}
sizes = {}
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
    for prefix, records in pool.map(query, string.ascii_uppercase + string.digits):
        sizes[prefix] = len(records)
        for record in records:
            assert record['gencodeVersion'] == 'v26' and record['genomeBuild'] == 'GRCh38/hg38'
            key = (record['gencodeId'], record['chromosome'])
            if key in genes:
                assert genes[key] == record
            genes[key] = record
coding = [r for r in genes.values() if r['geneType'].replace(' ', '_') == 'protein_coding']
auto = [r for r in coding if r['chromosome'] in {'chr' + str(i) for i in range(1, 23)}]
manifest = {r['gene'] for r in csv.DictReader(open('results_slide_prior_em/iteration_000/manifest.csv'))}
excluded = {r['gene'] for r in csv.DictReader(open('results_slide_prior_em/chromosome_exclusions.csv'))}
symbols = {r['geneSymbol'] for r in auto}
current = manifest - excluded
result = {
    'source': BASE,
    'parameters': {'gencodeVersion': 'v26', 'genomeBuild': 'GRCh38/hg38'},
    'method': 'Union geneSearch responses for A-Z and 0-9, deduplicate (gencodeId, chromosome), filter geneType protein coding and chr1-chr22',
    'reference_gene_records': len(genes),
    'protein_coding_records_by_chromosome': dict(collections.Counter(r['chromosome'] for r in coding)),
    'protein_coding_gene_ids_all_chromosomes': len({r['gencodeId'] for r in coding}),
    'autosomal_protein_coding_gene_ids': len({r['gencodeId'] for r in auto}),
    'autosomal_protein_coding_gene_symbols': len(symbols),
    'current_autosomal_scan_symbols': len(current),
    'shared_symbols': len(symbols & current),
    'reference_symbols_not_in_scan': len(symbols - current),
    'scan_symbols_not_in_v26_reference_protein_coding': len(current - symbols),
    'query_sizes': sizes,
}
Path('tmp/gtex_v8_api_gene_records.json').write_text(json.dumps(list(genes.values())))
Path('tmp/gtex_v8_protein_coding_counts.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
