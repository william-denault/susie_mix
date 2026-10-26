"""Read source cohorts and public annotations; write only separate audit files."""
import csv
import gzip
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'outputs/gtex_gene_cohort_audit_2026-10-26'
OUT.mkdir(parents=True, exist_ok=True)

def read_csv(path):
    with (ROOT/path).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))

def write_tsv(name, records):
    with (OUT/name).open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, list(records[0]), delimiter='\t')
        w.writeheader()
        w.writerows(records)

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''): h.update(block)
    return h.hexdigest()

def join(values):
    return '; '.join(sorted({str(v) for v in values if v not in (None, '')}))

def stable(gene_id):
    # Retain the PAR_Y suffix when removing only the numeric version.
    return re.sub(r'\.\d+', '', gene_id)

ref = json.loads((ROOT/'tmp/gtex_v8_api_gene_records.json').read_text())
new = json.loads((ROOT/'tmp/ensembl103_gene_records.json').read_text())
new_by_id = {r['gene_id']: r for r in new}
new_by_symbol = defaultdict(list)
for r in new: new_by_symbol[r['gene_name']].append(r)
ref_by_symbol = defaultdict(list)
for r in ref: ref_by_symbol[r['geneSymbol']].append(r)
expr_symbols = set((ROOT/'data/genes_all.txt').read_text().splitlines())
selected = set((ROOT/'data/genes_protein_coding.txt').read_text().splitlines())
manifest = {r['gene'] for r in read_csv('results_slide_prior_em/iteration_000/manifest.csv')}
excluded = {r['gene'] for r in read_csv('results_slide_prior_em/chromosome_exclusions.csv')}
scan = manifest - excluded
predicted = {r['gene_name'] for r in new if r['gene_type']=='protein_coding'} & expr_symbols
assert selected == manifest == predicted
assert len(scan) == 17617 and len(excluded) == 851

with gzip.open(ROOT/'tmp/gtex_v8_gene_median_tpm.gct.gz', 'rt') as f:
    next(f); dimensions = next(f).strip()
    expression = [{k:r[k] for k in ('Name','Description')} for r in csv.DictReader(f, delimiter='\t')]
assert len(expression)==56200
assert {r['Description'] for r in expression} == expr_symbols
expr_by_id = {stable(r['Name']):r for r in expression}
assert len(expr_by_id)==len(expression)

autosomes = {'chr'+str(i) for i in range(1,23)}
pc = [r for r in ref if r['chromosome'] in autosomes and r['geneType']=='protein coding']
pc_symbols = {r['geneSymbol'] for r in pc}
assert len(pc)==18901 and len(pc_symbols)==18884

all_rows=[]
for r in sorted(pc, key=lambda r:(r['geneSymbol'],r['gencodeId'])):
    stable_id=stable(r['gencodeId'])
    n=new_by_id.get(stable_id)
    e=expr_by_id.get(stable_id)
    in_scan=r['geneSymbol'] in scan
    if in_scan: reason='Selected symbol'
    elif r['geneSymbol'] not in expr_symbols: reason='Absent from GTEx expression list'
    elif n is None: reason='Stable ID absent from Ensembl 103'
    elif n['gene_type']!='protein_coding': reason='Reclassified as non-protein-coding'
    elif n['gene_name']!=r['geneSymbol']: reason='Renamed protein-coding gene'
    else: reason='Unresolved'
    if reason=='Renamed protein-coding gene':
        assert not e or e['Description']==r['geneSymbol']
    notes=[]
    if not e and r['geneSymbol'] in expr_symbols:
        notes.append('This ID is absent from the public expression rows; another ID has the same symbol.')
    if reason=='Renamed protein-coding gene' and n['gene_name'] in scan:
        assert r['geneSymbol']=='CTAGE5'
        notes.append('CTAGE5 becomes MIA2 in release 103; the selected GTEx MIA2 row is a different v26 ID. Review the locus/ID merge before adding.')
    all_rows.append({
        'gtex_v26_symbol':r['geneSymbol'], 'gencode_v26_id':r['gencodeId'],
        'stable_id':stable_id, 'chromosome':r['chromosome'],
        'in_current_scan_by_symbol':'yes' if in_scan else 'no',
        'in_local_expression_symbols':'yes' if r['geneSymbol'] in expr_symbols else 'no',
        'in_public_expression_by_id':'yes' if e else 'no',
        'public_expression_symbol':e['Description'] if e else '',
        'ensembl103_symbol':n['gene_name'] if n else '',
        'ensembl103_biotype':n['gene_type'] if n else '',
        'ensembl103_chromosome':n['chromosome'] if n else '',
        'reason':reason,
        'new_symbol_in_current_scan':'yes' if n and n['gene_name'] in scan else 'no',
        'note':' '.join(notes),
    })
missing=[r for r in all_rows if r['in_current_scan_by_symbol']=='no']
groups=defaultdict(list)
for r in missing: groups[r['gtex_v26_symbol']].append(r)
symbol_rows=[]
for symbol, rows in sorted(groups.items()):
    assert len({r['reason'] for r in rows})==1
    symbol_rows.append({
        'gtex_v26_symbol':symbol,
        'gencode_v26_ids':join(r['gencode_v26_id'] for r in rows),
        'chromosomes':join(r['chromosome'] for r in rows),
        'in_local_expression_symbols':rows[0]['in_local_expression_symbols'],
        'ensembl103_symbols':join(r['ensembl103_symbol'] for r in rows),
        'ensembl103_biotypes':join(r['ensembl103_biotype'] for r in rows),
        'reference_gene_ids':len(rows),
        'reason':rows[0]['reason'],
        'note':join(r['note'] for r in rows),
    })
counts=Counter(r['reason'] for r in symbol_rows)
assert counts=={'Absent from GTEx expression list':505,'Renamed protein-coding gene':679,
               'Reclassified as non-protein-coding':113,'Stable ID absent from Ensembl 103':37}
assert len(missing)==1338 and len(symbol_rows)==1334

extra_rows=[]
for symbol in sorted(scan-pc_symbols):
    old=ref_by_symbol[symbol]
    n=[r for r in new_by_symbol[symbol] if r['gene_type']=='protein_coding']
    assert old and n
    extra_rows.append({
        'current_symbol':symbol,'gencode_v26_ids':join(r['gencodeId'] for r in old),
        'gencode_v26_biotypes':join(r['geneType'] for r in old),
        'ensembl103_ids':join(r['gene_id'] for r in n),
        'ensembl103_biotypes':join(r['gene_type'] for r in n),
        'reason':'Protein-coding in Ensembl 103; different v26 biotype',
    })

pc_by_symbol=defaultdict(list)
for r in pc: pc_by_symbol[r['geneSymbol']].append(r)
duplicates=[]
for symbol, rows in sorted(pc_by_symbol.items()):
    if len(rows)>1:
        duplicates.append({'symbol':symbol,'reference_gene_ids':len(rows),'ids_minus_symbols':len(rows)-1,
                           'gencode_v26_ids':join(r['gencodeId'] for r in rows),
                           'in_current_scan':'yes' if symbol in scan else 'no'})
assert sum(r['ids_minus_symbols'] for r in duplicates)==17
assert 18901-17-len(symbol_rows)+len(extra_rows)==len(scan)

audit=read_csv('results_slide_prior_em/iteration_000/recovery/audit.csv')
audit_status=Counter(r['status'] for r in audit)
scan_status=Counter(r['status'] for r in audit if r['gene'] in scan)

source_paths=[
    'data/genes_all.txt','data/genes_protein_coding.txt',
    'results_slide_prior_em/iteration_000/manifest.csv',
    'results_slide_prior_em/chromosome_exclusions.csv',
    'results_slide_prior_em/iteration_000/recovery/audit.csv',
    'script/scan_tissue_attempt/get_gene_names.R',
    'script/scan_tissue_attempt/get_protein_coding_gene_names.R',
    'tmp/gtex_v8_api_gene_records.json','tmp/Homo_sapiens.GRCh38.103.chr.gtf.gz',
    'tmp/gtex_v8_gene_median_tpm.gct.gz',
]
sources={p:sha(ROOT/p) for p in source_paths}
urls={
    'GTEx reference gene API':'https://gtexportal.org/api/v2/reference/geneSearch',
    'Ensembl 103 reference':'https://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz',
    'GTEx v8 expression metadata':'https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz',
    'GTEx collapsed gene model rules':'https://github.com/broadinstitute/gtex-pipeline/blob/master/gene_model/README.md',
}
summary={
    'audit_date':'2026-10-26','gencode26_autosomal_protein_coding_ids':len(pc),
    'gencode26_autosomal_protein_coding_symbols':len(pc_symbols),
    'original_manifest_symbols':len(manifest),'excluded_X_Y_MT_symbols':len(excluded),
    'current_autosomal_scan_symbols':len(scan),'shared_reference_and_scan_symbols':len(scan&pc_symbols),
    'reference_symbols_absent_from_scan':len(symbol_rows),'reference_ids_under_absent_symbols':len(missing),
    'scan_symbols_outside_v26_protein_coding':len(extra_rows),
    'missing_symbol_reasons':dict(counts),'public_expression_rows':len(expression),
    'public_expression_symbols':len(expr_symbols),'reference_ids_absent_from_public_expression':sum(r['in_public_expression_by_id']=='no' for r in all_rows),
    'ensembl103_filter_reproduces_original_manifest_exactly':True,
    'public_expression_symbols_match_local_list_exactly':True,
    'synced_iteration000_statuses':dict(audit_status),'synced_iteration000_autosomal_statuses':dict(scan_status),
    'source_sha256':sources,'public_sources':urls,
    'limitations':[
        'No live RCC access; active job state was not queried.',
        'Official uncollapsed Ensembl 103 annotation used. The custom RCC GTF is not local; its symbol selection is exactly reproduced.',
        'GTEx public median TPM row metadata used. Its symbol set exactly matches the local list, but the RCC read-count matrix itself is not local.',
        'Absent expression rows are not evidence of low expression. Exact per-gene removal by model construction was not reconstructed.',
        'Stable IDs absent in release 103 may have been retired, merged, or replaced; their individual histories were not resolved.',
    ],
}
write_tsv('missing_gene_symbols.tsv',symbol_rows)
write_tsv('missing_gene_ids.tsv',missing)
write_tsv('renamed_protein_coding_genes.tsv',[r for r in missing if r['reason']=='Renamed protein-coding gene'])
write_tsv('current_extra_symbols.tsv',extra_rows)
write_tsv('duplicate_reference_symbols.tsv',duplicates)
write_tsv('all_reference_gene_ids.tsv',all_rows)
(OUT/'audit_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
(ROOT/'tmp/gtex_cohort_workbook_data.json').write_text(json.dumps(dict(summary=summary,missing=symbol_rows,extras=extra_rows,duplicates=duplicates)))
print(json.dumps({k:v for k,v in summary.items() if k not in ('source_sha256','public_sources','limitations')},indent=2))
