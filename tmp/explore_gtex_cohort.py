import csv, gzip, json, re
from pathlib import Path
from collections import Counter, defaultdict

root = Path(__file__).resolve().parents[1]
ref = json.loads((root/'tmp/gtex_v8_api_gene_records.json').read_text())
genes = []
with gzip.open(root/'tmp/Homo_sapiens.GRCh38.103.chr.gtf.gz','rt') as f:
    for line in f:
        if line.startswith('#'): continue
        fields=line.rstrip('\n').split('\t')
        if fields[2]!='gene': continue
        a=dict(re.findall(r'(\w+) "([^"]*)"',fields[8]))
        genes.append(dict(gene_id=a['gene_id'],gene_name=a.get('gene_name',''),gene_type=a.get('gene_biotype',''),chromosome=fields[0],start=int(fields[3]),end=int(fields[4]),strand=fields[6]))
(root/'tmp/ensembl103_gene_records.json').write_text(json.dumps(genes))
byid={r['gene_id']:r for r in genes}
expr=set((root/'data/genes_all.txt').read_text().splitlines())
manifest={r['gene'] for r in csv.DictReader((root/'results_slide_prior_em/iteration_000/manifest.csv').open())}
excluded={r['gene'] for r in csv.DictReader((root/'results_slide_prior_em/chromosome_exclusions.csv').open())}
scan=manifest-excluded
predicted={r['gene_name'] for r in genes if r['gene_type']=='protein_coding'} & expr
print('Ensembl103 genes',len(genes),'predicted selection',len(predicted),'actual',len(manifest),'predicted-only',len(predicted-manifest),'actual-only',len(manifest-predicted))
print('predicted-only',sorted(predicted-manifest))
print('actual-only',sorted(manifest-predicted))
autosomes={'chr'+str(i) for i in range(1,23)}
missing=[r for r in ref if r['geneType']=='protein coding' and r['chromosome'] in autosomes and r['geneSymbol'] not in scan]
counts=Counter(); examples=defaultdict(list)
for r in missing:
    newer=byid.get(r['gencodeId'].split('.')[0])
    if r['geneSymbol'] not in expr: reason='Absent expression symbol'
    elif not newer: reason='ID absent in Ensembl103'
    elif newer['gene_type']!='protein_coding': reason='No longer protein coding'
    elif newer['gene_name']!=r['geneSymbol']: reason='Renamed protein-coding ID'
    else: reason='Other'
    counts[reason]+=1
    if len(examples[reason])<20: examples[reason].append((r['geneSymbol'],newer and newer['gene_name'],newer and newer['gene_type']))
print('missing IDs',len(missing),'counts',counts)
print(json.dumps(examples,indent=2))
