"""Build an isolated RCC recovery distribution; never replace active scripts."""
import csv, hashlib, json, shutil, zipfile
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
NAME = 'renamed_gene_preflight_2026-10-26'
DEST = ROOT/'output'/NAME
DEST.mkdir(parents=True, exist_ok=True)
src = ROOT/'script/scan_tissue_attempt'
scripts = DEST/'script/scan_tissue_attempt'
scripts.mkdir(parents=True, exist_ok=True)
for name in ['em_utils.R','slide_prior_em_utils.R','workhorse_utils.R','get_gene_annotations.R']:
    shutil.copyfile(src/name, scripts/name)

for name in ['workhorse_em.R','workhorse_slide_prior_em.R']:
    if name == 'workhorse_slide_prior_em.R':
        with zipfile.ZipFile(ROOT/'output/slide_prior_shared_rcc_2026-10-23.zip') as z:
            text=z.read('script/scan_tissue_attempt/'+name).decode()
    else: text=(src/name).read_text()
    text=text.replace('\r\n','\n')
    text=text.replace('    target_gene = "GTF2H2",','    target_gene = "GTF2H2",\n    target_gene_id = NULL,',1)
    old='''  j <- which(
    gene_info$Description == target_gene
  )'''
    new='''  if (is.null(target_gene_id)) {
    j <- which(gene_info$Description == target_gene)
  } else {
    ids <- sub("\\\\.[0-9]+$", "", gene_info$Name)
    j <- which(ids == target_gene_id)
    if (length(j) != 1L || gene_info$Description[j] != target_gene)
      stop("Recovery requires exactly one expression row with the frozen stable ID and original name.")
  }'''
    assert text.count(old)==1, name
    text=text.replace(old,new)
    old='''  genes <- subset(
    genes,
    gene_name == target_gene
  )'''
    new='''  if (is.null(target_gene_id)) {
    genes <- subset(genes, gene_name == target_gene)
  } else {
    ids <- sub("\\\\.[0-9]+$", "", genes$ensembl)
    genes <- genes[ids == target_gene_id, , drop = FALSE]
    if (nrow(genes) != 1L) stop("Recovery requires one protein-coding annotation for the stable ID.")
  }'''
    assert text.count(old)==1, name
    text=text.replace(old,new).replace('--threads 2 --memory','--threads 1 --memory')
    if name=='workhorse_em.R': text=text.replace('weighted_fit_mix <- susie(', 'weighted_fit_mix <- susieR::susie(')
    (scripts/name).write_text(text, newline='\n')

data=ROOT/'outputs/gtex_gene_cohort_audit_2026-10-26/renamed_protein_coding_genes.tsv'
with data.open() as f: candidates=list(csv.DictReader(f, delimiter='\t'))
available=[r for r in candidates if r['in_public_expression_by_id']=='yes']
counts=Counter(r['gtex_v26_symbol'] for r in available)
annotation={r['gene_id']:r for r in json.loads((ROOT/'tmp/ensembl103_gene_records.json').read_text())}
assert len(available)==680 and len(counts)==679
rows=[]
for r in available:
    a=annotation[r['stable_id']]
    key=r['gtex_v26_symbol']
    if counts[key]>1: key+='__'+r['stable_id']
    rows.append(dict(gene=key,expression_symbol=r['gtex_v26_symbol'],expression_id=r['gencode_v26_id'],
        stable_id=r['stable_id'],annotation_symbol=a['gene_name'],chromosome=a['chromosome'],
        start=a['start'],end=a['end'],strand=a['strand'],note=r['note']))
with (DEST/'renamed_gene_map.csv').open('w', newline='') as f:
    writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
for f in ['preflight.R','check_bim.py','pilot.R','backup_full.py','preflight_job','pilot_job','backup_job','README.txt']:
    shutil.copyfile(Path(__file__).parent/f, DEST/f)
for d in ['tests']:
    (DEST/d).mkdir(exist_ok=True)
    for f in (Path(__file__).parent/d).glob('*'):
        if f.is_file(): shutil.copyfile(f, DEST/d/f.name)
files=sorted(p for p in DEST.rglob('*') if p.is_file() and p.name!='SHA256SUMS')
assert not any('runs' in p.relative_to(DEST).parts for p in files)
(DEST/'SHA256SUMS').write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.relative_to(DEST).as_posix()+'\n' for p in files))
archive=ROOT/'output'/'renamed_gene_preflight_rcc_2026-10-26.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
    for p in files+[DEST/'SHA256SUMS']: z.write(p, NAME+'/'+p.relative_to(DEST).as_posix())
print(archive)
print(hashlib.sha256(archive.read_bytes()).hexdigest())
