"""Build a fixed, audited 646-ID selection patch; preserve all delivered archives."""
import csv,hashlib,io,json,zipfile
from pathlib import Path
SRC=Path(__file__).resolve().parent;ROOT=SRC.parents[1]
BASE=ROOT/'output/slider_inplace_rcc_2026-10-29.zip'
assert hashlib.sha256(BASE.read_bytes()).hexdigest()=='030764542e08eceef4068ac0e227eb5f150a69b0519dd413a06a927d7548280c'
with zipfile.ZipFile(BASE) as z:
    old={'/'.join(n.split('/')[1:]):z.read(n) for n in z.namelist()}
new=old.copy();deferred=SRC.joinpath('deferred_genes.txt').read_text().splitlines()
assert len(deferred)==len(set(deferred))==34
rows=list(csv.DictReader(io.StringIO(old['renamed_gene_map.csv'].decode())))
assert len(rows)==680 and set(deferred)<=set(r['gene'] for r in rows)
kept=[r for r in rows if r['gene'] not in deferred]
excluded=[dict(r,reason='Deferred at user request: RCC annotation start/end differs from the audited reference locus.') for r in rows if r['gene'] in deferred]
assert len(kept)==646 and len(set(r['expression_symbol'] for r in kept))==645
def csv_bytes(rows):
    s=io.StringIO(newline='');w=csv.DictWriter(s,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows);return s.getvalue().encode()
new['renamed_gene_map.csv']=csv_bytes(kept);new['excluded_candidates.csv']=csv_bytes(excluded)
def replace(name,a,b):
    text=new[name].decode();assert text.count(a)==1,(name,a);new[name]=text.replace(a,b).encode()
replace('preflight.R','nrow(x) == 680L, length(unique(x$expression_symbol)) == 679L','nrow(x) == 646L, length(unique(x$expression_symbol)) == 645L')
replace('inplace.py',"genotype['candidate_ids']!=680","genotype['candidate_ids']!=646")
replace('inplace.py',"'prior_history.before_expansion.csv')","'prior_history.before_expansion.csv','excluded_candidates.csv')")
replace('initialize.R','nrow(map)!=680L','nrow(map)!=646L')
replace('initialize.R','nrow(parent$manifest)+680L','nrow(parent$manifest)+646L')
replace('initialize.R','added_expression_ids=680L','added_expression_ids=646L,deferred_expression_ids=34L')
replace('initialize.R','  # Archive only small provenance/metadata. No posterior files are copied.',
'''  excluded<-read.csv(file.path(bundle,"excluded_candidates.csv"),stringsAsFactors=FALSE)
  if(nrow(excluded)!=34L || anyDuplicated(excluded$gene) || anyDuplicated(excluded$stable_id) ||
     length(intersect(map$gene,excluded$gene)) || length(intersect(map$stable_id,excluded$stable_id))) stop("Invalid explicit deferred-gene list")
  # Archive only small provenance/metadata. No posterior files are copied.''')
replace('initialize.R','  save_expected(manifest,file.path(meta,"manifest.csv"),TRUE)',
'''  save_expected(excluded,file.path(meta,"excluded_candidates.csv"),TRUE)
  save_expected(manifest,file.path(meta,"manifest.csv"),TRUE)''')
new['README.md']=new['README.md'].decode().replace('680','646').replace('679','645').replace('19,148','19,114').replace('18,297','18,263').encode()
new['README.md']+=b'\n## Selection update, 2026-10-30\n\nThe user chose to defer the 34 candidates whose RCC annotation endpoints differ\nfrom the reference. The 646 passing IDs are retained unchanged. Their fixed\nmanifest applies from 006 through 040. The deferred IDs and reasons are preserved\nin `excluded_candidates.csv` and copied to the expansion metadata directory.\nThe original failed 680-row RCC report is backed up by the subset patch installer.\nNo original 005 gene is removed or refitted by this selection change.\n'
new['tests/test_inplace.R']=new['tests/test_inplace.R'].decode().replace('output/slider_inplace_2026-10-29','output/slider_inplace_646_2026-10-30').replace('683L','649L').replace('680','646').replace('tmp/slider_inplace_fixture_path.txt','tmp/slider_646_fixture_path.txt').encode()
new['tests/test_worker_and_cli.R']=new['tests/test_worker_and_cli.R'].decode().replace('tmp/slider_inplace_fixture_path.txt','tmp/slider_646_fixture_path.txt').encode()
def digest(b):return hashlib.sha256(b).hexdigest()
def write(p,b):p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(b)
new['SHA256SUMS']=''.join(digest(new[n])+'  '+n+'\n' for n in sorted(new) if n!='SHA256SUMS').encode()
candidate=ROOT/'output/slider_inplace_646_2026-10-30'
for n,b in new.items():write(candidate/n,b)
dest=ROOT/'output/slider_subset_patch_2026-10-30'
changed=[n for n in sorted(new) if old.get(n)!=new[n]]
for n in changed:write(dest/'replacement'/n,new[n])
config=dict(old_hashes={n:digest(b) for n,b in old.items()},new_hashes={n:digest(b) for n,b in new.items()},
    changed=changed,deferred_genes=deferred,kept_ids=[r['stable_id'] for r in kept],excluded_ids=[r['stable_id'] for r in excluded])
write(dest/'patch.json',(json.dumps(config,indent=2,sort_keys=True)+'\n').encode())
for n in ('install_subset.py','apply_and_resume.sh','README.txt'):
    write(dest/n,(SRC/n).read_text().replace('\r\n','\n').encode())
files=sorted(p for p in dest.rglob('*') if p.is_file() and p!=dest/'SHA256SUMS')
write(dest/'SHA256SUMS',''.join(digest(p.read_bytes())+'  '+p.relative_to(dest).as_posix()+'\n' for p in files).encode())
archive=ROOT/'output/slider_subset_patch_rcc_2026-10-30.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
    for p in files+[dest/'SHA256SUMS']:z.write(p,dest.name+'/'+p.relative_to(dest).as_posix())
print(archive);print('SHA256:',digest(archive.read_bytes()));print('Changed files:',changed)
