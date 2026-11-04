"""Independent, verified copy of one completed SuSiE-mix iteration and its code.

No hard links, deletion flags, source writes, or Slurm job modifications.
"""
import argparse, csv, hashlib, json, os, shutil, subprocess, time
from pathlib import Path

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''): h.update(block)
    return h.hexdigest()

def inventory(project, iteration):
    relative=['results_em/iteration_{:03d}'.format(iteration),'script/scan_tissue_attempt','data/temp_index']
    relative += ['results_em/prior_history.csv','data/genes_all.txt','data/genes_protein_coding.txt']
    relative += [str(p.relative_to(project)) for p in sorted((project/'job').glob('em_susie*')) if p.is_file()]
    for name in ['job/README_em.md','results_em/last_array_job_id.txt','results_em/last_continuation_job_id.txt']:
        if (project/name).is_file(): relative.append(name)
    files={}
    for name in relative:
        source=project/name
        if not source.exists(): raise ValueError('Required backup source missing: '+str(source))
        entries=sorted(source.rglob('*')) if source.is_dir() else [source]
        for p in entries:
            if p.is_symlink() and p.is_dir(): raise ValueError('Inspect directory symlink before backup: '+str(p))
            if p.is_file():
                s=p.stat();files[p.relative_to(project).as_posix()]=[s.st_size,s.st_mtime_ns]
    return relative,files

def check_completed(project,iteration):
    d=project/'results_em'/'iteration_{:03d}'.format(iteration)
    with (d/'manifest.csv').open() as f: genes=[r['gene'] for r in csv.DictReader(f)]
    rows=[]
    for p in sorted((d/'completed').glob('chunk_*.csv')):
        if not p.with_suffix('.done').exists(): raise ValueError('Chunk completion marker missing: '+str(p))
        with p.open() as f: rows.extend(csv.DictReader(f))
    if len(genes)!=len(set(genes)) or len(rows)!=len(genes) or {r['gene'] for r in rows}!=set(genes):
        raise ValueError('Iteration completion records do not cover its manifest exactly.')
    if any(int(r.get('n_nonconverged',0)) for r in rows): raise ValueError('Unconverged fits recorded in source.')
    if any(not (d/'results'/(gene+'.rds')).is_file() for gene in genes): raise ValueError('Result files missing.')
    return len(genes)

def backup(project,dest,iteration=18):
    project=project.resolve(strict=True); dest=dest.resolve()
    if dest==project or project in dest.parents: raise ValueError('Use a backup directory outside the source project.')
    if not shutil.which('rsync'): raise ValueError('rsync is required.')
    n=check_completed(project,iteration)
    roots,before=inventory(project,iteration)
    total=sum(x[0] for x in before.values())
    dest.parent.mkdir(parents=True,exist_ok=True)
    if shutil.disk_usage(dest.parent).free < total: raise ValueError('Insufficient filesystem free space for an independent copy; also check your RCC quota.')
    # Exclusive creation: an existing backup is never overwritten.
    dest.mkdir()
    info={'project':str(project),'iteration':iteration,'created_at_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
          'files':len(before),'bytes':total,'result_records':n,'status':'COPYING',
          'raw_gtex_inputs_copied':False,'hardlinks_used':False}
    (dest/'backup_info.json').write_text(json.dumps(info,indent=2)+'\n')
    for name in roots:
        source=project/name; target=dest/name
        target.parent.mkdir(parents=True,exist_ok=True)
        if source.is_dir():
            target.mkdir(exist_ok=True)
            args=['rsync','-aL','--info=progress2','--',str(source)+'/',str(target)+'/']
        else: args=['rsync','-aL','--',str(source),str(target)]
        subprocess.run(args,check=True)
    print('Verifying SHA-256 for every copied file.',flush=True)
    with (dest/'backup_manifest.jsonl').open('w') as output:
        for i,name in enumerate(sorted(before),1):
            source=project/name; target=dest/name
            if not target.is_file() or source.stat().st_size!=target.stat().st_size:
                raise ValueError('Copy missing or size mismatch: '+name)
            if os.path.samefile(source,target): raise ValueError('Backup is not independent: '+name)
            source_hash=digest(source)
            if source_hash!=digest(target): raise ValueError('Checksum mismatch: '+name)
            output.write(json.dumps({'path':name,'bytes':target.stat().st_size,'sha256':source_hash})+'\n')
            if i%1000==0: print('Verified {}/{} files'.format(i,len(before)),flush=True)
    if inventory(project,iteration)[1]!=before: raise ValueError('Sources changed during backup. Copy retained for inspection; not marked verified.')
    info['status']='VERIFIED'; (dest/'backup_info.json').write_text(json.dumps(info,indent=2)+'\n')
    (dest/'BACKUP_VERIFIED.txt').write_text('All {} copied files matched source SHA-256 checksums.\nNo source files were changed. Raw GTEx inputs remain at their original paths.\n'.format(len(before)))
    print('VERIFIED BACKUP: '+str(dest),flush=True)

def main():
    p=argparse.ArgumentParser(); p.add_argument('project',type=Path); p.add_argument('destination',type=Path)
    p.add_argument('--iteration',type=int,default=18); args=p.parse_args()
    backup(args.project,args.destination,args.iteration)

if __name__=='__main__': main()
