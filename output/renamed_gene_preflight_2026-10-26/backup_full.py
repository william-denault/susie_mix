"""Back up an entire live project; never write to the source or stop jobs."""
import argparse, hashlib, json, os, shutil, stat, subprocess, time
from pathlib import Path

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''): h.update(block)
    return h.hexdigest()

def stamp(path):
    s=path.lstat()
    if stat.S_ISLNK(s.st_mode): return dict(type='symlink',target=os.readlink(str(path)))
    if stat.S_ISDIR(s.st_mode): return dict(type='directory')
    if stat.S_ISREG(s.st_mode): return dict(type='file',size=s.st_size,mtime_ns=s.st_mtime_ns)
    return dict(type='other',mode=s.st_mode)

def inventory(project):
    files={}
    for folder,dirs,names in os.walk(str(project),followlinks=False):
        base=Path(folder)
        if base!=project: files[base.relative_to(project).as_posix()]=dict(type='directory')
        for name in list(dirs):
            p=base/name
            if p.is_symlink(): files[p.relative_to(project).as_posix()]=stamp(p); dirs.remove(name)
        for name in names:
            p=base/name
            try: files[p.relative_to(project).as_posix()]=stamp(p)
            except FileNotFoundError: pass
    return files

def verify(project,tree,metadata):
    before=inventory(project); problems=[]; verified=0
    cache={}; manifest_file=metadata/'file_checksums.jsonl'
    if manifest_file.exists():
        with manifest_file.open() as old:
            for line in old:
                try:
                    record=json.loads(line);cache[record['path']]=record
                except (ValueError,KeyError): pass  # interrupted final record
        with manifest_file.open('ab') as old: old.write(b'\n')
    # Append verified checkpoints so another Slurm job can resume a long hash pass.
    with manifest_file.open('a') as manifest:
        for i,(name,info) in enumerate(sorted(before.items()),1):
            source=project/name; target=tree/name; problem=None; sha=None; target_info=None
            try:
                target_info=stamp(target)
                if info['type']=='symlink':
                    if target_info!=info: problem='symlink differs'
                elif info['type']=='directory':
                    if target_info!=info: problem='directory missing or replaced'
                elif info['type']=='file':
                    if target_info['type']!='file' or target_info['size']!=info['size']: problem='type or size differs'
                    elif os.path.samefile(source,target): problem='not an independent copy'
                    else:
                        cached=cache.get(name,{})
                        if cached.get('source')==info and cached.get('target')==target_info and cached.get('sha256') and cached.get('issue') is None:
                            sha=cached['sha256']
                        else:
                            sha=digest(source)
                            if sha!=digest(target): problem='checksum differs'
                else: problem='special file requires separate verification'
                if stamp(source)!=info: problem='source changed during verification'
            except FileNotFoundError: problem='source or backup file disappeared/is missing'
            if problem: problems.append(dict(path=name,issue=problem))
            else: verified+=1
            manifest.write(json.dumps(dict(path=name,source=info,target=target_info,sha256=sha,issue=problem))+'\n')
            if i%100==0: manifest.flush()
            if i%1000==0: print('Checked {}/{} entries'.format(i,len(before)),flush=True)
    after=inventory(project)
    changed=sorted(k for k in set(before)|set(after) if before.get(k)!=after.get(k))
    extras=sorted(set(inventory(tree))-set(after))
    result=dict(source_entries=len(before),verified_entries=verified,problems=problems,
                source_changes_during_verification=changed,backup_only_entries=extras,
                exact_match=not problems and not changed and not extras)
    (metadata/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    return result

def backup(project,dest,resume=False):
    project=project.resolve(strict=True); dest=dest.resolve()
    if dest==project or project in dest.parents or dest in project.parents:
        raise ValueError('Use a separate sibling backup directory, outside the project.')
    if not shutil.which('rsync'): raise ValueError('rsync is required.')
    tree=dest/project.name; metadata=dest/'backup_metadata'
    if dest.exists():
        if not resume: raise ValueError('Destination exists. To resume an interrupted/live copy, use --resume with its exact path.')
        info=json.loads((metadata/'backup_info.json').read_text())
        if info['project']!=str(project) or info['tree']!=str(tree): raise ValueError('Backup belongs to another source.')
        if info['status']=='VERIFIED': raise ValueError('Preserve this verified backup; use a new destination.')
        if tree.is_symlink() or metadata.is_symlink(): raise ValueError('Backup paths must be real directories.')
    else:
        if resume: raise ValueError('No backup to resume here.')
        dest.mkdir(parents=True); tree.mkdir(); metadata.mkdir()
        info=dict(project=str(project),tree=str(tree),created_at_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
                  scope='entire project including hidden files',symlink_policy='preserve links; do not follow external targets')
    info['status']='COPYING'; (metadata/'backup_info.json').write_text(json.dumps(info,indent=2)+'\n')
    # No deletion, in-place updates, hard links to the source, or following external symlinks.
    command=['rsync','-a','--partial','--info=progress2','--',str(project)+'/',str(tree)+'/']
    print('Full project copy: {} -> {}'.format(project,tree),flush=True)
    result=subprocess.run(command,check=False)
    if result.returncode not in (0,24): raise RuntimeError('rsync failed with status {}. Partial backup retained.'.format(result.returncode))
    if result.returncode==24: print('Some live source files vanished during rsync; continuing with the second pass.',flush=True)
    print('Second pass to catch files written during the first pass.',flush=True)
    result=subprocess.run(command,check=False)
    if result.returncode not in (0,24): raise RuntimeError('rsync failed with status {}. Partial backup retained.'.format(result.returncode))
    info['status']='VERIFYING'; (metadata/'backup_info.json').write_text(json.dumps(info,indent=2)+'\n')
    result=verify(project,tree,metadata)
    info['status']='VERIFIED' if result['exact_match'] else 'LIVE_COPY_NEEDS_REVIEW'
    (metadata/'backup_info.json').write_text(json.dumps(info,indent=2)+'\n')
    marker=metadata/('BACKUP_VERIFIED.txt' if result['exact_match'] else 'LIVE_COPY_REPORT.txt')
    marker.write_text('Status: {}\nCopied project: {}\nVerified entries: {} / {}\nProblems: {}\nSource changes during verification: {}\nBackup-only entries: {}\nThis is a live file copy, not an atomic filesystem snapshot.\n'.format(
        info['status'],tree,result['verified_entries'],result['source_entries'],len(result['problems']),
        len(result['source_changes_during_verification']),len(result['backup_only_entries'])))
    print(marker.read_text(),flush=True)
    if not result['exact_match']: print('Inspect backup_metadata/verification.json. Resume when source activity settles; no jobs were stopped.',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('project',type=Path); p.add_argument('destination',type=Path)
    p.add_argument('--resume',action='store_true'); a=p.parse_args(); backup(a.project,a.destination,a.resume)
