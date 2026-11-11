"""Replace only unstarted online iterations; preserve every started plan row."""
import csv, io, re
from pathlib import Path
from support import read

BATCHES=12
RHO=0.5
TAG='rho05_12_batches_2026-11-06'

def csv_bytes(rows):
    buf=io.StringIO(newline='')
    writer=csv.DictWriter(buf,fieldnames=list(rows[0]),lineterminator='\n')
    writer.writeheader();writer.writerows(rows)
    return buf.getvalue().encode()

def prepare(branch,mb,state):
    cfg=read(branch/'config.json');old=mb.csvrows(branch/'stages.csv')
    checkpoint=mb.csvrows(state/'checkpoint.csv')[0]
    progress=int(checkpoint['stage']);published=int(mb.csvrows(branch/'progress.csv')[0]['stage'])
    genes=(state/'cohort.txt').read_text().splitlines()
    if int(checkpoint['genes'])!=18263 or len(genes)!=18263 or len(set(genes))!=len(genes):
        raise RuntimeError('Unexpected expanded gene cohort')
    if any(not re.fullmatch(r'[A-Za-z0-9_.-]+',g) or g in ('.','..') for g in genes):
        raise RuntimeError('Unsafe gene identifier')
    if [int(r['stage']) for r in old]!=list(range(1,60)) or not 9<=progress<54 or published not in (progress,progress-1):
        raise RuntimeError('Unexpected old plan/progress, or no future online pass remains')
    cohort=set(genes)
    for iteration in range(15,26):
        rows=[r for r in old if int(r['iteration'])==iteration]
        online=iteration<=20;count=9 if online else 1
        expected_phase='minibatch' if online else 'full'
        if len(rows)!=count or [int(r['batch']) for r in rows]!=list(range(1,count+1)):
            raise RuntimeError('Unexpected baseline batch numbering')
        selected=[]
        for row in rows:
            expected_pass=iteration-14 if online else iteration-20
            group='../../results_slide_prior_em/iteration_%03d'%iteration
            suffix='/batches/batch_%02d'%int(row['batch']) if online else '/batches/full_step'
            if row['phase']!=expected_phase or float(row['rho'])!=(.3 if online else 1.) or \
               int(row['pass'])!=expected_pass or row['group']!=group or row['results_dir']!=group+'/results' or row['directory']!=group+suffix:
                raise RuntimeError('Unexpected baseline stage definition')
            p=Path(row['genes_file'])
            if p.is_absolute() or '..' in p.parts or p.parts[0]!='plan':raise RuntimeError('Unexpected plan path')
            batch=(branch/p).read_text().splitlines()
            if not batch or len(set(batch))!=len(batch):raise RuntimeError('Empty/duplicate baseline batch')
            selected+=batch
        if len(selected)!=len(genes) or set(selected)!=cohort:raise RuntimeError('Baseline cohort differs across passes')

    # A prepared directory or submitted task counts as started even if no gene
    # result exists yet. Never repartition/relabel such work.
    started={int(r['iteration']) for r in old if int(r['stage'])<=progress}
    root=Path(cfg['source'])/'results_slide_prior_em'
    for d in root.glob('iteration_*'):
        if not re.fullmatch(r'iteration_[0-9]+',d.name):continue
        iteration=int(d.name[10:])
        if iteration<15:continue
        if d.is_symlink() or d.resolve().parent!=root.resolve() or \
           not (d/'ONLINE_OWNER.json').exists() or read(d/'ONLINE_OWNER.json').get('workflow_id')!=cfg['id']:
            raise RuntimeError('Unexpected ownership of iteration '+str(iteration))
        started.add(iteration)
    for p in (branch/'tasks').glob('task_*/task.json'):
        stage=int(read(p)['stage'])
        if not 1<=stage<=len(old):raise RuntimeError('Task outside baseline plan')
        started.add(int(old[stage-1]['iteration']))
    first=max(started)+1
    if not 16<=first<=20:raise RuntimeError('No unstarted online iteration remains; current run is unchanged')
    prefix=[r.copy() for r in old if int(r['iteration'])<first]
    if progress>len(prefix):raise RuntimeError('Cannot replace a committed stage')
    new=prefix[:];files={};sid=len(prefix)
    for iteration in range(first,21):
        if (root/('iteration_%03d'%iteration)).exists():raise RuntimeError('Future iteration already exists')
        epoch=iteration-14
        for batch,selected in enumerate(mb.balanced(genes,BATCHES,cfg['seed']+epoch),1):
            sid+=1;row=next(r.copy() for r in old if int(r['iteration'])==iteration)
            gf='plan/rho05_12/pass_%03d_batch_%02d.txt'%(epoch,batch)
            if (branch/gf).exists():raise RuntimeError('New plan file already exists outside this transaction')
            row.update(stage=str(sid),batch=str(batch),rho=str(RHO),genes_file=gf,
                       directory=row['group']+'/batches/batch_%02d'%batch)
            new.append(row);files[gf]=('\n'.join(selected)+'\n').encode()
    for original in old:
        if original['phase']=='full':
            sid+=1;row=original.copy();row['stage']=str(sid);new.append(row)
    if new[:len(prefix)]!=old[:len(prefix)] or [int(r['stage']) for r in new]!=list(range(1,len(new)+1)):
        raise RuntimeError('Started stage definitions changed')
    cfg.update(rho=RHO,batches=BATCHES,
        online_schedule=[dict(iteration=i,batches=9 if i<first else BATCHES,rho=.3 if i<first else RHO) for i in range(15,21)],
        schedule_change=dict(tag=TAG,start_iteration=first,previous_rho=.3,previous_batches=9,
                             policy='preserve_every_started_iteration',rho=RHO,batches=BATCHES))
    summary=dict(start_iteration=first,kept_through_iteration=first-1,checkpoint_stage=progress,
                 preserved_plan_stages=len(prefix),total_stages=len(new),rho=RHO,batches=BATCHES,
                 genes=len(genes),batch_sizes=[len(b) for b in mb.balanced(genes,BATCHES,cfg['seed']+first-14)],
                 final_iteration=25)
    return cfg,csv_bytes(new),files,summary
