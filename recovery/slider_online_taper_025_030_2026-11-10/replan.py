"""Preserve 015..025, retaper unstarted 026..028, finish with full 029..030."""
import csv, io, re
from pathlib import Path
from support import read

TAG='taper_025_030_2026-11-10'


def csv_bytes(rows):
    buf=io.StringIO(newline='');writer=csv.DictWriter(buf,fieldnames=list(rows[0]),lineterminator='\n')
    writer.writeheader();writer.writerows(rows)
    return buf.getvalue().encode()


def baseline_schedule():
    schedule=[]
    for i in range(15,39):
        online=i<=20 or 24<=i<=33
        groups=(9 if i<=16 else 18 if i<=20 else 61) if online else 1
        rho=(.3 if i<=16 else .5 if i<=20 else .8 if i<=25 else .5 if i<=28 else .3) if online else 1.
        pass_index=(i-14 if i<=20 else i-20 if i<=23 else i-17 if online else i-30)
        schedule.append(dict(iteration=i,phase='minibatch' if online else 'full',batches=groups,rho=rho,pass_index=pass_index))
    return schedule


def prepare(branch,mb,state):
    cfg=read(branch/'config.json');old=mb.csvrows(branch/'stages.csv')
    checkpoint=mb.csvrows(state/'checkpoint.csv')[0]
    progress=int(checkpoint['stage']);published=int(mb.csvrows(branch/'progress.csv')[0]['stage'])
    genes=(state/'cohort.txt').read_text().splitlines();cohort=set(genes)
    if int(checkpoint['genes'])!=18263 or len(genes)!=18263 or len(cohort)!=18263:
        raise RuntimeError('Unexpected expanded gene cohort')
    if any(not re.fullmatch(r'[A-Za-z0-9_.-]+',g) or g in ('.','..') for g in genes):
        raise RuntimeError('Unsafe gene identifier')
    if cfg.get('iteration_schedule')!=baseline_schedule() or len(old)!=708 or \
       [int(r['stage']) for r in old]!=list(range(1,709)) or not 93<=progress<=215 or published not in (progress,progress-1):
        raise RuntimeError('Expected the 708-stage baseline with 024/025 current or complete')
    for expected in baseline_schedule():
        i=expected['iteration'];rows=[r for r in old if int(r['iteration'])==i]
        if len(rows)!=expected['batches'] or [int(r['batch']) for r in rows]!=list(range(1,len(rows)+1)):
            raise RuntimeError('Unexpected baseline batch numbering')
        selected=[]
        for row in rows:
            group='../../results_slide_prior_em/iteration_%03d'%i
            suffix='/batches/batch_%02d'%int(row['batch']) if expected['phase']=='minibatch' else '/batches/full_step'
            if row['phase']!=expected['phase'] or int(row['pass'])!=expected['pass_index'] or float(row['rho'])!=expected['rho'] or \
               row['group']!=group or row['directory']!=group+suffix or row['results_dir']!=group+'/results':
                raise RuntimeError('Unexpected baseline stage definition')
            p=Path(row['genes_file'])
            if p.is_absolute() or '..' in p.parts or not p.parts or p.parts[0]!='plan':raise RuntimeError('Unexpected gene-list path')
            members=(branch/p).read_text().splitlines()
            if not members or len(set(members))!=len(members):raise RuntimeError('Empty/duplicate baseline batch')
            selected+=members
        if len(selected)!=18263 or set(selected)!=cohort:raise RuntimeError('Baseline cohort differs across scans')
    root=Path(cfg['source'])/'results_slide_prior_em'
    for d in root.glob('iteration_*'):
        if not re.fullmatch(r'iteration_[0-9]+',d.name):continue
        iteration=int(d.name[10:])
        if iteration>=26:raise RuntimeError('Iteration %03d already exists; no started work will be replaced'%iteration)
        if iteration>=15 and (d.is_symlink() or d.resolve().parent!=root.resolve() or not (d/'ONLINE_OWNER.json').exists() or read(d/'ONLINE_OWNER.json').get('workflow_id')!=cfg['id']):
            raise RuntimeError('Unexpected ownership of iteration '+str(iteration))
    for p in (branch/'tasks').glob('task_*/task.json'):
        if not 1<=int(read(p)['stage'])<=215:
            raise RuntimeError('A task for 026 or later already exists; current work is unchanged')
    prefix=[r.copy() for r in old if int(r['iteration'])<=25]
    if len(prefix)!=215:raise RuntimeError('Unexpected protected stage prefix')
    # Existing 026..028 gene files retain exactly their membership and ordering.
    new=prefix[:];files={};schedule=[r.copy() for r in baseline_schedule() if r['iteration']<=25]
    for i,rho in ((26,.8),(27,.5),(28,.3)):
        for row in old:
            if int(row['iteration'])==i:
                row=row.copy();row['rho']=str(rho);new.append(row)
        record=next(r.copy() for r in baseline_schedule() if r['iteration']==i)
        record['rho']=rho;schedule.append(record)
    # The old full-data list has the same complete cohort and frozen runtime order.
    full=next(r for r in old if int(r['iteration'])==34)
    full_genes=(branch/full['genes_file']).read_bytes()
    for i in (29,30):
        group='../../results_slide_prior_em/iteration_%03d'%i
        gf='plan/'+TAG+'/iteration_%03d_full.txt'%i
        if (branch/gf).exists():raise RuntimeError('New gene plan already exists outside this transaction')
        files[gf]=full_genes
        row=full.copy();row.update(stage=str(len(new)+1),iteration=str(i),group=group,
            directory=group+'/batches/full_step',results_dir=group+'/results',genes_file=gf,
            **{'pass':str(i-25)})
        new.append(row);schedule.append(dict(iteration=i,phase='full',batches=1,rho=1.,pass_index=i-25))
    if new[:215]!=old[:215] or len(new)!=400 or [int(r['stage']) for r in new]!=list(range(1,401)):
        raise RuntimeError('Protected prefix or new schedule changed unexpectedly')
    previous_schedule=cfg['iteration_schedule'];previous_change=cfg['schedule_change']
    cfg.update(passes=11,full_steps=5,iteration_schedule=schedule,
        online_schedule=[dict(iteration=r['iteration'],batches=r['batches'],rho=r['rho']) for r in schedule if r['phase']=='minibatch'],
        schedule_change=dict(tag=TAG,start_iteration=26,previous_change=previous_change,
            previous_schedule=previous_schedule,policy='preserve_every_started_iteration',final_iteration=30))
    # runtime_order_file and its bytes are deliberately unchanged.
    summary=dict(start_iteration=26,kept_through_iteration=25,checkpoint_stage=progress,
        preserved_plan_stages=215,total_stages=400,final_iteration=30,genes=18263,batches=61,
        batch_sizes=[300]*24+[299]*37,runtime_order_file=cfg['runtime_order_file'],
        schedule=[r for r in schedule if r['iteration']>=24])
    return cfg,csv_bytes(new),files,summary
