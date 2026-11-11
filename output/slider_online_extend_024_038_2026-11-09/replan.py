"""Preserve 015..023 exactly; append ten online scans and five full steps."""
import csv, io, re, sys
from pathlib import Path
from support import read
sys.path.insert(0, str(Path(__file__).resolve().parent / 'runtime'))
import runtime_order

TAG = 'extend_024_038_2026-11-09'
BATCHES = 61


def csv_bytes(rows):
    f = io.StringIO(newline='')
    w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator='\n')
    w.writeheader(); w.writerows(rows)
    return f.getvalue().encode()


def prepare(branch, mb, state):
    cfg = read(branch / 'config.json'); old = mb.csvrows(branch / 'stages.csv')
    checkpoint = mb.csvrows(state / 'checkpoint.csv')[0]
    progress = int(checkpoint['stage']); published = int(mb.csvrows(branch / 'progress.csv')[0]['stage'])
    genes = (state / 'cohort.txt').read_text().splitlines()
    if int(checkpoint['genes']) != 18263 or len(genes) != 18263 or len(set(genes)) != len(genes):
        raise RuntimeError('Unexpected expanded gene cohort')
    if any(not re.fullmatch(r'[A-Za-z0-9_.-]+', g) or g in ('.', '..') for g in genes):
        raise RuntimeError('Unsafe gene identifier')
    if len(old) != 95 or [int(r['stage']) for r in old] != list(range(1, 96)) or progress not in (92, 93) or published not in (progress, progress-1):
        raise RuntimeError('Expected completed 022 or 023 with no work after 023')
    cohort = set(genes)
    for iteration in range(15, 26):
        rows = [r for r in old if int(r['iteration']) == iteration]
        online = iteration <= 20
        count, rho = mb.online_parameters(cfg, iteration) if online else (1, 1.)
        if len(rows) != count or [int(r['batch']) for r in rows] != list(range(1, count+1)):
            raise RuntimeError('Unexpected baseline batch numbering')
        selected = []
        for row in rows:
            group = '../../results_slide_prior_em/iteration_%03d' % iteration
            suffix = '/batches/batch_%02d' % int(row['batch']) if online else '/batches/full_step'
            if row['phase'] != ('minibatch' if online else 'full') or float(row['rho']) != rho or int(row['pass']) != (iteration-14 if online else iteration-20) or row['group'] != group or row['results_dir'] != group+'/results' or row['directory'] != group+suffix:
                raise RuntimeError('Unexpected baseline stage definition')
            p = Path(row['genes_file'])
            if p.is_absolute() or '..' in p.parts or not p.parts or p.parts[0] != 'plan':
                raise RuntimeError('Unexpected plan path')
            batch = (branch/p).read_text().splitlines()
            if not batch or len(set(batch)) != len(batch):
                raise RuntimeError('Empty/duplicate baseline batch')
            selected += batch
        if len(selected) != len(genes) or set(selected) != cohort:
            raise RuntimeError('Baseline cohort differs across scans')
    root = Path(cfg['source']) / 'results_slide_prior_em'
    for d in root.glob('iteration_*'):
        if not re.fullmatch(r'iteration_[0-9]+', d.name):
            continue
        i = int(d.name[10:])
        if i >= 24:
            raise RuntimeError('Iteration %03d already exists: no started iteration will be replaced' % i)
        if i >= 15 and (d.is_symlink() or not (d/'ONLINE_OWNER.json').exists() or read(d/'ONLINE_OWNER.json').get('workflow_id') != cfg['id']):
            raise RuntimeError('Unexpected ownership of iteration '+str(i))
    for p in (branch/'tasks').glob('task_*/task.json'):
        stage = int(read(p)['stage'])
        if not 1 <= stage <= 93:
            raise RuntimeError('A task after 023 already exists; current work is unchanged')
    prefix = [r.copy() for r in old if int(r['iteration']) <= 23]
    if len(prefix) != 93:
        raise RuntimeError('Unexpected protected prefix')
    print('Reading saved worker timings to order future submissions; no result RDS files are loaded.',flush=True)
    observations = runtime_order.collect(branch, mb, genes, prefix)
    runtime_file = 'plan/'+TAG+'/runtime_estimates.csv'
    files = {runtime_file: csv_bytes(observations)}
    estimates = {r['gene']: r['seconds'] for r in observations}
    new = prefix[:]; sid = len(prefix)
    schedule = []
    for i in range(15, 24):
        r = [r for r in prefix if int(r['iteration']) == i]
        schedule.append(dict(iteration=i, phase=r[0]['phase'], batches=len(r), rho=float(r[0]['rho']), pass_index=int(r[0]['pass'])))
    for i in range(24, 39):
        online = i <= 33
        rho = .8 if i <= 25 else .5 if i <= 28 else .3 if online else 1.
        count = BATCHES if online else 1
        pass_index = i-17 if online else i-30
        schedule.append(dict(iteration=i, phase='minibatch' if online else 'full', batches=count, rho=rho, pass_index=pass_index))
        # Seed depends on iteration; ordering is applied after balanced assignment.
        batches = mb.balanced(genes, count, cfg['seed']+i-14) if online else [genes]
        for batch, selected in enumerate(batches, 1):
            sid += 1; group = '../../results_slide_prior_em/iteration_%03d' % i
            gf = 'plan/'+TAG+'/iteration_%03d_batch_%02d.txt' % (i, batch)
            if (branch/gf).exists():
                raise RuntimeError('New plan file exists outside this transaction')
            files[gf] = ('\n'.join(runtime_order.order(selected, estimates))+'\n').encode()
            row = dict(prefix[-1]); row.update(stage=str(sid), phase='minibatch' if online else 'full',
                **{'pass':str(pass_index)}, batch=str(batch), rho=str(rho), group=group, iteration=str(i),
                directory=group+('/batches/batch_%02d'%batch if online else '/batches/full_step'),
                genes_file=gf, results_dir=group+'/results')
            new.append(row)
    previous_change = cfg['schedule_change']; previous_schedule = cfg['online_schedule']
    cfg.update(passes=16, full_steps=8, batches=BATCHES, rho=.8, iteration_schedule=schedule,
        online_schedule=[dict(iteration=r['iteration'],batches=r['batches'],rho=r['rho']) for r in schedule if r['phase']=='minibatch'],
        runtime_order_file=runtime_file,
        schedule_change=dict(tag=TAG, start_iteration=24, previous_change=previous_change,
            previous_schedule=previous_schedule, policy='preserve_every_started_iteration', final_iteration=38))
    summary = dict(start_iteration=24, kept_through_iteration=23, checkpoint_stage=progress,
        preserved_plan_stages=93, total_stages=len(new), final_iteration=38, genes=len(genes),
        batches=BATCHES, batch_sizes=[len(b) for b in mb.balanced(genes,BATCHES,cfg['seed']+10)],
        observed_runtimes=sum(r['observed'] for r in observations), runtime_order_file=runtime_file,
        schedule=[r for r in schedule if r['iteration'] >= 24])
    if new[:93] != old[:93] or len(new) != 708:
        raise RuntimeError('Protected prefix or new schedule is inconsistent')
    return cfg, csv_bytes(new), files, summary
