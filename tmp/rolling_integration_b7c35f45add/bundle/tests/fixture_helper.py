"""Local test driver, not used by the RCC launchers."""
import contextlib, json, os, sys
from pathlib import Path
bundle=Path(sys.argv[1]).resolve();sys.path.insert(0,str(bundle));sys.dont_write_bytecode=True
import expansion,retention,storage
mode=sys.argv[2];project=Path(sys.argv[3]).resolve()
if mode=='code':
    expansion.make_code(project)
elif mode=='freeze':
    source=Path(sys.argv[4]).resolve()
    storage.atomic_json(project/'identity.json',dict(source=str(source),source_iteration=5,bundle=str(bundle),retention=2,target=40))
    files=list((project/'script').rglob('*.R'))+list((project/'job').iterdir())
    files += [project/n for n in ('renamed_gene_map.csv','verify_expanded.py','retention.py','retention.R','storage.py','maintenance.R','identity.json','results_slide_prior_em/seed_provenance.rds')]
    storage.atomic_json(project/'frozen_expansion.json',dict(hashes={p.relative_to(project).as_posix():storage.hashes(p)['sha256'] for p in files}))
    storage.atomic_json(project/'INITIALIZED.json',dict(test_fixture=True))
    storage.atomic_json(source/'results_slide_prior_em/storage_dependencies.json',[dict(branch=str(project),source_iteration=5)])
elif mode in ('interrupt','delete','finish','step'):
    current=int(sys.argv[4]);state=project/'results_slide_prior_em/rolling_cleanup'/('before_%03d'%current)
    if mode=='step':
        @contextlib.contextmanager
        def local_lock(*args): yield
        retention.barrier_lock=local_lock;retention.check_queue=lambda:None
        retention.step(project,current,project/'queue.txt')
    else:
        if not (state/'plan.json').exists(): retention.plan_files(project,state,current)
        if mode=='interrupt':
            def interrupt(p): raise RuntimeError('simulated interruption after unlink')
            try: retention.execute_plan(project,state,current,after_unlink=interrupt)
            except RuntimeError as e:
                assert str(e)=='simulated interruption after unlink'
                print('Expected interruption; retry will use the original durable plan.')
            else: raise AssertionError('No interruption')
        else:
            plan=retention.execute_plan(project,state,current)
            if mode=='finish': retention.archive_certificates(plan)
else: raise RuntimeError('Unknown fixture command')
