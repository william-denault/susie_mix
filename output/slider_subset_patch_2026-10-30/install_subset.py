"""Install the user-approved 646-ID selection before the expanded run starts."""
import argparse,csv,hashlib,io,json,os,subprocess,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
META='results_slide_prior_em/inplace_expansion_2026-10-29'
ISSUE='RCC annotation differs from the audited Ensembl 103 locus'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def atomic_bytes(p,data):
    p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_name(p.name+'.subset_tmp')
    with tmp.open('wb') as f:f.write(data);f.flush();os.fsync(f.fileno())
    os.replace(str(tmp),str(p))
def atomic_json(p,value):atomic_bytes(p,(json.dumps(value,indent=2,sort_keys=True)+'\n').encode())
def rows(p):
    with p.open(newline='') as f:return list(csv.DictReader(f))
def validate_report(original_map,report,config):
    reference={r['gene']:r for r in rows(original_map)};observed=rows(report)
    if len(observed)!=680 or len({r['gene'] for r in observed})!=680 or set(reference)!=set(r['gene'] for r in observed):
        raise RuntimeError('RCC report does not cover exactly the original 680 candidates')
    failed=[];passing=[]
    for r in observed:
        expected=reference[r['gene']]
        if any(r.get(k)!=v for k,v in expected.items()):raise RuntimeError('Report reference changed: '+r['gene'])
        if r['expression_matches']!='1' or r['annotation_matches']!='1':raise RuntimeError('Unexpected missing/ambiguous ID: '+r['gene'])
        for a,b in [('expression_id','actual_expression_id'),('expression_symbol','actual_expression_symbol'),('annotation_symbol','actual_annotation_symbol'),('chromosome','actual_chromosome')]:
            if r[a]!=r[b]:raise RuntimeError('Unexpected identity mismatch: '+r['gene'])
        endpoints=(int(r['actual_start']),int(r['actual_end']))
        equal=endpoints==(int(r['start']),int(r['end']))
        if r['status']=='PASS' and not r['issue'] and equal:passing.append(r['stable_id'])
        elif r['status']=='FAIL' and r['issue']==ISSUE and not equal:failed.append(r['gene'])
        else:raise RuntimeError('Unexpected mapping-check outcome: '+r['gene'])
    if set(failed)!=set(config['deferred_genes']) or len(failed)!=34 or set(passing)!=set(config['kept_ids']):
        raise RuntimeError('Selection differs from the 646 accepted and 34 explicitly deferred IDs')
    return dict(accepted_ids=646,deferred_ids=34,report_sha256=sha(report))

def verify_complete(target,hashes):
    for rel,digest in hashes.items():
        p=target/rel
        if not p.is_file() or sha(p)!=digest:raise RuntimeError('Unexpected bundle file: '+rel)
def verify_mixed(target,config):
    for rel,digest in config['new_hashes'].items():
        p=target/rel;old=config['old_hashes'].get(rel)
        if not p.exists() and old is None:continue
        if not p.is_file() or sha(p) not in {digest,old}:raise RuntimeError('Unrecognized/modified bundle file: '+rel)

def install_files(target,state,config,after_write=None):
    verify_mixed(target,config)
    # The bundle checksum is switched last. An interrupted install is resumed
    # only from old/new hashes recorded in the durable patch plan.
    order=[r for r in config['changed'] if r!='SHA256SUMS']+['SHA256SUMS']
    for rel in order:
        if (target/rel).is_file() and sha(target/rel)==config['new_hashes'][rel]:continue
        atomic_bytes(target/rel,(HERE/'replacement'/rel).read_bytes())
        if after_write:after_write(rel)
    verify_complete(target,config['new_hashes'])

def run(a):
    import fcntl
    for line in (HERE/'SHA256SUMS').read_text().splitlines():
        digest,rel=line.split('  ',1)
        if sha(HERE/rel)!=digest:raise RuntimeError('Patch checksum failed: '+rel)
    config=json.loads((HERE/'patch.json').read_text())
    project,target,state=[p.resolve() for p in (a.project,a.target,a.state)]
    if project==target or project in target.parents or project==state or project in state.parents or target==state or target in state.parents:
        raise RuntimeError('Use the original HOME bundle and a separate HOME patch-state directory')
    identity=dict(project=str(project),target=str(target),patch_sha256=sha(HERE/'patch.json'))
    state.mkdir(parents=True,exist_ok=True)
    with (state/'.install.lock').open('a') as lock:
        fcntl.flock(lock.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
        receipt=state/'COMPLETE.json'
        if receipt.exists():
            if json.loads(receipt.read_text())['identity']!=identity:raise RuntimeError('Foreign patch receipt')
            verify_complete(target,config['new_hashes']);print('Selection patch already installed; verified.');return
        # storage.py is unchanged and its hash is checked before importing it.
        if sha(target/'storage.py')!=config['old_hashes']['storage.py']:raise RuntimeError('Unknown storage helper')
        sys.path.insert(0,str(target));import storage
        with storage.maintenance_lock(project):
            live=[r for r in storage.queue() if r[1].startswith(('slide_','em_slide_','em_susie_slide'))]
            if live:raise RuntimeError('Slider jobs are still active; patch not installed: '+repr(live[:5]))
            if (project/META).exists() or (project/'slider_inplace_runtime_2026-10-29').exists() or (project/'results_slide_prior_em/iteration_006').exists():
                raise RuntimeError('Expanded setup has already created metadata/runtime/006; do not change its frozen cohort')
            if not (project/'results_slide_prior_em/iteration_005/recovery/parallel/COMPLETE.rds').is_file():
                raise RuntimeError('Original 005 has no completion certificate')
            storage.check_sources(project)
            plan=state/'plan.json'
            if plan.exists():
                saved=json.loads(plan.read_text())
                if saved['identity']!=identity:raise RuntimeError('Foreign/incompatible patch plan')
                verify_mixed(target,config)
            else:
                verify_complete(target,config['old_hashes'])
                report=target/'checks/mapping_checks.csv'
                summary=validate_report(target/'renamed_gene_map.csv',report,config)
                for rel in config['changed']:
                    if rel in config['old_hashes']:atomic_bytes(state/'original'/rel,(target/rel).read_bytes())
                atomic_bytes(state/'mapping_checks.before_subset.csv',report.read_bytes())
                saved=dict(identity=identity,summary=summary,created=time.time())
                atomic_json(plan,saved)
            validate_report(state/'original/renamed_gene_map.csv',state/'mapping_checks.before_subset.csv',config)
            install_files(target,state,config)
            storage.sync_dir(target);storage.sync_dir(state)
            atomic_json(receipt,dict(identity=identity,accepted=646,deferred=34,original_005_unchanged=True,finished=time.time()))
            print('Installed fixed selection: 646 added IDs; 34 deferred and recorded.')
            print('Original 005 results, priors and counts were not changed.')
            print('Backup and original RCC mapping report: '+str(state))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('project','target','state'):p.add_argument(n,type=Path)
    a=p.parse_args()
    try:run(a)
    except Exception as e:print('STOP: '+str(e)+'. No EM job was submitted by this installer.',file=sys.stderr);sys.exit(1)
