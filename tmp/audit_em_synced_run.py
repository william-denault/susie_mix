from pathlib import Path
import csv
import json
import re
from collections import Counter
from datetime import datetime, timezone

root = Path('results_em')
rows = list(csv.DictReader((root/'prior_history.csv').open(newline='', encoding='utf-8')))
priors = list(csv.DictReader((root/'priors.csv').open(newline='', encoding='utf-8')))
audit = list(csv.DictReader((root/'source_audit.csv').open(newline='', encoding='utf-8')))
manifest = list(csv.DictReader((root/'manifest.csv').open(newline='', encoding='utf-8')))
logs = list((root/'logs').glob('*.err'))
report = {
    'array_id': (root/'last_array_job_id.txt').read_text().strip(),
    'continuation_id': (root/'last_continuation_job_id.txt').read_text().strip(),
    'iterations': sorted({r['iteration'] for r in rows}),
    'tissues': len(rows),
    'methods': sorted({r['update_method'] for r in rows}),
    'history_matches_priors': rows == priors,
    'manifest_genes': len(manifest),
    'manifest_chunks': len({r['chunk'] for r in manifest}),
    'synced_error_logs': len(logs),
    'synced_completion_markers': len(list((root/'completed').glob('*.done'))),
    'local_layout': 'iteration files at results_em root; no iteration_001 subdirectory',
    'node_failures': [],
}
other_errors = []
gene_errors = []
for file in sorted(logs, key=lambda f: int(f.stem.rsplit('_', 1)[1])):
    chunk = int(file.stem.rsplit('_', 1)[1])
    text = file.read_text(errors='replace')
    for line in text.splitlines():
        if 'DUE TO NODE FAILURE' in line:
            detail = re.search(r'JOB (\d+) ON (\S+) CANCELLED AT (\S+) DUE TO NODE FAILURE', line)
            report['node_failures'].append({'chunk': chunk, 'job': detail[1], 'node': detail[2], 'time': detail[3]})
        if re.search(r'ERROR \[', line):
            gene_errors.append({'chunk': chunk, 'message': line})
        elif re.search(r'Error|ERROR|Execution halted|oom.kill|Out Of Memory|unconverged|CANCELLED|TIME LIMIT', line) and 'DUE TO NODE FAILURE' not in line and 'Error: All 67234533 variants in' not in line:
            other_errors.append({'chunk': chunk, 'message': line})
report['other_error_lines'] = other_errors
report['worker_gene_error_count'] = len(gene_errors)
report['worker_gene_error_classes'] = dict(Counter(re.sub(r'ERROR \[[^]]+\]: PLINK failed for \S+ on ', '', r['message']) for r in gene_errors))
report['source_audit_issues'] = dict(Counter(r['issue'] for r in audit))
report['source_gene_error_classes'] = dict(Counter(re.sub(r'PLINK failed for \S+ on ', '', r['message']) for r in audit))
source_error_genes = {r['gene'] for r in audit if r['issue']=='gene_error'}
report['worker_errors_absent_from_source_audit'] = sorted({re.search(r'ERROR \[([^]]+)\]', r['message'])[1] for r in gene_errors} - source_error_genes)
report['prior_checks'] = {
    'all_finite_positive_normalized': all(all(0 < float(r[c]) < 1 for c in ('pi_add','pi_rec','pi_dom')) and abs(sum(float(r[c]) for c in ('pi_add','pi_rec','pi_dom'))-1)<1e-8 for r in rows),
    'all_alpha_totals_match_active_components': all(abs(float(r['alpha_total'])-float(r['n_active_components']))<1e-6 for r in rows),
    'all_component_counts_reconcile': all(int(r['n_components'])==int(r['n_active_components'])+int(r['n_zero_variance_components']) for r in rows),
    'all_mstep_gains_nonnegative': all(float(r['mstep_q_gain'])>=0 for r in rows),
    'all_source_fits_converged': all(int(r['n_nonconverged'])==0 for r in rows),
    'all_source_elbos_complete': all(int(r['n_source_elbo'])==int(r['n_fits']) and r['source_elbo_sum'] != '' for r in rows),
}
report['total_fits'] = sum(int(r['n_fits']) for r in rows)
report['total_active_fits'] = sum(int(r['n_active_fits']) for r in rows)
report['total_components'] = sum(int(r['n_components']) for r in rows)
report['total_active_components'] = sum(int(r['n_active_components']) for r in rows)
report['pooled_active_coding_shares'] = {c: sum(float(r['alpha_'+c]) for r in rows)/report['total_active_components'] for c in ('add','rec','dom')}
report['tissue_prior_ranges'] = {c: [min(float(r['pi_'+c]) for r in rows), max(float(r['pi_'+c]) for r in rows)] for c in ('add','rec','dom')}
out_files = list((root/'logs').glob('*.out'))
report['logged_gene_starts'] = 0
report['logs_with_weighted_fits'] = 0
report['observed_remote_iteration_paths'] = []
timestamps = []
for file in out_files:
    text = file.read_text(errors='replace')
    starts = re.findall(r'\[([^]]+)\] iteration 1, chunk \d+: ', text)
    timestamps += starts
    report['logged_gene_starts'] += len(starts)
    report['logs_with_weighted_fits'] += bool('Running weighted SuSiE-mix' in text)
    report['observed_remote_iteration_paths'] += re.findall(r'Logging to ([^\r\n]+/iteration_\d+)/temp_plink', text)
report['observed_remote_iteration_paths'] = sorted(set(report['observed_remote_iteration_paths']))
report['gene_start_time_range_cluster_local'] = [min(timestamps), max(timestamps)]
report['latest_log_file_mtime_utc'] = datetime.fromtimestamp(max(f.stat().st_mtime for f in logs+out_files), timezone.utc).isoformat()
target = Path('tmp/em_run_49011158_audit.json')
target.write_text(json.dumps(report, indent=2), encoding='utf-8')
print(json.dumps(report, indent=2))
