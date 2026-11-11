"""Deterministic longest-observed-first ordering; never changes batch membership."""
import math
from pathlib import Path


def collect(branch, mb, genes, stages):
    # Use the latest successful observation per gene. Attempt + receipt times
    # measure the entire R worker, including input/output, without reading fits.
    wanted = set(genes)
    records = {}
    for row in reversed(stages):
        d = branch / row['directory']
        receipts = d / 'receipts'
        if not receipts.is_dir():
            continue
        starts = {}
        for p in (d / 'attempts').glob('*.json'):
            a = mb.read(p)
            g = a.get('gene')
            if g in wanted and int(a.get('stage', -1)) == int(row['stage']):
                starts.setdefault(g, []).append((float(a['started']), p.name))
        for g in sorted(wanted):
            p = receipts / (g + '.json')
            if not p.exists() or g not in starts:
                continue
            r = mb.read(p)
            if r.get('gene') != g:
                raise RuntimeError('Runtime receipt gene mismatch: ' + str(p))
            end = float(r['completed'])
            eligible = [(t, name) for t, name in starts[g] if math.isfinite(t) and t < end]
            if not eligible or not math.isfinite(end):
                continue
            start, name = max(eligible)
            elapsed = end - start
            if not 0 < elapsed <= 24 * 3600:
                continue
            records[g] = dict(gene=g, seconds=elapsed, source_iteration=int(row['iteration']),
                              source_stage=int(row['stage']), attempt=name,
                              receipt=(Path(row['directory']) / 'receipts' / p.name).as_posix())
        wanted.difference_update(records)
        if not wanted:
            break
    # Unknown genes use the observed median, preserving the seeded order on ties.
    values = sorted(r['seconds'] for r in records.values())
    n = len(values)
    fallback = (values[(n-1)//2] + values[n//2]) / 2 if n else 0.0
    rows = []
    for g in genes:
        r = records.get(g, dict(gene=g, seconds=fallback, source_iteration='', source_stage='', attempt='', receipt=''))
        r = dict(r, observed=g in records)
        rows.append(r)
    return rows


def order(genes, estimates):
    return sorted(genes, key=lambda g: -estimates.get(g, 0.0))


def pending_order(branch, mb, genes):
    cfg = mb.read(branch / 'config.json')
    rel = cfg.get('runtime_order_file')
    if not rel:
        return genes
    p = Path(rel)
    if p.is_absolute() or '..' in p.parts or p.parts[0] != 'plan':
        raise RuntimeError('Invalid runtime-order path')
    rows = mb.csvrows(branch / p)
    if len({r['gene'] for r in rows}) != len(rows):
        raise RuntimeError('Duplicate runtime-order gene')
    estimates = {r['gene']: float(r['seconds']) for r in rows}
    if any(not math.isfinite(v) or v < 0 for v in estimates.values()) or not set(genes) <= set(estimates):
        raise RuntimeError('Invalid/incomplete runtime estimates')
    return order(genes, estimates)
