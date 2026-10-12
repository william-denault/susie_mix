#!/usr/bin/env python3
"""Pinned GWAS downloads and streaming regional extraction. Python >=3.9, stdlib only.

Input regions are 1-based GRCh38 intervals from prepare_em_gwas_coloc(). No GWAS
P-value filtering. Original files, checksums and source definitions are retained.
"""
import argparse
import bisect
import collections
import contextlib
import csv
import gzip
import hashlib
import io
import json
import math
import pathlib
import sys
import time
import urllib.request
import zipfile

FIELDS = ['chromosome', 'position', 'effect_allele', 'other_allele', 'beta',
          'se', 'p', 'n', 'eaf', 'log_bf', 'info']
COMPLEMENT = str.maketrans('ACGT', 'TGCA')


def checksum(path, algorithm='sha256'):
    h = hashlib.new(algorithm)
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def atomic_json(path, data):
    tmp = pathlib.Path(str(path) + '.part')
    tmp.write_text(json.dumps(data, indent=2) + '\n', encoding='utf8')
    tmp.replace(path)


def download(source, cache):
    """A cache entry is immutable; validate it before reuse. Never trust .part files."""
    key = hashlib.sha256((source['url'] + source.get('post', '')).encode()).hexdigest()[:20]
    dest = cache / (key + '.data')
    receipt = cache / (key + '.json')
    cache.mkdir(parents=True, exist_ok=True)
    if dest.exists() and receipt.exists():
        old = json.loads(receipt.read_text())
        if checksum(dest) != old['sha256']:
            raise ValueError('Cached file changed: ' + str(dest))
        if source.get('md5') and checksum(dest, 'md5') != source['md5']:
            raise ValueError('Cached file disagrees with pinned source MD5')
        return dest, old
    tmp = pathlib.Path(str(dest) + '.part')
    request = urllib.request.Request(source['url'],
        data=source['post'].encode() if source.get('post') else None,
        headers={'User-Agent': 'susie-mix-coloc/1.0'})
    for attempt in range(3):
        try:
            with urllib.request.urlopen(request, timeout=120) as response, open(tmp, 'wb') as f:
                expected = response.headers.get('Content-Length')
                total = 0
                last = time.monotonic()
                for block in iter(lambda: response.read(1024 * 1024), b''):
                    f.write(block)
                    total += len(block)
                    if time.monotonic() - last > 30:
                        print(f"  downloaded {total / 1e6:.0f} MB", flush=True)
                        last = time.monotonic()
                if expected and total != int(expected):
                    raise IOError('Truncated download')
            if source.get('md5') and checksum(tmp, 'md5') != source['md5']:
                raise ValueError('Download disagrees with pinned source MD5')
            tmp.replace(dest)
            info = dict(url=source['url'], post=source.get('post'), bytes=total,
                        sha256=checksum(dest), md5=checksum(dest, 'md5'),
                        downloaded_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()))
            atomic_json(receipt, info)
            return dest, info
        except (OSError, ValueError):
            if attempt == 2:
                raise
            time.sleep(2 ** attempt)


class Regions:
    def __init__(self, rows):
        self.intervals = collections.defaultdict(list)
        for row in rows:
            self.intervals[str(row['chromosome']).removeprefix('chr')].append(
                (int(row['start']), int(row['end'])))
        for chrom, pairs in self.intervals.items():
            merged = []
            for lo, hi in sorted(pairs):
                if merged and lo <= merged[-1][1] + 1:
                    merged[-1] = (merged[-1][0], max(hi, merged[-1][1]))
                else:
                    merged.append((lo, hi))
            self.intervals[chrom] = merged
        self.starts = {c: [p[0] for p in ps] for c, ps in self.intervals.items()}

    def overlaps(self, chrom, lo, hi):
        pairs = self.intervals.get(chrom, [])
        i = bisect.bisect_right(self.starts.get(chrom, []), hi) - 1
        return i >= 0 and pairs[i][1] >= lo


class ChainMap:
    """Map SNV bases only; reject nonunique mappings and complement minus strands.

    UCSC chain t = source build, q = destination build. Keep only blocks that
    intersect requested destination regions; source buckets make streaming cheap.
    """
    def __init__(self, path, regions):
        self.bins = collections.defaultdict(list)
        candidate_bins = set()
        for tc, block in self.blocks(path):
            lo, hi, qc, qp, qs, strand = block
            qlo = qp if strand == '+' else qs - qp - (hi - lo)
            if regions.overlaps(qc, qlo + 1, qlo + hi - lo):
                candidate_bins.update((tc, b) for b in range(lo // 100000, (hi - 1) // 100000 + 1))
        # Include competing mappings OUTSIDE the requested target regions too.
        # Otherwise a multiply mapped source base could appear spuriously unique.
        for tc, block in self.blocks(path):
            lo, hi = block[:2]
            for bucket in range(lo // 100000, (hi - 1) // 100000 + 1):
                if (tc, bucket) in candidate_bins:
                    self.bins[(tc, bucket)].append(block)

    @staticmethod
    def blocks(path):
        with gzip.open(path, 'rt') as f:
            for line in f:
                x = line.split()
                if not x:
                    continue
                if x[0] == 'chain':
                    _, _, tc, _, ts, tp, _, qc, qs, strand, qp, _, _ = x
                    if ts != '+':
                        raise ValueError('Unsupported chain source strand')
                    tc, qc = tc.removeprefix('chr'), qc.removeprefix('chr')
                    tp, qp, qs = int(tp), int(qp), int(qs)
                    continue
                size = int(x[0])
                yield tc, (tp, tp + size, qc, qp, qs, strand)
                tp += size + (int(x[1]) if len(x) == 3 else 0)
                qp += size + (int(x[2]) if len(x) == 3 else 0)

    def map(self, chrom, position):
        pos = position - 1
        hits = []
        for lo, hi, qc, qp, qs, strand in self.bins.get((chrom, pos // 100000), []):
            if lo <= pos < hi:
                offset = qp + pos - lo
                hits.append((qc, offset + 1 if strand == '+' else qs - offset, strand))
        return hits[0] if len(hits) == 1 else None


@contextlib.contextmanager
def text_file(path, source):
    with contextlib.ExitStack() as stack:
        f = stack.enter_context(open(path, 'rb'))
        magic = f.read(2)
        f.seek(0)
        if magic == b'PK':
            archive = stack.enter_context(zipfile.ZipFile(f))
            member = source['archive_member']  # Never guess among archive members.
            f = stack.enter_context(archive.open(member))
        elif magic == b'\x1f\x8b':
            f = stack.enter_context(gzip.GzipFile(fileobj=f))
        yield io.TextIOWrapper(f, encoding='utf8', errors='replace')


def number(x):
    try:
        return float(x)
    except (ValueError, TypeError):
        return float('nan')


def normalize(source, path, output, regions, mapper=None):
    counts = collections.Counter()
    tmp = pathlib.Path(str(output) + '.part')
    columns = source['columns']
    with text_file(path, source) as f, gzip.open(tmp, 'wt', newline='') as out:
        # PGC has ## metadata followed by CHROM or #CHROM. Do not skip its header.
        header = next(line for line in f if line.strip() and not line.startswith('##')).split()
        missing = set(columns.values()) - set(header)
        if missing:
            raise ValueError('Changed source schema; missing ' + ','.join(sorted(missing)))
        indexes = {k: header.index(v) for k, v in columns.items()}
        writer = csv.DictWriter(out, fieldnames=FIELDS, delimiter='\t')
        writer.writeheader()
        for line in f:
            counts['source_rows'] += 1
            fields = line.split()
            if len(fields) != len(header):
                counts['malformed_rows'] += 1
                continue
            row = {key: fields[index] for key, index in indexes.items()}
            chrom = row['chromosome'].removeprefix('chr')
            pos = number(row['position'])
            if chrom not in map(str, range(1, 23)) or not math.isfinite(pos) or pos < 1 or pos != int(pos):
                counts['invalid_coordinate_or_non_autosome'] += 1
                continue
            ea, oa = row['effect_allele'].upper(), row['other_allele'].upper()
            if ea not in ('A', 'C', 'G', 'T') or oa not in ('A', 'C', 'G', 'T') or ea == oa:
                counts['not_biallelic_snv'] += 1
                continue
            pos = int(pos)
            if mapper:
                hit = mapper.map(chrom, pos)
                if hit is None:
                    counts['outside_target_or_unmappable'] += 1
                    continue
                chrom, pos, strand = hit
                if strand == '-':
                    ea, oa = ea.translate(COMPLEMENT), oa.translate(COMPLEMENT)
            if not regions.overlaps(chrom, pos, pos):
                continue
            values = {k: number(row.get(k)) for k in FIELDS[4:]}
            if 'log_bf' in columns:
                valid = math.isfinite(values['log_bf'])
            else:
                valid = math.isfinite(values['beta']) and math.isfinite(values['se']) and values['se'] > 0
            if not valid:
                counts['invalid_effect_in_region'] += 1
                continue
            if 'n' in columns:
                values['n'] *= source.get('n_multiplier', 1)
            else:
                values['n'] = number(source.get('n_fallback'))
            writer.writerow(dict(chromosome=chrom, position=pos, effect_allele=ea,
                                 other_allele=oa, **{k: v if math.isfinite(v) else 'NA' for k, v in values.items()}))
            counts['regional_rows'] += 1
    tmp.replace(output)
    return dict(counts)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--job-dir', required=True)
    ap.add_argument('--sources', default=str(pathlib.Path(__file__).with_name('coloc_gwas_sources.json')))
    ap.add_argument('--cache', required=True)
    ap.add_argument('--traits', help='Comma-separated trait IDs; default all 32')
    ap.add_argument('--chain-grch37', help='Optional local hg19ToHg38.over.chain.gz')
    ap.add_argument('--chain-grch36', help='Optional local hg18ToHg38.over.chain.gz')
    args = ap.parse_args(argv)
    job, cache = pathlib.Path(args.job_dir), pathlib.Path(args.cache)
    region_path = job / 'regions.csv'
    with region_path.open(encoding='utf8') as f:
        regions = Regions(list(csv.DictReader(f)))
    registry = json.loads(pathlib.Path(args.sources).read_text())
    sources = registry['sources']
    if args.traits:
        selected = set(args.traits.split(','))
        if selected - {s['trait_id'] for s in sources}:
            ap.error('Unknown trait ID')
        sources = [s for s in sources if s['trait_id'] in selected]
    dest = job / 'gwas'
    dest.mkdir(exist_ok=True)
    chain_maps = {}
    statuses = []
    for source in sources:
        trait = source['trait_id']
        print('Preparing ' + trait, flush=True)
        receipt_path = dest / (trait + '.json')
        output = dest / (trait + '.tsv.gz')
        local_chain = {'GRCh37': args.chain_grch37, 'GRCh36': args.chain_grch36}.get(source['build'])
        chain_hash = checksum(local_chain) if local_chain else 'automatic_UCSC_chain'
        fingerprint = hashlib.sha256((json.dumps(source, sort_keys=True) + checksum(region_path) + chain_hash).encode()).hexdigest()
        try:
            if receipt_path.exists() and output.exists():
                old = json.loads(receipt_path.read_text())
                if old.get('fingerprint') == fingerprint and old.get('output_sha256') == checksum(output):
                    statuses.append(dict(trait_id=trait, status='cached', message=''))
                    continue
                raise ValueError('Existing regional output has different settings/content; prepare a new job directory')
            mapper = None
            chain_receipt = None
            if source['build'] != 'GRCh38':
                build = {'GRCh36': 'hg18', 'GRCh37': 'hg19'}[source['build']]
                if build not in chain_maps:
                    if local_chain:
                        chain_path = pathlib.Path(local_chain).resolve()
                        chain_info = dict(url=chain_path.as_uri(), sha256=checksum(chain_path), md5=checksum(chain_path, 'md5'))
                    else:
                        chain_source = dict(url=f'https://hgdownload.soe.ucsc.edu/goldenPath/{build}/liftOver/{build}ToHg38.over.chain.gz')
                        chain_path, chain_info = download(chain_source, cache)
                    chain_maps[build] = (ChainMap(chain_path, regions), chain_info)
                mapper, chain_receipt = chain_maps[build]
            raw, receipt = download(source, cache)
            counts = normalize(source, raw, output, regions, mapper)
            atomic_json(receipt_path, dict(source=source, download=receipt, chain=chain_receipt,
                        fingerprint=fingerprint, output_sha256=checksum(output), counts=counts, build='GRCh38'))
            statuses.append(dict(trait_id=trait, status='ready', message=str(counts)))
        except Exception as error:
            statuses.append(dict(trait_id=trait, status='error', message=str(error)))
            print(f'  {trait}: {error}', file=sys.stderr, flush=True)
        atomic_json(dest / 'download_status.json', statuses)
    atomic_json(dest / 'download_status.json', statuses)
    return 1 if any(x['status'] == 'error' for x in statuses) else 0


if __name__ == '__main__':
    sys.exit(main())
