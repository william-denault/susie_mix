"""Full-region molecular summary-data adapters; no P-value/lead-SNP filtering."""
import collections
import contextlib
import csv
import gzip
import io
import math
import pathlib
import re
import tarfile


@contextlib.contextmanager
def bounded_writers(paths, fields, limit=32):
    """Append gzip members with at most limit descriptors for large gene panels."""
    for path in paths.values():
        with gzip.open(path, 'wt', newline='') as f:
            csv.DictWriter(f, fieldnames=fields, delimiter='\t').writeheader()
    opened = collections.OrderedDict()
    def writer(gene):
        if gene not in opened:
            if len(opened) >= limit:
                _, (handle, _) = opened.popitem(last=False)
                handle.close()
            handle = gzip.open(paths[gene], 'at', newline='')
            opened[gene] = (handle, csv.DictWriter(handle, fieldnames=fields, delimiter='\t'))
        opened.move_to_end(gene)
        return opened[gene][1]
    try:
        yield writer
    finally:
        for handle, _ in opened.values():
            handle.close()


@contextlib.contextmanager
def ukb_ppp_lines(path):
    """Stream tar members without extracting paths supplied by an archive."""
    with tarfile.open(path, 'r:*') as archive:
        def lines():
            first = None
            found = 0
            for member in archive:
                if not member.isfile() or not member.name.endswith('.gz'):
                    continue
                with archive.extractfile(member) as raw, gzip.GzipFile(fileobj=raw) as compressed:
                    text = io.TextIOWrapper(compressed)
                    header = next(text).strip()
                    fields = header.split()
                    if not {'CHROM', 'GENPOS', 'ALLELE0', 'ALLELE1', 'BETA', 'SE', 'N', 'TEST'} <= set(fields):
                        continue
                    if first is None:
                        first = header; yield header+'\n'
                    elif header != first:
                        raise ValueError('Inconsistent REGENIE headers within UKB-PPP archive')
                    found += 1
                    test = fields.index('TEST')
                    for line in text:
                        row = line.split()
                        if len(row) == len(fields) and row[test] == 'ADD':
                            yield line
            if not found:
                raise ValueError('No full REGENIE chromosome files in UKB-PPP archive')
        yield lines()


def normalize_eqtl_batch(sources, lbf_path, cs_path, dest):
    """Scan each Catalogue dataset once for all requested genes.

    CS membership identifies the SuSiE component only. Every SNP in its original
    fine-mapped region is exported, including weak and non-CS variants. Catalogue
    lbf_variable columns are natural-log SNP BFs; cs_log10bf is never used.
    """
    wanted = {s['molecular_trait_id']: s for s in sources}
    components = collections.defaultdict(dict)
    with gzip.open(cs_path, 'rt') as f:
        reader = csv.DictReader(f, delimiter='\t')
        if not {'molecular_trait_id', 'cs_id', 'cs_min_r2', 'region'} <= set(reader.fieldnames):
            raise ValueError('Changed Catalogue credible-set schema')
        for row in reader:
            gene = row['molecular_trait_id']
            if gene not in wanted:
                continue
            try:
                purity = float(row['cs_min_r2'])
            except ValueError:
                purity = float('nan')
            if row.get('cs_size') == '1':
                purity = 1.0
            if not math.isfinite(purity) or purity < .25:
                continue
            m = re.fullmatch(re.escape(gene)+r'_L(\d+)', row['cs_id'])
            if not m:
                raise ValueError('Cannot identify published SuSiE component: '+row['cs_id'])
            component = int(m[1])
            previous = components[gene].setdefault(component, row['region'])
            if previous != row['region']:
                raise ValueError('Multiple fine-mapped regions for a gene/component')
    counts = {s['trait_id']: dict(regional_rows=0, n_components=len(components[s['molecular_trait_id']])) for s in sources}
    fields = ['component','chromosome','position','effect_allele','other_allele','log_bf']
    # Temporary outputs are published only after an entire input scan succeeds.
    paths = {gene: dest/(source['trait_id']+'.tsv.gz.part') for gene, source in wanted.items()}
    with bounded_writers(paths, fields) as writer:
        with gzip.open(lbf_path, 'rt') as f:
            reader = csv.DictReader(f, delimiter='\t')
            if not {'molecular_trait_id','region','variant','chromosome','position'} <= set(reader.fieldnames):
                raise ValueError('Changed Catalogue log-BF schema')
            for gene, cs in components.items():
                if any('lbf_variable'+str(c) not in reader.fieldnames for c in cs):
                    raise ValueError('Credible-set component missing from log-BF file')
            for row in reader:
                gene = row['molecular_trait_id']
                if gene not in wanted or not components[gene]:
                    continue
                # Unsupported variants remain represented for the coverage denominator.
                bits = row['variant'].split('_')
                if len(bits) != 4:
                    raise ValueError('Invalid Catalogue variant ID')
                chrom, pos, ref, alt = bits
                chrom = chrom.removeprefix('chr')
                if chrom != row['chromosome'].removeprefix('chr') or int(pos) != int(row['position']):
                    raise ValueError('Catalogue variant coordinates disagree')
                for component, region in components[gene].items():
                    if row['region'] != region:
                        raise ValueError('Published BF region does not match the saved credible set')
                    value = float(row['lbf_variable'+str(component)])
                    if not math.isfinite(value):
                        raise ValueError('Nonfinite published log BF')
                    writer(gene).writerow(dict(component=component,chromosome=chrom,position=pos,
                        effect_allele=alt,other_allele=ref,log_bf=value))
                    counts[wanted[gene]['trait_id']]['regional_rows'] += 1
    for source in sources:
        path = dest/(source['trait_id']+'.tsv.gz')
        pathlib.Path(str(path)+'.part').replace(path)
    return counts
