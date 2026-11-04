"""Count variants in all candidate cis windows in one read-only BIM pass.

Uses O(number of query boundaries) RAM; no PLINK files or genotypes are changed.
Counts presence only, before donor-specific MAF/HWE/other fitting QC.
"""
import argparse, bisect, csv, itertools, json
from collections import defaultdict
from pathlib import Path

def count_windows(bim, queries):
    edges = defaultdict(set)
    for q in queries:
        edges[q['chromosome']].update((q['lo']-1, q['hi']))
    edges={chrom:sorted(values) for chrom,values in edges.items()}
    hist={chrom:[0]*(len(values)+1) for chrom,values in edges.items()}
    total=0
    with Path(bim).open() as f:
        for line in f:
            fields=line.split()
            if len(fields)!=6: raise ValueError('Invalid BIM row at line {}'.format(total+1))
            chrom=fields[0].removeprefix('chr') if hasattr(str,'removeprefix') else fields[0].replace('chr','',1)
            if chrom in edges:
                hist[chrom][bisect.bisect_left(edges[chrom],int(fields[3]))]+=1
            total+=1
            if total%10000000==0: print('Read {:,} BIM variants'.format(total),flush=True)
    cumulative={c:dict(zip(edges[c],itertools.accumulate(hist[c]))) for c in edges}
    counts=[cumulative[q['chromosome']][q['hi']]-cumulative[q['chromosome']][q['lo']-1] for q in queries]
    return counts,total

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('bundle',type=Path); parser.add_argument('bim',type=Path)
    args=parser.parse_args()
    with (args.bundle/'checks/mapping_checks.csv').open() as f: rows=list(csv.DictReader(f))
    if not rows or any(r['status']!='PASS' for r in rows): raise ValueError('Run and pass the metadata preflight first.')
    queries=[]
    for r in rows:
        tss=int(r['start'] if r['strand']=='+' else r['end'])
        queries.append(dict(chromosome=r['chromosome'],lo=max(0,tss-500000),hi=tss+500000))
    counts,total=count_windows(args.bim,queries)
    for r,q,n in zip(rows,queries,counts):
        r.update(cis_start=q['lo'],cis_end=q['hi'],bim_variants=n,genotype_window_status='PASS' if n else 'NO_VARIANTS')
    dest=args.bundle/'checks/genotype_windows.csv'
    with dest.open('w',newline='') as f:
        w=csv.DictWriter(f,list(rows[0])); w.writeheader();w.writerows(rows)
    result={'candidate_ids':len(rows),'windows_with_variants':sum(n>0 for n in counts),
            'windows_without_variants':sum(n==0 for n in counts),'bim_variants_read':total,
            'limitation':'Variant presence before SNP-only, donor-specific MAF, HWE, sample and tissue QC.'}
    (args.bundle/'checks/genotype_summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    if not all(counts): raise SystemExit('Some candidate cis windows have no variants; inspect genotype_windows.csv.')

if __name__=='__main__': main()
