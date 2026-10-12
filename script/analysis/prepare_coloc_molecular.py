#!/usr/bin/env python3
"""Freeze same-gene pQTL/eQTL comparisons before downloading association data.

Uses the frozen screening cases, not association results or H4, for selection.
UKB-PPP requires synapseclient and an existing login or SYNAPSE_AUTH_TOKEN.
"""
import argparse
import csv
import json
import os
import pathlib
import re

from download_coloc_gwas import atomic_json, checksum, download

ROOT = pathlib.Path(__file__).resolve().parent
HGNC = 'https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt'


def synapse_client():
    try:
        import synapseclient
    except ImportError as e:
        raise RuntimeError('UKB-PPP needs pip install synapseclient and a Synapse login') from e
    syn = synapseclient.Synapse()
    token = os.environ.get('SYNAPSE_AUTH_TOKEN')
    syn.login(**({'authToken': token} if token else {}), silent=True)
    return syn


def gene_ids(path, genes):
    """Exact approved symbols/Ensembl IDs only; ambiguous aliases never guessed."""
    result = {}
    with open(path, encoding='utf8') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if row['status'] != 'Approved' or not row['ensembl_gene_id']:
                continue
            for key in (row['symbol'], row['ensembl_gene_id']):
                if key in genes:
                    result.setdefault(key, set()).add(row['ensembl_gene_id'])
    return {k: next(iter(v)) for k, v in result.items() if len(v) == 1}


def eqtl_sources(cases, genes, datasets):
    sources, audit = [], []
    for case in cases:
        selected = [d for d in datasets if case['tissue'] in d['case_tissues']]
        if not selected:
            audit.append(dict(case_id=case['case_id'], group='eqtl_replication', status='no_prespecified_tissue_dataset'))
        for d in selected:
            gene = case['gene']
            if gene not in genes:
                audit.append(dict(case_id=case['case_id'], group='eqtl_replication', status='unmapped_gene', gene=gene))
                continue
            if 'gtex' in d['study'].lower():
                raise ValueError('GTEx cannot be an independent replication dataset')
            s = dict(d, trait_id=f"eqtl_{d['dataset_id']}_{genes[gene]}", target_gene=gene,
                     molecular_trait_id=genes[gene], trait=gene+' expression in '+d['external_tissue'],
                     format='eqtl_catalogue_susie', trait_group='eqtl_replication', build='GRCh38',
                     bf_scale='natural_log_absolute', snp_coverage='full_finemapped_region',
                     sample_overlap='distinct_cohort_design', trait_type='quant', columns={})
            sources.append(s)
            audit.append(dict(case_id=case['case_id'], group='eqtl_replication', status='selected', trait_id=s['trait_id']))
    # Same gene/tissue group may occur in multiple cases, but download each once.
    return list({s['trait_id']: s for s in sources}.values()), audit


def pqtl_sources(cases, children, syn):
    if not children:
        raise ValueError('Empty UKB-PPP file listing; release/access must be checked')
    genes = {c['gene'].upper(): c['gene'] for c in cases}
    sources = []
    for child in children:
        name = child['name']
        m = re.fullmatch(r'([^_]+)_([^_]+)_(OID\d+)_(v\d+)_(.+)\.tar', name, re.I)
        if not m or m[1].upper() not in genes:
            continue
        gene = genes[m[1].upper()]
        entity = syn.get(child['id'], downloadFile=False)
        version = entity['versionNumber']
        sources.append(dict(trait_id='pqtl_'+child['id'], trait=gene+' protein ('+m[3]+')',
            trait_group='cis_pqtl', target_gene=gene, external_tissue='plasma',
            study='UKB-PPP European discovery', consortium='UKB-PPP', ancestry='EUR',
            source_page='https://www.synapse.org/Synapse:syn51365303',
            sample_overlap='distinct_cohort_design', assay_id=m[3], uniprot=m[2],
            protein_assay_caveat='Protein-altering variants may affect assay binding; validate across platforms.',
            synapse_id=child['id'], synapse_version=version,
            url=f"synapse://{child['id']}/version/{version}", format='ukb_ppp_tar',
            archive_name=name, build='GRCh38', trait_type='quant', sd_y=1,
            primary_eligible=True, columns=dict(chromosome='CHROM', position='GENPOS',
                effect_allele='ALLELE1', other_allele='ALLELE0', beta='BETA', se='SE',
                n='N', eaf='A1FREQ'),
            notes='Full additive REGENIE statistics, no association P filter. N is per variant; European discovery subset, not total UKB-PPP.'))
    covered = {s['target_gene'] for s in sources}
    audit = [dict(case_id=c['case_id'], group='cis_pqtl', gene=c['gene'],
                  status='selected' if c['gene'] in covered else 'no_exact_protein_match') for c in cases]
    return sources, audit


def prepare(job, cache, include_eqtl=True, include_pqtl=True):
    receipt = job/'molecular_selection.json'
    if receipt.exists():
        old = json.loads(receipt.read_text())
        if old['sources_sha256'] != checksum(job/'sources.json') or old['cases_sha256'] != checksum(job/'cases.csv'):
            raise ValueError('Frozen molecular selection changed; prepare a new job')
        return old
    if list((job/'gwas').glob('*.json')) or list(job.glob('analysis_*')):
        raise ValueError('Add molecular sources only to a new job, before downloads/analysis')
    with (job/'cases.csv').open(encoding='utf8') as f:
        cases = list(csv.DictReader(f))
    registry = json.loads((job/'sources.json').read_text())
    sources, audit, provenance = [], [], {}
    if include_eqtl:
        try:
            path, provenance['hgnc'] = download(dict(url=HGNC), cache)
            genes = gene_ids(path, {c['gene'] for c in cases})
            config = json.loads((ROOT/'coloc_eqtl_datasets.json').read_text())
            add, rows = eqtl_sources(cases, genes, config['datasets'])
            sources += add; audit += rows
            provenance['eqtl_config'] = config
        except Exception as e:
            audit.append(dict(group='eqtl_replication', status='preparation_failed', detail=str(e)))
    if include_pqtl:
        try:
            syn = synapse_client()
            children = list(syn.getChildren('syn51365303', includeTypes=['file']))
            add, rows = pqtl_sources(cases, children, syn)
            sources += add; audit += rows
            provenance['ukb_ppp_folder'] = 'syn51365303'
        except Exception:
            # Credentials and service exception text must not enter manifests.
            audit.append(dict(group='cis_pqtl', status='preparation_failed',
                detail='UKB-PPP metadata unavailable. Check synapseclient, login and access to syn51365303; prepare a new job to retry.'))
    ids = [s['trait_id'] for s in registry['sources']+sources]
    if len(ids) != len(set(ids)):
        raise ValueError('Duplicate trait IDs')
    registry['sources'] += sources
    atomic_json(job/'sources.json', registry)
    info = dict(sources_sha256=checksum(job/'sources.json'), cases_sha256=checksum(job/'cases.csv'),
                provenance=provenance, selection_audit=audit, n_molecular_sources=len(sources))
    atomic_json(receipt, info)
    return info


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--job-dir', required=True); ap.add_argument('--cache', required=True)
    ap.add_argument('--no-eqtl', action='store_true'); ap.add_argument('--no-pqtl', action='store_true')
    args = ap.parse_args()
    info = prepare(pathlib.Path(args.job_dir), pathlib.Path(args.cache), not args.no_eqtl, not args.no_pqtl)
    print('Frozen molecular comparisons:', info['n_molecular_sources'])
    failed = [r for r in info['selection_audit'] if r['status']=='preparation_failed']
    for row in failed: print(row['group']+': '+row['detail'])
    return 1 if failed else 0


if __name__ == '__main__':
    raise SystemExit(main())
