import csv
import gzip
import io
import json
import pathlib
import sys
import tarfile
import tempfile
import unittest
from unittest import mock

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import download_coloc_gwas as d
import coloc_molecular_io as m
import prepare_coloc_molecular as p


class MolecularTests(unittest.TestCase):
    def test_18_verified_biomarker_headers(self):
        sources = {s['trait_id']: s for s in json.loads((ROOT/'coloc_biomarker_sources.json').read_text())['sources']}
        fixtures = json.loads((ROOT/'tests/biomarker_header_fixtures.json').read_text())
        self.assertEqual(len(sources), 18)
        self.assertEqual(set(sources), {f['trait_id'] for f in fixtures})
        regions = d.Regions([dict(chromosome=c,start=1,end=300000000) for c in range(1,23)])
        with tempfile.TemporaryDirectory() as tmp:
            tmp = pathlib.Path(tmp)
            for f in fixtures:
                with self.subTest(trait=f['trait_id']):
                    raw=tmp/'raw'; raw.write_text(f['header']+'\n'+f['row']+'\n')
                    counts=d.normalize(sources[f['trait_id']],raw,tmp/'out.gz',regions)
                    self.assertEqual(counts['source_rows'],1)
                    self.assertEqual(counts.get('malformed_rows',0),0)
                    self.assertEqual(sum(counts.get(k,0) for k in ('regional_rows','not_biallelic_snv','invalid_effect_in_region')),1)
                    if f['trait_id']=='tsh':
                        with gzip.open(tmp/'out.gz','rt') as out: row=next(csv.DictReader(out,delimiter='\t'))
                        self.assertEqual((row['chromosome'],row['position']),('10','100000625'))
                        self.assertAlmostEqual(float(row['beta']),-.0006)

    def test_published_bfs_include_non_cs_and_unsupported_variant_mass(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp=pathlib.Path(tmp)
            cs=tmp/'cs.gz'; bf=tmp/'bf.gz'
            with gzip.open(cs,'wt') as f:
                f.write('molecular_trait_id\tcs_id\tcs_min_r2\tregion\nG1\tG1_L2\t.8\tchr1:1-500\nG2\tG2_L1\t.9\tchr1:1-500\n')
            with gzip.open(bf,'wt') as f:
                f.write('molecular_trait_id\tregion\tvariant\tchromosome\tposition\tlbf_variable1\tlbf_variable2\n')
                f.write('G1\tchr1:1-500\tchr1_10_A_G\t1\t10\t99\t20\n')
                f.write('G1\tchr1:1-500\tchr1_20_A_G\t1\t20\t99\t-3\n')
                f.write('G1\tchr1:1-500\tchr1_30_A_AT\t1\t30\t99\t4\n')
                f.write('G2\tchr1:1-500\tchr1_10_A_G\t1\t10\t-2\t99\n')
            sources=[dict(trait_id='one',molecular_trait_id='G1'),dict(trait_id='two',molecular_trait_id='G2')]
            counts=m.normalize_eqtl_batch(sources,bf,cs,tmp)
            self.assertEqual(counts['one']['regional_rows'],3)
            with gzip.open(tmp/'one.tsv.gz','rt') as f: rows=list(csv.DictReader(f,delimiter='\t'))
            self.assertEqual([float(r['log_bf']) for r in rows],[20,-3,4])
            self.assertTrue(all(r['component']=='2' for r in rows))
            self.assertEqual(rows[-1]['effect_allele'],'AT')
            with gzip.open(tmp/'two.tsv.gz','rt') as f: row=next(csv.DictReader(f,delimiter='\t'))
            self.assertEqual(float(row['log_bf']),-2)
            # Actual downloader integration: both genes share one dataset/cache.
            job=tmp/'job';job.mkdir();(job/'regions.csv').write_text('chromosome,start,end\n1,1,100\n')
            registry=job/'sources.json'
            frozen=[dict(s,format='eqtl_catalogue_susie',url=bf.as_uri(),cs_url=cs.as_uri(),build='GRCh38') for s in sources]
            registry.write_text(json.dumps(dict(sources=frozen)))
            args=['--job-dir',str(job),'--sources',str(registry),'--cache',str(tmp/'cache')]
            self.assertEqual(d.main(args),0)
            self.assertEqual(d.main(args),0)
            self.assertEqual([s['status'] for s in json.loads((job/'gwas/download_status.json').read_text())],['cached','cached'])
            self.assertEqual(json.loads((job/'gwas/one.json').read_text())['counts']['regional_rows'],3)
            # Truncated components/changed schemas must fail before publishing.
            with gzip.open(cs,'wt') as f: f.write('molecular_trait_id\tcs_id\tcs_min_r2\tregion\nG1\tG1_L3\t.9\tchr1:1-500\n')
            with self.assertRaisesRegex(ValueError,'component missing'):
                m.normalize_eqtl_batch(sources,bf,cs,tmp)

    def test_ukb_tar_does_not_extract_paths_or_keep_other_tests(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp=pathlib.Path(tmp); archive=tmp/'protein.tar'
            header='CHROM GENPOS ALLELE0 ALLELE1 BETA SE N A1FREQ TEST\n'
            with tarfile.open(archive,'w') as tar:
                for name,pos in [('../outside.regenie.gz',10),('chr2.regenie.gz',20)]:
                    data=gzip.compress((header+f'1 {pos} A G .2 .1 1000 .2 ADD\n1 30 A G 99 .1 1000 .2 DOM\n').encode())
                    member=tarfile.TarInfo(name); member.size=len(data);tar.addfile(member,io.BytesIO(data))
            source=dict(format='ukb_ppp_tar',columns=dict(chromosome='CHROM',position='GENPOS',effect_allele='ALLELE1',
                other_allele='ALLELE0',beta='BETA',se='SE',n='N',eaf='A1FREQ'))
            counts=d.normalize(source,archive,tmp/'out.gz',d.Regions([dict(chromosome=1,start=1,end=100)]))
            self.assertEqual(counts['regional_rows'],2)
            self.assertFalse((tmp.parent/'outside.regenie.gz').exists())

    def test_gene_tissue_selection_and_no_gtex(self):
        cases=[dict(case_id='a',gene='GENE',tissue='Blood'),dict(case_id='b',gene='GENE',tissue='Liver')]
        dataset=dict(dataset_id='QTD1',study='Separate cohort',case_tissues=['Blood'],external_tissue='blood')
        sources,audit=p.eqtl_sources(cases,{'GENE':'ENSG1'},[dataset])
        self.assertEqual(len(sources),1);self.assertEqual(sources[0]['target_gene'],'GENE')
        self.assertEqual(audit[-1]['status'],'no_prespecified_tissue_dataset')
        with self.assertRaisesRegex(ValueError,'GTEx'):
            p.eqtl_sources(cases,{'GENE':'ENSG1'},[dict(dataset,study='GTEx_v8')])
        class Syn:
            def get(self, id, downloadFile): return dict(versionNumber=3)
        sources,audit=p.pqtl_sources(cases,[dict(id='syn1',name='GENE_P12345_OID123_v1_Neurology.tar'),
            dict(id='syn2',name='OTHER_P12345_OID999_v1_Neurology.tar')],Syn())
        self.assertEqual(len(sources),1); self.assertEqual(sources[0]['synapse_version'],3)
        self.assertEqual(sources[0]['columns']['effect_allele'],'ALLELE1')

    def test_molecular_selection_frozen_before_analysis(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp=pathlib.Path(tmp);(tmp/'gwas').mkdir()
            (tmp/'sources.json').write_text('{"sources": []}')
            (tmp/'cases.csv').write_text('case_id,gene,tissue\na,GENE,Blood\n')
            hgnc=tmp/'hgnc.tsv';hgnc.write_text('symbol\tstatus\tensembl_gene_id\nGENE\tApproved\tENSG000001\n')
            with mock.patch.object(p,'download',return_value=(hgnc,dict(sha256=d.checksum(hgnc)))):
                result=p.prepare(tmp,tmp/'cache',include_pqtl=False)
            self.assertEqual(result['n_molecular_sources'],1)
            self.assertEqual(p.prepare(tmp,tmp/'cache',include_pqtl=False),result)
            (tmp/'sources.json').write_text('{"sources": []}')
            with self.assertRaisesRegex(ValueError,'Frozen molecular selection changed'):
                p.prepare(tmp,tmp/'cache',include_pqtl=False)


if __name__=='__main__': unittest.main()
