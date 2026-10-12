"""python -m unittest discover -s script/analysis/tests -p test_download_coloc_gwas.py"""
import csv
import gzip
import importlib.util
import json
import pathlib
import tempfile
import unittest
import zipfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('download', ROOT / 'download_coloc_gwas.py')
d = importlib.util.module_from_spec(spec)
spec.loader.exec_module(d)


class DownloadTests(unittest.TestCase):
    def test_all_32_verified_headers(self):
        sources = {s['trait_id']: s for s in json.loads((ROOT / 'coloc_gwas_sources.json').read_text())['sources']}
        fixtures = json.loads(pathlib.Path(__file__).with_name('gwas_header_fixtures.json').read_text())
        self.assertEqual(set(sources), {f['trait_id'] for f in fixtures})
        self.assertEqual(len(sources), 32)
        regions = d.Regions([dict(chromosome=c, start=1, end=300000000) for c in range(1, 23)])
        with tempfile.TemporaryDirectory() as tmp:
            tmp = pathlib.Path(tmp)
            for fixture in fixtures:
                with self.subTest(trait=fixture['trait_id']):
                    source = sources[fixture['trait_id']]
                    text = fixture['header'] + '\n' + fixture['row'] + '\n'
                    raw = tmp / 'raw'
                    if 'archive_member' in source:
                        with zipfile.ZipFile(raw, 'w') as archive:
                            archive.writestr(source['archive_member'], text)
                    else:
                        with gzip.open(raw, 'wt') as f:
                            if fixture['trait_id'] in ['bipolar', 'schizophrenia']:
                                f.write('##source_metadata=ignored\n')
                            f.write(text)
                    out = tmp / 'out.gz'
                    counts = d.normalize(source, raw, out, regions)
                    self.assertEqual(counts['source_rows'], 1)
                    self.assertEqual(counts.get('malformed_rows', 0), 0)
                    # The first lipid row can be an indel and must be rejected.
                    self.assertEqual(counts.get('regional_rows', 0) + counts.get('not_biallelic_snv', 0) +
                                     counts.get('invalid_effect_in_region', 0), 1)
                    if fixture['trait_id'] == 'child_obesity':
                        # The first public row has African, but no European, results.
                        self.assertEqual(counts['invalid_effect_in_region'], 1)
                    if counts.get('regional_rows'):
                        with gzip.open(out, 'rt') as f:
                            result = list(csv.DictReader(f, delimiter='\t'))[0]
                        if fixture['trait_id'] == 'bipolar':
                            self.assertAlmostEqual(float(result['n']), 2*50981.48)
                        if fixture['trait_id'] == 'prostate_cancer':
                            self.assertEqual(result['effect_allele'], 'C')
                            self.assertAlmostEqual(float(result['beta']), -.3608)
                        if fixture['trait_id'] == 'pubertal_growth':
                            self.assertAlmostEqual(float(result['log_bf']), -2.8652)
                            self.assertEqual(result['beta'], 'NA')

    def test_chain_one_based_reverse_strand_and_boundaries(self):
        regions = d.Regions([dict(chromosome='2', start=101, end=120)])
        with tempfile.TemporaryDirectory() as tmp:
            path = pathlib.Path(tmp) / 'chain.gz'
            with gzip.open(path, 'wt') as f:
                f.write('chain 10 chr1 1000 + 10 20 chr2 1000 + 100 110 1\n10\n\n')
                f.write('chain 10 chr1 1000 + 30 40 chr2 1000 - 880 890 2\n10\n')
            mapper = d.ChainMap(path, regions)
            self.assertIsNone(mapper.map('1', 10))
            self.assertEqual(mapper.map('1', 11), ('2', 101, '+'))
            self.assertEqual(mapper.map('1', 20), ('2', 110, '+'))
            self.assertIsNone(mapper.map('1', 21))
            self.assertEqual(mapper.map('1', 31), ('2', 120, '-'))
            self.assertEqual(mapper.map('1', 40), ('2', 111, '-'))
            source = dict(columns=dict(chromosome='chr',position='pos',effect_allele='ea',other_allele='oa',beta='b',se='s'))
            raw = pathlib.Path(tmp) / 'raw'
            raw.write_text('chr pos ea oa b s\n1 31 A G .2 .1\n')
            output = pathlib.Path(tmp) / 'out.gz'
            d.normalize(source, raw, output, regions, mapper)
            with gzip.open(output, 'rt') as f:
                row = next(csv.DictReader(f, delimiter='\t'))
            self.assertEqual((row['effect_allele'], row['other_allele']), ('T', 'C'))
            self.assertEqual(float(row['beta']), .2)  # Complement is not an effect swap.
            with gzip.open(path, 'at') as f:
                f.write('\nchain 10 chr1 1000 + 10 20 chr3 1000 + 500 510 3\n10\n')
            ambiguous = d.ChainMap(path, regions)
            self.assertIsNone(ambiguous.map('1', 11)) # Competing mapping is outside target regions.

    def test_cache_checksum_and_incomplete_download(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = pathlib.Path(tmp)
            raw = tmp / 'source'
            raw.write_bytes(b'abc')
            source = dict(url=raw.as_uri(), md5=d.checksum(raw, 'md5'))
            cached, receipt = d.download(source, tmp / 'cache')
            same, again = d.download(source, tmp / 'cache')
            self.assertEqual((cached, receipt), (same, again))
            cached.write_bytes(b'changed')
            with self.assertRaisesRegex(ValueError, 'changed'):
                d.download(source, tmp / 'cache')

    def test_missing_column_fails_without_publishing_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            raw = pathlib.Path(tmp) / 'raw'
            raw.write_text('chr pos ea oa\n1 100 A G\n')
            output = pathlib.Path(tmp) / 'out.gz'
            with self.assertRaisesRegex(ValueError, 'schema'):
                d.normalize(dict(columns=dict(beta='BETA')), raw, output, d.Regions([]))
            self.assertFalse(output.exists())

    def test_cli_extraction_resume_and_failure_accounting(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = pathlib.Path(tmp)
            raw = tmp / 'source.tsv'
            raw.write_text('chr pos ea oa beta se\n1 100 A G .2 .1\n')
            source = dict(trait_id='good', url=raw.as_uri(), build='GRCh38',
                          columns=dict(chromosome='chr', position='pos', effect_allele='ea',
                                       other_allele='oa', beta='beta', se='se'))
            bad = dict(source, trait_id='bad', columns=dict(beta='missing'))
            registry = tmp / 'sources.json'
            registry.write_text(json.dumps(dict(sources=[source, bad])))
            (tmp / 'regions.csv').write_text('chromosome,start,end\n1,1,200\n')
            args = ['--job-dir', str(tmp), '--sources', str(registry), '--cache', str(tmp/'cache')]
            self.assertEqual(d.main(args), 1)
            statuses = json.loads((tmp/'gwas/download_status.json').read_text())
            self.assertEqual([s['status'] for s in statuses], ['ready', 'error'])
            self.assertTrue((tmp/'gwas/good.json').exists())
            self.assertFalse((tmp/'gwas/bad.json').exists())
            self.assertEqual(d.main(args + ['--traits', 'good']), 0)
            statuses = json.loads((tmp/'gwas/download_status.json').read_text())
            self.assertEqual(statuses[0]['status'], 'cached')


if __name__ == '__main__':
    unittest.main()
