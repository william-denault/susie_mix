"""Regression: Windows builds must produce Linux-readable checksum manifests."""
import contextlib
import hashlib
import importlib.util
import io
import tempfile
import unittest
import zipfile
from pathlib import Path


class PackagingTest(unittest.TestCase):
    def test_built_manifest_has_literal_lf_and_names_match_archive(self):
        source=Path(__file__).resolve().parents[1]/'build_bundle.py'
        spec=importlib.util.spec_from_file_location('online_builder',source)
        builder=importlib.util.module_from_spec(spec);spec.loader.exec_module(builder)
        scratch=Path('tmp').resolve();scratch.mkdir(exist_ok=True)
        temp=tempfile.TemporaryDirectory(prefix='online_package_',dir=str(scratch))
        root=Path(temp.name).resolve()
        try:
            self.assertIn(scratch,root.parents)
            builder.ROOT=root;builder.HERE=root/'authoring';builder.BASE=root/'baseline'
            builder.DEST=root/'output/bundle'
            builder.HERE.mkdir();(root/'output').mkdir()
            model=builder.BASE/'script/scan_tissue_attempt';model.mkdir(parents=True)
            (model/'example.R').write_bytes(b'x <- 1\n')
            for name in ('renamed_gene_map.csv','excluded_candidates.csv'):
                (builder.BASE/name).write_bytes(b'gene\nEXAMPLE\n')
            (builder.HERE/'start.sh').write_bytes(b'#!/bin/bash\r\nexit 0\r\n')
            (builder.HERE/'a file.txt').write_bytes(b'fixture\n')
            with contextlib.redirect_stdout(io.StringIO()):builder.build()
            archive=next((root/'output').glob('*.zip'))
            with zipfile.ZipFile(archive) as z:
                manifest=z.read('bundle/SHA256SUMS')
                self.assertNotIn(b'\r',manifest)
                self.assertTrue(manifest.endswith(b'\n'))
                # Deliberately do NOT use splitlines()/read_text(): those hid
                # the CRLF defect by stripping CR from Linux filenames.
                for line in manifest.split(b'\n')[:-1]:
                    digest,name=line.split(b'  ',1)
                    payload=z.read('bundle/'+name.decode('utf-8'))
                    self.assertEqual(hashlib.sha256(payload).hexdigest().encode(),digest)
                self.assertNotIn(b'\r',z.read('bundle/start.sh'))
        finally:
            assert scratch in root.parents
            temp.cleanup()


if __name__=='__main__':unittest.main()
