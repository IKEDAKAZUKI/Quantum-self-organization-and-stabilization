"""Portable input and failure-path tests for the C++ generator."""
import json
import subprocess
import tempfile
import unittest
from pathlib import Path
import numpy as np
from qtp.io import ROOT,binary_to_npz
from qtp.simulation import compile_engine


class GeneratorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp=tempfile.TemporaryDirectory(prefix='qtp generator tests ')
        cls.root=Path(cls.temp.name)
        cls.engine=compile_engine(cls.root/'build')
        cls.params=json.loads((ROOT/'data/primary/confirm_stripe_00_params.json').read_text())

    @classmethod
    def tearDownClass(cls):cls.temp.cleanup()

    def config(self,tag,**changes):
        folder=self.root/tag;folder.mkdir(exist_ok=True)
        np.savetxt(folder/'initial state.txt',np.c_[np.zeros(25,dtype=int),np.ones(25,dtype=int)],fmt='%d')
        p=self.params|dict(nx=5,ny=5,T=.01,save=.01,init='initial state.txt')|changes
        cfg=folder/'case.cfg';cfg.write_text('\n'.join(f'{k} {v}' for k,v in p.items() if isinstance(v,(int,float,str)))+'\n')
        return cfg,folder/'out'

    def run_case(self,cfg,prefix):
        return subprocess.run([str(self.engine),str(cfg),str(prefix)],cwd='/',text=True,capture_output=True,timeout=30)

    def test_relative_initial_path_and_different_working_directory(self):
        cfg,out=self.config('relative path')
        result=self.run_case(cfg,out);self.assertEqual(result.returncode,0,result.stderr)
        binary_to_npz(out.with_suffix('.bin'),out.with_suffix('.npz'),self.params['eps'])
        with np.load(out.with_suffix('.npz')) as z:
            np.testing.assert_array_equal(z['t'],[0,.01]);self.assertEqual(z['n'].shape,(2,5,5))

    def test_missing_config_reports_path(self):
        cfg=self.root/'missing.cfg';result=self.run_case(cfg,self.root/'missing-output')
        self.assertNotEqual(result.returncode,0);self.assertIn(str(cfg),result.stderr)

    def test_missing_initial_file_reports_path(self):
        cfg,out=self.config('missing-init',init='absent.txt');result=self.run_case(cfg,out)
        self.assertNotEqual(result.returncode,0);self.assertIn('absent.txt',result.stderr)

    def test_invalid_mode_rejected(self):
        cfg,out=self.config('bad-mode',mode=2);self.assertNotEqual(self.run_case(cfg,out).returncode,0)

    def test_negative_rate_rejected(self):
        cfg,out=self.config('bad-rate',eta=-1);self.assertNotEqual(self.run_case(cfg,out).returncode,0)

    def test_existing_output_not_overwritten(self):
        cfg,out=self.config('existing');path=out.with_suffix('.bin');path.write_bytes(b'read-only-reference')
        result=self.run_case(cfg,out);self.assertNotEqual(result.returncode,0)
        self.assertEqual(path.read_bytes(),b'read-only-reference')

if __name__=='__main__':unittest.main()
