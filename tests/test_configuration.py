"""Regression tests for explicit CFG/JSON selection and actual engine execution."""
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
import numpy as np

from qtp.io import ROOT, load_parameters, parse_config, simulation_parameters, load_trajectory
from qtp.simulation import compile_engine, simulate


class ConfigurationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='qtp configuration test ')
        self.root = Path(self.tmp.name)
        self.base = json.loads((ROOT/'data/primary/confirm_stripe_00_params.json').read_text())
        self.cfg = self.root/'case.cfg'
        self.json = self.root/'case_params.json'

    def tearDown(self):
        self.tmp.cleanup()

    def test_explicit_cfg_edit_overrides_conflicting_sibling_json(self):
        self.cfg.write_text('eta 999\ninit chosen initial.txt\n')
        self.json.write_text(json.dumps(self.base))
        prefix, p = load_parameters(self.cfg)
        self.assertEqual(p['eta'], 999)
        self.assertEqual(p['init'], 'chosen initial.txt')
        self.assertNotIn('D2', p)  # No silent merge of sibling JSON.
        self.assertEqual(prefix, self.root/'case')

    def test_explicit_json_remains_authoritative(self):
        self.cfg.write_text('eta 999\n')
        self.json.write_text(json.dumps(self.base))
        _, p = load_parameters(self.json)
        self.assertEqual(p['eta'], 1.0)
        _, by_prefix = load_parameters(self.root/'case')
        self.assertEqual(by_prefix['eta'], 1.0)

    def test_standalone_cfg_with_no_sibling_json(self):
        self.cfg.write_text('# example\n\nnx 5\nny 6\neta 2e-4\nseed 7\ndisseed 8\ninit custom state.txt\n')
        np.savetxt(self.root/'custom state.txt', np.zeros((30,2),int),fmt='%d')
        prefix,p,initial=simulation_parameters(self.cfg)
        self.assertEqual(p['nx'],5)
        self.assertEqual(p['eta'],2e-4)
        self.assertEqual(p['D4'],1.0)  # C++ default, not a sidecar value.
        self.assertEqual(initial,self.root/'custom state.txt')

    def test_missing_explicit_cfg_never_falls_back_to_json(self):
        self.json.write_text(json.dumps(self.base))
        with self.assertRaises(FileNotFoundError):
            load_parameters(self.cfg)

    def test_duplicate_cfg_keys_rejected(self):
        self.cfg.write_text('eta 1\neta 2\n')
        with self.assertRaisesRegex(ValueError,'duplicate'):
            load_parameters(self.cfg)

    def test_unknown_cfg_key_rejected(self):
        self.cfg.write_text('etta 1\n')
        with self.assertRaisesRegex(ValueError,'unknown'):
            parse_config(self.cfg)

    def test_invalid_cfg_numbers_rejected(self):
        for line in ('eta nan','eta inf','nx 5.5','seed -1','seed 18446744073709551616','eta 2 junk'):
            with self.subTest(line=line):
                self.cfg.write_text(line+'\n')
                with self.assertRaises(ValueError):
                    parse_config(self.cfg)

    def test_uint64_seed_is_not_rounded(self):
        self.cfg.write_text('seed 18446744073709551615\n')
        self.assertEqual(parse_config(self.cfg)['seed'],2**64-1)

    def test_duplicate_json_keys_rejected(self):
        self.json.write_text('{"eta":1,"eta":2}')
        with self.assertRaisesRegex(ValueError,'duplicate'):
            load_parameters(self.json)

    def test_explicit_init_path_is_resolved_beside_source(self):
        (self.root/'inputs with spaces').mkdir()
        initial=self.root/'inputs with spaces'/'actual.txt'
        np.savetxt(initial,np.ones((25,2),int),fmt='%d')
        p=self.base|dict(nx=5,ny=5,init='inputs with spaces/actual.txt')
        self.json.write_text(json.dumps(p))
        _,_,found=simulation_parameters(self.json)
        self.assertEqual(found,initial)

    def test_edited_cfg_cannot_relabel_an_old_recorded_npz(self):
        with self.assertRaisesRegex(ValueError,'recorded trajectory'):
            load_trajectory(self.cfg)

    def test_all_packaged_configs_agree_with_recorded_parameters(self):
        paths=list((ROOT/'data').glob('*/*.cfg'))
        self.assertEqual(len(paths),228)
        for path in paths:
            cfg=parse_config(path)
            p=json.loads(path.with_name(path.stem+'_params.json').read_text())
            for k in cfg.keys() & p.keys():
                self.assertEqual(cfg[k],p[k],(path.name,k))


class ConfigurationExecutionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp=tempfile.TemporaryDirectory(prefix='qtp config execution ')
        cls.root=Path(cls.tmp.name)
        cls.engine=compile_engine(cls.root/'build')

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def make_case(self,folder):
        folder.mkdir()
        initial=np.c_[np.zeros(25,int),np.ones(25,int)]
        initial[0]=[1,0]
        np.savetxt(folder/'chosen state.txt',initial,fmt='%d')
        cfg=folder/'case.cfg'
        cfg.write_text('nx 5\nny 5\neps 0.1\nU 1\nT 0.03\nsave 0.01\n'
                       'eta 999\nseed 7\ndisseed 8\ninit chosen state.txt\n')
        stale=json.loads((ROOT/'data/primary/confirm_stripe_00_params.json').read_text())
        (folder/'case_params.json').write_text(json.dumps(stale))
        return cfg

    def test_edited_cfg_changes_executed_settings_and_matches_direct_cpp(self):
        cfg=self.make_case(self.root/'input with spaces')
        out=self.root/'python output'
        simulate(cfg,out,engine=self.engine)
        p=json.loads((out/'case_params.json').read_text())
        self.assertEqual(p['eta'],999)
        self.assertEqual(p['nx'],5)
        self.assertEqual(p['configuration_source_format'],'cfg')
        self.assertAlmostEqual(p['E0'],.05**2/25)
        direct=self.root/'direct'/'case'
        result=subprocess.run([str(self.engine),str(cfg),str(direct)],cwd='/',capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertEqual((out/'case.bin').read_bytes(),direct.with_suffix('.bin').read_bytes())
        self.assertEqual((out/'case.csv').read_bytes(),direct.with_suffix('.csv').read_bytes())
        with np.load(out/'case.npz') as z:
            self.assertEqual(z['n'].shape,(4,5,5))
            self.assertAlmostEqual(z['t'][-1],.03)

    def test_initial_contrast_uses_phase_preparation_override(self):
        cfg=self.make_case(self.root/'phase override')
        cfg.write_text(cfg.read_text()+'phaseinit 1\n')
        out=self.root/'phase output'
        simulate(cfg,out,until=0,engine=self.engine)
        p=json.loads((out/'case_params.json').read_text())
        self.assertAlmostEqual(p['E0'],.1**2/25)
        with np.load(out/'case.npz') as z:
            self.assertTrue((z['s'][0]==1).all())

    def test_cli_accepts_cfg_from_an_unrelated_working_directory(self):
        cfg=self.make_case(self.root/'cli input')
        out=self.root/'cli output'
        result=subprocess.run(['python',str(ROOT/'run_simulation.py'),str(cfg),'--until','0',
                               '--output',str(out)],cwd='/',capture_output=True,text=True,timeout=60)
        self.assertEqual(result.returncode,0,result.stderr)
        p=json.loads((out/'case_params.json').read_text())
        self.assertEqual(p['eta'],999)
        self.assertEqual(p['T'],0)
        self.assertEqual(p['nx'],5)


if __name__=='__main__':
    unittest.main()
