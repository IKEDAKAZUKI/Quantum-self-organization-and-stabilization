"""Shared configuration validation and reproducible run identifiers."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
from qtp.io import ROOT, load_parameters, simulation_parameters, validate_parameter_keys
from qtp.simulation import compile_engine, simulate


class JsonContractTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='qtp json contract ')
        self.root = Path(self.tmp.name)
        self.path = self.root/'edited.json'
        self.base = dict(nx=5, ny=5, eps=.1, U=1, T=.02, save=.01,
                         seed=7, disseed=8, init='state.txt')
        np.savetxt(self.root/'state.txt', np.c_[np.zeros(25, int), np.ones(25, int)], fmt='%d')
    def tearDown(self):
        self.tmp.cleanup()
    def write(self, values):
        self.path.write_text(json.dumps(values))
        return self.path
    def test_unknown_json_key_rejected_before_defaults(self):
        self.write(self.base | {'etta':999})
        with self.assertRaisesRegex(ValueError, r"unknown.*etta.*eta"):
            simulation_parameters(self.path)
    def test_typo_rejected_when_valid_key_is_also_present(self):
        with self.assertRaisesRegex(ValueError, 'etta'):
            load_parameters(self.write(self.base | {'eta':1,'etta':999}))
    def test_unknown_settings_not_silently_accepted_as_metadata(self):
        with self.assertRaisesRegex(ValueError, 'unrecognised_parameter'):
            load_parameters(self.write(self.base | {'unrecognised_parameter':1}))
    def test_rejected_simulation_does_not_compile_or_write_results(self):
        self.write(self.base | {'etta':999})
        output = self.root/'rejected'
        with patch('qtp.simulation.compile_engine') as build:
            with self.assertRaises(ValueError):
                simulate(self.path,output)
            build.assert_not_called()
        self.assertFalse(output.exists())
    def test_cli_typo_fails_with_input_path_and_suggestion(self):
        self.write(self.base | {'etta':999})
        result=subprocess.run([sys.executable,str(ROOT/'run_simulation.py'),str(self.path),
                               '--output',str(self.root/'cli')],cwd='/',text=True,capture_output=True)
        self.assertNotEqual(result.returncode,0)
        self.assertIn('etta',result.stderr)
        self.assertIn('eta',result.stderr)
        self.assertIn(str(self.path),result.stderr)
        self.assertFalse((self.root/'cli').exists())
    def test_cfg_json_use_the_same_allowed_key_set(self):
        for ext in ('cfg','json'):
            with self.subTest(ext=ext):
                p=self.root/f'case.{ext}'
                p.write_text('etta 999\n' if ext=='cfg' else '{"etta":999}')
                with self.assertRaisesRegex(ValueError,'unknown'):
                    load_parameters(p)
    def test_metadata_and_source_id_are_allowed(self):
        p=self.write(self.base | dict(id='input',source_id='ancestor',parent='parent',family='stripe',
                                    group='test',block=0,rep=1,stream=2))
        _,values=load_parameters(p)
        self.assertEqual(values['source_id'],'ancestor')
    def test_nonfinite_json_constants_rejected(self):
        for value in ('NaN','Infinity','-Infinity'):
            with self.subTest(value=value):
                self.path.write_text('{"eta":'+value+'}')
                with self.assertRaisesRegex(ValueError,'non-finite'):
                    load_parameters(self.path)
    def test_invalid_json_top_level_rejected(self):
        for value in ('[]','null','1'):
            self.path.write_text(value)
            with self.assertRaisesRegex(ValueError,'JSON object'):
                load_parameters(self.path)
    def test_invalid_init_type_rejected(self):
        self.write(self.base | {'init': []})
        with self.assertRaisesRegex(ValueError,'path string'):
            simulation_parameters(self.path)
    def test_valid_json_changes_eta_without_changing_default_others(self):
        _,p,_=simulation_parameters(self.write(self.base | {'eta':999}))
        self.assertEqual(p['eta'],999)
        self.assertEqual(p['D4'],1)
    def test_all_distributed_json_keys_remain_supported(self):
        paths=list((ROOT/'data').glob('*/*_params.json'))
        self.assertEqual(len(paths),228)
        for path in paths:
            load_parameters(path)


class RunIdentityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp=tempfile.TemporaryDirectory(prefix='qtp identity ')
        cls.root=Path(cls.tmp.name)
        cls.engine=compile_engine(cls.root/'build')
    @classmethod
    def tearDownClass(cls): cls.tmp.cleanup()
    def run_case(self,name,extension,old_id=None):
        folder=self.root/name;folder.mkdir()
        np.savetxt(folder/'init.txt',np.c_[np.zeros(25,int),np.ones(25,int)],fmt='%d')
        p=dict(nx=5,ny=5,eps=.1,U=1,T=.01,save=.01,eta=999,seed=7,disseed=8,init='init.txt')
        if old_id is not None:p['id']=old_id
        source=folder/(name+extension)
        source.write_text(json.dumps(p) if extension=='.json' else '\n'.join(f'{k} {v}' for k,v in p.items()))
        output=folder/'output';report=simulate(source,output,engine=self.engine)
        recorded=json.loads((output/(name+'_params.json')).read_text())
        self.assertEqual(report['id'],name)
        self.assertEqual(recorded['id'],name)
        self.assertTrue((output/(name+'.npz')).is_file())
        _,cfg=load_parameters(output/(name+'.cfg'))
        self.assertEqual(cfg['id'],name)
        _,json_pars=load_parameters(output/(name+'_params.json'))
        self.assertEqual(json_pars['id'],name)
        return recorded,source,output
    def test_cfg_without_id_gets_output_id(self):
        p,_,_=self.run_case('edited','.cfg')
        self.assertNotIn('source_id',p)
    def test_json_copied_old_id_is_preserved_as_source_id(self):
        p,source,output=self.run_case('standalone','.json','older_run')
        self.assertEqual(p['source_id'],'older_run')
        self.assertEqual(json.loads(source.read_text())['id'],'older_run')
    def test_matching_id_does_not_create_false_ancestor(self):
        p,_,_=self.run_case('matching','.json','matching')
        self.assertNotIn('source_id',p)
    def test_custom_json_matches_direct_cpp(self):
        p,source,output=self.run_case('direct_json','.json','ancestor')
        prefix=output/'direct_json';direct=source.parent/'direct_cpp'
        result=subprocess.run([str(self.engine),str(prefix.with_suffix('.cfg')),str(direct)],
                              cwd='/',capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
        for ext in ('.bin','.csv'):
            self.assertEqual(prefix.with_suffix(ext).read_bytes(),direct.with_suffix(ext).read_bytes())

if __name__=='__main__': unittest.main()
