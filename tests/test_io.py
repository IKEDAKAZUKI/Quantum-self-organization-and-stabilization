import json
import tempfile
import unittest
from pathlib import Path
import numpy as np
from qtp.io import binary_to_npz,load_trajectory,ROOT
from qtp.observables import mean_coherence,occupation_contrast
from qtp.geometry import lifted_components


class DataTests(unittest.TestCase):
    def test_mode_one_is_not_coherent(self):
        s=np.ones((2,5,6),dtype=np.int8);s[:,0,0]=0
        np.testing.assert_array_equal(mean_coherence(s,1),np.zeros(2))
        np.testing.assert_allclose(mean_coherence(s,0),[29/30]*2)

    def test_contrast_decoder(self):
        n=np.array([[0,1],[-1,2]],dtype=np.int32);s=np.array([[1,0],[-1,0]])
        np.testing.assert_array_equal(occupation_contrast(n,s,.1),.1*np.array([[0,.5],[-1,1.5]]))

    def test_binary_round_trip(self):
        with tempfile.TemporaryDirectory() as td:
            td=Path(td);p,z=load_trajectory('confirm_stripe_00');m=p['nx']*p['ny']
            a=np.empty(2,dtype=np.dtype([('t','<f8'),('n','<i4',(m,)),('s','i1',(m,))]))
            a['t']=z.t[:2];a['n']=z.n[:2].reshape(2,m);a['s']=z.s[:2].reshape(2,m)
            with (td/'a.bin').open('wb') as f:np.array([p['nx'],p['ny'],m],dtype='<i4').tofile(f);a.tofile(f)
            binary_to_npz(td/'a.bin',td/'a.npz',p['eps'])
            with np.load(td/'a.npz') as back:
                for name in ['t','n','s']:np.testing.assert_array_equal(back[name],getattr(z,name)[:2])
            with self.assertRaises(FileExistsError):binary_to_npz(td/'a.bin',td/'a.npz',p['eps'])
            (td/'broken.bin').write_bytes((td/'a.bin').read_bytes()[:-1])
            with self.assertRaises(ValueError):binary_to_npz(td/'broken.bin',td/'b.npz',p['eps'])

    def test_configuration_paths_are_relative_and_present(self):
        files=list(ROOT.glob('data/*/*.cfg'));self.assertEqual(len(files),228)
        for cfg in files:
            initial=next(line.split(maxsplit=1)[1] for line in cfg.read_text().splitlines() if line.startswith('init '))
            self.assertFalse(Path(initial).is_absolute());self.assertTrue((cfg.parent/initial).is_file())

    def test_periodic_winding(self):
        m=np.zeros((20,24),bool);m[:,3:5]=True
        _,a,x,y=lifted_components(m);self.assertTrue(y[1]);self.assertFalse(x[1])
        m=np.zeros((20,24),bool);m[np.ix_([19,0,1],[23,0,1])]=True
        _,a,x,y=lifted_components(m);self.assertEqual(len(a),2);self.assertFalse(x[1] or y[1])

    def test_invalid_mode_and_labels(self):
        with self.assertRaises(ValueError):mean_coherence(np.ones((2,2)),2)
        with self.assertRaises(ValueError):occupation_contrast(np.ones((2,2)),np.full((2,2),9),.1)

if __name__=='__main__':unittest.main()
