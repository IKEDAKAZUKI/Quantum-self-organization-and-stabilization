"""Numerical regression tests for equal-rate limits and the reported parameters."""
import json
from decimal import Decimal, localcontext
from pathlib import Path
import unittest
import numpy as np
from scipy.integrate import solve_ivp
from qtp.bounds import comparison_constants,decay_convolution,formation_bound,maintenance_bound,source_moments
from qtp.io import ROOT


class BoundTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.p=json.loads((ROOT/'data/primary/confirm_stripe_00_params.json').read_text())

    def test_equal_rate_exact(self):
        t=np.array([0,1e-10,.1,1,10,1000])
        np.testing.assert_array_equal(decay_convolution(t,.2,.2),t*np.exp(-.2*t))

    def test_near_rates_high_precision(self):
        worst=0.
        for a,k in [(.2,.2),(.2,.2+1e-15),(.2+1e-12,.2),(.2,2.0),(2.0,.2),(0,0),(0,.2)]:
            for t in [0.,1e-10,.1,1.,10.,1000.]:
                with localcontext() as ctx:
                    ctx.prec=75
                    aa,kk,tt=map(Decimal.from_float,(a,k,t))
                    exact=tt*(-aa*tt).exp() if aa==kk else ((-kk*tt).exp()-(-aa*tt).exp())/(aa-kk)
                value=float(decay_convolution(t,a,k));reference=float(exact)
                self.assertTrue(np.isfinite(value))
                self.assertAlmostEqual(value,reference,delta=abs(reference)*8e-14+1e-300)
                if reference:worst=max(worst,abs(value-reference)/abs(reference))
        self.assertLess(worst,8e-14)

    def test_no_overflow_at_long_times(self):
        with np.errstate(all='raise'):
            values=decay_convolution([0.,1.,1e6],.00001,1000.)
        self.assertTrue(np.all(np.isfinite(values)))

    def test_maintenance_equal_and_near_rates(self):
        p=dict(self.p);a=comparison_constants(p)['alpha']
        for delta in [0.,1e-15,-1e-15,1e-9,-1e-9]:
            p['repair']=a-p['leak']+delta
            q=comparison_constants(p);t=np.array([0,.1,1,10,100])
            got=maintenance_bound(t,p,.13)
            self.assertTrue(np.all(np.isfinite(got)))
            conv=decay_convolution(t,q['alpha'],q['k'])
            expected=q['B_inf']+(q['Ue']**2-q['B_inf'])*np.exp(-q['alpha']*t)+q['c1']*(.13-q['ell_inf'])*conv
            np.testing.assert_array_equal(got,expected)
            self.assertAlmostEqual(float(got[0]),q['Ue']**2)

    def test_bounds_against_ode(self):
        p=dict(self.p);alpha=comparison_constants(p)['alpha']
        for delta in [0.,1e-12,.1]:
            p['repair']=alpha-p['leak']+delta;q=comparison_constants(p)
            for ell0,initial,fn in [(0.,p['E0'],lambda t:formation_bound(t,p)),
                                    (.1,q['Ue']**2,lambda t:maintenance_bound(t,p,.1))]:
                def ode(t,y):return [-q['alpha']*y[0]+q['c0']+q['c1']*y[1],p['leak']-q['k']*y[1]]
                times=np.linspace(0,20,21);sol=solve_ivp(ode,(0,20),[initial,ell0],t_eval=times,rtol=2e-12,atol=1e-13)
                np.testing.assert_allclose(fn(times),sol.y[0],rtol=2e-10,atol=2e-11)

    def test_scalar_and_vector_consistency(self):
        t=np.array([0.,.3,10.]);b=maintenance_bound(t,self.p,.01)
        for j,x in enumerate(t):self.assertEqual(float(maintenance_bound(x,self.p,.01)),b[j])

    def test_source_initial_conditions_and_ode(self):
        p=dict(self.p,r=.4,eta=.1,leak=.02,repair=.6)
        # a=k here, requiring the same convolution limit in the source solution.
        t=np.linspace(0,10,30);e,l=source_moments(t,p,e0=.3,ell0=.1)
        sol=solve_ivp(lambda t,y:[p['r']*(1-y[1])-(p['r']+2*p['eta']+p['leak'])*y[0],
                                    p['leak']-(p['leak']+p['repair'])*y[1]],
                      (0,10),[.3,.1],t_eval=t,rtol=2e-12,atol=1e-13)
        np.testing.assert_allclose(e,sol.y[0],rtol=1e-10,atol=2e-12)
        np.testing.assert_allclose(l,sol.y[1],rtol=1e-10,atol=2e-12)

    def test_invalid_inputs(self):
        for f in [lambda:formation_bound(-1,self.p),lambda:maintenance_bound(1,self.p,-.1),
                  lambda:maintenance_bound(1,self.p,1.1),lambda:decay_convolution(1,np.nan,1),
                  lambda:comparison_constants(self.p|{'D4':0}),lambda:comparison_constants(self.p|{'delta':.0001})]:
            with self.assertRaises(ValueError):f()

if __name__=='__main__':unittest.main()
