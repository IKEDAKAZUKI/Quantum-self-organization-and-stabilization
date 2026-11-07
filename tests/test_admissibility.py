"""Finite-matrix checks of the control moment conditions used in Section 3."""
import unittest
import numpy as np


class ControlMomentTests(unittest.TestCase):
    def test_number_projection_preserves_both_diagonal_observables(self):
        nmax=3
        occupations=[(i,j) for i in range(nmax) for j in range(nmax)]
        u=np.array([(i+j-1)/2-1 for i,j in occupations])
        E=np.diag(u*u)
        P=np.diag([int(i==j) for i,j in occupations])
        projectors=[np.diag(np.arange(len(u))==i).astype(float) for i in range(len(u))]
        for O in (E,P):
            np.testing.assert_array_equal(sum(K.T@O@K for K in projectors),O)

    def test_number_dephasing_has_zero_moment_contribution(self):
        n=np.diag([0.,1.,2.,3.])
        E=np.diag(np.array([-.5,0,.5,1.])**2)
        P=np.diag([1.,0.,0.,1.])
        for O in (E,P):
            dual=n@O@n-.5*(n@n@O+O@n@n)
            np.testing.assert_allclose(dual,0,atol=0)

    def test_preserving_mean_contrast_does_not_control_its_second_moment(self):
        # Symmetric replacement of the central level by +/-1 preserves the
        # first moment for every input but increases the second at the centre.
        u=np.diag([-1.,0.,1.]);E=u@u
        K=[]
        for target,source,amplitude in ((0,0,1),(2,2,1),(0,1,2**-.5),(2,1,2**-.5)):
            A=np.zeros((3,3));A[target,source]=amplitude;K.append(A)
        dual=lambda O:sum(A.T@O@A for A in K)
        np.testing.assert_allclose(dual(np.eye(3)),np.eye(3))
        np.testing.assert_allclose(dual(u),u)
        self.assertGreater((dual(E)-E)[1,1],.99)


if __name__=='__main__':
    unittest.main()
