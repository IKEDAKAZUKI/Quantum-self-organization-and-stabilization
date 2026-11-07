"""Distinguish entangled basis states from their unconditional mixtures."""
import unittest
import numpy as np


def partial_transpose_b(rho):
    return rho.reshape(2, 2, 2, 2).transpose(0, 3, 2, 1).reshape(4, 4)


def negativity(rho):
    return float(-np.minimum(np.linalg.eigvalsh(partial_transpose_b(rho)), 0).sum())


class MixtureWitnessTests(unittest.TestCase):
    def test_equal_phase_mixture_is_separable(self):
        left = np.array([0., 1., 0., 0.])
        right = np.array([0., 0., 1., 0.])
        plus = (left + right) / np.sqrt(2.)
        minus = (left - right) / np.sqrt(2.)
        rho = (np.outer(plus, plus) + np.outer(minus, minus)) / 2.
        expected = (np.outer(left, left) + np.outer(right, right)) / 2.
        np.testing.assert_allclose(rho, expected, atol=2e-16)
        self.assertAlmostEqual(negativity(rho), 0., places=15)
        self.assertAlmostEqual(negativity(np.outer(plus, plus)), 0.5, places=15)

    def test_imbalance_is_certified_by_local_moments(self):
        left = np.array([0., 1., 0., 0.])
        right = np.array([0., 0., 1., 0.])
        plus = (left + right) / np.sqrt(2.)
        minus = (left - right) / np.sqrt(2.)
        C = np.outer(left, right) + np.outer(right, left)
        for probability in (0.5, 0.6, 0.9, 1.0):
            rho = probability*np.outer(plus, plus)+(1-probability)*np.outer(minus, minus)
            e = float(np.trace(C @ rho))
            self.assertAlmostEqual(e, 2*probability-1, places=14)
            self.assertGreaterEqual(negativity(rho)+1e-15, max(e, 0)/4)


if __name__ == '__main__':
    unittest.main()
