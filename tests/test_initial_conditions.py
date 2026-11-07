"""Check the second moments and leakage used in formation bounds."""
from pathlib import Path
import json
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


class InitialConditionTests(unittest.TestCase):
    def test_rounded_initial_second_moments_match_recorded_parameters(self):
        checked = 0
        for group in ("primary", "preparations"):
            for path in sorted((ROOT / "data" / group).glob("*_params.json")):
                p = json.loads(path.read_text())
                if p["group"] not in ("confirmation", "fresh_formation"):
                    continue
                with np.load(path.with_name(p["id"] + ".npz"), allow_pickle=False) as z:
                    u = p["eps"] * (z["n"][0] - 0.5 * (z["s"][0] == 0))
                    self.assertFalse((z["s"][0] == 0).any())
                    self.assertAlmostEqual(float(np.mean(u * u)), p["E0"], places=15)
                checked += 1
        self.assertEqual(checked, 72)


if __name__ == "__main__":
    unittest.main()
