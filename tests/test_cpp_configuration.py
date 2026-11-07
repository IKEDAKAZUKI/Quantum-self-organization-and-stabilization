"""Exercise the C++ configuration contract directly and through Python."""
from pathlib import Path
import json
import subprocess
import tempfile
import unittest
import numpy as np

from qtp.io import ROOT, parse_config, simulation_parameters
from qtp.simulation import compile_engine, simulate


class CppConfigurationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory(prefix="qtp cpp config ")
        cls.root = Path(cls.temp.name)
        cls.engine = compile_engine(cls.root / "build")
        cls.serial = 0

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def make_case(self, changes=None, extra=""):
        type(self).serial += 1
        directory = self.root / str(type(self).serial)
        directory.mkdir()
        initial = np.c_[np.zeros(25, dtype=int), np.full(25, 9, dtype=int)]
        np.savetxt(directory / "initial state.txt", initial, fmt="%d")
        p = {"nx": 5, "ny": 5, "eps": 0.1, "U": 1, "seed": 7,
             "disseed": 8, "T": 0.02, "save": 0.01,
             "init": "initial state.txt", "phase_draw_order": 0}
        p.update(changes or {})
        cfg = directory / "case.cfg"
        cfg.write_text("\n".join(f"{k} {v}" for k, v in p.items()) + "\n" + extra)
        return cfg, directory / "out"

    def run_case(self, cfg, out):
        return subprocess.run([str(self.engine), str(cfg), str(out)], cwd="/",
                              text=True, capture_output=True, timeout=20)

    def reject(self, cfg, out, expected):
        result = self.run_case(cfg, out)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(expected, result.stderr)
        for suffix in (".bin", ".csv", "_engine.json"):
            self.assertFalse(Path(str(out) + suffix).exists())

    def test_unknown_key_rejected_before_output(self):
        cfg, out = self.make_case(extra="etta 999\n")
        self.reject(cfg, out, "unknown configuration key: etta")
        with self.assertRaisesRegex(ValueError, "unknown"):
            parse_config(cfg)

    def test_duplicate_key_rejected_before_output(self):
        cfg, out = self.make_case(extra="seed 9\n")
        self.reject(cfg, out, "duplicate configuration key: seed")
        with self.assertRaisesRegex(ValueError, "duplicate"):
            parse_config(cfg)

    def test_each_integer_field_rejects_fraction(self):
        for key in ("nx", "ny", "mode", "phaseinit", "phase_draw_order", "seed", "disseed"):
            with self.subTest(key=key):
                cfg, out = self.make_case({key: "0.5"})
                self.reject(cfg, out, "invalid integer parameter: " + key)
                with self.assertRaises(ValueError):
                    parse_config(cfg)

    def test_seed_range_and_trailing_text_rejected(self):
        for value in ("-1", "18446744073709551616", "7junk", "7 8", "NaN", "Inf", "1e-3"):
            with self.subTest(value=value):
                cfg, out = self.make_case({"seed": value})
                self.reject(cfg, out, "invalid integer parameter: seed")
                with self.assertRaises(ValueError):
                    parse_config(cfg)

    def test_numeric_field_rejects_nonfinite_and_suffix(self):
        for value in ("nan", "inf", "-inf", "1.0junk", "1 2"):
            with self.subTest(value=value):
                cfg, out = self.make_case({"eta": value})
                self.reject(cfg, out, "invalid numeric parameter: eta")
                with self.assertRaises(ValueError):
                    parse_config(cfg)

    def test_integral_decimal_and_scientific_values_match_exact_seeds(self):
        cfg, out = self.make_case({"nx": "5.0", "phase_draw_order": "1e0",
                                  "seed": "18446744073709551615.0",
                                  "disseed": "8e0"})
        result = self.run_case(cfg, out)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(parse_config(cfg)["seed"], 2**64 - 1)
        cfg2, out2 = self.make_case({"phase_draw_order": 1, "seed": 2**64 - 1})
        result2 = self.run_case(cfg2, out2)
        self.assertEqual(result2.returncode, 0, result2.stderr)
        self.assertEqual(out.with_suffix(".bin").read_bytes(), out2.with_suffix(".bin").read_bytes())
        self.assertEqual(out.with_suffix(".csv").read_bytes(), out2.with_suffix(".csv").read_bytes())

    def test_allowed_metadata_does_not_change_trajectory(self):
        metadata = {"family": "stripe", "group": "test", "id": "case", "parent": "parent",
                    "block": 0, "rep": 0, "stream": 0, "tf": 600, "E0": 0,
                    "initseed": 9, "sigma_init": 0.02, "engine_sha256": "abc",
                    "initial_sha256": "def", "configuration_source_sha256": "123",
                    "configuration_source_format": "cfg", "source_id": "original"}
        cfg, out = self.make_case(metadata)
        result = self.run_case(cfg, out)
        self.assertEqual(result.returncode, 0, result.stderr)
        cfg2, out2 = self.make_case()
        self.assertEqual(self.run_case(cfg2, out2).returncode, 0)
        self.assertEqual(out.with_suffix(".bin").read_bytes(), out2.with_suffix(".bin").read_bytes())

    def test_schedule_rejected_before_output(self):
        for change in ({"t1": -1}, {"t1": 2, "t2": 1}, {"eta_mid": -0.1}):
            with self.subTest(change=change):
                cfg, out = self.make_case(change)
                self.reject(cfg, out, "invalid dephasing schedule")

    def test_extra_initial_entries_rejected(self):
        cfg, out = self.make_case()
        initial = cfg.parent / "initial state.txt"
        initial.write_text(initial.read_text() + "0 1\n")
        self.reject(cfg, out, "extra entries in initial configuration")

    def test_python_and_cpp_execute_identical_valid_configurations(self):
        for phase_order in (0, 1):
            with self.subTest(phase_order=phase_order):
                cfg, out = self.make_case({"phase_draw_order": phase_order})
                result = self.run_case(cfg, out)
                self.assertEqual(result.returncode, 0, result.stderr)
                wrapped = cfg.parent / "python"
                simulate(cfg, wrapped, engine=self.engine)
                self.assertEqual(out.with_suffix(".bin").read_bytes(), (wrapped / "case.bin").read_bytes())
                self.assertEqual(out.with_suffix(".csv").read_bytes(), (wrapped / "case.csv").read_bytes())


if __name__ == "__main__":
    unittest.main()
