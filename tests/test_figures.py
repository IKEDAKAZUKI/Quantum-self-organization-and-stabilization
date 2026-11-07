"""Check canvas margins of the shared normalized-axis and colorbar strip."""
from pathlib import Path
import tempfile
import unittest
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from qtp.figures import style, colorbar_panel, heatmap_panel, save_panel, PANEL_LAYOUTS

class FigureLayoutTests(unittest.TestCase):
    def setUp(self):
        style()
        self.tmp=tempfile.TemporaryDirectory(prefix='qtp figure test ')
        self.root=Path(self.tmp.name)
    def tearDown(self):
        plt.close('all')
        self.tmp.cleanup()
    def test_normalized_axis_label_and_colorbar_inside_canvas(self):
        colorbar_panel(self.root/'colorbar.pdf')
        self.assertLessEqual(PANEL_LAYOUTS[-1]['canvas_overflow_points'],.1)
    def test_heatmap_title_inside_canvas(self):
        heatmap_panel(np.zeros((56,48)),self.root/'map.pdf','(a) $t=600$')
        self.assertLessEqual(PANEL_LAYOUTS[-1]['canvas_overflow_points'],.1)
    def test_clipping_is_reported_before_saving(self):
        f=plt.figure(figsize=(2,1));f.text(.5,1.0,'clipped',va='bottom')
        with self.assertRaisesRegex(ValueError,'clips visible'):
            save_panel(f,self.root/'bad.pdf')
        self.assertFalse((self.root/'bad.pdf').exists())
    def test_maintenance_figure_contains_four_panels_and_two_images(self):
        import json
        import fitz
        import pandas as pd
        from qtp.figures import maintenance
        root = Path(__file__).resolve().parents[1]
        summary = pd.read_csv(root / 'data/reference/geometry_confirmation_summary.csv')
        maintenance(self.root, self.root, summary)
        with fitz.open(self.root / 'figure3_maintenance.pdf') as doc:
            text = doc[0].get_text()
        for title in ('(a) Formed configuration', '(b) Quantum evolution',
                      '(c) Mean pattern intensity', '(d) Final intensity'):
            self.assertIn(title, text)
        self.assertNotIn('(e)', text)
        self.assertNotIn('Number-projected evolution', text)
        images = json.loads((self.root / 'figure3_images.json').read_text())
        self.assertEqual([r['id'] for r in images],
                         ['confirm_labyrinth_01', 'maint_labyrinth_01_0'])
        endpoints = pd.read_csv(self.root / 'figure3_endpoint_values.csv')
        self.assertEqual(len(endpoints), 6)
        self.assertEqual(set(endpoints['mode']), {0, 1})
        for row in endpoints.itertuples():
            ref = summary[(summary.group == 'maintenance') &
                          (summary.family == row.family) &
                          (summary['mode'] == row.mode) & (summary.t == 120)].iloc[0]
            self.assertAlmostEqual(row.ratio, ref.mean_G / ref.mean_B, places=14)
            self.assertAlmostEqual(row.ratio_SE, ref.SE_G / ref.mean_B, places=14)

if __name__=='__main__':unittest.main()

