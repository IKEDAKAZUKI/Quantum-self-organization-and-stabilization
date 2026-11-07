"""Render all four article figures from trajectories, analysis tables and equations.

Each Matplotlib canvas contains one scientific panel. PyMuPDF assembles the
panels rendered from recorded configurations, summaries and analytical curves.
"""
from __future__ import annotations
from pathlib import Path
import hashlib
import json
import tempfile
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.lines import Line2D
import fitz
from .io import ROOT, load_trajectory, load_parameters
from .bounds import source_moments
from .safety import output_directory
FAMILIES = ('stripe', 'spot', 'labyrinth')
COLORS = ('#1f77b4', '#ff7f0e', '#2ca02c')
MARKERS = ('o', 's', '^')
PANEL_LAYOUTS = []


def style():
    """Times-compatible serif typography and the article's plotting conventions."""
    plt.rcParams.update({'font.family': 'serif', 'font.serif': ['Tinos', 'Times New Roman', 'Nimbus Roman'], 'mathtext.fontset': 'stix', 'font.size': 10, 'axes.titlesize': 11, 'axes.labelsize': 11, 'axes.linewidth': 0.7, 'xtick.direction': 'out', 'ytick.direction': 'out', 'xtick.top': False, 'ytick.right': False, 'legend.frameon': False, 'pdf.fonttype': 42, 'savefig.dpi': 180})


def save_panel(fig, path):
    """Save a panel after checking the visible artists against its canvas."""
    fig.canvas.draw()
    bbox = fig.get_tightbbox(fig.canvas.get_renderer())
    width, height = fig.get_size_inches()
    excess = max(0.0, -bbox.x0, -bbox.y0, bbox.x1 - width, bbox.y1 - height)
    PANEL_LAYOUTS.append(dict(panel=Path(path).name, width_inches=float(width), height_inches=float(height), visible_bounds_inches=list(bbox.extents), canvas_overflow_points=float(excess * 72)))
    if excess * 72 > 0.1:
        plt.close(fig)
        raise ValueError(f'Panel {path.name} clips visible text/artists by {excess * 72:.3f} pt')
    fig.savefig(path)
    plt.close(fig)
    return path


def assemble(parts, width, height, out):
    with fitz.open() as doc:
        page = doc.new_page(width=width, height=height)
        for path, rect in parts:
            with fitz.open(path) as part:
                page.show_pdf_page(fitz.Rect(*rect), part, 0, keep_proportion=False)
        doc.save(out, garbage=4, deflate=True)
        page.get_pixmap(matrix=fitz.Matrix(2.2, 2.2), alpha=False).save(out.with_suffix('.png'))


def line_panel(size=(7, 2.1), pos=(0.11, 0.25, 0.85, 0.56)):
    fig = plt.figure(figsize=size)
    return (fig, fig.add_axes(pos))


def heatmap_panel(u, path, title, show_y=True):
    fig = plt.figure(figsize=(1.62, 1.83))
    ax = fig.add_axes([0.19, 0.17, 0.77, 0.7])
    ax.imshow(u, origin='lower', extent=(0, 1, 0, 1), aspect='auto', interpolation='bilinear', cmap='RdBu_r', vmin=-1, vmax=1)
    ax.set_xticks([0, 0.5, 1], ['0', '0.5', '1'])
    ax.set_yticks([0, 0.5, 1], ['0', '0.5', '1'] if show_y else [])
    ax.set_title(title, pad=5)
    return save_panel(fig, path)


def colorbar_panel(path, width=7):
    fig = plt.figure(figsize=(width, 0.9))
    ax = fig.add_axes([0.3, 0.5, 0.43, 0.17])
    cb = fig.colorbar(ScalarMappable(norm=Normalize(-1, 1), cmap='RdBu_r'), cax=ax, orientation='horizontal', ticks=[-1, 0, 1])
    cb.set_label('Occupation contrast $u$', labelpad=2)
    fig.text(0.52, 0.97, 'Normalized x', ha='center', va='top', fontsize=10)
    return save_panel(fig, path)


def formation(out, tmp):
    parts = []
    records = []
    selected = ['confirm_stripe_00', 'confirm_spot_00', 'confirm_labyrinth_01']
    for row, (family, run) in enumerate(zip(FAMILIES, selected)):
        p, z = load_trajectory(run)
        u = z.u
        for col, t in enumerate([0, 20, 60, 600]):
            idx = np.flatnonzero(np.isclose(z.t, t, rtol=0, atol=1e-12))
            if len(idx) != 1:
                raise ValueError(f'Missing formation snapshot {run} at {t}')
            image = u[int(idx[0])]
            path = tmp / f'form-{row}-{col}.pdf'
            heatmap_panel(image, path, f'({chr(97 + row * 4 + col)}) $t={t}$', col == 0)
            parts.append((path, (32 + col * 116, 12 + row * 134, 32 + (col + 1) * 116, 144 + row * 134)))
            records.append(dict(id=run, t=t, index=int(idx[0]), minimum=float(image.min()), maximum=float(image.max())))
        fig = plt.figure(figsize=(0.44, 1.85))
        fig.text(0.5, 0.5, family.title(), ha='center', va='center', rotation=90)
        path = save_panel(fig, tmp / f'row-{row}.pdf')
        parts.append((path, (0, 12 + row * 134, 30, 144 + row * 134)))
    parts.append((colorbar_panel(tmp / 'form-colorbar.pdf'), (0, 417, 504, 481.8)))
    assemble(parts, 504, 490, out / 'figure1_formation.pdf')
    return records


def resource(out, tmp, summary):
    parts = []
    _, p = load_parameters('confirm_stripe_00')
    t = np.linspace(0, 0.3, 301)
    e, ell = source_moments(t, p)
    fig, ax = line_panel()
    ax.plot(t, e, label='Exact $e(t)$', color=COLORS[0], lw=1.5)
    ax.plot(t, 2 * ell, label='PPT bound $2\\ell(t)$', color='.35', ls='--', lw=1.2)
    ax.set(xlim=(0, 0.3), ylim=(-0.02, 1.0), xlabel='Time $t$', ylabel='Local coherence')
    ax.set_xticks([0, 0.1, 0.2, 0.3])
    ax.set_title('(a) Local intermode coherence', loc='left', pad=8)
    ax.legend(loc='lower right', fontsize=9)
    parts.append((save_panel(fig, tmp / 'resource-a.pdf'), (0, 0, 504, 151.2)))
    exported = []
    for panel in (1, 2):
        fig, ax = line_panel(pos=(0.11, 0.25, 0.85, 0.55))
        for fam, color in zip(FAMILIES, COLORS):
            q = summary[(summary.group == 'confirmation') & (summary.family == fam) & (summary['mode'] == 0)].sort_values('t')
            x = q.t.to_numpy()
            y = q.transport_witness.to_numpy() if panel == 1 else (q.mean_G / q.mean_B).to_numpy()
            se = q.SE_transport_witness.to_numpy() if panel == 1 else (q.SE_G / q.mean_B).to_numpy()
            ax.plot(x, y, '.-', ms=3.8, lw=1.2, color=color, label=fam.title())
            ax.fill_between(x, y - se, y + se, color=color, alpha=0.16, linewidth=0)
            exported.extend((dict(panel=panel + 1, family=fam, t=float(a), value=float(b), SE=float(c)) for a, b, c in zip(x, y, se)))
        ax.set_xlim(0, 600)
        ax.set_xticks(np.arange(0, 601, 100))
        ax.set_xlabel('Time $t$')
        if panel == 1:
            ax.set_ylabel('$\\mathcal{W}_{\\mathrm{tr}}$')
            ax.set_ylim(-0.02, 0.75)
            ax.set_title('(b) Transport above the separable bound', loc='left', pad=22)
            ax.legend(loc='lower left', bbox_to_anchor=(0, 1.01), ncol=3, fontsize=9, borderaxespad=0)
        else:
            ax.set_ylabel('$\\mathcal{G}_g/B_{\\mathrm{sep}}$')
            ax.set_ylim(-0.25, 8.2)
            ax.axhline(1, color='.35', ls=':', lw=1)
            ax.set_title('(c) Pattern intensity relative to the separable bound', loc='left', pad=8)
        parts.append((save_panel(fig, tmp / f'resource-{panel}.pdf'), (0, panel * 151.2, 504, (panel + 1) * 151.2)))
    assemble(parts, 504, 453.6, out / 'figure2_resource.pdf')
    pd.DataFrame(exported).to_csv(out / 'figure2_values.csv', index=False)


def maintenance(out, tmp, summary):
    parts = []
    image_records = []
    image_cases = [('confirm_labyrinth_01', '(a) Formed configuration'), ('maint_labyrinth_01_0', '(b) Quantum evolution')]
    for col, (run, title) in enumerate(image_cases):
        p, z = load_trajectory(run)
        fig = plt.figure(figsize=(3.5, 2.8))
        ax = fig.add_axes([0.18, 0.16, 0.78, 0.72])
        ax.imshow(z.u[-1], origin='lower', extent=(0, 1, 0, 1), aspect='auto', interpolation='bilinear', cmap='RdBu_r', vmin=-1, vmax=1)
        ax.set_xticks([0, 0.5, 1], ['0', '0.5', '1'])
        ax.set_yticks([0, 0.5, 1], ['0', '0.5', '1'] if col == 0 else [])
        if col == 0:
            ax.set_ylabel('Normalized y', fontsize=10)
        ax.set_title(title, pad=8, fontsize=10.5)
        parts.append((save_panel(fig, tmp / f'maint-{col}.pdf'), (col * 252, 0, (col + 1) * 252, 201.6)))
        image_records.append(dict(id=run, time=float(z.t[-1]), mode=int(p['mode']), minimum=float(z.u[-1].min()), maximum=float(z.u[-1].max()), contrast_array_sha256=hashlib.sha256(z.u[-1].tobytes()).hexdigest()))
    parts.append((colorbar_panel(tmp / 'maint-cbar.pdf'), (0, 201.6, 504, 266.4)))
    fig, ax = line_panel((4.2, 2.8), (0.16, 0.24, 0.8, 0.57))
    for fam, color in zip(FAMILIES, COLORS):
        for mode in [0, 1]:
            q = summary[(summary.group == 'maintenance') & (summary.family == fam) & (summary['mode'] == mode)].sort_values('t')
            ax.plot(q.t, q.mean_G, color=color, lw=1.3, ls='-' if mode == 0 else '--', label=fam.title() if mode == 0 else None)
            ax.fill_between(q.t, q.mean_G - q.SE_G, q.mean_G + q.SE_G, color=color, alpha=0.13, linewidth=0)
    ax.set(xlim=(0, 120), ylim=(-0.005, 0.27), xlabel='Elapsed time $\\tau$', ylabel='$\\mathcal{G}_g$')
    ax.set_xticks([0, 30, 60, 90, 120])
    ax.set_yticks([0, 0.1, 0.2])
    fig.text(0.16, 0.94, '(c) Mean pattern intensity', fontsize=10.5)
    ax.legend(loc='lower left', bbox_to_anchor=(0, 1.04), ncol=3, fontsize=8.5, handlelength=1.1, columnspacing=0.7, handletextpad=0.35, borderaxespad=0)
    parts.append((save_panel(fig, tmp / 'maint-c.pdf'), (0, 266.4, 302.4, 468)))
    fig, ax = line_panel((2.8, 2.8), (0.25, 0.24, 0.7, 0.57))
    endrows = []
    for j, (fam, color, marker) in enumerate(zip(FAMILIES, COLORS, MARKERS)):
        for mode, offset in [(0, -0.08), (1, 0.08)]:
            q = summary[(summary.group == 'maintenance') & (summary.family == fam) & (summary['mode'] == mode) & (summary.t == 120)].iloc[0]
            y = float(q.mean_G / q.mean_B)
            se = float(q.SE_G / q.mean_B)
            ax.errorbar(j + offset, y, yerr=se, fmt=marker, color=color, ms=5, mfc=color if mode == 0 else 'white', capsize=3, lw=1.1)
            endrows.append(dict(family=fam, mode=mode, ratio=y, ratio_SE=se))
    ax.axhline(1, color='.35', ls=':', lw=1)
    ax.set(xlim=(-0.5, 2.5), ylim=(-0.35, 8.2), ylabel='$\\mathcal{G}_g/\\langle B_f\\rangle$')
    ax.set_xticks([0, 1, 2], ['Stripe', 'Spot', 'Lab.'], fontsize=9)
    ax.set_yticks([0, 2, 4, 6, 8])
    fig.text(0.25, 0.94, '(d) Final intensity', fontsize=10.5)
    handles = [Line2D([], [], marker='o', ls='', mfc='.3', mec='.3', ms=4, label='Quantum'), Line2D([], [], marker='o', ls='', mfc='white', mec='.3', ms=4, label='Projected')]
    ax.legend(handles=handles, loc='lower left', bbox_to_anchor=(-0.06, 1.04), ncol=2, fontsize=8.5, handlelength=0.6, handletextpad=0.4, columnspacing=0.8, borderaxespad=0)
    parts.append((save_panel(fig, tmp / 'maint-d.pdf'), (302.4, 266.4, 504, 468)))
    assemble(parts, 504, 468, out / 'figure3_maintenance.pdf')
    pd.DataFrame(endrows).to_csv(out / 'figure3_endpoint_values.csv', index=False)
    (out / 'figure3_images.json').write_text(json.dumps(image_records, indent=2) + '\n')


def control(out, tmp, cycles):
    parts = []
    growth = []
    e = np.linspace(0, 1, 301)
    fig, ax = line_panel((7, 2.5), (0.11, 0.22, 0.85, 0.57))
    for fam, color in zip(FAMILIES, COLORS):
        _, p = load_parameters('confirm_' + fam + ('_01' if fam == 'labyrinth' else '_00'))
        sigma = -p['delta'] + (p['D2'] * e + p['Dcl']) ** 2 / (4 * p['D4'])
        ax.plot(e, sigma, color=color, lw=1.4, label=fam.title())
        growth.extend((dict(family=fam, coherence=float(x), sigma_max=float(y)) for x, y in zip(e, sigma)))
    ax.axhline(0, color='.5', ls=':', lw=0.8)
    _, p = load_parameters('confirm_stripe_00')
    ell = p['leak'] / (p['leak'] + p['repair'])
    for eta in [30, 1]:
        es = p['r'] * (1 - ell) / (p['r'] + 2 * eta + p['leak'])
        ax.axvline(es, color='.5', ls=':', lw=0.9)
        ax.text(es + 0.012, 0.43 if eta == 30 else -0.25, f'$\\eta={eta}$', fontsize=9)
    ax.set(xlim=(0, 1), ylim=(-0.33, 0.55), xlabel='Mean local coherence $e$', ylabel='Growth rate $\\sigma_{\\max}$')
    ax.set_title('(a) Maximum linear growth rate', loc='left', pad=23)
    ax.legend(loc='lower left', bbox_to_anchor=(0, 1.02), ncol=3, fontsize=9, borderaxespad=0)
    parts.append((save_panel(fig, tmp / 'control-a.pdf'), (0, 0, 504, 180)))
    fig, ax = line_panel((7, 2.5), (0.11, 0.22, 0.85, 0.61))
    ax.axvspan(250, 320, color='.9', zorder=0)
    for fam, color in zip(FAMILIES, COLORS):
        q = cycles[cycles.family == fam].sort_values('t')
        ax.plot(q.t, q['mean'], color=color, lw=1.3)
        ax.fill_between(q.t, q['mean'] - q.SE, q['mean'] + q.SE, color=color, alpha=0.16, linewidth=0)
    ax.text(285, 0.24, '$\\eta=30$', ha='center', fontsize=9)
    ax.set(xlim=(0, 650), ylim=(-0.003, 0.27), xlabel='Time $t$', ylabel='Band power $\\mathcal{P}_{\\mathcal{B}}$')
    ax.set_xticks([0, 100, 200, 300, 400, 500, 600])
    ax.set_title('(b) Suppression and recovery of band power', loc='left', pad=8)
    parts.append((save_panel(fig, tmp / 'control-b.pdf'), (0, 180, 504, 360)))
    assemble(parts, 504, 360, out / 'figure4_control.pdf')
    pd.DataFrame(growth).to_csv(out / 'figure4_linear_values.csv', index=False)


def overview(out, tmp):
    """Display the three representative formed configurations at t=600."""
    parts = []
    selected = ('confirm_stripe_00', 'confirm_spot_00', 'confirm_labyrinth_01')
    for column, (family, run) in enumerate(zip(FAMILIES, selected)):
        _, trajectory = load_trajectory(run)
        indices = np.flatnonzero(np.isclose(trajectory.t, 600, rtol=0, atol=1e-12))
        if len(indices) != 1:
            raise ValueError(f'Missing overview snapshot for {run} at t=600')
        path = heatmap_panel(trajectory.u[int(indices[0])],
                             tmp / f'overview-{family}.pdf', family.title())
        parts.append((path, (168*column, 0, 168*(column+1), 189.6)))
    parts.append((colorbar_panel(tmp / 'overview-colorbar.pdf'),
                  (0, 189.6, 504, 254.4)))
    assemble(parts, 504, 254.4, out / 'research_overview.pdf')


def render_all(analysis: Path, output: Path) -> dict:
    """Render article figures using the output of analyze_trajectories.py."""
    output = output_directory(output)
    analysis = Path(analysis).resolve()
    for name in ['primary_summary.csv', 'dephasing_summary.csv']:
        if not (analysis / name).is_file():
            raise FileNotFoundError(f'Run analyze_trajectories.py first: missing {analysis / name}')
    PANEL_LAYOUTS.clear()
    style()
    primary = pd.read_csv(analysis / 'primary_summary.csv')
    cycles = pd.read_csv(analysis / 'dephasing_summary.csv')
    with tempfile.TemporaryDirectory(prefix='qtp-panels-') as td:
        tmp = Path(td)
        overview(output, tmp)
        images = formation(output, tmp)
        resource(output, tmp, primary)
        maintenance(output, tmp, primary)
        control(output, tmp, cycles)
    input_paths = [analysis / 'primary_summary.csv', analysis / 'dephasing_summary.csv']
    for run in ['confirm_stripe_00', 'confirm_spot_00', 'confirm_labyrinth_01', 'maint_labyrinth_01_0']:
        prefix, _ = load_parameters(run)
        input_paths.extend([prefix.with_suffix('.npz'), prefix.with_name(run + '_params.json')])
    provenance = dict(status='passed', figures=4, overview='research_overview', formation_configurations=images, input_sha256={str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in input_paths}, heatmap=dict(cmap='RdBu_r', vmin=-1, vmax=1, interpolation='bilinear'), source_curve='Analytical expectation from source_moments', growth_curve='Continuous-wavevector homogeneous maximum, Eq. (16)', panel_layouts=PANEL_LAYOUTS.copy())
    (output / 'figure_provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    return provenance
