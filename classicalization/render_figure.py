#!/usr/bin/env python3
"""Render the combined classicalization figure from saved simulation data.

Each scientific panel is drawn on its own Matplotlib canvas. The PDF panels are
assembled without rasterizing their text, axes, or curves.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import tempfile
import fitz
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from reference_time import completed_reference_time

ROOT = Path(__file__).resolve().parent
FAMILIES = ('stripe', 'spot', 'labyrinth')
COLORS = ('#1f77b4', '#ff7f0e', '#2ca02c')
WIDTH, HEIGHT = 7.55, 8.25


def style() -> None:
    plt.rcParams.update({
        'font.family': 'serif', 'font.serif': ['Tinos', 'Times New Roman', 'Nimbus Roman'],
        'mathtext.fontset': 'stix', 'font.size': 9.3, 'axes.labelsize': 9.7,
        'axes.titlesize': 10, 'xtick.labelsize': 8.6, 'ytick.labelsize': 8.6,
        'axes.linewidth': .7, 'xtick.direction': 'out', 'ytick.direction': 'out',
        'xtick.top': False, 'ytick.right': False, 'legend.frameon': False,
        'legend.fontsize': 8.5, 'pdf.fonttype': 42, 'ps.fonttype': 42,
        'savefig.dpi': 300,
    })


def subset(summary: pd.DataFrame, family: str, lo: float, hi: float) -> pd.DataFrame:
    q = summary[(summary.tag == 'main') & (summary.protocol == 'autonomous') &
                (summary.family == family) & summary.t.between(lo, hi)].sort_values('t')
    if q.empty:
        raise ValueError(f'Missing baseline data for {family}.')
    if not (q.power_count == 16).all():
        raise ValueError('The publication figure requires all 16 baseline realizations per family.')
    return q


def axes_labels(ax, ylabel: str) -> None:
    ax.set_xlabel('Time $t$', labelpad=2.5)
    ax.set_ylabel(ylabel, labelpad=3.5)
    ax.tick_params(length=3.1, width=.7, pad=2.1)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=ROOT)
    parser.add_argument('--output', type=Path, default=ROOT / 'figures')
    args = parser.parse_args()
    source, output = args.source.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    t_star, reference = completed_reference_time(source)
    summary = pd.read_csv(source / 'analysis/summary.csv')
    ends = pd.read_csv(source / 'analysis/endpoints.csv')
    if 'power_over_t_star_mean' not in summary:
        raise ValueError('Run analyze_study.py to compute the completion-referenced summaries.')
    stored_reference = json.loads((source / 'analysis/reference_time.json').read_text())
    if stored_reference['value'] != t_star:
        raise ValueError('The analysis reference time disagrees with the trajectory data.')
    style()
    selected = {}
    for family in FAMILIES:
        rows = ends[(ends.tag == 'main') & (ends.protocol == 'autonomous') &
                    (ends.family == family)].sort_values('power_final')
        selected[family] = rows.iloc[(len(rows)-1)//2].id
    times = (0, 40, 150, t_star, 1200)
    curve_rows, map_rows, artists = [], [], []
    parts = []

    def record(q, panel, **columns):
        for i, time in enumerate(q.t.to_numpy()):
            row = {'panel': panel, 'family': str(q.family.iloc[0]), 't': float(time)}
            row.update({key: float(np.asarray(value)[i]) for key, value in columns.items()})
            curve_rows.append(row)

    with tempfile.TemporaryDirectory(prefix='classicalization-panels-') as td:
        temporary = Path(td)
        def save(fig, name, rect):
            fig.canvas.draw()
            box = fig.get_tightbbox(fig.canvas.get_renderer())
            w, h = fig.get_size_inches()
            overflow = max(0., -box.x0, -box.y0, box.x1-w, box.y1-h)*72
            artists.append({'panel': name, 'canvas_overflow_points': float(overflow)})
            if overflow > .5:
                raise ValueError(f'{name} overflows its canvas by {overflow:.2f} points.')
            path = temporary / (name + '.pdf')
            fig.savefig(path, transparent=True)
            plt.close(fig)
            parts.append((path, rect))

        for index, panel in enumerate(('a', 'b', 'c')):
            fig = plt.figure(figsize=(2.48, 2.10))
            ax = fig.add_axes([.58/2.48, .35/2.10, 1.76/2.48, 1.33/2.10])
            if panel == 'a':
                for family, color in zip(FAMILIES, COLORS):
                    q = subset(summary, family, 0, 400)
                    y = q.activated_fraction_mean.to_numpy()
                    se = q.activated_fraction_sem.to_numpy()
                    ax.plot(q.t, y, color=color, lw=1.2, label=family.title())
                    ax.errorbar(q.t, y, yerr=se, fmt='none', ecolor=color,
                                elinewidth=.75, capsize=1.9, errorevery=16)
                    record(q, panel, mean=y, standard_error=se)
                ax.set(xlim=(0,400), ylim=(-.02,1.03), xticks=[0,100,200,300,400], yticks=[0,.5,1])
                axes_labels(ax, r'Population in $\vert1\rangle_{\mathrm{env}}$')
                ax.legend(loc='lower right', fontsize=8, handlelength=1.25,
                          borderaxespad=.35, labelspacing=.17)
            elif panel == 'b':
                q = subset(summary, 'labyrinth', 0, 400)
                ax.plot(q.t, q.ensemble_negativity_lower, color=COLORS[0], lw=1.25, label='Lower bound')
                ax.plot(q.t, q.negativity_upper_mean, color=COLORS[1], lw=1.25, ls='--', label='Upper bound')
                ax.set(xlim=(0,400), ylim=(-.012,.525), xticks=[0,100,200,300,400], yticks=[0,.25,.5])
                axes_labels(ax, 'Mean local negativity')
                ax.legend(loc='upper right', fontsize=8, handlelength=1.4,
                          borderaxespad=.35, labelspacing=.2)
                record(q, panel, lower_bound=q.ensemble_negativity_lower, upper_bound=q.negativity_upper_mean)
            else:
                lows, highs = [], []
                for family, color in zip(FAMILIES, COLORS):
                    q = subset(summary, family, t_star, 1200)
                    y = q.power_over_t_star_mean.to_numpy()
                    se = q.power_over_t_star_sem.to_numpy()
                    lows.append(np.min(y-se)); highs.append(np.max(y+se))
                    ax.plot(q.t, y, color=color, lw=1.2, label=family.title())
                    ax.fill_between(q.t,y-se,y+se,color=color,alpha=.14,linewidth=0)
                    ax.errorbar(q.t,y,yerr=se,fmt='none',ecolor=color,elinewidth=.65,capsize=1.6,errorevery=25)
                    record(q, panel, mean=y, standard_error=se)
                lower, upper = float(min(lows)), float(max(highs))
                gap = max((upper-lower)*.08, .001)
                ylim = (float(np.floor((lower-gap)/.005)*.005),
                        float(np.ceil((upper+gap)/.005)*.005))
                ax.set(xlim=(t_star,1200),ylim=ylim,xticks=[t_star,600,800,1000,1200])
                ax.set_xticklabels([r'$t_*$', '600','800','1000','1200'])
                ticks=np.arange(np.ceil(ylim[0]*100)/100,ylim[1]+1e-9,.01)
                ax.set_yticks(ticks,[f'{v:.2f}' for v in ticks])
                axes_labels(ax,r'$P_{\mathcal{B}}(t)/P_{\mathcal{B}}(t_*)$')
                ax.legend(loc='lower right',fontsize=8,handlelength=1.2,borderaxespad=.35,labelspacing=.17)
            save(fig, panel, ((0,2.48,4.95)[index],6.15,2.48,2.1))

        xs, ys = (1.02,2.32,3.62,4.92,6.22), (4.46,2.92,1.38)
        for row, family in enumerate(FAMILIES):
            run=selected[family]
            with np.load(source/'trajectories'/f'{run}.npz',allow_pickle=False) as z:
                for col,time in enumerate(times):
                    indices=np.flatnonzero(np.isclose(z['t'],time,rtol=0,atol=1e-11))
                    if len(indices)!=1: raise ValueError(f'{run}: missing saved time {time}.')
                    i=int(indices[0]);u=float(z['eps'])*(z['n'][i].astype(float)-.5*(z['s'][i]==0))
                    fig=plt.figure(figsize=(1.42,1.48))
                    ax=fig.add_axes([.25/1.42,.18/1.48,1.12/1.42,1.22/1.48])
                    ax.imshow(u,origin='lower',extent=(0,1,0,1),cmap='RdBu_r',vmin=-1,vmax=1,aspect='auto',interpolation='bilinear')
                    ax.set_xticks([0,.5,1],['0','0.5','1'])
                    ax.set_yticks([0,.5,1],['0','0.5','1'] if col==0 else [])
                    ax.tick_params(length=2.6,width=.65,pad=2,labelsize=8.5)
                    save(fig,f'map_{row}_{col}',(xs[col]-.25,ys[row]-.18,1.42,1.48))
                    map_rows.append({'id':run,'family':family,'t':float(time),'h_mean':float(z['h'][i].mean()),
                                     'contrast_sha256':hashlib.sha256(u.tobytes()).hexdigest()})

        fig=plt.figure(figsize=(WIDTH,HEIGHT))
        def label(x,y,text,**kwargs):
            fig.text(x/WIDTH,y/HEIGHT,text,**kwargs)
        for x,title in zip((.58,3.06,5.53),('(a) Environmental response','(b) Local entanglement','(c) Classical persistence')):
            label(x-.08,7.94,title,ha='left',va='center',fontsize=9.4)
        label(.48,6.00,'(d) Evolution of spatial configurations',ha='left',va='center',fontsize=9.4)
        for row,family in enumerate(FAMILIES):
            label(.64,ys[row]+.61,family.title(),ha='center',va='center',rotation=90,fontsize=9.3)
        for col,time in enumerate(times):
            title=r'$t=t_*$' if time==t_star else rf'$t={time:g}$'
            label(xs[col]+.56,5.81,title,ha='center',va='center',fontsize=9.7)
        label(.26,3.53,'Normalized y',ha='center',va='center',rotation=90,fontsize=9.5)
        label(4.18,1.12,'Normalized x',ha='center',va='center',fontsize=9.7)
        save(fig,'labels',(0,0,WIDTH,HEIGHT))
        fig=plt.figure(figsize=(3.4,.70))
        ax=fig.add_axes([.35/3.4,.43/.70,2.47/3.4,.12/.70])
        cb=fig.colorbar(ScalarMappable(norm=Normalize(-1,1),cmap='RdBu_r'),cax=ax,orientation='horizontal',ticks=[-1,0,1])
        cb.ax.tick_params(labelsize=9.1,length=2.6,pad=2,width=.7)
        cb.outline.set_linewidth(.7)
        cb.set_label('Occupation contrast $u$',labelpad=2.6,fontsize=9.7)
        save(fig,'colorbar',(2.30,.34,3.4,.70))
        pdf=output/'figure5_classicalization.pdf'
        with fitz.open() as document:
            page=document.new_page(width=WIDTH*72,height=HEIGHT*72)
            for path,(x,y,w,h) in parts:
                with fitz.open(path) as piece:
                    page.show_pdf_page(fitz.Rect(x*72,(HEIGHT-y-h)*72,(x+w)*72,(HEIGHT-y)*72),piece,0,keep_proportion=False)
            document.set_metadata({'title':'Autonomous classicalization','author':'Kazuki Ikeda',
                                   'subject':'Saved ensemble data with completion-referenced persistence'})
            document.save(pdf,garbage=4,deflate=True)
            page.get_pixmap(matrix=fitz.Matrix(300/72,300/72),alpha=False).save(output/'figure5_classicalization.png')
    pd.DataFrame(curve_rows).to_csv(output/'figure5_values.csv',index=False)
    info={'reference_time':reference,'representative_selection':'Lower median final band power, without morphology filtering',
          'runs':selected,'times':list(times),'heatmap':{'cmap':'RdBu_r','vmin':-1,'vmax':1,'interpolation':'bilinear'},
          'panel_c_ylim':list(ylim),'panel_c_mean_plus_minus_sem_range':[lower,upper],
          'maps':map_rows,'panel_layouts':artists,
          'source_sha256':{str(p.relative_to(source)):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [source/'analysis/summary.csv',source/'analysis/endpoints.csv',source/'analysis/reference_time.json']}}
    (output/'figure5_provenance.json').write_text(json.dumps(info,indent=2)+'\n')
    print(f'Figure saved to {pdf}; t_*={t_star:g}.')

if __name__=='__main__':
    main()
