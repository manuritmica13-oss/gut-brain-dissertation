#!/usr/bin/env python3
"""Plot shared taxa bar chart between MDD and PD.

Reads `results/shared_taxa_mdd_pd.tsv` and creates a horizontal grouped
bar chart of log2FC values (MDD vs PD), ordered by MDD log2FC.

Excludes taxa with abs(mdd_log2FC) > 10.
Bars: MDD = blue, PD = red. Bar edge colour indicates direction: depleted=teal, elevated=coral.
Adds vertical line at 0 and optional FDR significance markers if `mdd_padj`/`pd_padj`
columns are present in the TSV.

Saves output to `results/shared_taxa_barplot.png` at 300 dpi.
"""

from __future__ import annotations

import csv
from pathlib import Path
import math
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
IN_TSV = ROOT / 'results' / 'shared_taxa_mdd_pd.tsv'
OUT_PNG = ROOT / 'results' / 'shared_taxa_barplot.png'

if not IN_TSV.exists():
    print(f"Input not found: {IN_TSV}")
    sys.exit(1)

# Read TSV
rows = []
with IN_TSV.open('r', encoding='utf-8') as fh:
    reader = csv.DictReader(fh, delimiter='\t')
    for r in reader:
        rows.append(r)

if not rows:
    print('No rows in input TSV')
    sys.exit(0)

# parse values and optional padj columns
parsed = []
for r in rows:
    try:
        mdd = float(r.get('mdd_log2FC') or r.get('mdd_log2fc') or r.get('mdd') or 0.0)
    except Exception:
        mdd = None
    try:
        pd = float(r.get('pd_log2FC') or r.get('pd_log2fc') or r.get('pd') or 0.0)
    except Exception:
        pd = None
    direction = r.get('direction') or r.get('match_level') or ''
    match_level = r.get('match_level') or ''
    # optional padj columns
    try:
        mdd_padj = float(r.get('mdd_padj')) if r.get('mdd_padj') not in (None, '') else None
    except Exception:
        mdd_padj = None
    try:
        pd_padj = float(r.get('pd_padj')) if r.get('pd_padj') not in (None, '') else None
    except Exception:
        pd_padj = None
    parsed.append({'taxon': r.get('taxon'), 'mdd': mdd, 'pd': pd, 'direction': direction, 'match_level': match_level, 'mdd_padj': mdd_padj, 'pd_padj': pd_padj})

# deduplicate taxa: keep one record per species, prefer species-level match over genus-level
dedup = {}
for p in parsed:
    tax = p['taxon'] or ''
    # use species name as key; if tax contains underscore assume Genus_species
    species = tax
    if '_' in tax:
        species = tax
    else:
        # single token; treat as genus-level label
        species = tax
    key = species.strip()
    if key == '':
        continue
    if key in dedup:
        # prefer existing if it's species-level; otherwise prefer species-level
        prev = dedup[key]
        if prev.get('match_level') == 'species':
            continue
        if p.get('match_level') == 'species':
            dedup[key] = p
        else:
            # keep existing
            continue
    else:
        dedup[key] = p

# filter out unstable MDD estimates
filtered = [p for p in dedup.values() if p['mdd'] is not None and abs(p['mdd']) <= 10]
if not filtered:
    print('No taxa remain after filtering abs(mdd_log2FC) > 10')
    sys.exit(0)

# order by MDD log2FC descending
filtered.sort(key=lambda x: x['mdd'], reverse=True)

n = len(filtered)
ys = list(range(n))
taxa = [p['taxon'] for p in filtered]
mdd_vals = [p['mdd'] for p in filtered]
pd_vals = [p['pd'] for p in filtered]
dirs = [p['direction'] for p in filtered]
mdd_padjs = [p['mdd_padj'] for p in filtered]
pd_padjs = [p['pd_padj'] for p in filtered]

# plotting
fig, ax = plt.subplots(figsize=(14, max(4, n * 0.35)))
barh = 0.35
import numpy as np
y = np.arange(n)

# colours
color_mdd = '#1f77b4'  # blue
color_pd = '#d62728'   # red
dir_colors = {'depleted': '#008080', 'elevated': '#ff7f50'}  # teal, coral

# bars: MDD left, PD right (grouped horizontally)
ax.barh(y - barh/2, mdd_vals, height=barh, color=color_mdd, edgecolor=[dir_colors.get(d, 'black') for d in dirs], linewidth=1)
ax.barh(y + barh/2, pd_vals, height=barh, color=color_pd, edgecolor=[dir_colors.get(d, 'black') for d in dirs], linewidth=1)

# vertical line at 0
ax.axvline(0, color='grey', linewidth=0.8)

# significance markers
for i in range(n):
    # MDD
    if mdd_padjs[i] is not None and mdd_padjs[i] < 0.05:
        x = mdd_vals[i]
        ax.text(x + (0.2 if x >= 0 else -0.2), y[i] - barh/2, '*', va='center', ha=('left' if x>=0 else 'right'), color='black', fontsize=10)
    # PD
    if pd_padjs[i] is not None and pd_padjs[i] < 0.05:
        x = pd_vals[i]
        ax.text(x + (0.2 if x >= 0 else -0.2), y[i] + barh/2, '*', va='center', ha=('left' if x>=0 else 'right'), color='black', fontsize=10)

# y ticks — italicise species names (replace underscore with space)
yticklabels = []
for t in taxa:
    if t is None:
        yticklabels.append('')
        continue
    label = str(t).replace('_', ' ')
    # italicise using math text
    it_label = r"$\it{" + label + r"}$"
    yticklabels.append(it_label)

ax.set_yticks(y)
ax.set_yticklabels(yticklabels, fontsize=10)
ax.invert_yaxis()

# legend: create custom handles
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
handles = [Patch(facecolor=color_mdd, edgecolor='none', label='MDD (log2FC)'), Patch(facecolor=color_pd, edgecolor='none', label='PD (log2FC)')]
dir_handles = [Patch(facecolor='none', edgecolor=dir_colors['depleted'], label='depleted (both)'), Patch(facecolor='none', edgecolor=dir_colors['elevated'], label='elevated (both)')]
ax.legend(handles=handles + dir_handles, loc='lower right')

ax.set_xlabel('log2 fold change')
ax.set_title('Shared taxa: MDD vs PD')
plt.tight_layout()
OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(OUT_PNG, dpi=300)
print(f'Saved plot to {OUT_PNG}')
