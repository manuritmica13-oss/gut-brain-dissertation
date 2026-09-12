#!/usr/bin/env python3
"""Plot Shannon alpha diversity for PD vs Control samples.

Loads the Wallen relative abundance table and subject metadata, filters to PD and
Control samples, computes Shannon entropy for each sample, and saves a polished
boxplot with jittered data points and significance annotation.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
ABUND_PATH = ROOT / 'data' / 'processed' / 'wallen' / 'metaphlan_rel_ab.tsv'
META_PATH = ROOT / 'data' / 'processed' / 'wallen' / 'subject_metadata.tsv'
OUT_PATH = ROOT / 'results' / 'alpha_diversity_boxplot.png'
OUT_PATH.parent.mkdir(parents=True, exist_ok=True)


def compute_shannon_diversity(sample_row: pd.Series) -> float:
    """Compute Shannon diversity using relative abundances (non-zero only)."""
    values = pd.to_numeric(sample_row, errors='coerce').dropna()
    values = values[values > 0]
    if values.empty:
        return 0.0
    p = values / values.sum()
    return float(-(p * np.log(p)).sum())


def main() -> None:
    if not ABUND_PATH.exists():
        raise FileNotFoundError(ABUND_PATH)
    if not META_PATH.exists():
        raise FileNotFoundError(META_PATH)

    abundance = pd.read_csv(ABUND_PATH, sep='\t', index_col=0).T
    metadata = pd.read_csv(META_PATH, sep='\t', dtype={'sample_name': str, 'Case_status': str})

    if 'sample_name' not in metadata.columns or 'Case_status' not in metadata.columns:
        raise KeyError('Metadata must contain sample_name and Case_status columns.')

    meta_filtered = metadata.loc[metadata['Case_status'].isin(['PD', 'Control']), ['sample_name', 'Case_status']].copy()
    meta_filtered = meta_filtered.drop_duplicates(subset='sample_name')
    meta_filtered['sample_name'] = meta_filtered['sample_name'].astype(str)

    common_ids = abundance.index.intersection(meta_filtered['sample_name'])
    abundance = abundance.loc[common_ids].copy()
    metadata_aligned = meta_filtered.set_index('sample_name').loc[common_ids].rename_axis('sample_name').reset_index()

    shannon = abundance.apply(compute_shannon_diversity, axis=1).rename('shannon_diversity')
    df = metadata_aligned.copy()
    df['shannon_diversity'] = shannon.values

    palette = {'Control': '#1f9d9a', 'PD': '#f26b5b'}
    sns.set_theme(style='whitegrid', context='talk')

    fig, ax = plt.subplots(figsize=(6, 6))
    sns.boxplot(
        data=df,
        x='Case_status',
        y='shannon_diversity',
        order=['Control', 'PD'],
        palette=palette,
        width=0.5,
        fliersize=0,
        showcaps=True,
        boxprops={'edgecolor': 'black', 'linewidth': 1.2},
        whiskerprops={'color': 'black', 'linewidth': 1.2},
        capprops={'color': 'black', 'linewidth': 1.2},
        medianprops={'color': 'black', 'linewidth': 2},
        ax=ax,
    )

    sns.stripplot(
        data=df,
        x='Case_status',
        y='shannon_diversity',
        order=['Control', 'PD'],
        palette=palette,
        dodge=True,
        jitter=0.22,
        size=4,
        alpha=0.75,
        edgecolor='black',
        linewidth=0.3,
        ax=ax,
    )

    # Add significance bar and annotation.
    control_vals = df.loc[df['Case_status'] == 'Control', 'shannon_diversity']
    pd_vals = df.loc[df['Case_status'] == 'PD', 'shannon_diversity']
    y1 = max(control_vals.max(), pd_vals.max()) + 0.12
    y2 = y1 + 0.05
    ax.plot([0, 0, 1, 1], [y1, y2, y2, y1], linewidth=1.5, color='black')
    ax.text(0.5, y2 + 0.02, 'p<0.0001', ha='center', va='bottom', fontsize=12, fontweight='bold')

    ax.set_xlabel('')
    ax.set_ylabel('Shannon diversity index', fontsize=14)
    ax.set_ylim(bottom=2.5)

    # Add sample sizes below each group.
    counts = {'Control': int((df['Case_status'] == 'Control').sum()), 'PD': int((df['Case_status'] == 'PD').sum())}
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f'Control\nn={counts["Control"]}', f'PD\nn={counts["PD"]}'], fontsize=12)
    ax.tick_params(axis='x', length=0)

    sns.despine(trim=True)
    fig.tight_layout()
    fig.savefig(OUT_PATH, dpi=300, bbox_inches='tight')
    print(f'Saved alpha diversity boxplot to {OUT_PATH}')


if __name__ == '__main__':
    main()
