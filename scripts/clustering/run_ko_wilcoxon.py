#!/usr/bin/env python3
"""Run Wilcoxon rank-sum tests for key gut-brain KO genes.

This script:
1) loads the KO abundance table (rows = enzymes, columns = samples)
2) loads the subject metadata and keeps samples with known PD / Control status
3) transposes to samples x enzymes
4) compares PD vs Control abundance for IDO1 (K00463) and KMO (K00486)
5) reports mean abundance, fold change, and Wilcoxon rank-sum p-values
6) saves the summary to results/ko_wilcoxon_results.txt
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parents[2]
KO_PATH = ROOT / "data" / "processed" / "wallen" / "humann_ko_counts.tsv"
META_PATH = ROOT / "data" / "processed" / "wallen" / "subject_metadata.tsv"
OUT_PATH = ROOT / "results" / "ko_wilcoxon_results.txt"
OUT_PATH.parent.mkdir(parents=True, exist_ok=True)

TARGETS = {
    "IDO1": "K00463",
    "KMO": "K00486",
}


def load_ko_table(path: Path) -> pd.DataFrame:
    """Read KO abundances and return a samples x enzyme DataFrame."""
    if not path.exists():
        raise FileNotFoundError(path)

    table = pd.read_csv(path, sep='\t', index_col=0)
    table = table.T
    table.index = table.index.map(str)
    table.columns = table.columns.map(str)
    return table


def load_metadata(path: Path) -> pd.DataFrame:
    """Return metadata filtered to PD and Control samples."""
    if not path.exists():
        raise FileNotFoundError(path)

    meta = pd.read_csv(path, sep='\t', dtype={"sample_name": str, "Case_status": str})
    required = {"sample_name", "Case_status"}
    missing = required - set(meta.columns)
    if missing:
        raise KeyError(f"Missing metadata columns: {sorted(missing)}")

    meta = meta.loc[meta["Case_status"].isin(["PD", "Control"]), ["sample_name", "Case_status"]].copy()
    meta = meta.drop_duplicates(subset="sample_name")
    meta["sample_name"] = meta["sample_name"].map(str)
    return meta.set_index("sample_name")


def pick_enzyme_column(table: pd.DataFrame, accession: str) -> str:
    """Find the column corresponding to the requested KO accession."""
    matches = []
    for col in table.columns:
        if not isinstance(col, str):
            continue
        cleaned = col.strip()
        if cleaned == accession:
            matches.append(cleaned)
            continue
        if cleaned.startswith(f"{accession}:"):
            matches.append(cleaned)
            continue
        if ":" in cleaned and cleaned.split(":", 1)[0] == accession:
            matches.append(cleaned)
            continue
    if not matches:
        raise KeyError(f"Could not find KO {accession} in abundance table")
    return matches[0]


def safe_fold_change(pd_mean: float, control_mean: float) -> float:
    """Compute PD/control fold change while handling zero denominators."""
    if pd_mean == 0 and control_mean == 0:
        return 1.0
    if control_mean == 0:
        return float("inf") if pd_mean > 0 else 0.0
    return float(pd_mean / control_mean)


def run_wilcoxon(x: pd.Series, y: pd.Series) -> tuple[float, float, float]:
    """Run the Mann-Whitney U test and return statistic, p-value, and rank-biserial effect size."""
    if len(x) == 0 or len(y) == 0:
        return float("nan"), float("nan"), float("nan")

    # Guard against degenerate samples; scipy will often return p=1 for all-equal values,
    # but this explicit check keeps the output stable and interpretable.
    if x.nunique() <= 1 and y.nunique() <= 1 and x.mean() == y.mean():
        return 0.0, 1.0, 0.0

    stat, pval = mannwhitneyu(x, y, alternative="two-sided")
    n1 = len(x)
    n2 = len(y)
    if n1 > 0 and n2 > 0:
        effect_size = (2.0 * stat / (n1 * n2)) - 1.0
    else:
        effect_size = float("nan")
    return float(stat), float(pval), float(effect_size)


def main() -> None:
    ko_table = load_ko_table(KO_PATH)
    metadata = load_metadata(META_PATH)

    common_samples = ko_table.index.intersection(metadata.index)
    ko_table = ko_table.loc[common_samples].copy()
    metadata = metadata.loc[common_samples].copy()

    if ko_table.empty:
        raise ValueError("No overlapping samples remain between KO table and metadata")

    # Convert raw counts to copies per million (CPM) within each sample.
    row_totals = ko_table.sum(axis=1).replace(0, np.nan)
    ko_table = ko_table.div(row_totals, axis=0).fillna(0.0) * 1_000_000

    results = []
    for gene_name, accession in TARGETS.items():
        try:
            enzyme_col = pick_enzyme_column(ko_table, accession)
        except KeyError as exc:
            print(f"{gene_name} ({accession}): {exc}")
            results.append({
                "gene": gene_name,
                "accession": accession,
                "status": "not_found",
                "message": str(exc),
            })
            continue

        # Align abundance values to sample metadata and keep only known PD / Control labels.
        abundance = ko_table[[enzyme_col]].copy()
        abundance["Case_status"] = metadata["Case_status"].values
        abundance = abundance.dropna(subset=[enzyme_col])

        pd_vals = pd.to_numeric(
            abundance.loc[abundance["Case_status"] == "PD", enzyme_col],
            errors="coerce",
        ).dropna()
        ctrl_vals = pd.to_numeric(
            abundance.loc[abundance["Case_status"] == "Control", enzyme_col],
            errors="coerce",
        ).dropna()

        pd_mean = float(pd_vals.mean()) if not pd_vals.empty else float("nan")
        ctrl_mean = float(ctrl_vals.mean()) if not ctrl_vals.empty else float("nan")
        fold_change = safe_fold_change(pd_mean, ctrl_mean) if pd_mean == pd_mean and ctrl_mean == ctrl_mean else float("nan")
        statistic, p_value, effect_size = run_wilcoxon(pd_vals, ctrl_vals)

        record = {
            "gene": gene_name,
            "accession": accession,
            "column": enzyme_col,
            "n_pd": int(len(pd_vals)),
            "n_control": int(len(ctrl_vals)),
            "pd_mean": pd_mean,
            "control_mean": ctrl_mean,
            "fold_change_pd_vs_control": fold_change,
            "wilcoxon_statistic": statistic,
            "p_value": p_value,
            "mann_whitney_effect_size": effect_size,
        }
        results.append(record)

        print(f"{gene_name} ({accession})")
        print(f"  PD mean abundance: {pd_mean}")
        print(f"  Control mean abundance: {ctrl_mean}")
        print(f"  Fold change (PD / Control): {fold_change}")
        print(f"  Wilcoxon rank-sum statistic: {statistic}")
        print(f"  p-value: {p_value}")
        print(f"  Mann-Whitney effect size (rank-biserial): {effect_size}")
        print("")

    out_lines = [
        "KO Wilcoxon rank-sum results",
        "===========================",
        "",
    ]
    for record in results:
        if record.get("status") == "not_found":
            out_lines.append(f"{record['gene']} ({record['accession']}): {record['message']}")
            continue

        out_lines.append(
            f"{record['gene']} ({record['accession']})\n"
            f"  column: {record['column']}\n"
            f"  PD mean abundance (CPM): {record['pd_mean']}\n"
            f"  Control mean abundance (CPM): {record['control_mean']}\n"
            f"  Fold change (PD / Control): {record['fold_change_pd_vs_control']}\n"
            f"  Wilcoxon statistic: {record['wilcoxon_statistic']}\n"
            f"  p-value: {record['p_value']}\n"
            f"  Mann-Whitney effect size (rank-biserial): {record['mann_whitney_effect_size']}\n"
            f"  n_PD: {record['n_pd']}\n"
            f"  n_Control: {record['n_control']}"
        )
        out_lines.append("")

    OUT_PATH.write_text("\n".join(out_lines).rstrip() + "\n", encoding="utf-8")
    print(f"Saved KO Wilcoxon results to {OUT_PATH}")


if __name__ == "__main__":
    main()
