#!/usr/bin/env python3
"""Compare MDD and PD differential abundance results.

This script implements the requested steps without external DataFrame
dependencies: it uses openpyxl and the stdlib.
"""

import csv
import re
from pathlib import Path
import sys
from io import StringIO

import openpyxl

ROOT = Path(__file__).resolve().parents[2]
MDD_XLSX = ROOT / "data" / "raw" / "MDD_2026_zenodo" / "Table_S4_DAAResultsDESeq2.xlsx"
PD_TXT = ROOT / "results" / "phase2_wallen_results.txt"
OUT = ROOT / "results" / "shared_taxa_mdd_pd.tsv"
OUT.parent.mkdir(parents=True, exist_ok=True)


def normalize_key(s: str) -> str:
    if s is None:
        return ""
    return re.sub(r"[^a-z0-9]+", "", str(s).lower())


def clean_pd_taxon(t: str) -> str:
    if t is None:
        return ""
    s = str(t).strip()
    if '|' in s:
        parts = [p.strip() for p in s.split('|') if p.strip()]
        last = parts[-1]
        if '__' in last:
            return last.split('__', 1)[1].strip()
        return last
    if '__' in s:
        return s.split('__', 1)[1].strip()
    return s


def read_mdd():
    """Read MDD Excel sheet 'D_vs_C_alone_DeSEQ2' skipping first 2 rows.
    Returns dict keyed by normalized taxon -> {taxon, log2fc, padj} for padj<0.05.
    """
    if not MDD_XLSX.exists():
        print(f"MDD file not found: {MDD_XLSX}")
        sys.exit(1)
    wb = openpyxl.load_workbook(MDD_XLSX, read_only=True, data_only=True)
    sheet_name = 'D_vs_C_alone_DeSEQ2'
    ws = wb[sheet_name] if sheet_name in wb.sheetnames else wb.active

    rows = list(ws.iter_rows(min_row=3, values_only=True))
    if not rows:
        return {}
    header = [str(x).strip().lower() if x is not None else '' for x in rows[0]]

    def idx(colname):
        try:
            return header.index(colname.lower())
        except ValueError:
            return None

    tax_i = idx('taxon') or 0
    l2_i = idx('log2foldchange')
    padj_i = idx('padj')

    out = {}
    for r in rows[1:]:
        if not r:
            continue
        tax = r[tax_i] if tax_i is not None and tax_i < len(r) else None
        padj = r[padj_i] if padj_i is not None and padj_i < len(r) else None
        l2 = r[l2_i] if l2_i is not None and l2_i < len(r) else None
        try:
            padj_v = float(padj) if padj not in (None, '') else None
        except Exception:
            padj_v = None
        try:
            l2_v = float(l2) if l2 not in (None, '') else None
        except Exception:
            l2_v = None

        if padj_v is not None and padj_v < 0.05:
            key = normalize_key(str(tax))
            out[key] = {'taxon': str(tax).strip(), 'log2fc': l2_v, 'padj': padj_v}
    return out


def read_pd():
    """Parse the PD results file and return dict keyed by normalized lowest-level taxon."""
    if not PD_TXT.exists():
        print(f"PD results file not found: {PD_TXT}")
        sys.exit(1)
    txt = PD_TXT.read_text(encoding='utf-8')
    marker = 'Top 20 DESeq2 results (by padj):'
    i = txt.find(marker)
    if i == -1:
        print('PD results marker not found in phase2_wallen_results.txt')
        return {}
    rest = txt[i + len(marker):].lstrip('\n')
    lines = rest.splitlines()
    table_lines = []
    for ln in lines:
        if ln.strip() == '':
            break
        table_lines.append(ln)
    if not table_lines:
        return {}
    table_txt = '\n'.join(table_lines)
    reader = csv.DictReader(StringIO(table_txt), delimiter='\t')
    out = {}
    for row in reader:
        taxa_raw = row.get('taxa') or row.get('Taxa') or row.get('taxon')
        padj = row.get('padj') or row.get('padj ') or row.get('padj')
        l2 = row.get('log2FC') or row.get('log2fc') or row.get('log2')
        try:
            padj_v = float(padj) if padj not in (None, '') else None
        except Exception:
            padj_v = None
        try:
            l2_v = float(l2) if l2 not in (None, '') else None
        except Exception:
            l2_v = None

        if padj_v is not None and padj_v < 0.05:
            lowest = clean_pd_taxon(taxa_raw)
            key = normalize_key(lowest)
            out[key] = {'taxon': lowest, 'log2fc': l2_v, 'padj': padj_v}
    return out


def main():
    mdd = read_mdd()
    pd = read_pd()

    # prepare MDD lookup maps: species -> record, genus -> list of records
    mdd_species = {}
    mdd_genus = {}
    for k, v in mdd.items():
        # expect MDD taxon like Genus_species or Genus
        name = str(v['taxon'])
        norm_species = normalize_key(name)
        mdd_species[norm_species] = v
        genus = name.split('_')[0] if '_' in name else name
        norm_genus = normalize_key(genus)
        mdd_genus.setdefault(norm_genus, []).append(v)

    tested_mdd = len(mdd_species)
    tested_pd = len(pd)

    shared = []
    species_matches = 0
    genus_matches = 0

    for pk, pv in pd.items():
        pd_taxon = pv.get('taxon')
        if not pd_taxon:
            continue
        # extract species and genus from PD taxon
        pd_species = None
        pd_genus = None
        s = str(pd_taxon)
        # if input contains |, prefer the species part if present
        if '|' in s:
            parts = [p.strip() for p in s.split('|') if p.strip()]
            # look for species (prefix s__)
            species_part = next((p for p in parts if p.startswith('s__')), None)
            genus_part = next((p for p in parts if p.startswith('g__')), None)
            if species_part:
                pd_species = species_part.split('__', 1)[1]
            if genus_part:
                pd_genus = genus_part.split('__', 1)[1]
            # fallback: if no genus_part but species_part exists, derive genus
            if not pd_genus and pd_species and '_' in pd_species:
                pd_genus = pd_species.split('_', 1)[0]
        else:
            # single token like s__Genus_species or g__Genus
            if '__' in s:
                pref, val = s.split('__', 1)
                if pref == 's':
                    pd_species = val
                    if '_' in val:
                        pd_genus = val.split('_', 1)[0]
                elif pref == 'g':
                    pd_genus = val
            else:
                # unknown format, treat as a species-like string
                if '_' in s:
                    pd_species = s
                    pd_genus = s.split('_', 1)[0]
                else:
                    pd_genus = s

        matched = False
        # 1) exact species match (only if pd_species available)
        if pd_species:
            key = normalize_key(pd_species)
            if key in mdd_species:
                mv = mdd_species[key]
                if mv['log2fc'] is None or pv['log2fc'] is None:
                    matched = False
                else:
                    if mv['log2fc'] < 0 and pv['log2fc'] < 0:
                        direction = 'depleted'
                    elif mv['log2fc'] > 0 and pv['log2fc'] > 0:
                        direction = 'elevated'
                    else:
                        matched = False
                    if not matched:
                        # if direction set above, mark matched
                        pass
                    matched = ('direction' in locals() and direction in ('depleted', 'elevated'))
                    if matched:
                        species_matches += 1
                        shared.append({'taxon': mv['taxon'], 'mdd_log2FC': mv['log2fc'], 'pd_log2FC': pv['log2fc'], 'direction': direction, 'match_level': 'species'})
        # 2) genus-level match when species not available/matched
        if not matched and pd_genus:
            gkey = normalize_key(pd_genus)
            if gkey in mdd_genus:
                # pick any matching MDD taxon in that genus; prefer exact genus hit
                candidates = mdd_genus[gkey]
                # find a candidate with non-null log2fc
                chosen = None
                for c in candidates:
                    if c.get('log2fc') is not None:
                        chosen = c
                        break
                if chosen and chosen['log2fc'] is not None and pv.get('log2fc') is not None:
                    if chosen['log2fc'] < 0 and pv['log2fc'] < 0:
                        direction = 'depleted'
                    elif chosen['log2fc'] > 0 and pv['log2fc'] > 0:
                        direction = 'elevated'
                    else:
                        chosen = None
                if chosen:
                    genus_matches += 1
                    shared.append({'taxon': chosen['taxon'], 'mdd_log2FC': chosen['log2fc'], 'pd_log2FC': pv['log2fc'], 'direction': direction, 'match_level': 'genus'})

    # print summary
    print(f"Tested MDD taxa: {tested_mdd}")
    print(f"Tested PD taxa: {tested_pd}")
    print(f"Species-level matches: {species_matches}")
    print(f"Genus-level matches: {genus_matches}")

    if shared:
        print('Shared taxa examples:')
        for s in shared[:20]:
            print(f"- {s['taxon']}: {s['direction']} ({s['match_level']}) MDD {s['mdd_log2FC']}, PD {s['pd_log2FC']}")
    else:
        print('No shared taxa found.')

    # write TSV
    with OUT.open('w', encoding='utf-8', newline='') as fh:
        fieldnames = ['taxon', 'match_level', 'mdd_log2FC', 'pd_log2FC', 'direction']
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for s in shared:
            writer.writerow({k: s.get(k, '') for k in fieldnames})
    print(f"Saved {len(shared)} shared taxa to {OUT}")


if __name__ == '__main__':
    main()
