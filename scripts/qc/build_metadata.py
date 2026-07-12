#!/usr/bin/env python3
"""Build a combined metadata table from SRA run metadata files."""

from __future__ import annotations

import csv
import re
from collections import Counter
from pathlib import Path
from typing import Dict, List

import openpyxl

ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = ROOT / "data" / "raw"
OUTPUT_FILE = RAW_DIR / "metadata_combined.tsv"

STUDY_FOLDERS = [
    "PRJNA834801",
    "PRJEB55464",
    "PRJNA808166",
    "PRJNA588035",
    "PRJNA531273",
    "PRJNA762199",
    "PRJNA1083304",
    "PRJNA943232",
]

PD_STUDIES = {"PRJNA834801", "PRJEB55464", "PRJNA808166", "PRJNA588035", "PRJNA762199"}
MDD_STUDIES = {"PRJNA531273", "PRJNA1083304", "PRJNA943232"}

OUTPUT_COLUMNS = [
    "sample_id",
    "study",
    "disease",
    "disease_group",
    "disease_status",
    "age",
    "sex",
    "bmi",
    "medication_status",
]


def normalize_name(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", value.lower())


def get_value(row: Dict[str, str], aliases: List[str]) -> str:
    normalized_row = {normalize_name(k): v for k, v in row.items() if k is not None}
    for alias in aliases:
        key = normalize_name(alias)
        if key in normalized_row:
            value = (normalized_row[key] or "").strip()
            if value:
                return value
    return ""


def summarize_missingness(rows: List[Dict[str, str]]) -> List[Dict[str, object]]:
    if not rows:
        return []

    summaries = []
    for column in OUTPUT_COLUMNS:
        missing = sum(1 for row in rows if not (row.get(column) or "").strip())
        summaries.append(
            {
                "column": column,
                "missing_count": missing,
                "missing_fraction": round(missing / len(rows), 3),
            }
        )
    return summaries


def normalize_disease_status(value: str) -> str:
    if value is None:
        return "unknown"

    normalized = str(value).strip()
    if not normalized:
        return "unknown"

    lowered = normalized.lower()
    if any(token in lowered for token in ["control", "ctrl", "hc", "healthy", "health", "normal"]):
        return "Control"
    if any(token in lowered for token in ["case", "patient", "pd", "parkinson", "mdd", "depression", "depressive", "disease", "disorder", "illness"]):
        return "Case"
    if lowered in {"unknown", "na", "n/a", "none", "null"}:
        return "unknown"
    return normalized


def infer_disease_status(row: Dict[str, str]) -> str:
    text_fragments: List[str] = []
    for key, value in row.items():
        if not value:
            continue
        normalized_key = normalize_name(key)
        if any(token in normalized_key for token in [
            "samplename",
            "sample",
            "biosample",
            "subject",
            "run",
            "experiment",
            "description",
            "hostdisease",
            "subjectstatus",
            "disease",
            "affection",
        ]):
            text_fragments.append(str(value))

    combined_text = " ".join(text_fragments).lower()
    if re.search(r"\b(control|ctrl|hc|healthy|health|normal)\b", combined_text):
        return "control"
    if re.search(r"\b(patient|pd|parkinson|mdd|depression|depressive|case|disease|disorder|illness)\b", combined_text):
        return "case"
    return ""


def load_subject_metadata(path: Path) -> Dict[str, Dict[str, str]]:
    if not path.exists():
        return {}

    workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
    worksheet = workbook["subject_metadata"]
    rows = list(worksheet.rows)
    if not rows:
        return {}

    header = [cell.value for cell in rows[0]]
    header_map = {normalize_name(name): idx for idx, name in enumerate(header) if name is not None}
    metadata_by_sample: Dict[str, Dict[str, str]] = {}

    for row in rows[1:]:
        values = [cell.value for cell in row]
        sample_name = ""
        sample_idx = header_map.get(normalize_name("sample_name"))
        if sample_idx is not None and sample_idx < len(values):
            sample_name = str(values[sample_idx] or "").strip()
        if not sample_name:
            continue

        disease_status = ""
        sex = ""
        age = ""
        bmi = ""
        medication_status = ""

        for field_name, target in [
            ("Case_status", "disease_status"),
            ("Sex", "sex"),
            ("Age_at_collection", "age"),
            ("BMI", "bmi"),
            ("Depression_anxiety_mood_med", "medication_status"),
        ]:
            idx = header_map.get(normalize_name(field_name))
            if idx is None or idx >= len(values):
                continue
            value = values[idx]
            if value is None:
                continue
            if target == "disease_status":
                disease_status = str(value).strip()
            elif target == "sex":
                sex = str(value).strip()
            elif target == "age":
                age = str(value).strip()
            elif target == "bmi":
                bmi = str(value).strip()
            elif target == "medication_status":
                if str(value).strip() == "Y":
                    medication_status = "medicated"
                elif str(value).strip() == "N":
                    medication_status = "naive"

        metadata_by_sample[sample_name] = {
            "disease_status": disease_status,
            "sex": sex,
            "age": age,
            "bmi": bmi,
            "medication_status": medication_status,
        }

    return metadata_by_sample


def merge_subject_metadata(rows: List[Dict[str, str]], metadata_by_sample: Dict[str, Dict[str, str]]) -> None:
    for row in rows:
        if row.get("study") != "PRJNA834801":
            continue
        sample_id = (row.get("sample_id") or "").strip()
        if not sample_id:
            continue
        subject_info = metadata_by_sample.get(sample_id)
        if not subject_info:
            continue
        for field_name, value in subject_info.items():
            current_value = (row.get(field_name) or "").strip()
            if field_name == "disease_status" and current_value and current_value.lower() != "unknown":
                continue
            if field_name != "disease_status" and current_value:
                continue
            row[field_name] = value


def extract_pd_token(value: str) -> str:
    if not value:
        return ""
    normalized = re.sub(r"[^a-z0-9]+", "", str(value).lower())
    match = re.search(r"pd(\d+)", normalized)
    if match:
        return f"pd{match.group(1)}"
    return ""


def load_prjna588035_metadata(path: Path) -> Dict[str, Dict[str, str]]:
    if not path.exists():
        return {}

    workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
    worksheet = workbook.active
    rows = list(worksheet.rows)
    if not rows:
        return {}

    header = [cell.value for cell in rows[0]]
    header_map = {normalize_name(name): idx for idx, name in enumerate(header) if name is not None}
    metadata_by_sample: Dict[str, Dict[str, str]] = {}

    for row in rows[1:]:
        values = [cell.value for cell in row]
        sample_name = ""
        sample_idx = header_map.get(normalize_name("sample"))
        if sample_idx is not None and sample_idx < len(values):
            sample_name = str(values[sample_idx] or "").strip()
        if not sample_name:
            continue

        sex = ""
        age = ""
        sex_idx = header_map.get(normalize_name("gender"))
        if sex_idx is not None and sex_idx < len(values):
            sex = str(values[sex_idx] or "").strip()
        age_idx = header_map.get(normalize_name("age"))
        if age_idx is not None and age_idx < len(values):
            age = str(values[age_idx] or "").strip()

        sample_token = extract_pd_token(sample_name)
        if sample_token:
            metadata_by_sample[sample_token] = {"sex": sex, "age": age}

    return metadata_by_sample


def enrich_prjna588035_metadata(rows: List[Dict[str, str]], metadata_by_sample: Dict[str, Dict[str, str]]) -> None:
    for row in rows:
        if row.get("study") != "PRJNA588035":
            continue

        sample_token = extract_pd_token(row.get("sample_id") or "")
        if not sample_token:
            row["disease_status"] = "Control"
            continue

        metadata = metadata_by_sample.get(sample_token)
        if not metadata:
            row["disease_status"] = "Control"
            continue

        if not (row.get("sex") or "").strip() and metadata.get("sex"):
            row["sex"] = metadata["sex"]
        if not (row.get("age") or "").strip() and metadata.get("age"):
            row["age"] = metadata["age"]
        row["disease_status"] = "Case"


def enrich_prjeb55464_metadata(rows: List[Dict[str, str]], metadata_path: Path) -> None:
    if not metadata_path.exists():
        return

    with metadata_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            sample_id = (row.get("SampleName") or row.get("Sample") or row.get("Sample_Name") or "").strip()
            if not sample_id:
                continue
            library_name = (row.get("LibraryName") or "").strip()
            sex = (row.get("Sex") or "").strip()
            for target_row in rows:
                if target_row.get("study") != "PRJEB55464":
                    continue
                if (target_row.get("sample_id") or "").strip() != sample_id:
                    continue
                current_status = (target_row.get("disease_status") or "").strip()
                if not current_status or current_status.lower() == "unknown":
                    if "HC" in library_name:
                        target_row["disease_status"] = "Control"
                    elif "PD" in library_name:
                        target_row["disease_status"] = "Case"
                if not (target_row.get("sex") or "").strip() and sex:
                    target_row["sex"] = sex


def enrich_study_specific_metadata(rows: List[Dict[str, str]], study_id: str, metadata_path: Path) -> None:
    if not metadata_path.exists():
        return

    with metadata_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            sample_id = (row.get("SampleName") or row.get("Sample") or row.get("Sample_Name") or "").strip()
            if not sample_id:
                continue

            if study_id == "PRJNA943232":
                identifier = (row.get("SampleName") or "").strip()
                if not (identifier.startswith("HC") or identifier.startswith("MDD")):
                    continue
                if not (identifier.startswith("HC")):
                    label = "Case"
                else:
                    label = "Control"
            elif study_id == "PRJNA762199":
                identifier = (row.get("LibraryName") or "").strip()
                if identifier.startswith("HC_"):
                    label = "Control"
                elif identifier.startswith("PD_"):
                    label = "Case"
                else:
                    continue
            elif study_id == "PRJNA1083304":
                identifier = (row.get("SampleName") or row.get("LibraryName") or "").strip()
                if identifier.startswith("HC"):
                    label = "Control"
                elif identifier.startswith("MD") or identifier.startswith("SX"):
                    label = "Case"
                else:
                    continue
            else:
                continue

            for target_row in rows:
                if target_row.get("study") != study_id:
                    continue
                if (target_row.get("sample_id") or "").strip() != sample_id:
                    continue
                current_status = (target_row.get("disease_status") or "").strip()
                if not current_status or current_status.lower() == "unknown":
                    target_row["disease_status"] = label


def main() -> None:
    combined_rows: List[Dict[str, str]] = []
    case_control_counts: Counter = Counter()

    for study_id in STUDY_FOLDERS:
        metadata_path = RAW_DIR / study_id / "SraRunInfo.csv"
        if not metadata_path.exists():
            print(f"Skipping missing file: {metadata_path}")
            continue

        with metadata_path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                disease_group = "PD" if study_id in PD_STUDIES else "MDD" if study_id in MDD_STUDIES else ""
                inferred_status = infer_disease_status(row)
                if inferred_status:
                    case_control_counts[(study_id, inferred_status)] += 1
                initial_status = inferred_status or get_value(row, ["disease_status", "disease", "affection_status", "condition"])
                initial_status = normalize_disease_status(initial_status)
                combined_rows.append(
                    {
                        "sample_id": get_value(row, ["sample_id", "samplename", "sample_name", "sample", "biosample", "run"]),
                        "study": get_value(row, ["study", "bioproject", "sra_study", "study_accession"]) or study_id,
                        "disease": disease_group,
                        "disease_group": disease_group,
                        "disease_status": initial_status,
                        "age": get_value(row, ["age", "age_years", "subject_age", "age_at_sampling"]),
                        "sex": get_value(row, ["sex", "gender"]),
                        "bmi": get_value(row, ["bmi", "body_mass_index", "bodymassindex"]),
                        "medication_status": get_value(
                            row,
                            [
                                "medication_status",
                                "medication",
                                "medication_use",
                                "treatment_status",
                                "treatment",
                            ],
                        ),
                    }
                )

    subject_metadata_path = RAW_DIR / "PRJNA834801" / "source_data" / "Source_Data_24Oct2022.xlsx"
    subject_metadata = load_subject_metadata(subject_metadata_path)
    merge_subject_metadata(combined_rows, subject_metadata)

    prjeb55464_metadata_path = RAW_DIR / "PRJEB55464" / "SraRunInfo.csv"
    enrich_prjeb55464_metadata(combined_rows, prjeb55464_metadata_path)

    for study_id in ["PRJNA943232", "PRJNA762199", "PRJNA1083304"]:
        enrich_study_specific_metadata(combined_rows, study_id, RAW_DIR / study_id / "SraRunInfo.csv")

    prjna588035_metadata_path = None
    for candidate_path in [
        Path("/Users/manucoro/Downloads/Table_1.XLSX"),
        Path("/mnt/user-data/uploads/Table_1.XLSX"),
        RAW_DIR / "PRJNA588035" / "Table_1.XLSX",
    ]:
        if candidate_path.exists():
            prjna588035_metadata_path = candidate_path
            break

    if prjna588035_metadata_path is not None:
        prjna588035_metadata = load_prjna588035_metadata(prjna588035_metadata_path)
        enrich_prjna588035_metadata(combined_rows, prjna588035_metadata)
    else:
        print("PRJNA588035 supplementary metadata file not found; leaving disease_status as unknown/control based on sample naming.")

    for row in combined_rows:
        if row.get("study") in {"PRJNA531273", "PRJNA808166"}:
            row["disease_status"] = "unknown"
        else:
            row["disease_status"] = normalize_disease_status(row.get("disease_status", ""))

    with OUTPUT_FILE.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(combined_rows)

    counts_by_study = Counter(row["study"] for row in combined_rows)
    print("Metadata summary")
    print("=" * 40)
    print("Samples per study with disease group:")
    for study_id in sorted(counts_by_study):
        disease_group = "PD" if study_id in PD_STUDIES else "MDD" if study_id in MDD_STUDIES else ""
        print(f"- {study_id}: {counts_by_study[study_id]} samples, disease_group={disease_group}")

    prjna834801_rows = [row for row in combined_rows if row.get("study") == "PRJNA834801"]
    filled_disease_status = sum(1 for row in prjna834801_rows if (row.get("disease_status") or "").strip())
    missing_disease_status = len(prjna834801_rows) - filled_disease_status
    print(f"PRJNA834801 disease_status filled: {filled_disease_status}/{len(prjna834801_rows)}")
    print(f"PRJNA834801 disease_status missing: {missing_disease_status}/{len(prjna834801_rows)}")

    prjeb55464_rows = [row for row in combined_rows if row.get("study") == "PRJEB55464"]
    filled_prjeb55464_disease_status = sum(1 for row in prjeb55464_rows if (row.get("disease_status") or "").strip())
    print(f"PRJEB55464 disease_status filled: {filled_prjeb55464_disease_status}/{len(prjeb55464_rows)}")

    for study_id in ["PRJNA531273", "PRJNA588035", "PRJNA808166", "PRJNA943232", "PRJNA762199", "PRJNA1083304"]:
        study_rows = [row for row in combined_rows if row.get("study") == study_id]
        filled_count = sum(1 for row in study_rows if (row.get("disease_status") or "").strip())
        print(f"{study_id} disease_status filled: {filled_count}/{len(study_rows)}")

    print("\nCase/control labels per study:")
    for study_id in sorted({row["study"] for row in combined_rows}):
        case_count = case_control_counts.get((study_id, "case"), 0)
        control_count = case_control_counts.get((study_id, "control"), 0)
        print(f"- {study_id}: cases={case_count}, controls={control_count}")

    print("\nMissingness per column:")
    for summary in summarize_missingness(combined_rows):
        print(
            f"- {summary['column']}: {summary['missing_count']}/{len(combined_rows)} missing "
            f"({summary['missing_fraction']:.0%})"
        )

    print(f"\nCombined metadata saved to: {OUTPUT_FILE.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
