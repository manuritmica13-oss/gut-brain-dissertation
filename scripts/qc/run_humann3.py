#!/usr/bin/env python3
"""Run KneadData and HUMAnN3 on shotgun metagenomics FASTQ samples."""

import argparse
import gzip
import shutil
import subprocess
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run KneadData followed by HUMAnN3 on shotgun metagenomics FASTQ files."
    )
    parser.add_argument(
        "--input-dir",
        required=True,
        help="Directory containing shotgun FASTQ files for each sample.",
    )
    parser.add_argument(
        "--output-dir",
        default="data/processed/humann3",
        help="Directory to write KneadData and HUMAnN3 outputs.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Number of threads for KneadData and HUMAnN3.",
    )
    return parser.parse_args()


def find_fastq_pairs(input_dir):
    """Collect paired FASTQ paths from the input directory."""
    fastq_dir = Path(input_dir)
    fastq_files = sorted(fastq_dir.glob("*.fastq*"))
    if not fastq_files:
        raise FileNotFoundError(f"No FASTQ files found in {input_dir}")

    samples = {}
    for path in fastq_files:
        name = path.name
        if "_R1" in name or "_1" in name:
            sample_id = name.split("_R1")[0].split("_1")[0]
            samples.setdefault(sample_id, {})["forward"] = str(path.resolve())
        elif "_R2" in name or "_2" in name:
            sample_id = name.split("_R2")[0].split("_2")[0]
            samples.setdefault(sample_id, {})["reverse"] = str(path.resolve())

    paired_samples = []
    for sample_id, paths in samples.items():
        if "forward" in paths and "reverse" in paths:
            paired_samples.append((sample_id, paths["forward"], paths["reverse"]))

    if not paired_samples:
        raise ValueError("No paired FASTQ samples could be matched. Use *_R1/_R2 or *_1/*_2 naming.")

    return paired_samples


def run_command(command, description):
    """Run a shell command and raise if it fails."""
    print(f"\n[humann3] {description}")
    print(" ".join(command))
    subprocess.run(command, check=True)


def run_kneaddata(sample_id, forward, reverse, output_dir, threads):
    """Run KneadData on a paired sample and return the cleaned output file path."""
    kneaddata_out = output_dir / "kneaddata" / sample_id
    kneaddata_out.mkdir(parents=True, exist_ok=True)

    # Run KneadData to remove human contamination and produce cleaned paired reads.
    run_command(
        [
            "kneaddata",
            "--input",
            forward,
            "--input",
            reverse,
            "--output",
            str(kneaddata_out),
            "--threads",
            str(threads),
            "--verbose",
        ],
        f"Running KneadData for sample {sample_id}",
    )

    # KneadData outputs paired cleaned files with sample prefix in the sample directory.
    return kneaddata_out


def run_humann3(sample_id, kneaddata_dir, output_dir, threads):
    """Run HUMAnN3 on the cleaned sample reads."""
    humann_out = output_dir / "humann3" / sample_id
    humann_out.mkdir(parents=True, exist_ok=True)

    paired1 = next(kneaddata_dir.glob("*paired_1.fastq*"), None)
    paired2 = next(kneaddata_dir.glob("*paired_2.fastq*"), None)

    if paired1 and paired2:
        is_gz = paired1.suffix == ".gz" or paired2.suffix == ".gz"
        merged_name = f"{sample_id}_kneaddata_combined.fastq"
        if is_gz:
            merged_name += ".gz"
        merged_input = humann_out / merged_name

        # Concatenate the cleaned paired reads into one file for HUMAnN3.
        if is_gz:
            with gzip.open(merged_input, "wb") as out_handle:
                for part in (paired1, paired2):
                    with gzip.open(part, "rb") as in_handle:
                        shutil.copyfileobj(in_handle, out_handle)
        else:
            with open(merged_input, "wb") as out_handle:
                for part in (paired1, paired2):
                    with open(part, "rb") as in_handle:
                        shutil.copyfileobj(in_handle, out_handle)

        input_path = str(merged_input)
    else:
        single = next(kneaddata_dir.glob("*.fastq*"), None)
        if not single:
            raise FileNotFoundError(f"No cleaned FASTQ output found for sample {sample_id}")
        input_path = str(single)

    # Run HUMAnN3 to generate gene family and pathway abundance tables.
    run_command(
        [
            "humann",
            "--input",
            input_path,
            "--output",
            str(humann_out),
            "--threads",
            str(threads),
        ],
        f"Running HUMAnN3 for sample {sample_id}",
    )

    return humann_out


def merge_humann_outputs(humann_dir, output_dir):
    """Merge HUMAnN3 sample tables into a set of combined tables."""
    merged_dir = output_dir / "merged"
    merged_dir.mkdir(parents=True, exist_ok=True)

    # Merge gene families and pathway abundance outputs across all samples.
    for file_name in ["genefamilies.tsv", "pathabundance.tsv"]:
        output_file = merged_dir / f"humann_{file_name.replace('.tsv','')}_merged.tsv"
        command = [
            "humann_join_tables",
            "--input",
            str(humann_dir),
            "--output",
            str(output_file),
            "--file_name",
            file_name,
        ]
        run_command(command, f"Merging HUMAnN3 {file_name} tables")

    print(f"\nMerged HUMAnN3 tables written to: {merged_dir}")
    return merged_dir


def main():
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    paired_samples = find_fastq_pairs(args.input_dir)
    humann_output_dir = output_dir / "humann3"

    for sample_id, forward, reverse in paired_samples:
        kneaddata_dir = run_kneaddata(sample_id, forward, reverse, output_dir, args.threads)
        run_humann3(sample_id, kneaddata_dir, output_dir, args.threads)

    merge_humann_outputs(humann_output_dir, output_dir)

    print(f"\nShotgun metagenomics processing complete. Final results in {output_dir}")


if __name__ == "__main__":
    main()
