#!/usr/bin/env python3
"""Run a QIIME2 16S DADA2 pipeline for paired-end gut microbiome data."""

import argparse
import csv
import os
import subprocess
import sys
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Process paired-end 16S FASTQ reads with QIIME2 DADA2 and export results."
    )
    parser.add_argument(
        "--input-dir",
        required=True,
        help="Directory containing paired-end FASTQ files for each sample.",
    )
    parser.add_argument(
        "--metadata-file",
        required=True,
        help="Sample metadata TSV file for downstream QIIME2 analyses.",
    )
    parser.add_argument(
        "--output-dir",
        default="data/processed/qiime2",
        help="Directory to write QIIME2 outputs (default: data/processed/qiime2).",
    )
    return parser.parse_args()


def find_paired_fastqs(input_dir):
    """Scan the input directory and return a list of paired FASTQ samples."""
    fastq_dir = Path(input_dir)
    forward_reads = sorted(fastq_dir.glob("*R1*.fastq*"))
    reverse_reads = sorted(fastq_dir.glob("*R2*.fastq*"))

    if not forward_reads or not reverse_reads:
        raise FileNotFoundError("No paired-end FASTQ files found in input directory.")

    pairs = {}
    for fwd in forward_reads:
        sample_id = fwd.name.split("_")[0]
        pairs.setdefault(sample_id, {})["forward"] = str(fwd.resolve())
    for rev in reverse_reads:
        sample_id = rev.name.split("_")[0]
        pairs.setdefault(sample_id, {})["reverse"] = str(rev.resolve())

    paired_samples = [
        (sample_id, data["forward"], data["reverse"])
        for sample_id, data in pairs.items()
        if "forward" in data and "reverse" in data
    ]

    if not paired_samples:
        raise ValueError("No matching paired-end FASTQ samples found.")

    return paired_samples


def write_manifest(pairs, manifest_path):
    """Write a QIIME2 paired-end FASTQ manifest file."""
    with open(manifest_path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample-id", "forward-absolute-filepath", "reverse-absolute-filepath"])
        for sample_id, forward_path, reverse_path in pairs:
            writer.writerow([sample_id, forward_path, reverse_path])


def run_command(cmd, description):
    """Run a shell command and raise an error if it fails."""
    print(f"\n[qiime2] {description}")
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)


def main():
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    metadata_file = Path(args.metadata_file)
    if not metadata_file.exists():
        raise FileNotFoundError(f"Metadata file not found: {metadata_file}")

    paired_samples = find_paired_fastqs(args.input_dir)
    manifest_path = output_dir / "manifest.csv"
    write_manifest(paired_samples, manifest_path)

    demux_qza = output_dir / "demux.qza"
    demux_qzv = output_dir / "demux.qzv"
    table_qza = output_dir / "table.qza"
    rep_seqs_qza = output_dir / "rep-seqs.qza"
    stats_qza = output_dir / "denoising-stats.qza"
    filtered_table_qza = output_dir / "filtered-table.qza"
    export_table_dir = output_dir / "exported-feature-table"
    export_repseqs_dir = output_dir / "exported-rep-seqs"

    # Import paired-end FASTQ files into a Qiime2 artifact.
    run_command(
        [
            "qiime",
            "tools",
            "import",
            "--type",
            "SampleData[PairedEndSequencesWithQuality]",
            "--input-path",
            str(manifest_path),
            "--output-path",
            str(demux_qza),
            "--input-format",
            "PairedEndFastqManifestPhred33V2",
        ],
        "Importing paired-end FASTQ reads into QIIME2",
    )

    # Generate a demultiplexed summary visualization to inspect read quality.
    run_command(
        [
            "qiime",
            "demux",
            "summarize",
            "--i-data",
            str(demux_qza),
            "--o-visualization",
            str(demux_qzv),
        ],
        "Creating demux summary visualization",
    )

    # Denoise paired-end reads with DADA2 using default truncation settings.
    run_command(
        [
            "qiime",
            "dada2",
            "denoise-paired",
            "--i-demultiplexed-seqs",
            str(demux_qza),
            "--p-trunc-len-f",
            "0",
            "--p-trunc-len-r",
            "0",
            "--o-table",
            str(table_qza),
            "--o-representative-sequences",
            str(rep_seqs_qza),
            "--o-denoising-stats",
            str(stats_qza),
        ],
        "Running DADA2 denoising for paired-end reads",
    )

    # Filter out samples with fewer than 1000 total reads.
    run_command(
        [
            "qiime",
            "feature-table",
            "filter-samples",
            "--i-table",
            str(table_qza),
            "--p-min-frequency",
            "1000",
            "--o-filtered-table",
            str(filtered_table_qza),
        ],
        "Filtering low-depth samples with minimum 1000 reads",
    )

    # Export the filtered feature table and representative sequences for downstream use.
    export_table_dir.mkdir(exist_ok=True)
    export_repseqs_dir.mkdir(exist_ok=True)

    run_command(
        [
            "qiime",
            "tools",
            "export",
            "--input-path",
            str(filtered_table_qza),
            "--output-path",
            str(export_table_dir),
        ],
        "Exporting filtered feature table",
    )

    run_command(
        [
            "qiime",
            "tools",
            "export",
            "--input-path",
            str(rep_seqs_qza),
            "--output-path",
            str(export_repseqs_dir),
        ],
        "Exporting representative sequences",
    )

    print(f"\nQIIME2 16S processing complete. Outputs written to: {output_dir}")
    print(f"- Imported demux artifact: {demux_qza}")
    print(f"- Filtered feature table: {filtered_table_qza}")
    print(f"- Representative sequences: {rep_seqs_qza}")
    print(f"- Exported tables: {export_table_dir}")
    print(f"- Exported rep seqs: {export_repseqs_dir}")


if __name__ == "__main__":
    main()
