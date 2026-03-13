#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Dict, List, Optional


KEYMAP = {
    "Number of input reads": ("input_reads_n", "int"),
    "Uniquely mapped reads number": ("uniquely_mapped_reads_n", "int"),
    "Uniquely mapped reads %": ("uniquely_mapped_reads_pct", "pct"),
    "Number of reads mapped to multiple loci": ("mapped_multiple_loci_n", "int"),
    "% of reads mapped to multiple loci": ("mapped_multiple_loci_pct", "pct"),
    "Number of reads mapped to too many loci": ("mapped_too_many_loci_n", "int"),
    "% of reads mapped to too many loci": ("mapped_too_many_loci_pct", "pct"),
    "Number of reads unmapped: too many mismatches": ("unmapped_too_many_mismatches_n", "int"),
    "% of reads unmapped: too many mismatches": ("unmapped_too_many_mismatches_pct", "pct"),
    "Number of reads unmapped: too short": ("unmapped_too_short_n", "int"),
    "% of reads unmapped: too short": ("unmapped_too_short_pct", "pct"),
    "Number of reads unmapped: other": ("unmapped_other_n", "int"),
    "% of reads unmapped: other": ("unmapped_other_pct", "pct"),
    "Number of chimeric reads": ("chimeric_reads_n", "int"),
    "% of chimeric reads": ("chimeric_reads_pct", "pct"),
}

PERCENT_CATEGORIES = [
    ("uniquely_mapped_reads_pct", "Uniquely mapped"),
    ("mapped_multiple_loci_pct", "Mapped to multiple loci"),
    ("mapped_too_many_loci_pct", "Mapped to too many loci"),
    ("unmapped_too_many_mismatches_pct", "Unmapped: too many mismatches"),
    ("unmapped_too_short_pct", "Unmapped: too short"),
    ("unmapped_other_pct", "Unmapped: other"),
    ("chimeric_reads_pct", "Chimeric"),
]

COUNT_CATEGORIES = [
    ("uniquely_mapped_reads_n", "Uniquely mapped"),
    ("mapped_multiple_loci_n", "Mapped to multiple loci"),
    ("mapped_too_many_loci_n", "Mapped to too many loci"),
    ("unmapped_too_many_mismatches_n", "Unmapped: too many mismatches"),
    ("unmapped_too_short_n", "Unmapped: too short"),
    ("unmapped_other_n", "Unmapped: other"),
    ("chimeric_reads_n", "Chimeric"),
]


def parse_value(raw: str, kind: str):
    raw = raw.strip()
    if kind == "int":
        return int(raw.replace(",", ""))
    if kind == "pct":
        return float(raw.replace("%", "").strip())
    return raw


def parse_log_file(log_path: Path) -> Dict[str, object]:
    result: Dict[str, object] = {
        "sample": log_path.parent.name,
        "log_final_path": str(log_path.resolve()),
    }

    for _, (outcol, kind) in KEYMAP.items():
        result[outcol] = None

    with log_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if "|" not in line:
                continue
            left, right = line.split("|", 1)
            key = left.strip()
            value = right.strip()

            if key in KEYMAP:
                outcol, kind = KEYMAP[key]
                result[outcol] = parse_value(value, kind)

    input_reads = result.get("input_reads_n")
    if input_reads is not None:
        categorized = 0
        for col, _ in COUNT_CATEGORIES:
            val = result.get(col)
            if val is not None:
                categorized += int(val)
        result["categorized_reads_n"] = categorized
        result["uncategorized_reads_n"] = int(input_reads) - categorized
    else:
        result["categorized_reads_n"] = None
        result["uncategorized_reads_n"] = None

    return result


def write_table(rows: List[Dict[str, object]], out_tsv: Path, out_csv: Path) -> None:
    fieldnames = [
        "sample",
        "log_final_path",
        "input_reads_n",
        "uniquely_mapped_reads_n",
        "uniquely_mapped_reads_pct",
        "mapped_multiple_loci_n",
        "mapped_multiple_loci_pct",
        "mapped_too_many_loci_n",
        "mapped_too_many_loci_pct",
        "unmapped_too_many_mismatches_n",
        "unmapped_too_many_mismatches_pct",
        "unmapped_too_short_n",
        "unmapped_too_short_pct",
        "unmapped_other_n",
        "unmapped_other_pct",
        "chimeric_reads_n",
        "chimeric_reads_pct",
        "categorized_reads_n",
        "uncategorized_reads_n",
    ]

    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    with out_tsv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    with out_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_long_tables(rows: List[Dict[str, object]], outdir: Path) -> None:
    pct_tsv = outdir / "star_mapping_qc_percent_long.tsv"
    cnt_tsv = outdir / "star_mapping_qc_counts_long.tsv"

    with pct_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "category", "percent"])
        for row in rows:
            for col, label in PERCENT_CATEGORIES:
                writer.writerow([row["sample"], label, row.get(col)])

    with cnt_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "category", "count"])
        for row in rows:
            for col, label in COUNT_CATEGORIES:
                writer.writerow([row["sample"], label, row.get(col)])


def main():
    parser = argparse.ArgumentParser(
        description="Parse STAR Log.final.out files into a mapping QC summary table."
    )
    parser.add_argument(
        "--star-dir",
        required=True,
        help="Directory containing STAR sample subdirectories with Log.final.out files.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory for parsed tables.",
    )
    args = parser.parse_args()

    star_dir = Path(args.star_dir)
    outdir = Path(args.outdir)

    if not star_dir.is_dir():
        raise SystemExit(f"ERROR: STAR directory not found: {star_dir}")

    log_files = sorted(star_dir.glob("*/Log.final.out"))
    if not log_files:
        raise SystemExit(f"ERROR: no Log.final.out files found under {star_dir}")

    rows = [parse_log_file(p) for p in log_files]
    rows = sorted(rows, key=lambda x: x["sample"])

    out_tsv = outdir / "star_mapping_qc_summary.tsv"
    out_csv = outdir / "star_mapping_qc_summary.csv"

    write_table(rows, out_tsv, out_csv)
    write_long_tables(rows, outdir)

    print(f">>> Parsed {len(rows)} STAR logs")
    print(f">>> Summary TSV: {out_tsv}")
    print(f">>> Summary CSV: {out_csv}")


if __name__ == "__main__":
    main()
