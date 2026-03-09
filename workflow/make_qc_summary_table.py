#!/usr/bin/env python3
import argparse
import os
import re
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description="Build per-sample QC summary from seqkit raw/trimmed tables")
    p.add_argument("--raw", required=True, help="Path to seqkit_stats_raw.tsv")
    p.add_argument("--trimmed", required=True, help="Path to seqkit_stats_trimmed.tsv")
    p.add_argument("--out", required=True, help="Output TSV")
    return p.parse_args()


def strip_extensions(fn: str) -> str:
    for ext in [".fastq.gz", ".fq.gz", ".fastq", ".fq"]:
        if fn.endswith(ext):
            return fn[: -len(ext)]
    return fn


def parse_file_info(path: str):
    bn = os.path.basename(path)
    base = strip_extensions(bn)

    # remove ".trimmed" suffix if present
    is_trimmed = False
    if base.endswith(".trimmed"):
        base = base[:-8]
        is_trimmed = True

    # Illumina-style paired-end: sample_R1_001 / sample_R2_001
    m = re.match(r"(.+)_R([12])_001$", base)
    if m:
        sample = m.group(1)
        read = f"R{m.group(2)}"
        return sample, read

    # Common paired-end: sample_R1 / sample_R2
    m = re.match(r"(.+)_R([12])$", base)
    if m:
        sample = m.group(1)
        read = f"R{m.group(2)}"
        return sample, read

    # Alternative paired-end: sample_1 / sample_2
    m = re.match(r"(.+)_([12])$", base)
    if m:
        sample = m.group(1)
        read = f"R{m.group(2)}"
        return sample, read

    # otherwise treat as single-end
    return base, "SE"


def load_and_annotate(path: str, prefix: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        raise ValueError(f"Empty file: {path}")

    parsed = df["file"].apply(parse_file_info)
    df["sample"] = parsed.apply(lambda x: x[0])
    df["read"] = parsed.apply(lambda x: x[1])

    keep = [
        "sample", "read", "num_seqs", "sum_len", "min_len", "avg_len", "max_len",
        "Q20(%)", "Q30(%)", "AvgQual", "GC(%)", "sum_n"
    ]
    df = df[keep].copy()

    rename = {
        "num_seqs": f"{prefix}_num_seqs",
        "sum_len": f"{prefix}_sum_len",
        "min_len": f"{prefix}_min_len",
        "avg_len": f"{prefix}_avg_len",
        "max_len": f"{prefix}_max_len",
        "Q20(%)": f"{prefix}_Q20",
        "Q30(%)": f"{prefix}_Q30",
        "AvgQual": f"{prefix}_AvgQual",
        "GC(%)": f"{prefix}_GC",
        "sum_n": f"{prefix}_sum_n",
    }
    return df.rename(columns=rename)


def pct(numer, denom):
    return (numer / denom * 100).where(denom != 0, 0.0)


def main():
    args = parse_args()

    raw = load_and_annotate(args.raw, "raw")
    trimmed = load_and_annotate(args.trimmed, "trim")

    merged = raw.merge(trimmed, on=["sample", "read"], how="inner")

    if merged.empty:
        raise ValueError(
            "No matching sample/read pairs were found between raw and trimmed tables. "
            "Check filename parsing for raw vs trimmed FASTQs."
        )

    merged["reads_removed"] = merged["raw_num_seqs"] - merged["trim_num_seqs"]
    merged["reads_removed_pct"] = pct(merged["reads_removed"], merged["raw_num_seqs"])

    merged["bases_removed"] = merged["raw_sum_len"] - merged["trim_sum_len"]
    merged["bases_removed_pct"] = pct(merged["bases_removed"], merged["raw_sum_len"])

    merged["Q20_delta"] = merged["trim_Q20"] - merged["raw_Q20"]
    merged["Q30_delta"] = merged["trim_Q30"] - merged["raw_Q30"]
    merged["AvgQual_delta"] = merged["trim_AvgQual"] - merged["raw_AvgQual"]
    merged["GC_delta"] = merged["trim_GC"] - merged["raw_GC"]

    col_order = [
        "sample", "read",
        "raw_num_seqs", "trim_num_seqs", "reads_removed", "reads_removed_pct",
        "raw_sum_len", "trim_sum_len", "bases_removed", "bases_removed_pct",
        "raw_min_len", "trim_min_len",
        "raw_avg_len", "trim_avg_len",
        "raw_max_len", "trim_max_len",
        "raw_Q20", "trim_Q20", "Q20_delta",
        "raw_Q30", "trim_Q30", "Q30_delta",
        "raw_AvgQual", "trim_AvgQual", "AvgQual_delta",
        "raw_GC", "trim_GC", "GC_delta",
        "raw_sum_n", "trim_sum_n",
    ]
    merged = merged[col_order].copy()

    total_cols = {
        "raw_num_seqs", "trim_num_seqs", "reads_removed",
        "raw_sum_len", "trim_sum_len", "bases_removed",
        "raw_sum_n", "trim_sum_n"
    }

    summary = {}
    summary["sample"] = "Average/Total"
    summary["read"] = "ALL"

    for c in merged.columns:
        if c in {"sample", "read"}:
            continue
        if c in total_cols:
            summary[c] = merged[c].sum()
        else:
            summary[c] = merged[c].mean()

    out = pd.concat([merged, pd.DataFrame([summary])], ignore_index=True)

    int_like_cols = [
        "raw_num_seqs", "trim_num_seqs", "reads_removed",
        "raw_sum_len", "trim_sum_len", "bases_removed",
        "raw_min_len", "trim_min_len",
        "raw_max_len", "trim_max_len",
        "raw_sum_n", "trim_sum_n"
    ]
    for c in int_like_cols:
        if c in out.columns:
            out[c] = out[c].round(0).astype("Int64")

    float_cols_2 = [
        "reads_removed_pct", "bases_removed_pct",
        "raw_avg_len", "trim_avg_len",
        "raw_Q20", "trim_Q20", "Q20_delta",
        "raw_Q30", "trim_Q30", "Q30_delta",
        "raw_AvgQual", "trim_AvgQual", "AvgQual_delta",
        "raw_GC", "trim_GC", "GC_delta",
    ]
    for c in float_cols_2:
        if c in out.columns:
            out[c] = out[c].round(2)

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
