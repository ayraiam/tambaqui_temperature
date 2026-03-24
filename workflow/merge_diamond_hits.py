#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd

DIAMOND_COLS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
    "qlen", "slen", "qcovhsp"
]

def norm_text(x):
    if pd.isna(x):
        return None
    s = str(x).strip()
    if s == "" or s.lower() == "nan":
        return None
    return s

def acc_no_version(x):
    x = norm_text(x)
    if x is None:
        return None
    return x.split(".")[0]

def sanitize_columns(df: pd.DataFrame, prefix: str) -> pd.DataFrame:
    out = df.copy()
    out.columns = [
        prefix + c.strip().replace(" ", "_").replace("/", "_")
        for c in out.columns
    ]
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-annotation-tsv", required=True)
    ap.add_argument("--diamond-tsv", required=True)
    ap.add_argument("--target-metadata-tsv", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    source_ann = pd.read_csv(args.source_annotation_tsv, sep="\t")
    diamond = pd.read_csv(args.diamond_tsv, sep="\t", header=None, names=DIAMOND_COLS)
    target = pd.read_csv(args.target_metadata_tsv, sep="\t")

    if "Protein accession" not in target.columns:
        raise ValueError("Target metadata missing column: Protein accession")
    if "Gene ID" not in target.columns:
        raise ValueError("Target metadata missing column: Gene ID")

    diamond["qseqid"] = diamond["qseqid"].map(norm_text)
    diamond["sseqid"] = diamond["sseqid"].map(norm_text)
    diamond["qseqid_noversion"] = diamond["qseqid"].map(acc_no_version)
    diamond["sseqid_noversion"] = diamond["sseqid"].map(acc_no_version)

    diamond = diamond.sort_values(["qseqid", "bitscore", "evalue"], ascending=[True, False, True])
    best = diamond.drop_duplicates(subset=["qseqid"], keep="first").copy()

    target["target_protein_accession"] = target["Protein accession"].map(norm_text)
    target["target_protein_accession_noversion"] = target["target_protein_accession"].map(acc_no_version)
    target_small = sanitize_columns(target, "target_")

    best = best.merge(
        target_small,
        left_on="sseqid_noversion",
        right_on="target_target_protein_accession_noversion",
        how="left"
    )

    source_ann["source_protein_accession_noversion"] = source_ann["source_source_protein_accession"].map(acc_no_version)

    merged = source_ann.merge(
        best,
        left_on="source_protein_accession_noversion",
        right_on="qseqid_noversion",
        how="left"
    )

    out_tsv = outdir / "deseq_annotated_with_danio_hits.tsv"
    merged.to_csv(out_tsv, sep="\t", index=False)

    summary = pd.DataFrame({
        "metric": [
            "rows_source_annotation",
            "diamond_rows",
            "best_hits",
            "rows_with_target_hit",
            "rows_with_target_gene_id"
        ],
        "value": [
            len(source_ann),
            len(diamond),
            len(best),
            merged["sseqid"].notna().sum(),
            merged["target_Gene_ID"].notna().sum() if "target_Gene_ID" in merged.columns else 0
        ]
    })
    summary.to_csv(outdir / "merge_diamond_summary.tsv", sep="\t", index=False)

    print(f">>> Wrote merged annotation table: {out_tsv}")

if __name__ == "__main__":
    main()
