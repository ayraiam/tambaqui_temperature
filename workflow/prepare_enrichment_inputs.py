#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd
from Bio import SeqIO

def norm_text(x):
    if pd.isna(x):
        return None
    s = str(x).strip()
    if s == "" or s.lower() == "nan":
        return None
    if s.endswith(".0"):
        s = s[:-2]
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

def build_fasta_lookup(faa_path: Path):
    exact = {}
    noversion = {}
    for rec in SeqIO.parse(str(faa_path), "fasta"):
        rid = rec.id.strip()
        exact[rid] = str(rec.seq)
        noversion[acc_no_version(rid)] = str(rec.seq)
    return exact, noversion

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--deseq-tsv", required=True)
    ap.add_argument("--source-metadata-tsv", required=True)
    ap.add_argument("--source-protein-faa", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--deseq-join-col", default="Geneid")
    ap.add_argument("--source-join-col", default="Gene ID")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    deseq = pd.read_csv(args.deseq_tsv, sep="\t")
    source = pd.read_csv(args.source_metadata_tsv, sep="\t")

    if args.deseq_join_col not in deseq.columns:
        raise ValueError(f"DESeq table missing join col: {args.deseq_join_col}")
    if args.source_join_col not in source.columns:
        raise ValueError(f"Source metadata missing join col: {args.source_join_col}")
    if "Protein accession" not in source.columns:
        raise ValueError("Source metadata missing column: Protein accession")

    deseq["_join_key"] = deseq[args.deseq_join_col].map(norm_text)
    source["_join_key"] = source[args.source_join_col].map(norm_text)

    source_small = source.copy()
    source_small["source_protein_accession"] = source_small["Protein accession"].map(norm_text)
    source_small["source_protein_accession_noversion"] = source_small["source_protein_accession"].map(acc_no_version)
    source_small = sanitize_columns(source_small, "source_")

    merged = deseq.merge(
        source_small,
        left_on="_join_key",
        right_on="source__join_key",
        how="left"
    )

    exact_faa, noversion_faa = build_fasta_lookup(Path(args.source_protein_faa))

    def get_seq(acc):
        acc = norm_text(acc)
        if acc is None:
            return None
        if acc in exact_faa:
            return exact_faa[acc]
        acc_nv = acc_no_version(acc)
        return noversion_faa.get(acc_nv)

    merged["source_aa_sequence"] = merged["source_source_protein_accession"].map(get_seq)

    ann_tsv = outdir / "source_annotation_from_deseq.tsv"
    merged.to_csv(ann_tsv, sep="\t", index=False)

    query_rows = merged.loc[
        merged["source_source_protein_accession"].notna() &
        merged["source_aa_sequence"].notna(),
        ["source_source_protein_accession", "source_aa_sequence"]
    ].drop_duplicates()

    query_faa = outdir / "query_proteins.faa"
    with open(query_faa, "w") as fh:
        for _, row in query_rows.iterrows():
            fh.write(f">{row['source_source_protein_accession']}\n{row['source_aa_sequence']}\n")

    summary = pd.DataFrame({
        "metric": [
            "input_rows",
            "rows_with_source_metadata",
            "rows_with_protein_accession",
            "rows_with_sequence",
            "unique_query_proteins"
        ],
        "value": [
            len(merged),
            merged["source_source_Protein_accession"].notna().sum() if "source_source_Protein_accession" in merged.columns else merged["source_source_protein_accession"].notna().sum(),
            merged["source_source_protein_accession"].notna().sum(),
            merged["source_aa_sequence"].notna().sum(),
            len(query_rows)
        ]
    })
    summary.to_csv(outdir / "prepare_enrichment_summary.tsv", sep="\t", index=False)

    print(f">>> Wrote source annotation table: {ann_tsv}")
    print(f">>> Wrote query FASTA: {query_faa}")

if __name__ == "__main__":
    main()
