#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   bash workflow/run_rseqc.sh --bam <sample.bam> --bed <genes.bed12> --outdir results/rseqc/<sample>
#
# Notes:
# - BAM must be coordinate-sorted and indexed (.bai present or creatable)
# - BED must be BED12 (RSeQC format), not plain BED6

BAM=""
BED=""
OUTDIR=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam) BAM="$2"; shift 2 ;;
    --bed) BED="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    *) echo "Unknown arg: $1" >&2; exit 2 ;;
  esac
done

if [[ -z "${BAM}" || -z "${BED}" || -z "${OUTDIR}" ]]; then
  echo "ERROR: missing --bam/--bed/--outdir" >&2
  exit 2
fi

mkdir -p "${OUTDIR}"

# Ensure BAM index exists
if [[ ! -f "${BAM}.bai" && ! -f "${BAM%.bam}.bai" ]]; then
  samtools index "${BAM}"
fi

echo ">>> RSeQC: infer_experiment.py"
infer_experiment.py -i "${BAM}" -r "${BED}" > "${OUTDIR}/infer_experiment.txt" 2> "${OUTDIR}/infer_experiment.log" || true

echo ">>> RSeQC: geneBody_coverage.py"
geneBody_coverage.py -i "${BAM}" -r "${BED}" -o "${OUTDIR}/geneBody" \
  > "${OUTDIR}/geneBody_coverage.log" 2>&1 || true

echo ">>> Done. Outputs in: ${OUTDIR}"
