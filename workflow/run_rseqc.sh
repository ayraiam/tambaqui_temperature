#!/usr/bin/env bash
set -euo pipefail

# ==========================================================
# Script: workflow/run_rseqc.sh
# Purpose: Run RSeQC checks on one BAM or all BAMs in a directory
#
# Supports:
#   - strandedness inference with infer_experiment.py
#   - optional geneBody_coverage.py
#   - optional GTF -> BED12 conversion
#
# Examples:
#   bash workflow/run_rseqc.sh \
#     --bam-dir results/star \
#     --gtf reference/genes.gtf \
#     --make-bed12 \
#     --bed-out reference/genes.bed12 \
#     --outroot results/rseqc \
#     --infer-only
#
#   bash workflow/run_rseqc.sh \
#     --bam results/star/sample1/Aligned.sortedByCoord.out.bam \
#     --bed reference/genes.bed12 \
#     --outroot results/rseqc
# ==========================================================

BAM=""
BAM_DIR=""
BED=""
GTF=""
MAKE_BED12=0
BED_OUT=""
OUTROOT="results/rseqc"
INFER_ONLY=0
THREADS="${THREADS:-1}"

usage() {
  cat <<EOF
Usage:
  bash workflow/run_rseqc.sh [options]

Input BAM:
  --bam PATH            Single coordinate-sorted BAM
  --bam-dir DIR         Directory containing coordinate-sorted BAMs

Annotation:
  --bed PATH            BED12 file for RSeQC
  --gtf PATH            GTF file (used if BED12 needs to be made)
  --make-bed12          Convert GTF -> BED12
  --bed-out PATH        Output BED12 path (default: <gtf_basename>.bed12 next to GTF)

Output:
  --outroot DIR         Output root directory (default: results/rseqc)

Mode:
  --infer-only          Run only infer_experiment.py
  -h, --help            Show this help
EOF
  exit 0
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam) BAM="$2"; shift 2 ;;
    --bam-dir) BAM_DIR="$2"; shift 2 ;;
    --bed) BED="$2"; shift 2 ;;
    --gtf) GTF="$2"; shift 2 ;;
    --make-bed12) MAKE_BED12=1; shift 1 ;;
    --bed-out) BED_OUT="$2"; shift 2 ;;
    --outroot) OUTROOT="$2"; shift 2 ;;
    --infer-only) INFER_ONLY=1; shift 1 ;;
    -h|--help) usage ;;
    *) echo "Unknown arg: $1" >&2; exit 2 ;;
  esac
done

# -------------------------
# Validate BAM input
# -------------------------
if [[ -z "${BAM}" && -z "${BAM_DIR}" ]]; then
  echo "ERROR: provide --bam or --bam-dir" >&2
  exit 2
fi

if [[ -n "${BAM}" && -n "${BAM_DIR}" ]]; then
  echo "ERROR: use only one of --bam or --bam-dir" >&2
  exit 2
fi

# -------------------------
# Resolve/build BED12
# -------------------------
if [[ "${MAKE_BED12}" -eq 1 ]]; then
  if [[ -z "${GTF}" ]]; then
    echo "ERROR: --make-bed12 requires --gtf" >&2
    exit 2
  fi
  [[ -f "${GTF}" ]] || { echo "ERROR: GTF not found: ${GTF}" >&2; exit 2; }

  if [[ -z "${BED_OUT}" ]]; then
    BED_OUT="${GTF%.*}.bed12"
  fi

  mkdir -p "$(dirname "${BED_OUT}")"
  TMP_GP="${BED_OUT%.bed12}.genepred"

  command -v gtfToGenePred >/dev/null 2>&1 || {
    echo "ERROR: gtfToGenePred not found in PATH" >&2
    exit 2
  }
  command -v genePredToBed >/dev/null 2>&1 || {
    echo "ERROR: genePredToBed not found in PATH" >&2
    exit 2
  }

  echo ">>> Converting GTF -> GenePred"
  gtfToGenePred "${GTF}" "${TMP_GP}"

  echo ">>> Converting GenePred -> BED12"
  genePredToBed "${TMP_GP}" "${BED_OUT}"

  rm -f "${TMP_GP}"
  BED="${BED_OUT}"
fi

if [[ -z "${BED}" ]]; then
  echo "ERROR: provide --bed, or use --gtf --make-bed12" >&2
  exit 2
fi

[[ -f "${BED}" ]] || { echo "ERROR: BED12 not found: ${BED}" >&2; exit 2; }

mkdir -p "${OUTROOT}"

# -------------------------
# Gather BAMs
# -------------------------
declare -a BAMS=()

if [[ -n "${BAM}" ]]; then
  [[ -f "${BAM}" ]] || { echo "ERROR: BAM not found: ${BAM}" >&2; exit 2; }
  BAMS+=("${BAM}")
else
  while IFS= read -r -d '' f; do
    BAMS+=("$f")
  done < <(find "${BAM_DIR}" -type f -name "Aligned.sortedByCoord.out.bam" -print0 | sort -z)

  if [[ "${#BAMS[@]}" -eq 0 ]]; then
    echo "ERROR: no BAM files found in ${BAM_DIR}" >&2
    exit 2
  fi
fi

echo ">>> Found ${#BAMS[@]} BAM file(s)"

# -------------------------
# Run RSeQC
# -------------------------
for BAM_FILE in "${BAMS[@]}"; do
  BAM_BASE="$(basename "${BAM_FILE}")"
  SAMPLE="$(basename "$(dirname "${BAM_FILE}")")"

  # If BAMs are all in one folder, use BAM filename stem
  if [[ "${SAMPLE}" == "." || "${SAMPLE}" == "/" || -z "${SAMPLE}" ]]; then
    SAMPLE="${BAM_BASE%.bam}"
  fi
  if [[ "${SAMPLE}" == "star" || "${SAMPLE}" == "bam" || "${SAMPLE}" == "results" ]]; then
    SAMPLE="${BAM_BASE%.bam}"
  fi

  OUTDIR="${OUTROOT}/${SAMPLE}"
  mkdir -p "${OUTDIR}"

  echo "======================================"
  echo ">>> Processing: ${BAM_FILE}"
  echo ">>> Sample:     ${SAMPLE}"
  echo ">>> Output:     ${OUTDIR}"

  # Ensure BAM index exists
  if [[ ! -f "${BAM_FILE}.bai" && ! -f "${BAM_FILE%.bam}.bai" ]]; then
    echo ">>> Indexing BAM with samtools"
    samtools index "${BAM_FILE}"
  fi

  echo ">>> Running infer_experiment.py"
  infer_experiment.py -i "${BAM_FILE}" -r "${BED}" \
    > "${OUTDIR}/infer_experiment.txt" \
    2> "${OUTDIR}/infer_experiment.log" || true

  if [[ "${INFER_ONLY}" -eq 0 ]]; then
    echo ">>> Running geneBody_coverage.py"
    geneBody_coverage.py -i "${BAM_FILE}" -r "${BED}" -o "${OUTDIR}/geneBody" \
      > "${OUTDIR}/geneBody_coverage.log" 2>&1 || true
  fi
done

echo ">>> RSeQC finished."
echo ">>> Outputs written under: ${OUTROOT}"
