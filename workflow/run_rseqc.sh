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
#   - auto-check/install of missing tools in the active conda env
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

log() {
  echo ">>> $*"
}

die() {
  echo "ERROR: $*" >&2
  exit 2
}

have_cmd() {
  command -v "$1" >/dev/null 2>&1
}

ensure_conda_available() {
  if ! have_cmd conda; then
    die "conda not found in PATH. Activate the intended conda environment first."
  fi
}

ensure_active_conda_env() {
  if [[ -z "${CONDA_PREFIX:-}" || -z "${CONDA_DEFAULT_ENV:-}" ]]; then
    die "No active conda environment detected. Activate the target env first."
  fi
}

install_missing_tools() {
  local need_rseqc=0
  local need_samtools=0
  local need_ucsc=0

  if ! have_cmd infer_experiment.py; then
    need_rseqc=1
  fi

  if [[ "${INFER_ONLY}" -eq 0 ]] && ! have_cmd geneBody_coverage.py; then
    need_rseqc=1
  fi

  if ! have_cmd samtools; then
    need_samtools=1
  fi

  if [[ "${MAKE_BED12}" -eq 1 ]]; then
    if ! have_cmd gtfToGenePred || ! have_cmd genePredToBed; then
      need_ucsc=1
    fi
  fi

  if [[ "${need_rseqc}" -eq 0 && "${need_samtools}" -eq 0 && "${need_ucsc}" -eq 0 ]]; then
    log "All required tools already exist in the active environment: ${CONDA_DEFAULT_ENV}"
    return 0
  fi

  ensure_conda_available
  ensure_active_conda_env

  local -a pkgs=()
  [[ "${need_rseqc}" -eq 1 ]] && pkgs+=(rseqc)
  [[ "${need_samtools}" -eq 1 ]] && pkgs+=(samtools)
  if [[ "${need_ucsc}" -eq 1 ]]; then
    pkgs+=(ucsc-gtftogenepred ucsc-genepredtobed)
  fi

  log "Missing tools detected in env: ${CONDA_DEFAULT_ENV}"
  log "Installing required package(s): ${pkgs[*]}"

  conda install -y -c conda-forge -c bioconda "${pkgs[@]}"

  # Re-check after install
  have_cmd infer_experiment.py || die "infer_experiment.py still not found after installation"
  if [[ "${INFER_ONLY}" -eq 0 ]]; then
    have_cmd geneBody_coverage.py || die "geneBody_coverage.py still not found after installation"
  fi
  have_cmd samtools || die "samtools still not found after installation"
  if [[ "${MAKE_BED12}" -eq 1 ]]; then
    have_cmd gtfToGenePred || die "gtfToGenePred still not found after installation"
    have_cmd genePredToBed || die "genePredToBed still not found after installation"
  fi

  log "Tool installation/check completed successfully."
}

# -------------------------
# Validate BAM input
# -------------------------
if [[ -z "${BAM}" && -z "${BAM_DIR}" ]]; then
  die "provide --bam or --bam-dir"
fi

if [[ -n "${BAM}" && -n "${BAM_DIR}" ]]; then
  die "use only one of --bam or --bam-dir"
fi

# -------------------------
# Ensure tools are available
# -------------------------
install_missing_tools

# -------------------------
# Resolve/build BED12
# -------------------------
if [[ "${MAKE_BED12}" -eq 1 ]]; then
  if [[ -z "${GTF}" ]]; then
    die "--make-bed12 requires --gtf"
  fi
  [[ -f "${GTF}" ]] || die "GTF not found: ${GTF}"

  if [[ -z "${BED_OUT}" ]]; then
    BED_OUT="${GTF%.*}.bed12"
  fi

  mkdir -p "$(dirname "${BED_OUT}")"
  TMP_GP="${BED_OUT%.bed12}.genepred"

  log "Converting GTF -> GenePred"
  gtfToGenePred "${GTF}" "${TMP_GP}"

  log "Converting GenePred -> BED12"
  genePredToBed "${TMP_GP}" "${BED_OUT}"

  rm -f "${TMP_GP}"
  BED="${BED_OUT}"
fi

if [[ -z "${BED}" ]]; then
  die "provide --bed, or use --gtf --make-bed12"
fi

[[ -f "${BED}" ]] || die "BED12 not found: ${BED}"

mkdir -p "${OUTROOT}"

# -------------------------
# Gather BAMs
# -------------------------
declare -a BAMS=()

if [[ -n "${BAM}" ]]; then
  [[ -f "${BAM}" ]] || die "BAM not found: ${BAM}"
  BAMS+=("${BAM}")
else
  while IFS= read -r -d '' f; do
    BAMS+=("$f")
  done < <(find "${BAM_DIR}" -type f -name "Aligned.sortedByCoord.out.bam" -print0 | sort -z)

  if [[ "${#BAMS[@]}" -eq 0 ]]; then
    die "no BAM files found in ${BAM_DIR}"
  fi
fi

log "Found ${#BAMS[@]} BAM file(s)"

# -------------------------
# Run RSeQC
# -------------------------
for BAM_FILE in "${BAMS[@]}"; do
  BAM_BASE="$(basename "${BAM_FILE}")"
  SAMPLE="$(basename "$(dirname "${BAM_FILE}")")"

  if [[ "${SAMPLE}" == "." || "${SAMPLE}" == "/" || -z "${SAMPLE}" ]]; then
    SAMPLE="${BAM_BASE%.bam}"
  fi
  if [[ "${SAMPLE}" == "star" || "${SAMPLE}" == "bam" || "${SAMPLE}" == "results" ]]; then
    SAMPLE="${BAM_BASE%.bam}"
  fi

  OUTDIR="${OUTROOT}/${SAMPLE}"
  mkdir -p "${OUTDIR}"

  echo "======================================"
  log "Processing: ${BAM_FILE}"
  log "Sample:     ${SAMPLE}"
  log "Output:     ${OUTDIR}"

  if [[ ! -f "${BAM_FILE}.bai" && ! -f "${BAM_FILE%.bam}.bai" ]]; then
    log "Indexing BAM with samtools"
    samtools index "${BAM_FILE}"
  fi

  log "Running infer_experiment.py"
  infer_experiment.py -i "${BAM_FILE}" -r "${BED}" \
    > "${OUTDIR}/infer_experiment.txt" \
    2> "${OUTDIR}/infer_experiment.log" || true

  if [[ "${INFER_ONLY}" -eq 0 ]]; then
    log "Running geneBody_coverage.py"
    geneBody_coverage.py -i "${BAM_FILE}" -r "${BED}" -o "${OUTDIR}/geneBody" \
      > "${OUTDIR}/geneBody_coverage.log" 2>&1 || true
  fi
done

log "RSeQC finished."
log "Outputs written under: ${OUTROOT}"
