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
#   - dedicated Conda env creation/check for RSeQC tools
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

RSEQC_ENV_NAME="${RSEQC_ENV_NAME:-rseqc_env}"
RSEQC_ENV_FILE="${RSEQC_ENV_FILE:-envs/rseqc_env.yml}"

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

Environment:
  --env-name STR        Conda env name for RSeQC tools (default: rseqc_env)
  --env-file PATH       Conda YAML file for RSeQC env (default: envs/rseqc_env.yml)

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
    --env-name) RSEQC_ENV_NAME="$2"; shift 2 ;;
    --env-file) RSEQC_ENV_FILE="$2"; shift 2 ;;
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
  have_cmd conda || die "conda not found in PATH"
}

conda_env_exists() {
  conda env list | awk 'NR>2 {print $1}' | grep -Fxq "$1"
}

write_default_rseqc_env_file() {
  mkdir -p "$(dirname "${RSEQC_ENV_FILE}")"

  cat > "${RSEQC_ENV_FILE}" <<EOF
name: ${RSEQC_ENV_NAME}
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python=3.12
  - rseqc
  - samtools
  - ucsc-gtftogenepred
  - ucsc-genepredtobed
EOF

  log "Created default RSeQC environment file: ${RSEQC_ENV_FILE}"
}

create_rseqc_env_if_needed() {
  ensure_conda_available

  if conda_env_exists "${RSEQC_ENV_NAME}"; then
    log "RSeQC env already exists: ${RSEQC_ENV_NAME}"
    return 0
  fi

  if [[ ! -f "${RSEQC_ENV_FILE}" ]]; then
    log "RSeQC env YAML not found. Creating default one."
    write_default_rseqc_env_file
  fi

  log "Creating dedicated RSeQC env: ${RSEQC_ENV_NAME}"
  log "Using env file: ${RSEQC_ENV_FILE}"
  conda env create -f "${RSEQC_ENV_FILE}"
}

run_in_rseqc_env() {
  conda run -n "${RSEQC_ENV_NAME}" "$@"
}

check_tools_in_rseqc_env() {
  log "Checking required tools inside env: ${RSEQC_ENV_NAME}"

  conda run -n "${RSEQC_ENV_NAME}" which infer_experiment.py >/dev/null 2>&1 \
    || die "infer_experiment.py not found in env ${RSEQC_ENV_NAME}"

  conda run -n "${RSEQC_ENV_NAME}" which samtools >/dev/null 2>&1 \
    || die "samtools not found in env ${RSEQC_ENV_NAME}"

  if [[ "${MAKE_BED12}" -eq 1 ]]; then
    conda run -n "${RSEQC_ENV_NAME}" which gtfToGenePred >/dev/null 2>&1 \
      || die "gtfToGenePred not found in env ${RSEQC_ENV_NAME}"

    conda run -n "${RSEQC_ENV_NAME}" which genePredToBed >/dev/null 2>&1 \
      || die "genePredToBed not found in env ${RSEQC_ENV_NAME}"
  fi

  if [[ "${INFER_ONLY}" -eq 0 ]]; then
    conda run -n "${RSEQC_ENV_NAME}" which geneBody_coverage.py >/dev/null 2>&1 \
      || die "geneBody_coverage.py not found in env ${RSEQC_ENV_NAME}"
  fi
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
# Prepare dedicated env
# -------------------------
create_rseqc_env_if_needed
check_tools_in_rseqc_env

# -------------------------
# Resolve/build BED12
# -------------------------
if [[ "${MAKE_BED12}" -eq 1 ]]; then
  [[ -n "${GTF}" ]] || die "--make-bed12 requires --gtf"
  [[ -f "${GTF}" ]] || die "GTF not found: ${GTF}"

  if [[ -z "${BED_OUT}" ]]; then
    BED_OUT="${GTF%.*}.bed12"
  fi

  mkdir -p "$(dirname "${BED_OUT}")"
  TMP_GP="${BED_OUT%.bed12}.genepred"
  SANITIZED_GTF="${BED_OUT%.bed12}.sanitized.gtf"

  log "Sanitizing GTF: making unknown_transcript_* IDs unique per contig/strand"
  awk '
    BEGIN { FS=OFS="\t" }
    /^#/ { print; next }
    {
      attrs = $9
      if (attrs ~ /transcript_id "unknown_transcript_[^"]+"/) {
        strand = ($7 == "+" ? "plus" : ($7 == "-" ? "minus" : "unk"))
        match(attrs, /transcript_id "unknown_transcript_[^"]+"/)
        old = substr(attrs, RSTART, RLENGTH)
        tid = old
        sub(/^transcript_id "/, "", tid)
        sub(/"$/, "", tid)
        newtid = "transcript_id \"" tid "__" $1 "__" strand "\""
        sub(/transcript_id "unknown_transcript_[^"]+"/, newtid, attrs)
      }
      $9 = attrs
      print
    }
  ' "${GTF}" > "${SANITIZED_GTF}"

  log "Converting sanitized GTF -> GenePred"
  run_in_rseqc_env gtfToGenePred -ignoreGroupsWithoutExons -allErrors "${SANITIZED_GTF}" "${TMP_GP}"

  log "Converting GenePred -> BED12"
  run_in_rseqc_env genePredToBed "${TMP_GP}" "${BED_OUT}"

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
    run_in_rseqc_env samtools index "${BAM_FILE}"
  fi

  log "Running infer_experiment.py"
  run_in_rseqc_env infer_experiment.py -i "${BAM_FILE}" -r "${BED}" \
    > "${OUTDIR}/infer_experiment.txt" \
    2> "${OUTDIR}/infer_experiment.log" || true

  if [[ "${INFER_ONLY}" -eq 0 ]]; then
    log "Running geneBody_coverage.py"
    run_in_rseqc_env geneBody_coverage.py -i "${BAM_FILE}" -r "${BED}" -o "${OUTDIR}/geneBody" \
      > "${OUTDIR}/geneBody_coverage.log" 2>&1 || true
  fi
done

log "RSeQC finished."
log "Outputs written under: ${OUTROOT}"
