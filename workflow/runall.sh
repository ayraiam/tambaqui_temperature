#!/usr/bin/env bash
# ==========================================================
# Script: workflow/runall.sh
# Purpose: Run Illumina libsQC locally (optionally inside screen)
# ==========================================================
set -euo pipefail

ORIG_ARGS=("$@")

# Defaults
CPUS="8"
WDIR="$PWD"
FASTQ_DIR="data"
RESULTS="results"
RUN_LIBSQC=1
RAW_QC_ONLY=0
RUN_RSEQC=0
RSEQC_BED=""
RSEQC_BAM_DIR=""
SKIP_RAW_QC=0

# fastp knobs
FASTP_QUAL="20"
FASTP_LEN_MIN="30"
FASTP_TRIM_POLYG="${FASTP_TRIM_POLYG:-1}"
FASTP_CORRECTION="${FASTP_CORRECTION:-1}"

# screen options
USE_SCREEN=0
SCREEN_NAME="libsQC_illumina"

usage() {
  cat <<EOF
Usage: bash workflow/runall.sh [options]

General:
  --cpus INT            Threads to use (default: 8)
  --wd PATH             Working directory (default: current)
  --fastq-dir PATH      Input FASTQ dir (default: data)
  --results DIR         Output root (default: results)
  --raw-qc-only         Run only FastQC+MultiQC on raw reads (skip fastp + trimmed QC + seqkit)
  --rseqc              Run RSeQC checks (requires BAMs + BED12)
  --rseqc-bed PATH     Gene model in BED12 format for RSeQC
  --rseqc-bam-dir DIR  Directory containing coordinate-sorted BAMs

Stage control:
  --no-qc               Skip libsQC

fastp:
  --fastp-qual INT      phred cutoff (default: 20)
  --fastp-len-min INT   min length after trimming (default: 30)

screen:
  --screen              Run libsQC inside a detached screen session
  --screen-name STR     Screen session name (default: libsQC_illumina)

Notes:
  You can also toggle via environment variables:
    FASTP_TRIM_POLYG=0|1
    FASTP_CORRECTION=0|1
EOF
  exit 0
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --cpus) CPUS="$2"; shift 2 ;;
    --wd) WDIR="$2"; shift 2 ;;
    --fastq-dir) FASTQ_DIR="$2"; shift 2 ;;
    --results) RESULTS="$2"; shift 2 ;;
    --no-qc) RUN_LIBSQC=0; shift 1 ;;
    --fastp-qual) FASTP_QUAL="$2"; shift 2 ;;
    --fastp-len-min) FASTP_LEN_MIN="$2"; shift 2 ;;
    --screen) USE_SCREEN=1; shift 1 ;;
    --screen-name) SCREEN_NAME="$2"; shift 2 ;;
    --raw-qc-only) RAW_QC_ONLY=1; shift 1 ;;
    --rseqc) RUN_RSEQC=1; shift ;;
    --rseqc-bed) RSEQC_BED="$2"; shift 2 ;;
    --rseqc-bam-dir) RSEQC_BAM_DIR="$2"; shift 2 ;;
    --skip-raw-qc) SKIP_RAW_QC=1; shift ;;
    -h|--help) usage ;;
    *) echo "Unknown argument: $1"; usage ;;
  esac
done

mkdir -p logs metadata

TS="$(date +%Y%m%d_%H%M%S)"
OUT_LOG="logs/libsQC_${TS}.out"
ERR_LOG="logs/libsQC_${TS}.err"
INVOCATION_LOG="logs/invocation_${TS}.txt"

# Log invocation
{
  echo "=========================================="
  echo "Timestamp : $(date)"
  echo "Host      : $(hostname)"
  echo "PWD       : $PWD"
  echo "Command   : $0 ${ORIG_ARGS[*]}"
  echo "------------------------------------------"
  echo "CPUS      : ${CPUS}"
  echo "WDIR      : ${WDIR}"
  echo "FASTQ_DIR : ${FASTQ_DIR}"
  echo "RESULTS   : ${RESULTS}"
  echo "FASTP_QUAL: ${FASTP_QUAL}"
  echo "FASTP_LEN_MIN: ${FASTP_LEN_MIN}"
  echo "FASTP_TRIM_POLYG: ${FASTP_TRIM_POLYG}"
  echo "FASTP_CORRECTION: ${FASTP_CORRECTION}"
  echo "USE_SCREEN: ${USE_SCREEN}"
  echo "SCREEN_NAME: ${SCREEN_NAME}"
  echo "=========================================="
} > "$INVOCATION_LOG"

echo ">>> Invocation logged to: ${INVOCATION_LOG}"

if [[ "${RUN_LIBSQC}" -ne 1 ]]; then
  echo ">>> Skipping libsQC (--no-qc)"
  exit 0
fi

# Command we want to run (export env vars like Slurm --export used to)
RUN_CMD=$(
  cat <<EOF
set -euo pipefail
cd "$WDIR"
export THREADS="$CPUS"
export FASTQ_DIR="$FASTQ_DIR"
export RESULTS="$RESULTS"
export FASTP_QUAL="$FASTP_QUAL"
export FASTP_LEN_MIN="$FASTP_LEN_MIN"
export FASTP_TRIM_POLYG="$FASTP_TRIM_POLYG"
export FASTP_CORRECTION="$FASTP_CORRECTION"
export RAW_QC_ONLY="$RAW_QC_ONLY"
export RUN_RSEQC="$RUN_RSEQC"
export RSEQC_BED="$RSEQC_BED"
export RSEQC_BAM_DIR="$RSEQC_BAM_DIR"
export SKIP_RAW_QC="$SKIP_RAW_QC"
bash workflow/run_libsQC_illumina.sh
EOF
)

if [[ "$USE_SCREEN" -eq 1 ]]; then
  command -v screen >/dev/null 2>&1 || { echo "!!! screen not found in PATH"; exit 2; }

  # Start detached screen session running the pipeline, logging stdout/stderr
  screen -S "$SCREEN_NAME" -dm bash -lc "$RUN_CMD" \
    1> >(tee -a "$OUT_LOG") \
    2> >(tee -a "$ERR_LOG" >&2) || true

  echo ">>> Started in screen session: ${SCREEN_NAME}"
  echo "    Attach:  screen -r ${SCREEN_NAME}"
  echo "    Logs  :  ${OUT_LOG} / ${ERR_LOG}"
else
  echo ">>> Running libsQC in foreground..."
  bash -lc "$RUN_CMD" 1> >(tee -a "$OUT_LOG") 2> >(tee -a "$ERR_LOG" >&2)
  echo ">>> Done."
  echo "    Logs  :  ${OUT_LOG} / ${ERR_LOG}"
fi
