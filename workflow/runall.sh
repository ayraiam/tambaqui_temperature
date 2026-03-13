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
RSEQC_GTF=""
RSEQC_BAM=""
RSEQC_BAM_DIR=""
RSEQC_MAKE_BED12=0
RSEQC_BED_OUT=""
RSEQC_INFER_ONLY=1
RSEQC_GENE_BODY_ONLY=0
RSEQC_ENV_NAME="rseqc_env"
RSEQC_ENV_FILE="envs/rseqc_env.yml"
SKIP_RAW_QC=0
RUN_QC_SUMMARY_ONLY=0

# fastp knobs
FASTP_QUAL="20"
FASTP_LEN_MIN="30"
FASTP_TRIM_POLYG="${FASTP_TRIM_POLYG:-1}"
FASTP_CORRECTION="${FASTP_CORRECTION:-1}"

# screen options
USE_SCREEN=0
SCREEN_NAME="libsQC_illumina"

# STAR defaults
RUN_STAR=0
RUN_STAR_INDEX=0
RUN_FEATURECOUNTS=0
GENOME_FA=""
GTF=""
STAR_INDEX=""
READ_LENGTH="151"
STRANDNESS="0"
LIBRARY_TYPE="PE"
MAKE_BED12=0
BED12_OUT=""
COUNTS_DIR=""
TMP_DIR=""
STAR_TRIM_DIR=""

RUN_MAPPING_QC_VAR=0
MAPQC_STAR_DIR=""
MAPQC_OUTDIR="${WDIR}/results/FigMappingQCandVar"
MAPQC_ENV_NAME="mappingqc_var_env"
MAPQC_ENV_FILE="envs/mappingqc_var_env.yml"

usage() {
  cat <<EOF
Usage: bash workflow/runall.sh [options]

General:
  --cpus INT              Threads to use (default: 8)
  --wd PATH               Working directory (default: current)
  --fastq-dir PATH        Input FASTQ dir (default: data)
  --results DIR           Output root (default: results)
  --raw-qc-only           Run only FastQC+MultiQC on raw reads
  --qc-summary-only       Build only the final QC summary table
  --rseqc                 Run RSeQC as a separate stage
  --rseqc-bam PATH        Single coordinate-sorted BAM
  --rseqc-bam-dir DIR     Directory containing coordinate-sorted BAMs
  --rseqc-bed PATH        BED12 gene model for RSeQC
  --rseqc-gtf PATH        GTF annotation (used to build BED12 if requested)
  --rseqc-make-bed12      Convert GTF -> BED12 for RSeQC
  --rseqc-bed-out PATH    Output BED12 path
  --rseqc-full            Run infer_experiment.py + geneBody_coverage.py
  --gene-body-coverage    Run only geneBody_coverage.py
  --rseqc-env-name STR    Conda env name for RSeQC tools (default: rseqc_env)
  --rseqc-env-file PATH   Conda YAML file for RSeQC env

Stage control:
  --no-qc                 Skip libsQC

fastp:
  --fastp-qual INT        phred cutoff (default: 20)
  --fastp-len-min INT     min length after trimming (default: 30)

screen:
  --screen                Run libsQC inside detached screen session
  --screen-name STR       Screen session name (default: libsQC_illumina)

STAR mapping / counting:
    --star                  Run STAR mapping on trimmed reads
    --star-index            Build STAR genome index
    --featurecounts         Run featureCounts on existing STAR BAMs
    --genome-fa PATH        Genome FASTA file
    --gtf PATH              Annotation GTF file
    --star-index-dir PATH   Directory for STAR genome index
    --read-length INT       Read length (default: 151)
    --strandness 0|1|2      featureCounts strandedness (default: 0)
    --make-bed12            Create BED12 from GTF for RSeQC
    --bed12-out PATH        Output BED12 file path
    --library-type PE|SE    Library type for featureCounts (default: PE)
    --counts-dir PATH       Output directory for featureCounts files
    --tmp-dir PATH          Temporary directory for featureCounts

Making Mapping QC & variance Partition figure:
    --mapping-qc-var       Parse STAR Log.final.out files and plot mapping QC
    --mapqc-star-dir DIR   STAR directory containing sample subdirs with Log.final.out
    --mapqc-outdir DIR     Output dir for parser tables and plots
    --mapqc-env-name STR   Conda env name for mapping QC plots
    --mapqc-env-file PATH  Conda YAML file for mapping QC plots

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
    --rseqc-bam) RSEQC_BAM="$2"; shift 2 ;;
    --rseqc-bam-dir) RSEQC_BAM_DIR="$2"; shift 2 ;;
    --rseqc-bed) RSEQC_BED="$2"; shift 2 ;;
    --rseqc-gtf) RSEQC_GTF="$2"; shift 2 ;;
    --rseqc-make-bed12) RSEQC_MAKE_BED12=1; shift 1 ;;
    --rseqc-bed-out) RSEQC_BED_OUT="$2"; shift 2 ;;
    --rseqc-full) RSEQC_INFER_ONLY=0; RSEQC_GENE_BODY_ONLY=0; shift 1 ;;
    --gene-body-coverage) RUN_RSEQC=1; RSEQC_INFER_ONLY=0; RSEQC_GENE_BODY_ONLY=1; shift 1 ;;
    --rseqc-env-name) RSEQC_ENV_NAME="$2"; shift 2 ;;
    --rseqc-env-file) RSEQC_ENV_FILE="$2"; shift 2 ;;
    --skip-raw-qc) SKIP_RAW_QC=1; shift ;;
    --star) RUN_STAR=1; shift ;;
    --star-index) RUN_STAR_INDEX=1; shift ;;
    --featurecounts) RUN_FEATURECOUNTS=1; shift ;;
    --genome-fa) GENOME_FA="$2"; shift 2 ;;
    --gtf) GTF="$2"; shift 2 ;;
    --star-index-dir) STAR_INDEX="$2"; shift 2 ;;
    --trim-dir) STAR_TRIM_DIR="$2"; shift 2 ;;
    --read-length) READ_LENGTH="$2"; shift 2 ;;
    --strandness) STRANDNESS="$2"; shift 2 ;;
    --make-bed12) MAKE_BED12=1; shift ;;
    --bed12-out) BED12_OUT="$2"; shift 2 ;;
    --qc-summary-only) RUN_QC_SUMMARY_ONLY=1; shift 1 ;;
    --library-type) LIBRARY_TYPE="$2"; shift 2 ;;
    --counts-dir) COUNTS_DIR="$2"; shift 2 ;;
    --tmp-dir) TMP_DIR="$2"; shift 2 ;;
    --mapping-qc-var) RUN_MAPPING_QC_VAR=1; shift 1 ;;
    --mapqc-star-dir) MAPQC_STAR_DIR="$2"; shift 2 ;;
    --mapqc-outdir) MAPQC_OUTDIR="$2"; shift 2 ;;
    --mapqc-env-name) MAPQC_ENV_NAME="$2"; shift 2 ;;
    --mapqc-env-file) MAPQC_ENV_FILE="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "Unknown argument: $1"; usage ;;
  esac
done

# Resolve absolute results directory
if [[ "${RESULTS}" = /* ]]; then
  RESULTS_ABS="${RESULTS}"
else
  RESULTS_ABS="${WDIR}/${RESULTS}"
fi

mkdir -p logs metadata

TS="$(date +%Y%m%d_%H%M%S)"
OUT_LOG="logs/libsQC_${TS}.out"
ERR_LOG="logs/libsQC_${TS}.err"
INVOCATION_LOG="logs/invocation_${TS}.txt"

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
  echo "RUN_QC_SUMMARY_ONLY: ${RUN_QC_SUMMARY_ONLY}"
  echo "RUN_STAR: ${RUN_STAR}"
  echo "RUN_STAR_INDEX: ${RUN_STAR_INDEX}"
  echo "RUN_FEATURECOUNTS: ${RUN_FEATURECOUNTS}"
  echo "RUN_RSEQC: ${RUN_RSEQC}"
  echo "RSEQC_INFER_ONLY: ${RSEQC_INFER_ONLY}"
  echo "RSEQC_GENE_BODY_ONLY: ${RSEQC_GENE_BODY_ONLY}"
  echo "STRANDNESS: ${STRANDNESS}"
  echo "LIBRARY_TYPE: ${LIBRARY_TYPE}"
  echo "COUNTS_DIR: ${COUNTS_DIR}"
  echo "TMP_DIR: ${TMP_DIR}"
  echo "STAR_TRIM_DIR: ${STAR_TRIM_DIR}"
  echo "RUN_MAPPING_QC_VAR: ${RUN_MAPPING_QC_VAR}"
  echo "MAPQC_STAR_DIR: ${MAPQC_STAR_DIR}"
  echo "MAPQC_OUTDIR: ${MAPQC_OUTDIR}"
  echo "MAPQC_ENV_NAME: ${MAPQC_ENV_NAME}"
  echo "MAPQC_ENV_FILE: ${MAPQC_ENV_FILE}"
  echo "=========================================="
} > "$INVOCATION_LOG"

echo ">>> Invocation logged to: ${INVOCATION_LOG}"

if [[ "${RUN_LIBSQC}" -ne 1 ]]; then
  echo ">>> Skipping libsQC (--no-qc)"
else
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
export SKIP_RAW_QC="$SKIP_RAW_QC"
export RUN_STAR="$RUN_STAR"
export RUN_STAR_INDEX="$RUN_STAR_INDEX"
export RUN_FEATURECOUNTS="$RUN_FEATURECOUNTS"
export GENOME_FA="$GENOME_FA"
export GTF="$GTF"
export STAR_INDEX="$STAR_INDEX"
export READ_LENGTH="$READ_LENGTH"
export STRANDNESS="$STRANDNESS"
export MAKE_BED12="$MAKE_BED12"
export BED12_OUT="$BED12_OUT"
export RUN_QC_SUMMARY_ONLY="$RUN_QC_SUMMARY_ONLY"
export LIBRARY_TYPE="$LIBRARY_TYPE"
export COUNTS_DIR="$COUNTS_DIR"
export TMP_DIR="$TMP_DIR"
bash workflow/run_libsQC_illumina.sh
EOF
  )

  if [[ "$USE_SCREEN" -eq 1 ]]; then
    command -v screen >/dev/null 2>&1 || { echo "!!! screen not found in PATH"; exit 2; }
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
fi

# -------------------
# STAR / featureCounts stage
# -------------------
if [[ "${RUN_STAR_INDEX}" -eq 1 || "${RUN_STAR}" -eq 1 || "${RUN_FEATURECOUNTS}" -eq 1 || "${MAKE_BED12}" -eq 1 ]]; then
  STAR_ARGS=(
    --genome-fa "${GENOME_FA}"
    --gtf "${GTF}"
    --star-index-dir "${STAR_INDEX}"
    --threads "${CPUS}"
    --results "${RESULTS_ABS}"
    --trim-dir "${STAR_TRIM_DIR:-${RESULTS_ABS}/trimmed}"
    --read-length "${READ_LENGTH}"
    --strandness "${STRANDNESS}"
    --library-type "${LIBRARY_TYPE}"
  )

  [[ "${RUN_STAR_INDEX}" -eq 1 ]] && STAR_ARGS+=( --index )
  [[ "${RUN_STAR}" -eq 1 ]] && STAR_ARGS+=( --map )
  [[ "${RUN_FEATURECOUNTS}" -eq 1 ]] && STAR_ARGS+=( --counts )
  [[ "${MAKE_BED12}" -eq 1 ]] && STAR_ARGS+=( --make-bed12 )
  [[ -n "${BED12_OUT}" ]] && STAR_ARGS+=( --bed12-out "${BED12_OUT}" )
  [[ -n "${COUNTS_DIR}" ]] && STAR_ARGS+=( --counts-dir "${COUNTS_DIR}" )
  [[ -n "${TMP_DIR}" ]] && STAR_ARGS+=( --tmp-dir "${TMP_DIR}" )

  bash workflow/run_star.sh "${STAR_ARGS[@]}"
fi

# -------------------
# STAR mapping QC summary / plotting stage
# -------------------
if [[ "${RUN_MAPPING_QC_VAR}" -eq 1 ]]; then
  MAPQC_STAR_DIR_FINAL="${MAPQC_STAR_DIR:-${RESULTS_ABS}/star}"

  bash workflow/run_mapping_qc_var.sh \
    --star-dir "${MAPQC_STAR_DIR_FINAL}" \
    --outdir "${MAPQC_OUTDIR}" \
    --env-name "${MAPQC_ENV_NAME}" \
    --env-file "${MAPQC_ENV_FILE}"
fi

# -------------------
# RSeQC stage
# -------------------
if [[ "${RUN_RSEQC}" -eq 1 ]]; then
  RSEQC_ARGS=(
    --outroot "${RESULTS_ABS}/rseqc"
  )

  [[ -n "${RSEQC_BAM}" ]] && RSEQC_ARGS+=( --bam "${RSEQC_BAM}" )
  [[ -n "${RSEQC_BAM_DIR}" ]] && RSEQC_ARGS+=( --bam-dir "${RSEQC_BAM_DIR}" )
  [[ -n "${RSEQC_BED}" ]] && RSEQC_ARGS+=( --bed "${RSEQC_BED}" )
  [[ -n "${RSEQC_GTF}" ]] && RSEQC_ARGS+=( --gtf "${RSEQC_GTF}" )
  [[ "${RSEQC_MAKE_BED12}" -eq 1 ]] && RSEQC_ARGS+=( --make-bed12 )
  [[ -n "${RSEQC_BED_OUT}" ]] && RSEQC_ARGS+=( --bed-out "${RSEQC_BED_OUT}" )
  [[ "${RSEQC_INFER_ONLY}" -eq 1 ]] && RSEQC_ARGS+=( --infer-only )
  [[ "${RSEQC_GENE_BODY_ONLY}" -eq 1 ]] && RSEQC_ARGS+=( --gene-body-only )
  [[ -n "${RSEQC_ENV_NAME}" ]] && RSEQC_ARGS+=( --env-name "${RSEQC_ENV_NAME}" )
  [[ -n "${RSEQC_ENV_FILE}" ]] && RSEQC_ARGS+=( --env-file "${RSEQC_ENV_FILE}" )

  bash workflow/run_rseqc.sh "${RSEQC_ARGS[@]}"
fi
