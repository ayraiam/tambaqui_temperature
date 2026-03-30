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
RUN_VARPART=0
RUN_DESEQ2=0
DESEQ2_COUNTS_TSV=""
DESEQ2_METADATA_TSV=""
DESEQ2_OUTDIR=""
DESEQ2_ENV_NAME="deseq2_downstream_env"
DESEQ2_ENV_FILE="envs/deseq2_downstream_env.yml"

DESEQ2_SAMPLE_COL="sample"
DESEQ2_DESIGN="~ Condition"

DESEQ2_SUBSET_COLUMN=""
DESEQ2_SUBSET_VALUES=""

DESEQ2_REFERENCE_VARIABLE=""
DESEQ2_REFERENCE_LEVEL=""

DESEQ2_CONTRAST_VARIABLE=""
DESEQ2_CONTRAST_NUMERATOR=""
DESEQ2_CONTRAST_DENOMINATOR=""

DESEQ2_ALPHA="0.05"
DESEQ2_MIN_COUNT="10"
DESEQ2_MIN_SAMPLES="3"

DESEQ2_ANNOTATION_TSV=""
DESEQ2_ANNOTATION_ID_COL="Geneid"
DESEQ2_ANNOTATION_NAME_COL="GeneName"

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
MAPQC_OUTDIR=""
MAPQC_ENV_NAME="mappingqc_var_env"
MAPQC_ENV_FILE="envs/mappingqc_var_env.yml"
MAPQC_COUNTS_TSV=""
MAPQC_METADATA_TSV=""

RUN_ENRICH=0
ENRICH_DESEQ_TSV=""
ENRICH_OUTDIR=""
ENRICH_ENV_NAME="deseq2_downstream_env"
ENRICH_ENV_FILE="envs/deseq2_downstream_env.yml"

ENRICH_SOURCE_METADATA_TSV=""
ENRICH_SOURCE_PROTEIN_FAA=""
ENRICH_TARGET_METADATA_TSV=""
ENRICH_TARGET_PROTEIN_FAA=""
ENRICH_TARGET_DMND=""

ENRICH_DESEQ_JOIN_COL="Geneid"
ENRICH_SOURCE_JOIN_COL="Gene ID"

ENRICH_ALPHA="0.05"
ENRICH_EVALUE="1e-5"
ENRICH_MAX_TARGET_SEQS="1"
ENRICH_MODE="all"

ENRICH_NORMALIZED_COUNTS_TSV=""
ENRICH_METADATA_TSV=""
ENRICH_GSEA_GO_TSV=""
ENRICH_SAMPLE_COL="sample"
ENRICH_GROUP_COL="Condition"

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
  --variance-partition    Run variance partition analysis from featureCounts + metadata

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
    --mapqc-counts-tsv PATH  featureCounts.tsv file for library QC plots
    --mapqc-metadata-tsv PATH  metadata TSV for PCA and sample-distance heatmap

DESeq2 differential expression:
    --deseq2                    Run DESeq2 differential expression
    --deseq2-counts-tsv PATH    featureCounts.tsv file
    --deseq2-metadata-tsv PATH  metadata TSV file
    --deseq2-outdir DIR         output directory
    --deseq2-env-name STR       conda env name for DESeq2 stage
    --deseq2-env-file PATH      env YAML snapshot file

    --deseq2-sample-col STR         sample column in metadata
    --deseq2-design STR             design formula, e.g. "~ Condition"
    --deseq2-subset-column STR      metadata column used for subsetting
    --deseq2-subset-values STR      comma-separated values to keep
    --deseq2-reference-variable STR variable to relevel
    --deseq2-reference-level STR    reference level
    --deseq2-contrast-variable STR  contrast variable
    --deseq2-contrast-numerator STR contrast numerator
    --deseq2-contrast-denominator STR contrast denominator
    --deseq2-alpha FLOAT            adjusted p-value cutoff
    --deseq2-min-count INT          min count threshold
    --deseq2-min-samples INT        min samples passing min count
    --deseq2-annotation-tsv PATH    optional annotation TSV
    --deseq2-annotation-id-col STR  annotation gene ID column
    --deseq2-annotation-name-col STR annotation gene name
Enrichment annotation + ORA/GSEA:
    --enrich                          Run ortholog annotation + ORA + GSEA
    --enrich-deseq-tsv PATH           DESeq2 table to annotate
    --enrich-outdir DIR               Output root for enrichment
    --enrich-env-name STR             Conda env name
    --enrich-env-file PATH            Conda env snapshot path

    --enrich-source-metadata-tsv PATH Tambaqui ncbi_dataset.tsv
    --enrich-source-protein-faa PATH  Tambaqui protein.faa
    --enrich-target-metadata-tsv PATH Danio rerio ncbi_dataset.tsv
    --enrich-target-protein-faa PATH  Danio rerio protein.faa
    --enrich-target-dmnd PATH         Optional existing DIAMOND db

    --enrich-deseq-join-col STR       Join col in DESeq table [default: Geneid]
    --enrich-source-join-col STR      Join col in source metadata [default: Gene ID]

    --enrich-alpha FLOAT              padj cutoff for ORA [default: 0.05]
    --enrich-evalue FLOAT             DIAMOND evalue [default: 1e-5]
    --enrich-max-target-seqs INT      DIAMOND max target seqs [default: 1]

    --enrich-mode STR               all|prepare|diamond|merge|analysis [default: all]

    --enrich-normalized-counts-tsv PATH  DESeq2 normalized_counts.tsv for candidate mode
    --enrich-metadata-tsv PATH           DESeq2 metadata_used.tsv for candidate mode
    --enrich-gsea-go-tsv PATH            GSEA GO BP TSV for candidate mode
    --enrich-sample-col STR              sample column in metadata [default: sample]
    --enrich-group-col STR               group column in metadata [default: Condition]

    --enrich-mode STR               all|prepare|diamond|merge|analysis|candidates [default: all]
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
    --mapqc-counts-tsv) MAPQC_COUNTS_TSV="$2"; shift 2 ;;
    --mapqc-metadata-tsv) MAPQC_METADATA_TSV="$2"; shift 2 ;;
    --variance-partition) RUN_VARPART=1; shift 1 ;;
		--deseq2) RUN_DESEQ2=1; shift 1 ;;
    --deseq2-counts-tsv) DESEQ2_COUNTS_TSV="$2"; shift 2 ;;
    --deseq2-metadata-tsv) DESEQ2_METADATA_TSV="$2"; shift 2 ;;
    --deseq2-outdir) DESEQ2_OUTDIR="$2"; shift 2 ;;
    --deseq2-env-name) DESEQ2_ENV_NAME="$2"; shift 2 ;;
    --deseq2-env-file) DESEQ2_ENV_FILE="$2"; shift 2 ;;

    --deseq2-sample-col) DESEQ2_SAMPLE_COL="$2"; shift 2 ;;
    --deseq2-design) DESEQ2_DESIGN="$2"; shift 2 ;;

    --deseq2-subset-column) DESEQ2_SUBSET_COLUMN="$2"; shift 2 ;;
    --deseq2-subset-values) DESEQ2_SUBSET_VALUES="$2"; shift 2 ;;

    --deseq2-reference-variable) DESEQ2_REFERENCE_VARIABLE="$2"; shift 2 ;;
    --deseq2-reference-level) DESEQ2_REFERENCE_LEVEL="$2"; shift 2 ;;

    --deseq2-contrast-variable) DESEQ2_CONTRAST_VARIABLE="$2"; shift 2 ;;
    --deseq2-contrast-numerator) DESEQ2_CONTRAST_NUMERATOR="$2"; shift 2 ;;
    --deseq2-contrast-denominator) DESEQ2_CONTRAST_DENOMINATOR="$2"; shift 2 ;;

    --deseq2-alpha) DESEQ2_ALPHA="$2"; shift 2 ;;
    --deseq2-min-count) DESEQ2_MIN_COUNT="$2"; shift 2 ;;
    --deseq2-min-samples) DESEQ2_MIN_SAMPLES="$2"; shift 2 ;;

    --deseq2-annotation-tsv) DESEQ2_ANNOTATION_TSV="$2"; shift 2 ;;
    --deseq2-annotation-id-col) DESEQ2_ANNOTATION_ID_COL="$2"; shift 2 ;;
    --deseq2-annotation-name-col) DESEQ2_ANNOTATION_NAME_COL="$2"; shift 2 ;;
    --enrich) RUN_ENRICH=1; shift 1 ;;
    --enrich-deseq-tsv) ENRICH_DESEQ_TSV="$2"; shift 2 ;;
    --enrich-outdir) ENRICH_OUTDIR="$2"; shift 2 ;;
    --enrich-env-name) ENRICH_ENV_NAME="$2"; shift 2 ;;
    --enrich-env-file) ENRICH_ENV_FILE="$2"; shift 2 ;;

    --enrich-source-metadata-tsv) ENRICH_SOURCE_METADATA_TSV="$2"; shift 2 ;;
    --enrich-source-protein-faa) ENRICH_SOURCE_PROTEIN_FAA="$2"; shift 2 ;;
    --enrich-target-metadata-tsv) ENRICH_TARGET_METADATA_TSV="$2"; shift 2 ;;
    --enrich-target-protein-faa) ENRICH_TARGET_PROTEIN_FAA="$2"; shift 2 ;;
    --enrich-target-dmnd) ENRICH_TARGET_DMND="$2"; shift 2 ;;

    --enrich-deseq-join-col) ENRICH_DESEQ_JOIN_COL="$2"; shift 2 ;;
    --enrich-source-join-col) ENRICH_SOURCE_JOIN_COL="$2"; shift 2 ;;

    --enrich-alpha) ENRICH_ALPHA="$2"; shift 2 ;;
    --enrich-evalue) ENRICH_EVALUE="$2"; shift 2 ;;
    --enrich-max-target-seqs) ENRICH_MAX_TARGET_SEQS="$2"; shift 2 ;;
    --enrich-mode) ENRICH_MODE="$2"; shift 2 ;;
    --enrich-normalized-counts-tsv) ENRICH_NORMALIZED_COUNTS_TSV="$2"; shift 2 ;;
    --enrich-metadata-tsv) ENRICH_METADATA_TSV="$2"; shift 2 ;;
    --enrich-gsea-go-tsv) ENRICH_GSEA_GO_TSV="$2"; shift 2 ;;
    --enrich-sample-col) ENRICH_SAMPLE_COL="$2"; shift 2 ;;
    --enrich-group-col) ENRICH_GROUP_COL="$2"; shift 2 ;;
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

# Resolve default output dir for mapping QC plots
if [[ -z "${MAPQC_OUTDIR}" ]]; then
  MAPQC_OUTDIR="${WDIR}/results/FigMappingQCandVar"
fi

# Resolve default featureCounts file for mapping QC / sample QC figures
if [[ -z "${MAPQC_COUNTS_TSV}" ]]; then
  if [[ -n "${COUNTS_DIR}" ]]; then
    MAPQC_COUNTS_TSV="${COUNTS_DIR}/featureCounts.tsv"
  else
    MAPQC_COUNTS_TSV="${RESULTS_ABS}/counts/featureCounts.tsv"
  fi
fi

# Resolve default featureCounts file for DESeq2
if [[ -z "${DESEQ2_COUNTS_TSV}" ]]; then
  if [[ -n "${COUNTS_DIR}" ]]; then
    DESEQ2_COUNTS_TSV="${COUNTS_DIR}/featureCounts.tsv"
  else
    DESEQ2_COUNTS_TSV="${RESULTS_ABS}/counts/featureCounts.tsv"
  fi
fi

if [[ -z "${DESEQ2_OUTDIR}" ]]; then
  DESEQ2_OUTDIR="${RESULTS_ABS}/deseq2"
fi

if [[ -z "${ENRICH_DESEQ_TSV}" ]]; then
  ENRICH_DESEQ_TSV="${DESEQ2_OUTDIR}/$(basename "${DESEQ2_CONTRAST_VARIABLE}_${DESEQ2_CONTRAST_NUMERATOR}_vs_${DESEQ2_CONTRAST_DENOMINATOR}_all_genes.tsv")"
fi

if [[ -z "${ENRICH_OUTDIR}" ]]; then
  ENRICH_OUTDIR="${RESULTS_ABS}/enrichment"
fi

if [[ -z "${ENRICH_NORMALIZED_COUNTS_TSV}" ]]; then
  ENRICH_NORMALIZED_COUNTS_TSV="${DESEQ2_OUTDIR}/normalized_counts.tsv"
fi

if [[ -z "${ENRICH_METADATA_TSV}" ]]; then
  ENRICH_METADATA_TSV="${DESEQ2_OUTDIR}/metadata_used.tsv"
fi

if [[ -z "${ENRICH_GSEA_GO_TSV}" ]]; then
  ENRICH_GSEA_GO_TSV="${ENRICH_OUTDIR}/03_enrichment/gsea_go_bp.tsv"
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
  echo "MAPQC_COUNTS_TSV: ${MAPQC_COUNTS_TSV}"
  echo "MAPQC_METADATA_TSV: ${MAPQC_METADATA_TSV}"
  echo "RUN_VARPART: ${RUN_VARPART}"
	echo "RUN_DESEQ2: ${RUN_DESEQ2}"
  echo "DESEQ2_COUNTS_TSV: ${DESEQ2_COUNTS_TSV}"
  echo "DESEQ2_METADATA_TSV: ${DESEQ2_METADATA_TSV}"
  echo "DESEQ2_OUTDIR: ${DESEQ2_OUTDIR}"
  echo "DESEQ2_ENV_NAME: ${DESEQ2_ENV_NAME}"
  echo "DESEQ2_ENV_FILE: ${DESEQ2_ENV_FILE}"
  echo "DESEQ2_SAMPLE_COL: ${DESEQ2_SAMPLE_COL}"
  echo "DESEQ2_DESIGN: ${DESEQ2_DESIGN}"
  echo "DESEQ2_SUBSET_COLUMN: ${DESEQ2_SUBSET_COLUMN}"
  echo "DESEQ2_SUBSET_VALUES: ${DESEQ2_SUBSET_VALUES}"
  echo "DESEQ2_REFERENCE_VARIABLE: ${DESEQ2_REFERENCE_VARIABLE}"
  echo "DESEQ2_REFERENCE_LEVEL: ${DESEQ2_REFERENCE_LEVEL}"
  echo "DESEQ2_CONTRAST_VARIABLE: ${DESEQ2_CONTRAST_VARIABLE}"
  echo "DESEQ2_CONTRAST_NUMERATOR: ${DESEQ2_CONTRAST_NUMERATOR}"
  echo "DESEQ2_CONTRAST_DENOMINATOR: ${DESEQ2_CONTRAST_DENOMINATOR}"
  echo "DESEQ2_ALPHA: ${DESEQ2_ALPHA}"
  echo "DESEQ2_MIN_COUNT: ${DESEQ2_MIN_COUNT}"
  echo "DESEQ2_MIN_SAMPLES: ${DESEQ2_MIN_SAMPLES}"
  echo "DESEQ2_ANNOTATION_TSV: ${DESEQ2_ANNOTATION_TSV}"
  echo "DESEQ2_ANNOTATION_ID_COL: ${DESEQ2_ANNOTATION_ID_COL}"
  echo "DESEQ2_ANNOTATION_NAME_COL: ${DESEQ2_ANNOTATION_NAME_COL}"
  echo "RUN_ENRICH: ${RUN_ENRICH}"
  echo "ENRICH_DESEQ_TSV: ${ENRICH_DESEQ_TSV}"
  echo "ENRICH_OUTDIR: ${ENRICH_OUTDIR}"
  echo "ENRICH_ENV_NAME: ${ENRICH_ENV_NAME}"
  echo "ENRICH_ENV_FILE: ${ENRICH_ENV_FILE}"
  echo "ENRICH_SOURCE_METADATA_TSV: ${ENRICH_SOURCE_METADATA_TSV}"
  echo "ENRICH_SOURCE_PROTEIN_FAA: ${ENRICH_SOURCE_PROTEIN_FAA}"
  echo "ENRICH_TARGET_METADATA_TSV: ${ENRICH_TARGET_METADATA_TSV}"
  echo "ENRICH_TARGET_PROTEIN_FAA: ${ENRICH_TARGET_PROTEIN_FAA}"
  echo "ENRICH_TARGET_DMND: ${ENRICH_TARGET_DMND}"
  echo "ENRICH_DESEQ_JOIN_COL: ${ENRICH_DESEQ_JOIN_COL}"
  echo "ENRICH_SOURCE_JOIN_COL: ${ENRICH_SOURCE_JOIN_COL}"
  echo "ENRICH_ALPHA: ${ENRICH_ALPHA}"
  echo "ENRICH_EVALUE: ${ENRICH_EVALUE}"
  echo "ENRICH_MAX_TARGET_SEQS: ${ENRICH_MAX_TARGET_SEQS}"
  echo "ENRICH_MODE: ${ENRICH_MODE}"
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
# STAR mapping QC / PCA / variance partition stage
# -------------------
if [[ "${RUN_MAPPING_QC_VAR}" -eq 1 || "${RUN_VARPART}" -eq 1 ]]; then
  MAPQC_STAR_DIR_FINAL="${MAPQC_STAR_DIR:-${RESULTS_ABS}/star}"

  MAPQC_ARGS=(
    --star-dir "${MAPQC_STAR_DIR_FINAL}"
    --counts-tsv "${MAPQC_COUNTS_TSV}"
    --metadata-tsv "${MAPQC_METADATA_TSV}"
    --outdir "${MAPQC_OUTDIR}"
    --env-name "${MAPQC_ENV_NAME}"
    --env-file "${MAPQC_ENV_FILE}"
  )

  [[ "${RUN_MAPPING_QC_VAR}" -eq 1 ]] && MAPQC_ARGS+=( --mapping-qc )
  [[ "${RUN_VARPART}" -eq 1 ]] && MAPQC_ARGS+=( --variance-partition )

  bash workflow/run_mapping_qc_var.sh "${MAPQC_ARGS[@]}"
fi

# -------------------
# DESeq2 differential expression stage
# -------------------
if [[ "${RUN_DESEQ2}" -eq 1 ]]; then
  DESEQ2_TS="$(date +%Y%m%d_%H%M%S)"
  DESEQ2_OUT_LOG="logs/deseq2_${DESEQ2_TS}.out"
  DESEQ2_ERR_LOG="logs/deseq2_${DESEQ2_TS}.err"

  DESEQ2_ARGS=(
    --counts-tsv "${DESEQ2_COUNTS_TSV}"
    --metadata-tsv "${DESEQ2_METADATA_TSV}"
    --outdir "${DESEQ2_OUTDIR}"
    --env-name "${DESEQ2_ENV_NAME}"
    --env-file "${DESEQ2_ENV_FILE}"
    --sample-col "${DESEQ2_SAMPLE_COL}"
    --design "${DESEQ2_DESIGN}"
    --contrast-variable "${DESEQ2_CONTRAST_VARIABLE}"
    --contrast-numerator "${DESEQ2_CONTRAST_NUMERATOR}"
    --contrast-denominator "${DESEQ2_CONTRAST_DENOMINATOR}"
    --alpha "${DESEQ2_ALPHA}"
    --min-count "${DESEQ2_MIN_COUNT}"
    --min-samples "${DESEQ2_MIN_SAMPLES}"
  )

  [[ -n "${DESEQ2_SUBSET_COLUMN}" ]] && DESEQ2_ARGS+=( --subset-column "${DESEQ2_SUBSET_COLUMN}" )
  [[ -n "${DESEQ2_SUBSET_VALUES}" ]] && DESEQ2_ARGS+=( --subset-values "${DESEQ2_SUBSET_VALUES}" )
  [[ -n "${DESEQ2_REFERENCE_VARIABLE}" ]] && DESEQ2_ARGS+=( --reference-variable "${DESEQ2_REFERENCE_VARIABLE}" )
  [[ -n "${DESEQ2_REFERENCE_LEVEL}" ]] && DESEQ2_ARGS+=( --reference-level "${DESEQ2_REFERENCE_LEVEL}" )
  [[ -n "${DESEQ2_ANNOTATION_TSV}" ]] && DESEQ2_ARGS+=( --annotation-tsv "${DESEQ2_ANNOTATION_TSV}" )
  [[ -n "${DESEQ2_ANNOTATION_ID_COL}" ]] && DESEQ2_ARGS+=( --annotation-id-col "${DESEQ2_ANNOTATION_ID_COL}" )
  [[ -n "${DESEQ2_ANNOTATION_NAME_COL}" ]] && DESEQ2_ARGS+=( --annotation-name-col "${DESEQ2_ANNOTATION_NAME_COL}" )

  echo ">>> Running DESeq2 stage"
  bash workflow/run_deseq2.sh "${DESEQ2_ARGS[@]}" \
    1> >(tee -a "${DESEQ2_OUT_LOG}") \
    2> >(tee -a "${DESEQ2_ERR_LOG}" >&2)

  echo ">>> DESeq2 logs: ${DESEQ2_OUT_LOG} / ${DESEQ2_ERR_LOG}"
fi

# -------------------
# Enrichment annotation + ORA/GSEA stage
# -------------------
if [[ "${RUN_ENRICH}" -eq 1 ]]; then
  ENRICH_TS="$(date +%Y%m%d_%H%M%S)"
  ENRICH_OUT_LOG="logs/enrich_${ENRICH_TS}.out"
  ENRICH_ERR_LOG="logs/enrich_${ENRICH_TS}.err"

  ENRICH_ARGS=(
    --deseq-tsv "${ENRICH_DESEQ_TSV}"
    --outdir "${ENRICH_OUTDIR}"
    --env-name "${ENRICH_ENV_NAME}"
    --env-file "${ENRICH_ENV_FILE}"
    --source-metadata-tsv "${ENRICH_SOURCE_METADATA_TSV}"
    --source-protein-faa "${ENRICH_SOURCE_PROTEIN_FAA}"
    --target-metadata-tsv "${ENRICH_TARGET_METADATA_TSV}"
    --target-protein-faa "${ENRICH_TARGET_PROTEIN_FAA}"
    --deseq-join-col "${ENRICH_DESEQ_JOIN_COL}"
    --source-join-col "${ENRICH_SOURCE_JOIN_COL}"
    --alpha "${ENRICH_ALPHA}"
    --threads "${CPUS}"
    --evalue "${ENRICH_EVALUE}"
    --max-target-seqs "${ENRICH_MAX_TARGET_SEQS}"
    --mode "${ENRICH_MODE}"
  )

  [[ -n "${ENRICH_TARGET_DMND}" ]] && ENRICH_ARGS+=( --target-dmnd "${ENRICH_TARGET_DMND}" )

  echo ">>> Running enrichment stage"
  bash workflow/run_enrichment.sh "${ENRICH_ARGS[@]}" \
    1> >(tee -a "${ENRICH_OUT_LOG}") \
    2> >(tee -a "${ENRICH_ERR_LOG}" >&2)

  echo ">>> Enrichment logs: ${ENRICH_OUT_LOG} / ${ENRICH_ERR_LOG}"
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
