#!/usr/bin/env bash
set -euo pipefail

DESEQ_TSV=""
OUTDIR=""
ENV_NAME="deseq2_downstream_env"
ENV_FILE="envs/deseq2_downstream_env.yml"

SOURCE_METADATA_TSV=""
SOURCE_PROTEIN_FAA=""
TARGET_METADATA_TSV=""
TARGET_PROTEIN_FAA=""
TARGET_DMND=""

DESEQ_JOIN_COL="Geneid"
SOURCE_JOIN_COL="Gene ID"

ALPHA="0.05"
THREADS="8"
EVALUE="1e-5"
MAX_TARGET_SEQS="1"

MODE="all"

usage() {
  cat <<EOF
Usage:
  bash workflow/run_enrichment.sh [options]

Required:
  --deseq-tsv PATH
  --outdir DIR
  --source-metadata-tsv PATH
  --source-protein-faa PATH
  --target-metadata-tsv PATH
  --target-protein-faa PATH
  --mode STR                all|prepare|diamond|merge|analysis [default: all]

Optional:
  --target-dmnd PATH
  --env-name STR
  --env-file PATH
  --deseq-join-col STR      [default: Geneid]
  --source-join-col STR     [default: Gene ID]
  --alpha FLOAT             [default: 0.05]
  --threads INT             [default: 8]
  --evalue FLOAT            [default: 1e-5]
  --max-target-seqs INT     [default: 1]
EOF
  exit 0
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --deseq-tsv) DESEQ_TSV="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --env-name) ENV_NAME="$2"; shift 2 ;;
    --env-file) ENV_FILE="$2"; shift 2 ;;

    --source-metadata-tsv) SOURCE_METADATA_TSV="$2"; shift 2 ;;
    --source-protein-faa) SOURCE_PROTEIN_FAA="$2"; shift 2 ;;
    --target-metadata-tsv) TARGET_METADATA_TSV="$2"; shift 2 ;;
    --target-protein-faa) TARGET_PROTEIN_FAA="$2"; shift 2 ;;
    --target-dmnd) TARGET_DMND="$2"; shift 2 ;;

    --deseq-join-col) DESEQ_JOIN_COL="$2"; shift 2 ;;
    --source-join-col) SOURCE_JOIN_COL="$2"; shift 2 ;;

    --alpha) ALPHA="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --evalue) EVALUE="$2"; shift 2 ;;
    --max-target-seqs) MAX_TARGET_SEQS="$2"; shift 2 ;;
    --mode) MODE="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "ERROR: Unknown argument: $1" >&2; usage ;;
  esac
done

for f in "$DESEQ_TSV" "$SOURCE_METADATA_TSV" "$SOURCE_PROTEIN_FAA" "$TARGET_METADATA_TSV" "$TARGET_PROTEIN_FAA"; do
  [[ -n "$f" ]] || { echo "ERROR: missing required file argument" >&2; exit 2; }
  [[ -f "$f" ]] || { echo "ERROR: file not found: $f" >&2; exit 2; }
done

[[ -n "${OUTDIR}" ]] || { echo "ERROR: --outdir is required" >&2; exit 2; }

find_and_source_conda() {
  if command -v conda >/dev/null 2>&1; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
    return 0
  fi
  for CAND in "$HOME/mambaforge" "$HOME/miniforge3" "$HOME/miniconda3" "/opt/conda"; do
    if [[ -f "$CAND/etc/profile.d/conda.sh" ]]; then
      export PATH="$CAND/bin:$PATH"
      source "$CAND/etc/profile.d/conda.sh"
      command -v conda >/dev/null 2>&1 && return 0
    fi
  done
  return 1
}

have_cmd() {
  command -v "$1" >/dev/null 2>&1
}

conda_env_exists() {
  conda env list | awk 'NR>2 {print $1}' | grep -Fxq "$1"
}

ensure_channels() {
  conda config --remove-key channels >/dev/null 2>&1 || true
  conda config --add channels conda-forge
  conda config --add channels bioconda
  conda config --add channels defaults
  conda config --set channel_priority strict
  conda clean --index-cache -y >/dev/null 2>&1 || true
}

create_or_update_env() {
  ensure_channels
  local SOLVER="mamba"
  have_cmd mamba || SOLVER="conda"

  local PKGS=(
    python=3.11
    pandas
    biopython
    diamond
    r-base=4.3
    r-optparse
    r-readr
    r-dplyr
    r-tidyr
    r-tibble
    r-stringr
    r-ggplot2
    r-pheatmap
    r-rcolorbrewer
    bioconductor-deseq2
    bioconductor-apeglm
    bioconductor-clusterprofiler
    bioconductor-org.dr.eg.db
    bioconductor-enrichplot
  )

  if conda_env_exists "${ENV_NAME}"; then
    echo ">>> Env exists: ${ENV_NAME} (ensuring enrichment packages)"
    ${SOLVER} install -n "${ENV_NAME}" -y -c conda-forge -c bioconda "${PKGS[@]}"
  else
    echo ">>> Creating env: ${ENV_NAME}"
    ${SOLVER} create -n "${ENV_NAME}" -y -c conda-forge -c bioconda "${PKGS[@]}"
  fi
}

export_env_snapshot() {
  mkdir -p "$(dirname "${ENV_FILE}")"
  conda env export --name "${ENV_NAME}" > "${ENV_FILE}" || true
  echo ">>> Exported env snapshot to: ${ENV_FILE}"
}

check_tools() {
  conda run -n "${ENV_NAME}" python -c "import pandas, Bio" >/dev/null 2>&1 \
    || { echo "ERROR: Python enrichment packages unavailable" >&2; exit 2; }

  conda run -n "${ENV_NAME}" Rscript -e "library(clusterProfiler); library(org.Dr.eg.db); library(enrichplot)" >/dev/null 2>&1 \
    || { echo "ERROR: R enrichment packages unavailable" >&2; exit 2; }

  conda run -n "${ENV_NAME}" diamond version >/dev/null 2>&1 \
    || { echo "ERROR: diamond unavailable" >&2; exit 2; }
}

if ! find_and_source_conda; then
  echo "ERROR: conda not found in PATH" >&2
  exit 2
fi

create_or_update_env
export_env_snapshot
check_tools

mkdir -p "${OUTDIR}"
PREP_DIR="${OUTDIR}/01_prepare"
DIAMOND_DIR="${OUTDIR}/02_diamond"
ENRICH_DIR="${OUTDIR}/03_enrichment"
mkdir -p "${PREP_DIR}" "${DIAMOND_DIR}" "${ENRICH_DIR}"

echo ">>> Enrichment mode: ${MODE}"

if [[ "${MODE}" == "all" || "${MODE}" == "prepare" ]]; then
  echo ">>> Step 1: prepare source annotation + query FASTA"
  conda run -n "${ENV_NAME}" python workflow/prepare_enrichment_inputs.py \
    --deseq-tsv "${DESEQ_TSV}" \
    --source-metadata-tsv "${SOURCE_METADATA_TSV}" \
    --source-protein-faa "${SOURCE_PROTEIN_FAA}" \
    --outdir "${PREP_DIR}" \
    --deseq-join-col "${DESEQ_JOIN_COL}" \
    --source-join-col "${SOURCE_JOIN_COL}"
fi

if [[ "${MODE}" == "prepare" ]]; then
  echo ">>> Done."
  echo ">>> Prepare outputs: ${PREP_DIR}"
  exit 0
fi

if [[ -z "${TARGET_DMND}" ]]; then
  TARGET_DMND="${DIAMOND_DIR}/target_proteins.dmnd"
fi

if [[ "${MODE}" == "all" || "${MODE}" == "diamond" ]]; then
  if [[ ! -f "${TARGET_DMND}" ]]; then
    echo ">>> Step 2a: build DIAMOND database"
    conda run -n "${ENV_NAME}" diamond makedb \
      --in "${TARGET_PROTEIN_FAA}" \
      --db "${TARGET_DMND}"
  fi

  echo ">>> Step 2b: run DIAMOND blastp"
  conda run -n "${ENV_NAME}" diamond blastp \
    --db "${TARGET_DMND}" \
    --query "${PREP_DIR}/query_proteins.faa" \
    --out "${DIAMOND_DIR}/diamond_best_hits.tsv" \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovhsp \
    --evalue "${EVALUE}" \
    --max-target-seqs "${MAX_TARGET_SEQS}" \
    --threads "${THREADS}"
fi

if [[ "${MODE}" == "diamond" ]]; then
  echo ">>> Done."
  echo ">>> DIAMOND outputs: ${DIAMOND_DIR}"
  exit 0
fi

if [[ "${MODE}" == "all" || "${MODE}" == "merge" ]]; then
  echo ">>> Step 3: merge DIAMOND hits with zebrafish metadata"
  conda run -n "${ENV_NAME}" python workflow/merge_diamond_hits.py \
    --source-annotation-tsv "${PREP_DIR}/source_annotation_from_deseq.tsv" \
    --diamond-tsv "${DIAMOND_DIR}/diamond_best_hits.tsv" \
    --target-metadata-tsv "${TARGET_METADATA_TSV}" \
    --outdir "${DIAMOND_DIR}"
fi

if [[ "${MODE}" == "merge" ]]; then
  echo ">>> Done."
  echo ">>> Merge outputs: ${DIAMOND_DIR}"
  exit 0
fi

if [[ "${MODE}" == "all" || "${MODE}" == "analysis" ]]; then
  echo ">>> Step 4: run ORA + GSEA"
  conda run -n "${ENV_NAME}" Rscript workflow/run_enrichment_analysis.R \
    --annotated-tsv "${DIAMOND_DIR}/deseq_annotated_with_danio_hits.tsv" \
    --outdir "${ENRICH_DIR}" \
    --alpha "${ALPHA}"
fi

echo ">>> Done."
echo ">>> Enrichment root: ${OUTDIR}"
