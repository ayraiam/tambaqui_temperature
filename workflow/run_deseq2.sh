#!/usr/bin/env bash
set -euo pipefail

COUNTS_TSV=""
METADATA_TSV=""
OUTDIR=""
ENV_NAME="deseq2_downstream_env"
ENV_FILE="envs/deseq2_downstream_env.yml"

SAMPLE_COL="sample"
DESIGN="~ Condition"

SUBSET_COLUMN=""
SUBSET_VALUES=""

REFERENCE_VARIABLE=""
REFERENCE_LEVEL=""

CONTRAST_VARIABLE=""
CONTRAST_NUMERATOR=""
CONTRAST_DENOMINATOR=""

ALPHA="0.05"
MIN_COUNT="10"
MIN_SAMPLES="3"

ANNOTATION_TSV=""
ANNOTATION_ID_COL="Geneid"
ANNOTATION_NAME_COL="GeneName"

usage() {
  cat <<EOF
Usage:
  bash workflow/run_deseq2.sh [options]

Required:
  --counts-tsv PATH            featureCounts.tsv file
  --metadata-tsv PATH          metadata TSV file
  --outdir DIR                 output directory
  --contrast-variable STR      variable used in contrast (e.g. Condition)
  --contrast-numerator STR     numerator level (e.g. T2)
  --contrast-denominator STR   denominator/reference level (e.g. T1)

Optional:
  --env-name STR               conda env name [default: deseq2_downstream_env]
  --env-file PATH              env snapshot output [default: envs/deseq2_downstream_env.yml]

  --sample-col STR             sample column in metadata [default: sample]
  --design STR                 DESeq2 design formula [default: "~ Condition"]

  --subset-column STR          metadata column to subset on
  --subset-values STR          comma-separated values to keep

  --reference-variable STR     metadata variable to relevel
  --reference-level STR        reference level for that variable

  --alpha FLOAT                adjusted p-value cutoff [default: 0.05]
  --min-count INT              min count threshold for filtering [default: 10]
  --min-samples INT            min number of samples meeting min-count [default: 3]

  --annotation-tsv PATH        optional annotation table
  --annotation-id-col STR      gene ID column in annotation file [default: Geneid]
  --annotation-name-col STR    gene name column in annotation file [default: GeneName]

Examples:
  bash workflow/run_deseq2.sh \\
    --counts-tsv results/counts/featureCounts.tsv \\
    --metadata-tsv metadata/sample_metadata.tsv \\
    --outdir results/deseq2/T2_vs_T1 \\
    --subset-column Condition \\
    --subset-values T1,T2 \\
    --design "~ Condition" \\
    --reference-variable Condition \\
    --reference-level T1 \\
    --contrast-variable Condition \\
    --contrast-numerator T2 \\
    --contrast-denominator T1
EOF
  exit 0
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --counts-tsv) COUNTS_TSV="$2"; shift 2 ;;
    --metadata-tsv) METADATA_TSV="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --env-name) ENV_NAME="$2"; shift 2 ;;
    --env-file) ENV_FILE="$2"; shift 2 ;;

    --sample-col) SAMPLE_COL="$2"; shift 2 ;;
    --design) DESIGN="$2"; shift 2 ;;

    --subset-column) SUBSET_COLUMN="$2"; shift 2 ;;
    --subset-values) SUBSET_VALUES="$2"; shift 2 ;;

    --reference-variable) REFERENCE_VARIABLE="$2"; shift 2 ;;
    --reference-level) REFERENCE_LEVEL="$2"; shift 2 ;;

    --contrast-variable) CONTRAST_VARIABLE="$2"; shift 2 ;;
    --contrast-numerator) CONTRAST_NUMERATOR="$2"; shift 2 ;;
    --contrast-denominator) CONTRAST_DENOMINATOR="$2"; shift 2 ;;

    --alpha) ALPHA="$2"; shift 2 ;;
    --min-count) MIN_COUNT="$2"; shift 2 ;;
    --min-samples) MIN_SAMPLES="$2"; shift 2 ;;

    --annotation-tsv) ANNOTATION_TSV="$2"; shift 2 ;;
    --annotation-id-col) ANNOTATION_ID_COL="$2"; shift 2 ;;
    --annotation-name-col) ANNOTATION_NAME_COL="$2"; shift 2 ;;

    -h|--help) usage ;;
    *) echo "ERROR: Unknown argument: $1" >&2; usage ;;
  esac
done

[[ -n "${COUNTS_TSV}" ]] || { echo "ERROR: --counts-tsv is required" >&2; exit 2; }
[[ -f "${COUNTS_TSV}" ]] || { echo "ERROR: counts file not found: ${COUNTS_TSV}" >&2; exit 2; }

[[ -n "${METADATA_TSV}" ]] || { echo "ERROR: --metadata-tsv is required" >&2; exit 2; }
[[ -f "${METADATA_TSV}" ]] || { echo "ERROR: metadata file not found: ${METADATA_TSV}" >&2; exit 2; }

[[ -n "${OUTDIR}" ]] || { echo "ERROR: --outdir is required" >&2; exit 2; }

[[ -n "${CONTRAST_VARIABLE}" ]] || { echo "ERROR: --contrast-variable is required" >&2; exit 2; }
[[ -n "${CONTRAST_NUMERATOR}" ]] || { echo "ERROR: --contrast-numerator is required" >&2; exit 2; }
[[ -n "${CONTRAST_DENOMINATOR}" ]] || { echo "ERROR: --contrast-denominator is required" >&2; exit 2; }

if [[ -n "${ANNOTATION_TSV}" && ! -f "${ANNOTATION_TSV}" ]]; then
  echo "ERROR: annotation file not found: ${ANNOTATION_TSV}" >&2
  exit 2
fi

find_and_source_conda() {
  if command -v conda >/dev/null 2>&1; then
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    return 0
  fi

  for CAND in "$HOME/mambaforge" "$HOME/miniforge3" "$HOME/miniconda3" "/opt/conda"; do
    if [[ -f "$CAND/etc/profile.d/conda.sh" ]]; then
      export PATH="$CAND/bin:$PATH"
      # shellcheck disable=SC1091
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

create_env_if_needed() {
  if conda_env_exists "${ENV_NAME}"; then
    echo ">>> Env exists: ${ENV_NAME}"
    return 0
  fi

  echo ">>> Creating env: ${ENV_NAME}"
  ensure_channels

  local SOLVER="mamba"
  have_cmd mamba || SOLVER="conda"

  set +e
  ${SOLVER} create -n "${ENV_NAME}" -y -c conda-forge -c bioconda \
    python=3.11 \
    pandas \
    biopython \
    diamond \
    r-base=4.3 \
    r-optparse \
    r-readr \
    r-dplyr \
    r-tidyr \
    r-tibble \
    r-stringr \
    r-ggplot2 \
    r-pheatmap \
    r-rcolorbrewer \
    bioconductor-deseq2 \
    bioconductor-apeglm \
    bioconductor-clusterprofiler \
    bioconductor-org.dr.eg.db \
    bioconductor-enrichplot
  st=$?
  set -e

  if [[ $st -ne 0 ]]; then
    echo ">>> Strict solve failed; retrying with flexible channel priority"
    conda config --set channel_priority flexible
    ${SOLVER} create -n "${ENV_NAME}" -y -c conda-forge -c bioconda \
      python=3.11 \
      r-base=4.3 \
      r-optparse \
      r-readr \
      r-dplyr \
      r-tidyr \
      r-tibble \
      r-stringr \
      r-ggplot2 \
      r-pheatmap \
      r-rcolorbrewer \
      bioconductor-deseq2 \
      bioconductor-apeglm
    conda config --set channel_priority strict || true
  fi
}

export_env_snapshot() {
  mkdir -p "$(dirname "${ENV_FILE}")"
  conda env export --name "${ENV_NAME}" > "${ENV_FILE}" || true
  echo ">>> Exported env snapshot to: ${ENV_FILE}"
}

check_tools() {
  conda run -n "${ENV_NAME}" python -c "import pandas, Bio" >/dev/null 2>&1 \
    || { echo "ERROR: required Python packages not available in env ${ENV_NAME}" >&2; exit 2; }

  conda run -n "${ENV_NAME}" Rscript -e "library(optparse); library(readr); library(dplyr); library(tidyr); library(tibble); library(stringr); library(ggplot2); library(DESeq2); library(clusterProfiler); library(org.Dr.eg.db); library(enrichplot)" >/dev/null 2>&1 \
    || { echo "ERROR: required R packages not available in env ${ENV_NAME}" >&2; exit 2; }

  conda run -n "${ENV_NAME}" diamond version >/dev/null 2>&1 \
    || { echo "ERROR: diamond not available in env ${ENV_NAME}" >&2; exit 2; }
}

if ! find_and_source_conda; then
  echo "ERROR: conda not found in PATH" >&2
  exit 2
fi

create_env_if_needed
export_env_snapshot
check_tools

mkdir -p "${OUTDIR}"

CMD=(
  Rscript workflow/run_deseq2_analysis.R
  --counts-tsv "${COUNTS_TSV}"
  --metadata-tsv "${METADATA_TSV}"
  --outdir "${OUTDIR}"
  --sample-col "${SAMPLE_COL}"
  --design "${DESIGN}"
  --contrast-variable "${CONTRAST_VARIABLE}"
  --contrast-numerator "${CONTRAST_NUMERATOR}"
  --contrast-denominator "${CONTRAST_DENOMINATOR}"
  --alpha "${ALPHA}"
  --min-count "${MIN_COUNT}"
  --min-samples "${MIN_SAMPLES}"
)

[[ -n "${SUBSET_COLUMN}" ]] && CMD+=( --subset-column "${SUBSET_COLUMN}" )
[[ -n "${SUBSET_VALUES}" ]] && CMD+=( --subset-values "${SUBSET_VALUES}" )

[[ -n "${REFERENCE_VARIABLE}" ]] && CMD+=( --reference-variable "${REFERENCE_VARIABLE}" )
[[ -n "${REFERENCE_LEVEL}" ]] && CMD+=( --reference-level "${REFERENCE_LEVEL}" )

[[ -n "${ANNOTATION_TSV}" ]] && CMD+=( --annotation-tsv "${ANNOTATION_TSV}" )
[[ -n "${ANNOTATION_ID_COL}" ]] && CMD+=( --annotation-id-col "${ANNOTATION_ID_COL}" )
[[ -n "${ANNOTATION_NAME_COL}" ]] && CMD+=( --annotation-name-col "${ANNOTATION_NAME_COL}" )

echo ">>> Running DESeq2 analysis"
conda run -n "${ENV_NAME}" "${CMD[@]}"

echo ">>> Done."
echo ">>> Output dir: ${OUTDIR}"
