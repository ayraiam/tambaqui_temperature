#!/usr/bin/env bash
set -euo pipefail

STAR_DIR=""
COUNTS_TSV=""
OUTDIR=""
ENV_NAME="mappingqc_var_env"
ENV_FILE="envs/mappingqc_var_env.yml"

usage() {
  cat <<EOF
Usage:
  bash workflow/run_mapping_qc_var.sh [options]

Required:
  --star-dir DIR       STAR directory with sample subdirs containing Log.final.out
  --counts-tsv PATH    featureCounts.tsv file
  --outdir DIR         Output directory for parsed tables and plots

Optional:
  --env-name STR       Conda env name (default: mappingqc_var_env)
  --env-file PATH      YAML file to write/export env snapshot (default: envs/mappingqc_var_env.yml)
EOF
  exit 0
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --star-dir) STAR_DIR="$2"; shift 2 ;;
    --counts-tsv) COUNTS_TSV="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --env-name) ENV_NAME="$2"; shift 2 ;;
    --env-file) ENV_FILE="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "Unknown argument: $1" >&2; usage ;;
  esac
done

[[ -n "${STAR_DIR}" ]] || { echo "ERROR: --star-dir is required" >&2; exit 2; }
[[ -d "${STAR_DIR}" ]] || { echo "ERROR: STAR directory not found: ${STAR_DIR}" >&2; exit 2; }
[[ -n "${COUNTS_TSV}" ]] || { echo "ERROR: --counts-tsv is required" >&2; exit 2; }
[[ -f "${COUNTS_TSV}" ]] || { echo "ERROR: featureCounts file not found: ${COUNTS_TSV}" >&2; exit 2; }
[[ -n "${OUTDIR}" ]] || { echo "ERROR: --outdir is required" >&2; exit 2; }

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
    r-base=4.3 \
    r-ggbeeswarm \
    r-readr \
    r-dplyr \
    r-tidyr \
    r-ggplot2
  st=$?
  set -e

  if [[ $st -ne 0 ]]; then
    echo ">>> Strict solve failed; retrying with flexible channel priority"
    conda config --set channel_priority flexible
    ${SOLVER} create -n "${ENV_NAME}" -y -c conda-forge -c bioconda \
      python=3.11 \
      pandas \
      r-base=4.3 \
      r-ggbeeswarm \
      r-readr \
      r-dplyr \
      r-tidyr \
      r-ggplot2
    conda config --set channel_priority strict || true
  fi
}

export_env_snapshot() {
  mkdir -p "$(dirname "${ENV_FILE}")"
  conda env export --name "${ENV_NAME}" > "${ENV_FILE}" || true
  echo ">>> Exported env snapshot to: ${ENV_FILE}"
}

check_tools() {
  conda run -n "${ENV_NAME}" python -c "import pandas" >/dev/null 2>&1 \
    || { echo "ERROR: pandas not available in env ${ENV_NAME}" >&2; exit 2; }

  conda run -n "${ENV_NAME}" Rscript -e "library(readr); library(dplyr); library(tidyr); library(ggplot2); library(ggbeeswarm)" >/dev/null 2>&1 \
    || { echo "ERROR: required R packages not available in env ${ENV_NAME}" >&2; exit 2; }
}

if ! find_and_source_conda; then
  echo "ERROR: conda not found in PATH" >&2
  exit 2
fi

create_env_if_needed
export_env_snapshot
check_tools

mkdir -p "${OUTDIR}"

echo ">>> Parsing STAR Log.final.out files"
conda run -n "${ENV_NAME}" python workflow/parse_star_log_final.py \
  --star-dir "${STAR_DIR}" \
  --outdir "${OUTDIR}"

echo ">>> Plotting STAR mapping QC + featureCounts sample QC"
conda run -n "${ENV_NAME}" Rscript workflow/plot_star_mapping_qc.R \
  "${OUTDIR}/star_mapping_qc_summary.tsv" \
  "${OUTDIR}" \
  "${COUNTS_TSV}"

echo ">>> Done."
echo ">>> Output dir: ${OUTDIR}"
