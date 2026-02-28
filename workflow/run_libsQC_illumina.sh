#!/usr/bin/env bash
# ==========================================================
# Script: workflow/run_libsQC_illumina.sh
# Purpose: Illumina short-read QC + trimming using:
#          FastQC + MultiQC + fastp + seqkit (+ pigz, samtools)
# Inputs : data/*.fastq(.gz), data/*.fq(.gz)  (PE auto-detect)
# Outputs:
#   results/qc_raw/
#   results/multiqc_raw/
#   results/trimmed/
#   results/fastp/
#   results/qc_trimmed/
#   results/multiqc_trimmed/
#   results/summary/
#   metadata/fastq_meta.tsv
#   envs/libsQC_illumina.yml
# ==========================================================
set -euo pipefail

# --------------------------
# User-tunable via env vars
# --------------------------
THREADS="${THREADS:-8}"
RESULTS="${RESULTS:-results}"
FASTQ_DIR="${FASTQ_DIR:-data}"

# fastp defaults (override via env if desired)
FASTP_QUAL="${FASTP_QUAL:-20}"       # phred cutoff
FASTP_LEN_MIN="${FASTP_LEN_MIN:-30}" # min length after trimming
FASTP_TRIM_POLYG="${FASTP_TRIM_POLYG:-1}"  # NextSeq/NovaSeq polyG tail
FASTP_CORRECTION="${FASTP_CORRECTION:-1}"  # overlap correction for PE

RAW_QC_ONLY="${RAW_QC_ONLY:-0}"

# -----------------------------------------
# Conda bootstrap (robust)
# -----------------------------------------
find_and_source_conda() {
  # If conda already in PATH, just source conda.sh
  if command -v conda >/dev/null 2>&1; then
    # shellcheck disable=SC1090
    source "$(conda info --base)/etc/profile.d/conda.sh"
    return 0
  fi

  # Try common installations
  for CAND in "$HOME/mambaforge" "$HOME/miniforge3" "$HOME/miniconda3" "/opt/conda"; do
    if [[ -f "$CAND/etc/profile.d/conda.sh" ]]; then
      export PATH="$CAND/bin:$PATH"
      # shellcheck disable=SC1090
      source "$CAND/etc/profile.d/conda.sh"
      command -v conda >/dev/null 2>&1 && return 0
    fi
  done

  return 1
}

if ! find_and_source_conda; then
  echo "!!! conda not found."
  echo "    Looked in PATH and in:"
  echo "    \$HOME/mambaforge"
  echo "    \$HOME/miniforge3"
  echo "    \$HOME/miniconda3"
  echo "    /opt/conda"
  exit 2
fi

ensure_channels() {
  conda config --remove-key channels >/dev/null 2>&1 || true
  conda config --add channels conda-forge
  conda config --add channels bioconda
  conda config --add channels defaults
  conda config --set channel_priority strict
  conda clean --index-cache -y >/dev/null 2>&1 || true
}

create_env_libsQC_illumina() {
  local ENV="libsQC_illumina"

  if conda env list | grep -qE "^${ENV}\s"; then
    echo ">>> Env ${ENV} exists. Activating..."
    # shellcheck disable=SC1090
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${ENV}"
    return 0
  fi

  ensure_channels
  local SOLVER="mamba"
  command -v mamba >/dev/null 2>&1 || SOLVER="conda"

  echo ">>> Creating env ${ENV} with: fastqc multiqc fastp seqkit pigz samtools"
  set +e
  ${SOLVER} create -n "${ENV}" -y -c conda-forge -c bioconda \
    python=3.11 \
    fastqc=0.12.1 \
    multiqc=1.21 \
    fastp \
    seqkit \
    pigz \
    samtools \
    rseqc \
    bedtools \
    pysam
  st=$?
  set -e

  if [ $st -ne 0 ]; then
    echo "!!! Strict solve failed. Retrying with flexible priority..."
    conda config --set channel_priority flexible
    ${SOLVER} create -n "${ENV}" -y -c conda-forge -c bioconda \
      python=3.11 \
      fastqc=0.12.1 \
      multiqc=1.21 \
      fastp \
      seqkit \
      pigz \
      samtools
    conda config --set channel_priority strict || true
  fi

  # shellcheck disable=SC1090
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "${ENV}"
}

export_env() {
  mkdir -p envs
  conda env export --name libsQC_illumina > envs/libsQC_illumina.yml
  echo ">>> Exported env -> envs/libsQC_illumina.yml"
}

gather_fastqs() {
  shopt -s nullglob
  FASTQS=( "${FASTQ_DIR}"/*.fastq.gz "${FASTQ_DIR}"/*.fq.gz "${FASTQ_DIR}"/*.fastq "${FASTQ_DIR}"/*.fq )
  shopt -u nullglob
  if [ ${#FASTQS[@]} -eq 0 ]; then
    echo "!!! No FASTQs found in ${FASTQ_DIR}/"
    exit 1
  fi
  echo ">>> Found ${#FASTQS[@]} FASTQ files in ${FASTQ_DIR}/"
}

# Detect PE pairs by common patterns (R1/R2 or _1/_2). Others treated as SE.
detect_pairs() {
  R1=()
  R2=()
  SINGLES=()

  for f in "${FASTQS[@]}"; do
    bn="$(basename "$f")"

    if [[ "$bn" =~ (_R1[^0-9A-Za-z]|_R1\.|\.R1\.|_1[^0-9A-Za-z]|_1\.) ]]; then
      # candidate R1
      r2="$f"
      r2="${r2/_R1/_R2}"
      r2="${r2/_1./_2.}"
      r2="${r2/.R1./.R2.}"
      if [ -f "$r2" ]; then
        R1+=("$f"); R2+=("$r2")
      else
        SINGLES+=("$f")
      fi
    elif [[ "$bn" =~ (_R2[^0-9A-Za-z]|_R2\.|\.R2\.|_2[^0-9A-Za-z]|_2\.) ]]; then
      # skip; paired when we see R1
      :
    else
      SINGLES+=("$f")
    fi
  done

  echo ">>> PE pairs: ${#R1[@]} | SE/unpaired: ${#SINGLES[@]}"
}

build_manifest() {
  mkdir -p metadata
  out="metadata/fastq_meta.tsv"
  echo -e "type\tsample\tread\tfile" > "$out"

  if (( ${#SINGLES[@]} > 0 )); then
    for f in "${SINGLES[@]}"; do
      [[ -n "$f" ]] || continue
      bn="$(basename "$f")"
      sample="${bn%.fastq.gz}"; sample="${sample%.fq.gz}"; sample="${sample%.fastq}"; sample="${sample%.fq}"
      echo -e "SE\t${sample}\tNA\t${bn}" >> "$out"
    done
  fi

  if (( ${#R1[@]} > 0 )); then
    for i in "${!R1[@]}"; do
      f1="${R1[$i]}"; f2="${R2[$i]}"
      [[ -n "$f1" && -n "$f2" ]] || continue
      bn1="$(basename "$f1")"; bn2="$(basename "$f2")"
      sample="$bn1"
      sample="${sample%_R1*}"; sample="${sample%_1*}"; sample="${sample%.R1.*}"
      echo -e "PE\t${sample}\tR1\t${bn1}" >> "$out"
      echo -e "PE\t${sample}\tR2\t${bn2}" >> "$out"
    done
  fi

  echo ">>> Manifest -> ${out}"
}

run_fastqc_raw() {
  mkdir -p "${RESULTS}/qc_raw"
  echo ">>> FastQC raw..."
  fastqc -t "${THREADS}" -o "${RESULTS}/qc_raw" "${FASTQS[@]}"
}

run_multiqc_raw() {
  mkdir -p "${RESULTS}/multiqc_raw"
  multiqc -o "${RESULTS}/multiqc_raw" "${RESULTS}/qc_raw"
  echo ">>> MultiQC raw -> ${RESULTS}/multiqc_raw/multiqc_report.html"
}

run_fastp() {
  mkdir -p "${RESULTS}/trimmed" "${RESULTS}/fastp"
  echo ">>> fastp trimming..."
  echo ">>> DEBUG: #R1=${#R1[@]} #R2=${#R2[@]} #SINGLES=${#SINGLES[@]}"

  # SE/unpaired
if (( ${#SINGLES[@]} > 0 )); then
  for f in "${SINGLES[@]:-}"; do
    bn="$(basename "$f")"
    base="${bn%.fastq.gz}"; base="${base%.fq.gz}"; base="${base%.fastq}"; base="${base%.fq}"
    outF="${RESULTS}/trimmed/${base}.trimmed.fastq.gz"
    html="${RESULTS}/fastp/${base}.fastp.html"
    json="${RESULTS}/fastp/${base}.fastp.json"

    fastp \
      -w "${THREADS}" \
      -i "$f" \
      -o "$outF" \
      -q "${FASTP_QUAL}" \
      -l "${FASTP_LEN_MIN}" \
      $( [ "${FASTP_TRIM_POLYG}" = "1" ] && echo "--trim_poly_g" ) \
      --html "${html}" --json "${json}" \
      1>"${RESULTS}/fastp/${base}.stdout.log" \
      2>"${RESULTS}/fastp/${base}.stderr.log"
  done
fi

  # PE
if (( ${#R1[@]} > 0 )); then
  for i in "${!R1[@]}"; do
    f1="${R1[$i]}"; f2="${R2[$i]}"
    bn1="$(basename "$f1")"
    base="${bn1%.fastq.gz}"; base="${base%.fq.gz}"; base="${base%.fastq}"; base="${base%.fq}"
    base="${base%_R1*}"; base="${base%_1*}"; base="${base%.R1.*}"

    out1="${RESULTS}/trimmed/${base}_R1.trimmed.fastq.gz"
    out2="${RESULTS}/trimmed/${base}_R2.trimmed.fastq.gz"
    html="${RESULTS}/fastp/${base}.fastp.html"
    json="${RESULTS}/fastp/${base}.fastp.json"

    fastp \
      -w "${THREADS}" \
      --detect_adapter_for_pe \
      $( [ "${FASTP_CORRECTION}" = "1" ] && echo "--correction" ) \
      $( [ "${FASTP_TRIM_POLYG}" = "1" ] && echo "--trim_poly_g" ) \
      -i "$f1" -I "$f2" \
      -o "$out1" -O "$out2" \
      -q "${FASTP_QUAL}" \
      -l "${FASTP_LEN_MIN}" \
      --html "${html}" --json "${json}" \
      1>"${RESULTS}/fastp/${base}.stdout.log" \
      2>"${RESULTS}/fastp/${base}.stderr.log"
  done
fi

  echo ">>> fastp outputs -> ${RESULTS}/trimmed/ and ${RESULTS}/fastp/"
}

run_fastqc_trimmed() {
  mkdir -p "${RESULTS}/qc_trimmed"
  shopt -s nullglob
  TRIMMED=( "${RESULTS}/trimmed/"*.fastq.gz )
  shopt -u nullglob
  [ ${#TRIMMED[@]} -eq 0 ] && { echo "!!! No trimmed FASTQs found"; exit 1; }
  echo ">>> FastQC trimmed..."
  fastqc -t "${THREADS}" -o "${RESULTS}/qc_trimmed" "${TRIMMED[@]}"
}

run_multiqc_trimmed() {
  mkdir -p "${RESULTS}/multiqc_trimmed"
  multiqc -o "${RESULTS}/multiqc_trimmed" "${RESULTS}/qc_trimmed" "${RESULTS}/fastp"
  echo ">>> MultiQC trimmed -> ${RESULTS}/multiqc_trimmed/multiqc_report.html"
}

seqkit_stats() {
  mkdir -p "${RESULTS}/summary"
  echo ">>> seqkit stats (raw)..."
  seqkit stats -a -T "${FASTQS[@]}" > "${RESULTS}/summary/seqkit_stats_raw.tsv"

  shopt -s nullglob
  TRIMMED=( "${RESULTS}/trimmed/"*.fastq.gz )
  shopt -u nullglob
  if [ ${#TRIMMED[@]} -gt 0 ]; then
    echo ">>> seqkit stats (trimmed)..."
    seqkit stats -a -T "${TRIMMED[@]}" > "${RESULTS}/summary/seqkit_stats_trimmed.tsv"
  fi
}

# -------------------
# Main
# -------------------
create_env_libsQC_illumina
export_env

gather_fastqs
detect_pairs
build_manifest

SKIP_RAW_QC="${SKIP_RAW_QC:-0}"

if [[ "${SKIP_RAW_QC}" -ne 1 ]]; then
  run_fastqc_raw
  run_multiqc_raw
else
  echo ">>> SKIP_RAW_QC=1: skipping FastQC + MultiQC on raw reads"
fi

if [[ "${RAW_QC_ONLY}" -eq 1 ]]; then
  echo ">>> RAW_QC_ONLY=1: skipping fastp + trimmed QC + seqkit"
  echo ">>> DONE."
  echo "    Raw MultiQC : ${RESULTS}/multiqc_raw/multiqc_report.html"
  exit 0
fi

run_fastp
run_fastqc_trimmed
run_multiqc_trimmed
seqkit_stats

echo ">>> DONE."
echo "    Raw MultiQC     : ${RESULTS}/multiqc_raw/multiqc_report.html"
echo "    Trimmed MultiQC : ${RESULTS}/multiqc_trimmed/multiqc_report.html"

# -------------------------------------------------
# Optional: RSeQC checks (requires BAMs + BED12)
# -------------------------------------------------
RUN_RSEQC="${RUN_RSEQC:-0}"
RSEQC_BED="${RSEQC_BED:-}"
RSEQC_BAM_DIR="${RSEQC_BAM_DIR:-}"

if [[ "${RUN_RSEQC}" -eq 1 ]]; then
  if [[ -z "${RSEQC_BED}" || -z "${RSEQC_BAM_DIR}" ]]; then
    echo "!!! RUN_RSEQC=1 but missing RSEQC_BED or RSEQC_BAM_DIR"
    echo "    Provide: --rseqc --rseqc-bed <genes.bed12> --rseqc-bam-dir <bam_dir>"
    exit 2
  fi

  if [[ ! -d "${RSEQC_BAM_DIR}" ]]; then
    echo "!!! RSEQC_BAM_DIR does not exist: ${RSEQC_BAM_DIR}"
    exit 2
  fi

  mkdir -p "${RESULTS}/rseqc"

  shopt -s nullglob
  BAMS=( "${RSEQC_BAM_DIR}"/*.bam )
  shopt -u nullglob

  if [[ ${#BAMS[@]} -eq 0 ]]; then
    echo "!!! No BAMs found in ${RSEQC_BAM_DIR}/*.bam (skipping RSeQC)"
  else
    echo ">>> Running RSeQC on ${#BAMS[@]} BAMs..."
    for bam in "${BAMS[@]}"; do
      sample="$(basename "$bam" .bam)"
      bash workflow/run_rseqc.sh \
        --bam "$bam" \
        --bed "${RSEQC_BED}" \
        --outdir "${RESULTS}/rseqc/${sample}"
    done
  fi
fi
