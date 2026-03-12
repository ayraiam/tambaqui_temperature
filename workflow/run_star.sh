#!/usr/bin/env bash
set -euo pipefail

THREADS="${THREADS:-8}"
RESULTS="${RESULTS:-results}"
TRIM_DIR="${TRIM_DIR:-${RESULTS}/trimmed}"

GENOME_FA=""
GTF=""
STAR_INDEX=""

READ_LENGTH="${READ_LENGTH:-151}"
SJDB_OVERHANG="${SJDB_OVERHANG:-$((READ_LENGTH-1))}"

RUN_INDEX=0
RUN_MAP=0
RUN_COUNTS=0

STRANDNESS="${STRANDNESS:-0}"   # featureCounts -s: 0 unstranded, 1 stranded, 2 reverse

MAKE_BED12=0
BED12_OUT=""

usage() {
  cat <<EOF
Usage:
  bash workflow/run_star.sh [options]

Required:
  --genome-fa PATH
  --gtf PATH
  --star-index-dir PATH

Actions:
  --index                   Build STAR index
  --map                     Map trimmed FASTQs with STAR
  --counts                  Run featureCounts from existing STAR BAMs
  --make-bed12              Create BED12 from GTF for RSeQC

Optional:
  --threads INT             (default: \$THREADS)
  --results DIR             (default: results)
  --trim-dir DIR            (default: results/trimmed)
  --read-length INT         (default: 151)
  --sjdb-overhang INT       (default: read_length-1)
  --strandness 0|1|2        featureCounts strandedness (default: 0)
  --bed12-out PATH          Output BED12 path (default: reference/genes.bed12)

EOF
  exit 0
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --threads) THREADS="$2"; shift 2 ;;
    --results) RESULTS="$2"; shift 2 ;;
    --trim-dir) TRIM_DIR="$2"; shift 2 ;;
    --genome-fa) GENOME_FA="$2"; shift 2 ;;
    --gtf) GTF="$2"; shift 2 ;;
    --star-index-dir) STAR_INDEX="$2"; shift 2 ;;
    --read-length) READ_LENGTH="$2"; SJDB_OVERHANG=$((READ_LENGTH-1)); shift 2 ;;
    --sjdb-overhang) SJDB_OVERHANG="$2"; shift 2 ;;
    --strandness) STRANDNESS="$2"; shift 2 ;;
    --bed12-out) BED12_OUT="$2"; shift 2 ;;
    --index) RUN_INDEX=1; shift ;;
    --map) RUN_MAP=1; shift ;;
    --counts) RUN_COUNTS=1; shift ;;
    --make-bed12) MAKE_BED12=1; shift ;;
    -h|--help) usage ;;
    *) echo "Unknown argument: $1" >&2; usage ;;
  esac
done

if [[ -z "${GENOME_FA}" || -z "${GTF}" || -z "${STAR_INDEX}" ]]; then
  echo "ERROR: --genome-fa, --gtf, and --star-index-dir are required." >&2
  usage
fi

if [[ "${RUN_INDEX}" -eq 0 && "${RUN_MAP}" -eq 0 && "${RUN_COUNTS}" -eq 0 && "${MAKE_BED12}" -eq 0 ]]; then
  echo "ERROR: choose at least one action: --index and/or --map and/or --counts and/or --make-bed12" >&2
  usage
fi

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

if ! find_and_source_conda; then
  echo "!!! conda not found." >&2
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

create_env_STAR_map() {
  local ENV="STAR_map"

  if conda env list | grep -qE "^${ENV}\s"; then
    echo ">>> Env ${ENV} exists. Activating..."
    conda activate "${ENV}"
    return 0
  fi

  ensure_channels
  local SOLVER="mamba"
  command -v mamba >/dev/null 2>&1 || SOLVER="conda"

  echo ">>> Creating env ${ENV} (STAR mapping / counts / RSeQC)"
  set +e
  ${SOLVER} create -n "${ENV}" -y -c conda-forge -c bioconda \
    python=3.11 \
    star samtools subread rseqc bedtools pigz ucsc-gtftogenepred
  st=$?
  set -e

  if [[ $st -ne 0 ]]; then
    echo "!!! Strict solve failed. Retrying with flexible priority..."
    conda config --set channel_priority flexible
    ${SOLVER} create -n "${ENV}" -y -c conda-forge -c bioconda \
      python=3.11 \
      star samtools subread rseqc bedtools pigz ucsc-gtftogenepred
    conda config --set channel_priority strict || true
  fi

  conda activate "${ENV}"

  mkdir -p envs
  conda env export --name "${ENV}" > "envs/${ENV}.yml"
  echo ">>> Exported env -> envs/${ENV}.yml"
}

create_env_STAR_map

if [[ "${MAKE_BED12}" -eq 1 ]]; then
  [[ -n "${BED12_OUT}" ]] || BED12_OUT="reference/genes.bed12"
  mkdir -p "$(dirname "${BED12_OUT}")"

  echo ">>> Building BED12 from GTF (RSeQC input)"
  tmp_gp="$(mktemp).gp"

  gtfToGenePred "${GTF}" "${tmp_gp}"
  genePredToBed -ignoreGroupsWithoutExons "${tmp_gp}" "${BED12_OUT}" || true
  [[ -s "${BED12_OUT}" ]] || { echo "ERROR: BED12 conversion produced empty output: ${BED12_OUT}" >&2; exit 2; }

  rm -f "${tmp_gp}"
  echo ">>> BED12 written: ${BED12_OUT}"
fi

star_index_complete() {
  [[ -f "${STAR_INDEX}/Genome" ]] && \
  [[ -f "${STAR_INDEX}/SA" ]] && \
  [[ -f "${STAR_INDEX}/SAindex" ]] && \
  [[ -f "${STAR_INDEX}/chrName.txt" ]] && \
  [[ -f "${STAR_INDEX}/chrLength.txt" ]]
}

if [[ "${RUN_INDEX}" -eq 1 ]]; then
  mkdir -p "${STAR_INDEX}"

  if star_index_complete; then
    echo ">>> STAR index already present at ${STAR_INDEX} (skipping genomeGenerate)"
  else
    echo ">>> STAR genomeGenerate"
    STAR \
      --runThreadN "${THREADS}" \
      --runMode genomeGenerate \
      --genomeDir "${STAR_INDEX}" \
      --genomeFastaFiles "${GENOME_FA}" \
      --sjdbGTFfile "${GTF}" \
      --sjdbOverhang "${SJDB_OVERHANG}"
  fi
fi

STAR_OUT="${RESULTS}/star"
STAR_QC="${RESULTS}/star_qc"
COUNTS_DIR="${RESULTS}/counts"
mkdir -p "${STAR_OUT}" "${STAR_QC}" "${COUNTS_DIR}"

if [[ "${RUN_MAP}" -eq 1 ]]; then
  if ! star_index_complete; then
    echo "ERROR: STAR index not found/complete at ${STAR_INDEX}." >&2
    exit 2
  fi

  [[ -d "${TRIM_DIR}" ]] || { echo "ERROR: TRIM_DIR not found: ${TRIM_DIR}" >&2; exit 2; }

  shopt -s nullglob
  PE_R1=( "${TRIM_DIR}"/*_R1.trimmed.fastq.gz )
  SE_ALL=( "${TRIM_DIR}"/*.trimmed.fastq.gz )
  shopt -u nullglob

  SE=()
  for f in "${SE_ALL[@]:-}"; do
    bn="$(basename "$f")"
    [[ "$bn" == *_R1.trimmed.fastq.gz || "$bn" == *_R2.trimmed.fastq.gz ]] && continue
    SE+=( "$f" )
  done

  echo ">>> Mapping from TRIM_DIR: ${TRIM_DIR}"
  echo ">>> PE samples: ${#PE_R1[@]} | SE samples: ${#SE[@]}"

  if [[ ${#PE_R1[@]} -eq 0 && ${#SE[@]} -eq 0 ]]; then
    echo "ERROR: No trimmed FASTQs found in TRIM_DIR=${TRIM_DIR}" >&2
    exit 2
  fi

  STAR_COMMON=(
    --runThreadN "${THREADS}"
    --genomeDir "${STAR_INDEX}"
    --sjdbGTFfile "${GTF}"
    --readFilesCommand zcat
    --outSAMtype BAM SortedByCoordinate
    --outSAMattributes NH HI AS nM MD
    --quantMode GeneCounts
  )

  run_one() {
    local sample="$1"
    shift

    local outdir="${STAR_OUT}/${sample}"
    mkdir -p "${outdir}"

    local bam="${outdir}/Aligned.sortedByCoord.out.bam"
    if [[ -s "${bam}" ]]; then
      echo ">>> Skipping ${sample}: BAM already exists"
      return 0
    fi

    echo ">>> STAR: ${sample}"
    STAR "${STAR_COMMON[@]}" --outFileNamePrefix "${outdir}/" "$@"

    samtools index -@ "${THREADS}" "${bam}"

    local qc="${STAR_QC}/${sample}"
    mkdir -p "${qc}"
    samtools flagstat -@ "${THREADS}" "${bam}" > "${qc}/flagstat.txt"
    samtools stats    -@ "${THREADS}" "${bam}" > "${qc}/stats.txt"
    samtools idxstats "${bam}" > "${qc}/idxstats.txt"
  }

  if (( ${#PE_R1[@]} > 0 )); then
    for r1 in "${PE_R1[@]}"; do
      bn="$(basename "$r1")"
      sample="${bn%_R1.trimmed.fastq.gz}"
      r2="${TRIM_DIR}/${sample}_R2.trimmed.fastq.gz"
      [[ -f "${r2}" ]] || { echo "ERROR: missing R2 for ${sample}: ${r2}" >&2; exit 2; }
      run_one "${sample}" --readFilesIn "${r1}" "${r2}"
    done
  fi

  if (( ${#SE[@]} > 0 )); then
    for f in "${SE[@]}"; do
      bn="$(basename "$f")"
      sample="${bn%.trimmed.fastq.gz}"
      run_one "${sample}" --readFilesIn "${f}"
    done
  fi

  echo ">>> Mapping finished."
fi

if [[ "${RUN_COUNTS}" -eq 1 ]]; then
  shopt -s nullglob
  BAMS=( "${STAR_OUT}"/*/Aligned.sortedByCoord.out.bam )
  shopt -u nullglob
  [[ ${#BAMS[@]} -gt 0 ]] || { echo "ERROR: no BAMs found for featureCounts in ${STAR_OUT}" >&2; exit 2; }

  PE_COUNT=0
  SE_COUNT=0
  for bam in "${BAMS[@]}"; do
    sample="$(basename "$(dirname "$bam")")"
    if [[ -f "${TRIM_DIR}/${sample}_R1.trimmed.fastq.gz" && -f "${TRIM_DIR}/${sample}_R2.trimmed.fastq.gz" ]]; then
      ((PE_COUNT+=1))
    else
      ((SE_COUNT+=1))
    fi
  done

  if [[ "${PE_COUNT}" -gt 0 && "${SE_COUNT}" -gt 0 ]]; then
    echo "ERROR: mixed PE and SE libraries detected. Run featureCounts separately for each type." >&2
    exit 2
  fi

  if [[ "${PE_COUNT}" -gt 0 ]]; then
    echo ">>> featureCounts (PE, gene-level): -s ${STRANDNESS} -p --countReadPairs"
    featureCounts \
      -T "${THREADS}" \
      -a "${GTF}" \
      -o "${COUNTS_DIR}/featureCounts.tsv" \
      -g gene_id -t exon \
      -s "${STRANDNESS}" \
      -p --countReadPairs \
      "${BAMS[@]}" \
      1>"${COUNTS_DIR}/featureCounts.stdout.log" \
      2>"${COUNTS_DIR}/featureCounts.stderr.log"
  else
    echo ">>> featureCounts (SE, gene-level): -s ${STRANDNESS}"
    featureCounts \
      -T "${THREADS}" \
      -a "${GTF}" \
      -o "${COUNTS_DIR}/featureCounts.tsv" \
      -g gene_id -t exon \
      -s "${STRANDNESS}" \
      "${BAMS[@]}" \
      1>"${COUNTS_DIR}/featureCounts.stdout.log" \
      2>"${COUNTS_DIR}/featureCounts.stderr.log"
  fi

  echo ">>> featureCounts finished."
  echo "    Output: ${COUNTS_DIR}/featureCounts.tsv"
fi

echo ">>> DONE."
echo "    STAR BAMs : ${STAR_OUT}/<sample>/Aligned.sortedByCoord.out.bam"
echo "    STAR logs : ${STAR_OUT}/<sample>/Log.final.out"
