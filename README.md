<!-- LOGO -->
<p align="center">
  <img src="assets/logo-horizontal.png" alt="AY:RΔ Logo" width="160"/>
</p>

<h1 align="center" style="font-weight: normal;">Tambaqui_temperature</h1>

<p align="center">
  <img src="https://img.shields.io/badge/Pipeline-RNA--seq-blue">
  <img src="https://img.shields.io/badge/Conda-ready-green">
  <img src="https://img.shields.io/badge/STAR-supported-orange">
  <img src="https://img.shields.io/badge/Reproducible-workflow-purple">
  <img src="https://img.shields.io/badge/RSeQC-supported-yellow">
  <img src="https://img.shields.io/badge/featureCounts-ready-red">
  <img src="https://img.shields.io/badge/Subread-featureCounts-red">
  <img src="https://img.shields.io/badge/DESeq2-integrated-red">
  <img src="https://img.shields.io/badge/Functional%20Enrichment-ORA%20%2B%20GSEA-blue">
  <img src="https://img.shields.io/badge/DIAMOND-orthology-mapping-orange">
  <img src="https://img.shields.io/badge/clusterProfiler-supported-green">
</p>

<p align="center">
  <code>data and discovery in flow</code><br/>
  <a href="mailto:ayrabioinf@gmail.com">ayrabioinf@gmail.com</a> · 
  <a href="https://www.linkedin.com/company/aryaiam">LinkedIn</a>
</p>

---

<pre>
ABOUT
-----
Tambaqui_temperature is a reproducible, Conda-native RNA-seq pipeline
for transcriptomic analysis of temperature effects during tambaqui development.

It provides a structured, re-entrant framework for:

  1) Raw read quality control (FastQC + MultiQC)
  2) Adapter/quality trimming using fastp
  3) Post-trimming QC validation
  4) Read statistics reporting (SeqKit)
  5) Automatic QC summary table generation
  6) Genome alignment using STAR
  7) Gene-level quantification using featureCounts
  8) Optional RNA-seq diagnostics via RSeQC
  9) RNA-seq library strandedness inference via RSeQC
  10) Gene body coverage diagnostics
  11) Automatic GTF → BED12 conversion for RSeQC compatibility
  12) STAR alignment QC parsing (Log.final.out)
  13) Automatic generation of mapping QC summary tables
  14) Visualization of alignment composition across samples
  15) FeatureCounts-based sample QC (library size + detected genes)
  16) Exploratory transcriptomic QC (PCA + sample-to-sample distance heatmap)
  17) Variance partitioning analysis (Age vs Temperature effects)
  18) Modular execution of QC and statistical steps via CLI flags
  19) Automatic dependency detection + installation inside Conda environments
  20) Full provenance logging
  21) Differential gene expression analysis using DESeq2
  22) Statistical-first reporting without arbitrary log2FC thresholds
  23) Functional enrichment analysis via ORA and GSEA
  24) Orthology mapping via DIAMOND against Danio rerio
  25) Automated gene ID harmonization for non-model organisms
  26) Pathway-level biological interpretation of transcriptomic results  
  27) Candidate gene prioritization table integrating DESeq2 effect size,
      normalized expression, and GSEA leading-edge support from both
      GO Biological Process and KEGG
  28) Separate candidate-gene mode for downstream heatmap/qPCR target selection
  29) Combined pathway-aware candidate scoring based on repeated
      GSEA GO + KEGG core-enrichment support
  
The pipeline auto-detects:
  • Paired-end reads (R1/R2, _1/_2 patterns)
  • Single-end reads
  • Mixed layouts within the same directory

All stages are explicitly separated:
  - Raw QC
  - Trimming
  - Trimmed QC
  - Summary statistics
  - Optional alignment diagnostics

Designed for:
  • Local execution
  • Server execution
  • Screen-based background execution
  • Reproducible academic workflows
  • AY:RΔ production pipelines
</pre>

---

<pre>
STRUCTURE
---------
   /workflow/
   runall.sh
   run_libsQC_illumina.sh
   run_star.sh
   run_rseqc.sh
   run_mapping_qc_var.sh        - STAR Log.final.out parser + QC plotting
   parse_star_log_final.py      - STAR log parser
   plot_star_mapping_qc.R       - Mapping QC stacked barplots
   plot_pca_distance_qc.R      - PCA + sample-to-sample distance heatmap
   plot_variance_partition.R   - Variance partitioning (Age + Temperature)
   make_qc_summary_table.py
   run_deseq2.sh
   run_deseq2_analysis.R
   run_enrichment.sh
   run_enrichment_analysis.R
   prepare_enrichment_inputs.py
   merge_diamond_hits.py

 /envs/                          - Auto-exported Conda environments
 /logs/                          - Timestamped stdout/stderr logs
 /metadata/
   fastq_meta.tsv                - Auto-generated FASTQ manifest
 /data/                          - Input FASTQ(.gz) files
 /reference/
   genome.fa
   annotation.gtf
   protein.faa
   ncbi_dataset.tsv
   danio_rerio_protein.faa
   danio_rerio_ncbi_dataset.tsv
   star_index/                   - STAR genome index
  /results/
   qc_raw/
   multiqc_raw/
   trimmed/
   fastp/
   qc_trimmed/
   multiqc_trimmed/
   summary/
   star/
   star_qc/
   counts/
   rseqc/
   FigMappingQCandVar/          - STAR mapping QC tables + plots
  deseq2/
    enrichment/
     <contrast>/
       01_prepare/
       02_diamond/
       03_enrichment/
       04_candidates/
           candidate_gene_table.tsv
           candidate_gene_table_summary.tsv
    

 README.md
 LICENSE
 CITATION.cff
</pre>

---

<pre>
DESIGN PRINCIPLES
-----------------
 - Automatic Conda environment creation (strict channel priority)
 - Fully reproducible environment export (envs/libsQC_illumina.yml)
 - Raw reads never modified
 - Deterministic trimming (fastp)
 - Clear separation between raw QC and trimmed QC
 - Automatic QC summary table generation
 - Automatic PE/SE detection
 - Transparent logging (timestamped invocation logs)
 - Screen-compatible execution
 - Optional RNA-seq strandedness inference and coverage diagnostics
 - Fail-fast behavior (set -euo pipefail)
 - Explicit environment variable control
 - Re-entrant safe execution
 - STAR alignment integrated after trimming
 - Automatic paired-end detection for mapping
 - Alignment QC via samtools
 - Gene quantification via featureCounts
 - Safe re-entry (existing BAMs skipped)
 - Statistical significance prioritized over arbitrary fold-change thresholds
 - Full result tables retained for downstream interpretation
 - Designed for non-model organisms lacking curated pathway annotation
 - Orthology-based functional inference via DIAMOND alignment
</pre>

---

<pre>
PIPELINE STAGES
---------------

Stage 1 — Raw QC
  fastqc → results/qc_raw/
  multiqc → results/multiqc_raw/

Stage 2 — Trimming (fastp)
  • Adapter detection (PE auto-detect)
  • Phred filtering
  • Length filtering
  • Optional polyG trimming (NovaSeq/NextSeq)
  • Optional overlap correction (PE)

Outputs:
  results/trimmed/
  results/fastp/

Stage 3 — Trimmed QC
  fastqc → results/qc_trimmed/
  multiqc → results/multiqc_trimmed/

Stage 4 — Read statistics
  seqkit stats:
    results/summary/seqkit_stats_raw.tsv
    results/summary/seqkit_stats_trimmed.tsv

Stage 5 — QC summary table
  A consolidated QC report is generated automatically by comparing
  raw and trimmed read statistics.
    results/summary/qc_summary_table.tsv
  The table contains per-sample metrics including:
    raw read counts
    trimmed read counts
    reads removed
    percent reads removed
    raw vs trimmed total bases
    percent bases removed
    average read length change
    Q20/Q30 improvement
    GC content change

Stage 6 — Genome Alignment (STAR)
  STAR alignment of trimmed reads.

Outputs:
  results/star/<sample>/Aligned.sortedByCoord.out.bam
  results/star/<sample>/Log.final.out

Alignment QC:
  results/star_qc/<sample>/
    flagstat.txt
    stats.txt
    idxstats.txt

Stage 7 — Gene Quantification
  featureCounts gene-level quantification.

Output:
  results/counts/featureCounts.tsv

featureCounts settings

The pipeline uses fragment-level counting for paired-end RNA-seq libraries.

Default configuration:

  -t exon
  -g gene_id
  -p --countReadPairs
  -s <user specified>

Where:

  -t exon           counts exon features
  -g gene_id        summarizes reads at gene level
  -p                enables paired-end mode
  --countReadPairs  counts fragments instead of individual mates
  -s                library strandedness (determined via RSeQC)

Example for reverse-stranded libraries:

  -s 2 -p --countReadPairs -t exon -g gene_id

Stage 8 — RNA-seq Diagnostics (RSeQC)

Optional post-alignment diagnostics including:

  • Library strandedness inference
      infer_experiment.py

  • Gene body coverage assessment
      geneBody_coverage.py

The pipeline can automatically convert GTF annotations
to BED12 format required by RSeQC using UCSC utilities.

Outputs:
  results/rseqc/<sample>/
      infer_experiment.txt
      infer_experiment.log
      geneBody_coverage.*

Strandedness inference is typically run before gene counting
to determine the correct featureCounts parameter:

  -s 0   unstranded
  -s 1   forward-stranded
  -s 2   reverse-stranded

Stage 9 — STAR Mapping QC Visualization

The pipeline can parse STAR alignment logs (Log.final.out)
to generate a standardized mapping quality report across samples.

Metrics extracted include:

  • Uniquely mapped reads
  • Reads mapped to multiple loci
  • Reads mapped to too many loci
  • Unmapped reads (too many mismatches)
  • Unmapped reads (too short)
  • Unmapped reads (other)
  • Chimeric reads

Outputs:

  results/FigMappingQCandVar/
      star_mapping_qc_summary.tsv
      star_mapping_qc_summary.csv
      star_mapping_qc_stacked_percent.pdf
      star_mapping_qc_stacked_percent.png
      star_mapping_qc_stacked_counts.pdf
      star_mapping_qc_stacked_counts.png

Visualization:

A stacked barplot is generated showing the composition
of alignment outcomes for each sample, allowing rapid
identification of mapping issues or outlier libraries.

This stage uses:

  Python (pandas) for log parsing
  R (ggplot2) for visualization

and runs in a dedicated Conda environment:

  mappingqc_var_env

Stage 10 — FeatureCounts Sample QC + Exploratory Analysis

This stage extends QC beyond alignment metrics by leveraging
gene-level counts to assess sample quality and biological structure.

Inputs:
  results/counts/featureCounts.tsv
  metadata/sample_metadata.tsv

This stage generates:

1) Sample-level QC metrics:
   • Total assigned reads per sample (library size)
   • Number of detected genes (counts > 0)

Visualization:
   • Violin + boxplot + beeswarm plots
   • One point per sample
   • Custom color palette for individual samples

Outputs:
  results/FigMappingQCandVar/
      featurecounts_sample_qc.tsv
      featurecounts_library_size_violin_box_beeswarm.pdf/png
      featurecounts_detected_genes_violin_box_beeswarm.pdf/png

2) Principal Component Analysis (PCA):

Two PCA projections are generated:

  • Global PCA (all samples: C0 + T1 + T2)
  • Subset PCA (T1 vs T2 only)

These allow disentangling:

  • Developmental effects (C0 vs others)
  • Temperature effects (T1 vs T2)

Outputs:
  results/FigMappingQCandVar/
      pca_global_all_samples.pdf/png
      pca_global_all_samples_scores.tsv
      pca_t1_vs_t2.pdf/png
      pca_t1_vs_t2_scores.tsv

3) Sample-to-sample distance heatmap:

Pairwise distances are computed from variance-stabilized expression values.

Visualization:
  • Hierarchical clustering
  • Annotation by condition (and optional metadata)
  • Heatmap of sample similarity

Outputs:
  results/FigMappingQCandVar/
      sample_distance_matrix_all_samples.tsv
      sample_distance_heatmap_all_samples.pdf/png

Methodological notes:

  • Counts are transformed using DESeq2 variance stabilizing transformation (VST)
  • PCA is computed on VST-transformed data
  • Distance matrix uses Euclidean distance on VST values

This stage provides a critical bridge between:

  alignment QC → biological interpretation

and is recommended before differential expression analysis.

Stage 11 — Variance Partitioning Analysis

This stage quantifies the contribution of biological factors to global
gene expression variance across all samples.

Inputs:
  results/counts/featureCounts.tsv
  metadata/sample_metadata.tsv

Model:

  ~ Age + Temperature

Where:

  • Age reflects developmental stage (1 month vs 2 months)
  • Temperature reflects experimental exposure (pre, room, 37C)

Important design note:

  Temperature effects are evaluated within the same developmental stage
  (T1 vs T2), as temperature variation is only applied at 2 months.

Method:

  • Counts are normalized using edgeR
  • Lowly expressed genes are filtered (filterByExpr)
  • Expression values are transformed using voom (limma)
  • Variance is decomposed using variancePartition

Outputs:

  results/FigMappingQCandVar/
      variance_partition_gene_level.tsv
      variance_partition_summary.tsv
      variance_partition_violin.pdf/png

Outputs description:

1) Gene-level variance table:
   Fraction of variance explained per gene for:
     • Age
     • Temperature
     • Residual

2) Summary table:
   Median variance explained across all genes per factor

3) Visualization:
   • Violin plot of variance explained distributions
   • Boxplot overlay for median and dispersion

Interpretation:

  • Age captures global developmental effects
  • Temperature captures treatment-specific effects among same-age larvae
  • Residual captures unexplained or gene-specific variability

This stage provides a quantitative complement to PCA and clustering,
allowing formal decomposition of transcriptomic variance drivers.

Recommended usage:

  Run after exploratory QC and before differential expression analysis.

Stage 12 — Differential Gene Expression (DESeq2)

This stage performs statistical analysis of gene expression differences
between experimental conditions using DESeq2.

Inputs:
  results/counts/featureCounts.tsv
  metadata/sample_metadata.tsv

Key features:

  • Gene-level filtering based on expression:
      - Minimum count threshold (default: 10)
      - Minimum number of samples (default: 3)

  • Flexible experimental design via formula interface:
      ~ Condition

  • Subsetting of samples for pairwise comparisons:
      (e.g., T1 vs T2, C0 vs T1)

  • Proper factor releveling for contrast specification

  • Wald test for differential expression
  
  • Optional apeglm shrinkage for more stable log2 fold-change estimates

  • Multiple testing correction:
      Benjamini–Hochberg (FDR)

Outputs:

  results/deseq2/<comparison>/

      *_all_genes.tsv
          Full results table including:
            baseMean
            log2FoldChange
            lfcSE
            stat
            pvalue
            padj

      *_significant_genes.tsv
          Genes with:
            padj ≤ alpha (default: 0.05)

      *_ma_plot.pdf/png
      *_volcano_plot.pdf/png

Important methodological note:

  No arbitrary log2 fold-change cutoff is applied.

  All genes are retained in the output tables, allowing:

    • unbiased statistical interpretation
    • downstream effect size filtering if needed
    • full transparency of results

Recommended contrasts:

  • T1 vs C0 → developmental effect
  • T2 vs T1 → temperature effect (same age)
  • T2 vs C0 → combined effect

This stage represents the primary statistical analysis
of the transcriptomic response.

Stage 13 — Functional Enrichment Analysis (ORA + GSEA)

This stage identifies biological pathways and cellular processes
associated with differential gene expression results.

Because tambaqui is not a model organism, gene identifiers are first
mapped to orthologous genes in Danio rerio (zebrafish) using DIAMOND.

Workflow:

  1) Extract protein sequences for DE genes
  2) Align tambaqui proteins against zebrafish proteome using DIAMOND
  3) Retrieve ortholog gene identifiers
  4) Perform pathway enrichment analysis

Two complementary enrichment strategies are used:

ORA (Over-representation analysis)
  Tests whether significant DE genes are enriched for biological pathways.

GSEA (Gene set enrichment analysis)
  Uses the ranked gene list to detect coordinated pathway-level shifts.

Enrichment databases currently used in the pipeline:

  • Gene Ontology (GO)
      - Biological Process (ORA + GSEA)

  • KEGG pathways
      - ORA + GSEA
    
Methodological features:

  • Uses full DESeq2 results table (no log2FC filtering)
  • Automatically handles ID conversion
  • Deduplicates gene mappings for robust ranking
  • Compatible with non-model organisms

Inputs:

  results/deseq2/<contrast>/<contrast>_all_genes.tsv

  reference/protein.faa
      tambaqui protein sequences

  reference/ncbi_dataset.tsv
      tambaqui gene annotation table

  reference/danio_rerio_protein.faa
      zebrafish reference proteome

  reference/danio_rerio_ncbi_dataset.tsv
      zebrafish gene annotation table

Outputs:

  results/enrichment/<contrast>/

  01_prepare/
      source_annotation_from_deseq.tsv
      query_proteins.faa

  02_diamond/
      diamond_best_hits.tsv
      deseq_annotated_with_danio_hits.tsv

  03_enrichment/

      ora_go_bp_up.tsv
      ora_go_bp_down.tsv
      ora_kegg_up.tsv
      ora_kegg_down.tsv
      gsea_go_bp.tsv
      gsea_kegg.tsv
      gsea_ranked_gene_list.tsv
      enrichment_summary.tsv

      ora_go_bp_up_dotplot.png
      ora_go_bp_down_dotplot.png
      ora_kegg_up_dotplot.png
      ora_kegg_down_dotplot.png
      gsea_go_bp_dotplot.png
      gsea_kegg_dotplot.png

Interpretation:

  ORA identifies pathways over-represented among significant genes.

  GSEA identifies pathways with consistent directional expression shifts.

Together, these approaches provide complementary biological insight.

Stage 14 — Candidate Gene Prioritization

This optional downstream stage builds a compact candidate-gene table
for interpretation, heatmap design, and qPCR target selection.

It integrates:

  • DESeq2 effect size
      - log2FoldChange
      - abs_log2FoldChange
      - padj
      - baseMean

  • Gene expression level
      - mean normalized expression across all samples
      - mean normalized expression in each comparison group

  • Functional support from GSEA
      - number of significant GO Biological Process terms in which a gene
        appears in the GSEA leading edge (core_enrichment)
      - number of significant KEGG pathways in which a gene appears
        in the GSEA leading edge (core_enrichment)

This stage is designed to identify genes that are:

  • strongly differentially expressed
  • sufficiently expressed for robust interpretation
  • repeatedly implicated in enriched biological programs
  • supported across complementary pathway annotation systems

Inputs:

  results/deseq2/<contrast>/<contrast>_all_genes.tsv
  results/deseq2/<contrast>/normalized_counts.tsv
  results/deseq2/<contrast>/metadata_used.tsv
  results/enrichment/<contrast>/03_enrichment/gsea_go_bp.tsv
  results/enrichment/<contrast>/03_enrichment/gsea_kegg.tsv

Outputs:

  results/enrichment/<contrast>/04_candidates/
      candidate_gene_table.tsv
      candidate_gene_table_summary.tsv

Main columns in candidate_gene_table.tsv:

  • Geneid
  • GeneName
  • log2FoldChange
  • abs_log2FoldChange
  • padj
  • baseMean
  • mean_norm_all
  • mean_norm_<group1>
  • mean_norm_<group2>
  • gsea_go_core_term_count
  • gsea_kegg_core_term_count
  • candidate_score

Meaning of the GSEA support columns:

  • gsea_go_core_term_count
      Number of significant GO BP GSEA terms in which the gene appears
      in the leading edge

  • gsea_kegg_core_term_count
      Number of significant KEGG GSEA pathways in which the gene appears
      in the leading edge

Candidate score:

  candidate_score = abs_log2FoldChange *
                    log10(mean_norm_all + 1) *
                    (1 + log2(gsea_go_core_term_count +
                              gsea_kegg_core_term_count + 1))

Interpretation of candidate_score:

  Genes receive higher scores when they show:
    • larger expression changes
    • stronger overall expression
    • repeated support across enriched GO and/or KEGG leading-edge sets

Ranking strategy:

  Genes are prioritized using:
    1) total GSEA core support (GO + KEGG)
    2) absolute fold change
    3) adjusted p-value
    4) expression level

Interpretation:

  Genes with high fold change, robust normalized expression,
  and repeated GSEA leading-edge support across GO and KEGG
  are strong candidates for:
    • heatmap visualization
    • concise biological storytelling
    • qPCR validation

Important note:

  GSEA leading-edge support is counted using the mapped zebrafish
  gene identifier column (default: target_Gene_ID). Therefore,
  candidate-mode joins assume that the GSEA core_enrichment field
  and the mapped ortholog ID column use the same identifier system.
</pre>

---

<pre>
FASTP DEFAULTS
--------------
FASTP_QUAL        = 20
FASTP_LEN_MIN     = 30
FASTP_TRIM_POLYG  = 1   (recommended for NovaSeq/NextSeq)
FASTP_CORRECTION  = 1   (PE overlap correction)

Override via:
  --fastp-qual INT
  --fastp-len-min INT

Or environment variables:
  FASTP_TRIM_POLYG=0|1
  FASTP_CORRECTION=0|1
</pre>

---

<pre>
STEP CONTROL
------------

General:
  --cpus INT
  --fastq-dir PATH
  --results DIR
  --screen
  --screen-name STR
  --qc-summary-only      

QC Control:
  --no-qc
  --raw-qc-only
  --skip-raw-qc

RSeQC diagnostics:
  --rseqc                 Run RSeQC diagnostics
  --rseqc-bam PATH        Single coordinate-sorted BAM
  --rseqc-bam-dir DIR     Directory containing BAM files
  --rseqc-bed PATH        BED12 gene model
  --rseqc-gtf PATH        GTF annotation (can auto-convert to BED12)
  --rseqc-make-bed12      Convert GTF → BED12 automatically
  --rseqc-bed-out PATH    Output BED12 path
  --rseqc-full            Run full diagnostics (infer_experiment + geneBody_coverage)
  --gene-body-coverage   Run gene body coverage only

Default behavior runs only strandedness inference.
  
STAR Mapping:
  --star                   Run STAR mapping
  --star-index             Build STAR genome index
  --genome-fa PATH         Genome FASTA
  --gtf PATH               Gene annotation
  --star-index-dir PATH    STAR genome index directory
  --read-length INT        Read length (default: 151)

Quantification:
  --featurecounts          Run featureCounts using existing STAR BAMs
  --strandness 0|1|2       featureCounts strandedness

Variance Partitioning:
  --variance-partition     Run variance partitioning (Age + Temperature model)

DESeq2:
  --deseq2                    Run differential expression analysis
  --deseq2-counts-tsv PATH   featureCounts output
  --deseq2-metadata-tsv PATH Metadata file
  --deseq2-outdir PATH       Output directory
  --deseq2-design STR        Design formula (default: ~ Condition)
  --deseq2-subset-column STR
  --deseq2-subset-values STR
  --deseq2-reference-variable STR
  --deseq2-reference-level STR
  --deseq2-contrast-variable STR
  --deseq2-contrast-numerator STR
  --deseq2-contrast-denominator STR
  --deseq2-alpha FLOAT
  --deseq2-min-count INT
  --deseq2-min-samples INT  

Functional enrichment:
  --enrich
  --enrich-mode STR
      all|prepare|diamond|merge|analysis|candidates

  --enrich-deseq-tsv PATH
      DESeq2 results table

  --enrich-source-metadata-tsv PATH
      ncbi_dataset.tsv for tambaqui

  --enrich-source-protein-faa PATH
      protein sequences for tambaqui

  --enrich-target-metadata-tsv PATH
      ncbi_dataset.tsv for Danio rerio

  --enrich-target-protein-faa PATH
      zebrafish proteome FASTA

  --enrich-target-dmnd PATH
      optional prebuilt DIAMOND database

  --enrich-normalized-counts-tsv PATH
      DESeq2 normalized_counts.tsv for candidate mode

  --enrich-metadata-tsv PATH
      DESeq2 metadata_used.tsv for candidate mode

  --enrich-gsea-go-tsv PATH
      GSEA GO BP table for candidate mode

  --enrich-gsea-kegg-tsv PATH
      GSEA KEGG table for candidate mode

  --enrich-sample-col STR
      sample column in metadata [default: sample]

  --enrich-group-col STR
      grouping column in metadata [default: Condition]

  --enrich-outdir PATH
      output directory
</pre>

---

<pre>
EXAMPLES
--------

# 1) Default (raw QC → trim → trimmed QC → stats)
bash workflow/runall.sh

# 2) Raw QC only
bash workflow/runall.sh --raw-qc-only

# 3) Skip raw QC (re-run trimming only)
bash workflow/runall.sh --skip-raw-qc

# 4) Increase stringency
bash workflow/runall.sh --fastp-qual 25 --fastp-len-min 50

# 5) Run inside screen
bash workflow/runall.sh --screen

Attach later:
screen -r libsQC_illumina

# 6) Run RNA-seq diagnostics after alignment
bash workflow/runall.sh \
  --rseqc \
  --rseqc-bed genes.bed12 \
  --rseqc-bam-dir aligned_bams/

# 7) Custom FASTQ directory
bash workflow/runall.sh --fastq-dir data_batch2 --results results_batch2

# 8) Build STAR genome index
bash workflow/runall.sh \
  --star-index \
  --genome-fa reference/genome.fa \
  --gtf reference/annotation.gtf \
  --star-index-dir reference/star_index

# 9) Run STAR alignment
bash workflow/runall.sh \
  --star \
  --genome-fa reference/genome.fa \
  --gtf reference/annotation.gtf \
  --star-index-dir reference/star_index

# 10) Run gene counting from existing BAM files
bash workflow/runall.sh \
  --featurecounts \
  --strandness 2 \
  --genome-fa reference/genome.fa \
  --gtf reference/annotation.gtf \
  --star-index-dir reference/star_index
  
# 11) Build QC summary table only
bash workflow/runall.sh --qc-summary-only

# 12) Determine RNA-seq library strandedness
# (recommended before running featureCounts)

bash workflow/runall.sh \
  --no-qc \
  --rseqc \
  --rseqc-bam-dir results/star \
  --rseqc-gtf reference/annotation.gtf \
  --rseqc-make-bed12 \
  --rseqc-bed-out reference/annotation.bed12

Outputs:
  results/rseqc/<sample>/infer_experiment.txt

Inspect the fractions reported by infer_experiment.py
to determine the correct featureCounts strandedness mode.

# 13) Full mapping QC + featureCounts QC + PCA + clustering

bash workflow/runall.sh \
  --no-qc \
  --mapping-qc-var \
  --mapqc-star-dir /tmp/lab18/results/star \
  --mapqc-counts-tsv /home/lab18/results/counts/featureCounts.tsv \
  --mapqc-metadata-tsv /home/lab18/metadata/tambaqui_metadata.tsv \
  --mapqc-outdir /home/lab18/results/FigMappingQCandVar

This command generates:

  • STAR mapping QC stacked barplots
  • FeatureCounts sample QC plots (assigned reads, detected genes)
  • Global PCA (C0, T1, T2)
  • PCA restricted to T1 vs T2
  • Sample-to-sample distance heatmap

This represents the recommended exploratory QC workflow
prior to downstream statistical analysis.

# 14) Run variance partitioning only

bash workflow/runall.sh \
  --no-qc \
  --variance-partition \
  --mapqc-counts-tsv results/counts/featureCounts.tsv \
  --mapqc-metadata-tsv metadata/tambaqui_metadata.tsv

This command generates:

  • Variance partition violin plot
  • Gene-level variance decomposition table
  • Summary of variance explained by:
      - Age (development)
      - Temperature (treatment)
      - Residual variation

This analysis quantifies the relative contribution of biological
factors to global transcriptomic variation.

# 15) Run differential expression (T2 vs T1)

bash workflow/runall.sh \
  --no-qc \
  --deseq2 \
  --deseq2-counts-tsv results/counts/featureCounts.tsv \
  --deseq2-metadata-tsv metadata/sample_metadata.tsv \
  --deseq2-outdir results/deseq2/T2_vs_T1 \
  --deseq2-design "~ Condition" \
  --deseq2-subset-column Condition \
  --deseq2-subset-values T1,T2 \
  --deseq2-reference-variable Condition \
  --deseq2-reference-level T1 \
  --deseq2-contrast-variable Condition \
  --deseq2-contrast-numerator T2 \
  --deseq2-contrast-denominator T1

# 16) Prepare enrichment inputs only

bash workflow/runall.sh \
  --no-qc \
  --enrich \
  --enrich-mode prepare \
  --enrich-deseq-tsv results/deseq2/T2_vs_T1/condition_T2_vs_T1_all_genes.tsv \
  --enrich-source-metadata-tsv reference/ncbi_dataset.tsv \
  --enrich-source-protein-faa reference/protein.faa \
  --enrich-outdir results/enrichment/T2_vs_T1


# 17) Run full functional enrichment workflow

bash workflow/runall.sh \
  --no-qc \
  --enrich \
  --enrich-deseq-tsv results/deseq2/T2_vs_T1/condition_T2_vs_T1_all_genes.tsv \
  --enrich-source-metadata-tsv reference/ncbi_dataset.tsv \
  --enrich-source-protein-faa reference/protein.faa \
  --enrich-target-metadata-tsv reference/danio_rerio_ncbi_dataset.tsv \
  --enrich-target-protein-faa reference/danio_rerio_protein.faa \
  --enrich-outdir results/enrichment/T2_vs_T1

# 18) Build candidate gene table only
# (requires previous DESeq2 output + GSEA GO BP + GSEA KEGG outputs)

bash workflow/runall.sh \
  --no-qc \
  --enrich \
  --enrich-mode candidates \
  --enrich-deseq-tsv results/deseq2/T2_vs_T1/Condition_T2_vs_T1_all_genes.tsv \
  --enrich-normalized-counts-tsv results/deseq2/T2_vs_T1/normalized_counts.tsv \
  --enrich-metadata-tsv results/deseq2/T2_vs_T1/metadata_used.tsv \
  --enrich-gsea-go-tsv results/enrichment/T2_vs_T1/03_enrichment/gsea_go_bp.tsv \
  --enrich-gsea-kegg-tsv results/enrichment/T2_vs_T1/03_enrichment/gsea_kegg.tsv \
  --enrich-outdir results/enrichment/T2_vs_T1

# 19) Build candidate gene table using default inferred paths

bash workflow/runall.sh \
  --no-qc \
  --enrich \
  --enrich-mode candidates \
  --deseq2-outdir results/deseq2/T2_vs_T1 \
  --enrich-outdir results/enrichment/T2_vs_T1
</pre>

<pre>
DEPENDENCIES
------------
The pipeline automatically creates a dedicated Conda environment
for STAR alignment, featureCounts quantification, and RSeQC diagnostics.

Installed tools include:

  STAR
  subread (featureCounts)
  samtools
  rseqc
  ucsc-gtftogenepred
  ucsc-genepredtobed
  Python (pandas)
  R (ggplot2, dplyr, readr, tidyr)
Additional environments

mappingqc_var_env

Used for mapping QC, featureCounts QC, and exploratory
transcriptomic analysis (PCA and clustering).

Includes:

  Python
  pandas

  R
  ggplot2
  ggbeeswarm
  dplyr
  readr
  tidyr

  pheatmap
  RColorBrewer

  Bioconductor:
    DESeq2               (VST for PCA/clustering)
    limma                (voom transformation)
    edgeR                (normalization + filtering)
    variancePartition    (variance decomposition)  
    apeglm (optional for shrinkage)
    clusterProfiler
    enrichplot
    org.Dr.eg.db

Python:
  pandas
  biopython

Alignment:
  DIAMOND
  
If missing, the following packages will be installed automatically
into the active Conda environment:

  rseqc
  samtools
  ucsc-gtftogenepred
  ucsc-genepredtobed

This ensures that RSeQC diagnostics and BED12 conversion
work reproducibly across environments.
</pre>
---

<pre>
REPRODUCIBILITY
---------------
Each run logs:

  logs/invocation_YYYYMMDD_HHMMSS.txt
  logs/libsQC_YYYYMMDD_HHMMSS.out
  logs/libsQC_YYYYMMDD_HHMMSS.err

The Conda environment is exported to:

  envs/libsQC_illumina.yml

This ensures:
  • Exact tool versions
  • Exact dependency graph
  • Full reproducibility across machines
</pre>

---

<pre>
COMPATIBILITY
-------------
Tested with:
  • NovaSeq
  • NextSeq
  • Illumina Stranded Total RNA Prep (Ribo-Zero Plus)
  • Paired-end and single-end layouts
  • Mixed FASTQ directories
  • STAR ≥ 2.7
  • featureCounts (Subread)
  • samtools QC metrics
</pre>

---

<pre>
CITATION
--------
Lobo, I. (2026).
libsQC_illumina: A reproducible, Conda-native pipeline for
Illumina short-read quality control, trimming, and RNA-seq diagnostics.
AY:RΔ — data and discovery in flow.
</pre>

<p align="center"><sub>© 2026 AY:RΔ — data and discovery in flow</sub></p>
