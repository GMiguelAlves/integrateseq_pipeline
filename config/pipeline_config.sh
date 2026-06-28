#!/usr/bin/env bash
# Central configuration for the generic RNA-seq + ChIP-seq integration pipeline.
# Keep organism/project-specific paths here or in config/pipeline_config.local.sh.

set -euo pipefail

if [[ -n "${BASH_SOURCE[0]:-}" ]]; then
  export PROJECT_DIR="${PROJECT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
else
  export PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"
fi

export LOCAL_CONFIG="${LOCAL_CONFIG:-${PROJECT_DIR}/config/pipeline_config.local.sh}"
if [[ -f "${LOCAL_CONFIG}" ]]; then
  # shellcheck source=/dev/null
  source "${LOCAL_CONFIG}"
fi

# Source pipeline result directories. Override these in
# config/pipeline_config.local.sh or in the job environment.
export RNASEQ_RESULTS_DIR="${RNASEQ_RESULTS_DIR:-${PROJECT_DIR}/external/rnaseq_results}"
export CHIPSEQ_RESULTS_DIR="${CHIPSEQ_RESULTS_DIR:-${PROJECT_DIR}/external/chipseq_results}"
export INTEGRATION_OUTPUT_DIR="${INTEGRATION_OUTPUT_DIR:-${PROJECT_DIR}}"

# Organism-dependent files. Leave blank when unavailable.
export GENOME_FASTA="${GENOME_FASTA:-}"
export ANNOTATION_FILE="${ANNOTATION_FILE:-}"
export FUNCTIONAL_ANNOTATION="${FUNCTIONAL_ANNOTATION:-}"
export GENES_OF_INTEREST_FILE="${GENES_OF_INTEREST_FILE:-}"

# RNA-seq inputs.
export RNA_METADATA_FILE="${RNA_METADATA_FILE:-${RNASEQ_RESULTS_DIR}/025-parse/030-metadata_final/AllProjects_metadata_new.csv}"
export RNA_COUNTS_MATRIX="${RNA_COUNTS_MATRIX:-${RNASEQ_RESULTS_DIR}/050-quantification/counts_matrix.tsv}"
export RNA_NORMALIZED_MATRIX="${RNA_NORMALIZED_MATRIX:-${RNASEQ_RESULTS_DIR}/050-quantification/tpm_matrix.tsv}"
export RNA_DEG_RESULTS="${RNA_DEG_RESULTS:-${RNASEQ_RESULTS_DIR}/060-deg-analysis/DEGs_annotated_results.tsv}"
export RNA_GENE_CATALOG="${RNA_GENE_CATALOG:-${RNASEQ_RESULTS_DIR}/090-search-gene/results/tables/gene_catalog.tsv}"
export RNA_GENE_CATALOG_EXTRA="${RNA_GENE_CATALOG_EXTRA:-${RNASEQ_RESULTS_DIR}/../results/tables/gene_catalog.tsv}"
export RNA_EXPRESSION_CONTEXT="${RNA_EXPRESSION_CONTEXT:-${RNASEQ_RESULTS_DIR}/090-search-gene/results/tables/expression_summary_by_context.tsv}"
export RNA_WGCNA_HITS="${RNA_WGCNA_HITS:-${RNASEQ_RESULTS_DIR}/090-search-gene/results/tables/wgcna_hits.tsv}"
export RNA_MFUZZ_HITS="${RNA_MFUZZ_HITS:-${RNASEQ_RESULTS_DIR}/090-search-gene/results/tables/mfuzz_hits.tsv}"
export RNA_DTU_HITS="${RNA_DTU_HITS:-${RNASEQ_RESULTS_DIR}/090-search-gene/results/tables/dtu_hits.tsv}"
export RNA_SPLICING_HITS="${RNA_SPLICING_HITS:-${RNASEQ_RESULTS_DIR}/090-search-gene/results/tables/splicing_hits.tsv}"

# ChIP-seq inputs. Globs are allowed for peak tables/counts.
if [[ -z "${CHIP_METADATA_FILE:-}" ]]; then
  if [[ -s "${CHIPSEQ_RESULTS_DIR}/config/metadata.tsv" ]]; then
    export CHIP_METADATA_FILE="${CHIPSEQ_RESULTS_DIR}/config/metadata.tsv"
  else
    export CHIP_METADATA_FILE="${CHIPSEQ_RESULTS_DIR}/config/metadata_srp034587_no_input.tsv"
  fi
else
  export CHIP_METADATA_FILE
fi
export CHIP_ANNOTATED_PEAKS_GLOB="${CHIP_ANNOTATED_PEAKS_GLOB:-${CHIPSEQ_RESULTS_DIR}/090-peak-annotation/*.annotated.tsv*}"
export CHIP_PEAK_BED_GLOB="${CHIP_PEAK_BED_GLOB:-${CHIPSEQ_RESULTS_DIR}/110-consensus-peaks/groups/*.consensus.bed*}"
export CHIP_PEAK_COUNT_GLOB="${CHIP_PEAK_COUNT_GLOB:-${CHIPSEQ_RESULTS_DIR}/110-consensus-peaks/counts/*.counts.tsv*}"
export CHIP_DIFF_BINDING_FILE="${CHIP_DIFF_BINDING_FILE:-${CHIPSEQ_RESULTS_DIR}/120-differential-binding/differential_binding_results.tsv.gz}"
if [[ ! -s "${CHIP_DIFF_BINDING_FILE}" && -s "${CHIPSEQ_RESULTS_DIR}/120-differential-binding/differential_binding_results.tsv" ]]; then
  export CHIP_DIFF_BINDING_FILE="${CHIPSEQ_RESULTS_DIR}/120-differential-binding/differential_binding_results.tsv"
fi

# Column mapping and biological variables.
export GENE_ID_COLUMN="${GENE_ID_COLUMN:-gene_id}"
export GENE_NAME_COLUMN="${GENE_NAME_COLUMN:-gene_name}"
export TRANSCRIPT_ID_COLUMN="${TRANSCRIPT_ID_COLUMN:-transcript_id}"
export SAMPLE_ID_COLUMN="${SAMPLE_ID_COLUMN:-sample_id}"
export GROUP_COLUMNS="${GROUP_COLUMNS:-stage,sex,condition,treatment,tissue}"
export CONTRAST_ID_COLUMN="${CONTRAST_ID_COLUMN:-contrast_id}"
export MARK_COLUMN="${MARK_COLUMN:-mark_or_factor}"
export CONDITION_COLUMN="${CONDITION_COLUMN:-condition}"

# Integration thresholds.
export PEAK_GENE_WINDOW_BP="${PEAK_GENE_WINDOW_BP:-5000}"
export PROMOTER_UPSTREAM_BP="${PROMOTER_UPSTREAM_BP:-2000}"
export PROMOTER_DOWNSTREAM_BP="${PROMOTER_DOWNSTREAM_BP:-500}"
export DEG_PADJ_THRESHOLD="${DEG_PADJ_THRESHOLD:-0.05}"
export DEG_LOG2FC_THRESHOLD="${DEG_LOG2FC_THRESHOLD:-1}"
export DIFF_BINDING_PADJ_THRESHOLD="${DIFF_BINDING_PADJ_THRESHOLD:-0.05}"
export DIFF_BINDING_LOG2FC_THRESHOLD="${DIFF_BINDING_LOG2FC_THRESHOLD:-1}"
export TOP_CANDIDATES_N="${TOP_CANDIDATES_N:-100}"

# Runtime.
export RUN_MODE="${RUN_MODE:-slurm}"       # slurm or local
export THREADS="${THREADS:-8}"
export MEMORY="${MEMORY:-32G}"
export SLURM_TIME="${SLURM_TIME:-12:00:00}"
export SLURM_PARTITION="${SLURM_PARTITION:-compute}"
export SLURM_ACCOUNT="${SLURM_ACCOUNT:-}"
export PYTHON_BIN="${PYTHON_BIN:-python3}"
export RSCRIPT_BIN="${RSCRIPT_BIN:-Rscript}"
export ENV_BACKEND="${ENV_BACKEND:-none}"  # none, conda, apptainer, singularity
export CONDA_ENV="${CONDA_ENV:-integrateseq}"
export CONTAINER_IMAGE="${CONTAINER_IMAGE:-}"

# Safety.
export OVERWRITE="${OVERWRITE:-false}"
export CREATE_DONE_FILES="${CREATE_DONE_FILES:-true}"
