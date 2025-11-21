#!/bin/bash
#
# Portfolio/Educational Purpose Only
# -----------------------------------------------------------------------------
# This script is part of a bioinformatics portfolio demonstrating technical
# competencies in genomic data processing and pipeline engineering.
#
# It contains sanitized code derived from production workflows. All internal
# paths, keys, and proprietary data have been removed or replaced with
# generic placeholders.
#
# Disclaimer: This code is for demonstration purposes and is not intended
# for clinical use without validation.
# -----------------------------------------------------------------------------
#
# Script: run_hla_imputation_pipeline.sh
# Description: Orchestrates the end-to-end HLA imputation workflow.
#              Includes chromosome filtering, ID mapping, MHC extraction,
#              and parallelized imputation using SNP2HLA/Beagle.
#
# Usage: ./run_hla_imputation_pipeline.sh -b <batch_id> -i <input_dir> -o <output_dir>
#

set -e  # Exit immediately if a command exits with a non-zero status.
set -o pipefail # Return value of a pipeline is the status of the last command to exit with a non-zero status.

# --- Configuration Defaults ---
BATCH_ID="26"
INPUT_DIR="./data/raw"
OUTPUT_DIR="./output"
REF_PANEL="/path/to/reference/HM_CEU_REF"
SNP2HLA_DIR="/path/to/SNP2HLA"
PLINK_EXEC="plink"
JAVA_MEM="20000" # 20GB
THREADS=4

# --- Logging Helper ---
log() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] [INFO] $1"
}

error() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] [ERROR] $1" >&2
    exit 1
}

# --- Argument Parsing ---
while getopts "b:i:o:r:m:t:h" opt; do
  case $opt in
    b) BATCH_ID="$OPTARG" ;;
    i) INPUT_DIR="$OPTARG" ;;
    o) OUTPUT_DIR="$OPTARG" ;;
    r) REF_PANEL="$OPTARG" ;;
    m) JAVA_MEM="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    h) 
       echo "Usage: $0 -b <batch_id> -i <input_dir> -o <output_dir> [-r <ref_panel>] [-m <java_mem>] [-t <threads>]"
       exit 0
       ;;
    \?) echo "Invalid option -$OPTARG" >&2; exit 1 ;;
  esac
done

# --- Derived Paths ---
BATCH_NAME="b${BATCH_ID}"
WORKDIR="${OUTPUT_DIR}/${BATCH_NAME}_hla_imputation"
FILTERED_PREFIX="${WORKDIR}/${BATCH_NAME}_filtered"
RENAMED_PREFIX="${WORKDIR}/${BATCH_NAME}_filtered_renamed"
MHC_PREFIX="${WORKDIR}/${BATCH_NAME}_mhc"
AX_RS_MAP_FILE="${INPUT_DIR}/AX_rs_mapping.tsv"

# --- Main Execution ---

log "Starting HLA Imputation Pipeline for Batch: ${BATCH_ID}"
log "Input Directory: ${INPUT_DIR}"
log "Output Directory: ${WORKDIR}"

# 1. Environment Setup
log "Step 1: Environment Setup"
mkdir -p "${WORKDIR}/batches"
cd "${WORKDIR}" || error "Failed to access working directory."

# Check dependencies
command -v $PLINK_EXEC >/dev/null 2>&1 || error "PLINK executable '$PLINK_EXEC' not found."
if [ ! -d "$SNP2HLA_DIR" ]; then
    log "Warning: SNP2HLA directory not found at default location. Assuming tools are in path or local."
fi

# Clean previous run artifacts
rm -f "${FILTERED_PREFIX}".* "${RENAMED_PREFIX}".* "${MHC_PREFIX}".* plink.log

# 2. Filter Chromosomes
log "Step 2: Filtering Chromosomes (1-22, X, Y, MT)"
# Note: In a real scenario, ensure input files exist.
# $PLINK_EXEC --bfile "${INPUT_DIR}/UKBBAffy_SAX_${BATCH_NAME}" ...
# For portfolio demonstration, we simulate the command or assume inputs exist.

if [ -f "${INPUT_DIR}/UKBBAffy_SAX_${BATCH_NAME}.bed" ]; then
    $PLINK_EXEC --bfile "${INPUT_DIR}/UKBBAffy_SAX_${BATCH_NAME}" \
        --allow-extra-chr --chr 1-22,X,Y,MT \
        --make-bed --out "${FILTERED_PREFIX}" > /dev/null
else
    log "Simulation: Input files not found, skipping actual PLINK execution for portfolio demo."
    touch "${FILTERED_PREFIX}.bed" "${FILTERED_PREFIX}.bim" "${FILTERED_PREFIX}.fam"
fi

# 3. Rename IDs (AX -> rs)
log "Step 3: Renaming Marker IDs (AX -> rs)"
if [ -f "${AX_RS_MAP_FILE}" ] && [ -f "${FILTERED_PREFIX}.bim" ]; then
    awk -v mapfile="$AX_RS_MAP_FILE" '
        BEGIN{
            FS=OFS="\t";
            while ((getline < mapfile) > 0) { map[$1] = $2 }
        }
        {
            if ($2 in map) { $2 = map[$2]; }
            print $0
        }' "${FILTERED_PREFIX}.bim" > "${RENAMED_PREFIX}.bim"
    
    cp "${FILTERED_PREFIX}.bed" "${RENAMED_PREFIX}.bed"
    cp "${FILTERED_PREFIX}.fam" "${RENAMED_PREFIX}.fam"
else
    log "Simulation: Mapping file or input missing. Skipping rename step."
    touch "${RENAMED_PREFIX}.bed" "${RENAMED_PREFIX}.bim" "${RENAMED_PREFIX}.fam"
fi

# 4. Extract MHC Region
log "Step 4: Extracting MHC Region (Chr6: 26Mb - 34Mb)"
if [ -f "${RENAMED_PREFIX}.bed" ]; then
    $PLINK_EXEC --bfile "${RENAMED_PREFIX}" \
        --chr 6 --from-bp 26000000 --to-bp 34000000 \
        --make-bed --out "${MHC_PREFIX}" > /dev/null
else
    log "Simulation: Skipping MHC extraction."
    touch "${MHC_PREFIX}.bed" "${MHC_PREFIX}.bim" "${MHC_PREFIX}.fam"
fi

# 5. Prepare SNP2HLA
log "Step 5: Preparing SNP2HLA Environment"
# Copying necessary jars and scripts (Mocking this for portfolio)
# cp "${SNP2HLA_DIR}/SNP2HLA.csh" .
# sed -i 's/--nonfounder/--nonfounders/g' SNP2HLA.csh

# 6. Batch Processing
log "Step 6: Splitting and Running Imputation Batches"
if [ -f "${MHC_PREFIX}.fam" ]; then
    # Split into smaller chunks for parallel processing
    split -l 1000 "${MHC_PREFIX}.fam" batches/batch_
    
    BATCH_FILES=$(ls batches/batch_* 2>/dev/null || true)
    
    if [ -z "$BATCH_FILES" ]; then
        log "No batches generated (likely empty input in demo mode)."
    else
        for batch_fam in $BATCH_FILES; do
            sub_batch=$(basename "$batch_fam")
            log "Processing sub-batch: $sub_batch"
            
            # Extract sub-batch
            $PLINK_EXEC --bfile "${MHC_PREFIX}" --keep "$batch_fam" --make-bed --out "batches/${sub_batch}" >/dev/null
            
            # Run SNP2HLA (Mock command)
            # ./SNP2HLA.csh "batches/${sub_batch}" "${REF_PANEL}" "batches/${sub_batch}_imputed" "$PLINK_EXEC" "$JAVA_MEM"
            
            log "Finished sub-batch: $sub_batch"
        done
    fi
else
    log "Simulation: Skipping batch splitting."
fi

# 7. Summary
log "Pipeline Completed Successfully."
log "Results stored in: ${WORKDIR}/batches/"

exit 0
