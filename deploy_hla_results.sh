#!/bin/bash

# ==============================================================================
# Script Name: deploy_hla_results.sh
# Author:      Ugur Tuna
# Context:     Developed during tenure at NIHR BioResource (Cambridge).
# Disclaimer:  Sanitized version for educational/portfolio use.
#
# Description: Safely deploys verified HLA imputation results to the production
#              directory. Includes automated backup of existing files and 
#              post-copy verification.
# Usage:       ./deploy_hla_results.sh
# ==============================================================================

set -e # Exit on error

echo "======================================================="
echo "--- Deploying Final HLA Output Files ---"
echo "======================================================="
echo "Current Time: $(date)"
echo ""

# --- Configuration ---
BASE_DIR="${BIOBANK_DATA_ROOT:-/path/to/biobank/processed/data}"
FINAL_DIR="${BIOBANK_FINAL_DIR:-/path/to/biobank/release/HLA}"
BATCH_LIST=(22 25 26 27 28)
BACKUP_SUFFIX=$(date +%Y%m%d_%H%M%S)
BACKUP_DIR="${FINAL_DIR}/BACKUP_${BACKUP_SUFFIX}"

# --- Create Backup Directory ---
echo "Creating backup directory: ${BACKUP_DIR}"
mkdir -p "${BACKUP_DIR}"

# --- Loop Through Batches ---
for BATCH_ID in "${BATCH_LIST[@]}"; do
    echo ""
    echo "--- Processing Batch ${BATCH_ID} ---"

    WORK_DIR="${BASE_DIR}/batch_${BATCH_ID}_hla_revised"
    FINAL_PREFIX="${FINAL_DIR}/UKBBAffy_SAX_b${BATCH_ID}-HLA"
    MERGED_PREFIX="${WORK_DIR}/b${BATCH_ID}_imputed_merged"

    # Check if new merged files exist
    if [ ! -f "${MERGED_PREFIX}.bed" ] || [ ! -f "${MERGED_PREFIX}.bim" ] || [ ! -f "${MERGED_PREFIX}.fam" ]; then
        echo "  ERROR: Newly merged files not found at ${MERGED_PREFIX}.* - Skipping replacement for Batch ${BATCH_ID}."
        continue
    fi

    # --- Step 1: Backup Existing Files ---
    echo "  Step 1: Backing up existing files..."
    # Use find to move only existing files, avoid errors if some don't exist
    find "${FINAL_DIR}" -maxdepth 1 -type f -name "UKBBAffy_SAX_b${BATCH_ID}-HLA.*" -exec mv {} "${BACKUP_DIR}/" \;
    echo "  Existing files moved to backup."

    # --- Step 2: Copy New Files ---
    echo "  Step 2: Copying new files to production..."
    cp "${MERGED_PREFIX}.bed" "${FINAL_PREFIX}.bed"
    cp "${MERGED_PREFIX}.bim" "${FINAL_PREFIX}.bim"
    cp "${MERGED_PREFIX}.fam" "${FINAL_PREFIX}.fam"

    # --- Step 3: Verify Copy ---
    echo "  Step 3: Verifying copied files..."
    if [ -f "${FINAL_PREFIX}.bed" ] && [ -f "${FINAL_PREFIX}.bim" ] && [ -f "${FINAL_PREFIX}.fam" ]; then
        echo "  OK: New BED/BIM/FAM files successfully copied."
        ls -l "${FINAL_PREFIX}.bed" "${FINAL_PREFIX}.bim" "${FINAL_PREFIX}.fam"
    else
        echo "  ERROR: Failed to verify all copied files for Batch ${BATCH_ID}!"
        exit 1
    fi
    echo "--- Replacement for Batch ${BATCH_ID} Finished ---"

done

echo ""
echo "======================================================="
echo "--- Deployment Complete ---"
echo "Old files backed up in: ${BACKUP_DIR}"
echo "New files copied to: ${FINAL_DIR}"
echo "======================================================="

exit 0
