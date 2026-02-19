#!/bin/bash

# ==============================================================================
# Script Name: check_hla_imputation_status.sh
# Author:      Ugur Tuna
# Context:     Developed during tenure at NIHR BioResource (Cambridge).
# Disclaimer:  Sanitized version for educational/portfolio use.
#
# Description: Verifies the integrity of HLA imputation outputs across multiple
#              batches. Focuses on file existence and content validation rather
#              than potentially flawed log summaries.
# Usage:       ./check_hla_imputation_status.sh
# ==============================================================================

# set -e # Allow checking all batches even if one fails

echo "======================================================="
echo "--- Checking HLA Imputation Status (File Verification) ---"
echo "======================================================="
echo "Current Time: $(date)"
echo ""

# --- Configuration ---
# Base directory for processed data (Sanitized)
BASE_DIR="${BIOBANK_DATA_ROOT:-/path/to/biobank/processed/data}"
BATCH_LIST=(22 25 26 27)
OVERALL_STATUS=0 # 0=All OK, 1=Some Failed
BATCH_SIZE=1000 # Used to estimate expected batch count

# --- Helper Function ---
check_file_exists_and_non_empty() {
    local filepath=$1
    local filename=$(basename "$filepath")
    if [ ! -s "$filepath" ]; then # -s checks if file exists and is not empty
        echo "  ERROR: File missing or empty: $filename"
        return 1
    fi
    return 0
}

# --- Loop Through Batches ---
for BATCH_ID in "${BATCH_LIST[@]}"; do
    echo "--- Checking Batch ${BATCH_ID} ---"
    WORK_DIR="${BASE_DIR}/batch_${BATCH_ID}_hla_revised"
    MHC_PREFIX="b${BATCH_ID}_mhc"
    BATCH_STATUS="SUCCESS" # Assume success initially
    ERRORS=()

    # 1. Check Working Directory and MHC File
    if [ ! -d "${WORK_DIR}" ]; then
        ERRORS+=("Working directory ${WORK_DIR} not found.")
        BATCH_STATUS="FAILURE (No Work Dir)"
    elif [ ! -f "${WORK_DIR}/${MHC_PREFIX}.fam" ]; then
         ERRORS+=("MHC sample file ${MHC_PREFIX}.fam not found in ${WORK_DIR}.")
         BATCH_STATUS="FAILURE (No MHC FAM)"
    else
         # 2. Determine Expected Sub-batch Count
         # We use a subshell to avoid changing the main script's directory
         TOTAL_SAMPLES=$(wc -l < "${WORK_DIR}/${MHC_PREFIX}.fam")
         EXPECTED_BATCH_COUNT=$(( (TOTAL_SAMPLES + BATCH_SIZE - 1) / BATCH_SIZE )) # Ceiling division
         
         echo "Working directory: ${WORK_DIR}"
         echo "Expecting ${EXPECTED_BATCH_COUNT} sub-batches based on ${TOTAL_SAMPLES} samples."

         # 3. Check Output Files in batches/ directory
         BATCHES_DIR="${WORK_DIR}/batches"
         if [ ! -d "$BATCHES_DIR" ]; then
            ERRORS+=("Directory 'batches/' not found.")
            BATCH_STATUS="FAILURE (No Batches Dir)"
         else
            # Find actual imputed dosage files
            # Use find for safer file counting
            ACTUAL_COUNT=$(find "$BATCHES_DIR" -maxdepth 1 -name "batch_*_imputed.dosage" | wc -l)
            
            echo "Found ${ACTUAL_COUNT} sub-batch dosage files."

            if [[ "$ACTUAL_COUNT" -eq 0 ]]; then
                ERRORS+=("No imputed dosage files found.")
                BATCH_STATUS="FAILURE (No Output Files)"
            elif [[ "$ACTUAL_COUNT" -ne "$EXPECTED_BATCH_COUNT" ]]; then
                ERRORS+=("Mismatch: Expected ${EXPECTED_BATCH_COUNT} sub-batches, but found ${ACTUAL_COUNT} dosage files.")
                BATCH_STATUS="PARTIAL FAILURE (Batch Count Mismatch)"
            fi

            # Check individual output files for each found batch
            FILES_OK_COUNT=0
            ALL_SUBBATCHES_CHECKED_OK=1
            
            # Iterate over found dosage files
            for dosage_file in "$BATCHES_DIR"/batch_*_imputed.dosage; do
                if [ ! -f "$dosage_file" ]; then continue; fi

                prefix=${dosage_file%.dosage}
                filename_prefix=$(basename "$prefix")
                SUB_BATCH_OK=1
                
                # Verify all required file types exist
                check_file_exists_and_non_empty "${prefix}.bed" || SUB_BATCH_OK=0
                check_file_exists_and_non_empty "${prefix}.bim" || SUB_BATCH_OK=0
                check_file_exists_and_non_empty "${prefix}.fam" || SUB_BATCH_OK=0
                check_file_exists_and_non_empty "${prefix}.dosage" || SUB_BATCH_OK=0
                check_file_exists_and_non_empty "${prefix}.bgl.r2" || SUB_BATCH_OK=0
                check_file_exists_and_non_empty "${prefix}.bgl.log" || SUB_BATCH_OK=0

                if [ "$SUB_BATCH_OK" -eq 1 ]; then
                    # Content Validation
                    if ! grep -q -m 1 '^6\s\+HLA_' "${prefix}.bim"; then
                        echo "  ERROR: No HLA alleles found in ${filename_prefix}.bim"
                        SUB_BATCH_OK=0
                    fi
                    if ! grep -q -m 1 '^HLA_' "${prefix}.bgl.r2"; then
                        echo "  ERROR: No HLA R2 values found in ${filename_prefix}.bgl.r2"
                        SUB_BATCH_OK=0
                    fi
                    # Check beagle log for completion message
                    if ! tail -n 5 "${prefix}.bgl.log" | grep -q -E "finished|End time"; then
                        echo "  WARNING: Beagle completion message not found at end of ${filename_prefix}.bgl.log"
                    fi
                fi

                if [ "$SUB_BATCH_OK" -eq 0 ]; then
                    ERRORS+=("Sub-batch ${filename_prefix} failed file checks.")
                    BATCH_STATUS=${BATCH_STATUS/SUCCESS/PARTIAL FAILURE}
                    ALL_SUBBATCHES_CHECKED_OK=0
                else
                    FILES_OK_COUNT=$((FILES_OK_COUNT + 1))
                fi
            done

            # Refine status based on file checks
            if [[ "$BATCH_STATUS" != "FAILURE"* && "$ALL_SUBBATCHES_CHECKED_OK" -eq 1 && "$ACTUAL_COUNT" -eq "$EXPECTED_BATCH_COUNT" ]]; then
                 BATCH_STATUS="SUCCESS"
                 echo "All ${FILES_OK_COUNT} expected sub-batches appear valid."
            elif [[ "$BATCH_STATUS" != "FAILURE"* && "$FILES_OK_COUNT" -gt 0 ]]; then
                 BATCH_STATUS="PARTIAL FAILURE"
            else
                 BATCH_STATUS="FAILURE"
            fi
        fi
    fi

    # 4. Report Final Status for Batch
    echo "Batch ${BATCH_ID} Status: ${BATCH_STATUS}"
    if [ "${#ERRORS[@]}" -gt 0 ]; then
        echo "  Detected Issues:"
        for error in "${ERRORS[@]}"; do
            echo "  - $error"
        done
        OVERALL_STATUS=1
    fi
    echo "---------------------------"
    echo ""

done

# --- Overall Summary ---
echo "======================================================="
if [ "$OVERALL_STATUS" -eq 0 ]; then
    echo "--- Overall Status: ALL BATCHES PASSED verification! ---"
else
    echo "--- Overall Status: ISSUES DETECTED. Review errors above. ---"
fi
echo "======================================================="

exit $OVERALL_STATUS
