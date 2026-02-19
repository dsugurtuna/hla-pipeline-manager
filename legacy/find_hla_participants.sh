#!/bin/bash

# ==============================================================================
# Script Name: find_hla_participants.sh
# Author:      Ugur Tuna
# Context:     Developed during tenure at NIHR BioResource (Cambridge).
# Disclaimer:  Sanitized version for educational/portfolio use.
#
# Description: Automates the identification of participants carrying a specific
#              HLA allele (e.g., HLA_DRB1_0401) from imputed dosage files.
# Usage:       ./find_hla_participants.sh <allele_name> <participant_list.txt>
# ==============================================================================

set -e

# --- Configuration ---
DOSAGE_DIR="${BIOBANK_DOSAGE_DIR:-./data/dosages}"
OUTPUT_FILE="hla_positive_participants.txt"

# --- Validation ---
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <allele_name> <participant_list.txt>"
    echo "Example: $0 HLA_DRB1_0401 ra_participants.txt"
    exit 1
fi

TARGET_ALLELE=$1
PARTICIPANT_LIST=$2

if [ ! -f "$PARTICIPANT_LIST" ]; then
    echo "Error: Participant list file '$PARTICIPANT_LIST' not found."
    exit 1
fi

echo "========================================================"
echo "HLA Participant Finder"
echo "Target Allele: $TARGET_ALLELE"
echo "Input List:    $PARTICIPANT_LIST"
echo "========================================================"

# --- Step 1: Locate Dosage Files ---
# Find the first dosage file to read the header
FIRST_DOSAGE_FILE=$(find "$DOSAGE_DIR" -name "*.dosage" | head -n 1)

if [ -z "$FIRST_DOSAGE_FILE" ]; then
    echo "Error: No dosage files found in $DOSAGE_DIR"
    exit 1
fi

echo "Reading header from: $(basename "$FIRST_DOSAGE_FILE")"

# --- Step 2: Find Column Index ---
# Automatically find the column number from the header line
COLUMN_NUMBER=$(head -n 1 "${FIRST_DOSAGE_FILE}" | awk -v allele="${TARGET_ALLELE}" -F'\t' '{for(i=1;i<=NF;i++) if($i==allele) {print i; exit}}')

if [[ -z "${COLUMN_NUMBER}" ]]; then
    echo "ERROR: Could not find allele '${TARGET_ALLELE}' in the header."
    echo "Please check the allele name spelling."
    exit 1
fi

echo "--> Allele '${TARGET_ALLELE}' found at column: ${COLUMN_NUMBER}"

# --- Step 3: Filter Participants ---
echo "Filtering dosage files..."

# This awk command loads the participant list into a hash map, then processes dosage files.
# It checks if a person is on the list AND has a dosage > 0 for the target allele.
awk -v COL="${COLUMN_NUMBER}" -F'\t' '
    BEGIN {
        while( (getline < "'${PARTICIPANT_LIST}'") > 0 ) {
            target_pids[$1] = 1
        }
        close("'${PARTICIPANT_LIST}'")
    }
    FNR > 1 {
        # Check if ID (column 1) is in our target list AND dosage > 0
        if ( ($1 in target_pids) && ($COL > 0) ) {
            print $1
        }
    }
' "$DOSAGE_DIR"/*.dosage > "${OUTPUT_FILE}"

# --- Final Report ---
FINAL_COUNT=$(wc -l < "${OUTPUT_FILE}")
echo "========================================================"
echo "Search Complete."
echo "Found ${FINAL_COUNT} participants matching criteria."
echo "Results saved to: ${OUTPUT_FILE}"
echo "========================================================"
