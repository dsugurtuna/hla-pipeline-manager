#!/bin/bash
#
# Portfolio/Educational Purpose Only
# -----------------------------------------------------------------------------
# Script: locate_array_mapping_file.sh
# Description: A utility script to recursively search through raw Affymetrix
#              metadata directories to locate the correct AX-to-rs ID mapping
#              file. Essential for the initialization phase of the pipeline.
#

# --- Configuration ---
# The specific Array ID we know exists in the target dataset
KNOWN_AX_ID="AX-13216142"

# The root directory of the raw Affymetrix data
SEARCH_DIR="./data/raw/Affy-BIOAXIOMAX"

echo "################################################################"
echo "### Searching for the UKBBv1.0 (AX-to-rs) mapping file"
echo "################################################################"
echo "Target ID: ${KNOWN_AX_ID}"
echo "Search Path: ${SEARCH_DIR}"
echo "----------------------------------------------------------------"

# Check if search directory exists (Simulation)
if [ ! -d "${SEARCH_DIR}" ]; then
    echo "[Simulation] Search directory not found locally."
    echo "Creating mock directory and file for demonstration..."
    mkdir -p "${SEARCH_DIR}"
    echo "AX-13216142,rs123456" > "${SEARCH_DIR}/UKBB_Mapping_v1.csv"
fi

echo "--- Searching for '${KNOWN_AX_ID}' in all .csv, .txt, and .annot files ---"
echo "--- (This may take a moment)..."
echo

# Recursive grep to find the file containing the ID
# -r = recursive
# -w = whole word match
# -l = print filename only (cleaner output)
grep -r -w -l "$KNOWN_AX_ID" "$SEARCH_DIR" --include="*.csv" --include="*.txt" --include="*.annot"

echo
echo "################################################################"
echo "### Search complete."
echo "### The file(s) listed above contain the target ID and are likely"
echo "### the correct mapping files for the pipeline."
echo "################################################################"
