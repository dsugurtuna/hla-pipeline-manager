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
# Script: generate_clinical_hla_report.sh
# Description: Collates HLA imputation results from multiple batches into a
#              final clinical-grade CSV report. Demonstrates complex text
#              parsing (AWK) of dosage files (.dos) and integration of 
#              external data sources (e.g., CookHLA).
#
# Usage: ./generate_clinical_hla_report.sh
#

set -e
echo "--- Starting Final HLA Report Generation ---"

# --- Configuration ---
# In a production environment, these would be passed as arguments or env vars
BASE_PATH="./data/processed"
OUTPUT_DIR="./output/reports"
FINAL_REPORT="${OUTPUT_DIR}/Final_HLA_Report_DELIVERABLE.csv"

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

# --- Participant Data Definitions (Sanitized) ---
# Associative array mapping ParticipantID to metadata string.
# Format: [ParticipantID]="VariantID;BatchID;FamilyID_ArrayID"
# Note: Real IDs have been replaced with generic placeholders.
declare -A P_DATA
P_DATA=(
    ["P001_REF"]="V_REF_001;b18;Sample_Ref_01_Array_H03"
    ["P002_TEST"]="V_TEST_002;b06;Sample_Test_02_Array_B02"
    ["P003_TEST"]="V_TEST_003;b09;Sample_Test_03_Array_H07"
    ["P004_TEST"]="V_TEST_004;b11;Sample_Test_04_Array_E07"
    ["P005_TEST"]="V_TEST_005;b04;Sample_Test_05_Array_E02"
)

# --- Create final report with header ---
echo "V_ID,Participant_ID,Data_Source,Batch,HLA_A,HLA_B,HLA_C,HLA_DPA1,HLA_DPB1,HLA_DQA1,HLA_DQB1,HLA_DRB1" > "$FINAL_REPORT"

# --- 1. External Method Integration (Example: CookHLA) ---
# Demonstrates how to merge results from different calling algorithms (e.g., for validation)
echo "Processing External Reference Sample (P001_REF)..."

# Simulating the parsing of an external tool's output format
# In the original script, this parsed a specific .alleles file format
# COOKHLA_RESULTS=$(grep "$ARRAY_ID_B18" "$COOKHLA_FILE" || true)

# For portfolio demonstration, we inject the parsed result directly to show the data structure
# This represents a sample where we trust an external caller (CookHLA) over the pipeline
echo "V_REF_001,P001_REF,CookHLA_b18,b18,HLA-A*01:01/HLA-A*02:01,HLA-B*07:02/HLA-B*08:01,HLA-C*07:02/HLA-C*07:02,HLA-DPA1*01:03/HLA-DPA1*01:03,HLA-DPB1*04:01/HLA-DPB1*04:01,HLA-DQA1*01:02/HLA-DQA1*05:01,HLA-DQB1*02:01/HLA-DQB1*06:02,HLA-DRB1*03:01/HLA-DRB1*15:01" >> "$FINAL_REPORT"

# --- 2. Process All Imputed Batches (SNP2HLA) ---
for PID in "${!P_DATA[@]}"; do
    if [[ "$PID" == "P001_REF" ]]; then continue; fi # Skip the one handled above
    
    IFS=';' read -r VID BATCH FAM_ID <<< "${P_DATA[$PID]}"
    echo "Parsing results for $VID ($PID) from batch $BATCH..."
    
    # Construct paths dynamically based on batch ID
    WORK_DIR="${BASE_PATH}/imputation_batch_${BATCH}"
    
    R2_FILE="$WORK_DIR/${BATCH}_imputed.bgl.r2"
    DOS_FILE="$WORK_DIR/${BATCH}_imputed.dos"
    FAM_FILE="$WORK_DIR/${BATCH}_subset_mhc.fam"

    # Check if files exist (Mocking logic for portfolio if files are missing)
    if [ ! -f "$DOS_FILE" ] || [ ! -f "$R2_FILE" ]; then
        echo "   [Simulated] Input files not found locally. Generating mock entry for demonstration."
        # Generate a plausible mock entry to show what the output looks like
        echo "$VID,$PID,Imputed_SNP2HLA,$BATCH,HLA-A*02:01/HLA-A*03:01,HLA-B*07:02/HLA-B*44:02,HLA-C*05:01/HLA-C*07:02,HLA-DPA1*01:03/HLA-DPA1*02:01,HLA-DPB1*04:01/HLA-DPB1*04:02,HLA-DQA1*01:02/HLA-DQA1*05:05,HLA-DQB1*03:01/HLA-DQB1*06:02,HLA-DRB1*04:01/HLA-DRB1*11:01" >> "$FINAL_REPORT"
        continue
    fi

    # --- Core Parsing Logic ---
    # This AWK script is the heart of the report generation.
    # It reads three files simultaneously:
    # 1. R2 File: To check imputation quality (R2 > 0.3)
    # 2. FAM File: To find the column index of the specific individual
    # 3. DOS File: To read the dosage probabilities for that individual
    
    PARSED_RESULTS=$(awk -v target_sample="$FAM_ID" '
        # File 1: R2 Quality Metrics
        ARGIND==1 { 
            if ($1 ~ /^HLA_/) r2[$1] = $2; 
            next 
        }
        
        # File 2: FAM File (Sample Indexing)
        ARGIND==2 { 
            if ($2 == target_sample) { target_col = FNR + 3; } 
            next; 
        }
        
        # File 3: Dosage File (Genotype Probabilities)
        ARGIND==3 {
            marker = $1;
            # Filter: Must be HLA marker, valid R2, and R2 > 0.3 (QC threshold)
            if (marker ~ /^HLA_/ && r2[marker] != "NaN" && r2[marker] > 0.3) {
                split(marker, parts, "_");
                gene = parts[2];
                
                # Parse 4-digit resolution (e.g., 0101 -> 01:01)
                if (length(parts[3]) == 4 && parts[3] ~ /^[0-9]+$/) {
                    type = substr(parts[3],1,2) ":" substr(parts[3],3,4);
                    # Accumulate dosage for this specific allele type
                    dosages[gene, type] += $(target_col);
                }
            }
        }
        
        # Final Aggregation
        END {
            loci="A,B,C,DPA1,DPB1,DQA1,DQB1,DRB1";
            split(loci, la, ",");
            
            for (l_idx in la) {
                gene = la[l_idx];
                best1_type = ""; best1_dos = -1;
                best2_type = ""; best2_dos = -1;
                
                # Find the two alleles with the highest dosages for this gene
                for (type_dos_key in dosages) {
                    split(type_dos_key, parts, SUBSEP);
                    if (parts[1] == gene) {
                        current_dos = dosages[gene, parts[2]];
                        if (current_dos > best1_dos) {
                            best2_dos = best1_dos; best2_type = best1_type;
                            best1_dos = current_dos; best1_type = parts[2];
                        } else if (current_dos > best2_dos) {
                            best2_dos = current_dos; best2_type = parts[2];
                        }
                    }
                }
                
                # Determine Genotype Call based on dosages
                if (best1_type == "") { 
                    hla[gene] = "-/-";
                } else if (best1_dos > 1.5) { 
                    # Homozygous call (Dosage near 2.0)
                    hla[gene] = "HLA-" gene "*" best1_type "/HLA-" gene "*" best1_type;
                } else if (best2_type != "") { 
                    # Heterozygous call
                    hla[gene] = "HLA-" gene "*" best1_type "/HLA-" gene "*" best2_type;
                } else { 
                    # Ambiguous/Single call (Low confidence on second allele)
                    hla[gene] = "HLA-" gene "*" best1_type "/-"; 
                }
            }
            
            # Output formatted CSV string
            printf "%s,%s,%s,%s,%s,%s,%s,%s",
                hla["A"],hla["B"],hla["C"],hla["DPA1"],hla["DPB1"],
                hla["DQA1"],hla["DQB1"],hla["DRB1"];
        }
    ' "$R2_FILE" "$FAM_FILE" "$DOS_FILE")
    
    echo "$VID,$PID,Imputed_SNP2HLA,$BATCH,$PARSED_RESULTS" >> "$FINAL_REPORT"
done

# --- Final Notification ---
echo "============================================================"
echo "TASK COMPLETE. All participants processed."
echo "Final report is located at: $FINAL_REPORT"
echo "============================================================"
echo ""
# Display the first few lines of the report
head -n 5 "$FINAL_REPORT"
