# HLA Pipeline Manager 🧬

[![Bash](https://img.shields.io/badge/Language-Bash-blue.svg)](https://www.gnu.org/software/bash/)
[![Bioinformatics](https://img.shields.io/badge/Domain-Bioinformatics-green.svg)]()
[![Portfolio](https://img.shields.io/badge/Status-Portfolio_Project-purple.svg)]()

**End-to-End Orchestration for High-Throughput HLA Imputation Pipelines.**

> **Note:** This repository contains sanitized versions of scripts developed during my tenure at **NIHR BioResource**. They are presented here for **educational and portfolio purposes only** to demonstrate proficiency in pipeline orchestration, data deployment, and cohort selection. No real patient data or internal infrastructure paths are included.

The **HLA Pipeline Manager** is a suite of operational scripts designed to manage the lifecycle of large-scale HLA imputation tasks. It handles everything from verifying the integrity of thousands of sub-batches to safely deploying final datasets and querying them for specific immunogenetic profiles.

---

## 🌟 Use Case: Operational Excellence

Running imputation on 50,000+ samples involves generating terabytes of data across thousands of small files. A Data Manager needs tools to:
1.  **Verify**: Did all 5,000 sub-jobs finish successfully? (`check_hla_imputation_status.sh`)
2.  **Deploy**: How do we move new data to production without breaking existing access? (`deploy_hla_results.sh`)
3.  **Query**: How do we quickly find patients with a specific risk allele for a new study? (`find_hla_participants.sh`)

## 🚀 Key Features

*   **✅ Automated Verification**: Checks file existence, non-zero size, and internal content (e.g., "Are there actually HLA alleles in this file?") for thousands of batches in seconds.
*   **⚡ HPC Optimized**: Includes `execute_batch_imputation.sh`, a production-ready script for Slurm/SGE environments that handles module loading, legacy tool patching (SNP2HLA/Plink 1.9), and parallel execution.
*   **🛡️ Safe Deployment**: Implements a "Backup-First" strategy. Existing production data is automatically archived before new data is copied, ensuring zero data loss.
*   **🔍 Rapid Cohort Selection**: Uses optimized `awk` processing to filter massive dosage files for specific alleles (e.g., `HLA_DRB1*04:01`) against a target participant list.
*   **📉 Error Reporting**: Provides detailed, actionable error logs (e.g., "Batch 26, Sub-batch 12 failed: No R2 values found").

## 📂 Repository Structure

```text
.
├── check_hla_imputation_status.sh  # 🕵️ QC: Verifies batch integrity
├── deploy_hla_results.sh           # 🚀 Ops: Safely deploys to production
├── find_hla_participants.sh        # 🔍 Query: Finds allele carriers
├── execute_batch_imputation.sh     # ⚡ HPC: Production execution script
├── locate_array_mapping_file.sh    # 🗺️ Setup: Finds metadata mapping
└── README.md                       # 📖 Documentation
```

## 🛠️ Usage

### 1. Verify a Run
Check if your imputation run completed successfully:
```bash
./check_hla_imputation_status.sh
```
*Output:*
> Batch 26 Status: SUCCESS
> All 120 expected sub-batches appear valid.

### 2. Deploy Results
Move verified files to the release directory (with auto-backup):
```bash
./deploy_hla_results.sh
```

### 3. Find Participants
Find all participants in your study list who carry `HLA_DRB1_0401`:
```bash
./find_hla_participants.sh HLA_DRB1_0401 my_study_ids.txt
```

## 🤝 Contributing
Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for details.

---
*Developed to streamline immunogenetics operations and ensure data reliability.*
