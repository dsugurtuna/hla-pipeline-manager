# HLA Pipeline Manager

[![CI](https://github.com/dsugurtuna/hla-pipeline-manager/actions/workflows/ci.yml/badge.svg)](https://github.com/dsugurtuna/hla-pipeline-manager/actions)
[![Python](https://img.shields.io/badge/Python-3.9%2B-blue)](https://www.python.org/)
[![Portfolio](https://img.shields.io/badge/Status-Portfolio_Project-purple.svg)]()

End-to-end HLA imputation pipeline orchestration — batch execution, verification, safe deployment, and clinical report generation.

> **Portfolio disclaimer:** This repository contains sanitised, generalised versions of tooling developed at NIHR BioResource. No real participant data or internal paths are included.

---

## Overview

Running HLA imputation on 50,000+ samples produces thousands of sub-batch output files. This toolkit manages the full lifecycle:

- **Batch executor** — chromosome filtering, AX→rs ID renaming, MHC-region extraction, sub-batch splitting, SNP2HLA orchestration.
- **Imputation verifier** — checks .bed/.bim/.fam/.dosage/.bgl.r2/.bgl.log presence, HLA marker counts, and Beagle completion.
- **Safe deployer** — backup-first strategy with automatic timestamped archiving before production copy.
- **Clinical reporter** — genotype calling from dosage data with homozygous/heterozygous/negative classification and CSV export.

## Repository Structure

```text
.
├── src/hla_pipeline/               Python package
│   ├── __init__.py
│   ├── executor.py                 Batch imputation executor
│   ├── verifier.py                 Output verification engine
│   ├── deployer.py                 Safe deployment with backup
│   └── reporter.py                 Clinical genotype reporter
├── tests/
│   ├── test_executor.py
│   ├── test_verifier.py
│   ├── test_deployer.py
│   └── test_reporter.py
├── legacy/                         Original shell scripts
├── .github/workflows/ci.yml
├── pyproject.toml
├── Dockerfile
└── Makefile
```

## Quick Start

```bash
pip install -e ".[dev]"
```

### Python API

```python
from hla_pipeline import BatchExecutor, ImputationVerifier, ResultDeployer, ClinicalReporter

# Execute batch imputation
executor = BatchExecutor()
result = executor.execute_batch("B001", "batch.fam", "/work/")
print(f"{result.sub_batches} sub-batches, {result.success_rate:.0%} success")

# Verify outputs
verifier = ImputationVerifier()
report = verifier.verify_batch("/output/batch_001/")
print(verifier.format_report(report))

# Deploy safely
deployer = ResultDeployer("/production/hla/")
report = deployer.deploy("/staging/batch_001/")
print(f"Deployed {report.deployment_count} files, backup at {report.backup_dir}")

# Generate clinical report
reporter = ClinicalReporter()
report = reporter.generate_report("dosage.raw", "HLA_DRB1_0101")
reporter.export_csv(report, "clinical_report.csv")
```

## Testing

```bash
make test   # or: pytest tests/ -v
```

## Jira Provenance

- **HLA imputation pipeline** — full SNP2HLA/CookHLA orchestration across Cambridge HPC.
- **Batch verification** — validating completeness of thousands of sub-batch outputs.
- **Clinical reporting** — genotype calling and deliverable CSV generation for clinical teams.

## Licence

MIT
