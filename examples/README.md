# Synthetic example data

Every file here is made up for the quickstart. Sample IDs (SYN01-SYN12), dosages and r2 values are invented and do not come from any real cohort.

- `synthetic.fam`: 12 samples, used by `plan`.
- `synthetic_hla.raw`: PLINK `--recode A` layout. `HLA_DRB1_0401_P` counts the present allele; `HLA_DRB1_1501_A` counts the absent allele, so the reader flips it.
- `batch/`: two SNP2HLA-style sub-batch outputs. `sub_batch_002_imputed` has no Beagle log on purpose, so `verify` reports it as incomplete.
