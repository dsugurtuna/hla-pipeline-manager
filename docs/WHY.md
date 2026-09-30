# Why it's built this way

## The problem

HLA imputation of a large cohort is thousands of small SNP2HLA runs, and any one of them can fail quietly. The expensive mistakes are releasing an incomplete batch and overwriting a good release without a way back.

## Design choices

**Why separate planning from running?** Because a plan can be read, reviewed and tested without PLINK, Java or a reference panel, and because "what would run" and "what did run" must never be confused. An earlier version returned a 100% success rate for work it had not done.

**Why an injectable runner instead of calling `subprocess` directly?** Because the interesting logic is failure handling: stop a sub-batch at its first failed command, stop everything if a shared step fails. A fake runner tests that in milliseconds.

**Why split sub-batches from the input `.fam`?** Because extracting a region with `--chr/--from-bp/--to-bp` removes variants, never samples. The split is known before any tool runs, so the whole plan can be printed up front.

**Why write `.keep` files rather than `.fam` files for the split?** Because PLINK writes `<out>.fam`. If the keep file had the same name, PLINK would overwrite its own input.

**Why pass `--update-name <file> 2 1` explicitly?** Because the default column order (new ID in column 2, old in column 1) is easy to get backwards, and writing it out makes the command self-documenting.

**Why orient every dosage to the "present" allele?** Because SNP2HLA codes HLA alleles as presence (`P`) or absence (`A`), and both PLINK `--recode A` and SNP2HLA's `.dosage` count whichever allele comes first. If that happens to be `A`, carriers look like non-carriers. The readers check the counted allele and flip `x -> 2 - x` for `HLA_` markers only, because for ordinary SNPs `A` means adenine.

**Why verify by SHA-256 rather than "the file exists"?** Because a truncated copy exists too. A checksum compares what arrived with what was sent.

**Why is deployment not a CLI command?** Because it is the only operation that overwrites release data. Planning, verifying and reporting read or write local outputs; deployment has to be called deliberately from code, backs up first and checks every copy.

**Why keep `legacy/` at all?** Because it shows where the package came from and what behaviour it had to preserve (the verifier's checks are a direct translation). It is excluded from linting and not maintained.

## Questions worth asking

**"How do you know an imputed allele is right, not just present?"**
This package does not. It checks that outputs are structurally complete and reports each marker's Beagle r2, which estimates imputation quality from the posterior probabilities. Accuracy needs samples with lab HLA typing to compare against; that is a separate validation step and should be done per reference panel and per ancestry group.

**"Why hard calls at 0.5 and 1.5? Doesn't that throw away uncertainty?"**
Yes. The thresholds turn an expected allele count into a carrier list, which is what a recall study or look-up needs. For association analysis you should use the dosages directly. The report keeps the dosage and the r2 next to each call so a reader can see borderline cases.

**"What happens with strand or build mismatches between the array and the reference panel?"**
SNP2HLA itself aligns alleles to the reference, flips what it can and removes ambiguous A/T and C/G SNPs with discordant frequencies. This package does not repeat that. It does make the region and the ID renaming explicit, and the positions assume the same genome build as the reference panel; a build mismatch would need a liftover before planning.

## What's next

- Write the plan out as a Slurm array job, one task per sub-batch.
- Record tool versions and input checksums in a manifest next to the outputs.
- Add CookHLA output support and compare the two on the same synthetic input.
