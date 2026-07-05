**daniel-mr-exposure-selected** — a Codex skill for two-sample Mendelian Randomization.

**What it does:** Runs a full MR pipeline using your pre-cleaned exposure SNP list (`clean_confounder07.csv`) as exposure instruments, matched against FinnGen outcomt GWAS data as outcome.

**Pipeline:**
1. Reads your curated IVs
2. Matches them against FinnGen by rsID (handles GRCh37/GRCh38 build mismatch)
3. Harmonises exposure–outcome alleles
4. Runs MR: IVW, MR-Egger, weighted median
5. Sensitivity: heterogeneity, pleiotropy, MR-PRESSO, leave-one-out
6. Generates 5 figures: scatter, forest, funnel, LOO, custom bioForest

**No API dependencies** — entirely local. Just needs your CSV files and R packages.

**Install:** `$skill-installer https://github.com/9for1ikelee/daniel-mr-exposure-selected`
