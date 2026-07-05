---
name: daniel-mr-exposure-selected
description: Run two-sample Mendelian Randomization (GALNT2 pQTL as exposure, FinnGen NSCLC as outcome) using user-selected exposure IVs from a pre-cleaned confounder file. Use when the user has already curated their own clean_confounder07.csv SNP list and needs to run the full MR pipeline matching against FinnGen GWAS data, including harmonisation, MR analysis (IVW, MR-Egger, weighted median), sensitivity tests (heterogeneity, pleiotropy, MR-PRESSO, leave-one-out), and visualization (scatter, forest, funnel, leave-one-out, custom bioForest plot).
---

# Daniel MR - Exposure Selected

Two-sample Mendelian Randomization using a user-provided, pre-cleaned confounder SNP list as exposure instruments, matched against FinnGen R12 GWAS summary statistics as outcome.

## Prerequisites

- R packages: TwoSampleMR, MendelianRandomization, MRPRESSO, data.table, dplyr, tidyr, vroom, pacman
- Exposure file: clean_confounder07.csv (user-curated IVs after confounder screening)
- Outcome file: finngen_R12_C3_LUNG_NONSMALL_EXALLC.gz (FinnGen R12 NSCLC GWAS)

## Workflow

Run the bundled script scripts/run_mr.R (same as mr_analysis.R in the working directory).
Steps: 1) Read exposure IVs, 2) Match FinnGen by rsID, 3) Format via TwoSampleMR, 4) Align genome builds (GRCh37 vs GRCh38), 5) Harmonise, 6) MR: IVW, MR-Egger, weighted median, 7) Sensitivity: heterogeneity, pleiotropy, MR-PRESSO, LOO, 8) 5 figures.

## Output Files

All output goes to results/ and figures/ subdirectories: MR_results.csv, heterogeneity.csv, pleiotropy.csv, MR_PRESSO_outlier.csv, loo_results.csv, singlesnp_results.csv, outcome_iv.csv, scatter_plot.pdf, forest_plot.pdf, funnel_plot.pdf, leaveoneout_plot.pdf, forest_plot_custom.pdf.

## Notes

- Genome build mismatch handled via rsID-based position alignment
- FinnGen rsids may contain comma-separated multiple rsIDs
- MR-PRESSO: 1000 distributions, significance threshold 0.05
- No API dependencies needed — all processing is local
