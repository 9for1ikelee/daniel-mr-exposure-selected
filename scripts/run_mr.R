# ============================================================================
# 两样本孟德尔随机化：GALNT2 → NSCLC
# 暴露：GALNT2 pQTL (PROT-a-1171, GRCh37)
# 结局：NSCLC (FinnGen R12 C3_LUNG_NONSMALL_EXALLC, GRCh38)
# ============================================================================

rm(list=ls())
options(ieugwasr_jwt = Sys.getenv("IEUGWASR_JWT"))
if(!require("pacman")) install.packages("pacman", update=F, ask=F)
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
library("pacman")
p_load(data.table, tidyr, dplyr, purrr, readr, vroom,
       MendelianRandomization, MRPRESSO)
p_load_gh("mrcieu/gwasglue", "MRCIEU/TwoSampleMR")
library(TwoSampleMR)
library(MendelianRandomization)

setwd("/Users/muu/Desktop/ai-test")

cat("\n========================================\n")
cat("  GALNT2 -> NSCLC MR Analysis\n")
cat("========================================\n\n")

# ====== Step 1: 读取用户筛选的暴露IV ======
cat("[1/5] Reading user cleaned exposure IVs...\n")
expo_user <- fread("clean_confounder07.csv")
cat("  User IVs:", nrow(expo_user), "\n")

# ====== Step 2: 读取 FinnGen 结局并匹配 ======
cat("[2/5] Matching with FinnGen outcome...\n")
finngen <- vroom::vroom("finngen_R12_C3_LUNG_NONSMALL_EXALLC.gz",
                        show_col_types=FALSE)
outcome_raw <- finngen %>%
  rename(SNP="rsids", CHR="#chrom", BP="pos",
         effect_allele="alt", other_allele="ref",
         P="pval", EAF="af_alt", BETA="beta", SE="sebeta") %>%
  select(SNP, CHR, BP, effect_allele, other_allele,
         P, EAF, BETA, SE) %>%
  mutate(P=as.numeric(P))
cat("  FinnGen variants:", nrow(outcome_raw), "\n")

# 匹配用户IVs (处理逗号分隔rsID)
iv_snps <- expo_user$SNP
fg_rsids <- strsplit(as.character(outcome_raw$SNP), ",", fixed=TRUE)
match_idx <- which(sapply(fg_rsids, function(x) any(x %in% iv_snps)))
outcome_iv <- outcome_raw[match_idx, ]
cat("  Matched outcome SNPs:", nrow(outcome_iv), "\n")
write.csv(outcome_iv, "results/outcome_iv.csv", row.names=FALSE)

# ====== Step 3: 格式化为 MR 输入 ======
cat("[3/5] Formatting for MR...\n")
expo_MR <- read_exposure_data("clean_confounder07.csv",
  sep=",", snp_col="SNP", beta_col="beta.exposure",
  se_col="se.exposure", effect_allele_col="effect_allele.exposure",
  other_allele_col="other_allele.exposure",
  eaf_col="eaf.exposure", pval_col="pval.exposure",
  samplesize_col="samplesize.exposure", clump=F)
expo_MR$exposure <- "GALNT2"

outco_MR <- read_outcome_data(
  snps=expo_MR$SNP,
  filename="results/outcome_iv.csv",
  sep=",", snp_col="SNP", beta_col="BETA", se_col="SE",
  effect_allele_col="effect_allele", other_allele_col="other_allele",
  eaf_col="EAF", pval_col="P", chr_col="CHR", pos_col="BP")
outco_MR$outcome <- "NSCLC"

# 基因组版本对齐 (GRCh37 vs GRCh38)
pos_map <- expo_user[, .(SNP, chr.exposure, pos.exposure)]
outco_MR <- merge(outco_MR, pos_map, by="SNP", all.x=TRUE)
outco_MR$chr.outcome  <- outco_MR$chr.exposure
outco_MR$pos.outcome  <- outco_MR$pos.exposure
extra_cols <- grep("exposure$", colnames(outco_MR), value=TRUE)
for(col in extra_cols) outco_MR[[col]] <- NULL

cat("  Exposure:", nrow(expo_MR), " | Outcome:", nrow(outco_MR), "\n")

# ====== Step 4: Harmonise + MR ======
cat("[4/5] Harmonising & MR analysis...\n")
dat <- harmonise_data(exposure_dat=expo_MR, outcome_dat=outco_MR)
dat <- subset(dat, mr_keep==TRUE)
cat("  Harmonised:", nrow(dat), "\n")

# MR 主分析
mrResult <- mr(dat, method_list=c("mr_ivw", "mr_egger_regression",
                                    "mr_weighted_median"))
mrTab <- generate_odds_ratios(mrResult)
cat("\n--- MR Results ---\n")
print(mrTab, row.names=FALSE)
write.csv(mrTab, "results/MR_results.csv", row.names=FALSE)

# 异质性检验
heterTab <- mr_heterogeneity(dat)
cat("\n--- Heterogeneity ---\n")
print(heterTab, row.names=FALSE)
write.csv(heterTab, "results/heterogeneity.csv", row.names=FALSE)

# 多效性检验
pleioTab <- mr_pleiotropy_test(dat)
cat("\n--- Pleiotropy ---\n")
print(pleioTab, row.names=FALSE)
write.csv(pleioTab, "results/pleiotropy.csv", row.names=FALSE)

# MR-PRESSO 离群值检测
cat("\n--- MR-PRESSO ---\n")
tryCatch({
  presso <- run_mr_presso(dat, NbDistribution=1000, SignifThreshold=0.05)
  cat("Global Test P:", presso[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue, "\n")
  write.csv(presso[[1]]$`MR-PRESSO results`$`Outlier Test`,
            "results/MR_PRESSO_outlier.csv", row.names=FALSE)
}, error=function(e) cat("  Error:", e$message, "\n"))

# 留一法
loo_res <- mr_leaveoneout(dat)
write.csv(loo_res, "results/loo_results.csv", row.names=FALSE)

# 单SNP分析
res_single <- mr_singlesnp(dat)
write.csv(res_single, "results/singlesnp_results.csv", row.names=FALSE)

# ====== Step 5: 可视化 ======
cat("[5/5] Figures...\n")

pdf(file="figures/scatter_plot.pdf", width=9, height=8)
mr_scatter_plot(mrResult, dat)
dev.off()

pdf(file="figures/forest_plot.pdf", width=8, height=7)
mr_forest_plot(res_single)
dev.off()

pdf(file="figures/funnel_plot.pdf", width=8, height=7)
mr_funnel_plot(singlesnp_results=res_single)
dev.off()

pdf(file="figures/leaveoneout_plot.pdf", width=8, height=7)
mr_leaveoneout_plot(leaveoneout_results=loo_res)
dev.off()

# 脚本9：bioForest 自定义森林图
bioForest <- function(inputFile=NULL, forestFile=NULL, forestCol=NULL) {
  biofsci <- read.csv(inputFile, header=T, sep=",", check.names=F)
  row.names(biofsci) <- biofsci$method
  biofsci <- biofsci[biofsci$pval < 1, ]
  method <- rownames(biofsci)
  or <- sprintf("%.3f", biofsci$"or")
  orLow  <- sprintf("%.3f", biofsci$"or_lci95")
  orHigh <- sprintf("%.3f", biofsci$"or_uci95")
  OR <- paste0(or, "(", orLow, "-", orHigh, ")")
  pVal <- ifelse(biofsci$pval < 0.001, "<0.001", sprintf("%.3f", biofsci$pval))
  pdf(file=forestFile, width=7, height=4.6)
  n <- nrow(biofsci); nRow <- n + 1; ylim <- c(1, nRow)
  layout(matrix(c(1,2), nc=2), width=c(3.5, 2))
  par(mar=c(4, 2.5, 2, 1))
  plot(1, xlim=c(0, 3), ylim=ylim, type="n", axes=F, xlab="", ylab="")
  text(0, n:1, method, adj=0, cex=0.8)
  text(1.9, n:1, pVal, adj=1, cex=0.8)
  text(1.9, n+1, "pvalue", cex=1, font=2, adj=1)
  text(3.1, n:1, OR, adj=1, cex=0.8)
  text(2.7, n+1, "OR", cex=1, font=2, adj=1)
  par(mar=c(4, 1, 2, 1), mgp=c(2, 0.5, 0))
  xlim <- c(min(as.numeric(orLow)*0.975, as.numeric(orHigh)*0.975, 0.9),
            max(as.numeric(orLow), as.numeric(orHigh))*1.025)
  plot(1, xlim=xlim, ylim=ylim, type="n", axes=F, ylab="", xaxs="i", xlab="OR")
  arrows(as.numeric(orLow), n:1, as.numeric(orHigh), n:1,
         angle=90, code=3, length=0.05, col="darkgreen", lwd=3)
  abline(v=1, col="grey34", lty=2, lwd=2)
  boxcolor <- ifelse(as.numeric(or) > 1, forestCol, forestCol)
  points(as.numeric(or), n:1, pch=16, col=boxcolor, cex=2)
  axis(1)
  dev.off()
}
bioForest(inputFile="results/MR_results.csv",
          forestFile="figures/forest_plot_custom.pdf",
          forestCol="darkred")

cat("\n========================================\n")
cat("  Analysis Complete!\n")
cat("========================================\n")
