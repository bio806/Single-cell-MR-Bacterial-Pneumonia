library(dplyr)
library(tidyr)
library(TwoSampleMR)
library(coloc)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(ggpubr)

stopifnot(exists("exposure_dat"), exists("outcome_list"), is.list(outcome_list))

# 1. Exposure filtering and harmonization
exposure_dat <- exposure_dat %>%
  mutate(MAF.exposure = pmin(eaf.exposure, 1 - eaf.exposure)) %>%
  filter(is.na(MAF.exposure) | MAF.exposure > 0.01)

prepare_pair <- function(exposure_dat, outcome_dat) {
  h <- harmonise_data(exposure_dat, outcome_dat, action = 3) %>% filter(mr_keep)
  if ("gene_cell" %in% names(h)) {
    h$id.exposure <- gsub(" ", "_", h$gene_cell)
    h$exposure <- h$id.exposure
  }
  h$F_statistic <- (h$beta.exposure / h$se.exposure)^2
  h %>% filter(is.finite(F_statistic), F_statistic >= 10)
}

# 2. MR-PRESSO: only eligible exposure-outcome sets (>=4 instruments)
remove_presso_outliers <- function(dat, alpha = 0.05) {
  key <- interaction(dat$id.exposure, dat$id.outcome, drop = TRUE)
  parts <- split(dat, key)
  filtered <- lapply(parts, function(x) {
    if (nrow(x) < 4) return(x)
    fit <- tryCatch(TwoSampleMR::run_mr_presso(x, NbDistribution = max(1000, ceiling(nrow(x) / 1000) * 1000)), error = function(e) NULL)
    if (is.null(fit)) return(x)
    outliers <- fit[[1]][["MR-PRESSO results"]][["Outlier Test"]]
    if (is.null(outliers) || nrow(outliers) != nrow(x)) return(x)
    p <- suppressWarnings(as.numeric(gsub("<", "", as.character(outliers[, "Pvalue"]))))
    if (all(is.na(p))) return(x)
    x[is.na(p) | p >= alpha, , drop = FALSE]
  })
  bind_rows(filtered)
}

instrument_data <- lapply(outcome_list, function(outcome_dat) {
  remove_presso_outliers(prepare_pair(exposure_dat, outcome_dat))
})
mr_data <- bind_rows(instrument_data, .id = "dataset")

# 3. Instrument strength summary (unique exposure-SNP combinations)
f_summary <- mr_data %>% distinct(id.exposure, SNP, .keep_all = TRUE) %>%
  summarise(n_instruments = n(), median_F = median(F_statistic, na.rm = TRUE),
            Q1 = unname(quantile(F_statistic, 0.25, na.rm = TRUE)),
            Q3 = unname(quantile(F_statistic, 0.75, na.rm = TRUE)),
            min_F = min(F_statistic, na.rm = TRUE), max_F = max(F_statistic, na.rm = TRUE))
print(f_summary)

# 4. MR estimates and sensitivity analyses
mr_results <- mr(mr_data, method_list = c("mr_wald_ratio", "mr_ivw"))
mr_or <- generate_odds_ratios(mr_results) %>%
  filter(method %in% c("Wald ratio", "Inverse variance weighted")) %>%
  group_by(id.outcome) %>% mutate(FDR = p.adjust(pval, method = "BH")) %>% ungroup()
heterogeneity <- mr_heterogeneity(mr_data)
pleiotropy <- tryCatch(mr_pleiotropy_test(mr_data), error = function(e) NULL)
directionality <- tryCatch(directionality_test(mr_data), error = function(e) NULL)

# Significant same-direction associations in two outcomes
consistent_mr <- mr_or %>% filter(pval < 0.05) %>%
  group_by(id.exposure) %>%
  filter(n_distinct(id.outcome) == length(outcome_list), all(b > 0) | all(b < 0)) %>%
  ungroup()

# 5. Colocalization: regional summary data; requires full regional input, not clumped IVs
# One row per variant in each dataset. Use the same genome build and allele harmonization.
coloc_one <- function(eqtls, gwas, label, flank = 1e5) {
  x <- eqtls %>% filter(exposure == label) %>%
    filter(chrom == first(chr_gene), pos > first(start) - flank, pos < first(end) + flank) %>%
    distinct(SNP, .keep_all = TRUE)
  y <- gwas %>% distinct(SNP, .keep_all = TRUE)
  shared <- intersect(x$SNP, y$SNP)
  x <- x %>% filter(SNP %in% shared) %>% arrange(SNP)
  y <- y %>% filter(SNP %in% shared) %>% arrange(SNP)
  if (length(shared) < 2) return(NULL)
  fit <- coloc::coloc.abf(
    dataset1 = list(pvalues = y$P, snp = y$SNP, type = "cc", s = y$number_cases[1] / y$samplesize[1], N = y$samplesize[1], MAF = y$MAF),
    dataset2 = list(pvalues = x$P, snp = x$SNP, type = "quant", N = x$samplesize[1], MAF = x$MAF))
  s <- fit$summary