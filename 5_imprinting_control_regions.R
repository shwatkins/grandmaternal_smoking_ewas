# test whether hits are located in imprinting control regions

library(readxl)
# https://www.nature.com/articles/s41467-019-09301-y#Sec24
# supplementary data 1
imprinting_cpgs <- read_xlsx("/path/to/41467_2019_9301_MOESM4_ESM.xlsx")
ICR_cpgs <- as.character(imprinting_cpgs$CpG)

library(meffil)
probeDetails450k <- meffil.featureset("450k")
probeDetails450k <- probeDetails450k[!is.na(probeDetails450k$position),]
probeDetailsepic <- meffil.featureset("epic")
probeDetailsepic <- probeDetailsepic[!is.na(probeDetailsepic$position),]

ICRs_epic <- probeDetailsepic[probeDetailsepic$name %in% ICR_cpgs,]
dim(ICRs_epic)
ICRs_450k <- probeDetails450k[probeDetails450k$name %in% ICR_cpgs,]
dim(ICRs_450k)

# EWAS setup

ewas_ICR_results_epic <- function(x){
  ICR_results_p <- x$p.value
  colnames(ICR_results_p) <- c("pval_none", "pval_all", "pval_sva")
  head(ICR_results_p)
  ICR_results_e <- x$coefficient
  colnames(ICR_results_e) <- c("effect_none", "effect_all", "effect_sva")
  head(ICR_results_e)
  ICR_results <- merge(x=ICR_results_p, y=ICR_results_e, by.x="row.names", by.y="row.names")
  head(ICR_results)
  message(print("epic probe details dim:"))
  message(print(dim(ICRs_epic)))
  message(print(head(ICRs_epic)))
  ICR <- as.character(ICRs_epic$name)
  ICR <- ICR[ICR %in% ICR_results$Row.names]
  length(ICR)
  ICR_results <- ICR_results[ICR_results$Row.names %in% ICR,]
  dim(ICR_results)
  ICR_results <- ICR_results[order(ICR_results$pval_all),]
  ICR_cutoff <- 0.05/nrow(ICR_results)
  message(print("ICR cutoff epic:"))
  message(print(ICR_cutoff))
  ICR_results_corrected <- ICR_results[ICR_results$pval_all < ICR_cutoff,]
  dim(ICR_results_corrected)
  head(ICR_results_corrected)
  return(ICR_results_corrected)
}


ewas_ICR_results_450k <- function(x){
  ICR_results_p <- x$p.value
  colnames(ICR_results_p) <- c("pval_none", "pval_all", "pval_sva")
  head(ICR_results_p)
  ICR_results_e <- x$coefficient
  colnames(ICR_results_e) <- c("effect_none", "effect_all", "effect_sva")
  head(ICR_results_e)
  ICR_results <- merge(x=ICR_results_p, y=ICR_results_e, by.x="row.names", by.y="row.names")
  head(ICR_results)
  message(print("450k probe details dim:"))
  message(print(dim(ICRs_450k)))
  message(print(head(ICRs_450k)))
  ICR <- as.character(ICRs_450k$name)
  ICR <- ICR[ICR %in% ICR_results$Row.names]
  length(ICR)
  ICR_results <- ICR_results[ICR_results$Row.names %in% ICR,]
  dim(ICR_results)
  ICR_results <- ICR_results[order(ICR_results$pval_all),]
  ICR_cutoff <- 0.05/nrow(ICR_results)
  ICR_cutoff
  message(print("ICR cutoff 450k:"))
  message(print(ICR_cutoff))
  ICR_results_corrected <- ICR_results[ICR_results$pval_all < ICR_cutoff,]
  dim(ICR_results_corrected)
  head(ICR_results_corrected)
  return(ICR_results_corrected)
}

# EPIC MGM
load("/path/to/mgm_smoke_ewas.Robj")
ICR_epic_MGM <- ewas_ICR_results_epic(ewas.mgm)
dim(ICR_epic_MGM)
head(ICR_epic_MGM)

