# replication analysis in all individuals when MGM smokes

library(meffil)
probedetails_epic <- meffil.featureset("epic")
rownames(probedetails_epic) <- probedetails_epic$name
probedetails_epic <- na.omit(probedetails_epic)

probedetails_450 <- meffil.featureset("450k")
rownames(probedetails_450) <- probedetails_450$name
probedetails_450 <- na.omit(probedetails_450)

# replication (EPIC adjusted for plate, 450k for slide)
# load EPIC and 450k EWAS results
load("path/to/mgm_smoke_ewas.Robj")
epic_all <- ewas.mgm
load("path/to/mgm_smoke_ewas_450k.Robj")
four50_all <- ewas.mgm

# extract p values
epic_pval <- epic_all$p.value
colnames(epic_pval) <- c("pval_none", "pval_all", "pval_sva")
four50_pval <- four50_all$p.value
colnames(four50_pval) <- c("pval_none", "pval_all", "pval_sva")

# extract coefficients
epic_effect <- epic_all$coefficient
colnames(epic_effect) <- c("effect_none", "effect_all", "effect_sva")
four50_effect <- four50_all$coefficient
colnames(four50_effect) <- c("effect_none", "effect_all", "effect_sva")

# get sites common to EPIC and 450k
epic_sites <- as.character(rownames(epic_pval))
four50_sites <- as.character(rownames(four50_pval))
common_sites <- intersect(epic_sites, four50_sites)
length(common_sites)

# merge p value and effct sizes for epic
epic_full <- merge(x=epic_pval, y=epic_effect, by.x="row.names", by.y="row.names")
head(epic_full)
epic_full <- epic_full[order(epic_full$pval_all),]
head(epic_full)
# write out top 25 sites and all sites for meta-analysis (as 
# meta analysis needs a txt file)
epic_full_25 <- epic_full[1:25,]
epic_full_25$chr <- probedetails_epic$chromosome[match(epic_full_25$Row.names, probedetails_epic$name)]
epic_full_25$position <- probedetails_epic$position[match(epic_full_25$Row.names, probedetails_epic$name)]
write.csv(epic_full_25, file="~/transgenerational_inheritance/redo/top_25/epic_all_mgm.csv")
write.csv(epic_full, file="~/transgenerational_inheritance/redo/meta_analysis/epic_all_mgm.csv", quote = F, row.names = F)
all_info <- epic_all$analyses$all$table
all_info <- all_info[common_sites,]
dim(all_info)
all_info$cpg <- as.character(rownames(all_info))
write.csv(all_info, file="~/transgenerational_inheritance/redo/meta_analysis/epic_all_mgm_allinfo.csv", quote = F, row.names = F)

# repeat for the 450k
four50_full <- merge(x=four50_pval, y=four50_effect, by.x="row.names", by.y="row.names")
head(four50_full)
four50_full <- four50_full[order(four50_full$pval_all),]
head(four50_full)
four50_full$Row.names <- as.character(four50_full$Row.names)
four50_full_25 <- four50_full[1:25,]
four50_full_25$chr <- probedetails_450$chromosome[match(four50_full_25$Row.names, probedetails_450$name)]
four50_full_25$position <- probedetails_450$position[match(four50_full_25$Row.names, probedetails_450$name)]
write.csv(four50_full_25, file="~/transgenerational_inheritance/redo/top_25/four50_all_mgm.csv")
write.csv(four50_full, file="~/transgenerational_inheritance/redo/meta_analysis/four50_all_mgm.csv", quote = F, row.names = F)
all_info <- four50_all$analyses$all$table
all_info <- all_info[common_sites,]
dim(all_info)
all_info$cpg <- as.character(rownames(all_info))
write.csv(all_info, file="~/transgenerational_inheritance/redo/meta_analysis/four50_all_mgm_allinfo.csv", quote = F, row.names = F)

epic_top25 <- as.character(epic_full$Row.names[1:25])
length(epic_top25)
epic_top25
four50_in_EPIC <- four50_full[four50_full$Row.names %in% epic_top25,]
head(four50_in_EPIC)
dim(four50_in_EPIC)
# merge 450k and epic data
all_results <- merge(epic_full,four50_full, by.x="Row.names", by.y="Row.names")
dim(all_results)
all_results <- all_results[order(all_results$pval_all.x),]
head(all_results)
# get 450k results for the top 25 epic sites
all_results25 <- all_results[1:25,]
all_results25[all_results25$pval_all.y < 0.05/nrow(all_results25),]
all_results25

# correlate effect sizes between epic and 450k for top 10,25,50,100,200 sites
top_range <- c(10,25,50,100,200)
cors_top_range <- function(x){
  all_results_i <- all_results[1:x,]
  cor_xy <- cor.test(all_results_i$effect_all.x, all_results_i$effect_all.y, alternative = c("two.sided"))
  return(cor_xy)
}

top_cors <- lapply(top_range, FUN = cors_top_range)

top_cors

# run binomial test for direction of effect
cors_binom_test <- function(y){
  all_results_i <- all_results[1:y,]
  all_results_i <- all_results_i[all_results_i$pval_all.y < 0.05,]
  all_results_i <- all_results_i[sign(all_results_i$effect_all.x)==sign(all_results_i$effect_all.y),]
  binomial.results <- binom.test(x = nrow(all_results_i), n = y, p = 0.05/2, alternative = c("greater"), conf.level = 0.95)
  return(binomial.results)
}
binomial.results.list <- lapply(top_range, FUN = cors_binom_test)
binomial.results.list
