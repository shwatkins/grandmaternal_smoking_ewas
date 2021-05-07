# meta-analysis post processing 
library(data.table)

# MGM all
# load output from METAL
mgm_all_ma <- fread("~/transgenerational_inheritance/redo/meta_analysis/meta_analysis_metal_mgm_all1.txt")
head(mgm_all_ma)
mgm_all_ma <- as.data.frame(mgm_all_ma)
dim(mgm_all_ma)
# rename columns
test <- colnames(mgm_all_ma)
test[6] <- "p_value"
colnames(mgm_all_ma) <- test
mgm_all_ma <- mgm_all_ma[order(mgm_all_ma$p_value),]
mgm_all_ma$dataset <- "mgm_all_ma"
head(mgm_all_ma)
# set p-value threshold (for 450k as this is an intersection between EPIC and 450k)
threshold <- 2.4e-7
# subset output to genome-wide significant sites
mgm_all_ma_cutoff <- mgm_all_ma[mgm_all_ma$p_value < threshold,]
mgm_all_ma_cutoff