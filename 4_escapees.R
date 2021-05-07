# test escapee sites

# load escapees
library(data.table)
# load escapee sites from Dietmann et al 2020 biorxiv paper
# https://www.biorxiv.org/content/10.1101/2020.03.19.998930v1.supplementary-material
# media-3.xls is Supplementary table 1
escapees <- fread("/path/to/media-3.xls")
escapees <- as.data.frame(escapees)
# take cols 1:3 to create genomic range for each escapee region
escapees_ranges <- escapees[,c(1,3)]
escapees_ranges$chr <- sapply(strsplit(escapees_ranges$HyperMR.hg19, ":"),'[', 1)
escapees_ranges$start <- sapply(strsplit(escapees_ranges$HyperMR.hg19, ":"),'[', 2)
escapees_ranges$stop <- sapply(strsplit(escapees_ranges$HyperMR.hg19, ":"),'[', 3)

# load details of sites on 450k and epic
library(meffil)
probeDetails450k <- meffil.featureset("450k")
probeDetails450k <- probeDetails450k[!is.na(probeDetails450k$position),]
probeDetailsepic <- meffil.featureset("epic")
probeDetailsepic <- probeDetailsepic[!is.na(probeDetailsepic$position),]

# create granges object of escapee regions, sites on 450k and sites on EPIC
library(GenomicRanges)
escapees_granges <- makeGRangesFromDataFrame(escapees_ranges, ignore.strand = T,
                                             seqnames.field = "chr", start.field = "start",
                                             end.field = "stop")
# 450k
probeDetails450k$stop <- probeDetails450k$position
granges_450k <- makeGRangesFromDataFrame(probeDetails450k, ignore.strand = T,
                                         seqnames.field = "chromosome", start.field = "position",
                                         end.field = "stop", keep.extra.columns = T)
# find the overlaps of escapee regions and sites on the 450k
escapees_450k <- findOverlaps(escapees_granges, granges_450k)
escapees_450k <- intersect(escapees_granges, granges_450k)
escapees_450k_df <- as.data.frame(escapees_450k)
escapees_450k_deets <- merge(probeDetails450k, escapees_450k_df, 
                             by.x = c("chromosome", "position"),
                             by.y = c("seqnames", "start"))
head(escapees_450k_deets)
dim(escapees_450k_deets)

# epic
probeDetailsepic$stop <- probeDetailsepic$position
granges_epic <- makeGRangesFromDataFrame(probeDetailsepic, ignore.strand = T,
                                         seqnames.field = "chromosome", start.field = "position",
                                         end.field = "stop", keep.extra.columns = T)
# find the overlaps of escapee regions and sites on epic
escapees_epic <- findOverlaps(escapees_granges, granges_epic)
escapees_epic <- intersect(escapees_granges, granges_epic)
escapees_epic_df <- as.data.frame(escapees_epic)
escapees_epic_deets <- merge(probeDetailsepic, escapees_epic_df, 
                             by.x = c("chromosome", "position"),
                             by.y = c("seqnames", "start"))
head(escapees_epic_deets)
dim(escapees_epic_deets)

# EWAS setup
# epic
ewas_escapees_results_epic <- function(x){
  escapees_results_p <- x$p.value
  colnames(escapees_results_p) <- c("pval_none", "pval_all", "pval_sva")
  escapees_results_e <- x$coefficient
  colnames(escapees_results_e) <- c("effect_none", "effect_all", "effect_sva")
  escapees_results <- merge(x=escapees_results_p, y=escapees_results_e, by.x="row.names", by.y="row.names")
  escapees <- probeDetailsepic[probeDetailsepic$name %in% escapees_epic_deets$name,]
  message(print("epic probe details dim:"))
  message(print(dim(escapees)))
  message(print(head(escapees)))
  escapees <- as.character(escapees$name)
  escapees <- escapees[escapees %in% escapees_results$Row.names]
  escapees_results <- escapees_results[escapees_results$Row.names %in% escapees,]
  dim(escapees_results)
  escapees_results <- escapees_results[order(escapees_results$pval_all),]
  escapees_cutoff <- 0.05/nrow(escapees_results)
  message(print("escapees cutoff epic:"))
  message(print(escapees_cutoff))
  escapees_results_corrected <- escapees_results[escapees_results$pval_all < escapees_cutoff,]
  return(escapees_results_corrected)
}

# 450k
ewas_escapees_results_450k <- function(x){
  escapees_results_p <- x$p.value
  colnames(escapees_results_p) <- c("pval_none", "pval_all", "pval_sva")
  escapees_results_e <- x$coefficient
  colnames(escapees_results_e) <- c("effect_none", "effect_all", "effect_sva")
  escapees_results <- merge(x=escapees_results_p, y=escapees_results_e, by.x="row.names", by.y="row.names")
  escapees <- probeDetails450k[probeDetails450k$name %in% escapees_450k_deets$name,]
  message(print(head(escapees)))
  escapees <- as.character(escapees$name)
  escapees <- escapees[escapees %in% escapees_results$Row.names]
  message(length(escapees))
  escapees_results <- escapees_results[escapees_results$Row.names %in% escapees,]
  escapees_results <- escapees_results[order(escapees_results$pval_all),]
  escapees_cutoff <- 0.05/nrow(escapees_results)
  message(print("escapees cutoff 450k:"))
  message(print(escapees_cutoff))
  escapees_results_corrected <- escapees_results[escapees_results$pval_all < escapees_cutoff,]
  dim(escapees_results_corrected)
  head(escapees_results_corrected)
  return(escapees_results_corrected)
}

# EPIC all MGM
load("/path/to/mgm_smoke_ewas.Robj")
escapees_epic_MGM <- ewas_escapees_results_epic(ewas.mgm)
dim(escapees_epic_MGM)
head(escapees_epic_MGM)

## 450k

# 450k females MGM
load("/path/to/mgm_smoke_ewas_450k.Robj")
escapees_450_females_MGM <- ewas_escapees_results_450k(ewas.mgm)
dim(escapees_450_females_MGM)
head(escapees_450_females_MGM)



