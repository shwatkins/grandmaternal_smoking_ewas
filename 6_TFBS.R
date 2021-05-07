library(data.table)
library(LOLA)
library(simpleCache)
library(reshape2)
library(GenomicAlignments)
library(Rsamtools)
library(ggplot2)
library(biovizBase)
library(meffil)
library(dplyr)
library("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
library("BSgenome.Hsapiens.UCSC.hg19")
theme_set(theme_bw())

# create list of all probes in the analysis:
load("/path/to/mgm_smoke_ewas.Robj")
epic_pval <- as.data.frame(ewas.mgm$p.value)
colnames(epic_pval) <- c("pval_none", "pval_all", "pval_sva")
epic_pval <- epic_pval[order(epic_pval$pval_all),]

allcgs <- as.character(rownames(epic_pval))

intergen_top25 <- epic_pval[1:25,]
head(intergen_top25)

##load cell type conversion and colors
cellType_conversions=fread("~/CellTypes.tsv",drop="collection")
colors=fread("~/color.tsv")

setwd("path")

##load regiondb
# regiondb can be downloaded from lola
regionDB <- loadRegionDB("/resources/regions/LOLACore/hg19", collection="encode_tfbs")
head(regionDB)

retaincpg<-allcgs

#get 450k locations into a data table
Illumina450 <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Locations
Illumina450_dt=as.data.table(Illumina450)
Illumina450_dt[,cpgID:=row.names(Illumina450),]
Illumina450_dt <- Illumina450_dt[Illumina450_dt$cpgID%in%retaincpg,]
Illumina450_dt[,cpgstart_pre:=ifelse(strand=="-",pos-100,pos-99),]
Illumina450_dt[,cpgend_pre:=ifelse(strand=="-",pos+100,pos+101),]
data <- Illumina450_dt[Illumina450_dt$cpgID%in%rownames(intergen_top25),]
dim(data)

#collapse overlaps
gr_range = with(Illumina450_dt,GRanges(seqnames=chr,ranges=IRanges(cpgstart_pre,cpgend_pre)))
gr_cpg = with(Illumina450_dt,GRanges(seqnames=chr,ranges=IRanges(pos,pos)))

overlap=as.data.table(findOverlaps(gr_cpg, gr_range))
overlap_red=overlap[,list(subjectHit=min(subjectHits),NsubjectHits=.N),by=queryHits]

Illumina450_dt[,cpgstart:=start(gr_range[overlap_red$subjectHit])]
Illumina450_dt[,cpgend:=end(gr_range[overlap_red$subjectHit])]
Illumina450_dt[,NsubjectHits:=overlap_red$NsubjectHits]

head(data)
dim(data)
Illumina450_sub=Illumina450_dt[,c("cpgID","cpgstart","cpgend","pos"),with=FALSE]
setnames(Illumina450_sub,c("cpgID","pos"),c("cpg","ill_pos"))
setnames(data,c("cpgID"),c("cpg"))
head(Illumina450_sub)
head(data)
data=merge(data,Illumina450_sub,by="cpg",all.x=TRUE)
#check pos (should be all true)
table(data[,ill_pos==pos,])

####run with external background for CpGs

hg19_Illumina450_gr=with(Illumina450_dt, GRanges(seqnames = Rle(chr), IRanges(start=cpgstart, end=cpgend),strand=Rle(strand),ID=cpgID))
seq_Illumina450=getSeq(BSgenome.Hsapiens.UCSC.hg19,hg19_Illumina450_gr)

# add GC, CpG frequency of probes to the data table. Also add wheher or not the
# probe is in the list of top correlated probes.
Illumina450_dt[,GC_freq:=letterFrequency(seq_Illumina450, "CG", as.prob=T),]
Illumina450_dt[,CpG_freq:=dinucleotideFrequency(seq_Illumina450, step=2, as.prob=T)[,"CG"],]
Illumina450_dt[,isIntergenHit:=ifelse(cpgID%in%rownames(intergen_top25),TRUE,FALSE),]
Illumina450_dt<-Illumina450_dt[Illumina450_dt$cpgID%in%retaincpg,]

F7_cpg_gr=unique(with(Illumina450_dt,GRanges(seqnames = Rle(chr), IRanges(start=cpgstart_pre, end=cpgend_pre),strand=Rle("*"))))

#plot CG and CpG frequency for GoDMC cpgs and background 
jpeg("GCfreq_tfbs_enrichment_mgm_epic_all.jpg", width = 7, height = 5, units = 'in', res=300)
ggplot(Illumina450_dt,aes(x=GC_freq,col=isIntergenHit))+geom_density()
dev.off()
jpeg("CpGfreq_tfbs_enrichment_mgm_epic_all.jpg", width = 7, height = 5, units = 'in', res=300)
ggplot(Illumina450_dt,aes(x=CpG_freq,col=isIntergenHit))+geom_density()
dev.off()

background <- Illumina450_dt
background <- unique(with(Illumina450_dt,GRanges(seqnames = Rle(chr), IRanges(start=cpgstart, end=cpgend),strand=Rle("*"))))
intergen_hits <- unique(with(data,GRanges(seqnames = Rle(chr), IRanges(start=cpgstart, end=cpgend),strand=Rle("*"))))
intergen_hits_bed <- as.data.frame(intergen_hits)
head(intergen_hits_bed)
intergen_hits_bed <- intergen_hits_bed[,c(1:3)]
dim(intergen_hits_bed)
write.table(intergen_hits_bed, file="intergen_hits_bed_mgm_epic_all.txt", quote = F, row.names = F, col.names = F)

background_bed <- as.data.frame(background)
background_bed <- background_bed[,c(1:3)]
dim(background_bed)
write.table(background_bed, file="background_bed_mgm_epic_all.txt", quote = F, row.names = F, col.names = F)
random_bed <- background_bed[sample(nrow(background_bed), 250), ]
dim(random_bed)
write.table(random_bed, file="random_bed_mgm_epic_all.txt", quote = F, row.names = F, col.names = F)

random_25 <- as.data.frame(background)
random_25 <- random_25[sample(nrow(random_25), 25), ]
random_25 <- with(random_25,GRanges(seqnames = Rle(seqnames), IRanges(start=start, end=end),strand=Rle("*")))


dim(Illumina450_dt)
dim(background)
lola_res0_matched=runLOLA(intergen_hits, background, regionDB, cores=1)
lola_res0_matched$logOddsRatio<-log(lola_res0_matched$oddsRatio)
save(lola_res0_matched,file="lola_encodetfbs_intergen_mgm_epic_all.Rdata")
names(lola_res0_matched)
head(lola_res0_matched)
dat_sort <- lola_res0_matched[order(lola_res0_matched$oddsRatio, decreasing = T),]
dat_sort$pvalue <- 10^-dat_sort$pValueLog
head(dat_sort)
dat_sort <- dat_sort[order(dat_sort$pvalue),]
dat_sort_pval <- dat_sort[dat_sort$pvalue < 0.05,]
dat_sort_pval

lola_res0_matched$logOddsRatio
lola_res0_matched$tf <- sapply(strsplit(lola_res0_matched$antibody,"_"), `[`, 1)

pdf(paste0("intergen_encodetfbs_mgm_epic_all.pdf"),height=12,width=24)
pl3=ggplot(lola_res0_matched,aes(x=tf,y=logOddsRatio))+
  geom_hline(yintercept=0, linetype="dotted")+
  geom_jitter(width = 0.2, aes(colour=cellType,size=pValueLog))+
  facet_wrap(~userSet,scale="free_x",ncol=1)+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="bottom")+
  scale_size(range=c(1,8))+
  guides(fill = guide_legend(ncol=20))+
  theme(legend.text=element_text(size=12)) +
  labs(y="Odds ratio (log scale)",x="Encode segmentation") +
  scale_fill_brewer(type="qual")

print(pl3)
dev.off()

# change to cell types

cellType_conversions=fread("~/CellTypes.tsv",drop="collection")

process_LOLA = function (LOLA_res, collections=c("codex","encode_tfbs"),cellType_conversions){
  
  LOLA_res=LOLA_res[!is.na(userSet)]
  
  LOLA_res=LOLA_res[collection%in%collections]
  #changed form exp() to 10^ to accomodate change in LOLA
  
  LOLA_res[,p.adjust:=p.adjust(10^(-pValueLog),method="BY"),by=userSet]
  LOLA_res[,mlog10p.adjust:=-log10(p.adjust),]
  
  ##standardize cellTypes
  LOLA_res[cellType==""|is.na(cellType),cellType:="Not defined",]
  LOLA_res=merge(LOLA_res,cellType_conversions,by="cellType",all.x=TRUE)
  ##correct wrong annotation
  LOLA_res[description=="T-cell acute lymphoblastic leukaemia (T-ALL) cell line.",c("Lineage1","Lineage","cellType_corr"):=list(Lineage1="Lymphoid",Lineage="Lymphoid",cellType_corr="T lymphocyte"),]
  LOLA_res[,lineage_count_allstate:=length(unique(filename[!is.na(filename)])),by=c("Lineage","cellState")]
  LOLA_res[,lineage_count_all:=length(unique(filename[!is.na(filename)])),by=c("Lineage")]
  
  #standardize antibodies
  LOLA_res[,target:=toupper(sub("-","",unlist(lapply(antibody,function(x){spl=unlist(strsplit(x,"_|eGFP-"));spl[spl!=""][1]})))),]
  
  return(LOLA_res)  
  
}

locResults=process_LOLA(LOLA_res=lola_res0_matched,cellType_conversions=cellType_conversions)
theme_set(theme_bw())
pdf(paste0("intergen_encodetfbs_tissues_mgm_epic_all.pdf"),height=12,width=24)
pl3=ggplot(locResults,aes(x=target,y=logOddsRatio))+
  geom_hline(yintercept=0, linetype="dotted")+
  geom_jitter(width = 0.2, aes(colour=Tissue,size=pValueLog))+
  facet_wrap(~userSet,scale="free_x",ncol=1)+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="bottom")+
  scale_size(range=c(1,8))+
  guides(fill = guide_legend(ncol=20))+
  theme(legend.text=element_text(size=14), 
        legend.title=element_text(size=16, face = "bold"),
        axis.title=element_text(size=16)) +
  labs(y="Odds ratio (log scale)",x="Transcription factor") +
  scale_fill_brewer(type="qual")

print(pl3)
dev.off()

save(locResults, file = "mgm_epic_all_locresults.Rdata")

# now with TFs at p<0.05
dat_sort_pval$tf <- sapply(strsplit(dat_sort_pval$antibody,"_"), `[`, 1)
locResults2=process_LOLA(LOLA_res=dat_sort_pval,cellType_conversions=cellType_conversions)
save(locResults2, file="intergen_pval_table_mgm_epic_all.Rdata")
oneresult <- locResults2[,1]
theme_set(theme_bw())
pdf(paste0("intergen_pvalcutoff_encodetfbs_tissues_mgm_epic_all.pdf"),height=12,width=24)
pl3=ggplot(locResults2,aes(x=target,y=logOddsRatio))+
  geom_hline(yintercept=0, linetype="dotted")+
  geom_jitter(width = 0.2, aes(colour=Tissue,size=pValueLog))+
  facet_wrap(~userSet,scale="free_x",ncol=1)+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="bottom")+
  scale_size(range=c(1,8))+
  guides(fill = guide_legend(ncol=20))+
  theme(legend.text=element_text(size=14), 
        legend.title=element_text(size=16, face = "bold"),
        axis.title=element_text(size=16)) +
  labs(y="Odds ratio (log scale)",x="Transcription factor") +
  scale_fill_brewer(type="qual")

print(pl3)
dev.off()

