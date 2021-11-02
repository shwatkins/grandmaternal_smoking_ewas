## EWAS for all individuals when MGM smokes
# meth is methylation matrix
# ewas_covars is a dataframe containing the variable of interest and all covariates

# remove X and Y chromosomes
# create list of 450k probes inc gene and position
probeDetails <- meffil.featureset("epic")
rownames(probeDetails) <- probeDetails$name
probeDetails <- na.omit(probeDetails)
rownames(probeDetails) <- probeDetails$name

xchr <- probeDetails[probeDetails$chromosome=="chrX",]
ychr <- probeDetails[probeDetails$chromosome=="chrY",]
xchr <- as.character(xchr$name)
ychr <- as.character(ychr$name)
XY <- c(xchr, ychr)
length(XY)
dat_cpgs <- as.character(rownames(meth))
XY <- XY[XY %in% dat_cpgs]
XY <- na.omit(XY)
meth <- as.data.frame(meth)[! rownames(meth) %in% XY,]
meth <- as.matrix(meth)

# set covariates for EWAS
covars_for_ewas <- c("sex", "age", "plate", "Bcell", "CD4T", "CD8T", "Eos", "Mono", "Neu", "NK")

# run and save EWAS
ewas.mgm <- meffil.ewas(meth, variable=ewas_covars$mgm.smoke, covariates=ewas_covars[,covars_for_ewas], isva=F, random.seed=23) 
save(ewas.mgm, file="path/to/mgm_smoke_ewas.Robj")

# This function sets your default model and your significance threshold. It will show statistics for all four models for all CpGs that pass significance threshold with the default model.
ewas.parameters <- meffil.ewas.parameters(sig.threshold=9e-8,  ## EWAS p-value threshold
                                          max.plots=10, ## plot at most 100 CpG sites
                                          qq.inflation.method="median",  ## measure inflation using median
                                          model="all") ## select default EWAS model; 


ewas.summary<-meffil.ewas.summary(ewas.mgm,meth,parameters=ewas.parameters)                              

meffil.ewas.report(ewas.summary, output.file=paste("path/to/mgm_smoke.ewas.report.html",sep=""))
