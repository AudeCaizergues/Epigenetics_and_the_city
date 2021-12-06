#!/usr/bin/Rscript

####################################################################################################
#######################   DMR BETWEEN HABITATS AND SEXES     #######################################
####################################################################################################

#### Libraries ######################################################################################
library(devtools)
library(methylKit)
library (graphics)
library (tools)
#####################################################################################################

#####################################################################################################
#######  BARCELONA ##################################################################################
#####################################################################################################

### BUILD OBJECT #####
setwd("/PATH/TO/06_results/methylkit/")
file.list.bar = list("bar10_S64_CpG.txt",
                  "bar11_S75_CpG.txt",
                  "bar12_S25_CpG.txt",
                  "bar13_S42_CpG.txt",
                  "bar14_S3_CpG.txt",
                  "bar15_S70_CpG.txt",
                  "bar16_S10_CpG.txt",
                  "bar17_S54_CpG.txt",
                  "bar18_S52_CpG.txt",
                  "bar19b_S12_CpG.txt",
                  "bar1_S46_CpG.txt",
                  "bar20_S60_CpG.txt",
                  "bar2_S23_CpG.txt",
                  "bar3_S29_CpG.txt",
                  "bar4_S63_CpG.txt",
                  "bar5_S45_CpG.txt",
                  "bar6_S7_CpG.txt",
                  "bar7_S16_CpG.txt",
                  "bar8_S59_CpG.txt",
                  "bar9_S35_CpG.txt")                  

sample.id.bar=list("bar10_S64","bar11_S75","bar12_S25","bar13_S42","bar14_S3","bar15_S70","bar16_S10","bar17_S54","bar18_S52","bar19b_S12",
               "bar1_S46","bar20_S60","bar2_S23","bar3_S29","bar4_S63","bar5_S45","bar6_S7","bar7_S16","bar8_S59","bar9_S35")


### URBAN EFFECT ####
#rur=1 urb=2
myobj_bar=methRead(file.list.bar,
                   sample.id.bar,
                   assembly="unknown",
                   treatment=c(2,2,2,2,2,2,2,1,2,1,1,2,1,1,1,1,1,1,1,2),
                   context="CpG",
                   mincov=3)

setwd("/PATH/TO/06_results/methylkit/methylkit_bar/")
#### Normalize
normalize.bar_allow_missing=normalizeCoverage(myobj_bar)
#test cov10
filtered.myobj.10x.bar_allow_missing=filterByCoverage(normalize.bar_allow_missing,lo.count=10,lo.perc=NULL,
                                   hi.count=NULL,hi.perc=99.9)
save(filtered.myobj.10x.bar_allow_missing, file = "filtered.myobj.10x.bar_allow_missing.rda")
#### Merging samples
meth_bar_allow_missing=unite(myobj_bar, min.per.group=9L ,destrand=TRUE)
save(meth_bar_allow_missing, file = "meth_bar_allow_missing.rda")
meth_filtered.myobj.10x.bar_allow_missing=unite(filtered.myobj.10x.bar_allow_missing, min.per.group=9L,destrand=TRUE)
save(meth_filtered.myobj.10x.bar_allow_missing, file = "meth_filtered.myobj.10x.bar_allow_missing.rda")

#### TILES 1000BP #################################################################################
### Cutting tiles
clusterSamples(meth_bar_allow_missing, dist="correlation", method="ward", plot=TRUE)
hist(OSBDiff$diff.meth)

hist(meth_filtered.myobj.10x.bar_allow_missing$diff.meth)
filtered.myobj_bar_allow_missing=filterByCoverage(myobj_bar,lo.count=10,lo.perc=NULL,hi.count=NULL,hi.perc=99.9)
normalized.filtered.myobj_bar_allow_missing=normalizeCoverage(filtered.myobj_bar_allow_missing)
tiles_bar_allow_missing=tileMethylCounts(normalized.filtered.myobj_bar_allow_missing,win.size=1000,step.size=1000,cov.bases=10)
meth_tiles_bar_allow_missing=unite(tiles_bar_allow_missing, min.per.group=9L, destrand=TRUE)
save(meth_tiles_bar_allow_missing, file = "meth_tiles_bar_allow_missing.rda")

###
#load("meth_tiles_bar_allow_missing.rda")
covariates=data.frame(sex=c("female","female","female","male","male","male","male","male","female","female","male","female","male","male","male","female","female","female","female","female"))
myDiff_bar_allow_missing=calculateDiffMeth(meth_tiles_bar_allow_missing, covariates=covariates)
write.table(myDiff_bar_allow_missing, file="dmr_tiles_bar_allow_missing_env.tsv", sep='\t')
my.diffMeth_10_bar_allow_missing = getMethylDiff(myDiff_bar_allow_missing,difference=10,qvalue=0.05)
write.table(my.diffMeth_10_bar_allow_missing, file="dmr_tiles_bar_allow_missing_10pourcent_env.tsv", sep='\t', quote=FALSE)

###########

### SEX EFFECT ####
#male=1 female=2
setwd("/PATH/TO/06_results/methylkit/")

myobj_bar_allow_missing_sex=methRead(file.list.bar,
                   sample.id.bar,
                   assembly="unknown",
                   treatment=c(2,2,2,1,1,1,1,1,2,2,1,2,1,1,1,2,2,2,2,2),
                   context="CpG",
                   mincov=3)

setwd("/PATH/TO/06_results/methylkit/methylkit_bar/")
#### Normalize
normalize.bar_allow_missing_sex=normalizeCoverage(myobj_bar_allow_missing_sex)
#test 10
filtered.myobj.10x.bar_allow_missing_sex=filterByCoverage(normalize.bar_allow_missing_sex,lo.count=10,lo.perc=NULL,
                                        hi.count=NULL,hi.perc=99.9)
save(filtered.myobj.10x.bar_allow_missing_sex, file = "filtered.myobj.10x.bar_allow_missing_sex.rda")
#### Merging samples
meth_bar_allow_missing_sex=unite(myobj_bar_allow_missing_sex, destrand=TRUE)
save(meth_bar_allow_missing_sex_sex, file = "meth_bar_allow_missing_sex.rda")
meth_filtered.myobj.10x.bar_allow_missing_sex=unite(filtered.myobj.10x.bar_allow_missing_sex, min.per.group=9L, destrand=TRUE)
save(meth_filtered.myobj.10x.bar_allow_missing_sex, file = "meth_filtered.myobj.10x.bar_allow_missing_sex.rda")

#### TILES 1000BP #################################################################################
### Cutting tiles
clusterSamples(meth_bar_allow_missing_sex, dist="correlation", method="ward", plot=TRUE)

hist(meth_filtered.myobj.10x.bar_allow_missing_sex$diff.meth)
filtered.myobj_bar_allow_missing_sex=filterByCoverage(myobj_bar_allow_missing_sex,lo.count=10,lo.perc=NULL,hi.count=NULL,hi.perc=99.9)
normalized.filtered.myobj_bar_allow_missing_sex=normalizeCoverage(filtered.myobj_bar_allow_missing_sex)
tiles_bar_allow_missing_sex=tileMethylCounts(normalized.filtered.myobj_bar_allow_missing_sex,win.size=1000,step.size=1000,cov.bases=10)
meth_tiles_bar_allow_missing_sex=unite(tiles_bar_allow_missing_sex, min.per.group=9L, destrand=TRUE)
save(meth_tiles_bar_allow_missing_sex, file = "meth_tiles_bar_allow_missing_sex.rda")

###
#load("meth_tiles_bar_allow_missing.rda")

covariates2=data.frame(habitat=c("urb","urb","urb","urb","urb","urb","urb","rur","urb","rur","rur","urb","rur","rur","rur","rur","rur","rur","rur","urb"))
myDiff_bar_allow_missing_sex=calculateDiffMeth(meth_tiles_bar_allow_missing_sex, covariates=covariates2)
write.table(myDiff_bar_allow_missing_sex, file="dmr_tiles_bar_allow_missing_sex.tsv", sep='\t')
my.diffMeth_10_bar_allow_missing_sex = getMethylDiff(myDiff_bar_allow_missing_sex,difference=10,qvalue=0.05)
write.table(my.diffMeth_10_bar_allow_missing_sex, file="dmr_tiles_bar_allow_missing_10pourcent_sex.tsv", sep='\t', quote=FALSE)

###########


#####################################################################################################
#######  MONTPELLIER ################################################################################
#####################################################################################################

### BUILD OBJECT #####
setwd("/PATH/TO/06_results/methylkit/")
file.list.mtp = list("mtp10_S19_CpG.txt",
                     "mtp16b_S33_CpG.txt",
                     "mtp17_S24_CpG.txt",
                     "mtp18_S66_CpG.txt",
                     "mtp19_S51_CpG.txt",
                     "mtp1_S62_CpG.txt",
                     "mtp20_S49_CpG.txt",
                     "mtp21_S6_CpG.txt",
                     "mtp22_S21_CpG.txt",
                     "mtp23_S37_CpG.txt",
                     "mtp24_S32_CpG.txt",
                     "mtp25_S38_CpG.txt",
                     "mtp2_S14_CpG.txt",
                     "mtp3_S34_CpG.txt",
                     "mtp4_S58_CpG.txt",
                     "mtp5_S18_CpG.txt",
                     "mtp6_S78_CpG.txt",
                     "mtp7_S44_CpG.txt",
                     "mtp8_S27_CpG.txt",
                     "mtp9_S43_CpG.txt")                  

sample.id.mtp=list("mtp10_S19","mtp16b_S33","mtp17_S24","mtp18_S66","mtp19_S51","mtp1_S62","mtp20_S49","mtp21_S6","mtp22_S21","mtp23_S37","mtp24_S32","mtp25_S38",
                   "mtp2_S14","mtp3_S34","mtp4_S58","mtp5_S18","mtp6_S78","mtp7_S44","mtp8_S27","mtp9_S43")

### URBAN EFFECT ####
#rur=1 urb=2
myobj_mtp=methRead(file.list.mtp,
                   sample.id.mtp,
                   assembly="unknown",
                   treatment=c(2,1,1,1,1,1,1,2,2,2,2,2,1,1,1,1,2,2,2,2),
                   context="CpG",
                   mincov=3)

setwd("PATH/TO/06_results/methylkit/methylkit_mtp/")

#### Normalize
normalize.mtp_allow_missing=normalizeCoverage(myobj_mtp)
#test 10
filtered.myobj.10x.mtp_allow_missing=filterByCoverage(normalize.mtp_allow_missing,lo.count=10,lo.perc=NULL,
                                        hi.count=NULL,hi.perc=99.9)
save(filtered.myobj.10x.mtp_allow_missing, file = "filtered.myobj.10x.mtp_allow_missing.rda")
#### Merging samples
meth_mtp_allow_missing=unite(myobj_mtp, min.per.group=9L, destrand=TRUE)
save(meth_mtp_allow_missing, file = "meth_mtp_allow_missing.rda")
meth_filtered.myobj.10x.mtp_allow_missing=unite(filtered.myobj.10x.mtp_allow_missing, min.per.group=9L, destrand=TRUE)
save(meth_filtered.myobj.10x.mtp_allow_missing, file = "meth_filtered.myobj.10x.mtp_allow_missing.rda")

#### TILES 1000BP #################################################################################
### Cutting tiles
clusterSamples(meth_mtp_allow_missing, dist="correlation", method="ward", plot=TRUE)
hist(OSBDiff$diff.meth)

hist(meth_filtered.myobj.10x.mtp_allow_missing$diff.meth)
filtered.myobj_mtp_allow_missing=filterByCoverage(myobj_mtp,lo.count=10,lo.perc=NULL,hi.count=NULL,hi.perc=99.9)
normalized.filtered.myobj_mtp_allow_missing=normalizeCoverage(filtered.myobj_mtp_allow_missing)
tiles_mtp_allow_missing=tileMethylCounts(normalized.filtered.myobj_mtp_allow_missing,win.size=1000,step.size=1000,cov.bases=10)
meth_tiles_mtp_allow_missing=unite(tiles_mtp_allow_missing, min.per.group=9L, destrand=TRUE)
save(meth_tiles_mtp_allow_missing, file = "meth_tiles_mtp_allow_missing.rda")

###
#load("meth_tiles_mtp_allow_missing.rda")
covariates=data.frame(sex=c("male","female","female","female","male","male","male","female","male","male","female","male","male","female","male","female","female","male","female","female"))
myDiff_mtp_allow_missing=calculateDiffMeth(meth_tiles_mtp_allow_missing, covariates=covariates)
write.table(myDiff_mtp_allow_missing, file="dmr_tiles_mtp_allow_missing_env.tsv", sep='\t')
my.diffMeth_10_mtp_allow_missing = getMethylDiff(myDiff_mtp_allow_missing,difference=10,qvalue=0.05)
write.table(my.diffMeth_10_mtp_allow_missing, file="dmr_tiles_mtp_allow_missing_10pourcent_env.tsv", sep='\t', quote=FALSE)

###########

### SEX EFFECT ####
#male=1 female=2
setwd("/PATH/TO/06_results/methylkit/")

myobj_mtp_allow_missing_sex=methRead(file.list.mtp,
                       sample.id.mtp,
                       assembly="unknown",
                       treatment=c(1,2,2,2,1,1,1,2,1,1,2,1,1,2,1,2,2,1,2,2),
                       context="CpG",
                       mincov=3)

setwd("/PATH/TO/06_results/methylkit/methylkit_mtp/")
#### Normalize
normalize.mtp_allow_missing_sex=normalizeCoverage(myobj_mtp_allow_missing_sex)
#test 10
filtered.myobj.10x.mtp_allow_missing_sex=filterByCoverage(normalize.mtp_allow_missing_sex,lo.count=10,lo.perc=NULL,
                                            hi.count=NULL,hi.perc=99.9)
save(filtered.myobj.10x.mtp_allow_missing_sex, file = "filtered.myobj.10x.mtp_allow_missing_sex.rda")

#### Merging samples
meth_mtp_allow_missing_sex=unite(myobj_mtp_allow_missing_sex, min.per.group=9L, destrand=TRUE)
save(meth_mtp_allow_missing_sex, file = "meth_mtp_allow_missing_sex.rda")
meth_filtered.myobj.10x.mtp_allow_missing_sex=unite(filtered.myobj.10x.mtp_allow_missing_sex, min.per.group=9L, destrand=TRUE)
save(meth_filtered.myobj.10x.mtp_allow_missing_sex, file = "meth_filtered.myobj.10x.mtp_allow_missing_sex.rda")

#### TILES 1000BP #################################################################################
### Cutting tiles
clusterSamples(meth_mtp_allow_missing_sex, dist="correlation", method="ward", plot=TRUE)
hist(meth_filtered.myobj.10x.mtp_allow_missing_sex$diff.meth)
filtered.myobj_mtp_allow_missing_sex=filterByCoverage(myobj_mtp_allow_missing_sex,lo.count=10,lo.perc=NULL,hi.count=NULL,hi.perc=99.9)
normalized.filtered.myobj_mtp_allow_missing_sex=normalizeCoverage(filtered.myobj_mtp_allow_missing_sex)
tiles_mtp_allow_missing_sex=tileMethylCounts(normalized.filtered.myobj_mtp_allow_missing_sex,win.size=1000,step.size=1000,cov.bases=10)
meth_tiles_mtp_allow_missing_sex=unite(tiles_mtp_allow_missing_sex, min.per.group=9L, destrand=TRUE)
save(meth_tiles_mtp_allow_missing_sex, file = "meth_tiles_mtp_allow_missing_sex.rda")

###
#load("meth_tiles_mtp_allow_missing.rda")
covariates2=data.frame(habitat=c("urb","rur","rur","rur","rur","rur","rur","urb","urb","urb","urb","urb","rur","rur","rur","rur","urb","urb","urb","urb"))
myDiff_mtp_allow_missing_sex=calculateDiffMeth(meth_tiles_mtp_allow_missing_sex, covariates=covariates2)
write.table(myDiff_mtp_allow_missing_sex, file="dmr_tiles_mtp_allow_missing_sex.tsv", sep='\t')
my.diffMeth_10_mtp_allow_missing_sex = getMethylDiff(myDiff_mtp_allow_missing_sex,difference=10,qvalue=0.05)
write.table(my.diffMeth_10_mtp_allow_missing_sex, file="dmr_tiles_mtp_allow_missing_10pourcent_sex.tsv", sep='\t', quote=FALSE)

###########


#####################################################################################################
#######  WARSAW    ################################################################################
#####################################################################################################

### BUILD OBJECT #####
setwd("/PATH/TO/06_results/methylkit/")
file.list.var = list("var10_S8_CpG.txt",
                     "var11_S80_CpG.txt",
                     "var12_S67_CpG.txt",
                     "var13_S74_CpG.txt",
                     "var14_S81_CpG.txt",
                     "var15_S72_CpG.txt",
                     "var16_S50_CpG.txt",
                     "var17_S53_CpG.txt",
                     "var18_S39_CpG.txt",
                     "var19_S76_CpG.txt",
                     "var1_S79_CpG.txt",
                     "var20_S61_CpG.txt",
                     "var2_S77_CpG.txt",
                     "var3_S71_CpG.txt",
                     "var4_S68_CpG.txt",
                     "var5_S69_CpG.txt",
                     "var6_S65_CpG.txt",
                     "var7b_S30_CpG.txt",
                     "var8_S31_CpG.txt",
                     "var9_S41_CpG.txt")                  

sample.id.var=list("var10_S8","var11_S80","var12_S67","var13_S74","var14_S81","var15_S72","var16_S50","var17_S53","var18_S39","var19_S76",
                   "var1_S79","var20_S61","var2_S77","var3_S71","var4_S68","var5_S69","var6_S65","var7b_S30","var8_S31","var9_S41")

### URBAN EFFECT ####
#rur=1 urb=2
myobj_var_allow_missing=methRead(file.list.var,
                   sample.id.var,
                   assembly="unknown",
                   treatment=c(1,2,2,2,2,2,2,2,2,2,1,2,1,1,1,1,1,1,1,1),
                   context="CpG",
                   mincov=3)

setwd("/PATH/TO/06_results/methylkit/methylkit_var/")
#### Normalize
normalize.var_allow_missing=normalizeCoverage(myobj_var_allow_missing)
#test 10
filtered.myobj.10x.var_allow_missing=filterByCoverage(normalize.var_allow_missing,lo.count=10,lo.perc=NULL,
                                        hi.count=NULL,hi.perc=99.9)
save(filtered.myobj.10x.var_allow_missing, file = "filtered.myobj.10x.var_allow_missing.rda")
#### Merging samples
meth_var_allow_missing=unite(myobj_var_allow_missing, min.per.group=9L, destrand=TRUE)
save(meth_var_allow_missing, file = "meth_var_allow_missing.rda")
meth_filtered.myobj.10x.var_allow_missing=unite(filtered.myobj.10x.var_allow_missing, min.per.group=9L, destrand=TRUE)
save(meth_filtered.myobj.10x.var_allow_missing, file = "meth_filtered.myobj.10x.var_allow_missing.rda")

#### TILES 1000BP #################################################################################
### Cutting tiles
clusterSamples(meth_var_allow_missing, dist="correlation", method="ward", plot=TRUE)
hist(OSBDiff$diff.meth)

hist(meth_filtered.myobj.10x.var_allow_missing$diff.meth)
filtered.myobj_var_allow_missing=filterByCoverage(myobj_var_allow_missing,lo.count=10,lo.perc=NULL,hi.count=NULL,hi.perc=99.9)
normalized.filtered.myobj_var_allow_missing=normalizeCoverage(filtered.myobj_var_allow_missing)
tiles_var_allow_missing=tileMethylCounts(normalized.filtered.myobj_var_allow_missing,win.size=1000,step.size=1000,cov.bases=10)
meth_tiles_var_allow_missing=unite(tiles_var_allow_missing, min.per.group=9L, destrand=TRUE)
save(meth_tiles_var_allow_missing, file = "meth_tiles_var_allow_missing.rda")

###
#load("meth_tiles_var_allow_missing.rda")
covar_allow_missingiates=data.frame(sex=c("female","female", "female","female","male","male","male","male","male","female","female","female","female","female","male","male","male","male","male","female"))
myDiff_var_allow_missing=calculateDiffMeth(meth_tiles_var_allow_missing, covar_allow_missingiates=covar_allow_missingiates)
write.table(myDiff_var_allow_missing, file="dmr_tiles_var_allow_missing_env.tsv", sep='\t')
my.diffMeth_10_var_allow_missing = getMethylDiff(myDiff_var_allow_missing,difference=10,qvalue=0.05)
write.table(my.diffMeth_10_var_allow_missing, file="dmr_tiles_var_allow_missing_10pourcent_env.tsv", sep='\t', quote=FALSE)

###########

### SEX EFFECT ####
#male=1 female=2
setwd("/PATH/TO/06_results/methylkit/")

myobj_var_allow_missing_sex=methRead(file.list.var,
                       sample.id.var,
                       assembly="unknown",
                       treatment=c(2,2,2,2,1,1,1,1,1,2,2,2,2,2,1,1,1,1,1,2),
                       context="CpG",
                       mincov=3)

setwd("/PATH/TO/06_results/methylkit/methylkit_var/")
#### Normalize
normalize.var_allow_missing_sex=normalizeCoverage(myobj_var_allow_missing_sex)
#test 10
filtered.myobj.10x.var_allow_missing_sex=filterByCoverage(normalize.var_allow_missing_sex,lo.count=10,lo.perc=NULL,
                                            hi.count=NULL,hi.perc=99.9)
save(filtered.myobj.10x.var_allow_missing_sex, file = "filtered.myobj.10x.var_allow_missing_sex.rda")
#### Merging samples
meth_var_allow_missing_sex=unite(myobj_var_allow_missing_sex, min.per.group=9L, destrand=TRUE)
save(meth_var_allow_missing_sex, file = "meth_var_allow_missing_sex.rda")
meth_filtered.myobj.10x.var_allow_missing_sex=unite(filtered.myobj.10x.var_allow_missing_sex, min.per.group=9L, destrand=TRUE)
save(meth_filtered.myobj.10x.var_allow_missing_sex, file = "meth_filtered.myobj.10x.var_allow_missing_sex.rda")

#### TILES 1000BP #################################################################################
### Cutting tiles
clusterSamples(meth_var_allow_missing_sex, dist="correlation", method="ward", plot=TRUE)

hist(meth_filtered.myobj.10x.var_allow_missing_sex$diff.meth)
filtered.myobj_var_allow_missing_sex=filterByCoverage(myobj_var_allow_missing_sex,lo.count=10,lo.perc=NULL,hi.count=NULL,hi.perc=99.9)
normalized.filtered.myobj_var_allow_missing_sex=normalizeCoverage(filtered.myobj_var_allow_missing_sex)
tiles_var_allow_missing_sex=tileMethylCounts(normalized.filtered.myobj_var_allow_missing_sex,win.size=1000,step.size=1000,cov.bases=10)
meth_tiles_var_allow_missing_sex=unite(tiles_var_allow_missing_sex, min.per.group=9L, destrand=TRUE)
save(meth_tiles_var_allow_missing_sex, file = "meth_tiles_var_allow_missing_sex.rda")

###
#load("meth_tiles_var_allow_missing.rda")

covar_allow_missingiates2=data.frame(habitat=c("rur","urb","urb","urb","urb","urb","urb","urb","urb","urb","rur","urb","rur","rur","rur","rur","rur","rur","rur","rur"))
myDiff_var_allow_missing_sex=calculateDiffMeth(meth_tiles_var_allow_missing_sex, covar_allow_missingiates=covar_allow_missingiates2)
write.table(myDiff_var_allow_missing_sex, file="dmr_tiles_var_allow_missing_sex.tsv", sep='\t')
my.diffMeth_10_var_allow_missing_sex = getMethylDiff(myDiff_var_allow_missing_sex,difference=10,qvalue=0.05)
write.table(my.diffMeth_10_var_allow_missing_sex, file="dmr_tiles_var_allow_missing_10pourcent_sex.tsv", sep='\t', quote=FALSE)

###########





