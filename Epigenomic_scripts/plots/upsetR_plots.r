########################################################################
#############            UpsetR diagram urban         ##################
########################################################################

### LIBRARIES ##########################################################
library(anchors)
library(UpSetR)
library(VennDiagram)
library(tidyverse)
library(hrbrthemes)
library(tm)
library(proustr)
########################################################################


### UPSET DIAGRAM A PARTIR DES DMR VILLE PAR VILLE #####################
setwd("/PATH/TO/06_results/methylkit/methylkit_bar")
dmr_bar_env_signif<-read.table("dmr_tiles_bar_allow_missing_10pourcent_env.tsv",header=TRUE)
setwd("/PATH/TO/06_results/methylkit/methylkit_mtp")
dmr_mtp_env_signif<-read.table("dmr_tiles_mtp_allow_missing_10pourcent_env.tsv",header=TRUE)
setwd("/PATH/TO/06_results/methylkit/methylkit_var")
dmr_var_env_signif<-read.table("dmr_tiles_var_allow_missing_10pourcent_env.tsv",header=TRUE)

setwd("/PATH/TO/06_results/upsetR_venn_diagrams")

# dmr_bar_env_signif$DMR <- paste(dmr_bar_env_signif$chr, "-", dmr_bar_env_signif$start)
# dmr_mtp_env_signif$DMR <- paste(dmr_mtp_env_signif$chr, "-", dmr_mtp_env_signif$start)
# dmr_var_env_signif$DMR <- paste(dmr_var_env_signif$chr, "-", dmr_var_env_signif$start)

merg1<-merge(dmr_bar_env_signif,dmr_mtp_env_signif, by=c("chr","start","end"),all=TRUE, suffixes=c("_bar","_mtp"))
merg2<-merge(merg1,dmr_var_env_signif, by=c("chr","start","end"),all=TRUE)
matrice_DMR<-merg2[,c("chr", "start", "end","meth.diff_bar", "meth.diff_mtp","meth.diff")] # var n'a pas de suffixe
colnames(matrice_DMR)<-c("chr", "start", "end", "meth.diff_bar","meth.diff_mtp", "meth.diff_var")

matrice_DMR$meth.diff_bar[!is.na(matrice_DMR$meth.diff_bar)]<-1
matrice_DMR$meth.diff_bar[is.na(matrice_DMR$meth.diff_bar)]<-0
matrice_DMR$meth.diff_mtp[!is.na(matrice_DMR$meth.diff_mtp)]<-1
matrice_DMR$meth.diff_mtp[is.na(matrice_DMR$meth.diff_mtp)]<-0
matrice_DMR$meth.diff_var[!is.na(matrice_DMR$meth.diff_var)]<-1
matrice_DMR$meth.diff_var[is.na(matrice_DMR$meth.diff_var)]<-0

write.table(matrice_DMR, file="matrix_DMR_10percent_all_presence-absence.txt", sep='\t', quote=FALSE)

upset<-upset(matrice_DMR, sets = c("meth.diff_bar","meth.diff_mtp","meth.diff_var"), keep.order= TRUE, number.angles = 30, point.size = 3.5, line.size = 2, 
                             mainbar.y.label = "Number of common DMR \n between urban and rural habitat \n min 10% methylation difference", sets.x.label = "Nb DMR per model", 
                             text.scale = c(1.3, 1.3, 1, 1, 1.3, 1.3))
upset


matrice_DMR$sum<-(matrice_DMR$meth.diff_bar+matrice_DMR$meth.diff_mtp+matrice_DMR$meth.diff_var)
matrice_DMR$sum

intersection<-subset(matrice_DMR, matrice_DMR$sum<=2)
write.table(intersection, file="matrix_DMR_10percent_intersection2ou3villes.txt", sep='\t', quote=FALSE,row.names=FALSE)

intersection2<-subset(matrice_DMR, matrice_DMR$sum==2)
write.table(intersection2, file="matrix_DMR_10percent_intersection2villes.txt", sep='\t', quote=FALSE,row.names=FALSE)

intersection3<-subset(matrice_DMR, matrice_DMR$sum==3)
write.table(intersection3, file="matrix_DMR_10percent_intersection3villes.txt", sep='\t', quote=FALSE,row.names=FALSE)

########################################################################
