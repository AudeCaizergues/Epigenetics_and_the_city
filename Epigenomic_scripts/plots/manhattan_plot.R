#!/usr/bin/Rscript
################################################################################
###########       Script manhattan plot methylation        #####################
################################################################################

### libraries ##################################################################
library(qqman)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(ggExtra)
################################################################################

##### fabrication de la matrice pour le manhattan RDA ##########################
## FILE LOADING ####
setwd("/PATH/TO/06_results/rda")

bar.urb9prov<-read.table("bar9_S35_CpG.txt_filtered10cov.txt",header=T)
bar.urb9<-bar.urb9prov[,c("chr","base","freqT")]

bar.urb10prov<-read.table("bar10_S64_CpG.txt_filtered10cov.txt",header=T)
bar.urb10<-bar.urb10prov[,c("chr","base","freqT")]

bar.urb11prov<-read.table("bar11_S75_CpG.txt_filtered10cov.txt",header=T)
bar.urb11<-bar.urb11prov[,c("chr","base","freqT")]

bar.urb12prov<-read.table("bar12_S25_CpG.txt_filtered10cov.txt",header=T)
bar.urb12<-bar.urb12prov[,c("chr","base","freqT")]

bar.urb13prov<-read.table("bar13_S42_CpG.txt_filtered10cov.txt",header=T)
bar.urb13<-bar.urb13prov[,c("chr","base","freqT")]

bar.urb14prov<-read.table("bar14_S3_CpG.txt_filtered10cov.txt",header=T)
bar.urb14<-bar.urb14prov[,c("chr","base","freqT")]

bar.urb15prov<-read.table("bar15_S70_CpG.txt_filtered10cov.txt",header=T)
bar.urb15<-bar.urb15prov[,c("chr","base","freqT")]

bar.urb16prov<-read.table("bar16_S10_CpG.txt_filtered10cov.txt",header=T)
bar.urb16<-bar.urb16prov[,c("chr","base","freqT")]

bar.urb18prov<-read.table("bar18_S52_CpG.txt_filtered10cov.txt",header=T)
bar.urb18<-bar.urb18prov[,c("chr","base","freqT")]

bar.urb20prov<-read.table("bar20_S60_CpG.txt_filtered10cov.txt",header=T)
bar.urb20<-bar.urb20prov[,c("chr","base","freqT")]

bar.rur17prov<-read.table("bar17_S54_CpG.txt_filtered10cov.txt",header=T)
bar.rur17<-bar.rur17prov[,c("chr","base","freqT")]

bar.rur19prov<-read.table("bar19b_S12_CpG.txt_filtered10cov.txt",header=T)
bar.rur19<-bar.rur19prov[,c("chr","base","freqT")]

bar.rur1prov<-read.table("bar1_S46_CpG.txt_filtered10cov.txt",header=T)
bar.rur1<-bar.rur1prov[,c("chr","base","freqT")]

bar.rur2prov<-read.table("bar2_S23_CpG.txt_filtered10cov.txt",header=T)
bar.rur2<-bar.rur2prov[,c("chr","base","freqT")]

bar.rur3prov<-read.table("bar3_S29_CpG.txt_filtered10cov.txt",header=T)
bar.rur3<-bar.rur3prov[,c("chr","base","freqT")]

bar.rur4prov<-read.table("bar4_S63_CpG.txt_filtered10cov.txt",header=T)
bar.rur4<-bar.rur4prov[,c("chr","base","freqT")]

bar.rur5prov<-read.table("bar5_S45_CpG.txt_filtered10cov.txt",header=T)
bar.rur5<-bar.rur5prov[,c("chr","base","freqT")]

bar.rur6prov<-read.table("bar6_S7_CpG.txt_filtered10cov.txt",header=T)
bar.rur6<-bar.rur6prov[,c("chr","base","freqT")]

bar.rur7prov<-read.table("bar7_S16_CpG.txt_filtered10cov.txt",header=T)
bar.rur7<-bar.rur7prov[,c("chr","base","freqT")]

bar.rur8prov<-read.table("bar8_S59_CpG.txt_filtered10cov.txt",header=T)
bar.rur8<-bar.rur8prov[,c("chr","base","freqT")]

mtp.urb6prov<-read.table("mtp6_S78_CpG.txt_filtered10cov.txt",header=T)
mtp.urb6<-mtp.urb6prov[,c("chr","base","freqT")]

mtp.urb7prov<-read.table("mtp7_S44_CpG.txt_filtered10cov.txt",header=T)
mtp.urb7<-mtp.urb7prov[,c("chr","base","freqT")]

mtp.urb8prov<-read.table("mtp8_S27_CpG.txt_filtered10cov.txt",header=T)
mtp.urb8<-mtp.urb8prov[,c("chr","base","freqT")]

mtp.urb9prov<-read.table("mtp9_S43_CpG.txt_filtered10cov.txt",header=T)
mtp.urb9<-mtp.urb9prov[,c("chr","base","freqT")]

mtp.urb10prov<-read.table("mtp10_S19_CpG.txt_filtered10cov.txt",header=T)
mtp.urb10<-mtp.urb10prov[,c("chr","base","freqT")]

mtp.urb21prov<-read.table("mtp21_S6_CpG.txt_filtered10cov.txt",header=T)
mtp.urb21<-mtp.urb21prov[,c("chr","base","freqT")]

mtp.urb22prov<-read.table("mtp22_S21_CpG.txt_filtered10cov.txt",header=T)
mtp.urb22<-mtp.urb22prov[,c("chr","base","freqT")]

mtp.urb23prov<-read.table("mtp23_S37_CpG.txt_filtered10cov.txt",header=T)
mtp.urb23<-mtp.urb23prov[,c("chr","base","freqT")]

mtp.urb24prov<-read.table("mtp24_S32_CpG.txt_filtered10cov.txt",header=T)
mtp.urb24<-mtp.urb24prov[,c("chr","base","freqT")]

mtp.urb25prov<-read.table("mtp25_S38_CpG.txt_filtered10cov.txt",header=T)
mtp.urb25<-mtp.urb25prov[,c("chr","base","freqT")]

mtp.rur16prov<-read.table("mtp16b_S33_CpG.txt_filtered10cov.txt",header=T)
mtp.rur16<-mtp.rur16prov[,c("chr","base","freqT")]

mtp.rur17prov<-read.table("mtp17_S24_CpG.txt_filtered10cov.txt",header=T)
mtp.rur17<-mtp.rur17prov[,c("chr","base","freqT")]

mtp.rur18prov<-read.table("mtp18_S66_CpG.txt_filtered10cov.txt",header=T)
mtp.rur18<-mtp.rur18prov[,c("chr","base","freqT")]

mtp.rur19prov<-read.table("mtp19_S51_CpG.txt_filtered10cov.txt",header=T)
mtp.rur19<-mtp.rur19prov[,c("chr","base","freqT")]

mtp.rur20prov<-read.table("mtp20_S49_CpG.txt_filtered10cov.txt",header=T)
mtp.rur20<-mtp.rur20prov[,c("chr","base","freqT")]

mtp.rur1prov<-read.table("mtp1_S62_CpG.txt_filtered10cov.txt",header=T)
mtp.rur1<-mtp.rur1prov[,c("chr","base","freqT")]

mtp.rur2prov<-read.table("mtp2_S14_CpG.txt_filtered10cov.txt",header=T)
mtp.rur2<-mtp.rur2prov[,c("chr","base","freqT")]

mtp.rur3prov<-read.table("mtp3_S34_CpG.txt_filtered10cov.txt",header=T)
mtp.rur3<-mtp.rur3prov[,c("chr","base","freqT")]

mtp.rur4prov<-read.table("mtp4_S58_CpG.txt_filtered10cov.txt",header=T)
mtp.rur4<-mtp.rur4prov[,c("chr","base","freqT")]

mtp.rur5prov<-read.table("mtp5_S18_CpG.txt_filtered10cov.txt",header=T)
mtp.rur5<-mtp.rur5prov[,c("chr","base","freqT")]

var.urb11prov<-read.table("var11_S80_CpG.txt_filtered10cov.txt",header=T)
var.urb11<-var.urb11prov[,c("chr","base","freqT")]

var.urb12prov<-read.table("var12_S67_CpG.txt_filtered10cov.txt",header=T)
var.urb12<-var.urb12prov[,c("chr","base","freqT")]

var.urb13prov<-read.table("var13_S74_CpG.txt_filtered10cov.txt",header=T)
var.urb13<-var.urb13prov[,c("chr","base","freqT")]

var.urb14prov<-read.table("var14_S81_CpG.txt_filtered10cov.txt",header=T)
var.urb14<-var.urb14prov[,c("chr","base","freqT")]

var.urb15prov<-read.table("var15_S72_CpG.txt_filtered10cov.txt",header=T)
var.urb15<-var.urb15prov[,c("chr","base","freqT")]

var.urb16prov<-read.table("var16_S50_CpG.txt_filtered10cov.txt",header=T)
var.urb16<-var.urb16prov[,c("chr","base","freqT")]

var.urb17prov<-read.table("var17_S53_CpG.txt_filtered10cov.txt",header=T)
var.urb17<-var.urb17prov[,c("chr","base","freqT")]

var.urb18prov<-read.table("var18_S39_CpG.txt_filtered10cov.txt",header=T)
var.urb18<-var.urb18prov[,c("chr","base","freqT")]

var.urb19prov<-read.table("var19_S76_CpG.txt_filtered10cov.txt",header=T)
var.urb19<-var.urb19prov[,c("chr","base","freqT")]

var.urb20prov<-read.table("var20_S61_CpG.txt_filtered10cov.txt",header=T)
var.urb20<-var.urb20prov[,c("chr","base","freqT")]

var.rur1prov<-read.table("var1_S79_CpG.txt_filtered10cov.txt",header=T)
var.rur1<-var.rur1prov[,c("chr","base","freqT")]

var.rur2prov<-read.table("var2_S77_CpG.txt_filtered10cov.txt",header=T)
var.rur2<-var.rur2prov[,c("chr","base","freqT")]

var.rur3prov<-read.table("var3_S71_CpG.txt_filtered10cov.txt",header=T)
var.rur3<-var.rur3prov[,c("chr","base","freqT")]

var.rur4prov<-read.table("var4_S68_CpG.txt_filtered10cov.txt",header=T)
var.rur4<-var.rur4prov[,c("chr","base","freqT")]

var.rur5prov<-read.table("var5_S69_CpG.txt_filtered10cov.txt",header=T)
var.rur5<-var.rur5prov[,c("chr","base","freqT")]

var.rur6prov<-read.table("var6_S65_CpG.txt_filtered10cov.txt",header=T)
var.rur6<-var.rur6prov[,c("chr","base","freqT")]

var.rur7prov<-read.table("var7b_S30_CpG.txt_filtered10cov.txt",header=T)
var.rur7<-var.rur7prov[,c("chr","base","freqT")]

var.rur8prov<-read.table("var8_S31_CpG.txt_filtered10cov.txt",header=T)
var.rur8<-var.rur8prov[,c("chr","base","freqT")]

var.rur9prov<-read.table("var9_S41_CpG.txt_filtered10cov.txt",header=T)
var.rur9<-var.rur9prov[,c("chr","base","freqT")]

var.rur10prov<-read.table("var10_S8_CpG.txt_filtered10cov.txt",header=T)
var.rur10<-var.rur10prov[,c("chr","base","freqT")]
#############

### TOTAL###
#fusion de tous les individus avec uniquement les chrBases qui ont sont chez tout le monde
merge1<-merge(bar.urb9,bar.urb10, by=c("chr","base"),all=F)
merge2<-merge(merge1,bar.urb11, by=c("chr","base"),all=F)
merge3<-merge(merge2,bar.urb12, by=c("chr","base"),all=F)
merge4<-merge(merge3,bar.urb13, by=c("chr","base"),all=F)
merge5<-merge(merge4,bar.urb14, by=c("chr","base"),all=F)
merge6<-merge(merge5,bar.urb15, by=c("chr","base"),all=F)
merge7<-merge(merge6,bar.urb16, by=c("chr","base"),all=F)
merge8<-merge(merge7,bar.urb18, by=c("chr","base"),all=F)
merge9<-merge(merge8,bar.urb20, by=c("chr","base"),all=F)
merge10<-merge(merge9,bar.rur17, by=c("chr","base"),all=F)
merge11<-merge(merge10,bar.rur19, by=c("chr","base"),all=F)
merge12<-merge(merge11,bar.rur1, by=c("chr","base"),all=F)
merge13<-merge(merge12,bar.rur2, by=c("chr","base"),all=F)
merge14<-merge(merge13,bar.rur3, by=c("chr","base"),all=F)
merge15<-merge(merge14,bar.rur4, by=c("chr","base"),all=F)
merge16<-merge(merge15,bar.rur5, by=c("chr","base"),all=F)
merge17<-merge(merge16,bar.rur6, by=c("chr","base"),all=F)
merge18<-merge(merge17,bar.rur7, by=c("chr","base"),all=F)
merge19<-merge(merge18,bar.rur8, by=c("chr","base"),all=F)
merge20<-merge(merge19,mtp.urb6, by=c("chr","base"),all=F)
merge21<-merge(merge20,mtp.urb7, by=c("chr","base"),all=F)
merge22<-merge(merge21,mtp.urb8, by=c("chr","base"),all=F)
merge23<-merge(merge22,mtp.urb9, by=c("chr","base"),all=F)
merge24<-merge(merge23,mtp.urb10, by=c("chr","base"),all=F)
merge30<-merge(merge24,mtp.urb21, by=c("chr","base"),all=F)
merge31<-merge(merge30,mtp.urb22, by=c("chr","base"),all=F)
merge32<-merge(merge31,mtp.urb23, by=c("chr","base"),all=F)
merge33<-merge(merge32,mtp.urb24, by=c("chr","base"),all=F)
merge34<-merge(merge33,mtp.urb25, by=c("chr","base"),all=F)
merge40<-merge(merge34,mtp.rur16, by=c("chr","base"),all=F)
merge41<-merge(merge40,mtp.rur17, by=c("chr","base"),all=F)
merge42<-merge(merge41,mtp.rur18, by=c("chr","base"),all=F)
merge43<-merge(merge42,mtp.rur19, by=c("chr","base"),all=F)
merge44<-merge(merge43,mtp.rur20, by=c("chr","base"),all=F)
merge45<-merge(merge44,mtp.rur1, by=c("chr","base"),all=F)
merge46<-merge(merge45,mtp.rur2, by=c("chr","base"),all=F)
merge47<-merge(merge46,mtp.rur3, by=c("chr","base"),all=F)
merge48<-merge(merge47,mtp.rur4, by=c("chr","base"),all=F)
merge49<-merge(merge48,mtp.rur5, by=c("chr","base"),all=F)
merge50<-merge(merge49,var.urb11, by=c("chr","base"),all=F)
merge51<-merge(merge50,var.urb12, by=c("chr","base"),all=F)
merge52<-merge(merge51,var.urb13, by=c("chr","base"),all=F)
merge53<-merge(merge52,var.urb14, by=c("chr","base"),all=F)
merge54<-merge(merge53,var.urb15, by=c("chr","base"),all=F)
merge55<-merge(merge54,var.urb16, by=c("chr","base"),all=F)
merge56<-merge(merge55,var.urb17, by=c("chr","base"),all=F)
merge57<-merge(merge56,var.urb18, by=c("chr","base"),all=F)
merge58<-merge(merge57,var.urb19, by=c("chr","base"),all=F)
merge59<-merge(merge58,var.urb20, by=c("chr","base"),all=F)
merge60<-merge(merge59,var.rur1, by=c("chr","base"),all=F)
merge61<-merge(merge60,var.rur2, by=c("chr","base"),all=F)
merge62<-merge(merge61,var.rur3, by=c("chr","base"),all=F)
merge62<-merge(merge61,var.rur3, by=c("chr","base"),all=F)
merge63<-merge(merge62,var.rur4, by=c("chr","base"),all=F)
merge64<-merge(merge63,var.rur5, by=c("chr","base"),all=F)
merge65<-merge(merge64,var.rur6, by=c("chr","base"),all=F)
merge66<-merge(merge65,var.rur7, by=c("chr","base"),all=F)
merge67<-merge(merge66,var.rur8, by=c("chr","base"),all=F)
merge68<-merge(merge67,var.rur9, by=c("chr","base"),all=F)
merge69<-merge(merge68,var.rur10, by=c("chr","base"),all=F)
colnames(merge69)<-c("chr","pos","bar9.X","bar10.X","bar11.X","bar12.X","bar13.X"
                     ,"bar14.X","bar15.X","bar16.X","bar18.X","bar20.X"
                     ,"bar17.X","bar19.X","bar1.X","bar2.X","bar3.X"
                     ,"bar4.X","bar5.X","bar6.X","bar7.X","bar8.X"
                     ,"mtp6.X","mtp7.X","mtp8.X","mtp9.X","mtp10.X"
                     ,"mtp21.X","mtp22.X","mtp23.X","mtp24.X","mtp25.X"
                     ,"mtp16.X","mtp17.X","mtp18.X","mtp19.X","mtp20.X"
                     ,"mtp1.X","mtp2.X","mtp3.X","mtp4.X","mtp5.X"
                     ,"var11.X","var12.X","var13.X","var14.X","var15.X"
                     ,"var16.X","var17.X","var18.X","var19.X","var20.X"
                     ,"var1.X","var2.X","var3.X","var4.X","var5.X"
                     ,"var6.X","var7.X","var8.X","var9.X","var10.X")
merge69$chr<-as.character(merge69$chr)

#changement des noms de chromosomes
for (i in 1:nrow(merge69)) {
  print(i)
  if(merge69[i,1]=="NC_031768.1"){ 
    merge69[i,1]<-'01'
  }else if(merge69[i,1]=="NC_031769.1") {
    merge69[i,1]<-'02'
  }else if(merge69[i,1]=="NC_031770.1") {
    merge69[i,1]<-'03'
  }else if(merge69[i,1]=="NC_031771.1") {
    merge69[i,1]<-'04'
  }else if(merge69[i,1]=="NC_031772.1") {
    merge69[i,1]<-'04A'
  }else if(merge69[i,1]=="NC_031773.1") { 
    merge69[i,1]<-'01A'
  }else if(merge69[i,1]=="NC_031774.1") {
    merge69[i,1]<-'05'
  }else if(merge69[i,1]=="NC_031775.1") {
    merge69[i,1]<-'06'
  }else if(merge69[i,1]=="NC_031776.1") {
    merge69[i,1]<-'07'
  }else if(merge69[i,1]=="NC_031777.1") {
    merge69[i,1]<-'08'
  }else if(merge69[i,1]=="NC_031778.1") {
    merge69[i,1]<-'09'
  }else if(merge69[i,1]=="NC_031779.1") {
    merge69[i,1]<-'10'
  }else if(merge69[i,1]=="NC_031780.1") {
    merge69[i,1]<-'11'
  }else if(merge69[i,1]=="NC_031781.1") {
    merge69[i,1]<-'12'
  }else if(merge69[i,1]=="NC_031782.1") {
    merge69[i,1]<-'13'
  }else if(merge69[i,1]=="NC_031783.1") {
    merge69[i,1]<-'14'
  }else if(merge69[i,1]=="NC_031784.1") {
    merge69[i,1]<-'15'
  }else if(merge69[i,1]=="NC_031785.1") {
    merge69[i,1]<-'17'
  }else if(merge69[i,1]=="NC_031786.1") {
    merge69[i,1]<-'18'
  }else if(merge69[i,1]=="NC_031787.1") {
    merge69[i,1]<-'19'
  }else if(merge69[i,1]=="NC_031788.1") {
    merge69[i,1]<-'20'
  }else if(merge69[i,1]=="NC_031789.1") {
    merge69[i,1]<-'21'
  }else if(merge69[i,1]=="NC_031790.1") {
    merge69[i,1]<-'22'
  }else if(merge69[i,1]=="NC_031791.1") {
    merge69[i,1]<-'23'
  }else if(merge69[i,1]=="NC_031792.1") {
    merge69[i,1]<-'24'
  }else if(merge69[i,1]=="NC_031793.1") {
    merge69[i,1]<-'25LG1'
  }else if(merge69[i,1]=="NC_031794.1") {
    merge69[i,1]<-'25LG2'
  }else if(merge69[i,1]=="NC_031795.1") {
    merge69[i,1]<-'26'
  }else if(merge69[i,1]=="NC_031796.1") {
    merge69[i,1]<-'27'
  }else if(merge69[i,1]=="NC_031797.1") {
    merge69[i,1]<-'28'
  }else if(merge69[i,1]=="NC_031798.1") {
    merge69[i,1]<-'LGE22'
  }else if(merge69[i,1]=="NC_031799.1") {
    merge69[i,1]<-'Z'
  }else if(is.na(merge69[i,1])==TRUE) {
    merge69[i,1]<-'unplaced'
  }else {merge69[i,1]<-'unplaced'
  }
}

setwd("/media/aude/data_local/aude/RRBS_methylkit/06_results/manhattan_plot")
write.csv(merge69,file="merge_meth_par_position_par_indiv_pourmanhattan.txt", row.names = F, quote=F)
###############"

##### BARCELONE #####
#fusion de tous les individus avec uniquement les chrBases qui ont sont chez tout le monde
merge1<-merge(bar.urb9,bar.urb10, by=c("chr","base"),all=F)
merge2<-merge(merge1,bar.urb11, by=c("chr","base"),all=F)
merge3<-merge(merge2,bar.urb12, by=c("chr","base"),all=F)
merge4<-merge(merge3,bar.urb13, by=c("chr","base"),all=F)
merge5<-merge(merge4,bar.urb14, by=c("chr","base"),all=F)
merge6<-merge(merge5,bar.urb15, by=c("chr","base"),all=F)
merge7<-merge(merge6,bar.urb16, by=c("chr","base"),all=F)
merge8<-merge(merge7,bar.urb18, by=c("chr","base"),all=F)
merge9<-merge(merge8,bar.urb20, by=c("chr","base"),all=F)
merge10<-merge(merge9,bar.rur17, by=c("chr","base"),all=F)
merge11<-merge(merge10,bar.rur19, by=c("chr","base"),all=F)
merge12<-merge(merge11,bar.rur1, by=c("chr","base"),all=F)
merge13<-merge(merge12,bar.rur2, by=c("chr","base"),all=F)
merge14<-merge(merge13,bar.rur3, by=c("chr","base"),all=F)
merge15<-merge(merge14,bar.rur4, by=c("chr","base"),all=F)
merge16<-merge(merge15,bar.rur5, by=c("chr","base"),all=F)
merge17<-merge(merge16,bar.rur6, by=c("chr","base"),all=F)
merge18<-merge(merge17,bar.rur7, by=c("chr","base"),all=F)
mergebar<-merge(merge18,bar.rur8, by=c("chr","base"),all=F)
colnames(mergebar)<-c("chr","pos","bar9.X","bar10.X","bar11.X","bar12.X","bar13.X"
                      ,"bar14.X","bar15.X","bar16.X","bar18.X","bar20.X"
                      ,"bar17.X","bar19.X","bar1.X","bar2.X","bar3.X"
                      ,"bar4.X","bar5.X","bar6.X","bar7.X","bar8.X")
mergebar$chr<-as.character(mergebar$chr)


#changement des noms de chromosomes
for (i in 1:nrow(mergebar)) {
  print(i)
  if(mergebar[i,1]=="NC_031768.1"){ 
    mergebar[i,1]<-'01'
  }else if(mergebar[i,1]=="NC_031769.1") {
    mergebar[i,1]<-'02'
  }else if(mergebar[i,1]=="NC_031770.1") {
    mergebar[i,1]<-'03'
  }else if(mergebar[i,1]=="NC_031771.1") {
    mergebar[i,1]<-'04'
  }else if(mergebar[i,1]=="NC_031772.1") {
    mergebar[i,1]<-'04A'
  }else if(mergebar[i,1]=="NC_031773.1") { 
    mergebar[i,1]<-'01A'
  }else if(mergebar[i,1]=="NC_031774.1") {
    mergebar[i,1]<-'05'
  }else if(mergebar[i,1]=="NC_031775.1") {
    mergebar[i,1]<-'06'
  }else if(mergebar[i,1]=="NC_031776.1") {
    mergebar[i,1]<-'07'
  }else if(mergebar[i,1]=="NC_031777.1") {
    mergebar[i,1]<-'08'
  }else if(mergebar[i,1]=="NC_031778.1") {
    mergebar[i,1]<-'09'
  }else if(mergebar[i,1]=="NC_031779.1") {
    mergebar[i,1]<-'10'
  }else if(mergebar[i,1]=="NC_031780.1") {
    mergebar[i,1]<-'11'
  }else if(mergebar[i,1]=="NC_031781.1") {
    mergebar[i,1]<-'12'
  }else if(mergebar[i,1]=="NC_031782.1") {
    mergebar[i,1]<-'13'
  }else if(mergebar[i,1]=="NC_031783.1") {
    mergebar[i,1]<-'14'
  }else if(mergebar[i,1]=="NC_031784.1") {
    mergebar[i,1]<-'15'
  }else if(mergebar[i,1]=="NC_031785.1") {
    mergebar[i,1]<-'17'
  }else if(mergebar[i,1]=="NC_031786.1") {
    mergebar[i,1]<-'18'
  }else if(mergebar[i,1]=="NC_031787.1") {
    mergebar[i,1]<-'19'
  }else if(mergebar[i,1]=="NC_031788.1") {
    mergebar[i,1]<-'20'
  }else if(mergebar[i,1]=="NC_031789.1") {
    mergebar[i,1]<-'21'
  }else if(mergebar[i,1]=="NC_031790.1") {
    mergebar[i,1]<-'22'
  }else if(mergebar[i,1]=="NC_031791.1") {
    mergebar[i,1]<-'23'
  }else if(mergebar[i,1]=="NC_031792.1") {
    mergebar[i,1]<-'24'
  }else if(mergebar[i,1]=="NC_031793.1") {
    mergebar[i,1]<-'25LG1'
  }else if(mergebar[i,1]=="NC_031794.1") {
    mergebar[i,1]<-'25LG2'
  }else if(mergebar[i,1]=="NC_031795.1") {
    mergebar[i,1]<-'26'
  }else if(mergebar[i,1]=="NC_031796.1") {
    mergebar[i,1]<-'27'
  }else if(mergebar[i,1]=="NC_031797.1") {
    mergebar[i,1]<-'28'
  }else if(mergebar[i,1]=="NC_031798.1") {
    mergebar[i,1]<-'LGE22'
  }else if(mergebar[i,1]=="NC_031799.1") {
    mergebar[i,1]<-'Z'
  }else if(is.na(mergebar[i,1])==TRUE) {
    mergebar[i,1]<-'unplaced'
  }else {mergebar[i,1]<-'unplaced'
  }
}

setwd("/media/aude/data_local/aude/RRBS_methylkit/06_results/manhattan_plot")
write.csv(mergebar,file="merge_meth_par_position_par_indiv_pourmanhattan_bar.txt", row.names = F, quote=F)
########

#### MONTPELLIER ####
merge21<-merge(mtp.urb6,mtp.urb7, by=c("chr","base"),all=F)
merge22<-merge(merge21,mtp.urb8, by=c("chr","base"),all=F)
merge23<-merge(merge22,mtp.urb9, by=c("chr","base"),all=F)
merge24<-merge(merge23,mtp.urb10, by=c("chr","base"),all=F)
merge30<-merge(merge24,mtp.urb21, by=c("chr","base"),all=F)
merge31<-merge(merge30,mtp.urb22, by=c("chr","base"),all=F)
merge32<-merge(merge31,mtp.urb23, by=c("chr","base"),all=F)
merge33<-merge(merge32,mtp.urb24, by=c("chr","base"),all=F)
merge34<-merge(merge33,mtp.urb25, by=c("chr","base"),all=F)
merge40<-merge(merge34,mtp.rur16, by=c("chr","base"),all=F)
merge41<-merge(merge40,mtp.rur17, by=c("chr","base"),all=F)
merge42<-merge(merge41,mtp.rur18, by=c("chr","base"),all=F)
merge43<-merge(merge42,mtp.rur19, by=c("chr","base"),all=F)
merge44<-merge(merge43,mtp.rur20, by=c("chr","base"),all=F)
merge45<-merge(merge44,mtp.rur1, by=c("chr","base"),all=F)
merge46<-merge(merge45,mtp.rur2, by=c("chr","base"),all=F)
merge47<-merge(merge46,mtp.rur3, by=c("chr","base"),all=F)
merge48<-merge(merge47,mtp.rur4, by=c("chr","base"),all=F)
mergemtp<-merge(merge48,mtp.rur5, by=c("chr","base"),all=F)

colnames(mergemtp)<-c("chr","pos",
                      "mtp6.X","mtp7.X","mtp8.X","mtp9.X","mtp10.X"
                      ,"mtp21.X","mtp22.X","mtp23.X","mtp24.X","mtp25.X"
                      ,"mtp16.X","mtp17.X","mtp18.X","mtp19.X","mtp20.X"
                      ,"mtp1.X","mtp2.X","mtp3.X","mtp4.X","mtp5.X")
mergemtp$chr<-as.character(mergemtp$chr)


#changement des noms de chromosomes
for (i in 1:nrow(mergemtp)) {
  print(i)
  if(mergemtp[i,1]=="NC_031768.1"){ 
    mergemtp[i,1]<-'01'
  }else if(mergemtp[i,1]=="NC_031769.1") {
    mergemtp[i,1]<-'02'
  }else if(mergemtp[i,1]=="NC_031770.1") {
    mergemtp[i,1]<-'03'
  }else if(mergemtp[i,1]=="NC_031771.1") {
    mergemtp[i,1]<-'04'
  }else if(mergemtp[i,1]=="NC_031772.1") {
    mergemtp[i,1]<-'04A'
  }else if(mergemtp[i,1]=="NC_031773.1") { 
    mergemtp[i,1]<-'01A'
  }else if(mergemtp[i,1]=="NC_031774.1") {
    mergemtp[i,1]<-'05'
  }else if(mergemtp[i,1]=="NC_031775.1") {
    mergemtp[i,1]<-'06'
  }else if(mergemtp[i,1]=="NC_031776.1") {
    mergemtp[i,1]<-'07'
  }else if(mergemtp[i,1]=="NC_031777.1") {
    mergemtp[i,1]<-'08'
  }else if(mergemtp[i,1]=="NC_031778.1") {
    mergemtp[i,1]<-'09'
  }else if(mergemtp[i,1]=="NC_031779.1") {
    mergemtp[i,1]<-'10'
  }else if(mergemtp[i,1]=="NC_031780.1") {
    mergemtp[i,1]<-'11'
  }else if(mergemtp[i,1]=="NC_031781.1") {
    mergemtp[i,1]<-'12'
  }else if(mergemtp[i,1]=="NC_031782.1") {
    mergemtp[i,1]<-'13'
  }else if(mergemtp[i,1]=="NC_031783.1") {
    mergemtp[i,1]<-'14'
  }else if(mergemtp[i,1]=="NC_031784.1") {
    mergemtp[i,1]<-'15'
  }else if(mergemtp[i,1]=="NC_031785.1") {
    mergemtp[i,1]<-'17'
  }else if(mergemtp[i,1]=="NC_031786.1") {
    mergemtp[i,1]<-'18'
  }else if(mergemtp[i,1]=="NC_031787.1") {
    mergemtp[i,1]<-'19'
  }else if(mergemtp[i,1]=="NC_031788.1") {
    mergemtp[i,1]<-'20'
  }else if(mergemtp[i,1]=="NC_031789.1") {
    mergemtp[i,1]<-'21'
  }else if(mergemtp[i,1]=="NC_031790.1") {
    mergemtp[i,1]<-'22'
  }else if(mergemtp[i,1]=="NC_031791.1") {
    mergemtp[i,1]<-'23'
  }else if(mergemtp[i,1]=="NC_031792.1") {
    mergemtp[i,1]<-'24'
  }else if(mergemtp[i,1]=="NC_031793.1") {
    mergemtp[i,1]<-'25LG1'
  }else if(mergemtp[i,1]=="NC_031794.1") {
    mergemtp[i,1]<-'25LG2'
  }else if(mergemtp[i,1]=="NC_031795.1") {
    mergemtp[i,1]<-'26'
  }else if(mergemtp[i,1]=="NC_031796.1") {
    mergemtp[i,1]<-'27'
  }else if(mergemtp[i,1]=="NC_031797.1") {
    mergemtp[i,1]<-'28'
  }else if(mergemtp[i,1]=="NC_031798.1") {
    mergemtp[i,1]<-'LGE22'
  }else if(mergemtp[i,1]=="NC_031799.1") {
    mergemtp[i,1]<-'Z'
  }else if(is.na(mergemtp[i,1])==TRUE) {
    mergemtp[i,1]<-'unplaced'
  }else {mergemtp[i,1]<-'unplaced'
  }
}

setwd("/media/aude/data_local/aude/RRBS_methylkit/06_results/manhattan_plot")
write.csv(mergemtp,file="merge_meth_par_position_par_indiv_pourmanhattan_mtp.txt", row.names = F, quote=F)
######

#### VARSOVIE ####
merge50<-merge(merge49,var.urb11, by=c("chr","base"),all=F)
merge51<-merge(merge50,var.urb12, by=c("chr","base"),all=F)
merge52<-merge(merge51,var.urb13, by=c("chr","base"),all=F)
merge53<-merge(merge52,var.urb14, by=c("chr","base"),all=F)
merge54<-merge(merge53,var.urb15, by=c("chr","base"),all=F)
merge55<-merge(merge54,var.urb16, by=c("chr","base"),all=F)
merge56<-merge(merge55,var.urb17, by=c("chr","base"),all=F)
merge57<-merge(merge56,var.urb18, by=c("chr","base"),all=F)
merge58<-merge(merge57,var.urb19, by=c("chr","base"),all=F)
merge59<-merge(merge58,var.urb20, by=c("chr","base"),all=F)
merge60<-merge(merge59,var.rur1, by=c("chr","base"),all=F)
merge61<-merge(merge60,var.rur2, by=c("chr","base"),all=F)
merge62<-merge(merge61,var.rur3, by=c("chr","base"),all=F)
merge62<-merge(merge61,var.rur3, by=c("chr","base"),all=F)
merge63<-merge(merge62,var.rur4, by=c("chr","base"),all=F)
merge64<-merge(merge63,var.rur5, by=c("chr","base"),all=F)
merge65<-merge(merge64,var.rur6, by=c("chr","base"),all=F)
merge66<-merge(merge65,var.rur7, by=c("chr","base"),all=F)
merge67<-merge(merge66,var.rur8, by=c("chr","base"),all=F)
merge68<-merge(merge67,var.rur9, by=c("chr","base"),all=F)
mergevar<-merge(merge68,var.rur10, by=c("chr","base"),all=F)
colnames(mergevar)<-c("chr","pos"
                      ,"var11.X","var12.X","var13.X","var14.X","var15.X"
                      ,"var16.X","var17.X","var18.X","var19.X","var20.X"
                      ,"var1.X","var2.X","var3.X","var4.X","var5.X"
                      ,"var6.X","var7.X","var8.X","var9.X","var10.X")
mergevar$chr<-as.character(mergevar$chr)


#changement des noms de chromosomes
for (i in 1:nrow(mergevar)) {
  print(i)
  if(mergevar[i,1]=="NC_031768.1"){ 
    mergevar[i,1]<-'01'
  }else if(mergevar[i,1]=="NC_031769.1") {
    mergevar[i,1]<-'02'
  }else if(mergevar[i,1]=="NC_031770.1") {
    mergevar[i,1]<-'03'
  }else if(mergevar[i,1]=="NC_031771.1") {
    mergevar[i,1]<-'04'
  }else if(mergevar[i,1]=="NC_031772.1") {
    mergevar[i,1]<-'04A'
  }else if(mergevar[i,1]=="NC_031773.1") { 
    mergevar[i,1]<-'01A'
  }else if(mergevar[i,1]=="NC_031774.1") {
    mergevar[i,1]<-'05'
  }else if(mergevar[i,1]=="NC_031775.1") {
    mergevar[i,1]<-'06'
  }else if(mergevar[i,1]=="NC_031776.1") {
    mergevar[i,1]<-'07'
  }else if(mergevar[i,1]=="NC_031777.1") {
    mergevar[i,1]<-'08'
  }else if(mergevar[i,1]=="NC_031778.1") {
    mergevar[i,1]<-'09'
  }else if(mergevar[i,1]=="NC_031779.1") {
    mergevar[i,1]<-'10'
  }else if(mergevar[i,1]=="NC_031780.1") {
    mergevar[i,1]<-'11'
  }else if(mergevar[i,1]=="NC_031781.1") {
    mergevar[i,1]<-'12'
  }else if(mergevar[i,1]=="NC_031782.1") {
    mergevar[i,1]<-'13'
  }else if(mergevar[i,1]=="NC_031783.1") {
    mergevar[i,1]<-'14'
  }else if(mergevar[i,1]=="NC_031784.1") {
    mergevar[i,1]<-'15'
  }else if(mergevar[i,1]=="NC_031785.1") {
    mergevar[i,1]<-'17'
  }else if(mergevar[i,1]=="NC_031786.1") {
    mergevar[i,1]<-'18'
  }else if(mergevar[i,1]=="NC_031787.1") {
    mergevar[i,1]<-'19'
  }else if(mergevar[i,1]=="NC_031788.1") {
    mergevar[i,1]<-'20'
  }else if(mergevar[i,1]=="NC_031789.1") {
    mergevar[i,1]<-'21'
  }else if(mergevar[i,1]=="NC_031790.1") {
    mergevar[i,1]<-'22'
  }else if(mergevar[i,1]=="NC_031791.1") {
    mergevar[i,1]<-'23'
  }else if(mergevar[i,1]=="NC_031792.1") {
    mergevar[i,1]<-'24'
  }else if(mergevar[i,1]=="NC_031793.1") {
    mergevar[i,1]<-'25LG1'
  }else if(mergevar[i,1]=="NC_031794.1") {
    mergevar[i,1]<-'25LG2'
  }else if(mergevar[i,1]=="NC_031795.1") {
    mergevar[i,1]<-'26'
  }else if(mergevar[i,1]=="NC_031796.1") {
    mergevar[i,1]<-'27'
  }else if(mergevar[i,1]=="NC_031797.1") {
    mergevar[i,1]<-'28'
  }else if(mergevar[i,1]=="NC_031798.1") {
    mergevar[i,1]<-'LGE22'
  }else if(mergevar[i,1]=="NC_031799.1") {
    mergevar[i,1]<-'Z'
  }else if(is.na(mergevar[i,1])==TRUE) {
    mergevar[i,1]<-'unplaced'
  }else {mergevar[i,1]<-'unplaced'
  }
}

setwd("/media/aude/data_local/aude/RRBS_methylkit/06_results/manhattan_plot")
write.csv(mergevar,file="merge_meth_par_position_par_indiv_pourmanhattan_var.txt", row.names = F, quote=F)
######

#########################################################################################

#### PREPARATION DATA ###################################################################
#on calcule la methylation moyenne a chaque position par population
data<-read.csv("merge_meth_par_position_par_indiv_pourmanhattan.txt")
data.means.milieu<-data.frame(chr=data[,1],pos=data[,2], 
                              mean.bar.urb=rowMeans(data[,3:12]), 
                              mean.bar.rur=rowMeans(data[,13:22]),
                              mean.mtp.urb=rowMeans(data[,23:32]), 
                              mean.mtp.rur=rowMeans(data[,33:42]),
                              mean.var.urb=rowMeans(data[,43:52]), 
                              mean.var.rur=rowMeans(data[,53:62]))   

# on calcule les differences de methylation moyennes entre rur et urb a chaque position par ville
data.means.milieu$diff.bar<-(data.means.milieu$mean.bar.rur-data.means.milieu$mean.bar.urb)
data.means.milieu$diff.mtp<-(data.means.milieu$mean.mtp.rur-data.means.milieu$mean.mtp.urb)
data.means.milieu$diff.var<-(data.means.milieu$mean.var.rur-data.means.milieu$mean.var.urb)
data.means.milieu$diff.rururb<-((data.means.milieu$diff.bar+data.means.milieu$diff.mtp+data.means.milieu$diff.var)/3)

data.means.milieu$diff.rururb<-as.numeric(data.means.milieu$diff.rururb)
data.means.milieu$pos<-as.numeric(data.means.milieu$pos)

#### PREPARE MANHATTAN PLOT #############################################################
### compute cumulative position of SNP
don <- data.means.milieu %>% 
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(pos)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(data.means.milieu, ., by=c("chr"="chr")) %>%
  # Add a cumulative position of each SNP
  arrange(chr, pos) %>%
  mutate( BPcum=pos+tot)

### prepare axis to show only chromosome name
axisdf = don %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

#### MANHATTAN PLOT GLOBAL ########################################################################
manhattan.glob<-ggplot(don, aes(x=BPcum, y=diff.rururb)) +
  # Show all points
  geom_point( aes(color=as.factor(chr)), alpha=0.2, size=0.8) +
  scale_color_manual(values = rep(c("grey", "black"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center) +
  #scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  labs(x = "Position on chromosomes", y="Change in % methylation in urban versus rural") +
  ylim(-40,40)+
  ggtitle("Global")+
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text = element_text(angle=90, hjust=1)
  )

manhattan.glob


freq.glob <-ggplot()+
  aes(data.means.milieu$diff.rururb) +
  geom_histogram(binwidth=0.5, colour="black", fill="black")+
  coord_flip()+
  theme_bw()+
  xlim(-40,40)+
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text = element_text(angle=90, hjust=1),
    axis.text.y = element_blank(),
    axis.title.y=element_blank(),
    axis.ticks.y=element_blank()
  )

freq.glob

freq.glob2<-ggdensity(data.means.milieu,data.means.milieu$diff.rururb)

ggMarginal(manhattan.glob, type= "density", margins = "y",size=6)


#### MANHATTAN PLOT GLOBAL 5% de difference minimum###############################
data.means.milieu.5pourcent<-subset(data.means.milieu,data.means.milieu$diff.rururb>=5 | data.means.milieu$diff.rururb<=-5)
### compute cumulative position of SNP
don <- data.means.milieu.5pourcent %>% 
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(pos)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(data.means.milieu.5pourcent, ., by=c("chr"="chr")) %>%
  # Add a cumulative position of each SNP
  arrange(chr, pos) %>%
  mutate( BPcum=pos+tot)

### prepare axis to show only chromosome name
axisdf = don %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

### plot with ggplot
ggplot(don, aes(x=BPcum, y=diff.rururb)) +
  
  # Show all points
  geom_point( aes(color=as.factor(chr)), alpha=0.8, size=0.8) +
  scale_color_manual(values = rep(c("grey", "black"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
  #scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  labs(x = "Chromosome", y="Change in % methylation in urban versus rural \n 5% minimum") +
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text = element_text(angle=90, hjust=1)
    
  )
#### MANHATTAN PLOT PAR VILLE #####################################################################
### BARCELONE
par(mfrow=c(1,1))
ggplot(don, aes(x=BPcum, y=diff.bar)) +
  # Show all points
  geom_point( aes(color=as.factor(chr)), alpha=0.8, size=0.8) +
  scale_color_manual(values = rep(c("indianred4", "indianred1"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
  #scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  labs(x = "Position on chromosomes", y="Change in % methylation in urban versus rural") +
  ylim(-40,40)+
  ggtitle("Barcelone")+
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text = element_text(angle=90, hjust=1)
  )

### MONTPELLIER 
ggplot(don, aes(x=BPcum, y=diff.mtp)) +
  
  # Show all points
  geom_point( aes(color=as.factor(chr)), alpha=0.8, size=0.8) +
  scale_color_manual(values = rep(c("chartreuse4", "olivedrab3"), 22 )) +
  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
  #scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  labs(x = "Position on chromosomes", y="Change in % methylation in urban versus rural") +
  ylim(-40,40)+
  ggtitle("Montpellier")+
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text = element_text(angle=90, hjust=1)
  )

### VARSOVIE
ggplot(don, aes(x=BPcum, y=diff.var)) +
  # Show all points
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=0.8) +
  scale_color_manual(values = rep(c("steelblue4", "steelblue1"), 22 )) +
  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
  #scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  labs(x = "Position on chromosomes", y="Change in % methylation in urban versus rural") +
  ylim(-40,40)+
  ggtitle("Varsovie")+
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text = element_text(angle=90, hjust=1)
  )



### BOXPLOT METH MOY/ VILLE / CHROM ##################################################
# subset BARCELONE
data.bar<-data[,(1:22)]
names(data.bar) # les 10 premiers sont urbains, les 10 derniers ruraux
data.means.bar<-data.frame(chr=data.bar[,1],pos=data.bar[,2], 
                           mean.urb=rowMeans(data.bar[,3:12]), 
                           mean.rur=rowMeans(data.bar[,13:22]))   
data.means.bar$ville<-c(rep("bar",length(data.means.bar$pos)))
# subset MONTPELLIER
data.mtp<-data[,c(1,2,23:52)]
names(data.mtp) # les 20 premiers sont urbains, les 10 derniers ruraux
data.means.mtp<-data.frame(chr=data.mtp[,1],pos=data.mtp[,2], 
                           mean.urb=rowMeans(data.mtp[,3:22]), 
                           mean.rur=rowMeans(data.mtp[,23:32])) 
data.means.mtp$ville<-c(rep("mtp",length(data.means.mtp$pos)))
# subset VARSOVIE
data.var<-data[,c(1,2,53:72)]
names(data.var) # les 10 premiers sont urbains, les 10 derniers ruraux
data.means.var<-data.frame(chr=data.var[,1],pos=data.var[,2], 
                           mean.urb=rowMeans(data.var[,3:12]), 
                           mean.rur=rowMeans(data.var[,13:22]))
data.means.var$ville<-c(rep("var",length(data.means.var$pos)))

# fusion des 3 tables
mean_pour_boxplot<-rbind(data.means.bar,data.means.mtp,data.means.var)
mean_pour_boxplot$diff.rururb<-(mean_pour_boxplot$mean.rur-mean_pour_boxplot$mean.urb)*100

# boxplot général
ggplot(mean_pour_boxplot, aes(x=chr, y=diff.rururb, fill=ville)) + 
  geom_boxplot()

### BOXPLOT pour les diff >10% ou <10%
meanpos_pour_boxplot<-subset(mean_pour_boxplot, mean_pour_boxplot$diff.rururb>=10)
meanneg_pour_boxplot<-subset(mean_pour_boxplot, mean_pour_boxplot$diff.rururb<=-10)

ggplot(meanpos_pour_boxplot, aes(x=chr, y=diff.rururb, fill=ville)) + 
  labs(x = "Chromosome", y="Change in % methylation in urban versus rural \n 10% minimum") +
  geom_boxplot()

ggplot(meanneg_pour_boxplot, aes(x=chr, y=diff.rururb, fill=ville)) + 
  labs(x = "Chromosome", y="Change in % methylation in urban versus rural \n 10% minimum") +
  geom_boxplot()
