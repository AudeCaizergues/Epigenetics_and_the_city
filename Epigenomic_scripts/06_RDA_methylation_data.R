#!/usr/bin/Rscript
### PURPOSE OF THE SCRIPT: Performing a RDA on methylation data

# A GOOD TUTORIAL TO HELP UNDERSTAND RDA ANALYSES
#https://rstudio-pubs-static.s3.amazonaws.com/64619_2f93b223a318410bbf999d092ecf05ec.html

########################################################################
#############	   	 RDA on methylation data 		  ##################
########################################################################

ls()
rm(list=ls())
ls()

### set working directory
setwd("PATH/TO/FILTERED_FILES/06_results/rda")

##### BUILDING MATRIX FOR RDA ANALYSIS #################################
### Input files from 05_filtre 05_coverage_filtration_for_RDA.sh
### freqT = proportion of methylation per base

bar.urb9prov<-read.table("bar9_S35_CpG.txt_filtered10cov.txt",header=T)
bar.urb9<-bar.urb9prov[,c("chrBase","freqT")]
bar.urb10prov<-read.table("bar10_S64_CpG.txt_filtered10cov.txt",header=T)
bar.urb10<-bar.urb10prov[,c("chrBase","freqT")]
bar.urb11prov<-read.table("bar11_S75_CpG.txt_filtered10cov.txt",header=T)
bar.urb11<-bar.urb11prov[,c("chrBase","freqT")]
bar.urb12prov<-read.table("bar12_S25_CpG.txt_filtered10cov.txt",header=T)
bar.urb12<-bar.urb12prov[,c("chrBase","freqT")]
bar.urb13prov<-read.table("bar13_S42_CpG.txt_filtered10cov.txt",header=T)
bar.urb13<-bar.urb13prov[,c("chrBase","freqT")]
bar.urb14prov<-read.table("bar14_S3_CpG.txt_filtered10cov.txt",header=T)
bar.urb14<-bar.urb14prov[,c("chrBase","freqT")]
bar.urb15prov<-read.table("bar15_S70_CpG.txt_filtered10cov.txt",header=T)
bar.urb15<-bar.urb15prov[,c("chrBase","freqT")]
bar.urb16prov<-read.table("bar16_S10_CpG.txt_filtered10cov.txt",header=T)
bar.urb16<-bar.urb16prov[,c("chrBase","freqT")]
bar.urb18prov<-read.table("bar18_S52_CpG.txt_filtered10cov.txt",header=T)
bar.urb18<-bar.urb18prov[,c("chrBase","freqT")]
bar.urb20prov<-read.table("bar20_S60_CpG.txt_filtered10cov.txt",header=T)
bar.urb20<-bar.urb20prov[,c("chrBase","freqT")]
bar.rur17prov<-read.table("bar17_S54_CpG.txt_filtered10cov.txt",header=T)
bar.rur17<-bar.rur17prov[,c("chrBase","freqT")]
bar.rur19prov<-read.table("bar19b_S12_CpG.txt_filtered10cov.txt",header=T)
bar.rur19<-bar.rur19prov[,c("chrBase","freqT")]
bar.rur1prov<-read.table("bar1_S46_CpG.txt_filtered10cov.txt",header=T)
bar.rur1<-bar.rur1prov[,c("chrBase","freqT")]
bar.rur2prov<-read.table("bar2_S23_CpG.txt_filtered10cov.txt",header=T)
bar.rur2<-bar.rur2prov[,c("chrBase","freqT")]
bar.rur3prov<-read.table("bar3_S29_CpG.txt_filtered10cov.txt",header=T)
bar.rur3<-bar.rur3prov[,c("chrBase","freqT")]
bar.rur4prov<-read.table("bar4_S63_CpG.txt_filtered10cov.txt",header=T)
bar.rur4<-bar.rur4prov[,c("chrBase","freqT")]
bar.rur5prov<-read.table("bar5_S45_CpG.txt_filtered10cov.txt",header=T)
bar.rur5<-bar.rur5prov[,c("chrBase","freqT")]
bar.rur6prov<-read.table("bar6_S7_CpG.txt_filtered10cov.txt",header=T)
bar.rur6<-bar.rur6prov[,c("chrBase","freqT")]
bar.rur7prov<-read.table("bar7_S16_CpG.txt_filtered10cov.txt",header=T)
bar.rur7<-bar.rur7prov[,c("chrBase","freqT")]
bar.rur8prov<-read.table("bar8_S59_CpG.txt_filtered10cov.txt",header=T)
bar.rur8<-bar.rur8prov[,c("chrBase","freqT")]
mtp.urb6prov<-read.table("mtp6_S78_CpG.txt_filtered10cov.txt",header=T)
mtp.urb6<-mtp.urb6prov[,c("chrBase","freqT")]
mtp.urb7prov<-read.table("mtp7_S44_CpG.txt_filtered10cov.txt",header=T)
mtp.urb7<-mtp.urb7prov[,c("chrBase","freqT")]
mtp.urb8prov<-read.table("mtp8_S27_CpG.txt_filtered10cov.txt",header=T)
mtp.urb8<-mtp.urb8prov[,c("chrBase","freqT")]
mtp.urb9prov<-read.table("mtp9_S43_CpG.txt_filtered10cov.txt",header=T)
mtp.urb9<-mtp.urb9prov[,c("chrBase","freqT")]
mtp.urb10prov<-read.table("mtp10_S19_CpG.txt_filtered10cov.txt",header=T)
mtp.urb10<-mtp.urb10prov[,c("chrBase","freqT")]
mtp.urb21prov<-read.table("mtp21_S6_CpG.txt_filtered10cov.txt",header=T)
mtp.urb21<-mtp.urb21prov[,c("chrBase","freqT")]
mtp.urb22prov<-read.table("mtp22_S21_CpG.txt_filtered10cov.txt",header=T)
mtp.urb22<-mtp.urb22prov[,c("chrBase","freqT")]
mtp.urb23prov<-read.table("mtp23_S37_CpG.txt_filtered10cov.txt",header=T)
mtp.urb23<-mtp.urb23prov[,c("chrBase","freqT")]
mtp.urb24prov<-read.table("mtp24_S32_CpG.txt_filtered10cov.txt",header=T)
mtp.urb24<-mtp.urb24prov[,c("chrBase","freqT")]
mtp.urb25prov<-read.table("mtp25_S38_CpG.txt_filtered10cov.txt",header=T)
mtp.urb25<-mtp.urb25prov[,c("chrBase","freqT")]
mtp.rur16prov<-read.table("mtp16b_S33_CpG.txt_filtered10cov.txt",header=T)
mtp.rur16<-mtp.rur16prov[,c("chrBase","freqT")]
mtp.rur17prov<-read.table("mtp17_S24_CpG.txt_filtered10cov.txt",header=T)
mtp.rur17<-mtp.rur17prov[,c("chrBase","freqT")]
mtp.rur18prov<-read.table("mtp18_S66_CpG.txt_filtered10cov.txt",header=T)
mtp.rur18<-mtp.rur18prov[,c("chrBase","freqT")]
mtp.rur19prov<-read.table("mtp19_S51_CpG.txt_filtered10cov.txt",header=T)
mtp.rur19<-mtp.rur19prov[,c("chrBase","freqT")]
mtp.rur20prov<-read.table("mtp20_S49_CpG.txt_filtered10cov.txt",header=T)
mtp.rur20<-mtp.rur20prov[,c("chrBase","freqT")]
mtp.rur1prov<-read.table("mtp1_S62_CpG.txt_filtered10cov.txt",header=T)
mtp.rur1<-mtp.rur1prov[,c("chrBase","freqT")]
mtp.rur2prov<-read.table("mtp2_S14_CpG.txt_filtered10cov.txt",header=T)
mtp.rur2<-mtp.rur2prov[,c("chrBase","freqT")]
mtp.rur3prov<-read.table("mtp3_S34_CpG.txt_filtered10cov.txt",header=T)
mtp.rur3<-mtp.rur3prov[,c("chrBase","freqT")]
mtp.rur4prov<-read.table("mtp4_S58_CpG.txt_filtered10cov.txt",header=T)
mtp.rur4<-mtp.rur4prov[,c("chrBase","freqT")]
mtp.rur5prov<-read.table("mtp5_S18_CpG.txt_filtered10cov.txt",header=T)
mtp.rur5<-mtp.rur5prov[,c("chrBase","freqT")]
var.urb11prov<-read.table("var11_S80_CpG.txt_filtered10cov.txt",header=T)
var.urb11<-var.urb11prov[,c("chrBase","freqT")]
var.urb12prov<-read.table("var12_S67_CpG.txt_filtered10cov.txt",header=T)
var.urb12<-var.urb12prov[,c("chrBase","freqT")]
var.urb13prov<-read.table("var13_S74_CpG.txt_filtered10cov.txt",header=T)
var.urb13<-var.urb13prov[,c("chrBase","freqT")]
var.urb14prov<-read.table("var14_S81_CpG.txt_filtered10cov.txt",header=T)
var.urb14<-var.urb14prov[,c("chrBase","freqT")]
var.urb15prov<-read.table("var15_S72_CpG.txt_filtered10cov.txt",header=T)
var.urb15<-var.urb15prov[,c("chrBase","freqT")]
var.urb16prov<-read.table("var16_S50_CpG.txt_filtered10cov.txt",header=T)
var.urb16<-var.urb16prov[,c("chrBase","freqT")]
var.urb17prov<-read.table("var17_S53_CpG.txt_filtered10cov.txt",header=T)
var.urb17<-var.urb17prov[,c("chrBase","freqT")]
var.urb18prov<-read.table("var18_S39_CpG.txt_filtered10cov.txt",header=T)
var.urb18<-var.urb18prov[,c("chrBase","freqT")]
var.urb19prov<-read.table("var19_S76_CpG.txt_filtered10cov.txt",header=T)
var.urb19<-var.urb19prov[,c("chrBase","freqT")]
var.urb20prov<-read.table("var20_S61_CpG.txt_filtered10cov.txt",header=T)
var.urb20<-var.urb20prov[,c("chrBase","freqT")]
var.rur1prov<-read.table("var1_S79_CpG.txt_filtered10cov.txt",header=T)
var.rur1<-var.rur1prov[,c("chrBase","freqT")]
var.rur2prov<-read.table("var2_S77_CpG.txt_filtered10cov.txt",header=T)
var.rur2<-var.rur2prov[,c("chrBase","freqT")]
var.rur3prov<-read.table("var3_S71_CpG.txt_filtered10cov.txt",header=T)
var.rur3<-var.rur3prov[,c("chrBase","freqT")]
var.rur4prov<-read.table("var4_S68_CpG.txt_filtered10cov.txt",header=T)
var.rur4<-var.rur4prov[,c("chrBase","freqT")]
var.rur5prov<-read.table("var5_S69_CpG.txt_filtered10cov.txt",header=T)
var.rur5<-var.rur5prov[,c("chrBase","freqT")]
var.rur6prov<-read.table("var6_S65_CpG.txt_filtered10cov.txt",header=T)
var.rur6<-var.rur6prov[,c("chrBase","freqT")]
var.rur7prov<-read.table("var7b_S30_CpG.txt_filtered10cov.txt",header=T)
var.rur7<-var.rur7prov[,c("chrBase","freqT")]
var.rur8prov<-read.table("var8_S31_CpG.txt_filtered10cov.txt",header=T)
var.rur8<-var.rur8prov[,c("chrBase","freqT")]
var.rur9prov<-read.table("var9_S41_CpG.txt_filtered10cov.txt",header=T)
var.rur9<-var.rur9prov[,c("chrBase","freqT")]
var.rur10prov<-read.table("var10_S8_CpG.txt_filtered10cov.txt",header=T)
var.rur10<-var.rur10prov[,c("chrBase","freqT")]

## Merging by positions that are present in every individual
merge1<-merge(bar.urb9,bar.urb10, by=c('chrBase'),all=F)
merge2<-merge(merge1,bar.urb11, by=c('chrBase'),all=F)
merge3<-merge(merge2,bar.urb12, by=c('chrBase'),all=F)
merge4<-merge(merge3,bar.urb13, by=c('chrBase'),all=F)
merge5<-merge(merge4,bar.urb14, by=c('chrBase'),all=F)
merge6<-merge(merge5,bar.urb15, by=c('chrBase'),all=F)
merge7<-merge(merge6,bar.urb16, by=c('chrBase'),all=F)
merge8<-merge(merge7,bar.urb18, by=c('chrBase'),all=F)
merge9<-merge(merge8,bar.urb20, by=c('chrBase'),all=F)
merge10<-merge(merge9,bar.rur17, by=c('chrBase'),all=F)
merge11<-merge(merge10,bar.rur19, by=c('chrBase'),all=F)
merge12<-merge(merge11,bar.rur1, by=c('chrBase'),all=F)
merge13<-merge(merge12,bar.rur2, by=c('chrBase'),all=F)
merge14<-merge(merge13,bar.rur3, by=c('chrBase'),all=F)
merge15<-merge(merge14,bar.rur4, by=c('chrBase'),all=F)
merge16<-merge(merge15,bar.rur5, by=c('chrBase'),all=F)
merge17<-merge(merge16,bar.rur6, by=c('chrBase'),all=F)
merge18<-merge(merge17,bar.rur7, by=c('chrBase'),all=F)
merge19<-merge(merge18,bar.rur8, by=c('chrBase'),all=F)
merge20<-merge(merge19,mtp.urb6, by=c('chrBase'),all=F)
merge21<-merge(merge20,mtp.urb7, by=c('chrBase'),all=F)
merge22<-merge(merge21,mtp.urb8, by=c('chrBase'),all=F)
merge23<-merge(merge22,mtp.urb9, by=c('chrBase'),all=F)
merge24<-merge(merge23,mtp.urb10, by=c('chrBase'),all=F)
merge30<-merge(merge24,mtp.urb21, by=c('chrBase'),all=F)
merge31<-merge(merge30,mtp.urb22, by=c('chrBase'),all=F)
merge32<-merge(merge31,mtp.urb23, by=c('chrBase'),all=F)
merge33<-merge(merge32,mtp.urb24, by=c('chrBase'),all=F)
merge34<-merge(merge33,mtp.urb25, by=c('chrBase'),all=F)
merge40<-merge(merge34,mtp.rur16, by=c('chrBase'),all=F)
merge41<-merge(merge40,mtp.rur17, by=c('chrBase'),all=F)
merge42<-merge(merge41,mtp.rur18, by=c('chrBase'),all=F)
merge43<-merge(merge42,mtp.rur19, by=c('chrBase'),all=F)
merge44<-merge(merge43,mtp.rur20, by=c('chrBase'),all=F)
merge45<-merge(merge44,mtp.rur1, by=c('chrBase'),all=F)
merge46<-merge(merge45,mtp.rur2, by=c('chrBase'),all=F)
merge47<-merge(merge46,mtp.rur3, by=c('chrBase'),all=F)
merge48<-merge(merge47,mtp.rur4, by=c('chrBase'),all=F)
merge49<-merge(merge48,mtp.rur5, by=c('chrBase'),all=F)
merge50<-merge(merge49,var.urb11, by=c('chrBase'),all=F)
merge51<-merge(merge50,var.urb12, by=c('chrBase'),all=F)
merge52<-merge(merge51,var.urb13, by=c('chrBase'),all=F)
merge53<-merge(merge52,var.urb14, by=c('chrBase'),all=F)
merge54<-merge(merge53,var.urb15, by=c('chrBase'),all=F)
merge55<-merge(merge54,var.urb16, by=c('chrBase'),all=F)
merge56<-merge(merge55,var.urb17, by=c('chrBase'),all=F)
merge57<-merge(merge56,var.urb18, by=c('chrBase'),all=F)
merge58<-merge(merge57,var.urb19, by=c('chrBase'),all=F)
merge59<-merge(merge58,var.urb20, by=c('chrBase'),all=F)
merge60<-merge(merge59,var.rur1, by=c('chrBase'),all=F)
merge61<-merge(merge60,var.rur2, by=c('chrBase'),all=F)
merge62<-merge(merge61,var.rur3, by=c('chrBase'),all=F)
merge62<-merge(merge61,var.rur3, by=c('chrBase'),all=F)
merge63<-merge(merge62,var.rur4, by=c('chrBase'),all=F)
merge64<-merge(merge63,var.rur5, by=c('chrBase'),all=F)
merge65<-merge(merge64,var.rur6, by=c('chrBase'),all=F)
merge66<-merge(merge65,var.rur7, by=c('chrBase'),all=F)
merge67<-merge(merge66,var.rur8, by=c('chrBase'),all=F)
merge68<-merge(merge67,var.rur9, by=c('chrBase'),all=F)
merge69<-merge(merge68,var.rur10, by=c('chrBase'),all=F)
colnames(merge69)<-c("position","bar9.X","bar10.X","bar11.X","bar12.X","bar13.X"
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

write.csv(merge69,file="merge_meth_par_position_par_indiv.txt", row.names = F, quote=F)
fichierOK<-read.table("merge_meth_par_position_par_indiv.txt",sep=",", header=T)

fichiertest<-merge69[1:10,1:10]
datatest<-t(fichiertest)
colnames(datatest)<-datatest[1,]

data<-t(fichierOK)
colnames(data)<-data[1,]
data2<-data[2:nrow(data),]

write.table(data2,file="merge_meth_par_position_par_indiv_transposee.txt",sep=",", row.names = F, quote=F)


##### RDA ##############################################################
library(vegan)
## donnees methylation
dataOK<-read.table("merge_meth_par_position_par_indiv_transposee.txt",sep=",", header=T)

## explainatory variables 
data_detail=read.csv("indiv_detail_env.csv",header=TRUE, sep=",")
colnames(data_detail)<-c("indiv","latitude","milieu","sexe","code","milieu.sexe","ville")
data_detail$log.lat<-log(data_detail$latitude)
names(data_detail)
str(data_detail)

#### full RDA ##########################################################
rda1 <- rda(dataOK ~ ville+milieu+sexe,data_detail)
rda1
a1=RsquareAdj(rda1)
a2=summary(eigenvals(rda1, model = "constrained"))
screeplot(rda1)
a3 <- anova.cca(rda1, parallel=getOption("mc.cores"), permutations = how(nperm=99))
a4=anova(rda1, step=100)
a5=anova(rda1, by="margin", step=100)
a6=head(summary(rda1))
a1
a2
a3
a4
a5
a6
plot(rda1)

## Constrained for the habitat
cond.rda <- rda(dataOK ~ milieu + Condition(ville + sexe), data_detail)
cond.rda
b1=RsquareAdj(cond.rda)
b1
b2=summary(eigenvals(cond.rda, model = "constrained"))
b2
screeplot(cond.rda)
b3 <- anova.cca(cond.rda, parallel=getOption("mc.cores"), permutations = how(nperm=99))
b3
b4=anova(cond.rda, step=100)
b4
b5=anova(cond.rda, by="margin", step=100)
b5
b6=head(summary(cond.rda))
b6

### Constrained for the city
cond.rda <- rda(dataOK ~ ville + Condition(milieu + sexe), data_detail)
cond.rda
c1<-RsquareAdj(cond.rda)
c2<-summary(eigenvals(cond.rda, model = "constrained"))
screeplot(cond.rda)
c3<-anova.cca(cond.rda, parallel=getOption("mc.cores"), permutations = how(nperm=99))
c4<-anova(cond.rda, step=100)
c5<-anova(cond.rda, by="margin", step=100)
c6<-head(summary(cond.rda))

### Constrained for the sex
cond.rda <- rda(dataOK ~ sexe + Condition(milieu + ville), data_detail)
cond.rda
d1<-RsquareAdj(cond.rda)
d2<-summary(eigenvals(cond.rda, model = "constrained"))
screeplot(cond.rda)
d3<-anova.cca(cond.rda, parallel=getOption("mc.cores"), permutations = how(nperm=99))
d4<-anova(cond.rda, step=100)
d5<-anova(cond.rda, by="margin", step=100)
d6<-head(summary(cond.rda))


###### graphs #####################################
bg <- c("#E66101","#5E3C99","#FDB863") 
bg2=c(2,1,17,16) #shapes

pdf(file = "rda_1_2.pdf", width = 6, height = 5)
plot(rda1,scaling=3,  type="n",xlim=c(-40,40),ylim=c(-30,30),xlab="RDA1",ylab="RDA2")
#points(rda1, display="species", pch=1,  cex=0.25, col="gray32", scaling=3)#SNPs
points(rda1, cex=0.75, display="sites", pch=bg2[data_detail$milieu.sexe], scaling=3, col=bg[as.factor(data_detail$latitude)])#tits
text(rda1, scaling=3, display="bp", col="black",labels=c("Montpellier","Warsaw","Urban","Male"),  cex=1,xlab="RDA1",ylab="RDA2")#predictors
legend("topleft", legend=levels(data_detail$milieu.sexe), bty="n", pch=bg2, cex=1)
dev.off()

pdf(file = "rda_1_3.pdf", width = 6, height = 5)
plot(rda1,scaling=3,  type="n",xlim=c(-40,40),ylim=c(-30,30), choice=c(1,3),xlab="RDA1",ylab="RDA3")
#points(rda1, display="species", pch=1, cex=0.25, col="gray32", scaling=3, choice=c(1,3))#SNPs
points(rda1, cex=0.75, display="sites", pch=bg2[data_detail$milieu.sexe], scaling=3, col=bg[as.factor(data_detail$latitude)], choice=c(1,3))#tits
text(rda1, scaling=3, display="bp", col="black", cex=1,labels=c("Montpellier","Warsaw","Urban","Male"),choice=c(1,3),xlab="RDA1",ylab="RDA3")#predictors
legend("topright", legend=levels(data_detail$ville), bty="n", col=bg, pch=15, cex=1)
dev.off()

plot(rda1,scaling=3,  type="n",xlim=c(-40,40),ylim=c(-30,30), choice=c(1,3),xlab="RDA1",ylab="RDA3")
#points(rda1, display="species", pch=1, cex=0.25, col="gray32", scaling=3, choice=c(1,3))#SNPs
points(rda1, cex=0.75, display="sites", pch=bg2[data_detail$milieu.sexe], scaling=1, col=bg[as.factor(data_detail$latitude)], choice=c(1,3))#tits
text(rda1, scaling=3, display="bp", col="black", cex=1,labels=c("Montpellier","Warsaw","Urban","Male"),choice=c(1,3),xlab="RDA1",ylab="RDA2")#predictors
legend("topright", legend=levels(data_detail$ville), bty="n", col=bg, pch=15, cex=1)

par(mfrow=c(1,1))

plot(cond.rda,scaling=1,  type="n",ylim=c(-10,10))
#points(rda1, display="species", pch=1,  cex=0.25, col="gray32", scaling=3)#SNPs
points(cond.rda, cex=0.75, display="sites", pch=bg2[data_detail$milieu.sexe], scaling=3, col=bg[as.factor(data_detail$latitude)])#tits
text(cond.rda, scaling=3, display="bp", col="purple", cex=1)#predictors
legend("topleft", legend=levels(data_detail$milieu.sexe), bty="n", pch=bg2, cex=1)


## loadings RDA
load.rda.con <- scores(cond.rda, choices=1, display="species")
hist(load.rda.con[,1], main="Loadings on RDA1")
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}
cand1 <- outliers(load.rda.con[,1],3) # 1840
ncand <- length(cand1)
ncand
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
colnames(cand1) <- c("axis","snp","loading")
cand1$snp <- as.character(cand1$snp)
write.table(cand1, "loadings-RDAconstrained-outliers.txt")


