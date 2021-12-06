#!/usr/bin/Rscript
### PURPOSE OF THE SCRIPT: Performing a RDA on methylation data without Z

########################################################################
############    RDA on methylation data without       ##################
############			 Z chromosome      		      ##################
########################################################################

ls()
rm(list=ls())
ls()


library(vegan)
setwd("/PATH/TO/06_results/rda/without_Z")
fichierOK<-read.table("merge_meth_par_position_par_indiv_sans_Z.txt",sep=",", header=T)
data<-t(fichierOK)
colnames(data)<-data[1,]
data2<-data[2:nrow(data),]
write.table(data2,file="merge_meth_par_position_par_indiv_transposee_sans_Z.txt",sep=",", row.names = F, quote=F)
data_detail=read.table("indiv_detail_env.csv",header=TRUE, sep=",")
colnames(data_detail)<-c("indiv","latitude","milieu","sexe","code","milieu.sexe","ville")
str(data_detail)
dataOK<-read.table("merge_meth_par_position_par_indiv_transposee_sans_Z.txt",sep=",", header=T)

#### full RDA normale ##################################################
rda1 <- rda(dataOK ~ ville+milieu+sexe,data_detail)
rda1
a1<-RsquareAdj(rda1)
a2<-summary(eigenvals(rda1, model = "constrained"))
screeplot(rda1)
a3 <- anova.cca(rda1, parallel=getOption("mc.cores"), permutations = how(nperm=99))
a4<-anova(rda1, step=100)
a5<-anova(rda1, by="margin", step=100)
a6<-head(summary(rda1))
a1
a2
a3
a4
a5
a6
plot(rda1)


##RDA2
## Constrained for habitat
cond.rda <- rda(dataOK ~ milieu + Condition(ville + sexe), data_detail)
cond.rda
b1<-RsquareAdj(cond.rda)
b2<-summary(eigenvals(cond.rda, model = "constrained"))
screeplot(cond.rda)
b3<-anova.cca(cond.rda, parallel=getOption("mc.cores"), permutations = how(nperm=99))
b4<-anova(cond.rda, step=100)
b5<-anova(cond.rda, by="margin", step=100)
b6<-head(summary(cond.rda))

## Constrained for city
cond.rda <- rda(dataOK ~ ville + Condition(milieu + sexe), data_detail)
cond.rda
c1<-RsquareAdj(cond.rda)
c2<-summary(eigenvals(cond.rda, model = "constrained"))
screeplot(cond.rda)
c3<-anova.cca(cond.rda, parallel=getOption("mc.cores"), permutations = how(nperm=99))
c4<-anova(cond.rda, step=100)
c5<-anova(cond.rda, by="margin", step=100)
c6<-head(summary(cond.rda))

## Constrained for sex
cond.rda <- rda(dataOK ~ sexe + Condition(milieu + ville), data_detail)
cond.rda
d1<-RsquareAdj(cond.rda)
d2<-summary(eigenvals(cond.rda, model = "constrained"))
screeplot(cond.rda)
d3<-anova.cca(cond.rda, parallel=getOption("mc.cores"), permutations = how(nperm=99))
d4<-anova(cond.rda, step=100)
d5<-anova(cond.rda, by="margin", step=100)
d6<-head(summary(cond.rda))
