########################################################################
#############             PLOT MEAN METH              ##################
########################################################################

### libraries ##########################################################
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(dplyr)
theme_set(theme_pubr())
########################################################################


###############################################
#############   DATA WITH Z   #################
###############################################
setwd("/PATH/TO/06_results/mean_meth_plot")
data<-read.table("mean_meth_per_indiv_avecZ")
colnames(data)<-c("IND","meth")
data_order <- data[order(data$IND),] #order 
env <- read.csv("indiv_detail_env.csv", header=T) #load env
str(env) # Look at the structure of env
env$individual <- as.character(env$INDV) # Make individual names characters (not factors)
env2 <- env[order(env$INDV),] #order env
identical(data_order$IND, env2$INDV)
data_complet<-cbind(data_order,env2)
bar<-subset(data_complet,data_complet$ville=="Barcelona")
mtp<-subset(data_complet,data_complet$ville=="Montpellier")
var<-subset(data_complet,data_complet$ville=="Warsaw")


## PLOT GLOBAL
plot <-ggplot(data_complet, aes(x = meth)) + 
  geom_density(aes(y=..count..,fill = ville, color=ville, linetype=ENV), alpha = 0.4)+
  scale_fill_manual(values=c("#e66101","#e66101","#5e3c99","#5e3c99","#fdb863","#fdb863"))+
  scale_color_manual(values=c("#e66101","#e66101","#5e3c99","#5e3c99","#fdb863","#fdb863"))
plot

## PLOT CITY ######
mean_bar <- bar %>% 
  group_by(ville) %>%
  summarise(grp.mean = mean(meth))
mean_bar

plot_bar <-ggplot(bar, aes(x = meth)) + 
  geom_density(aes(y=..count..,fill=ville,color=ville, linetype=ENV), alpha = 0.4)+
  scale_fill_manual(values=c("#e66101","#e66101"))+
  scale_color_manual(values=c("#e66101","#e66101"))+
  xlim(71,80)+
  ylim(0,11)+
  labs(x="Mean meth per individual (%)",y=("Number of individuals"))+
  geom_vline(aes(xintercept = grp.mean, color=ville,linetype=ville),data=mean_bar, size = 0.6)+
  scale_linetype_manual(values=c("solid","dashed","solid","dashed"))
plot_bar


mean_mtp <- mtp %>% 
  group_by(ville) %>%
  summarise(grp.mean = mean(meth))
mean_mtp

plot_mtp <-ggplot(mtp, aes(x = meth)) + 
  geom_density(aes(y=..count..,fill=ville,color=ville, linetype=ENV), alpha = 0.4)+
  scale_fill_manual(values=c("#5e3c99","#5e3c99"))+
  scale_color_manual(values=c("#5e3c99","#5e3c99"))+
  xlim(71,80)+
  ylim(0,11)+
  labs(x="Mean meth per individual (%)",y=("Number of individuals"))+
  geom_vline(aes(xintercept = grp.mean, color=ville,linetype=ville),data=mean_mtp, size = 0.6)+
  scale_linetype_manual(values=c("solid","dashed","solid","dashed"))
plot_mtp

mean_var <- var %>% 
  group_by(ville) %>%
  summarise(grp.mean = mean(meth))
mean_var

plot_var <-ggplot(var, aes(x = meth)) + 
  geom_density(aes(y=..count..,fill=ville,color=ville, linetype=ENV), alpha = 0.4)+
  scale_fill_manual(values=c("#fdb863","#fdb863"))+
  scale_color_manual(values=c("#fdb863","#fdb863"))+
  xlim(71,80)+
  ylim(0,11)+
  labs(x="Mean meth per individual (%)",y=("Number of individuals"))+
  geom_vline(aes(xintercept = grp.mean, color=ville,linetype=ville),data=mean_var, size = 0.6)+
  scale_linetype_manual(values=c("solid","dashed","solid","dashed"))
plot_var

grid.arrange(plot_bar, plot_mtp, plot_var, ncol=1, nrow =3)

#########################

### PLOT SEX #######
mean_sex <- data_complet %>% 
  group_by(sexe) %>%
  summarise(grp.mean = mean(meth))
mean_sex

plot_sex <-ggplot(data_complet, aes(x = meth)) + 
  geom_density(aes(y=..count..,fill=sexe,color=sexe), alpha = 0.4)+
  scale_fill_manual(values=c("#f1b6da","#2c7bb6"))+
  scale_color_manual(values=c("#f1b6da","#2c7bb6"))+
  xlim(71.5,78.5)+
  ylim(0,20)+
  labs(x="Mean meth per individual (%)",y=("Number of individuals"))+
  geom_vline(aes(xintercept = grp.mean, color=sexe),data=mean_sex, size = 0.6)
plot_sex

#########################

### PLOT CITY ######
mean_ville <- data_complet %>% 
  group_by(ville) %>%
  summarise(grp.mean = mean(meth))
mean_ville

plot_ville <-ggplot(data_complet, aes(x = meth)) + 
  geom_density(aes(y=..count..,fill=ville,color=ville), alpha = 0.4)+
  scale_fill_manual(values=c("#e66101","#5e3c99","#fdb863"))+
  scale_color_manual(values=c("#e66101","#5e3c99","#fdb863"))+
  xlim(71.5,78.5)+
  ylim(0,11)+
  labs(x="Mean meth per individual (%)",y=("Number of individuals"))+
  geom_vline(aes(xintercept = grp.mean, color=ville),data=mean_ville, size = 0.6)
plot_ville


#########################


### PLOT HABITAT #######
mean_env <- data_complet %>% 
  group_by(ENV) %>%
  summarise(grp.mean = mean(meth))
mean_env

plot_env <-ggplot(data_complet, aes(x = meth)) + 
  geom_density(aes(y=..count..,fill=ENV,color=ENV), alpha = 0.4)+
  scale_fill_manual(values=c("#4dac26","grey69"))+
  scale_color_manual(values=c("#4dac26","grey69"))+
  xlim(71.5,78.5)+
  ylim(0,20)+
  labs(x="Mean meth per individual (%)",y=("Number of individuals"))+
  geom_vline(aes(xintercept = grp.mean, color=ENV),data=mean_env, size = 0.6)
plot_env

#########################


grid.arrange(plot_env, plot_sex, plot_ville, ncol=1, nrow =3)

###########################################################
###########################################################

