rm(list=ls())
require(mgcv);require(data.table);require(dplyr);require(ggplot2);require(tidyr);require(lme4);require(lmerTest);require(parallel)
library(effects);require(umap)
require(Rfast);require(ider)
library("dendextend");library(ggdendro);require(ggsci);require(viridis)
require("Rdimtools")

discStatsData <- data.table(Cohort="Discovery cohort",read.csv(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/Figure4/DiscoveryCohort_Stats_FibroblastERBBligandexpressionincreasesmyCAFs.csv"))
validStatsData <- data.table(Cohort="Validation cohort",read.csv(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/Figure4/ValidationCohort_Stats_FibroblastERBBligandexpressionincreasesmyCAFs.csv"))

discLU<-data.table(emtlevel=sort( unique(discStatsData$emtlevel)) , EMTlevel= formatC(  1:length(unique(discStatsData$emtlevel))   , width = 2, format = "d", flag = "0")  )
discLU[,PlotID:=paste0("D",EMTlevel)]
validLU<-data.table(emtlevel=sort( unique(validStatsData$emtlevel))[-1] , EMTlevel=formatC( 1:length(unique(validStatsData$emtlevel)[-1] )   , width = 2, format = "d", flag = "0"))
validLU[,PlotID:=paste0("D",EMTlevel)]


plotall <- rbind( merge(discLU,discStatsData, by="emtlevel"), merge(validLU,validStatsData, by="emtlevel"))


ggplot(plotall[],aes(y= lnvalmu ,x=(as.numeric(EMTlevel)),group=interaction(Cohort,EMTlevel)    ,fill= as.numeric(EMTlevel) ))+
  theme_classic(base_size=26)+theme(aspect.ratio=1)+ # or aspect.ratio=0.71
  labs(x="Fibroblast myCAF differentiation level",y="ERBB ligand signaling (scaled)")+
  geom_boxplot(position = position_dodge(.7),outlier.color=NA) +
  stat_boxplot(aes(col=Cohort),geom = "errorbar", position = position_dodge(.7),width =0.5,size=2) +
  geom_point(data=plotall[abs(lnvalmu)<2],position=position_dodge(.7),size=0.4) +
  scale_fill_gradient(guide=F,low =pal_aaas()(6)[3] ,high =pal_aaas()(6)[1]) + 
  scale_color_manual(values=c("slategray","grey")) + 
  theme(legend.position = "none")+
  lims(y=c(-0.5,2))+
  scale_x_continuous(breaks = 1:10)
#ggsave("~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Discovery and Validation combined Fibroblast ERBBligandexpression increases in myCAFs.png",height=8,width=8,dpi = 320)


load("/Users/jason/Jason Griffiths Dropbox/CancerMacrophageCPM_FELINE1/scRNAseq_FelineCohort1and2_CellMetaData.RData")

plotallSummar <- data.table(merge(plotall[], CellMetaData%>%dplyr::select(Cell.ID,Patient.Study.ID,Day),by="Cell.ID")%>%group_by(emtlevel,EMTlevel,PlotID,Cohort,Timepoint,Day,
                                                                                                                               Patient.Study.ID)%>%
                            dplyr::summarise(lnvalmu1=mean(lnvalmu),lnvalmu2=median(lnvalmu)))



ggplot(plotallSummar,aes(y= lnvalmu1 ,x=(as.numeric(EMTlevel)) ))+
  theme_classic(base_size=26)+theme(aspect.ratio=1)+ 
  labs(x="Fibroblast myCAF differentiation level",y="ERBB ligand signaling (scaled)")+
  geom_boxplot(aes(group=interaction(Cohort,EMTlevel)    ,fill= as.numeric(EMTlevel) ),position = position_dodge(.7),outlier.color=NA) +
  stat_boxplot(aes(group=interaction(Cohort,EMTlevel)    ,fill= as.numeric(EMTlevel),col=Cohort),geom = "errorbar", position = position_dodge(.7),width =0.75,size=2) +
  geom_point(data=plotallSummar[],aes(group=interaction(Cohort,EMTlevel)  ,col=Cohort  ,fill= as.numeric(EMTlevel) ),pch=21,position=position_dodge(.7),size=2.5) +
  scale_fill_gradient(guide=F,low =pal_aaas()(6)[3] ,high =pal_aaas()(6)[1]) + 
  scale_color_manual(values=c("slategray","grey")) + 
  theme(legend.position = "none")+
  scale_x_continuous(breaks = 1:10)

#ggsave("~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Discovery and Validation combined Fibroblast ERBBligandexpression increases in myCAFs sampleAverage.pdf",height=7.5,width=7.5,dpi = 320)

#write.csv(plotallSummar[Cohort=="Discovery cohort"]%>%dplyr::select(-c(lnvalmu2)),file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/Figure4/DiscoveryCohort_Stats_FibroblastERBBligandexpressionincreasesmyCAFssampleAverage.csv")
#write.csv(plotallSummar[Cohort=="Validation cohort"]%>%dplyr::select(-c(lnvalmu2)),file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/Figure4/ValidationCohort_Stats_FibroblastERBBligandexpressionincreasesmyCAFssampleAverage.csv")

summary(lm(lnvalmu1~emtlevel,plotallSummar[Cohort=="Discovery cohort"]))
summary(lm(lnvalmu1~emtlevel,plotallSummar[Cohort=="Validation cohort"]))

plotallSummar[,EMTlevel2:=(as.numeric(EMTlevel)>7)]
ggplot(plotallSummar,aes(y= lnvalmu1 ,x=(as.numeric(EMTlevel2)) ))+
  theme_classic(base_size=26)+theme(aspect.ratio=1)+ 
  labs(x="Fibroblast myCAF differentiation level",y="ERBB ligand signaling (scaled)")+
  geom_boxplot(aes(group=interaction(Cohort,EMTlevel2)    ,fill= as.numeric(EMTlevel2) ),position = position_dodge(.7),outlier.color=NA) +
  stat_boxplot(aes(group=interaction(Cohort,EMTlevel2)    ,fill= as.numeric(EMTlevel2),col=Cohort),geom = "errorbar", position = position_dodge(.7),width =0.75,size=2) +
  geom_point(data=plotallSummar[],aes(group=interaction(Cohort,EMTlevel2)  ,col=Cohort  ,fill= as.numeric(EMTlevel2) ),pch=21,position=position_dodge(.7),size=2.5) +
  scale_fill_gradient(guide=F,low =pal_aaas()(6)[3] ,high =pal_aaas()(6)[1]) + 
  scale_color_manual(values=c("slategray","grey")) + 
  theme(legend.position = "none")+
  scale_x_continuous(breaks = 1:10)

ggplot(plotallSummar,aes(y= lnvalmu1 , x=Cohort ))+
  theme_classic(base_size=26)+theme(aspect.ratio=1)+ # or aspect.ratio=0.71
  labs(x="Cohort",y="ERBB ligand signaling (scaled)")+
  geom_boxplot(aes(group=interaction(Cohort,EMTlevel2)    ,fill= as.numeric(EMTlevel2) ),position = position_dodge(.7),width=0.5,outlier.color=NA) +
  stat_boxplot(aes(group=interaction(Cohort,EMTlevel2)    ,fill= as.numeric(EMTlevel2)),geom = "errorbar", position = position_dodge(.7),width =0.75) +
  geom_point(data=plotallSummar[],aes(group=interaction(Cohort,EMTlevel2)  ,fill= as.numeric(EMTlevel2) ), pch=21,position=position_dodge(.7),size=2.5) +
  scale_fill_gradient(guide=F,low =pal_aaas()(6)[3] ,high =pal_aaas()(6)[1]) + 
  scale_x_discrete(labels=c("Discovery","Validation"))
#Fibroblast myCAF differentiation
#ggsave("~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Discovery and Validation Fibroblast differentiation increases ERBBligandexpression of fibrobasts.pdf",height=8,width=8,dpi = 320)


ggplot(plotall[Cohort=="Validation cohort"],aes(y= lnvalmu ,x=emtlevel,group=interaction(Cohort,emtlevel)          #,fill=emtlevel)
   ))+
  theme_classic(base_size=26)+theme(aspect.ratio=0.71)+
  labs(x="Fibroblast myCAF phenotype \n (TGFb stimulated)",y="ERBB ligand signaling \n (scaled)")+
  geom_boxplot() +stat_boxplot(geom = "errorbar",
                               width = 0.5) +geom_point(size=0.4) +
  scale_fill_gradient(low =pal_aaas()(6)[3] ,high =pal_aaas()(6)[1])+theme(legend.position = "none")+
  facet_wrap(~Cohort,ncol=1, scales="free")
#ggsave("~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Validation Fibroblast ERBBligandexpression increases in myCAFs.pdf",height=8,width=8,dpi = 320)
