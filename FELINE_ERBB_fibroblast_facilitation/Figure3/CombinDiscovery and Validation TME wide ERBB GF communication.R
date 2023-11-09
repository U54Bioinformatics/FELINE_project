rm(list=ls())
require(mgcv);require(data.table);require(dplyr);require(ggplot2);require(tidyr);require(lme4);require(lmerTest);require(parallel)
library(effects);require(umap)
require(Rfast);require(ider)
library("dendextend");library(ggdendro);require(ggsci);require(viridis)
require("Rdimtools")

discData <- data.table(read.csv(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/Figure3/DiscoveryCohortTMEwideERBBcommunicationwithtoCancerPrePost.csv"))
validData <- data.table(read.csv(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/Figure3/ValidationCohortTMEwideERBBcommunicationwithtoCancerPrePost.csv"))
discData[,grpx:=paste(Patient.Study.ID,TimePoint,sep="AT")]
validData[,grpx:=paste(Patient.Study.ID,TimePoint,sep="AT")]

d_wide<-data.table(Cohort="Discovery cohort", discData %>% dplyr::select(-c(Patient.Study.ID,TimePoint,ARM,X,treat,CommunicationScore)) %>% spread(grpx,ScaledlnCommunicationScore,fill=0) %>%gather(grpx,y,-CommunicationReceptor))
v_wide<-data.table(Cohort="Validation cohort", validData %>% dplyr::select(-c(Patient.Study.ID,TimePoint,ARM,X,treat,CommunicationScore)) %>% spread(grpx,ScaledlnCommunicationScore,fill=0) %>%gather(grpx,y,-CommunicationReceptor))

plotall <- rbind(d_wide,v_wide)
plotall[,c("Patient.Study.ID","TimePoint"):= tstrsplit(grpx,"AT",fixed=T)]
plotall$TimePoint <- factor(plotall$TimePoint , levels= rev( unique(plotall$TimePoint) ) )

ggplot(plotall, 
       aes(x=CommunicationReceptor ,y= y , fill=TimePoint,group=interaction(CommunicationReceptor,TimePoint) ))+
  geom_boxplot(outlier.color=NA)+
  geom_point(position=position_dodge(width=0.9))+
  theme_classic(base_size=26)+
  scale_fill_manual(name="Timepoint",values= c("yellow","red")) +
  labs(x="Cancer growth factor receptor", y="Growth factor signaling to cancer cells \n from across the TME (relative to baseline average)")+ theme(aspect.ratio=1,legend.position="top")+
  facet_wrap(~Cohort,ncol=1)+stat_boxplot(geom = "errorbar" )+
  theme( strip.text = element_text(face="bold" ))
#ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Discovery and validation cohort2 Growth factor signalling to cancer cells increased relative to baseline ERBB onlyB.pdf", height=12.25, width = 10, dpi=320)

ggplot(plotall, 
       aes(x=CommunicationReceptor ,y= y , fill=TimePoint,group=interaction(CommunicationReceptor,TimePoint) ))+
  geom_boxplot(outlier.color=NA)+
  geom_point(position=position_dodge(width=0.9))+
  theme_classic(base_size=26)+
  scale_fill_manual(name="Timepoint",values= c("yellow","red")) +
  labs(x="Cancer growth factor receptor", y="Tumor-wide growth factor signaling to cancer cells \n (relative to baseline average)")+ theme(aspect.ratio=1,legend.position="top")+
  facet_wrap(~Cohort,ncol=1)+stat_boxplot(geom = "errorbar" )+
  theme( strip.text = element_text(face="bold" ))
#ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Discovery and validation cohort2 Growth factor signalling to cancer cells increased relative to baseline ERBB onlyC.pdf", height=12.25, width = 10, dpi=320)


