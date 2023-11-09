rm(list=ls())
require(mgcv);require(data.table);require(dplyr);require(ggplot2);require(tidyr);require(lme4);require(lmerTest);require(parallel)
library(effects);require(umap)
require(Rfast);require(ider)
library("dendextend");library(ggdendro);require(ggsci);require(viridis)
require("Rdimtools")


discData <- data.table(Cohort="Discovery cohort", read.csv(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/Figure3/DiscoveryCohortFibroblastERBBcommunicationwithtoCancerPrePost.csv"))
validData <- data.table(Cohort="Validation cohort", read.csv(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/Figure3/ValidationCohortFibroblastERBBcommunicationwithtoCancerPrePost.csv"))
plotall <- rbind(discData,validData)
plotall$TimePoint <- factor(plotall$TimePoint , levels= rev( unique(plotall$TimePoint) ) )


ggplot(plotall, 
       aes(x=CommunicationReceptor ,y= logFibroblasttoCancerCommunicationScore , fill=TimePoint,group=interaction(CommunicationReceptor,TimePoint) ))+
  geom_boxplot(outlier.color=NA)+
  geom_point(position=position_dodge(width=0.9))+
  theme_classic(base_size=26)+
  scale_fill_manual(name="Timepoint",values=c("yellow","red")) +
  labs(x="Cancer growth factor receptor", y="Fibroblast-cancer cell growth factor signaling \n (relative to baseline average)")+ theme(aspect.ratio=1,legend.position="top")+
  facet_wrap(~Cohort,ncol=1)+stat_boxplot(geom = "errorbar" )+
  theme( strip.text = element_text(face="bold" ))

#ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Discovery and validation cohort2 FibroblastERBBcommunicationwithtoCancerPrePostB.pdf", height=12.25, width = 10, dpi=320)


