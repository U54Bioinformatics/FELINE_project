rm(list=ls())
require(mgcv);require(data.table);require(dplyr);require(ggplot2);require(tidyr);require(lme4);require(lmerTest);require(parallel)
discData <- data.table(Cohort="Discovery cohort", read.csv(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/Figure4/Fig4C/DiscoveryCohortFibroblastdifferentiationsignalsfromcancerandnoncancer.csv"))
validData <- data.table(Cohort="Validation cohort", read.csv(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/Figure4/Fig4C/ValidationCohortFibroblastdifferentiationsignalsfromcancerandnoncancer.csv"))

plotall <- rbind( discData, validData)

ggplot(plotall,  
       aes(y=scalelnTransductionMu,x=HPMR.Ligand ,
           fill=Malignance ) ) +
  theme_classic(base_size=26)+
  theme(aspect.ratio=1)+
  geom_boxplot(aes(group=interaction(HPMR.Ligand,Malignance)))+
  geom_jitter(pch=21,col="black",size=2.5, position = position_dodge(width=0.5)) +
  scale_fill_manual(name="Cell type",values=c(pal_aaas("default")(6)[2],"grey"))+
  labs(y="myCAF differentiation stimulus \n received by fibroblasts"
       ,x="Signaling ligand") +
  facet_wrap(~Cohort, ncol=2,scale="free") +
  theme(legend.position = "none")
#ggsave("~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Discovery and Validation combined Fibroblast activating signals by cell type.pdf",height=9,width=11,dpi = 320)
