rm(list=ls())
require(ggplot2)
require(data.table)
require(dplyr)
require(tidyr)
require(ggsci)
require(lmerTest)
require(mgcv) 

# Define sourvce data location and read data 
SourceDataLoc <- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript Feline ERBB facilitation/Molecular Systems Biology/Revision and submission folder/Source Data/"
Intermediateloc <- paste0(SourceDataLoc,"Figure 5/")
FfreqDat <- data.table( read.csv( file=paste0(Intermediateloc,
                                        "SourceData_Figure_5C_Cancer growth increase with fibroblast freq DMSO and fulv.csv")))


collabs<-(ggsci::pal_aaas()(6)[c(2,3,1)])

ggplot(FfreqDat, aes(y= RGR-RGRalone, x= 100*(1-InitialCancerFraction)))+
  geom_smooth(method="gam",col="black",fill="black", formula=y~s(x,k=3,bs="cr"),alpha=0.2,se=T)+
  stat_boxplot(aes(fill=CompositionLab,col=CompositionLab, group= interaction(InitialCancerFraction,CancerCellName,ResistanceLab,Fulvestrant_nM,CompositionLab)),geom="errorbar",position=position_dodge(width=8))+
  geom_boxplot(col="black",aes(fill=CompositionLab,col=CompositionLab, group= interaction(InitialCancerFraction,CancerCellName,ResistanceLab,Fulvestrant_nM,CompositionLab)),
               width=6,position=position_dodge(width=8))+
  geom_point(aes(fill=CompositionLab,col=CompositionLab),pch=21,col="black",
             position=position_dodge(width=8))+
  theme_classic(base_size=26)+
  facet_grid(wrapLab~.)+theme(aspect.ratio=1)+
  labs(y="Cancer growth rate \n (compared to monoculture)",x="Fibroblast coculture frequency \n(% of fixed cancer population)")+
  scale_color_manual(name="",values=collabs[-2])+
  scale_fill_manual(name="",values=collabs[-2])+
  scale_x_continuous(breaks=(100*unique(1-FfreqDat$InitialCancerFraction))  )+
  theme(aspect.ratio=1,legend.position = "bottom")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
  theme(legend.position="none")

#ggsave(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Fibroblast_Experiments/coculture/20230223_Feng Afatinib Fulvestrant CAMA1 MCF7 T47D/Modeling_JG/MUTUALISM FINALIZE part3 cancer growth rate improves under more frequent fibroblast coculture DMSO and fulv treatment.pdf", width=9, height=12, dpi=320)
#ggsave("~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/MUTUALISM FINALIZE part3 cancer growth rate improves under more frequent fibroblast coculture DMSO and fulv treatment.pdf",height=9,width=12,dpi = 320)

# Stats
# Gam models predicting cancer growth rates as a function of initial fibroblast frequency under DMSO or fulv 
dmod<-gam( I(RGR-RGRalone)~ s(I(100*(1-InitialCancerFraction)),k=3),data=FfreqDat[wrapLab=="DMSO"])
fmod<-gam( I(RGR-RGRalone)~ s(I(100*(1-InitialCancerFraction)),k=3),data=FfreqDat[wrapLab!="DMSO"])
summary(dmod)
summary(fmod)

# linear mixed model to get effect sizes of coculture with fibroblasts
FfreqDat[,CancerCellNameSensitivity:= as.factor(paste0(CancerCellName,Sensitivity))]
dmodlmer<-lmer( I(RGR-RGRalone)~ I(100*(1-InitialCancerFraction))+
                  (-1+InitialCancerFraction|CancerCellNameSensitivity ),data=FfreqDat[wrapLab=="DMSO"])
fmodlmer<-lmer( I(RGR-RGRalone)~ I(100*(1-InitialCancerFraction))+
                  (-1+InitialCancerFraction|CancerCellNameSensitivity ),data=FfreqDat[wrapLab!="DMSO"])

summary(dmodlmer)
summary(fmodlmer)
###### 