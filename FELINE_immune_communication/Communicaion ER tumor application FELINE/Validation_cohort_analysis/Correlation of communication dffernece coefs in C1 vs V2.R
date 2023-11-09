rm(list=ls())   
require(rTensor); require(nnTensor); library(abind)
require(data.table); require(dplyr); require(ggplot2); require(tidyr)
require(igraph)
require(ggsci)
load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Cohort2 When cells communicate_communication changes AllArms/Cohort2 trendsby response AllArms.RData")
C2_assessment <- assessment
load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/When cells communicate_communication changes AllArms/trendsby response AllArms.RData")
C1_assessment <- assessment

joinAssess <- merge(
  C1_assessment[Treat=="CombinationRibo"]%>%dplyr::select(Model, Pair.Name,LigandPhenoCelltype, ReceptorPhenoCelltype, Treat,coefficient,Estimate,pval),
  C2_assessment[Treat=="CombinationRibo"]%>%dplyr::select(Model, Pair.Name,LigandPhenoCelltype, ReceptorPhenoCelltype, Treat,coefficient,Estimate,pval) ,
  by=c("Model", "Pair.Name","LigandPhenoCelltype", "ReceptorPhenoCelltype", "Treat","coefficient"))

setnames( joinAssess, old=c("Estimate.x","Estimate.y","pval.x" , "pval.y"), new=c("EstimateC1","EstimateC2","pval.C1" , "pval.C2"))

summary( lm(EstimateC1~ EstimateC2, joinAssess[Model=="two way anova"][LigandPhenoCelltype!="B cells"][ReceptorPhenoCelltype!="B cells"] ))

joinAssess[ReceptorPhenoCelltype=="Macrophages",ReceptorPhenoCelltype:="Myeloid cells"]
joinAssess[ReceptorPhenoCelltype=="Normal epithelial cells",ReceptorPhenoCelltype:="Diploid epithelial cells"]
ggplot(joinAssess[][Model=="two way anova"][LigandPhenoCelltype!="B cells"][ReceptorPhenoCelltype!="B cells"] , aes(EstimateC1, EstimateC2, col=ReceptorPhenoCelltype,fill=ReceptorPhenoCelltype)  ) +  
  theme_classic(base_size = 22) + 
  geom_point(size=0.15,alpha=0.8) + geom_smooth(method="lm") + 
  theme(aspect.ratio= 1) +
  scale_x_continuous(name="Discovery cohort \n communication difference between \n resistant and sensitive tumors",breaks=5*(-3:3),limits=c(-7,12)) + 
  scale_y_continuous(name="Validation cohort \n communication difference between \n resistant and sensitive tumors",breaks=5*(-3:3),limits=c(-7,12)) +
  scale_colour_npg(name="Communication \n receiver \n cell type")+
  scale_fill_npg(name="Communication \n receiver \n cell type")
#ggsave(filename = "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cohort 2 results/Cohort 1 vs 2 communication difference between resistant and sensitive tumors.png")
paperfile<- "/Users/jason/Dropbox/FELINE Project (1)/Manuscript  Feline immune communication/Figures Communication Project/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Validation IL communicaiton Ribo and letrozole Myeloid to T cell post treatment activation communications Validation C2.png"),height=10,width=10)
ggsave(paste0(paperfile,"Ribo and Letrozole SI Validation communicaiton consistency across cohorts.png"),height=10,width=10)


ggplot(joinAssess[Model=="two way anova"][LigandPhenoCelltype!="B cells"][ReceptorPhenoCelltype!="B cells"] , aes(EstimateC1, EstimateC2, col=ReceptorPhenoCelltype,fill=ReceptorPhenoCelltype)  ) +  
  theme_classic() + 
  geom_point(size=0.25,alpha=0.8) + geom_smooth(size=1.5,method="lm") + 
  theme(aspect.ratio= 1) +
  scale_x_continuous(name="Discovery cohort \n communication difference between \n resistant and sensitive tumors",breaks=5*(-3:3),limits=c(-7,12)) + 
  scale_y_continuous(name="Validation cohort \n communication difference between \n resistant and sensitive tumors",breaks=5*(-3:3),limits=c(-7,12)) +
  scale_colour_npg(name="Communication receiver \n cell type")+
  scale_fill_npg(name="Communication receiver \n cell type") +
  theme(axis.text=element_blank(), legend.position = "none", axis.title = element_blank())
ggsave(filename = "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cohort 2 results/BLANK Cohort 1 vs 2 communication difference between resistant and sensitive tumors.png")


ggplot(joinAssess[Model=="two way anova"][LigandPhenoCelltype!="B cells"][ReceptorPhenoCelltype!="B cells"] , aes(EstimateC1, EstimateC2, col=LigandPhenoCelltype,fill=LigandPhenoCelltype)  ) +  theme_classic() + geom_point(size=0.5,alpha=0.8) + geom_smooth(method="lm") + 
  theme(aspect.ratio= 1) +
  scale_x_continuous(name="Discovery cohort \n communication difference between \n resistant and sensitive tumors",breaks=5*(-3:3),limits=c(-7,12)) + 
  scale_y_continuous(name="Validation cohort \n communication difference between \n resistant and sensitive tumors",breaks=5*(-3:3),limits=c(-7,12)) +
  scale_colour_npg()+
  scale_fill_npg()+facet_wrap(~LigandPhenoCelltype)


