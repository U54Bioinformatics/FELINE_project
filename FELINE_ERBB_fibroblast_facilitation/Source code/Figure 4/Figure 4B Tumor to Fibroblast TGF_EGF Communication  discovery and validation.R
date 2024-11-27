rm(list=ls())
library(glmnet);require(dplyr);require(data.table);require(tidyr)
require(ggraph);require(tidygraph)
require(scales);require(ggsci); require(lmerTest)
# Define data location
SourceDataLoc <- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript Feline ERBB facilitation/Molecular Systems Biology/Revision and submission folder/Source Data/"
Intermediateloc <- paste0(SourceDataLoc,"Figure 4/")

discovery_LRcomms <- data.table( read.csv(file= paste0(Intermediateloc, "SourceData_Figure_4AB_DiscoveryTumorComms.csv" ) ) ) 

# Treatment labels
discovery_LRcomms[,Treat:="CombinationRibo"]
discovery_LRcomms[ARM=="A",Treat:="LetrozoleAlone"]
# Scaling and transformation of signals
discovery_LRcomms[,scaleTransduction:=scale(TransductionMu,center=F),by=c("Pair.Name")] #"ReceptorPhenoCelltype",
discovery_LRcomms[!is.finite(TransductionMu),scaleTransduction:=0]
discovery_LRcomms[,scalelnTransductionMu:=scale(log(1+TransductionMu)),by=Pair.Name]
discovery_LRcomms[, c("HPMR.Ligand", "HPMR.Receptor") := tstrsplit(Pair.Name, "_", fixed=TRUE)]


# Specify relevant LR pairs
pathwaysLetrozole<-c( "CDH1_EGFR",
                      "TGFB2_TGFBR3",
                      "TGFB1_ITGAV",
                      "TGFB1_TGFBR1",
                      "TGFB2_TGFBR1",
                      "INHBA_TGFBR3",
                      "HBEGF_EGFR",
                      "GNAI2_EGFR",
                      "HBEGF_ERBB2",
                      "HBEGF_CD9",
                      "HBEGF_CD44" )

#  Gather LR pairs using typical ligands for the focal receptor (TGF for TGFB and EGF for EGFR)
part1<- discovery_LRcomms[(grepl("TGFB1_",Pair.Name)|grepl("TGFB2_",Pair.Name)|grepl("TGFB3_",Pair.Name))
  &grepl(c("TGFBR"),Pair.Name)][ReceptorPhenoCelltype=="Fibroblasts"]
part1[,Signal:="TGF"]
part2<- discovery_LRcomms[grepl("HBEGF_",Pair.Name)][ReceptorPhenoCelltype=="Fibroblasts"]
part2[,Signal:="HBEGF"]
# Join data
partboth <- rbind(part1,part2)[!LigandPhenoCelltype%in%c("B cells")]
partboth[LigandPhenoCelltype%in%c("CD4+ T cells","CD8+ T cells" ),LigandPhenoCelltype:="T cells"]
partboth[LigandPhenoCelltype%in%c("Pericytes","Endothelial cells" ),LigandPhenoCelltype:="Endothelial/Pericytes cells"]

ggplot(partboth[  ],#[  LigandPhenoCelltype%in%c("Cancer cells","Endothelial cells","Fibroblasts","Macrophages","Normal epothelial cells","T cells")],  
       aes(y=scalelnTransductionMu,x=(log(1+Day)) ,col=LigandPhenoCelltype ,fill=LigandPhenoCelltype) ) + theme_classic(base_size=26)+
  facet_grid(Signal~dynamic_class3)+ 
  geom_smooth(aes(linetype=LigandPhenoCelltype!="Cancer cells"),method="lm",se=T,alpha=0.3)+ geom_point(pch=21,col="black",size=2.5) +
  scale_color_viridis_d(name="Signaling \n cell type", option="B",direction = 1)+
  scale_fill_viridis_d( name="Signaling \n cell type", option="B",direction = 1)+
  scale_linetype_discrete( name="Malignance", labels=c("Cancer","Non-cancer"))+theme(aspect.ratio=1)+
  labs(y="Fibroblasts activation communications",x="Day")+scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))
#ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Fibroblast activating signals by cell type and tumor response over time.pdf",width=14, height=14)

# Group cancer and non-cancer cell types into two classes for presentation
partboth[,Malignance:="Cancer"]
partboth[LigandPhenoCelltype!="Cancer cells",Malignance:="Non-cancer"]
partboth[,TumorResponse:="Growing tumor"]
partboth[dynamic_class3!="Non-response",TumorResponse:="Shrinking tumor"]

ggplot(partboth[  Day==180]%>%group_by(Patient.Study.ID,TumorResponse,HPMR.Ligand,Signal,Malignance,dynamic_class3)%>%summarise(scalelnTransductionMu=mean(scalelnTransductionMu)),#[  LigandPhenoCelltype%in%c("Cancer cells","Endothelial cells","Fibroblasts","Macrophages","Normal epothelial cells","T cells")],  
       aes(y=scalelnTransductionMu,x=HPMR.Ligand ,#col=Malignance ,
           #shape=dynamic_class3,
           fill=Malignance ) ) +
  theme_classic(base_size=26)+
  #facet_grid(~Signal)+ 
  theme(aspect.ratio=1)+#legend.position="none")+
  #geom_smooth(aes(linetype=LigandPhenoCelltype!="Cancer cells"),method="lm",se=T,alpha=0.3)+
  geom_boxplot(aes(group=interaction(HPMR.Ligand,Malignance)))+
  geom_jitter(pch=21,col="black",size=2.5, position = position_dodge(width=0.5)) +
  scale_fill_manual(name="Cell type",values=c(pal_aaas("default")(6)[2],"grey"))+
  labs(y="Fibroblasts differentiation \n communication strength",x="Signaling ligand")

#  scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))
#ggsave(file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Fibroblast differentiation signals by cell type and tumor response over time.pdf",width=8.5, height=8.5)
#ggsave(file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Fibroblast differentiation signals from cancer and noncancerRed.pdf",width=10, height=10)

discData<-data.table( partboth[  Day==180]%>%group_by(Patient.Study.ID,TumorResponse,HPMR.Ligand,Signal,Malignance,dynamic_class3)%>%dplyr::summarise(scalelnTransductionMu=mean(scalelnTransductionMu)) )
#write.csv(discData,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/DiscoveryCohortFibroblastdifferentiationsignalsfromcancerandnoncancer.csv")

# Linear model analysis
discData$Malignance <- factor(discData$Malignance , levels=c("Non-cancer","Cancer"))
summary( lm(scalelnTransductionMu~ -1+HPMR.Ligand:Malignance,data=discData) )








## Replicate analysis for validation cohort
rm(list=ls())
SourceDataLoc <- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript Feline ERBB facilitation/Molecular Systems Biology/Revision and submission folder/Source Data/"
Intermediateloc <- paste0(SourceDataLoc,"Figure 4/")

validation_LRcomms <- data.table( read.csv(file= paste0(Intermediateloc, "SourceData_Figure_4AB_ValidationTumorComms.csv" ) ) ) 
validation_LRcomms[,Treat:="CombinationRibo"]
validation_LRcomms[ARM=="A",Treat:="LetrozoleAlone"]
validation_LRcomms[,scaleTransduction:=scale(TransductionMu,center=F),by=c("Pair.Name")] #"ReceptorPhenoCelltype",
validation_LRcomms[!is.finite(TransductionMu),scaleTransduction:=0]
validation_LRcomms[,scalelnTransductionMu:=scale(log(1+TransductionMu)),by=Pair.Name]
validation_LRcomms[,scaleTransduction2:=exp(scale(log(TransductionMu),center=F)),by=c("Pair.Name")] 
validation_LRcomms[, c("HPMR.Ligand", "HPMR.Receptor") := tstrsplit(Pair.Name, "_", fixed=TRUE)]

pathwaysLetrozole<-c( #"CDH1_EGFR",
  "TGFB2_TGFBR3",
  "TGFB1_ITGAV",
  "TGFB1_TGFBR1",
  "TGFB2_TGFBR1"
  #,
  #"INHBA_TGFBR3",
  #"HBEGF_EGFR",
  #"GNAI2_EGFR",
  #"HBEGF_ERBB2",
  #"HBEGF_CD9",
  #"HBEGF_CD44" 
)

part1<- validation_LRcomms[
][(grepl("TGFB1_",Pair.Name)|grepl("TGFB2_",Pair.Name)|grepl("TGFB3_",Pair.Name))
  &grepl(c("TGFBR"),Pair.Name)][ReceptorPhenoCelltype=="Fibroblasts"]
part1[,Signal:="TGF"]
part2<- validation_LRcomms[
][grepl("HBEGF_",Pair.Name)][ReceptorPhenoCelltype=="Fibroblasts"]
part2[,Signal:="HBEGF"]
partboth <- rbind(part1,part2)[!LigandPhenoCelltype%in%c("B cells")]
partboth[LigandPhenoCelltype%in%c("CD4+ T cells","CD8+ T cells" ),LigandPhenoCelltype:="T cells"]
partboth[LigandPhenoCelltype%in%c("Pericytes","Endothelial cells" ),LigandPhenoCelltype:="Endothelial/Pericytes cells"]
#partboth[LigandPhenoCelltype%in%c("Adipocytes","Macrophages" ),LigandPhenoCelltype:="Macrophages/Adipocytes"]
ggplot(partboth[  ],#[  LigandPhenoCelltype%in%c("Cancer cells","Endothelial cells","Fibroblasts","Macrophages","Normal epothelial cells","T cells")],  
       aes(y=scalelnTransductionMu,x=(log(1+Day)) ,col=LigandPhenoCelltype ,fill=LigandPhenoCelltype) ) + theme_classic(base_size=26)+
  facet_grid(Signal~dynamic_class3)+ 
  geom_smooth(aes(linetype=LigandPhenoCelltype!="Cancer cells"),method="lm",se=T,alpha=0.3)+ geom_point(pch=21,col="black",size=2.5) +
  scale_color_viridis_d(name="Signaling \n cell type", option="B",direction = 1)+
  scale_fill_viridis_d( name="Signaling \n cell type", option="B",direction = 1)+
  scale_linetype_discrete( name="Malignance", labels=c("Cancer","Non-cancer"))+theme(aspect.ratio=1)+
  labs(y="Fibroblasts activation communications",x="Day")+scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))+facet_wrap(~Pair.Name)
#ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Fibroblast activating signals by cell type and tumor response over time.pdf",width=14, height=14)

partboth[,Malignance:="Cancer"]
partboth[LigandPhenoCelltype!="Cancer cells",Malignance:="Non-cancer"]
partboth[,TumorResponse:="Growing tumor"]
partboth[dynamic_class3!="Non-response",TumorResponse:="Shrinking tumor"]

ggplot(partboth[][  Day==180]%>%group_by(Patient.Study.ID,TumorResponse,HPMR.Ligand,Signal,Malignance,dynamic_class3)%>%summarise(scalelnTransductionMu=mean(scalelnTransductionMu)),#[  LigandPhenoCelltype%in%c("Cancer cells","Endothelial cells","Fibroblasts","Macrophages","Normal epothelial cells","T cells")],  
       aes(y=scalelnTransductionMu,x=HPMR.Ligand ,#col=Malignance ,
           #shape=dynamic_class3,
           fill=Malignance ) ) +
  theme_classic(base_size=26)+
  #facet_grid(~Signal)+ 
  theme(aspect.ratio=1)+#legend.position="none")+
  #geom_smooth(aes(linetype=LigandPhenoCelltype!="Cancer cells"),method="lm",se=T,alpha=0.3)+
  geom_boxplot(aes(group=interaction(HPMR.Ligand,Malignance)))+
  geom_jitter(pch=21,col="black",size=2.5, position = position_dodge(width=0.5)) +
  scale_fill_manual(name="Cell type",values=c(pal_aaas("default")(6)[2],"grey"))+
  labs(y="Fibroblasts differentiation \n communication strength",x="Signaling ligand")

partboth[,LigandPhenoCelltype2:=LigandPhenoCelltype]
partboth[LigandPhenoCelltype2=="Normal epithelial cells",LigandPhenoCelltype2:="Cancer cells"]
partboth1 <- data.table(partboth%>%group_by(Patient.Study.ID,TumorResponse,HPMR.Ligand,Signal,LigandPhenoCelltype2,dynamic_class3,Pair.Name,Day)%>%summarise(scalelnTransductionMu=sum(scalelnTransductionMu)))

partboth1[,Malignance:="Cancer"]
partboth1[LigandPhenoCelltype2!="Cancer cells",Malignance:="Non-cancer"]

ggplot(partboth1[  Day==180]%>%group_by(Patient.Study.ID,TumorResponse,HPMR.Ligand,Signal,LigandPhenoCelltype2,dynamic_class3)%>%summarise(scalelnTransductionMu=mean(scalelnTransductionMu)),#[  LigandPhenoCelltype%in%c("Cancer cells","Endothelial cells","Fibroblasts","Macrophages","Normal epothelial cells","T cells")],  
       aes(y=scalelnTransductionMu,x=HPMR.Ligand ,#col=Malignance ,
           #shape=dynamic_class3,
           fill=LigandPhenoCelltype2 ) ) +
  theme_classic(base_size=26)+
  #facet_grid(~Signal)+ 
  theme(aspect.ratio=1)+#legend.position="none")+
  #geom_smooth(aes(linetype=LigandPhenoCelltype!="Cancer cells"),method="lm",se=T,alpha=0.3)+
  geom_boxplot(aes(group=interaction(HPMR.Ligand,LigandPhenoCelltype2)))+
  geom_jitter(pch=21,col="black",size=2.5, position = position_dodge(width=0.5))


ggplot(partboth1[  Day==180]%>%group_by(Patient.Study.ID,TumorResponse,HPMR.Ligand,Signal,Malignance,dynamic_class3)%>%summarise(scalelnTransductionMu=mean(scalelnTransductionMu)),#[  LigandPhenoCelltype%in%c("Cancer cells","Endothelial cells","Fibroblasts","Macrophages","Normal epothelial cells","T cells")],  
       aes(y=scalelnTransductionMu,x=HPMR.Ligand ,#col=Malignance ,
           #shape=dynamic_class3,
           fill=Malignance ) ) +
  theme_classic(base_size=26)+
  #facet_grid(~Signal)+ 
  theme(aspect.ratio=1)+#legend.position="none")+
  #geom_smooth(aes(linetype=LigandPhenoCelltype!="Cancer cells"),method="lm",se=T,alpha=0.3)+
  geom_boxplot(aes(group=interaction(HPMR.Ligand,Malignance)))+
  geom_jitter(pch=21,col="black",size=2.5, position = position_dodge(width=0.5))+
  scale_fill_manual(name="Cell type",values=c(pal_aaas("default")(6)[2],"grey"))+
  labs(y="Fibroblasts differentiation \n communication strength",x="Signaling ligand")


#ggsave(file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Fibroblast differentiation signals by cell type and tumor response over time.pdf",width=8.5, height=8.5)
#ggsave(file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Validation Fibroblast differentiation signals from cancer and noncancerRed.pdf",width=10, height=10)


validData<-data.table( partboth1[  Day==180]%>%group_by(Patient.Study.ID,TumorResponse,HPMR.Ligand,Signal,Malignance,dynamic_class3)%>%dplyr::summarise(scalelnTransductionMu=mean(scalelnTransductionMu)) )
#write.csv(validData,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/ValidationCohortFibroblastdifferentiationsignalsfromcancerandnoncancer.csv")

validData$Malignance <- factor(validData$Malignance , levels=c("Non-cancer","Cancer"))
summary( lm(scalelnTransductionMu~ -1+HPMR.Ligand:Malignance,data=validData) )


