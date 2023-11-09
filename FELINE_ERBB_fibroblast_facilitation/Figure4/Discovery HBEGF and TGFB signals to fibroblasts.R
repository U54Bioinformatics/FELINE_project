rm(list=ls())
library(glmnet);require(dplyr);require(data.table);require(tidyr)
require(ggraph); require(tidygraph);
require(scales);require(ggsci); require(lmerTest)
load( "~/Dropbox/Cancer_pheno_evo/data/FELINE2/LigandReceptor/Filtered_Human-2015-Ramilowski-LR-pairs.RData" )
load(file ="~/Dropbox/Cancer_pheno_evo/data/FELINE2/Communication output merged ALL ARMS/PopulationCommunicationMerged ALL ARMS.RData" )#CCI,allphenotypes, uu,perIndiv,
LRgenelist <- unique( data.table(merge(CCI[],LRpairsFiltered,by="Pair.Name"))$HPMR.Receptor) 
rm(list= c("LRpairs", "LRpairsFiltered"))

CCI[,Treat:="CombinationRibo"]
CCI[ARM=="A",Treat:="LetrozoleAlone"]
CCI[,scaleTransduction:=scale(TransductionMu,center=F),by=c("Pair.Name")] 
CCI[!is.finite(TransductionMu),scaleTransduction:=0]
CCI[,scalelnTransductionMu:=scale(log(1+TransductionMu)),by=Pair.Name]
CCI[, c("HPMR.Ligand", "HPMR.Receptor") := tstrsplit(Pair.Name, "_", fixed=TRUE)]

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


CCI[,scaleTransduction2:=exp(scale(log(TransductionMu),center=F)),by=c("Pair.Name")] 

load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/GeneOntology/growthFactorReceptors.RData" )


part1<- CCI[
  ][(grepl("TGFB1_",Pair.Name)|grepl("TGFB2_",Pair.Name)|grepl("TGFB3_",Pair.Name))
    &grepl(c("TGFBR"),Pair.Name)][ReceptorPhenoCelltype=="Fibroblasts"][]
part1[,Signal:="TGF"]
part2<- CCI[
  ][grepl("HBEGF_",Pair.Name)][ReceptorPhenoCelltype=="Fibroblasts"][]
part2[,Signal:="HBEGF"]
partboth <- rbind(part1,part2)[!LigandPhenoCelltype%in%c("B cells")]
partboth[LigandPhenoCelltype%in%c("CD4+ T cells","CD8+ T cells" ),LigandPhenoCelltype:="T cells"]
partboth[LigandPhenoCelltype%in%c("Pericytes","Endothelial cells" ),LigandPhenoCelltype:="Endothelial/Pericytes cells"]
ggplot(partboth[  ],
       aes(y=scalelnTransductionMu,x=(log(1+Day)) ,col=LigandPhenoCelltype ,fill=LigandPhenoCelltype) ) + theme_classic(base_size=26)+
  facet_grid(Signal~dynamic_class3)+ 
  geom_smooth(aes(linetype=LigandPhenoCelltype!="Cancer cells"),method="lm",se=T,alpha=0.3)+ geom_point(pch=21,col="black",size=2.5) +
  scale_color_viridis_d(name="Signaling \n cell type", option="B",direction = 1)+
  scale_fill_viridis_d( name="Signaling \n cell type", option="B",direction = 1)+
  scale_linetype_discrete( name="Malignance", labels=c("Cancer","Non-cancer"))+theme(aspect.ratio=1)+
  labs(y="Fibroblasts activation communications",x="Day")+scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))
ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Fibroblast activating signals by cell type and tumor response over time.pdf",width=14, height=14)

partboth[,Malignance:="Cancer"]
partboth[LigandPhenoCelltype!="Cancer cells",Malignance:="Non-cancer"]
partboth[,TumorResponse:="Growing tumor"]
partboth[dynamic_class3!="Non-response",TumorResponse:="Shrinking tumor"]

ggplot(partboth[  Day==180]%>%group_by(Patient.Study.ID,TumorResponse,HPMR.Ligand,Signal,Malignance,dynamic_class3)%>%summarise(scalelnTransductionMu=mean(scalelnTransductionMu)),#[  LigandPhenoCelltype%in%c("Cancer cells","Endothelial cells","Fibroblasts","Macrophages","Normal epothelial cells","T cells")],  
       aes(y=scalelnTransductionMu,x=HPMR.Ligand ,
           fill=Malignance ) ) +
  theme_classic(base_size=26)+
  theme(aspect.ratio=1)+
  geom_boxplot(aes(group=interaction(HPMR.Ligand,Malignance)))+
  geom_jitter(pch=21,col="black",size=2.5, position = position_dodge(width=0.5)) +
  scale_fill_manual(name="Cell type",values=c(pal_aaas("default")(6)[2],"grey"))+
  labs(y="Fibroblasts differentiation \n communication strength",x="Signaling ligand")
#ggsave(file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Fibroblast differentiation signals from cancer and noncancerRed.pdf",width=10, height=10)


discData<-data.table( partboth[  Day==180]%>%group_by(Patient.Study.ID,TumorResponse,HPMR.Ligand,Signal,Malignance,dynamic_class3)%>%dplyr::summarise(scalelnTransductionMu=mean(scalelnTransductionMu)) )
write.csv(discData,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/DiscoveryCohortFibroblastdifferentiationsignalsfromcancerandnoncancer.csv")

discData$Malignance <- factor(discData$Malignance , levels=c("Non-cancer","Cancer"))
summary( lm(scalelnTransductionMu~ -1+HPMR.Ligand:Malignance,data=discData) )

