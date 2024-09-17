rm(list=ls())   
require(rTensor); require(nnTensor); library(abind)
require(data.table); require(dplyr); require(ggplot2); require(tidyr)
require(igraph)



### Load Ligand Receptor database list of Ramilowski et al 2015
LRpairsFiltered <- data.table(read.csv( file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/SI data/FigureS4/SourceData_Figure S4 LigandReceptorPairsRamilowski2015.csv"))
LRgenelist <- unique( c(LRpairsFiltered$HPMR.Receptor, LRpairsFiltered$HPMR.Ligand) )

### Load GO database list of GF receptors
#growthFactorReceptors2 <- read.csv(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/SI data/FigureS4/Figure S4 GeneOntologyGrowthFactorReceptors.csv")$x

### Load some signalling data for downstream analysis
savloc<-"/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/Figure3/"
CCI_all <-  read.csv(file=paste0(savloc, "SourceData_Figure3_CellCummunicationIndexTumorWideSingnaling.csv"))
CCI_discovery <- data.table(CCI_all )[Cohort=="Discovery"]
CCI_validation <-  data.table(CCI_all )[Cohort=="Validation"]

CCI <- CCI_discovery

#extended.signiftable<- data.table(read.csv( file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Final communication files/Cancer to macropage communication supervised analysis/Version2 Cancer to macropage M2 stimulating communication supervised analysis.csv"))
extendelistC1<-c("CLCF1_CRLF1","CLCF1_IL6ST","CSF1_CSF1R","FGF17_FGFR1","FGF18_FGFR1",  
 "FGF1_FGFR1" ,   "FGF23_FGFR1" ,  "FGF2_CD44",     "FGF2_FGFR1" ,   "FGF2_NRP1"  ,  
 "FGF2_SDC1"  ,   "FGF2_SDC2"  ,   "FGF2_SDC3" ,    "FGF2_SDC4"  ,   "FGF7_FGFR1" ,  
  "GDF9_ACVR2A" ,  "GDF9_BMPR1A" ,  "GDF9_BMPR1B" ,  "GDF9_BMPR2" ,   "GDF9_TGFBR1"  ,
  "IL11_IL11RA",   "IL11_IL6ST"  ,  "IL12A_IL12RB1", "IL12A_IL12RB2", "IL5_CSF2RB"  , 
 "MADCAM1_CD44" , "MADCAM1_ITGA4" ,"PLAU_IGF2R"  ,  "PLAU_LRP1"  ,   "PLAU_LRP2" ,   
 "PLAU_PLAUR")

corruptsigslong<-CCI[Pair.Name%in%extendelistC1][ReceptorPhenoCelltype=="Macrophages"]#[Treat=="CombinationRibo"]
corruptsigs <- data.table(corruptsigslong%>%select(Pair.Name, Treat,LigandPhenoCelltype,scalelnTransductionMu, Patient.Study.ID,dynamic_class3,Day)%>% spread(Pair.Name,scalelnTransductionMu,fill=0))[order(dynamic_class3,)]

macrosignals<-CCI[Pair.Name%in%extendelistC1][ReceptorPhenoCelltype=="Macrophages"][LigandPhenoCelltype!="B cells"]#[Treat=="CombinationRibo"]
macrosignals[,broadCat:=LigandPhenoCelltype]
macrosignals[grep("T cell",LigandPhenoCelltype) ,broadCat:="T cells"]

macrosignals$broadCat <- factor(macrosignals$broadCat, levels=c("Cancer cells", "Normal epithelial cells","Pericytes","Endothelial cells", "Macrophages"  ,"T cells",
                                                                "Fibroblasts","Adipocytes"          
                                                                 ) )
avmacrosignals <- data.table(macrosignals[]%>%group_by(broadCat,Treat,dynamic_class3)%>%summarise(m=mean(scalelnTransductionMu),
                                                                                            m2=median(scaleTransduction),
                                                                             ucl=mean(scalelnTransductionMu)+sd(scalelnTransductionMu)/sqrt(n()), mean(scalelnTransductionMu)+sd(scalelnTransductionMu)/sqrt(n()),
                                                                             lcl=mean(scalelnTransductionMu)-sd(scalelnTransductionMu)/sqrt(n()),
                                                                             ))
avmacrosignals2 <- data.table(macrosignals[]%>%group_by(Day,broadCat,Treat,dynamic_class3)%>%summarise(m=mean(scalelnTransductionMu),
                                                                                            m2=median(scaleTransduction),
                                                                                            m3=mean(scaleTransduction),
                                                                                            ucl=mean(scalelnTransductionMu)+sd(scalelnTransductionMu)/sqrt(n()),
                                                                                            lcl=mean(scalelnTransductionMu)-sd(scalelnTransductionMu)/sqrt(n()),
))

avmacrosignals2[,"Tumorresponse":="Sensitive"]
avmacrosignals2[dynamic_class3=="Non-response","Tumorresponse":="Resistant"]
avmacrosignals2[,"Treatlabel":="Letrozole alone"]
avmacrosignals2[Treat=="CombinationRibo","Treatlabel":="Combination ribociclib"]

avmacrosignals2[broadCat=="Macrophages","broadCat":="Myeloid cells"]
avmacrosignals2[broadCat=="Normal epithelial cells","broadCat":="Diploid epithelial cells"]

ggplot( avmacrosignals2[Day==0],aes(y=m2,x=Tumorresponse,fill=broadCat))+geom_bar(stat="identity")+ theme_classic(base_size=14)+
  labs(y="Communications stimulating \n M2 myeloid  differentiation",x="Tumor response")+
  scale_fill_discrete(name="Sending \n cell \n type")+theme(aspect.ratio=1)+facet_wrap(~Treatlabel)
#ggsave(file="~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/cancer to macrophage communications/Version2b Ribociclib and Letrozole supervised list of cancer signals to macrophages.pdf")
pthfig<-"/Users/jason/Dropbox/FELINE Project (1)/Manuscript  Feline immune communication/Figures Communication Project/Paper figures Cancer immune communication/"
ggsave(file=paste0(pthfig,"Ribo and Letrozole SI cancer contribution of M2 differentiation signals.png"))

resout <- avmacrosignals2[Day==0]%>%dplyr::select( broadCat , Treat , Treatlabel,Day,dynamic_class3,Tumorresponse,m2)
                                                
ggplot( resout,aes(y=m2,x=Tumorresponse,fill=broadCat))+geom_bar(stat="identity")+ theme_classic(base_size=14)+
  labs(y="Communications stimulating \n M2 myeloid  differentiation",x="Tumor response")+
  scale_fill_discrete(name="Sending \n cell \n type")+theme(aspect.ratio=1)+facet_wrap(~Treatlabel)

savloc <- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/SI data/FigureS13/"
write.csv(resout, file=paste0(savloc,"SourceData_FigureS13_CancerCellMajorContributorOfM2DifferentiationSignals.csv"))




#0.66328706
avmacrosignals2[Day==0][Treatlabel=="Combination ribociclib"][Tumorresponse=="Resistant"][broadCat=="Cancer cells"]$m2
avmacrosignals2[Day==0][Treatlabel=="Combination ribociclib"][Tumorresponse!="Resistant"]$m2%>%sum()



