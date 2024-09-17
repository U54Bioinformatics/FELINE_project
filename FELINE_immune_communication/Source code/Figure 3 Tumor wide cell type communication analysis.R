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

CCI <-  read.csv(file=paste0(savloc, "SourceData_Figure3_CellCummunicationIndexTumorWideSingnaling.csv"))
CCI_discovery <- data.table(CCI)[Cohort=="Discovery"]
CCI_validation <-  data.table(CCI)[Cohort=="Validation"]


#### Discovery cohort analysis
# calculate the average strength of communication in the data
summaryCCI_discovery <- data.table( CCI_discovery[][order(Treat,LigandPhenoCelltype, ReceptorPhenoCelltype ) ] %>% 
                            group_by(Treat,LigandPhenoCelltype, ReceptorPhenoCelltype, dynamic_class3, Day ) %>% 
                            #dplyr::summarise(meanSig= mean(scalelnTransductionMu) ) )
                            dplyr::summarise(meanSig=median(log(1+scaleTransduction)))  ) 

# Exclude the rarely found B cell type
plt_discovery <- summaryCCI_discovery[LigandPhenoCelltype!="B cells"&ReceptorPhenoCelltype!="B cells"]
# Adjust labeling for plotting
plt_discovery[LigandPhenoCelltype=="Normal epithelial cells",LigandPhenoCelltype:= "Diploid epithelial cells"]
plt_discovery[ReceptorPhenoCelltype=="Normal epithelial cells",ReceptorPhenoCelltype:= "Diploid epithelial cells"]
plt_discovery[LigandPhenoCelltype=="Macrophages",LigandPhenoCelltype:= "Myeloid cells"]
plt_discovery[ReceptorPhenoCelltype=="Macrophages",ReceptorPhenoCelltype:= "Myeloid cells"]

# order cell names
plt_discovery$LigandPhenoCelltype <-factor(plt_discovery$LigandPhenoCelltype, 
                                  levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Diploid epithelial cells","Cancer cells","B cells","Myeloid cells","CD8+ T cells","CD4+ T cells"))
plt_discovery$ReceptorPhenoCelltype <-factor(plt_discovery$ReceptorPhenoCelltype, 
                                    levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Diploid epithelial cells","Cancer cells","B cells","Myeloid cells","CD8+ T cells","CD4+ T cells") )
# label treatments
plt_discovery[,Treatmentlab:= "Combination ribociclib"]
plt_discovery[Treat=="LetrozoleAlone",Treatmentlab:= "Letrozole alone"]

ggplot(plt_discovery, aes(y= log(exp(meanSig)-1), x= LigandPhenoCelltype, col= LigandPhenoCelltype,fill=LigandPhenoCelltype,group=interaction(LigandPhenoCelltype)))+geom_violin(scale="width",alpha=0.5)+ geom_point()+
  theme_classic(base_size=26)+theme(axis.text.x=element_text(angle=90,vjust=0.5))+
  theme(aspect.ratio=1,legend.position="none") +labs(y="Contribution to \n tumor wide communication (mean)",x= "Cell type") +
  facet_wrap(~Treatmentlab,nrow=1)+ 
  scale_y_continuous(breaks=log(c(0.001,0.01,0.1,1)),labels=c(0.001,0.01,0.1,1))
ggplot(plt_discovery, aes(y= log(exp(meanSig)-1), x= LigandPhenoCelltype, col= LigandPhenoCelltype,fill=LigandPhenoCelltype,group=interaction(LigandPhenoCelltype)))+geom_violin(scale="width",alpha=0.5)+ geom_point()+
  theme_classic(base_size=26)+theme(axis.text.x=element_text(angle=90,vjust=0.5))+
  theme(aspect.ratio=1,legend.position="none") +labs(y="Contribution to tumor wide communication (mean)",x= "Cell type") +
  facet_wrap(~Treatmentlab,nrow=1)+
  theme(axis.title=element_blank(), axis.text.x = element_blank(), axis.text=element_blank(),strip.text = element_blank(),legend.title=element_blank(),legend.text=element_blank() ,legend.position = "none")+ 
  scale_y_continuous(breaks=log(c(0.001,0.01,0.1,1)),labels=c(0.001,0.01,0.1,1))

plt_discovery[, Ligandiscancer:=F]
plt_discovery[LigandPhenoCelltype=="Cancer cells", Ligandiscancer:=T]
summary(lm( log(exp(meanSig)-1)~ Ligandiscancer , plt_discovery))

plt_discovery[, Receptoriscancer:=F]
plt_discovery[ReceptorPhenoCelltype=="Cancer cells", Receptoriscancer:=T]
summary(lm( log(exp(meanSig)-1)~ Receptoriscancer , plt_discovery))


ggplot(plt_discovery, aes(y= log(exp(meanSig)-1), x= ReceptorPhenoCelltype, col= ReceptorPhenoCelltype,fill=ReceptorPhenoCelltype,group=interaction(ReceptorPhenoCelltype)))+geom_violin(col=NA,scale="width",alpha=0.5)+ geom_point()+
  theme_classic(base_size=26)+
  theme(axis.text.x=element_text(angle=90,vjust=0.5)) +
  theme(aspect.ratio=1,legend.position="none") +labs(y="Strength of communicaiton \n received (mean)",x= "Cell type") +
  facet_wrap(~Treatmentlab)+ 
  scale_y_continuous(breaks=log(c(0.001,0.01,0.1,1)),labels=c(0.001,0.01,0.1,1))

ggplot(plt_discovery, aes(y= log(exp(meanSig)-1), x= ReceptorPhenoCelltype, col= ReceptorPhenoCelltype,fill=ReceptorPhenoCelltype,group=interaction(ReceptorPhenoCelltype)))+geom_violin(col=NA,scale="width",alpha=0.5)+ geom_point()+theme_classic(base_size=16)+
  theme(axis.text.x=element_text(angle=90,vjust=0.5)) +
  theme(aspect.ratio=1,legend.position="none") +labs(y="Strength of communicaiton received (mean)",x= "Cell type") +
  facet_wrap(~Treatmentlab)+
  theme(axis.title=element_blank(), axis.text.x = element_blank(), axis.text=element_blank(),strip.text = element_blank(),legend.title=element_blank(),legend.text=element_blank() ,legend.position = "none")







## Duplicate analysis for Validation cohort
# calculate the average strength of communication in the actual data
summaryCCI_validation <- data.table( CCI_validation[
  order(LigandPhenoCelltype, ReceptorPhenoCelltype ) ] %>% 
    group_by(Treat,LigandPhenoCelltype, ReceptorPhenoCelltype, dynamic_class3, Day ) %>% 
    #dplyr::summarise(meanSig= mean(scalelnTransductionMu) ) )
    dplyr::summarise(meanSig=median(log(1+scaleTransduction)))  ) 


plt_validation <- summaryCCI_validation[LigandPhenoCelltype!="B cells"&ReceptorPhenoCelltype!="B cells"]

plt_validation[LigandPhenoCelltype=="Normal epithelial cells",LigandPhenoCelltype:= "Diploid epithelial cells"]
plt_validation[ReceptorPhenoCelltype=="Normal epithelial cells",ReceptorPhenoCelltype:= "Diploid epithelial cells"]
plt_validation[LigandPhenoCelltype=="Macrophages",LigandPhenoCelltype:= "Myeloid cells"]
plt_validation[ReceptorPhenoCelltype=="Macrophages",ReceptorPhenoCelltype:= "Myeloid cells"]

plt_validation$LigandPhenoCelltype <-factor(plt_validation$LigandPhenoCelltype, 
                                  levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Diploid epithelial cells","Cancer cells","B cells","Myeloid cells","CD8+ T cells","CD4+ T cells"))
plt_validation$ReceptorPhenoCelltype <-factor(plt_validation$ReceptorPhenoCelltype, 
                                    levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Diploid epithelial cells","Cancer cells","B cells","Myeloid cells","CD8+ T cells","CD4+ T cells") )
plt_validation[,Treatmentlab:= "Combination ribociclib"]
plt_validation[Treat=="LetrozoleAlone",Treatmentlab:= "Letrozole alone"]


ggplot(plt_validation,
       aes(y= log(exp(meanSig)-1), x= LigandPhenoCelltype, col= LigandPhenoCelltype,fill=LigandPhenoCelltype,group=interaction(LigandPhenoCelltype)))+geom_violin(scale="width",alpha=0.5)+ geom_point()+
  theme_classic(base_size=26)+theme(axis.text.x=element_text(angle=90,vjust=0.5))+
  theme(aspect.ratio=1,legend.position="none") +labs(y="Contribution to \n tumor wide communication (mean)",x= "Cell type") +
  facet_wrap(~Treatmentlab,nrow=1)+ 
  scale_y_continuous(breaks=log(c(0.001,0.01,0.1,1)),labels=c(0.001,0.01,0.1,1))

ggplot(plt_validation, aes(y= log(exp(meanSig)-1), x= ReceptorPhenoCelltype, col= ReceptorPhenoCelltype,fill=ReceptorPhenoCelltype,group=interaction(ReceptorPhenoCelltype)))+geom_violin(col=NA,scale="width",alpha=0.5)+ geom_point()+
  theme_classic(base_size=26)+
  theme(axis.text.x=element_text(angle=90,vjust=0.5)) +
  theme(aspect.ratio=1,legend.position="none") +labs(y="Strength of communicaiton \n received (mean)",x= "Cell type") +
  facet_wrap(~Treatmentlab)+ 
  scale_y_continuous(breaks=log(c(0.001,0.01,0.1,1)),labels=c(0.001,0.01,0.1,1))

#save(plt_validation, file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/temporary/whocommValidation.RData")


plt1<- rbind(data.table(Cohort="Discovery",plt_discovery)%>%dplyr::select(-c(Ligandiscancer, Receptoriscancer)),data.table(Cohort="Validation",plt_validation))

ggplot(plt1, aes(y= log(exp(meanSig)-1), x= LigandPhenoCelltype, col= LigandPhenoCelltype,fill=LigandPhenoCelltype,group=interaction(LigandPhenoCelltype)))+
  geom_boxplot(alpha=0.5,outlier.color =NA)+ 
  stat_boxplot(geom="errorbar")+
  geom_jitter(aes(shape=Cohort),size=1.5,height=0,width=0.2)+
  theme_classic(base_size=26)+theme(axis.text.x=element_text(angle=90,vjust=0.5))+
  theme(aspect.ratio=1,legend.position="none") +labs(y="Contribution to \n tumor wide communication \n (mean)",x= "Cell type") +
  facet_wrap(~Treatmentlab,nrow=1)+ 
  scale_y_continuous(breaks=log(c(0.001,0.01,0.1,1)),labels=c(0.001,0.01,0.1,1))

ggplot(plt1, aes(y= log(exp(meanSig)-1), x= ReceptorPhenoCelltype, col= ReceptorPhenoCelltype,fill=ReceptorPhenoCelltype,group=interaction(ReceptorPhenoCelltype)))+
  geom_boxplot(alpha=0.5,outlier.color =NA)+
  stat_boxplot(geom="errorbar")+
  geom_jitter(aes(shape=Cohort),size=1.5,height=0,width=0.2)+
  theme_classic(base_size=26)+
  theme(axis.text.x=element_text(angle=90,vjust=0.5)) +
  theme(aspect.ratio=1,legend.position="none") +labs(y="Strength of communicaiton \n received (mean)",x= "Cell type") +
  facet_wrap(~Treatmentlab)+ 
  scale_y_continuous(breaks=log(c(0.001,0.01,0.1,1)),labels=c(0.001,0.01,0.1,1))

plt1[, Ligandiscancer:=F]
plt1[LigandPhenoCelltype=="Cancer cells", Ligandiscancer:=T]
summary(lm( log(exp(meanSig)-1)~  Ligandiscancer*Cohort , plt1[Treatmentlab=="Combination ribociclib"]))
summary(lm( log(exp(meanSig)-1)~  Ligandiscancer*Cohort , plt1[Treatmentlab=="Letrozole alone"]))

summary(lm( log(exp(meanSig)-1)~  Ligandiscancer , plt1[Treatmentlab=="Combination ribociclib"]))
summary(lm( log(exp(meanSig)-1)~  Ligandiscancer , plt1[Treatmentlab=="Letrozole alone"]))

plt1[, Receptoriscancer:=F]
plt1[ReceptorPhenoCelltype=="Cancer cells", Receptoriscancer:=T]
summary(lm( log(exp(meanSig)-1)~ Receptoriscancer*Cohort , plt1[Treatmentlab=="Combination ribociclib"]))
summary(lm( log(exp(meanSig)-1)~ Receptoriscancer*Cohort , plt1[Treatmentlab=="Letrozole alone"]))
summary(lm( log(exp(meanSig)-1)~  Receptoriscancer , plt1[Treatmentlab=="Combination ribociclib"]))
summary(lm( log(exp(meanSig)-1)~  Receptoriscancer , plt1[Treatmentlab=="Letrozole alone"]))

setwd(savloc)
#write.csv(plt1, file="SourceData_Figure3_Discovery and Validation cell type signal contribution_Output.csv")









# SI verificatiom analyses
CCI_all<-data.table(CCI)
# FigS5) select cell type communications known to be specifically recieved by a given cell type and verify that this is observed in the results
# ERBB  
growthFactorReceptors2 <- unique(CCI_discovery[grepl("ERBB",R)]$R) #growthFactorReceptors2 <- c("EGFR","ERBB2","ERBB3","ERBB4")
CCI_discovery[,isGrowthFactor:=F]
CCI_discovery[R%in%growthFactorReceptors2, isGrowthFactor:=T]
# calculate the average strength of communication in the data
summaryCCI <- data.table( CCI_discovery[][order(Treat,LigandPhenoCelltype, ReceptorPhenoCelltype ) ] %>% 
                            group_by(Treat,LigandPhenoCelltype, ReceptorPhenoCelltype, dynamic_class3, Day,isGrowthFactor ) %>% 
                            #dplyr::summarise(meanSig= mean(scalelnTransductionMu) ) )
                            dplyr::summarise(meanSig=mean(log(1+scaleTransduction)))  )

plt1 <- summaryCCI[LigandPhenoCelltype!="B cells"&ReceptorPhenoCelltype!="B cells"]
plt1[LigandPhenoCelltype=="Normal epithelial cells",LigandPhenoCelltype:= "Diploid epithelial cells"]
plt1[ReceptorPhenoCelltype=="Normal epithelial cells",ReceptorPhenoCelltype:= "Diploid epithelial cells"]
plt1[LigandPhenoCelltype=="Macrophages",LigandPhenoCelltype:= "Myeloid cells"]
plt1[ReceptorPhenoCelltype=="Macrophages",ReceptorPhenoCelltype:= "Myeloid cells"]
plt1$LigandPhenoCelltype <-factor(plt1$LigandPhenoCelltype, 
                                  levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Diploid epithelial cells","Cancer cells","B cells","Myeloid cells","CD8+ T cells","CD4+ T cells"))
plt1$ReceptorPhenoCelltype <-factor(plt1$ReceptorPhenoCelltype, 
                                    levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Diploid epithelial cells","Cancer cells","B cells","Myeloid cells","CD8+ T cells","CD4+ T cells") )

plt1[, Receptoriscancer:=F]
plt1[ReceptorPhenoCelltype=="Cancer cells", Receptoriscancer:=T]
summary(lm( log(exp(meanSig)-1)~ Receptoriscancer , plt1[isGrowthFactor==T]))



ggplot(plt1[isGrowthFactor==T], aes(y= log(exp(meanSig)-1), x= ReceptorPhenoCelltype, col= ReceptorPhenoCelltype,fill=ReceptorPhenoCelltype,group=interaction(ReceptorPhenoCelltype)))+geom_violin(col=NA,scale="width",alpha=0.5)+ geom_point()+theme_classic(base_size=16)+
  theme(axis.text.x=element_text(angle=90,vjust=0.5)) +
  theme(aspect.ratio=1,legend.position="none") +labs(y="Strength of ERBB growth factor \n communication received (mean)",x= "Cell type") +
  facet_wrap(~Treat)+  scale_y_continuous(breaks=log(c(0.001,0.01,0.1,1)),labels=c(0.001,0.01,0.1,1))
#ggsave(filename = "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Strength communication cell types AllArms/ERBB family specific Cell type strength of communicaiton received by treatment.png")

ERBBout<-plt1[isGrowthFactor==T]




#macrophage colony-stimulating factor receptor
growthFactorReceptors2 <- c("CSF1R","CSF2RA")#,"FGFR1","FGFR2","FGFR3","FGFR4","FGFRL1" ,"TGFBR1","TGFBR2" ,"TGFBR3"   )
CCI_discovery[,isGrowthFactor:=F]
CCI_discovery[R%in%growthFactorReceptors2, isGrowthFactor:=T]
summaryCCI <- data.table( CCI_discovery[][order(Treat,LigandPhenoCelltype, ReceptorPhenoCelltype ) ] %>% 
                            group_by(Treat,LigandPhenoCelltype, ReceptorPhenoCelltype, dynamic_class3, Day,isGrowthFactor ) %>% 
                            dplyr::summarise(meanSig=mean(log(1+scaleTransduction)))  )

plt1 <- summaryCCI[LigandPhenoCelltype!="B cells"&ReceptorPhenoCelltype!="B cells"]
plt1[LigandPhenoCelltype=="Normal epithelial cells",LigandPhenoCelltype:= "Diploid epithelial cells"]
plt1[ReceptorPhenoCelltype=="Normal epithelial cells",ReceptorPhenoCelltype:= "Diploid epithelial cells"]
plt1[LigandPhenoCelltype=="Macrophages",LigandPhenoCelltype:= "Myeloid cells"]
plt1[ReceptorPhenoCelltype=="Macrophages",ReceptorPhenoCelltype:= "Myeloid cells"]
plt1$LigandPhenoCelltype <-factor(plt1$LigandPhenoCelltype, 
                                  levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Diploid epithelial cells","Cancer cells","B cells","Myeloid cells","CD8+ T cells","CD4+ T cells"))
plt1$ReceptorPhenoCelltype <-factor(plt1$ReceptorPhenoCelltype, 
                                    levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Diploid epithelial cells","Cancer cells","B cells","Myeloid cells","CD8+ T cells","CD4+ T cells") )



ggplot(plt1[isGrowthFactor==T], aes(y= log(exp(meanSig)-1), x= ReceptorPhenoCelltype, col= ReceptorPhenoCelltype,fill=ReceptorPhenoCelltype,group=interaction(ReceptorPhenoCelltype)))+geom_violin(col=NA,scale="width",alpha=0.5)+ geom_point()+theme_classic(base_size=16)+
  theme(axis.text.x=element_text(angle=90,vjust=0.5)) +
  theme(aspect.ratio=1,legend.position="none") +labs(y="Strength of macrophage colony stimulating factor \n communication received (mean)",x= "Cell type") +
  facet_wrap(~Treat)+
  scale_y_continuous(breaks=log(c(0.001,0.01,0.1,1)),labels=c(0.001,0.01,0.1,1))
#ggsave(filename = "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Strength communication cell types AllArms/Macrophage CSF specific Cell type strength of communicaiton received by treatment.png")

CSFRout<-plt1[isGrowthFactor==T]


###CCR5
growthFactorReceptors2 <- c("CCR5") 
CCI_discovery[,isGrowthFactor:=F]
CCI_discovery[R%in%growthFactorReceptors2, isGrowthFactor:=T]
summaryCCI <- data.table( CCI_discovery[][order(Treat,LigandPhenoCelltype, ReceptorPhenoCelltype ) ] %>% 
                            group_by(Treat,LigandPhenoCelltype, ReceptorPhenoCelltype, dynamic_class3, Day,isGrowthFactor ) %>% 
                            dplyr::summarise(meanSig=mean(log(1+scaleTransduction)))  )

plt1 <- summaryCCI[LigandPhenoCelltype!="B cells"&ReceptorPhenoCelltype!="B cells"]
plt1[LigandPhenoCelltype=="Normal epithelial cells",LigandPhenoCelltype:= "Diploid epithelial cells"]
plt1[ReceptorPhenoCelltype=="Normal epithelial cells",ReceptorPhenoCelltype:= "Diploid epithelial cells"]
plt1[LigandPhenoCelltype=="Macrophages",LigandPhenoCelltype:= "Myeloid cells"]
plt1[ReceptorPhenoCelltype=="Macrophages",ReceptorPhenoCelltype:= "Myeloid cells"]
plt1$LigandPhenoCelltype <-factor(plt1$LigandPhenoCelltype, 
                                  levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Diploid epithelial cells","Cancer cells","B cells","Myeloid cells","CD8+ T cells","CD4+ T cells"))
plt1$ReceptorPhenoCelltype <-factor(plt1$ReceptorPhenoCelltype, 
                                    levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Diploid epithelial cells","Cancer cells","B cells","Myeloid cells","CD8+ T cells","CD4+ T cells") )



ggplot(plt1[isGrowthFactor==T], aes(y= log(exp(meanSig)-1), x= ReceptorPhenoCelltype, col= ReceptorPhenoCelltype,fill=ReceptorPhenoCelltype,group=interaction(ReceptorPhenoCelltype)))+geom_violin(col=NA,scale="width",alpha=0.5)+ geom_point()+theme_classic(base_size=16)+
  theme(axis.text.x=element_text(angle=90,vjust=0.5)) +
  theme(aspect.ratio=1,legend.position="none") +labs(y="Strength of CCR5 \n communication received (mean)",x= "Cell type") +
  facet_wrap(~Treat)+
  scale_y_continuous(breaks=log(c(0.001,0.01,0.1,1)),labels=c(0.001,0.01,0.1,1))

#ggsave(filename = "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Strength communication cell types AllArms/CCR5 specific Cell type strength of communicaiton received by treatment.png")

CCRout<-plt1[isGrowthFactor==T]



## VEGF receptor : FLT1
growthFactorReceptors2 <- c("FLT1")#,"FGFR1","FGFR2","FGFR3","FGFR4","FGFRL1" ,"TGFBR1","TGFBR2" ,"TGFBR3"   )

CCI_discovery[,isGrowthFactor:=F]
CCI_discovery[R%in%growthFactorReceptors2, isGrowthFactor:=T]
# calculate the average strength of communication in the data
summaryCCI <- data.table( CCI_discovery[][order(Treat,LigandPhenoCelltype, ReceptorPhenoCelltype ) ] %>% 
                            group_by(Treat,LigandPhenoCelltype, ReceptorPhenoCelltype, dynamic_class3, Day,isGrowthFactor ) %>% 
                            #dplyr::summarise(meanSig= mean(scalelnTransductionMu) ) )
                            dplyr::summarise(meanSig=mean(log(1+scaleTransduction)))  )

plt1 <- summaryCCI[LigandPhenoCelltype!="B cells"&ReceptorPhenoCelltype!="B cells"]
plt1[LigandPhenoCelltype=="Normal epithelial cells",LigandPhenoCelltype:= "Diploid epithelial cells"]
plt1[ReceptorPhenoCelltype=="Normal epithelial cells",ReceptorPhenoCelltype:= "Diploid epithelial cells"]
plt1[LigandPhenoCelltype=="Macrophages",LigandPhenoCelltype:= "Myeloid cells"]
plt1[ReceptorPhenoCelltype=="Macrophages",ReceptorPhenoCelltype:= "Myeloid cells"]
plt1$LigandPhenoCelltype <-factor(plt1$LigandPhenoCelltype, 
                                  levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Diploid epithelial cells","Cancer cells","B cells","Myeloid cells","CD8+ T cells","CD4+ T cells"))
plt1$ReceptorPhenoCelltype <-factor(plt1$ReceptorPhenoCelltype, 
                                    levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Diploid epithelial cells","Cancer cells","B cells","Myeloid cells","CD8+ T cells","CD4+ T cells") )



ggplot(plt1[isGrowthFactor==T], aes(y= log(exp(meanSig)-1), x= ReceptorPhenoCelltype, col= ReceptorPhenoCelltype,fill=ReceptorPhenoCelltype,group=interaction(ReceptorPhenoCelltype)))+geom_violin(col=NA,scale="width",alpha=0.5)+ geom_point()+theme_classic(base_size=16)+
  theme(axis.text.x=element_text(angle=90,vjust=0.5)) +
  theme(aspect.ratio=1,legend.position="none") +labs(y="Strength of VEGF \n communication received (mean)",x= "Cell type") +
  facet_wrap(~Treat)+
  scale_y_continuous(breaks=log(c(0.001,0.01,0.1,1)),labels=c(0.001,0.01,0.1,1))
#ggsave(filename = "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Strength communication cell types AllArms/VEGF specific Cell type strength of communicaiton received by treatment.png")

VEGout<-plt1[isGrowthFactor==T]


## FGFR1
growthFactorReceptors2 <- c("FGFR")
CCI_discovery[,isGrowthFactor:=F]
CCI_discovery[grep("FGFR1",R), isGrowthFactor:=T]
# calculate the average strength of communication in the data
summaryCCI <- data.table( CCI_discovery[][order(Treat,LigandPhenoCelltype, ReceptorPhenoCelltype ) ] %>% 
                            group_by(Treat,LigandPhenoCelltype, ReceptorPhenoCelltype, dynamic_class3, Day,isGrowthFactor ) %>% 
                            #dplyr::summarise(meanSig= mean(scalelnTransductionMu) ) )
                            dplyr::summarise(meanSig=mean(log(1+scaleTransduction)))  )

plt1 <- summaryCCI[LigandPhenoCelltype!="B cells"&ReceptorPhenoCelltype!="B cells"]
plt1[LigandPhenoCelltype=="Normal epithelial cells",LigandPhenoCelltype:= "Diploid epithelial cells"]
plt1[ReceptorPhenoCelltype=="Normal epithelial cells",ReceptorPhenoCelltype:= "Diploid epithelial cells"]
plt1[LigandPhenoCelltype=="Macrophages",LigandPhenoCelltype:= "Myeloid cells"]
plt1[ReceptorPhenoCelltype=="Macrophages",ReceptorPhenoCelltype:= "Myeloid cells"]
plt1$LigandPhenoCelltype <-factor(plt1$LigandPhenoCelltype, 
                                  levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Diploid epithelial cells","Cancer cells","B cells","Myeloid cells","CD8+ T cells","CD4+ T cells"))
plt1$ReceptorPhenoCelltype <-factor(plt1$ReceptorPhenoCelltype, 
                                    levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Diploid epithelial cells","Cancer cells","B cells","Myeloid cells","CD8+ T cells","CD4+ T cells") )



ggplot(plt1[isGrowthFactor==T], aes(y= log(exp(meanSig)-1), x= ReceptorPhenoCelltype, col= ReceptorPhenoCelltype,fill=ReceptorPhenoCelltype,group=interaction(ReceptorPhenoCelltype)))+geom_violin(col=NA,scale="width",alpha=0.5)+ geom_point()+theme_classic(base_size=16)+
  theme(axis.text.x=element_text(angle=90,vjust=0.5)) +
  theme(aspect.ratio=1,legend.position="none") +labs(y="Strength of FGFR1 \n communication received (mean)",x= "Cell type") +
  facet_wrap(~Treat)+
  scale_y_continuous(breaks=log(c(0.001,0.01,0.1,1)),labels=c(0.001,0.01,0.1,1))
#ggsave(filename = "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Strength communication cell types AllArms/FGFR1 specific Cell type strength of communicaiton received by treatment.png")

FGFout<-plt1[isGrowthFactor==T]

outtable<-rbind(ERBBout%>%dplyr::select(-Receptoriscancer),
CSFRout,
CCRout,
VEGout,
FGFout)

savloc<- "/Users/jason/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/SI data/FigureS5/"
write.csv(outtable,file=paste0(savloc,"Figure S5 Recovering known cell type communications.csv"))
