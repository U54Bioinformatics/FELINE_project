rm(list=ls())   
require(rTensor); require(nnTensor); library(abind)
require(data.table); require(dplyr); require(ggplot2); require(tidyr)
require(igraph)

### Phenotype data
load(file="~/Dropbox/Cancer_pheno_evo/data/FELINE2/PhenotypesOfAllCellTypesAllArms/PhenotypesOfAllCellTypesAllArms.RData") #allphenotypes,UMAPlocs ,UMAPfiles,umapDImRedloc,umapDImRedfiles,nCellTypes,
# Encode the subtypes that have been analysed using UMAP
tmp <- data.table( allphenotypes %>% group_by(key_) %>% slice(1) ) %>% dplyr::select(Celltype , Celltype_subtype) %>% unique
tmp[ , PhenoCelltype:= Celltype]
tmp[Celltype_subtype %in% c("CD4+ T cells", "Tregs"), PhenoCelltype:= "CD4+ T cells"]
tmp[Celltype_subtype %in% c("CD8+ T cells", "NK cells"), PhenoCelltype:= "CD8+ T cells"]
tmp[Celltype_subtype %in% c("Macrophages", "DC", "Monocytes"), PhenoCelltype:= "Macrophages"]
tmp[Celltype_subtype %in% c("Vas-Endo", "Lym-Endo","Endothelial cells"), PhenoCelltype:= "Endothelial cells"]
tmp[Celltype_subtype %in% c("Cancer cells"), PhenoCelltype:= "Cancer cells"]
tmp[Celltype_subtype %in% c("Normal epithelial cells"), PhenoCelltype:= "Normal epithelial cells"]
allphenotypes <- merge(allphenotypes %>% dplyr::select(-c( paste0("V", 1:5), "nCount_RNA", "nFeature_RNA", "Infercnv_CNA", "Platform","Timepoint","Sample_p_t","prop_change", "file_string", "day_fact" , "Ribo","TreatLab", "Burden_t0" ,"BurdenTracked" ,"Day0", "DayLastBurd", "Dose_lab","TreatCode","TreatCodeOrd","dynamic_class2","rgr_A","rgr_B")), tmp, by= c("Celltype", "Celltype_subtype") )

### Unique clusters of cells (ALL) and their umap discretization level
uu <- unique( allphenotypes %>% dplyr::select( c("Celltype","PhenoCelltype", "key_", paste0("Disc_V", 1:5) ) ) ) 

### Load Ligand Receptor database list of Ramilowski et al 2015
load( "~/Dropbox/Cancer_pheno_evo/data/FELINE2/LigandReceptor/Filtered_Human-2015-Ramilowski-LR-pairs.RData")
LRgenelist <- unique( c(LRpairsFiltered$HPMR.Receptor, LRpairsFiltered$HPMR.Ligand) )

### Load some signalling data for downstream analysis
savelocCCI <- "~/Dropbox/Cancer_pheno_evo/data/FELINE2/CommunicationOutputAllArms/"

# Are we studying the per indiv or per cell type level of signalling?
perIndiv=FALSE

# Names of LR data to analyze
filenamesCCI <- list.files(savelocCCI)#%>%length

# First let's do this for the population level
load(file ="~/Dropbox/Cancer_pheno_evo/data/FELINE2/Communication output merged ALL ARMS/PopulationCommunicationMerged ALL ARMS.RData" )#CCI,allphenotypes, uu,perIndiv,
CCI[,Treat:="CombinationRibo"]
CCI[ARM=="A",Treat:="LetrozoleAlone"]
CCI[,scaleTransduction:=scale(TransductionMu,center=F),by=c("Pair.Name")] #"ReceptorPhenoCelltype",
CCI[!is.finite(TransductionMu),scaleTransduction:=0]
CCI[,scalelnTransductionMu:=scale(log(1+TransductionMu)),by=Pair.Name]
CCI[,c("L","R"):= tstrsplit(Pair.Name,"_",fixed=T)]

load(file="~/Dropbox/Cancer_pheno_evo/data/FELINE2/GeneOntology/growthFactorReceptors.RData" )
growthFactorReceptors2<- c("EGFR","ERBB2","ERBB3","ERBB4","FGFR1","FGFR2","FGFR3","FGFR4","FGFRL1" ,"TGFBR1","TGFBR2" ,"TGFBR3"  )
growthFactorReceptors2<- growthFactorReceptors #c("ERBB4" )
CCI[,isGrowthFactor:=F]
CCI[R%in%growthFactorReceptors2,isGrowthFactor:=T]

####
# calculate the average strength of communication in the actual data
summaryCCI <- data.table( CCI[][order(Treat,LigandPhenoCelltype, ReceptorPhenoCelltype ) ] %>% 
                            group_by(Treat,LigandPhenoCelltype, ReceptorPhenoCelltype, dynamic_class3, Day ) %>% 
                            dplyr::summarise(meanSig=median(log(1+scaleTransduction)))  ) 

plt1 <- summaryCCI[LigandPhenoCelltype!="B cells"&ReceptorPhenoCelltype!="B cells"]
plt1[LigandPhenoCelltype=="Normal epithelial cells",LigandPhenoCelltype:= "Diploid epithelial cells"]
plt1[ReceptorPhenoCelltype=="Normal epithelial cells",ReceptorPhenoCelltype:= "Diploid epithelial cells"]
plt1[LigandPhenoCelltype=="Macrophages",LigandPhenoCelltype:= "Myeloid cells"]
plt1[ReceptorPhenoCelltype=="Macrophages",ReceptorPhenoCelltype:= "Myeloid cells"]



plt1$LigandPhenoCelltype <-factor(plt1$LigandPhenoCelltype, 
         levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Diploid epithelial cells","Cancer cells","B cells","Myeloid cells","CD8+ T cells","CD4+ T cells"))
plt1$ReceptorPhenoCelltype <-factor(plt1$ReceptorPhenoCelltype, 
         levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Diploid epithelial cells","Cancer cells","B cells","Myeloid cells","CD8+ T cells","CD4+ T cells") )
plt1[,Treatmentlab:= "Combination ribociclib"]
plt1[Treat=="LetrozoleAlone",Treatmentlab:= "Letrozole alone"]

ggplot(plt1, aes(y= log(exp(meanSig)-1), x= LigandPhenoCelltype, col= LigandPhenoCelltype,fill=LigandPhenoCelltype,group=interaction(LigandPhenoCelltype)))+geom_violin(scale="width",alpha=0.5)+ geom_point()+
  theme_classic(base_size=26)+theme(axis.text.x=element_text(angle=90,vjust=0.5))+
  theme(aspect.ratio=1,legend.position="none") +labs(y="Contribution to \n tumor wide communication (mean)",x= "Cell type") +
  facet_wrap(~Treatmentlab,nrow=1)+ 
  scale_y_continuous(breaks=log(c(0.001,0.01,0.1,1)),labels=c(0.001,0.01,0.1,1))
#ggsave(filename = "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Strength communication cell types AllArms/Cell type contribution to tumor wide communication by treatment.png")
paperfile<- "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
paperfile<- "~/Dropbox/FELINE Project (1)/Manuscript  Feline immune communication/Figures Communication Project/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Ribo and Letrozole Cell type contribution to tumor wide communication.png"),height=10,width=10)
ggplot(plt1, aes(y= log(exp(meanSig)-1), x= LigandPhenoCelltype, col= LigandPhenoCelltype,fill=LigandPhenoCelltype,group=interaction(LigandPhenoCelltype)))+geom_violin(scale="width",alpha=0.5)+ geom_point()+
  theme_classic(base_size=26)+theme(axis.text.x=element_text(angle=90,vjust=0.5))+
  theme(aspect.ratio=1,legend.position="none") +labs(y="Contribution to tumor wide communication (mean)",x= "Cell type") +
  facet_wrap(~Treatmentlab,nrow=1)+
  theme(axis.title=element_blank(), axis.text.x = element_blank(), axis.text=element_blank(),strip.text = element_blank(),legend.title=element_blank(),legend.text=element_blank() ,legend.position = "none")+ 
  scale_y_continuous(breaks=log(c(0.001,0.01,0.1,1)),labels=c(0.001,0.01,0.1,1))
#ggsave(paste0(paperfile,"BLANK Ribo and Letrozole Cell type contribution to tumor wide communication.png"),height=10,width=10)

plt1[, Ligandiscancer:=F]
plt1[LigandPhenoCelltype=="Cancer cells", Ligandiscancer:=T]
summary(lm( log(exp(meanSig)-1)~ Ligandiscancer , plt1))

plt1[, Receptoriscancer:=F]
plt1[ReceptorPhenoCelltype=="Cancer cells", Receptoriscancer:=T]
summary(lm( log(exp(meanSig)-1)~ Receptoriscancer , plt1))


ggplot(plt1, aes(y= log(exp(meanSig)-1), x= ReceptorPhenoCelltype, col= ReceptorPhenoCelltype,fill=ReceptorPhenoCelltype,group=interaction(ReceptorPhenoCelltype)))+geom_violin(col=NA,scale="width",alpha=0.5)+ geom_point()+
  theme_classic(base_size=26)+
  theme(axis.text.x=element_text(angle=90,vjust=0.5)) +
  theme(aspect.ratio=1,legend.position="none") +labs(y="Strength of communicaiton \n received (mean)",x= "Cell type") +
  facet_wrap(~Treatmentlab)+ 
  scale_y_continuous(breaks=log(c(0.001,0.01,0.1,1)),labels=c(0.001,0.01,0.1,1))
#ggsave(filename = "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Strength communication cell types AllArms/Cell type strength of communicaiton received by treatment.png")
ggsave(paste0(paperfile,"Ribo and Letrozole Cell type strength of communicaiton received.png"),height=10,width=10)

ggplot(plt1, aes(y= log(exp(meanSig)-1), x= ReceptorPhenoCelltype, col= ReceptorPhenoCelltype,fill=ReceptorPhenoCelltype,group=interaction(ReceptorPhenoCelltype)))+geom_violin(col=NA,scale="width",alpha=0.5)+ geom_point()+theme_classic(base_size=16)+
  theme(axis.text.x=element_text(angle=90,vjust=0.5)) +
  theme(aspect.ratio=1,legend.position="none") +labs(y="Strength of communicaiton received (mean)",x= "Cell type") +
  facet_wrap(~Treatmentlab)+
  theme(axis.title=element_blank(), axis.text.x = element_blank(), axis.text=element_blank(),strip.text = element_blank(),legend.title=element_blank(),legend.text=element_blank() ,legend.position = "none")
#ggsave(paste0(paperfile,"BLANK Ribo and Letrozole Cell type strength of communicaiton received.png"),height=10,width=10)





sort(unique(CCI$R))
growthFactorReceptors2 <- c("EGFR","ERBB2","ERBB3","ERBB4")
growthFactorReceptors2 <- c("ERBB2","ERBB3","ERBB4")
CCI[,isGrowthFactor:=F]
CCI[R%in%growthFactorReceptors2, isGrowthFactor:=T]
# calculate the average strength of communication in the data
summaryCCI <- data.table( CCI[][order(Treat,LigandPhenoCelltype, ReceptorPhenoCelltype ) ] %>% 
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

ggplot(plt1[isGrowthFactor==T], aes(y= log(exp(meanSig)-1), x= LigandPhenoCelltype, col= LigandPhenoCelltype,fill=LigandPhenoCelltype,group=interaction(LigandPhenoCelltype)))+geom_violin(scale="width",alpha=0.5)+ geom_point()+theme_classic(base_size=16)+theme(axis.text.x=element_text(angle=90,vjust=0.5))+
  theme(aspect.ratio=1,legend.position="none") +labs(y="Contribution to tumor wide \n ERBB growth factor communication (mean)",x= "Cell type") +
  facet_wrap(~Treat)+  scale_y_continuous(breaks=log(c(0.001,0.01,0.1,1)),labels=c(0.001,0.01,0.1,1))

ggsave(filename = "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Strength communication cell types AllArms/ERBB family specific Cell type contribution to tumor wide communication by treatment.png")




ggplot(plt1[isGrowthFactor==T], aes(y= log(exp(meanSig)-1), x= ReceptorPhenoCelltype, col= ReceptorPhenoCelltype,fill=ReceptorPhenoCelltype,group=interaction(ReceptorPhenoCelltype)))+geom_violin(col=NA,scale="width",alpha=0.5)+ geom_point()+theme_classic(base_size=16)+
  theme(axis.text.x=element_text(angle=90,vjust=0.5)) +
  theme(aspect.ratio=1,legend.position="none") +labs(y="Strength of ERBB growth factor \n communication received (mean)",x= "Cell type") +
  facet_wrap(~Treat)+  scale_y_continuous(breaks=log(c(0.001,0.01,0.1,1)),labels=c(0.001,0.01,0.1,1))
ggsave(filename = "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Strength communication cell types AllArms/ERBB family specific Cell type strength of communicaiton received by treatment.png")






#macrophage colony-stimulating factor receptor
growthFactorReceptors2 <- c("CSF1R","CSF2RA")#,"FGFR1","FGFR2","FGFR3","FGFR4","FGFRL1" ,"TGFBR1","TGFBR2" ,"TGFBR3"   )
CCI[,isGrowthFactor:=F]
CCI[R%in%growthFactorReceptors2, isGrowthFactor:=T]
summaryCCI <- data.table( CCI[][order(Treat,LigandPhenoCelltype, ReceptorPhenoCelltype ) ] %>% 
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

ggplot(plt1[isGrowthFactor==T], aes(y= log(exp(meanSig)-1), x= LigandPhenoCelltype, col= LigandPhenoCelltype,fill=LigandPhenoCelltype,group=interaction(LigandPhenoCelltype)))+geom_violin(scale="width",alpha=0.5)+ geom_point()+theme_classic(base_size=16)+theme(axis.text.x=element_text(angle=90,vjust=0.5))+
  theme(aspect.ratio=1,legend.position="none") +labs(y="Contribution to tumor wide \n macrophage colony stimulating factor communication (mean)",x= "Cell type") +
  facet_wrap(~Treat)+
  scale_y_continuous(breaks=log(c(0.001,0.01,0.1,1)),labels=c(0.001,0.01,0.1,1))
ggsave(filename = "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Strength communication cell types AllArms/Macrophage CSF specific Cell type contribution to tumor wide communication by treatment.png")


ggplot(plt1[isGrowthFactor==T], aes(y= log(exp(meanSig)-1), x= ReceptorPhenoCelltype, col= ReceptorPhenoCelltype,fill=ReceptorPhenoCelltype,group=interaction(ReceptorPhenoCelltype)))+geom_violin(col=NA,scale="width",alpha=0.5)+ geom_point()+theme_classic(base_size=16)+
  theme(axis.text.x=element_text(angle=90,vjust=0.5)) +
  theme(aspect.ratio=1,legend.position="none") +labs(y="Strength of macrophage colony stimulating factor \n communication received (mean)",x= "Cell type") +
  facet_wrap(~Treat)+
  scale_y_continuous(breaks=log(c(0.001,0.01,0.1,1)),labels=c(0.001,0.01,0.1,1))
ggsave(filename = "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Strength communication cell types AllArms/Macrophage CSF specific Cell type strength of communicaiton received by treatment.png")




###CCR5
growthFactorReceptors2 <- c("CCR5") 
CCI[,isGrowthFactor:=F]
CCI[R%in%growthFactorReceptors2, isGrowthFactor:=T]
summaryCCI <- data.table( CCI[][order(Treat,LigandPhenoCelltype, ReceptorPhenoCelltype ) ] %>% 
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

ggplot(plt1[isGrowthFactor==T], aes(y= log(exp(meanSig)-1), x= LigandPhenoCelltype, col= LigandPhenoCelltype,fill=LigandPhenoCelltype,group=interaction(LigandPhenoCelltype)))+geom_violin(scale="width",alpha=0.5)+ geom_point()+theme_classic(base_size=16)+theme(axis.text.x=element_text(angle=90,vjust=0.5))+
  theme(aspect.ratio=1,legend.position="none") +labs(y="Contribution to tumor wide \n CCR5 communication (mean)",x= "Cell type") +
  facet_wrap(~Treat)+
  scale_y_continuous(breaks=log(c(0.001,0.01,0.1,1)),labels=c(0.001,0.01,0.1,1))

ggsave(filename = "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Strength communication cell types AllArms/CCR5 specific Cell type contribution to tumor wide communication by treatment.png")


ggplot(plt1[isGrowthFactor==T], aes(y= log(exp(meanSig)-1), x= ReceptorPhenoCelltype, col= ReceptorPhenoCelltype,fill=ReceptorPhenoCelltype,group=interaction(ReceptorPhenoCelltype)))+geom_violin(col=NA,scale="width",alpha=0.5)+ geom_point()+theme_classic(base_size=16)+
  theme(axis.text.x=element_text(angle=90,vjust=0.5)) +
  theme(aspect.ratio=1,legend.position="none") +labs(y="Strength of CCR5 \n communication received (mean)",x= "Cell type") +
  facet_wrap(~Treat)+
  scale_y_continuous(breaks=log(c(0.001,0.01,0.1,1)),labels=c(0.001,0.01,0.1,1))

ggsave(filename = "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Strength communication cell types AllArms/CCR5 specific Cell type strength of communicaiton received by treatment.png")




## VEGF receptor : FLT1
growthFactorReceptors2 <- c("FLT1")#,"FGFR1","FGFR2","FGFR3","FGFR4","FGFRL1" ,"TGFBR1","TGFBR2" ,"TGFBR3"   )

CCI[,isGrowthFactor:=F]
CCI[R%in%growthFactorReceptors2, isGrowthFactor:=T]
# calculate the average strength of communication in the data
summaryCCI <- data.table( CCI[][order(Treat,LigandPhenoCelltype, ReceptorPhenoCelltype ) ] %>% 
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

ggplot(plt1[isGrowthFactor==T], aes(y= log(exp(meanSig)-1), x= LigandPhenoCelltype, col= LigandPhenoCelltype,fill=LigandPhenoCelltype,group=interaction(LigandPhenoCelltype)))+geom_violin(scale="width",alpha=0.5)+ geom_point()+theme_classic(base_size=16)+theme(axis.text.x=element_text(angle=90,vjust=0.5))+
  theme(aspect.ratio=1,legend.position="none") +labs(y="Contribution to tumor wide \n VEGF growth factor communication (mean)",x= "Cell type") +
  facet_wrap(~Treat)+
  scale_y_continuous(breaks=log(c(0.001,0.01,0.1,1)),labels=c(0.001,0.01,0.1,1))

ggsave(filename = "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Strength communication cell types AllArms/VEGF specific Cell type contribution to tumor wide communication by treatment.png")


ggplot(plt1[isGrowthFactor==T], aes(y= log(exp(meanSig)-1), x= ReceptorPhenoCelltype, col= ReceptorPhenoCelltype,fill=ReceptorPhenoCelltype,group=interaction(ReceptorPhenoCelltype)))+geom_violin(col=NA,scale="width",alpha=0.5)+ geom_point()+theme_classic(base_size=16)+
  theme(axis.text.x=element_text(angle=90,vjust=0.5)) +
  theme(aspect.ratio=1,legend.position="none") +labs(y="Strength of VEGF \n communication received (mean)",x= "Cell type") +
  facet_wrap(~Treat)+
  scale_y_continuous(breaks=log(c(0.001,0.01,0.1,1)),labels=c(0.001,0.01,0.1,1))
ggsave(filename = "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Strength communication cell types AllArms/VEGF specific Cell type strength of communicaiton received by treatment.png")

## FGFR1
growthFactorReceptors2 <- c("FGFR")
CCI[,isGrowthFactor:=F]
CCI[grep("FGFR1",R), isGrowthFactor:=T]
# calculate the average strength of communication in the data
summaryCCI <- data.table( CCI[][order(Treat,LigandPhenoCelltype, ReceptorPhenoCelltype ) ] %>% 
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

ggplot(plt1[isGrowthFactor==T], aes(y= log(exp(meanSig)-1), x= LigandPhenoCelltype, col= LigandPhenoCelltype,fill=LigandPhenoCelltype,group=interaction(LigandPhenoCelltype)))+geom_violin(scale="width",alpha=0.5)+ geom_point()+theme_classic(base_size=16)+theme(axis.text.x=element_text(angle=90,vjust=0.5))+
  theme(aspect.ratio=1,legend.position="none") +labs(y="Contribution to tumor wide \n FGFR1 growth factor communication (mean)",x= "Cell type") +
  facet_wrap(~Treat)+
  scale_y_continuous(breaks=log(c(0.001,0.01,0.1,1)),labels=c(0.001,0.01,0.1,1))
ggsave(filename = "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Strength communication cell types AllArms/FGFR1 specific Cell type contribution to tumor wide communication by treatment.png")


ggplot(plt1[isGrowthFactor==T], aes(y= log(exp(meanSig)-1), x= ReceptorPhenoCelltype, col= ReceptorPhenoCelltype,fill=ReceptorPhenoCelltype,group=interaction(ReceptorPhenoCelltype)))+geom_violin(col=NA,scale="width",alpha=0.5)+ geom_point()+theme_classic(base_size=16)+
  theme(axis.text.x=element_text(angle=90,vjust=0.5)) +
  theme(aspect.ratio=1,legend.position="none") +labs(y="Strength of FGFR1 \n communication received (mean)",x= "Cell type") +
  facet_wrap(~Treat)+
  scale_y_continuous(breaks=log(c(0.001,0.01,0.1,1)),labels=c(0.001,0.01,0.1,1))
ggsave(filename = "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Strength communication cell types AllArms/FGFR1 specific Cell type strength of communicaiton received by treatment.png")
