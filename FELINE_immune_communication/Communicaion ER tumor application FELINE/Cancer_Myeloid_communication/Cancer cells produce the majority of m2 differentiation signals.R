rm(list=ls())   
require(rTensor); require(nnTensor); library(abind)
require(data.table); require(dplyr); require(ggplot2); require(tidyr)
require(igraph)

### Phenotype data
load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/PhenotypesOfAllCellTypesAllArms/PhenotypesOfAllCellTypesAllArms.RData") #allphenotypes,UMAPlocs ,UMAPfiles,umapDImRedloc,umapDImRedfiles,nCellTypes,
# Encode the subtypes that have been analysed using UMAP
tmp <- data.table(allphenotypes %>% group_by(key_) %>% slice(1)) %>% dplyr::select(Celltype , Celltype_subtype) %>% unique
tmp[ , PhenoCelltype:= Celltype]
tmp[Celltype_subtype %in% c("CD4+ T cells", "Tregs"), PhenoCelltype:= "CD4+ T cells"]
tmp[Celltype_subtype %in% c("CD8+ T cells", "NK cells"), PhenoCelltype:= "CD8+ T cells"]
tmp[Celltype_subtype %in% c("Macrophages", "DC", "Monocytes"), PhenoCelltype:= "Macrophages"]
tmp[Celltype_subtype %in% c("Vas-Endo", "Lym-Endo","Endothelial cells"), PhenoCelltype:= "Endothelial cells"]
tmp[Celltype_subtype %in% c("Cancer cells"), PhenoCelltype:= "Cancer cells"]
tmp[Celltype_subtype %in% c("Normal epithelial cells"), PhenoCelltype:= "Normal epithelial cells"]
allphenotypes <- merge(allphenotypes %>% dplyr::select(-c( paste0("V", 1:5), "nCount_RNA", "nFeature_RNA", "Infercnv_CNA", "Platform","Timepoint","Sample_p_t","prop_change", "file_string", "day_fact" , "Ribo","TreatLab", "Burden_t0" ,"BurdenTracked" ,"Day0", "DayLastBurd", "Dose_lab","TreatCode","TreatCodeOrd","dynamic_class2","rgr_A","rgr_B")), tmp, by= c("Celltype", "Celltype_subtype") )
#allphenotypes[Patient.Study.ID=="001-138"&Day==14][Celltype=="Macrophages"]

### Unique clusters of cells (ALL) and their umap discretization level
uu <- unique( allphenotypes %>% dplyr::select( c("Celltype","PhenoCelltype", "key_", paste0("Disc_V", 1:5) ) ) )  #"Celltype_subtype",

### Load Ligand Receptor database list of Ramilowski et al 2015
load( "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/LigandReceptor/Filtered_Human-2015-Ramilowski-LR-pairs.RData")
LRgenelist <- unique( c(LRpairsFiltered$HPMR.Receptor, LRpairsFiltered$HPMR.Ligand) )

### Load some signalling data for downstream analysis
savelocCCI <- "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/CommunicationOutputAllArms/"

# Are we studying the per indiv or per cell type level of signalling?
perIndiv=FALSE

# Names of LR data to analyze
filenamesCCI <- list.files(savelocCCI)#%>%length

# First let's do this for the population level
#summarytable <- read.csv("/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/TensorAnalysisOutput/Ligand Receptor list fro NTD model of population cellcell communication_AB.csv")
load(file ="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Communication output merged ALL ARMS/PopulationCommunicationMerged ALL ARMS.RData" )#CCI,allphenotypes, uu,perIndiv,
CCI[,Treat:="CombinationRibo"]
CCI[ARM=="A",Treat:="LetrozoleAlone"]
CCI[,scaleTransduction:=scale(TransductionMu,center=F),by=c("Pair.Name")] #"ReceptorPhenoCelltype",
CCI[!is.finite(TransductionMu),scaleTransduction:=0]
CCI[,scalelnTransductionMu:=scale(log(1+TransductionMu)),by=Pair.Name]

extended.signiftable<- data.table(read.csv( file="~/Dropbox/Cancer_pheno_evo/data/FELINE2/Final communication files/Cancer to macropage communication supervised analysis/Version2 Cancer to macropage M2 stimulating communication supervised analysis.csv"))
#extended.signiftable<- data.table(read.csv( file="~/Dropbox/Cancer_pheno_evo/data/FELINE2/Final communication files/Ribo Extended Communication differences between resistant and sensitive tumors at day0.csv"))
extendelistC1<-as.character(extended.signiftable$Pair.Name)

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


0.66328706
avmacrosignals2[Day==0][Treatlabel=="Combination ribociclib"][Tumorresponse=="Resistant"][broadCat=="Cancer cells"]$m2
avmacrosignals2[Day==0][Treatlabel=="Combination ribociclib"][Tumorresponse!="Resistant"]$m2%>%sum()

ggplot(macrosignals,
       aes(x=broadCat, y=scalelnTransductionMu, col=dynamic_class3, fill= dynamic_class3))+
  geom_violin(scale="width", alpha=0.5) + geom_point(position=position_dodge(width=0.9))+
  theme_classic()+theme(aspect.ratio=1, axis.text.x = element_text(angle=90))+
  geom_errorbar(data=avmacrosignals,aes(y=1,ymax=ucl,ymin=lcl),col="black",width=0.4,position=position_dodge(width=0.9))+
  geom_point(data=avmacrosignals,aes(y=m),col="black",width=0.4,position=position_dodge(width=0.9)) +
  facet_wrap(~Treat)

macrosignals[,is.Cancer:=F]
macrosignals[broadCat=="Cancer cells",is.Cancer:=T]
ggplot(macrosignals,
       aes(x=Pair.Name, y=scalelnTransductionMu, col=is.Cancer, fill= is.Cancer))+
  geom_violin(scale="width", alpha=0.5) + geom_point(position=position_dodge(width=0.9))+
  theme_classic()+theme(aspect.ratio=1, axis.text.x = element_text(angle=90))+
  geom_errorbar(data=avmacrosignals,aes(y=1,ymax=ucl,ymin=lcl),col="black",width=0.4,position=position_dodge(width=0.9))+
  geom_point(data=avmacrosignals,aes(y=m),col="black",width=0.4,position=position_dodge(width=0.9))



pp<- data.frame( corruptsigs[Day!=180] %>% select(-c(Patient.Study.ID,dynamic_class3,Day))  )
rownames(pp)<- paste0("samp",1:nrow(pp))
rownam <- data.frame( dynamic_class3= corruptsigs[Day!=180] $dynamic_class3 )
rownames(rownam)<- paste0("samp",1:nrow(pp))

pheatmap::pheatmap(  data.frame(pp),annotation_row=rownam)

pp<-data.frame(corruptsigs[Day!=180][]%>%dplyr::select(-c(dynamic_class3,Day, Patient.Study.ID)))
annotation_row <-data.frame(corruptsigs[Day!=180]%>%dplyr::select(dynamic_class3, Day,Patient.Study.ID)) #orig.ident
rownames(annotation_row) <- paste0(annotation_row$Patient.Study.ID, annotation_row$Day);annotation_row$Patient.Study.ID <-NULL;annotation_row$Day <-NULL
names(annotation_row) <- "Tumor response"
require(ggsci)

annotation_col = list(
  "Tumor response" = rev(pal_npg("nrc")(2) ))
names(annotation_col[[1]] ) <- c("Response","Non-response")
rownames(pp) <-rownames(annotation_row)
pheatmap::pheatmap((pp),cluster_rows=F,annotation_row = annotation_row, annotation_colors = annotation_col,scale="none",show_rownames=F,annotation_legend=F,border_color=NA)
pheatmap::pheatmap((pp),cluster_rows=F,annotation_row = annotation_row, annotation_colors = annotation_col,scale="column",show_rownames=F,annotation_legend=F,border_color=NA)
#pheatmap::pheatmap(t(pp),cluster_cols=F,annotation_col = annotation_row, annotation_colors = annotation_col,scale="none",show_colnames=F,annotation_legend=F,border_color=NA)


pheatmap::pheatmap(sqrt( corruptsigs%>%dplyr::select(-c(dynamic_class3,Day, Patient.Study.ID)) ),annotation_row = annotation_row)







