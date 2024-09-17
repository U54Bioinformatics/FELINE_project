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

# Select cancer - macrophage communications for communication pathways linked to M2 differentiation.
corruptsigs0 <- CCI[Day!=180][Treat=="CombinationRibo"][LigandPhenoCelltype=="Cancer cells"][ReceptorPhenoCelltype=="Macrophages"]
nsm1<- unique(corruptsigs0[Day!=180]$Pair.Name)
extractThese0 <- c(grep("FGF2_",nsm1),#grep("FGF1_",nsm1),
                   grep("IL11_",nsm1),grep("IL12A_",nsm1),
  grep("CLCF1_",nsm1),grep("MADCAM1_",nsm1),
  #grep("SEMA3",nsm1),grep("SEMA4",nsm1),
  grep("PLAU_",nsm1),grep("GDF9_",nsm1))
extractThese1 <- c(grep("_IL4",nsm1),grep("_IL13",nsm1),grep("_IL10",nsm1)         )
extractThese2 <- c(grep("_CSF1R",nsm1),grep("_CSF2RB",nsm1))#,grep("_MERTK",nsm1),grep("GAS_",nsm1) )
extractThese3 <- c(grep("_FGFR1",nsm1),grep("IL6_",nsm1) )
selectList<- nsm1[c(extractThese0,extractThese1,extractThese2,extractThese3)]
selectList<-selectList[selectList!="NCAM1_FGFR1"]
Supervisedcorruptsigs <- corruptsigs0[Pair.Name%in%selectList] 
#Supervisedcorruptsigs <- corruptsigs[Pair.Name%in%nsm1[extractThese]] 

corruptsigstmp <- data.table(Supervisedcorruptsigs%>%dplyr::select(Pair.Name, scalelnTransductionMu, Patient.Study.ID, dynamic_class3, Day)%>% spread(Pair.Name,scalelnTransductionMu,fill=0))[order(dynamic_class3,Patient.Study.ID,Day)]

clustmatResponder<- data.frame( corruptsigstmp[][dynamic_class3=="Response"]%>% dplyr::select(-c(Patient.Study.ID,dynamic_class3,Day))  )
clustmatNonResponder<- data.frame( corruptsigstmp[][dynamic_class3=="Non-response"]%>% dplyr::select(-c(Patient.Study.ID,dynamic_class3,Day))  )
TumorOrder <- c(
unique(corruptsigstmp[][dynamic_class3=="Non-response"]$Patient.Study.ID[hclust(dist(clustmatNonResponder))$order]),
unique(corruptsigstmp[][dynamic_class3=="Response"]$Patient.Study.ID   [hclust(dist(clustmatResponder))$order]))

Supervisedcorruptsigs$Patient.Study.ID <- factor(Supervisedcorruptsigs$Patient.Study.ID , levels=TumorOrder)
  
corruptsigs <- data.table(Supervisedcorruptsigs%>%dplyr::select(Pair.Name, scalelnTransductionMu, Patient.Study.ID, dynamic_class3, Day)%>% spread(Pair.Name,scalelnTransductionMu,fill=NA))[order(dynamic_class3,Patient.Study.ID,Day)]

# Plot as heatmap
pp<- data.frame( corruptsigs%>% dplyr::select(-c(Patient.Study.ID,dynamic_class3,Day))  )
pp<- pp[,colMeans(!is.na(pp))>0.3]

rownames(pp)<- paste0("samp",1:nrow(pp))
rownam <- data.frame( dynamic_class3= corruptsigs$dynamic_class3 )
rownames(rownam)<- paste0("samp",1:nrow(pp))
wellMeasured <- rowMeans(!is.na(pp))>0.33

summary(lm(CSF1_CSF1R~dynamic_class3,data= corruptsigs[Day!=180][wellMeasured,][]))
ggplot(corruptsigs[Day!=180][wellMeasured,][],aes(y=CSF1_CSF1R,x=dynamic_class3,col=Day) ) +geom_point()


pp<- pp[wellMeasured,]
rownam <- rownam[wellMeasured,]

annotation_row <- data.frame(corruptsigs[Day!=180][wellMeasured,]%>%dplyr::select(dynamic_class3, Day,Patient.Study.ID)) #orig.ident
rownames(annotation_row) <- paste0(annotation_row$Patient.Study.ID, annotation_row$Day);annotation_row$Patient.Study.ID <-NULL;annotation_row$Day <-NULL
names(annotation_row) <- "Tumor response"
require(ggsci)

annotation_col = list( "Tumor response" = rev(pal_npg("nrc")(2) ))
names(annotation_col[[1]] ) <- c("Response","Non-response")
rownames(pp) <-rownames(annotation_row)
brk <- sum(annotation_row$'Tumor response'=="Non-response")

fillseq<- seq(-2.35,4.45,length=101)
#pheatmap::pheatmap((pp),cluster_rows=F,annotation_row = annotation_row, annotation_colors = annotation_col,scale="none",show_rownames=F,annotation_legend=F,border_color=NA)
pheatmap::pheatmap((pp), gaps_row=brk, cluster_rows=F, cluster_cols=F, breaks=fillseq, na_col="white", annotation_row = annotation_row, annotation_colors = annotation_col, scale="none",show_rownames=F,annotation_legend=F,border_color=NA)
#pheatmap::pheatmap((pp),gaps_row=brk,cluster_rows=F,cluster_cols=F,na_col="white",annotation_row = annotation_row, annotation_colors = annotation_col,scale="column",show_rownames=F,annotation_legend=F,border_color=NA)
#save 5x 6 plot (file="~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/cancer to macrophage communications/Version2b Letrozole supervised list of cancer signals to macrophages.pdf")

blankrowannot <- annotation_row
names( blankrowannot)<- c(" ")

names( annotation_col)<- names( blankrowannot)
names(annotation_col[[1]] ) <- c("Response","Non-response")

p1 <- pheatmap::pheatmap((pp), gaps_row=brk, cluster_rows=F, cluster_cols=F, breaks=fillseq, na_col="white", 
                   annotation_row = blankrowannot, annotation_colors = annotation_col, 
                   scale= "none", show_rownames=F, annotation_legend=F, border_color=NA,labels_col=rep("",100),
                   labels_row= rep("",100) ,
                   legend_breaks = -4:2,
                   legend_labels = rep("",7) )

paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"

#save_pheatmap(p1,filename= paste0(paperfile,"BLANK Ribo and Letrozole Ribo heatmap of cancer signals to macrophages.pdf"),height=480,width=480)

#BLANK Ribo and Letrozole RIBO heatmap of cancer signals to macrophages.pdf




# look up p values of the difference in communication between resistant and sensitive tumors
load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/When cells communicate_communication changes AllArms/trendsby response AllArms.RData")
assessmentSupervised <- assessment[Pair.Name%in%names(pp)][Treat=="CombinationRibo"] [Model=="one way anova_DayStart_Response"][LigandPhenoCelltype=="Cancer cells"][ReceptorPhenoCelltype=="Macrophages"][coefficient!="(Intercept)"][order(Pair.Name)]
assessmentSupervised$adjust.p.val <- NULL
assessmentSupervised$i <- NULL
assessmentSupervised$adj.r.squared <- NULL
#assessmentSupervised[Estimate>0]
write.csv(assessmentSupervised,file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Final communication files/Cancer to macropage communication supervised analysis/Version2 Cancer to macropage M2 stimulating communication supervised analysis.csv")

assessmentSupervised <- assessment[Pair.Name%in%names(pp)][Treat!="CombinationRibo"] [Model=="one way anova_DayStart_Response"][LigandPhenoCelltype=="Cancer cells"][ReceptorPhenoCelltype=="Macrophages"][coefficient!="(Intercept)"][order(Pair.Name)]
assessmentSupervised$adjust.p.val <- NULL
assessmentSupervised$i <- NULL
assessmentSupervised$adj.r.squared <- NULL
write.csv(assessmentSupervised,file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Final communication files/Cancer to macropage communication supervised analysis/Version2Letrozole Cancer to macropage M2 stimulating communication supervised analysis.csv")
#selectedLR<- read.csv(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Final communication files/Cancer to macropage communication supervised analysis/Version2Letrozole Cancer to macropage M2 stimulating communication supervised analysis.csv") %>%select(Pair.Name, LigandPhenoCelltype, ReceptorPhenoCelltype)
write.csv(assessmentSupervised,file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Final communication files/Cancer to macropage communication supervised analysis/Version2Letrozole Cancer to macropage M2 stimulating communication supervised analysis.csv")

# Select cancer - macrophage communications for communication pathways linked to M2 differentiation.
corruptsigs<-CCI[Day!=180][Treat!="CombinationRibo"][LigandPhenoCelltype=="Cancer cells"][ReceptorPhenoCelltype=="Macrophages"]
nsm1<- unique(corruptsigs[Day!=180]$Pair.Name)
Supervisedcorruptsigs <- corruptsigs[Pair.Name%in%assessmentSupervised$Pair.Name] 
corruptsigs <- data.table(Supervisedcorruptsigs%>%dplyr::select(Pair.Name, scalelnTransductionMu, Patient.Study.ID,dynamic_class3,Day)%>% spread(Pair.Name,scalelnTransductionMu,fill=NA))[order(dynamic_class3,Patient.Study.ID,Day)]

# Plot as heatmap
pp <- data.frame( corruptsigs%>% dplyr::select(-c(Patient.Study.ID,dynamic_class3,Day))  )
#pp<- pp[,colMeans(!is.na(pp))>0.3]
rownames(pp)<- paste0("samp",1:nrow(pp))
rownam <- data.frame( dynamic_class3= corruptsigs$dynamic_class3 )
rownames(rownam)<- paste0("samp",1:nrow(pp))

annotation_row <-data.frame(corruptsigs[Day!=180]%>%dplyr::select(dynamic_class3, Day,Patient.Study.ID)) #orig.ident
rownames(annotation_row) <- paste0(annotation_row$Patient.Study.ID, annotation_row$Day);annotation_row$Patient.Study.ID <-NULL;annotation_row$Day <-NULL
names(annotation_row) <- "Tumor response"
require(ggsci)

annotation_col = list( "Tumor response" = rev(pal_npg("nrc")(2) ))
names(annotation_col[[1]] ) <- c("Response","Non-response")
rownames(pp) <-rownames(annotation_row)
brk <- sum(annotation_row$'Tumor response'=="Non-response")
#pheatmap::pheatmap((pp),cluster_rows=F,annotation_row = annotation_row, annotation_colors = annotation_col,scale="none",show_rownames=F,annotation_legend=F,border_color=NA)
pheatmap::pheatmap((pp),gaps_row=brk,breaks= fillseq,cluster_rows=F,cluster_cols=F,na_col="white",annotation_row = annotation_row, annotation_colors = annotation_col,scale="none",show_rownames=F,annotation_legend=F,border_color=NA)
#save 5 x 6 plot (file="~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/cancer to macrophage communications/Version2b Letrozole supervised list of cancer signals to macrophages.png.pdf")

#/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/Ribo and Letrozole LETROZOLE heatmap of cancer signals to macrophages
blankrowannot <- annotation_row
names( blankrowannot)<- c(" ")

names( annotation_col)<- names( blankrowannot)
names(annotation_col[[1]] ) <- c("Response","Non-response")

p1 <- pheatmap::pheatmap((pp), gaps_row=brk, cluster_rows=F, cluster_cols=F, breaks=fillseq, na_col="white", 
                         annotation_row = blankrowannot, annotation_colors = annotation_col, 
                         scale= "none", show_rownames=F, annotation_legend=F, border_color=NA,labels_col=rep("",100),
                         labels_row= rep("",100) ,
                         legend_breaks = -4:2,
                         legend_labels = rep("",7) )

paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"

#save_pheatmap(p1,filename= paste0(paperfile,"BLANK Ribo and Letrozole LETROZOLE heatmap of cancer signals to macrophages.pdf"),height=480,width=480)

# BLANK Ribo and Letrozole LETROZOLE heatmap of cancer signals to macrophages.pdf





# extended.signiftable<- data.table(read.csv( file="~/Dropbox/Cancer_pheno_evo/data/FELINE2/Final communication files/Ribo Extended Communication differences between resistant and sensitive tumors at day0.csv"))
# extendelistC1<-as.character(extended.signiftable$Pair.Name)
# 
# corruptsigslong<-CCI[Pair.Name%in%extendelistC1][Treat=="CombinationRibo"][LigandPhenoCelltype=="Cancer cells"][ReceptorPhenoCelltype=="Macrophages"]
# corruptsigs <- data.table(corruptsigslong%>%select(Pair.Name, scalelnTransductionMu, Patient.Study.ID,dynamic_class3,Day)%>% spread(Pair.Name,scalelnTransductionMu,fill=0))[order(dynamic_class3,)]
# 
# ggplot(gather(corruptsigs,Pair.Name,scalelnTransductionMu, ADIPOQ_ADIPOR1:THBS2_ITGA6 ),
#        aes(y=Pair.Name,x=Patient.Study.ID,fill=scalelnTransductionMu ))+geom_tile()+
#   facet_wrap(~dynamic_class3,scales="free_x")
# 
# pp<- data.frame( corruptsigs[Day!=180] %>% select(-c(Patient.Study.ID,dynamic_class3,Day))  )
# rownames(pp)<- paste0("samp",1:nrow(pp))
# rownam <- data.frame( dynamic_class3= corruptsigs[Day!=180] $dynamic_class3 )
# rownames(rownam)<- paste0("samp",1:nrow(pp))
# 
# pheatmap::pheatmap(  data.frame(pp),annotation_row=rownam)
# 
# pp <- data.frame(corruptsigs[Day!=180][]%>%dplyr::select(-c(dynamic_class3,Day, Patient.Study.ID)))
# annotation_row <-data.frame(corruptsigs[Day!=180]%>%dplyr::select(dynamic_class3, Day,Patient.Study.ID)) #orig.ident
# rownames(annotation_row) <- paste0(annotation_row$Patient.Study.ID, annotation_row$Day);annotation_row$Patient.Study.ID <-NULL;annotation_row$Day <-NULL
# names(annotation_row) <- "Tumor response"
# require(ggsci)
# 
# annotation_col = list(
#   "Tumor response" = rev(pal_npg("nrc")(2) ))
# names(annotation_col[[1]] ) <- c("Response","Non-response")
# rownames(pp) <-rownames(annotation_row)
# pheatmap::pheatmap((pp),cluster_rows=F,annotation_row = annotation_row, annotation_colors = annotation_col,scale="none",show_rownames=F,annotation_legend=F,border_color=NA)
# pheatmap::pheatmap((pp),cluster_rows=F,annotation_row = annotation_row, annotation_colors = annotation_col,scale="column",show_rownames=F,annotation_legend=F,border_color=NA)
# #pheatmap::pheatmap(t(pp),cluster_cols=F,annotation_col = annotation_row, annotation_colors = annotation_col,scale="none",show_colnames=F,annotation_legend=F,border_color=NA)
# 
# 
# pheatmap::pheatmap(sqrt( corruptsigs%>%dplyr::select(-c(dynamic_class3,Day, Patient.Study.ID)) ),annotation_row = annotation_row)
# 
# 
# 
# 
# 
# 
# 
# 
# 

