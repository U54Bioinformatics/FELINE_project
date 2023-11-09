rm(list=ls())   
require(rTensor); require(nnTensor); library(abind)
require(data.table); require(dplyr); require(ggplot2); require(tidyr)
require(igraph)
require(ggsci)
### Phenotype data
load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/PhenotypesOfAllCellTypesAllArmsCohort2/PhenotypesOfAllCellTypesAllArmsCohort2.RData")
# Encode the subtypes that have been analysed using UMAP
tmp <- data.table(allphenotypes %>% group_by(key_) %>% slice(1)) %>% dplyr::select(Celltype , Celltype_subtype) %>% unique
tmp[ , PhenoCelltype:= Celltype]
tmp[Celltype_subtype %in% c("CD4+ T cells", "Tregs"), PhenoCelltype:= "CD4+ T cells"]
tmp[Celltype_subtype %in% c("CD8+ T cells", "NK cells"), PhenoCelltype:= "CD8+ T cells"]
tmp[Celltype_subtype %in% c("Macrophages", "DC", "Monocytes"), PhenoCelltype:= "Macrophages"]
tmp[Celltype_subtype %in% c("Vas-Endo", "Lym-Endo","Endothelial cells"), PhenoCelltype:= "Endothelial cells"]
tmp[Celltype_subtype %in% c("Cancer cells"), PhenoCelltype:= "Cancer cells"]
tmp[Celltype_subtype %in% c("Normal epithelial cells"), PhenoCelltype:= "Normal epithelial cells"]
allphenotypes <- merge(allphenotypes %>%dplyr::select(-c( paste0("V", 1:5), "nCount_RNA", "nFeature_RNA",  "Platform","Sample_p_t", "file_string", "day_fact")), tmp, by= c("Celltype", "Celltype_subtype") )
allphenotypes[Patient.Study.ID=="001-171_853686",Patient.Study.ID:="001-171"]
load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE/Feline clinical input.RData")
allphenotypes <- merge(allphenotypes, unique(Clin_resp_dd_classAdd%>%dplyr::select(Patient.Study.ID,prop_change,dynamic_class)), by="Patient.Study.ID")

### Unique clusters of cells (ALL) and their umap discretization level
uu <- unique( allphenotypes %>% dplyr::select( c("Celltype","PhenoCelltype", "key_", paste0("Disc_V", 1:5) ) ) )  #"Celltype_subtype",

### Load Ligand Receptor database list of Ramilowski et al 2015
load( "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/LigandReceptor/Filtered_Human-2015-Ramilowski-LR-pairs.RData")
LRgenelist <- unique( c(LRpairsFiltered$HPMR.Receptor, LRpairsFiltered$HPMR.Ligand) )

### Load some signalling data for downstream analysis
savelocCCI <- "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/CommunicationOutputAllArmsCohort2/"

# Are we studying the per indiv or per cell type level of signalling?
perIndiv=FALSE

# Names of LR data to analyze
filenamesCCI <- list.files(savelocCCI)#%>%length

# First let's do this for the population level
#summarytable <- read.csv("/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/TensorAnalysisOutput/Ligand Receptor list fro NTD model of population cellcell communication_AB.csv")
load(file ="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/cohort2 Communication output merged/PopulationCommunicationMerged ALL ARMS cohort2.RData" )#CCI,allphenotypes, uu,perIndiv,
CCI[,Treat:="CombinationRibo"]
CCI[ARM=="A",Treat:="LetrozoleAlone"]

CCI[,scaleTransduction:=scale(TransductionMu,center=F),by=c("Pair.Name")] #"ReceptorPhenoCelltype",
CCI[!is.finite(TransductionMu),scaleTransduction:=0]
CCI[Patient.Study.ID=="001-171_853686",Patient.Study.ID:="001-171"]
CCI2<- CCI <-merge(CCI, unique(Clin_resp_dd_classAdd%>%dplyr::select(Patient.Study.ID,prop_change,dynamic_class)), by="Patient.Study.ID")
CCI[ prop_change < (2/3) ,dynamic_class3:="Response"]
CCI2[ prop_change < (2/3) ,dynamic_class3:="Response"]

CCI[,scalelnTransductionMu:=scale(log(1+TransductionMu)),by=Pair.Name]


# Select cancer - macrophage communications for communication pathways linked to M2 differentiation.
corruptsigs<-CCI[][Day!=180][Treat=="CombinationRibo"][LigandPhenoCelltype=="Cancer cells"][ReceptorPhenoCelltype=="Macrophages"]
nsm1<- unique(corruptsigs[Day!=180]$Pair.Name)
extractThese <- c(grep("CSF1_CSF1R",nsm1) ,grep("FGF2_",nsm1),grep("FGF1_",nsm1),grep("IL11_",nsm1),grep("IL12A_",nsm1),
                  grep("CLCF1_",nsm1),grep("MADCAM1_",nsm1),grep("SEMA3",nsm1),grep("SEMA4",nsm1),grep("PLAU_",nsm1),grep("GDF9_",nsm1))

extractThese <- which(nsm1 %in% c("CLCF1_CRLF1","CLCF1_IL6ST","CLCF1_LIFR", "CSF1_CSF1R", "FGF1_FGFR1","FGF17_FGFR1","FGF18_FGFR1","FGF2_CD44","FGF2_FGFR1","FGF2_NRP1","FGF2_SDC1","FGF2_SDC2","FGF2_SDC3","FGF2_SDC4","FGF23_FGFR1","FGF7_FGFR1",
                  "GDF9_ACVR2A", "GDF9_BMPR1A" ,"GDF9_BMPR1B", "GDF9_BMPR2" , "GDF9_TGFBR1","IL11_IL6ST","IL11_IL11RA","IL12A_IL12RB1","IL12A_IL12RB2","IL5_CSF2RB","MADCAM1_CD44", "MADCAM1_ITGA4",
                  "PLAU_IGF2R","PLAU_LRP1" ,"PLAU_LRP2" ,"PLAU_PLAUR") )
#nsm1[grep("FGF17",nsm1)]
  
Supervisedcorruptsigs <- corruptsigs[Pair.Name%in%nsm1[extractThese]] 
corruptsigs <- data.table(Supervisedcorruptsigs%>%dplyr::select(Pair.Name, scalelnTransductionMu, Patient.Study.ID,dynamic_class3,Day,prop_change)%>% spread(Pair.Name,scalelnTransductionMu,fill=NA))[order(dynamic_class3,-prop_change,Patient.Study.ID,Day)]

# Plot as heatmap
pp <- data.frame( corruptsigs%>% dplyr::select(-c(Patient.Study.ID,dynamic_class3,Day))  )
pp<- pp[,colMeans(!is.na(pp))>0.3]
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
#pheatmap::pheatmap((pp),cluster_rows=F,annotation_row = annotation_row, annotation_colors = annotation_col,scale="none",show_rownames=F,annotation_legend=F,border_color=NA)
pheatmap::pheatmap((pp),cluster_rows=F,cluster_cols=F,na_col="white",annotation_row = annotation_row, annotation_colors = annotation_col,scale="none",show_rownames=F,annotation_legend=F,border_color=NA)
pheatmap::pheatmap((pp),cluster_rows=F,cluster_cols=F,na_col="white",annotation_row = annotation_row, annotation_colors = annotation_col,scale="column",show_rownames=F,annotation_legend=F,border_color=NA)

summary(lm(CSF1_CSF1R~dynamic_class3,data= corruptsigs[Day==0]))
ggplot(corruptsigs[Day==0],aes(y=CSF1_CSF1R,x=dynamic_class3,col=Day) ) +geom_point()

# look up p values of the difference in communication between resistant and sensitive tumors
#load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/When cells communicate_communication changes AllArms/trendsby response AllArms.RData")
load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Cohort2 When cells communicate_communication changes AllArms/Cohort2 trendsby response AllArms.RData")
assessmentSupervised <- assessment[Pair.Name%in%names(pp)][Treat=="CombinationRibo"] [Model=="one way anova_DayStart_Response"][LigandPhenoCelltype=="Cancer cells"][ReceptorPhenoCelltype=="Macrophages"][coefficient!="(Intercept)"][order(Pair.Name)]
assessmentSupervised$adjust.p.val <- NULL
assessmentSupervised$i <- NULL
assessmentSupervised$adj.r.squared <- NULL
assessmentSupervised[,"Significantly validates from C1":=F]
assessmentSupervised[,"Effect direction validates from C1":=F]
assessmentSupervised[pval<0.05&Estimate<0,"Significantly validates from C1":=T]
assessmentSupervised[Estimate<0,"Effect direction validates from C1":=T]
#assessmentSupervised[Estimate>0]
write.csv(assessmentSupervised,file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Final communication files/Cancer to macropage communication supervised analysis/Cohort 2 Cancer to macropage M2 stimulating communication supervised analysis.csv")
pheatmap::pheatmap((pp[,assessmentSupervised[pval<0.05]$Pair.Name]),cluster_rows=F,cluster_cols=F,na_col="white",annotation_row = annotation_row, annotation_colors = annotation_col,scale="none",show_rownames=F,annotation_legend=F,border_color=NA)

assessmentSupervised <- assessment[Pair.Name%in%names(pp)][Treat!="CombinationRibo"] [Model=="one way anova_DayStart_Response"][LigandPhenoCelltype=="Cancer cells"][ReceptorPhenoCelltype=="Macrophages"][coefficient!="(Intercept)"][order(Pair.Name)]
assessmentSupervised$adjust.p.val <- NULL
assessmentSupervised$i <- NULL
assessmentSupervised$adj.r.squared <- NULL
assessmentSupervised[,"Significantly validates from C1":=F]
assessmentSupervised[,"Effect direction validates from C1":=F]
assessmentSupervised[pval<0.05&Estimate<0,"Significantly validates from C1":=T]
assessmentSupervised[Estimate<0,"Effect direction validates from C1":=T]
write.csv(assessmentSupervised,file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Final communication files/Cancer to macropage communication supervised analysis/Cohort 2 Letrozole Cancer to macropage M2 stimulating communication supervised analysis.csv")

# Select cancer - macrophage communications for communication pathways linked to M2 differentiation.
corruptsigs<-CCI[Day!=180][Treat!="CombinationRibo"][LigandPhenoCelltype=="Cancer cells"][ReceptorPhenoCelltype=="Macrophages"]
nsm1<- unique(corruptsigs[Day!=180]$Pair.Name)
extractThese <- c(grep("FGF2_",nsm1),grep("FGF1_",nsm1),grep("IL11_",nsm1),grep("IL12A_",nsm1),
                  grep("CLCF1_",nsm1),grep("MADCAM1_",nsm1),grep("SEMA3",nsm1),grep("SEMA4",nsm1),grep("PLAU_",nsm1),grep("GDF9_",nsm1))
Supervisedcorruptsigs <- corruptsigs[Pair.Name%in%assessmentSupervised$Pair.Name] 
corruptsigs <- data.table(Supervisedcorruptsigs%>%select(Pair.Name, scalelnTransductionMu, Patient.Study.ID,dynamic_class3,Day)%>% spread(Pair.Name,scalelnTransductionMu,fill=NA))[order(dynamic_class3,Patient.Study.ID,Day)]

# Plot as heatmap
pp <- data.frame( corruptsigs%>% select(-c(Patient.Study.ID,dynamic_class3,Day))  )
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
#pheatmap::pheatmap((pp),cluster_rows=F,annotation_row = annotation_row, annotation_colors = annotation_col,scale="none",show_rownames=F,annotation_legend=F,border_color=NA)
pheatmap::pheatmap((pp),cluster_rows=F,cluster_cols=F,na_col="white",annotation_row = annotation_row, annotation_colors = annotation_col,scale="none",show_rownames=F,annotation_legend=F,border_color=NA)


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

