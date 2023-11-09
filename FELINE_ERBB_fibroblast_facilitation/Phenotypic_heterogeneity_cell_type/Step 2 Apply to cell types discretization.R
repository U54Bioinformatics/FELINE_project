rm(list=ls())
require(data.table);require(dplyr);require(ggplot2);library(MASS);require(arules)
# Load clinical data
load(file="~/Dropbox/Cancer_pheno_evo/data/FELINE/Feline clinical input.RData")

# load counts per million read data for a specific patient's samples and gen Ligand+Receptor genes
UMAPlocs <- "~/Dropbox/Cancer_pheno_evo/data/FELINE2/PhenotypesAllArms/"

UMAPfiles <- c("AdipocytesAdipocytes.RData",
               "B_cellsPlasma cells.RData",
               "Endothelial_cellsVas-Endo.RData",
               "FibroblastsFibroblasts.RData",
               "MacrophagesMacrophages.RData",
               "Cancer_cellsCancer cells.RData",
               "Normal_epithelial_cellsNormal epithelial cells.RData", 
               "PericytesPericytes.RData" ,
               "T_cellsCD4+ T cells.RData",                       
               "T_cellsCD8+ T cells.RData")#list.files(UMAPlocs)
nCellTypes <- length(UMAPfiles)


for(x in 1: length(UMAPfiles)){
  # Load umap output data 
  load(paste0(UMAPlocs,UMAPfiles[x]))
  
  #Number of phenotype components
  NN <-umap_data_tt$config$n_components
  PhenoVars<-paste0("V",1:NN)
  UMAPonly <- data.table( u_dat%>%dplyr::select(Cell.ID:TreatCodeOrd,PhenoVars ) %>%dplyr::select(-c(classification_hclust, classification_hclust2, Type, TypeM, TypeM2,Arm.x,Arm.y,Class,      Response )))
  # discretize each phenotype
  # Adip nbreaks=6
  # Bcell nbreaks=6
  # Endo nbreaks=6
  if(NN>3){nbreaks=5}else{nbreaks=5}

  discretedd <- data.table(discretizeDF(UMAPonly%>%dplyr::select(PhenoVars),default = list(method = "interval",breaks = nbreaks,labels = FALSE)))
  setnames(discretedd,old=colnames(discretedd),new=paste0("Disc_",colnames(discretedd)))  
  disc_levels <- expand.grid(data.table(sapply(1:NN,function(x){ 1:nbreaks })))
  setnames(disc_levels,old=colnames(disc_levels),new=colnames(discretedd))
  freqtable<-table(discretedd)
  median(table(discretedd))
  hist( freqtable[freqtable>0],breaks=20)
  median(freqtable[freqtable>0])
  discretedd[, key_ := paste0(cell_types_all,"_",Subtype[1] ,"_",do.call(paste, c(.SD, sep = "_"))), .SDcols = names(discretedd)]
  UMAPonlydiscretized <- data.table( UMAPonly,discretedd)  
  outloc<-"~/Dropbox/Cancer_pheno_evo/data/FELINE2/UMAP reduced dimensions ALL ARMS/"
  save(UMAPonlydiscretized,disc_levels,freqtable,file=paste0(outloc,"UMAPonlydiscret_",UMAPfiles[x]))
  print(x/length(UMAPfiles))
}