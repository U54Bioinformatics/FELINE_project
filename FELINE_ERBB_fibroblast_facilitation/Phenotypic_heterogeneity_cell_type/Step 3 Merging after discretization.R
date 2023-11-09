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
               "T_cellsCD8+ T cells.RData")

umapDImRedloc<-"~/Dropbox/Cancer_pheno_evo/data/FELINE2/UMAP reduced dimensions ALL ARMS/"
umapDImRedfiles <- list.files(umapDImRedloc)
nCellTypes <- length(umapDImRedfiles)

allphenotypes <- rbindlist(lapply(1:nCellTypes,function(x){
  load(file=paste0(umapDImRedloc,"UMAPonlydiscret_",UMAPfiles[x]))
  #UMAPonlydiscretized,disc_levels,freqtable,
  if(is.null(UMAPonlydiscretized$V3)){ UMAPonlydiscretized$V3<-NA}
  if(is.null(UMAPonlydiscretized$V4)){ UMAPonlydiscretized$V4<-NA}
  if(is.null(UMAPonlydiscretized$V5)){ UMAPonlydiscretized$V5<-NA}
  if(is.null(UMAPonlydiscretized$V6)){ UMAPonlydiscretized$V6<-NA}
  if(is.null(UMAPonlydiscretized$V7)){ UMAPonlydiscretized$V7<-NA}
  
  if(is.null(UMAPonlydiscretized$Disc_V3)){ UMAPonlydiscretized$Disc_V3<-NA}
  if(is.null(UMAPonlydiscretized$Disc_V4)){ UMAPonlydiscretized$Disc_V4<-NA}
  if(is.null(UMAPonlydiscretized$Disc_V5)){ UMAPonlydiscretized$Disc_V5<-NA}
  if(is.null(UMAPonlydiscretized$Disc_V6)){ UMAPonlydiscretized$Disc_V6<-NA}
  if(is.null(UMAPonlydiscretized$Disc_V7)){ UMAPonlydiscretized$Disc_V7<-NA}
  
  yy<-UMAPonlydiscretized%>%dplyr::select(Cell.ID:TreatCodeOrd, key_,paste0("V",1:6),paste0("Disc_V",1:6))
  return(yy)
}))


allphenotypes%>%group_by(Celltype_subtype)%>% dplyr::summarise(n=n() )
allphenotypes%>%group_by(Celltype_subtype)%>% dplyr::summarise(n=sum(!is.na(V7)) )
allphenotypes[,V6:=NULL]
allphenotypes[,Disc_V6:=NULL]

save(allphenotypes,UMAPlocs ,UMAPfiles,umapDImRedloc,umapDImRedfiles,nCellTypes,file="~/Dropbox/Cancer_pheno_evo/data/FELINE2/PhenotypesOfAllCellTypesAllArms/PhenotypesOfAllCellTypesAllArms.RData")


