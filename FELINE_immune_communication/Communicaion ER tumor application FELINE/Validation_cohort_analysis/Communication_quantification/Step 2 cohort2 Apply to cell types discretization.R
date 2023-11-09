rm(list=ls())
require(data.table)
require(dplyr)
require(ggplot2)
library(MASS)
require(arules)
# Load clinical data
#save(Clin_resp_dd,Clin_resp_dd_class,Clin_resp_dd_classAdd,patientCode_LU,res_file_46,response_code_dd,response_dd,file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE/Feline clinical input.RData")
load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE/Feline clinical input.RData")
#gsea_path<-"~/Dropbox/FELINE Project/Data_analysis/scRNA/05_ssGSEA_score/Signature_c2_hallmark/results/"   #"~/Dropbox/FELINE/Data_share/Modeling_Data/pathway_seperated_files/Data_gene_count_per_celltype_model_zinbwave_ssGSEA/FEL001043/"

# where are the umap phenotype landscapes
UMAPlocs <- "~/Dropbox/Cancer_pheno_evo/data/FELINE2/PhenotypesAllArmsCohort2/"

UMAPfiles <- c("UpdatedRevisednewAdipocytesAdipocytes.RData",
"UpdatedRevisednewB_cellsB cells.RData",
"UpdatedRevisednewCancer_cellsCancer cells.RData",                      
"UpdatedRevisednewEndothelial_cellsLym-Endo.RData" ,                    
 "UpdatedRevisednewFibroblastsFibroblasts.RData" ,                       
 "UpdatedRevisednewMacrophagesDC.RData",                                 
 "UpdatedRevisednewNormal_epithelial_cellsNormal epithelial cells.RData",
 "UpdatedRevisednewPericytesPericytes.RData"     ,                       
 "UpdatedRevisednewT_cellsCD4+ T cells.RData"  ,                         
 "UpdatedRevisednewT_cellsCD8+ T cells.RData" )


nCellTypes <- length(UMAPfiles)

# x=4
# # Load umap output data 
# load(paste0(UMAPlocs,UMAPfiles[x]))
# # identify genes highly correlated with a phenotype
# informative_ssgsea <- colnames(corVall)[colSums( abs(corVall)>0.4 )>0]
# outloc<-"/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/InformativeSSGSEAforUMAP/"
# save(informative_ssgsea,file=paste0(outloc,"informSSGSEA",UMAPfiles[x] ) )

#x=6
for(x in 1: length(UMAPfiles)){
  # Load umap output data 
  load(paste0(UMAPlocs,UMAPfiles[x]))
  
  #u_dat[1:10,1:50]
  #Number of phenotype components
  NN <-umap_data_tt$config$n_components
  PhenoVars<-paste0("V",1:NN)
  UMAPonly <- data.table( u_dat%>%dplyr::select(Cell.ID:file_string,PhenoVars ))
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
  #hist( freqtable[freqtable>0],breaks=20)
  #median(freqtable[freqtable>0])
  discretedd[, key_ := paste0(cell_types_all,"_",Subtype[1] ,"_",do.call(paste, c(.SD, sep = "_"))), .SDcols = names(discretedd)]
  UMAPonlydiscretized <- data.table( UMAPonly,discretedd)
  
  #ggplot(UMAPonlydiscretized,aes(V1,V2,col=interaction(Disc_V1,Disc_V2)))+geom_point()
  #ggplot(UMAPonlydiscretized,aes(V4,V3,col=interaction(Disc_V4,Disc_V3)))+geom_point()
  #ggplot(UMAPonlydiscretized,aes(y=V3,x=Disc_V3,col=interaction(Disc_V3),fill=interaction(Disc_V3)))+geom_violin()
  
  #ggplot(UMAPonlydiscretized, aes(Disc_V3, y = 1, fill = Disc_V3)) +
  #  geom_col(position = 'stack',  show.legend = FALSE)+scale_fill_viridis()+theme_classic(base_size=18)+labs(y="Number of cells", x="Decile")+scale_x_continuous(breaks=1:10)+theme(aspect.ratio=1)
  
  
  
  outloc<-"/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Cohort2UMAP reduced dimensions ALL ARMS/"
  save(UMAPonlydiscretized,disc_levels,freqtable,file=paste0(outloc,"UMAPonlydiscret_",UMAPfiles[x]))
  print(x/length(UMAPfiles))
}



#UMAPinformSSGSEA<-data.table( u_dat%>%dplyr::select(Cell.ID:TreatCodeOrd,paste0("V",1:NN) ,informative_ssgsea) )
# colnames(corVall[abs(corVall)>0.4])# hist(corVall[abs(corVall)>0.4])
