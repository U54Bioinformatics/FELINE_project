rm(list=ls())
require(mgcv);require(data.table);require(dplyr);require(ggplot2);require(tidyr);require(lme4);require(lmerTest);require(parallel)
library(effects);require(umap);require(Rfast);require(ider);library("dendextend");library(ggdendro);require(ggsci);require(viridis);require("Rdimtools")
# Load clinical data
load(file= "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE/Feline clinical input.RData" )
# File path to ssgsea data
gsea_path <- "~/Dropbox/FELINE Project/FELINE Cohort 2/DNANexusCopy/ssgsea/"          #gsea_path <- "~/Dropbox/FELINE Project/Data_analysis/scRNA/05_ssGSEA_score/Signature_c2_hallmark/results/"   #"~/Dropbox/FELINE/Data_share/Modeling_Data/pathway_seperated_files/Data_gene_count_per_celltype_model_zinbwave_ssGSEA/FEL001043/"
gsea_pathC1  <- "~/Dropbox/FELINE Project/Data_analysis/scRNA/05_ssGSEA_score/Signature_c2_hallmark/results/"   #"~/Dropbox/FELINE/Data_share/Modeling_Data/pathway_seperated_files/Data_gene_count_per_celltype_model_zinbwave_ssGSEA/FEL001043/"
# Load cohort 2 metadata
load( file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Cohort2Metadata/Cohort2Metadata.RData") #metadd,cohort1metadd,annotation.file,compdataLU
# metadd[,ARM:="B"]
# what and where are the cohort 1 umap files? 
# UMAPlocs <- "~/Dropbox/Cancer_pheno_evo/data/FELINE2/PhenotypesAllArms/"
# UMAPfiles <- c("AdipocytesAdipocytes.RData",
#                "B_cellsPlasma cells.RData",
#                "Endothelial_cellsVas-Endo.RData",
#                "FibroblastsFibroblasts.RData",
#                "MacrophagesMacrophages.RData",
#                "Cancer_cellsCancer cells.RData",
#                "Normal_epithelial_cellsNormal epithelial cells.RData", 
#                "PericytesPericytes.RData" ,
#                "T_cellsCD4+ T cells.RData",                       
#                "T_cellsCD8+ T cells.RData")#list.files(UMAPlocs)
# nCellTypes <- length(UMAPfiles)

# data to use: Cell type, timepoint and treatment
cell_types_all <-c("T_cells");   #### CHANGE HERE compdataLU
Subtype <- c("CD4+ T cells","Tregs")
DAY=c(0,14,180);     ARMS <-c("A","B","C")
metadd[file_string == cell_types_all][ Celltype_subtype %in% Subtype ]$Celltype_subtype%>%unique()
## Load normalized gsea scores
ssgsealist <- lapply(cell_types_all , function (cell_type_i){
  # cell types in cohort 2 that link to cohort 1 cell type
  CELL_Subtype <- unique( metadd[file_string == cell_type_i]$Celltype_subtype )  
  CELL_SubtypeC2 <- unique( metadd[file_string == cell_type_i]$Annotation2 ) 
  
  c1ssgsea<- data.table(readRDS(paste0(gsea_pathC1,"FEL001046_",cell_type_i,"_scRNA.zinbwave.normalized.ssGSEA_scores.RDS") ))$'Gene Set'
  # for each cohort 2 cell type, load ssgsea data, ensure metadata cell type annotations match labels of ssgsea file. Select ssgsea for cells in the metadata file
  gsealist <- lapply(1:length(CELL_SubtypeC2),function(j){
    # Subset metadata for that cell type
    FULL_cell_type_meta_dd_CT <- metadd[ Annotation2 %in% CELL_SubtypeC2[j] ][ Day %in% DAY ][ ARM %in% ARMS ] #metadd[Cell.ID %in% colnames(cell_type_gsea_adj)[-1] ]$Celltype%>%table()
    # Read HQ & LQ cell subtype ssgsea file
    cell_type_gsea_adj_CT <- unique( data.table( fread(file=paste0(gsea_path,"text_filesAll/feline2_hq_and_hqlq_C2_", CELL_SubtypeC2[j], "_ssgsea.txt.gz" ))) ,by="Gene_set")  [Gene_set%in% c1ssgsea,c("Gene_set", FULL_cell_type_meta_dd_CT$Cell.ID), with=FALSE]   
    setnames( cell_type_gsea_adj_CT, old="Gene_set", new="Gene_Set" )     #cell_type_gsea_adj[1:5,1:5]
    return(cell_type_gsea_adj_CT) #cell_type_gsea_adj_CT[1:10,1:10]
  })
  # Join list of cell type specific datasets together if more that one C2 cell type for the C1 cell type
  if( length(CELL_SubtypeC2) >1 ){
    cell_type_gsea_result <- Reduce(function(...){ merge(... , all=FALSE ) }, gsealist)
  }else{
    cell_type_gsea_result <- gsealist[[1]]
  }
  return(cell_type_gsea_result)
})


## Filter for intersecting pathways across cell types'ssgsea data if jointly analysing multiple cell types
common_ssgsea <- Reduce( intersect, lapply(1:length(ssgsealist), function(x){ ssgsealist[[x]]$Gene_Set }) )  #ssgsealist[[1]]$Gene_Set 
n_ssgsea <- length(common_ssgsea)
for(i in 1:length(ssgsealist)){ssgsealist[[i]] <- ssgsealist[[i]][Gene_Set %in% common_ssgsea] }
SSGSEA <- Reduce(function(...){ merge(... , all=FALSE ,by="Gene_Set") }, ssgsealist)    #SSGSEA <- ssgsealist[[1]]#cbind( ssgsealist[[1]],ssgsealist[[2]][,-1])
rm(list=c("ssgsealist","CELL_Subtype", "cell_type_gsea_adj" , "cell_type_gsea_adj_CT", "cell_type_i" , "FULL_cell_type_meta_dd",  "cell_type_meta_dd",  "FULL_cell_type_meta_dd_CT"  ))

## Downsample cells to train umap if large numbers of cells in some tumor samples to prevent over-representation bias
downsam.metadd <- data.table( metadd[ Day %in% DAY][ ARM %in% ARMS ][file_string%in%cell_types_all][ Celltype_subtype %in% Subtype ][ order( -nCount_RNA ) ] %>% 
                                group_by(Sample, orig.ident, Sample_p_t, Celltype, Celltype_subtype, Platform, Day, day_fact,  ARM) %>%
                                slice( 1:200 ) ) [!is.na(Cell.ID)]
#Geneomic data for subset of downsampled cells
subcel_gen_sub <- SSGSEA [ ,][, c("Gene_Set", downsam.metadd$Cell.ID ), with=FALSE ] 
rm(list="downsam.metadd")

##### Umap nonlinear dimension reduction
# subsetted data transposed for input into umap
ss_dd <- t( as.matrix( subcel_gen_sub , rownames= "Gene_Set" )  ) 
rm(list="subcel_gen_sub")

NN <- 3  ### CHANGE TO MATCH Cohort 1
set.seed(123); umap_data_tt <- umap(ss_dd, n_components=NN, n_neighbors=20)     #cor(umap_data_tt$layout) #pairs(umap_data_tt$layout)

# Reshape all ssgsea data (train and out of sample cells)
fullcel_gen_sub <- SSGSEA [Gene_Set %in% colnames(ss_dd) ,][, names(SSGSEA) %in% c("Gene_Set", metadd[Day %in% DAY][ARM %in% ARMS][file_string%in%cell_types_all][ Celltype_subtype %in% Subtype ]$Cell.ID ), with=FALSE ] 
rm(list="ss_dd")
full_dd_transp <- t( as.matrix(fullcel_gen_sub, rownames="Gene_Set" ) ) ;#colnames(full_dd_transp) <- EXAM_genes
# get full ssgsea data set in a format to merge with umap data output
all_1 <- SSGSEA [, names(SSGSEA) %in% c("Gene_Set", metadd[Day %in% DAY][ARM %in% ARMS][file_string%in%cell_types_all][ Celltype_subtype %in% Subtype ]$Cell.ID), with=FALSE ]
allgene_dd_transp <- t(as.matrix( all_1, rownames="Gene_Set" )) ;#colnames(allgene_dd_transp) <- gene_matrix$Gene.ID
rm(list=c("all_1","fullcel_gen_sub"))
# Checks #SSGSEA[Gene_Set=="ABBUD_LIF_SIGNALING_1_DN"][,FEL002_C25288_104770_C36_R02] #allgene_dd_transp["FEL002_C25288_104770_C36_R02","ABBUD_LIF_SIGNALING_1_DN"] #full_dd_transp["FEL002_C25288_104770_C36_R02","ABBUD_LIF_SIGNALING_1_DN"] #ss_dd["FEL002_C25288_104770_C36_R02","ABBUD_LIF_SIGNALING_1_DN"]

#Combine all metadd, umap and ssgsea data together. Predict out of sample data if umap trained by downsampling
if(nrow(umap_data_tt$layout) == nrow(metadd[Day %in% DAY][ARM %in% ARMS][file_string%in%cell_types_all][ Celltype_subtype %in% Subtype ]) ){
  u_dat <- merge(metadd[Day %in% DAY][ARM %in% ARMS][file_string%in%cell_types_all][ Celltype_subtype %in% Subtype ],
                 data.table(Cell.ID=rownames(full_dd_transp), umap_data_tt$layout, allgene_dd_transp )
                 ,by="Cell.ID") 
}else{
  u_dat <- merge(metadd[Day %in% DAY][ARM %in% ARMS][file_string%in%cell_types_all][ Celltype_subtype %in% Subtype ],
                 data.table(Cell.ID=rownames(full_dd_transp), predict(umap_data_tt,data= full_dd_transp  ), allgene_dd_transp )
                 ,by="Cell.ID") 
}

# Save output
Subtype <- unique(metadd[file_string%in%cell_types_all][ Celltype_subtype %in% Subtype ]$Celltype_subtype)
save(umap_data_tt,u_dat,SSGSEA,DAY,cell_types_all,ARMS,common_ssgsea,n_ssgsea,Subtype,
     file=paste0("~/Dropbox/Cancer_pheno_evo/data/FELINE2/PhenotypesAllArmsCohort2/UpdatedRevisednew",cell_types_all,Subtype[1],".RData"))
# 
# names(u_dat) [grepl("REACTOME", names(u_dat))]
# ggplot(u_dat[ARM!="A"][abs(V1)<2],aes(V1,V4,col=HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION))+theme_classic()+geom_point(size=1.5,alpha=0.5)+facet_grid(Day~dynamic_class3)
# 
# cor(u_dat[ARM!="A"]$ REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION , u_dat[ARM!="A"] %>%select(V1,V2,V3,V4,V5))
# cordd<-data.table(t(cor(u_dat[ARM!="A"] %>%dplyr::select(V1,V2,V3,V4,V5) , u_dat[ARM!="A"]%>%dplyr::select(common_ssgsea) )),keep.rownames = T)
# cordd[order(-abs(V1))]
# cordd[order(-abs(V2))]
# cordd[order(-abs(V3))]
# cordd[order(-abs(V4))]
# cordd[order(-abs(V5))]
# 
# ggplot(u_dat[ARM!="A"][abs(V1)<2], aes(V1,BIOCARTA_HER2_PATHWAY ) ) + theme_classic()+geom_point(size=1.5,alpha=0.5)
# ggplot(u_dat[ARM!="A"][abs(V2)<2], aes(V2,BIOCARTA_HER2_PATHWAY ) ) + theme_classic()+geom_point(size=1.5,alpha=0.5)
