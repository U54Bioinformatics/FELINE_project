rm(list=ls())
require(mgcv);require(data.table);require(dplyr);require(ggplot2);require(tidyr);require(lme4);require(lmerTest);require(parallel)
library(effects);require(umap)
require(Rfast);require(ider)
library("dendextend");library(ggdendro);require(ggsci);require(viridis)
require("Rdimtools")

# Load clinical data
load(file= "~/Dropbox/Cancer_pheno_evo/data/FELINE/Feline clinical input.RData" )

# Load metadata for a specific cell type
gsea_path <- "~/Dropbox/FELINE Project/Data_analysis/scRNA/05_ssGSEA_score/Signature_c2_hallmark/results/"   

cell_types_all <-c("B_cells")
DAY=c(0,14,180)
ARMS <-c("A","B","C")
metadd <- rbindlist( lapply(cell_types_all , function (cell_type_i){
  #load cell metadata
  annotation.file <- paste0( "~/Dropbox/FELINE Project/Data_analysis/scRNA/01_metadata/results/FEL001046_meta_for_each_celltype/FEL001046_", cell_type_i, "_scRNA.metadata.txt")
  cell_type_meta_dd <- data.table(  fread(annotation.file))[Celltype_subtype!="Low-quality cells"]
  setnames( cell_type_meta_dd,old="Sample",new="Sample_p_t")
  cell_type_meta_dd[, c("Sample", "Timepoint") := tstrsplit(Sample_p_t, "_", fixed=TRUE)]
  cell_type_meta_dd[,Day:=0]  ; cell_type_meta_dd[grepl("_M",Sample_p_t),Day:=14]   ; cell_type_meta_dd[grepl("_E",Sample_p_t),Day:=180]   
  cell_type_meta_dd[,day_fact:=as.factor(Day)]
  cell_type_meta_dd[,file_string:=cell_type_i]
  # merge response scores
  FULL_cell_type_meta_dd <- merge(cell_type_meta_dd,response_code_dd,by="Sample")      
}) )[Platform!="ICELL8"]


# Load normalized gsea scores
ssgsealist <- lapply(cell_types_all , function (cell_type_i){
  cell_type_gsea_adj <- data.table( readRDS(file=paste0(gsea_path, "FEL001046_", cell_type_i, "_scRNA.zinbwave.normalized.ssGSEA_scores.RDS" ))) 
  setnames( cell_type_gsea_adj, old="Gene Set", new="Gene_Set" )     
  CELL_Subtype <- unique( metadd[file_string == cell_type_i]$Celltype_subtype )  
  FULL_cell_type_meta_dd_CT <- metadd[ Celltype_subtype %in% CELL_Subtype ][ Day %in% DAY ][ ARM %in% ARMS ]
  cell_type_gsea_adj_CT <- cell_type_gsea_adj[ ,c("Gene_Set", FULL_cell_type_meta_dd_CT$Cell.ID), with=FALSE]  
  return(cell_type_gsea_adj_CT)
})

# can look for intersecting genes if jointly analysing multiple cell types
common_ssgsea <- ssgsealist[[1]]$Gene_Set 
n_ssgsea <- length(common_ssgsea)
ssgsealist[[1]] <- ssgsealist[[1]][Gene_Set %in% common_ssgsea]#ssgsealist[[2]] <- ssgsealist[[2]][Gene_Set%in%common_ssgsea]
SSGSEA <- ssgsealist[[1]]#cbind( ssgsealist[[1]],ssgsealist[[2]][,-1])
rm(list=c("ssgsealist","CELL_Subtype", "cell_type_gsea_adj" , "cell_type_gsea_adj_CT", "cell_type_i" , "FULL_cell_type_meta_dd",  "cell_type_meta_dd",  "FULL_cell_type_meta_dd_CT"  ))

## downsample cells to train umap if large numbers of cells in some samples to prevent overrepresentation
downsam.metadd <- data.table( metadd[ Day %in% DAY][ ARM %in% ARMS ][ order( -nCount_RNA ) ] %>% 
                                group_by(Sample, orig.ident, Sample_p_t, Celltype, Celltype_subtype, Platform, Day, day_fact, Patient.Study.ID, ARM) %>%
                                slice( 1:1000 ) ) [!is.na(Cell.ID)]
subcel_gen_sub <- SSGSEA [ ,][, c("Gene_Set", downsam.metadd$Cell.ID ), with=FALSE ] 
# subset data and put into correct matrix input format for dimensionality estimation and will need transposing for umap
dim_est_input <- as.matrix( subcel_gen_sub , rownames= "Gene_Set" )


## Perform Intrinsic dimension estimation
dim_estout1 <- est.boxcount(dim_est_input , nlevel = 30, cut = c(0.05, 0.95))
dim_estout1.df <- data.table(data.frame(dim_estout1),delta_eps=diff(log(dim_estout1$r))[1])
dim_estout1.df[,x:=log(1/dim_estout1$r)]
dim_estout1.df[,y:=log(dim_estout1$Nr)]
dim_estout1.df[,delta_y:=c(diff(y),NA)]
dim_estout1.df[delta_y<0,delta_y:=0]
dim_estout1.df[,D:=- delta_y/delta_eps] 

overall_D <- round( median(na.omit(dim_estout1.df[D>0])$D) )
ggplot(dim_estout1.df,aes(y= - delta_y/delta_eps,x=r )) + 
  geom_point() +  
  geom_line(aes(y= D,x=r)) +
  geom_hline(aes(yintercept = overall_D))

##### Umap nonlinear dimension reduction
# subsetted data transposed for input into umap
ss_dd <- t( dim_est_input ) ; rm(list="dim_est_input")
NN <- 2
set.seed(123); umap_data_tt <- umap(ss_dd, n_components=NN, n_neighbors=20) # v good
umap_MOD <- umap_data_tt

cor(umap_data_tt$layout)
pairs(umap_data_tt$layout)

# Reshape all ssgsea data
fullcel_gen_sub <- SSGSEA [Gene_Set %in% colnames(ss_dd) ,][, names(SSGSEA) %in% c("Gene_Set", metadd[Day %in% DAY][ARM %in% ARMS]$Cell.ID ), with=FALSE ] 
full_dd_transp <- t( as.matrix(fullcel_gen_sub, rownames="Gene_Set" ) ) ;#colnames(full_dd_transp) <- EXAM_genes
# get full ssgsea data set in a format to merge with umap data output
all_1 <- SSGSEA [, names(SSGSEA) %in% c("Gene_Set", metadd[Day %in% DAY][ARM %in% ARMS]$Cell.ID), with=FALSE ]
allgene_dd_transp <- t(as.matrix( all_1, rownames="Gene_Set" )) ;

# Checks
#SSGSEA[Gene_Set=="ABBUD_LIF_SIGNALING_1_DN"][,FEL002_C25288_104770_C36_R02]
#all_1[Gene_Set=="ABBUD_LIF_SIGNALING_1_DN"][,FEL002_C25288_104770_C36_R02]
#allgene_dd_transp["FEL002_C25288_104770_C36_R02","ABBUD_LIF_SIGNALING_1_DN"]
#full_dd_transp["FEL002_C25288_104770_C36_R02","ABBUD_LIF_SIGNALING_1_DN"]
#ss_dd["FEL002_C25288_104770_C36_R02","ABBUD_LIF_SIGNALING_1_DN"]

# look at pca of dimension reduced space
pc_rot <- prcomp(umap_data_tt$layout)$x[,1:(ncol(umap_data_tt$layout)-1)]
if(nrow(umap_data_tt$layout) == nrow(metadd[Day %in% DAY][ARM %in% ARMS]) ){
  u_dat <- merge(metadd[Day %in% DAY][ARM %in% ARMS],
                 data.table(#Cell.ID=rownames(ss_dd), umap_data_tt$layout, ss_dd)
                   Cell.ID=rownames(full_dd_transp), umap_data_tt$layout, allgene_dd_transp )
                 ,by="Cell.ID") 
}else{
  u_dat <- merge(metadd[Day %in% DAY][ARM %in% ARMS],
                 data.table(#Cell.ID=rownames(ss_dd), umap_data_tt$layout, ss_dd)
                   Cell.ID=rownames(full_dd_transp), predict(umap_data_tt,data= full_dd_transp  ), allgene_dd_transp )
                 ,by="Cell.ID") 
}


methodcor= "spearman"
corVall <-cor(u_dat[]%>%dplyr::select(paste0("V",1:NN)),u_dat[]%>%dplyr::select(colnames(ss_dd)),method=methodcor)#corV1<-cor(u_dat$V1,full_dd_transp)
Subtype<-metadd[1]$Celltype_subtype
save(umap_data_tt,u_dat,SSGSEA,DAY,cell_types_all,ARMS,common_ssgsea,n_ssgsea,corVall,Subtype,
     file=paste0("~/Dropbox/Cancer_pheno_evo/data/FELINE2/PhenotypesAllArms/new",cell_types_all,Subtype[1],".RData"))

corResponse_Vs <- cor(log(u_dat[u_dat$V2<5][]$prop_change),u_dat[u_dat$V2<5][]%>%dplyr::select(V1:V2))
methodcor= "pearson"
methodcor= "spearman"
corV1 <-cor(u_dat[u_dat$V2<5][]$V1,u_dat[u_dat$V2<5]%>%dplyr::select(colnames(ss_dd)),method=methodcor)#corV1<-cor(u_dat$V1,full_dd_transp)
corV1ord<-corV1[,order(-abs(corV1))]
corV2<-cor(u_dat[u_dat$V2<5][]$V2,u_dat[u_dat$V2<5]%>%dplyr::select(colnames(ss_dd)),method=methodcor)#corV2<-cor(u_dat$V2,full_dd_transp[,])
corV2ord<-corV2[,order(-abs(corV2))]
nn=20
corV1ord[1:nn]
corV2ord[1:nn]


