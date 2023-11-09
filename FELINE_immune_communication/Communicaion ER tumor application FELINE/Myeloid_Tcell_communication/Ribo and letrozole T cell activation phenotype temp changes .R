rm(list=ls())
require(mgcv);require(data.table);require(dplyr);require(ggplot2);require(tidyr);require(lme4);require(lmerTest);require(parallel)
library(effects);require(umap)
require(Rfast);require(ider)
library("dendextend");library(ggdendro);require(ggsci);require(viridis)
require("Rdimtools")

# Load clinical data
load(file= "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE/Feline clinical input.RData" )

# Load metadata for a specific cell type
#gsea_path <- "~/Dropbox/FELINE Project/Data_analysis/scRNA/05_ssGSEA_score/Signature_c2_hallmark/results/"   #"~/Dropbox/FELINE/Data_share/Modeling_Data/pathway_seperated_files/Data_gene_count_per_celltype_model_zinbwave_ssGSEA/FEL001043/"
gsea_path <- "~/Dropbox/FELINE Project (1)/Data_analysis/scRNA/05_ssGSEA_score/Signature_c7/results/"
  
  
  
cell_types_all <-c("T_cells")
Subtype <- c("CD8+ T cells","NK cells")
DAY=c(0,14,180)
ARMS <-c("A","B","C")
metadd <- rbindlist( lapply(cell_types_all , function (cell_type_i){
  #load cell metadata#file.exists("~/Dropbox/FELINE Project/Data_analysis/scRNA/01_metadata/results/FEL001046_meta_for_each_celltype/FEL001046_Cancer_cells_scRNA.metadata.txt")
  annotation.file <- paste0("~/Dropbox/FELINE Project (1)/Data_analysis/scRNA/01_metadata/results/FEL001046_meta_for_each_celltype/FEL001046_",cell_type_i,"_scRNA.metadata.txt")
  cell_type_meta_dd <- data.table(  fread(annotation.file))[Celltype_subtype!="Low-quality cells"]
  setnames( cell_type_meta_dd,old="Sample",new="Sample_p_t")
  cell_type_meta_dd[, c("Sample", "Timepoint") := tstrsplit(Sample_p_t, "_", fixed=TRUE)]
  cell_type_meta_dd[,Day:=0]  ; cell_type_meta_dd[grepl("_M",Sample_p_t),Day:=14]   ; cell_type_meta_dd[grepl("_E",Sample_p_t),Day:=180]   
  cell_type_meta_dd[,day_fact:=as.factor(Day)]
  cell_type_meta_dd[,file_string:=cell_type_i]
  # merge response scores
  FULL_cell_type_meta_dd <- merge(cell_type_meta_dd,response_code_dd,by="Sample")        #FULL_cell_type_meta_dd[Celltype_subtype=="Fibroblasts",Celltype_subtype:="CAF-S1"]
  #print(unique(FULL_cell_type_meta_dd$Celltype_subtype))
}) )[Platform!="ICELL8"][Celltype_subtype%in%Subtype]


# Load normalized gsea scores
# ssgsealist <- lapply(cell_types_all , function (cell_type_i){
#   
#   cell_type_gsea_adj <- data.table( readRDS(file=paste0(gsea_path,"FEL001046_",cell_type_i,"_scRNA.zinbwave.normalized.ssGSEA_scores.RDS" ))) 
#   setnames( cell_type_gsea_adj,old="Gene Set",new="Gene_Set") #cell_type_gsea_adj[1:5,1:5]
#   CELL_Subtype<- unique(metadd[file_string==cell_type_i]$Celltype_subtype)  #"Cancer cells"
#   FULL_cell_type_meta_dd_CT <- metadd[Celltype_subtype==CELL_Subtype][Day==DAY][ARM%in%ARMS]
#   cell_type_gsea_adj_CT <- cell_type_gsea_adj[,c("Gene_Set",FULL_cell_type_meta_dd_CT$Cell.ID),with=FALSE]  
#   return(cell_type_gsea_adj_CT)
# })
#ssgsea.FEL001046.T_cells.zbw.txt
# Load normalized gsea scores
ssgsealist <- lapply(cell_types_all , function (cell_type_i){
  cell_type_gsea_adj <- data.table( fread(file=paste0(gsea_path, "ssgsea.FEL001046.", cell_type_i, ".zbw.txt" ))) 
  
  #cell_type_gsea_adj <- data.table( readRDS(file=paste0(gsea_path, "FEL001046_", cell_type_i, "_scRNA.zinbwave.normalized.ssGSEA_scores.RDS" ))) 
  setnames( cell_type_gsea_adj, old="Gene Set", new="Gene_Set" )     #cell_type_gsea_adj[1:5,1:5]
  CELL_Subtype <- unique( metadd[file_string == cell_type_i]$Celltype_subtype )  
  FULL_cell_type_meta_dd_CT <- metadd[ Celltype_subtype %in% CELL_Subtype ][ Day %in% DAY ][ ARM %in% ARMS ]
  cell_type_gsea_adj_CT <- cell_type_gsea_adj[ ,c("Gene_Set", FULL_cell_type_meta_dd_CT$Cell.ID), with=FALSE]  
  return(cell_type_gsea_adj_CT)
})
ssgsealist[[1]][1:10,1:10]
# can look for intersecting genes if jointly analysing multiple cell types
common_ssgsea <- ssgsealist[[1]]$Gene_Set # intersect(ssgsealist[[1]]$Gene_Set, ssgsealist[[2]]$Gene_Set)
n_ssgsea <- length(common_ssgsea)
ssgsealist[[1]] <- ssgsealist[[1]][Gene_Set %in% common_ssgsea]#ssgsealist[[2]] <- ssgsealist[[2]][Gene_Set%in%common_ssgsea]
SSGSEA <- ssgsealist[[1]]#cbind( ssgsealist[[1]],ssgsealist[[2]][,-1])
rm(list=c("ssgsealist","CELL_Subtype", "cell_type_gsea_adj" , "cell_type_gsea_adj_CT", "cell_type_i" , "FULL_cell_type_meta_dd",  "cell_type_meta_dd",  "FULL_cell_type_meta_dd_CT"  ))

SSGSEA[1:10,1:10]

# get full ssgsea data set in a format to merge with umap data output
all_1 <- SSGSEA [, names(SSGSEA) %in% c("Gene_Set", metadd[Day %in% DAY][ARM %in% ARMS]$Cell.ID), with=FALSE ]
allgene_dd_transp <- t(as.matrix( all_1, rownames="Gene_Set" )) ;#colnames(allgene_dd_transp) <- gene_matrix$Gene.ID
allgene_dd_transp<-data.table(allgene_dd_transp,keep.rownames = T)
setnames(allgene_dd_transp,old="rn",new="Cell.ID")

metagsea<-merge(metadd,allgene_dd_transp,by="Cell.ID")
metagsea[,Treatmentlab:= "Combination ribociclib"]
metagsea[ARM=="A",Treatmentlab:= "Letrozole alone"]

names(metagsea)[grepl("NK_CELL",names(metagsea))]
names(metagsea)[grepl("NK_CELL",names(metagsea))&grepl("CD8",names(metagsea)) ]
names(metagsea)[grepl("EFFEC",names(metagsea))&grepl("CD8",names(metagsea)) ]

TcellPhenotytpedata <- metagsea%>%select(Cell.ID,Day,Sample_p_t,dynamic_class3,ARM,Treatmentlab,GSE22886_NAIVE_CD8_TCELL_VS_NKCELL_UP)%>%mutate%>%mutate(Tcellactivation=-GSE22886_NAIVE_CD8_TCELL_VS_NKCELL_UP )%>%select(-GSE22886_NAIVE_CD8_TCELL_VS_NKCELL_UP)
#save( TcellPhenotytpedata ,file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/T cell ssgsea phenotype C1/T cell ssgsea phenotype C1.RData")
load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/T cell ssgsea phenotype C1/T cell ssgsea phenotype C1.RData")

ggplot( TcellPhenotytpedata[] , aes(y=Tcellactivation, x=log(1+Day), col=dynamic_class3, fill=dynamic_class3, group=interaction(dynamic_class3) ))+
  theme_classic(base_size=18)+facet_wrap(~Treatmentlab,ncol=1)+
  geom_violin(alpha=0.6,aes(group=interaction(dynamic_class3, Day)),position=position_dodge() )+#geom_smooth(method="lm")+
  geom_point(position =position_jitterdodge(dodge.width=2.3,jitter.width=0.0)  )+
  theme(aspect.ratio=1)+
  scale_x_continuous(breaks=log(1+c(0,14,180)) ,labels=c(0,14,180))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant","Sensitive"))+
  scale_color_npg(name="Tumor response", labels=c("Resistant","Sensitive"))+
  labs(y="CD8 T cell  cytotoxicity phenotype \n (GSE22886 NK vs naive ssGSEA signature)", x="Day") +
  geom_smooth(method="gam", formula=y~s(x,k=3),se=F)

summar<- TcellPhenotytpedata%>%group_by(Treatmentlab,Day,dynamic_class3,Patient.Study.ID)%>%summarise(Tcellactivation=mean(Tcellactivation))

ggplot( TcellPhenotytpedata[] , aes(y=Tcellactivation, x=log(1+Day), col=dynamic_class3, fill=dynamic_class3, group=interaction(dynamic_class3) ))+
  theme_classic(base_size=26)+facet_wrap(~Treatmentlab,nrow=2)+
  stat_boxplot(geom ='errorbar',position=position_dodge(),aes(group=interaction(dynamic_class3, Day)) ) + 
  geom_boxplot(alpha=0.6,outlier.shape=NA,aes(group=interaction(dynamic_class3, Day)),position=position_dodge() )+#geom_smooth(method="lm")+
  geom_point(position =position_jitterdodge(dodge.width=2.3,jitter.width=0.0)  )+
  theme(aspect.ratio=1)+scale_x_continuous(breaks=log(1+c(0,14,180)) ,labels=c(0,14,180))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant","Sensitive"))+
  scale_color_npg(name="Tumor response", labels=c("Resistant","Sensitive"))+
  labs(y="CD8+ T cell  cytotoxicity phenotype \n (GSE22886 NK vs naive ssGSEA signature)", x="Day") 
paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Ribo and Letrozole Activated vs Naive T cell cytotoxicity phenotype over time boxplot.png"),height=10,width=10)

ggplot( TcellPhenotytpedata[] , aes(y=Tcellactivation, x=log(1+Day), col=dynamic_class3, fill=dynamic_class3, group=interaction(dynamic_class3) ))+
  theme_classic(base_size=26)+facet_wrap(~Treatmentlab,nrow=2)+
  stat_boxplot(geom ='errorbar',position=position_dodge(),aes(group=interaction(dynamic_class3, Day)) ) + 
  geom_boxplot(alpha=0.6,outlier.shape=NA, aes(group=interaction(dynamic_class3, Day)),position=position_dodge() )+#geom_smooth(method="lm")+
  geom_point(position =position_jitterdodge(dodge.width=2.3,jitter.width=0.0)  )+
  theme(aspect.ratio=1)+scale_x_continuous(breaks=log(1+c(0,14,180)) ,labels=c(0,14,180))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant","Sensitive"))+
  scale_color_npg(name="Tumor response", labels=c("Resistant","Sensitive"))+
  labs(y="CD8+ T cell  cytotoxicity phenotype \n (GSE22886 NK vs naive ssGSEA signature)", x="Day") +
  theme(axis.title=element_blank(),  axis.text=element_blank(),strip.text = element_blank(),legend.title=element_blank(),legend.text=element_blank() ,legend.position = "none")

paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"BLANK Ribo and Letrozole Activated vs Naive T cell cytotoxicity phenotype over time boxplot.png"),height=10,width=10)

TcellPhenotytpedata[,Daylab:= paste0("Day ", Day)]
TcellPhenotytpedata[,TumorResponselab:= "Resistant"]
TcellPhenotytpedata[dynamic_class3=="Response",TumorResponselab:= "Sensitive"]

ggplot( TcellPhenotytpedata[] , aes(y=Tcellactivation, x=TumorResponselab, col=dynamic_class3, fill=dynamic_class3, group=interaction(dynamic_class3) ))+
  theme_classic(base_size=18)+facet_wrap(Treatmentlab~Daylab)+
  geom_boxplot(alpha=0.6,aes(group=interaction(dynamic_class3, Day)),position=position_dodge() )+#geom_smooth(method="lm")+
  geom_point(position =position_jitterdodge(dodge.width=2.3,jitter.width=0.0)  )+
  theme(aspect.ratio=1)+#scale_x_continuous(breaks=log(1+c(0,14,180)) ,labels=c(0,14,180))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant","Sensitive"))+
  scale_color_npg(name="Tumor response", labels=c("Resistant","Sensitive"))+
  labs(y="CD8 T cell  cytotoxicity phenotype \n (GSE22886 NK vs naive ssGSEA signature)", x="Day") 


# ggplot( TcellPhenotytpedata[] , aes(y=Tcellactivation, x=log(1+Day), col=dynamic_class3, fill=dynamic_class3, group=interaction(dynamic_class3) ))+
#   theme_classic(base_size=18)+facet_wrap(~Treatmentlab)+
#   geom_boxplot(alpha=0.6,aes(group=interaction(dynamic_class3, Day)),position=position_dodge() )+#geom_smooth(method="lm")+
#   geom_point(position =position_jitterdodge(dodge.width=2.3,jitter.width=0.0)  )+
#   theme(aspect.ratio=1)+scale_x_continuous(breaks=log(1+c(0,14,180)) ,labels=c(0,14,180))+
#   scale_fill_npg(name="Tumor response", labels=c("Resistant","Sensitive"))+
#   scale_color_npg(name="Tumor response", labels=c("Resistant","Sensitive"))+
#   labs(y="CD8 T cell  cytotoxicity phenotype \n (GSE22886 NK vs naive ssGSEA signature)", x="Day") 

ggplot( TcellPhenotytpedata[] , aes(y=Tcellactivation, x=log(1+Day), col=dynamic_class3, fill=dynamic_class3, group=interaction(dynamic_class3) ))+
  theme_classic(base_size=28)+facet_wrap(~Treatmentlab,ncol=1)+
  geom_violin(alpha=0.6,aes(group=interaction(dynamic_class3, Day)),position=position_dodge() )+#geom_smooth(method="lm")+
  geom_point(position =position_jitterdodge(dodge.width=2.3,jitter.width=0.0)  )+
  theme(aspect.ratio=1)+scale_x_continuous(breaks=log(1+c(0,14,180)) ,labels=c(0,14,180))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant","Sensitive"))+
  scale_color_npg(name="Tumor response", labels=c("Resistant","Sensitive"))+
  labs(y="CD8 T cell  cytotoxicity phenotype \n (GSE22886 NK vs naive ssGSEA signature)", x="Day") #+
  #geom_smooth(method="gam", formula=y~s(x,k=3),se=F)

ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Tcell communication phenotypes ALLArms/Ribo and Letrozole NK vs naive cell cytotoxicity phenotype of T cells over time.png",width=10)

paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Ribo and Letrozole Activated vs Naive T cell cytotoxicity phenotype over time.png"),height=10,width=10)

ggplot( TcellPhenotytpedata[] , aes(y=Tcellactivation, x=log(1+Day), col=dynamic_class3, fill=dynamic_class3, group=interaction(dynamic_class3) ))+
  theme_classic(base_size=28)+facet_wrap(~Treatmentlab,ncol=1)+
  geom_violin(alpha=0.6,aes(group=interaction(dynamic_class3, Day)),position=position_dodge() )+#geom_smooth(method="lm")+
  geom_point(position =position_jitterdodge(dodge.width=2.3,jitter.width=0.0)  )+
  theme(aspect.ratio=1)+scale_x_continuous(breaks=log(1+c(0,14,180)) ,labels=c(0,14,180))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant","Sensitive"))+
  scale_color_npg(name="Tumor response", labels=c("Resistant","Sensitive"))+
  labs(y="CD8 T cell  cytotoxicity phenotype \n (GSE22886 NK vs naive ssGSEA signature)", x="Day") +
  #geom_smooth(method="gam", formula=y~s(x,k=3),se=F)+
  theme(axis.title=element_blank(),  axis.text=element_blank(),strip.text = element_blank(),legend.title=element_blank(),legend.text=element_blank() ,legend.position = "none")

paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"BLANK Ribo and Letrozole Activated vs Naive T cell cytotoxicity phenotype over time.png"),height=10,width=10)







ggplot( TcellPhenotytpedata[] , aes(y=Tcellactivation, x=log(1+Day), col=dynamic_class3, fill=dynamic_class3, group=interaction(dynamic_class3) ))+
  theme_classic(base_size=18)+facet_wrap(~Treatmentlab)+
  geom_boxplot(alpha=0.6,aes(group=interaction(dynamic_class3, Day)),position=position_dodge() )+#geom_smooth(method="lm")+
  geom_point(position =position_jitterdodge(dodge.width=2.3,jitter.width=0.0)  )+
  theme(aspect.ratio=1)+scale_x_continuous(breaks=log(1+c(0,14,180)) ,labels=c(0,14,180))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant","Sensitive"))+
  scale_color_npg(name="Tumor response", labels=c("Resistant","Sensitive"))+
  labs(y="CD8 T cell  cytotoxicity phenotype \n (GSE22886 NK vs naive ssGSEA signature)", x="Day") +
  geom_smooth(method="gam", formula=y~s(x,k=3),se=F)



ggplot( metagsea[] , aes(y=GSE27786_CD8_TCELL_VS_NKTCELL_UP, x=log(1+Day), fill=dynamic_class3, group=interaction(dynamic_class3, Day) ))+
  theme_classic(base_size=18)+
  geom_violin(position=position_dodge(),scale="width" )+#geom_smooth(method="lm")+
  geom_point(position =position_jitterdodge(dodge.width=2.3,jitter.width=0.0)  )+
  theme(aspect.ratio=1)+scale_x_continuous(breaks=log(1+c(0,14,180)) ,labels=c(0,14,180))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant","Sensitive"))+
  labs(y="CD8 T cell  cytotoxicity phenotype \n (GSE22886 NK vs naive ssGSEA signature)", x="Day") +
  facet_wrap(~Treatmentlab)
