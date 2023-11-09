rm(list=ls())
require(data.table)
require(dplyr)
require(ggplot2)
require(tidyr)
require(ggsci)
### Load clinical data
#load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE/Feline clinical input.RData") #save(Clin_resp_dd,Clin_resp_dd_class,Clin_resp_dd_classAdd,patientCode_LU,res_file_46,response_code_dd,response_dd,file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE/Feline clinical input.RData") #gsea_path<-"~/Dropbox/FELINE Project/Data_analysis/scRNA/05_ssGSEA_score/Signature_c2_hallmark/results/"   #"~/Dropbox/FELINE/Data_share/Modeling_Data/pathway_seperated_files/Data_gene_count_per_celltype_model_zinbwave_ssGSEA/FEL001043/"

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
allphenotypes <- merge(allphenotypes %>%dplyr::select(-c( paste0("V", 1:5), "nCount_RNA", "nFeature_RNA", "Infercnv_CNA", "Platform","Timepoint","Sample_p_t","prop_change", "file_string", "day_fact" , "Ribo","TreatLab", "Burden_t0" ,"BurdenTracked" ,"Day0", "DayLastBurd", "Dose_lab","TreatCode","TreatCodeOrd","dynamic_class2","rgr_A","rgr_B")), tmp, by= c("Celltype", "Celltype_subtype") )


### Unique clusters of cells (ALL) and their umap discretization level
uu <- unique( allphenotypes %>% dplyr::select( c("Celltype","PhenoCelltype", "key_", paste0("Disc_V", 1:5) ) ) )  #"Celltype_subtype",


### Load Ligand Receptor database list of Ramilowski et al 2015
#load( "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/LigandReceptor/Filtered_Human-2015-Ramilowski-LR-pairs.RData")
#LRgenelist <- unique( c(LRpairsFiltered$HPMR.Receptor, LRpairsFiltered$HPMR.Ligand) )


### Load counts per million read data for a specific patient's samples and gen Ligand+Receptor genes
CPMlocs <- "/Users/jason/Dropbox/FELINE Project (1)/Data_analysis/FELINE_data_folder/scRNA_count_CPM/output/"
CPMfiles <- list.files(CPMlocs)[ grep("CPM", list.files(CPMlocs) ) ]
n10Xpats <- length(CPMfiles)
# settings
shouldscale <- FALSE ## should cpm be scaled
macrocelldd<-allphenotypes[Celltype=="Macrophages"]
#saveloc <- "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/CommunicationOutputAllArms/"

## Select a patient index to load cpm data
fulldd<-rbindlist(lapply(1:length(CPMfiles) ,function(i){
  ## Load gene expression of all macrophages
  cat("patient index    ");cat(i);cat("    loading data    ") #i= 2 #patient index
  cpm_i <- data.table( fread( paste0(CPMlocs, CPMfiles[i]) ) )   # load full gene expression
  t_cpm_i <- cpm_i[, data.table(t(.SD), keep.rownames=TRUE), .SDcols=-"Gene.ID"]
  colnames(t_cpm_i) <- c("Cell.ID",cpm_i$Gene.ID)
  y<- merge(macrocelldd,t_cpm_i,by="Cell.ID")
  return(y)  
}))  
fulldd[1:10,1:50]

##  Perform umap on log expression of genes with non zero coverage in >5% of cells
umap_in  <- fulldd %>% dplyr::select(-c(colnames(macrocelldd)))
gene_summary <- data.table( GeneName= colnames(umap_in),Mean= colMeans(umap_in), Variance= colVars( as.matrix(umap_in) ),coverage= colMeans(umap_in>0), nCellsExpressed= colSums(umap_in>0) )
#gene_summary[Mean>0][coverage>0.05]$Mean %>% log()%>%hist()
#gene_summary[Mean>0][coverage>0.05]$Variance %>% log()%>%hist()
#hist( colMeans(umap_in[,1:10]) )
#hist( colVars(umap_in[,1:10]) )
umap_in2 <-  umap_in[,gene_summary[Mean>0][coverage>0.05]$GeneName ,with=FALSE ]
umap_mod2 <- umap::umap(log(1+umap_in2),n_components=4)
u_dat <- cbind(fulldd,umap_mod2$layout)

## Gather results
#umap_mod <- umap::umap(umap_in,n_components=2)
#cor(umap_mod2$layout)
#pairs(umap_mod2$layout)
u_dat <- cbind(fulldd,umap_mod2$layout)
u_dat[,MyloidStatus:= "Macrophage"]
u_dat[Celltype_subtype=="DC",MyloidStatus:= "Dendritic cell"]
names(u_dat ) <- gsub("-","_", names(u_dat))

## Load saved data
load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/UMAP genes macrophage only ALLARMS/UMAP genes macrophage only ALLARMS.RData")


M1 <- lmer(V2~dynamic_class3+(1|Patient.Study.ID) , data=u_dat[Day!=180][ARM!="A"][!(Celltype_subtype!="DC"&V1>2&V2< -1)])
summary(M1)
M1b <- lmer(V2~dynamic_class3+(1|Patient.Study.ID) , data=u_dat[Day!=180][ARM=="A"][!(Celltype_subtype!="DC"&V1>2&V2< -1)])
summary(M1b)





# visualize
sz<-1.25
ggplot( u_dat[Day!=180][ARM!="A"][!(Celltype_subtype!="DC"&V1>2&V2< -1)],aes(-V1,V2,col=  dynamic_class3)) +  
  geom_point(size=2.5)+scale_color_npg(name="Tumor response")+
  theme_classic()+
  theme(aspect.ratio=1,axis.text=element_blank(),axis.title=element_blank(),legend.position="none")
ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Macrophage phenotypes AllArms/BLANK Macrophage Ribociclib TumorResponse UMAP AllArms Cohort1 pts1_5.png")
paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
ggsave(paste0(paperfile,"BLANK Ribo and Letrozole RIBO Early treatment myeloid phenotype UMAP.png"))

ggplot( u_dat[Day==0][ARM=="A"][!(Celltype_subtype!="DC"&V1>2&V2< -1)],aes(-V1,V2,col=  dynamic_class3)) +  
  geom_point(size=2.5)+scale_color_npg(name="Tumor response")+
  theme_classic()+
  theme(aspect.ratio=1,axis.text=element_blank(),axis.title=element_blank(),legend.position="none")
ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Macrophage phenotypes AllArms/BLANK Macrophage Letrozole TumorResponse UMAP AllArms Cohort1 pts1_5.png")
ggsave(paste0(paperfile,"BLANK Ribo and Letrozole LETROZOLE Early treatment myeloid phenotype UMAP.png"))

ggplot( u_dat[!(Celltype_subtype!="DC"&V1>2&V2< -1)],aes(-V1,V2,col=  Celltype_subtype)) +  
  geom_point(size=sz)+#scale_color_ucscgb(name="Myleoid cell types")+
  theme_classic()+
  theme(aspect.ratio=1,axis.text=element_blank(),axis.title=element_blank(),legend.position="none")
ggsave(paste0(paperfile,"BLANK Ribo and Letrozole LETROZOLE myeloid classification UMAP.png"))


ggplot( u_dat[!(Celltype_subtype!="DC"&V1>2&V2< -1)] ,aes(-V1,V2,col=  log(1+CD36))) +  
  geom_point(size=sz)+scale_color_viridis_c(name="Myleoid cell types")+
  theme_classic()+
  theme(aspect.ratio=1,axis.text=element_blank(),axis.title=element_blank(),legend.text=element_blank(),legend.title=element_blank(),
        legend.key.size = unit(1.25, 'cm'))
ggsave(paste0(paperfile,"BLANK Ribo and Letrozole LETROZOLE myeloid M2 differentiation UMAP.png"))


corgenes[order(-V2)]



ggplot( u_dat[ARM!="A"][!(Celltype_subtype!="DC"&V1>2&V2< -1)],aes(-V1,V2,col=  dynamic_class3)) +  
  geom_point(size=1.5)+scale_color_npg(name="Tumor response",labels=c("Resistant","Sensitive"))+labs(x="UMAP 1",y="UMAP 2")+
  theme_classic()+
  theme(aspect.ratio=1)
ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Macrophage phenotypes AllArms/Macrophage Ribociclib TumorResponse UMAP AllArms Cohort1 pts1_5.png")

ggplot( u_dat[ARM=="A"][!(Celltype_subtype!="DC"&V1>2&V2< -1)],aes(-V1,V2,col=  dynamic_class3)) +  
  geom_point(size=1.5)+scale_color_npg(name="Tumor response",labels=c("Resistant","Sensitive"))+labs(x="UMAP 1",y="UMAP 2")+
  theme_classic()+
  theme(aspect.ratio=1)
ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Macrophage phenotypes AllArms/Macrophage Letrozole TumorResponse UMAP AllArms Cohort1 pts1_5.png")



ggplot( u_dat[ARM=="A"][!(Celltype_subtype!="DC"&V1>2&V2< -1)][Patient.Study.ID!="2972-005-202" ],aes(-V1,V2,col=  dynamic_class3)) +  
  geom_point(size=1.5)+
  #scale_color_npg(name="Tumor response",labels=c("Resistant","Sensitive"))+labs(x="UMAP 1",y="UMAP 2")+
  theme_classic()+
  theme(aspect.ratio=1)

corgenes[order(V2)][1:10]
ggplot( u_dat[ARM=="A"][!(Celltype_subtype!="DC"&V1>2&V2< -1)][],aes(-V1,V2,col=  log(1+ESR1))) +  
  geom_point(size=1.5)+
  #scale_color_npg(name="Tumor response",labels=c("Resistant","Sensitive"))+labs(x="UMAP 1",y="UMAP 2")+
  theme_classic()+
  theme(aspect.ratio=1)


ggplot( u_dat[!(Celltype_subtype!="DC"&V1>2&V2< -1)],aes(-V1,V2,col=  Celltype_subtype)) +  
  geom_point(size=sz)+#scale_color_ucscgb(name="Myleoid cell types")+
  theme_classic()+
  theme(aspect.ratio=1,axis.text=element_blank(),axis.title=element_blank(),legend.position="none")
#ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Macrophage phenotypes AllArms/BLANK MacrophageClassificationUMAP AllArms Cohort1.png")
ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Macrophage phenotypes AllArms/BLANK MacrophageClassificationUMAP AllArms Cohort1 pts1_5.png")

ggplot( u_dat[!(Celltype_subtype!="DC"&V1>2&V2< -1)],aes(-V1,V2,col=  Celltype_subtype)) +  
  geom_point(size=sz)+ labs(x="UMAP 1",y="UMAP 2")+
  scale_color_discrete(name="Myleoid cell types")+
  theme_classic()+
  theme(aspect.ratio=1)
ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Macrophage phenotypes AllArms/MacrophageClassificationUMAP AllArms Cohort1 pts1_5.png")


ggplot( u_dat[!(Celltype_subtype!="DC"&V1>2&V2< -1)] ,aes(-V1,V2,col=  log(1+CD36))) +  
  geom_point(size=sz)+scale_color_viridis_c(name="Myleoid cell types")+
  theme_classic()+
  theme(aspect.ratio=1,axis.text=element_blank(),axis.title=element_blank(),legend.text=element_blank(),legend.title=element_blank(),
        legend.key.size = unit(1.25, 'cm'))
ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Macrophage phenotypes AllArms/BLANK MacrophagePhenotypeM1vsM2 UMAP AllArms Cohort1 pts1_5.png")

ggplot( u_dat[!(Celltype_subtype!="DC"&V1>2&V2< -1)] ,aes(-V1,V2,col=  log(1+CD36))) +  
  geom_point(size=sz)+scale_color_viridis_c(name="CD36 \n expression \n (log(1+x))")+
  theme_classic()+ labs(x="UMAP 1",y="UMAP 2")+
  theme(aspect.ratio=1,
        legend.key.size = unit(1.25, 'cm'))
ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Macrophage phenotypes AllArms/MacrophagePhenotypeM1vsM2 UMAP AllArms Cohort1 pts1_5.png")


# 
# 
# 
# 
# # ggplot( u_dat ,aes(-V1,V2,col=  log(1+CYP27A1))) +  
# #   geom_point(size=1.5)+scale_color_viridis_c(name="Myleoid cell types")+
# #   theme_classic()+
# #   theme(aspect.ratio=1,axis.text=element_blank(),axis.title=element_blank(),legend.text=element_blank(),legend.title=element_blank(),
# #         legend.key.size = unit(1.25, 'cm'))
# # ggplot( u_dat ,aes(-V1,V2,col=  log(1+DHRS9))) +  
# #   geom_point(size=1.5)+scale_color_viridis_c(name="Myleoid cell types")+
# #   theme_classic()+
# #   theme(aspect.ratio=1,axis.text=element_blank(),axis.title=element_blank(),legend.text=element_blank(),legend.title=element_blank(),
# #         legend.key.size = unit(1.25, 'cm'))
# 
# ggplot( u_dat[ARM!="A"],aes(-V1,V2,col=  dynamic_class3)) +  
#   geom_point(size=1.5)+scale_color_npg(name="Tumor response")+
#   theme_classic()+
#   theme(aspect.ratio=1)#axis.text=element_blank(),axis.title=element_blank(),legend.position="none"
# ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Macrophage phenotypes AllArms/MacrophageGeneUmap by Response all days.png")
# 
# 
# 
# ggplot( u_dat[ARM!="A"],aes(-V1,V2,col=  dynamic_class3)) +  
#   geom_point(size=1.5)+scale_color_npg(name="Tumor response")+
#   theme_classic()+
#   theme(aspect.ratio=1,axis.text=element_blank(),axis.title=element_blank(),legend.position="none")
# ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Macrophage phenotypes AllArms/BLANK MacrophageGeneUmap by Response.png")
# 
# ggplot( u_dat[ARM!="A"],aes(-V1,V2,col=  Celltype_subtype)) +  
#   geom_point(size=1.5)+#scale_color_ucscgb(name="Myleoid cell types")+
#   theme_classic()+
#   theme(aspect.ratio=1,axis.text=element_blank(),axis.title=element_blank(),legend.position="none")
# ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Macrophage phenotypes AllArms/BLANK MacrophageGeneUmap by Immune classifier.png")
# 
# 
# ggplot( u_dat[ARM!="A"],aes(-V1,V2,col=  log(1+CD36))) +  
#   geom_point(size=1.5)+scale_color_viridis_c(name="Myleoid cell types")+
#   theme_classic()+
#   theme(aspect.ratio=1,axis.text=element_blank(),axis.title=element_blank(),legend.text=element_blank(),legend.title=element_blank(),
#         legend.key.size = unit(1.25, 'cm'))
# ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Macrophage phenotypes AllArms/BLANK MacrophageGeneUmap by M2 marker CD36 expression.png")
# 
# ggplot( u_dat[ARM!="A"],aes(-V1,V2,col=  log(1+LRP1))) +  
#   geom_point(size=1.5)+scale_color_viridis_c(name="Myleoid cell types")+
#   theme_classic()+
#   theme(aspect.ratio=1,axis.text=element_blank(),axis.title=element_blank(),legend.text=element_blank(),legend.title=element_blank(),
#         legend.key.size = unit(1.25, 'cm'))
# ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Macrophage phenotypes AllArms/BLANK MacrophageGeneUmap by M2 marker LRP1 expression.png")
# 
# ggplot( u_dat[ARM!="A"],aes(-V1,V2,col=  log(1+NRP1))) +  
#   geom_point(size=1.5)+scale_color_viridis_c(name="Myleoid cell types")+
#   theme_classic()+
#   theme(aspect.ratio=1,axis.text=element_blank(),axis.title=element_blank(),legend.text=element_blank(),legend.title=element_blank(),
#         legend.key.size = unit(1.25, 'cm'))
# ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Macrophage phenotypes AllArms/BLANK MacrophageGeneUmap by M2 marker NRP1 expression.png")
# 
# ggplot( u_dat[ARM!="A"],aes(-V1,V2,col=  log(1+SDC2))) +  
#   geom_point(size=1.5)+scale_color_viridis_c(name="Myleoid cell types")+
#   theme_classic()+
#   theme(aspect.ratio=1,axis.text=element_blank(),axis.title=element_blank(),legend.text=element_blank(),legend.title=element_blank(),
#         legend.key.size = unit(1.25, 'cm'))
# #theme(aspect.ratio=1,axis.text=element_blank(),axis.title=element_blank(),legend.position="none")
# ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Macrophage phenotypes AllArms/BLANK MacrophageGeneUmap by M2 marker SDC2 expression.png")
# 
# ggplot( u_dat[ARM!="A"],aes(-V1,V2,col=  log(1+IQGAP2))) +  
#   geom_point(size=1.5)+scale_color_viridis_c(name="Myleoid cell types")+
#   theme_classic()+
#   theme(aspect.ratio=1,axis.text=element_blank(),axis.title=element_blank(),legend.text=element_blank(),legend.title=element_blank(),
#         legend.key.size = unit(1.25, 'cm'))
# #theme(aspect.ratio=1,axis.text=element_blank(),axis.title=element_blank(),legend.position="none")
# ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Macrophage phenotypes AllArms/BLANK MacrophageGeneUmap by M2 marker IQGAP2 expression.png")
# 
# 
# 
# ggplot( u_dat[ARM!="A"][],aes(-V1,V2,col=  log(1+CCR5))) +  
#   geom_point(col="grey",size=0.25)+scale_color_viridis_c(name="Myleoid cell types")+
#   theme_classic()+
#   theme(aspect.ratio=1,axis.text=element_blank(),axis.title=element_blank(),legend.text=element_blank(),legend.title=element_blank(),
#         legend.key.size = unit(1.25, 'cm'))+
#   geom_point(data= u_dat[ARM!="A"][CCR5>0],size=1.5)
# ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Macrophage phenotypes AllArms/BLANK MacrophageGeneUmap by M2 marker CCR5 expression.png")
# 
# 
# 
# ggplot( u_dat[ARM!="A"],aes(-V1,V2,col=  log(1+CD36))) +  
#   geom_point(size=1.5)+scale_color_viridis_c(name="Myleoid cell types")+
#   theme_classic()+
#   theme(aspect.ratio=1, legend.key.size = unit(1.25, 'cm'))
# ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Macrophage phenotypes AllArms/MacrophageGeneUmap by M2 marker CD36 expression.png")
# 
# ggplot( u_dat[ARM!="A"],aes(-V1,V2,col=  log(1+LRP1))) +  
#   geom_point(size=1.5)+scale_color_viridis_c(name="Myleoid cell types")+
#   theme_classic()+
#   theme(aspect.ratio=1,  legend.key.size = unit(1.25, 'cm'))
# ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Macrophage phenotypes AllArms/MacrophageGeneUmap by M2 marker LRP1 expression.png")
# 
# ggplot( u_dat[ARM!="A"],aes(-V1,V2,col=  log(1+NRP1))) +  
#   geom_point(size=1.5)+scale_color_viridis_c(name="Myleoid cell types")+
#   theme_classic()+
#   theme(aspect.ratio=1,  legend.key.size = unit(1.25, 'cm'))
# ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Macrophage phenotypes AllArms/MacrophageGeneUmap by M2 marker NRP1 expression.png")
# 
# ggplot( u_dat[],aes(-V1,V2,col=  log(1+SDC2))) +  
#   geom_point(size=1.5)+scale_color_viridis_c(name="Myleoid cell types")+
#   theme_classic()+
#   theme(aspect.ratio=1,  legend.key.size = unit(1.25, 'cm'))
# #theme(aspect.ratio=1,axis.text=element_blank(),axis.title=element_blank(),legend.position="none")
# ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Macrophage phenotypes AllArms/MacrophageGeneUmap by M2 marker SDC2 expression.png")
# 
# ggplot( u_dat[],aes(-V1,V2,col=  log(1+IQGAP2))) +  
#   geom_point(size=1.5)+scale_color_viridis_c(name="Myleoid cell types")+
#   theme_classic()+ 
#   theme(aspect.ratio=1,  legend.key.size = unit(1.25, 'cm'))
# #theme(aspect.ratio=1,axis.text=element_blank(),axis.title=element_blank(),legend.position="none")
# ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Macrophage phenotypes AllArms/MacrophageGeneUmap by M2 marker IQGAP2 expression.png")
# 
# ggplot( u_dat[], aes(-V1, V2, col=  log(1+CCR5))) +  
#   geom_point(col= "grey", size= 0.25)+ geom_point(data= u_dat[ARM!="A"][CCR5>0], size= 1.5)+ 
#   scale_color_viridis_c(name= "Expression \n (log (1+x))")+
#   theme_classic()+
#   theme(aspect.ratio= 1, legend.key.size= unit(1.25, 'cm'))
# ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Macrophage phenotypes AllArms/MacrophageGeneUmap by M2 marker CCR5 expression.png")
# 
# 
# ggplot( u_dat[], aes(-V1, V2, col=  log(1+ALK))) +  
#   geom_point(col= "grey", size= 0.25)+ geom_point(data= u_dat[ARM!="A"][ALK>0], size= 1.5)+ 
#   scale_color_viridis_c(name= "Expression \n (log (1+x))")+
#   theme_classic()+
#   theme(aspect.ratio= 1, legend.key.size= unit(1.25, 'cm'))
# ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Macrophage phenotypes AllArms/MacrophageGeneUmap by M2 marker ALK expression.png")
# 
# 
# 
# ggplot( u_dat[], aes(-V1, V2, col=  log(1+CSF1R))) +  
#   geom_point(col= "grey", size= 0.25)+ geom_point(data= u_dat[ARM!="A"][CSF1R>0], size= 1.5)+ 
#   scale_color_viridis_c(name= "Expression \n (log (1+x))")+
#   theme_classic()+
#   theme(aspect.ratio= 1, legend.key.size= unit(1.25, 'cm'))
# ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Macrophage phenotypes AllArms/MacrophageGeneUmap by M2 marker CSF1R expression.png")
# 
# 
# ggplot( u_dat[], aes(-V1, V2, col=  log(1+CSF1))) +  
#   geom_point(col= "grey", size= 0.25)+ geom_point(data= u_dat[ARM!="A"][CSF1>0], size= 1.5)+ 
#   scale_color_viridis_c(name= "Expression \n (log (1+x))")+
#   theme_classic()+
#   theme(aspect.ratio= 1, legend.key.size= unit(1.25, 'cm'))
# ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Macrophage phenotypes AllArms/MacrophageGeneUmap by M2 marker CSF1 expression.png")
# 
# 
# 
# 
# 
# ggplot( u_dat[ARM!="A"],aes(-V1,V2,col=  Celltype_subtype)) +  
#   geom_point(size=1.5)+scale_color_discrete(name="Myleoid cell types")+
#   theme_classic()+
#   theme(aspect.ratio=1) #axis.text=element_blank(),axis.title=element_blank(),legend.position="none"
# ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Macrophage phenotypes AllArms/MacrophageGeneUmap by Immune classifier annotation.png")
# 
# 
# 
# 
# 
# 
# 
# ggplot( u_dat[ARM!="A"],aes(y=V2,x=log(1+Day) ,col=dynamic_class3,fill=  dynamic_class3)) +  geom_violin(aes(group=interaction(Day,dynamic_class3)),alpha=0.6)+theme_classic()+
#   geom_point(col="black",shape=21, position = position_jitterdodge(dodge.width=2.3,jitter.width=0.1) )+scale_color_npg(name="Tumor response")+scale_fill_npg(name="Tumor response")+ylab("Macrophage polarisation (high =M2 resolving, low=M1 inflammatory)")#+facet_wrap(~Day)
# ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Macrophage phenotypes AllArms/Macrophageimmunesuppressivestate by Response and Day.png")
# 
# 
# ggplot( u_dat[ARM!="A"],aes(y=V2,x=log(1+Day) ,col=dynamic_class3,fill=  dynamic_class3)) +  
#   geom_violin(aes(group=interaction(Day,dynamic_class3)),alpha=0.6)+theme_classic()+
#   geom_point(col="black",shape=21, position = position_jitterdodge(dodge.width=2.3,jitter.width=0.1) )+
#   scale_color_npg(name="Tumor response")+scale_fill_npg(name="Tumor response")+
#   theme(aspect.ratio=1,axis.text=element_blank(),axis.title=element_blank(),legend.position="none")
# 
# ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Macrophage phenotypes AllArms/BLANK Macrophageimmunesuppressivestate by Response and Day.png")
# 
# 
# 
# ggplot( u_dat[ARM!="A"],aes(y=V2,x=log(1+Day) ,col=dynamic_class3,fill=  dynamic_class3,group=interaction(Day,dynamic_class3)))+
#   theme_classic(base_size=18)+
#   geom_violin(position=position_dodge() )+#geom_smooth(method="lm")+
#   geom_point(col="black",position =position_jitterdodge(dodge.width=2.3,jitter.width=0.0)  )+
#   theme(aspect.ratio=1) + scale_x_continuous(breaks=log(1+c(0,14,180)) , labels=c(0,14,180))+
#   scale_fill_npg(name="Tumor response", labels=c("Resistant","Sensitive"))+
#   theme(axis.title = element_blank(),axis.text = element_blank(), legend.position = "none")
# ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Macrophage phenotypes AllArms/BLANK Macrophageimmunesuppressivestate by Response and Day2.png")
