rm(list=ls())
require(data.table)
require(dplyr)
require(ggplot2)
require(tidyr)
require(ggsci)
require(umap)

### Load clinical data#load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE/Feline clinical input.RData") #save(Clin_resp_dd,Clin_resp_dd_class,Clin_resp_dd_classAdd,patientCode_LU,res_file_46,response_code_dd,response_dd,file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE/Feline clinical input.RData") #gsea_path<-"~/Dropbox/FELINE Project/Data_analysis/scRNA/05_ssGSEA_score/Signature_c2_hallmark/results/"   #"~/Dropbox/FELINE/Data_share/Modeling_Data/pathway_seperated_files/Data_gene_count_per_celltype_model_zinbwave_ssGSEA/FEL001043/"
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

load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Cancer gene expression wrangling/Cancer summary gene expression.RData") #gene_summary
#gene_summary[coverage>0.01]
expressedGenes <- gene_summary[][coverage>0.01]$GeneName

load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Gene expression landscapes/cancer gene expression landscape/cancer gene expression landscape.RData") # u_dat 
#intersect( unique(u_dat[Day==0]$Patient.Study.ID) , unique(u_dat[Day==180]$Patient.Study.ID) )
u_dattrim <- data.table( u_dat[Day %in% c(0,180)] %>% select(c(Cell.ID:PhenoCelltype, one_of(expressedGenes))) )
u_dattrim[,hasDay0:=any(Day==0), by=Patient.Study.ID]
u_dattrim[,hasDay180:=any(Day==180), by=Patient.Study.ID]
u_dattrim[,hasStartEnd:=hasDay0*hasDay180]
u_dattrim[,ncellsSampled:=length(unique(Cell.ID)), by=c("Patient.Study.ID","Day")]
u_dattrim[,ncellsSampledDay0:=sum(Day==0), by=Patient.Study.ID]
u_dattrim[,ncellsSampledDay180:=sum(Day==180), by=Patient.Study.ID]
rm(list="u_dat")
#u_dattrim[1:10,1:30]

lu1 <- unique( data.table(u_dattrim%>% dplyr::select(Patient.Study.ID,Day) ) )

summar_all <- rbindlist(lapply(1: nrow(lu1) , function(i){
  meta_i <- data.table( u_dattrim[Patient.Study.ID==lu1[i]$Patient.Study.ID][Day ==lu1[i]$Day][1,] %>% select(c(Celltype:dynamic_class3,PhenoCelltype)) )
  mean_sqrt_gene_i <-  data.table( t(t(
    colMeans( sqrt(u_dattrim[Patient.Study.ID==lu1[i]$Patient.Study.ID][Day ==lu1[i]$Day]  %>% dplyr::select(c( one_of(expressedGenes)))) )
  )),keep.rownames = T)
  setnames( mean_sqrt_gene_i , old=c("rn","V1"), new=c("gene","mean_sqrtCPM"))
  summar_i <- cbind(meta_i,mean_sqrt_gene_i)  
  print(i)
  return(summar_i)
}))

summar_all[,TreatmentLb2:= "endocrine"]
summar_all[ARM!="A", TreatmentLb2:= "ribocombo"]

GenesToTest <- unique(summar_all$gene)
# Perform pathway analysis 
pathanalysis0 <- rbindlist(lapply( 1:length(GenesToTest), function(pp){
  cat(pp)
  res <- tryCatch({
    m1<- lm(mean_sqrtCPM~ as.factor(Day)*TreatmentLb2 ,    summar_all[gene==GenesToTest[pp]])
    data.table(gene=GenesToTest[pp],  coef( summary(m1) ), keep.rownames = T)
  },
  error= function(x){
    return(NULL)
  })
  return(res)  
} ))
setnames(pathanalysis0, old= c( "Std. Error", "t value", "Pr(>|t|)" ), new= c( "Std.Error", "tval", "pval" ))

#pathanalysis0 <- data.table( pathanalysis0%>%group_by(rn) %>%mutate(adj.pval2= p.adjust(p= pval, method= "fdr" ) ) )
pathanalysis0$adj.pval <- p.adjust(p= pathanalysis0$pval, method= "fdr" )#,n=length(unique(resultout$Pair.Name)))
#pathanalysis0[rn!= "(Intercept)"][adj.pval<0.05][][order(adj.pval)]

pathanalysis0[rn=="as.factor(Day)180"][adj.pval<0.05][order(-abs(Estimate))][1:10] #pathanalysis0[rn=="as.factor(Day)180:TreatmentLb2ribocombo"][adj.pval<0.05][order(-abs(Estimate))]
pathanalysis0[rn=="as.factor(Day)180"][ pval<0.01][order(gene)]
pathanalysis0[rn=="as.factor(Day)180"][ adj.pval<0.05][order(gene)]

toppaths <- pathanalysis0[rn=="as.factor(Day)180"][adj.pval<0.05][order(-abs(Estimate))][1:25]
plotthese <- toppaths$gene
pdat1 <- summar_all[gene%in%plotthese ] 

#save( summar_all,u_dattrim , file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Gene expression landscapes/cancer gene expression landscape/cancer gene expression raw and summarised.RData")
#save( pathanalysis0 , file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Gene expression landscapes/cancer gene expression landscape/Endocrine Letrozole strt vs end cancer gene expression summarised.RData")

pdat1$gene <- factor(pdat1$gene ,  levels= toppaths[order((Estimate))]$gene   )
pdat1$TreatmentLb2 <- factor(pdat1$TreatmentLb2 ,  levels= sort( pdat1$TreatmentLb2%>%unique() ),labels=c("Endocrine alone","Ribociclib combination")   )

ggplot( pdat1 %>% group_by(gene)%>%mutate( scale_mean_sqrtCPM=scale(mean_sqrtCPM) ), aes( y=scale_mean_sqrtCPM, x=gene,col=as.factor(Day), fill=as.factor(Day) ,group=interaction(Day,gene) ))+
  theme_classic(base_size=20)+
  geom_violin(position = position_dodge(width=1),alpha=0.6,scale="width")+geom_point(position = position_dodge(width=1))+
  scale_color_lancet(name="Day")+scale_fill_lancet(name="Day")+facet_wrap(~TreatmentLb2)+coord_flip() +theme(aspect.ratio=2)+
  labs(x="Gene",y="Mean cancer cell expression (scaled mean sqrt CPM)") +
  geom_hline(aes(yintercept = 0), linetype=2)

ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Tumor Mean changes in cancer cell expression in resistant cells.png", height=10, width = 10)


signifChanges <- pathanalysis0[rn=="as.factor(Day)180"][adj.pval < 0.05][order(-abs(Estimate))]
summ_wide <- summar_all[Treatment=="letrozole"][ gene %in% signifChanges$gene ] %>% spread(gene, mean_sqrtCPM)
set.seed(12345)
umap_mod <- umap( summ_wide %>% dplyr::select(one_of(signifChanges$gene )) ,n_components=2 ,n_neighbors=10 )
summ_wide_umap <- cbind(summ_wide, umap_mod$layout)
ggplot( summ_wide_umap , aes( V1, V2 , col= Day ) ) + theme_classic(base_size=18) + theme(aspect.ratio=1) + geom_point(size= 4) + #facet_wrap(~TreatmentLb2) +
  labs(y="Umap B",x="Umap A")
ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/UMAP of signaif Tumor Mean changes in cancer cell expression in resistant cells.png", height=4, width = 4)

corumap <- data.table( t(cor( summ_wide_umap%>%dplyr::select(V1,V2) , summ_wide%>%dplyr::select(one_of(signifChanges$gene )) )) , keep.rownames = T)
corumap[order(-abs(V1))][1:10]
corumap[order(-abs(V2))][1:20]
corumap[order(-(V2))][V2>0.5]
ggplot( summ_wide_umap , aes( x=V1, y= V2 , col= TRIB1 ) ) + theme_classic() + geom_point( size= 4) + labs(y="Umap B",x="Umap A") + scale_color_viridis_c()
ggplot( summ_wide_umap , aes( x=V1, y= V2 , col= DUSP1 ) ) + theme_classic() + geom_point( size= 4) + labs(y="Umap B",x="Umap A") + scale_color_viridis_c()
ggplot( summ_wide_umap , aes( x=V1, y= V2 , col= JUND ) ) + theme_classic() + geom_point( size= 4) + labs(y="Umap B",x="Umap A") + scale_color_viridis_c()


ggplot( summ_wide_umap , aes( x=V1, y= V2 , col= Day ) ) + theme_classic() + geom_point( size= 4) + labs(y="Umap B",x="Umap A")
ggplot( summ_wide_umap , aes( y=V1, x= V2 , col= EGR1 ) ) + theme_classic() + geom_point( size= 4) + labs(y="Umap B",x="Umap A") + scale_color_viridis_c()
ggplot( summ_wide_umap , aes( y=V1, x= V2 , col= FOS ) ) + theme_classic() + geom_point( size= 4) + labs(y="Umap B",x="Umap A") + scale_color_viridis_c()
ggplot( summ_wide_umap , aes( y=V1, x= V2 , col= C3orf67 ) ) + theme_classic() + geom_point( size= 4) + labs(y="Umap B",x="Umap A") + scale_color_viridis_c()
ggplot( summ_wide_umap , aes( y=V1, x= V2 , col= BHLHE40 ) ) + theme_classic() + geom_point( size= 4) + labs(y="Umap B",x="Umap A") + scale_color_viridis_c()


ggplot( summ_wide_umap , aes( y= FOS, x= Day  ) ) + geom_point( size= 4 ) 
ggplot( summ_wide_umap , aes( y= EGR1, x= Day  ) ) + geom_point( size= 4 ) 
ggplot( summ_wide_umap , aes( y= SOX9, x= Day  ) ) + geom_point( size= 4 ) 
ggplot( summ_wide , aes( y= EZH2, x= Day , group=Day ) ) + theme_classic() + geom_boxplot(fill="red", alpha=0.2) + geom_point( size= 4 )  +facet_wrap(~Treatment)+theme(aspect.ratio=1)
ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/EZH2 Tumor Mean changes in cancer cell expression in resistant cells.png", height=4, width = 4)


ggplot( u_dattrim[] , aes(y=sqrt(EZH2), x=Day ,group= Day)  ) + theme_classic() + 
  geom_boxplot(fill="red", alpha=0.2) + 
  #geom_violin(fill="red", alpha=0.2) + 
  geom_point() +facet_wrap(~Treatment)
ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/EZH2 Singel cell changes in cancer cell expression in resistant cells.png", height=4, width = 4)


ggplot( u_dattrim[] , aes(y=sqrt(FOS), x=Day ,group= Day)  ) + theme_classic() + 
  geom_boxplot(fill="red", alpha=0.2) + 
  #geom_violin(fill="red", alpha=0.2) + 
  geom_point() +facet_wrap(~Treatment)
ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/FOS Single cell changes in cancer cell expression in resistant cells.png", height=4, width = 4)

ggplot( u_dattrim[] , aes(y=sqrt(COLEC12), x=Day ,group= Day)  ) + theme_classic() + 
  geom_boxplot(fill="red", alpha=0.2) + 
  #geom_violin(fill="red", alpha=0.2) + 
  geom_point() +facet_wrap(~Treatment)
ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/COLEC12 Single cell changes in cancer cell expression in resistant cells.png", height=4, width = 4)

ggplot( u_dattrim[] , aes(y=sqrt(JUND), x=Day ,group= Day)  ) + theme_classic() + 
  geom_boxplot(fill="red", alpha=0.2) + 
  #geom_violin(fill="red", alpha=0.2) + 
  geom_point() +facet_wrap(~Treatment)
ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/JUND Single cell changes in cancer cell expression in resistant cells.png", height=4, width = 4)
BHLHE40
SOX9

ggplot( u_dattrim[] , aes(y=sqrt(DUSP1), x=Day ,group= Day)  ) + theme_classic() + 
  geom_boxplot(fill="red", alpha=0.2) + 
  #geom_violin(fill="red", alpha=0.2) + 
  geom_point() +facet_wrap(~Treatment)
ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/DUSP1 Single cell changes in cancer cell expression in resistant cells.png", height=4, width = 4)
ggplot( u_dattrim[] , aes(y=sqrt(SOX9), x=Day ,group= Day)  ) + theme_classic() + 
  geom_boxplot(fill="red", alpha=0.2) + 
  #geom_violin(fill="red", alpha=0.2) + 
  geom_point() +facet_wrap(~Treatment)
ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/SOX9 Single cell changes in cancer cell expression in resistant cells.png", height=4, width = 4)

ggplot( u_dattrim[] , aes(y=sqrt(BHLHE40), x=Day ,group= Day)  ) + theme_classic() + 
  geom_boxplot(fill="red", alpha=0.2) + 
  #geom_violin(fill="red", alpha=0.2) + 
  geom_point() +facet_wrap(~Treatment)
ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/BHLHE40 Single cell changes in cancer cell expression in resistant cells.png", height=4, width = 4)



ggplot( u_dattrim[Treatment=="letrozole"] , aes(y=sqrt(FOS),  x=Day ,group= Day)  ) + geom_boxplot(fill="red", alpha=0.2) + geom_point()
ggplot( u_dattrim[Treatment=="letrozole"] , aes(y=sqrt(SOX9), x=Day ,group= Day)  ) + geom_boxplot(fill="red", alpha=0.2) + geom_point()
ggplot( u_dattrim[Treatment=="letrozole"] , aes(y=sqrt(JUN),  x=Day ,group= Day)  ) + geom_boxplot(fill="red", alpha=0.2) + geom_point()

ggplot( u_dattrim , aes(y=sqrt(FOS), x=sqrt(JUN) ,col= Day)  ) +  geom_point()



require(glmnet)

signifChanges <- pathanalysis0[rn=="as.factor(Day)180"][pval < 10.01][order(-abs(Estimate))]
summ_wide <- summar_all[ gene %in% signifChanges$gene ] %>% spread(gene, mean_sqrtCPM)
summ_wide[, response:=0]
summ_wide[Day==180, response:=1]

#glmnet_binom <- glmnet(x= as.matrix(summ_wide %>% dplyr::select( signifChanges$gene ) ), y= summ_wide$response,family= "binomial")
#plot(glmnet_binom)

lassoIterative <- rbindlist( lapply(1:30,function(ii){
  glmnet_binomB <- cv.glmnet(x= as.matrix(summ_wide %>% dplyr::select( signifChanges$gene ) ), y= summ_wide$response,family= "binomial",type.measure="deviance")
  #plot(glmnet_binomB)
  coef_est <- coef(glmnet_binomB,s="lambda.min")
  coefnonzero<-data.table( t(t(coef_est [which(coef_est!=0),])) ,keep.rownames = T)
  coefvals <- data.table(repli=ii, coefnonzero[order(-abs(V1))] )
  return(coefvals)
}) )

summarlassoIterative <- data.table( lassoIterative%>%group_by(rn )%>%summarise(ndetections = length(repli) , meancoef=sum(V1)/ length(unique(lassoIterative$repli))))
plot(ndetections ~     meancoef, summarlassoIterative)
summarlassoIterative[ndetections>0.5*max(ndetections)]

signifChanges[1:20]
ggplot( summar_all[TreatmentLb2=="endocrine"][gene=="NFKBIA"] , aes( y=mean_sqrtCPM, x=Day) )+geom_point()
ggplot( summar_all[TreatmentLb2=="endocrine"][gene=="NFKBIA"] , aes( y=mean_sqrtCPM, x=Day) )+geom_point()
ggplot( summar_all[TreatmentLb2=="endocrine"][gene=="HDAC7"] , aes( y=mean_sqrtCPM, x=Day) )+geom_point()
ggplot( summar_all[TreatmentLb2=="endocrine"][gene=="UROD"] , aes( y=mean_sqrtCPM, x=Day) )+geom_point()

ggplot( summar_all[gene=="FOS"] , aes( y=mean_sqrtCPM, x=Day) )+geom_point()
ggplot( summar_all[gene=="COLEC12"] , aes( y=mean_sqrtCPM, x=Day) )+geom_point()
ggplot( summar_all[gene=="JUND"] , aes( y=mean_sqrtCPM, x=Day) )+geom_point()
ggplot( summar_all[gene=="BHLHE40"] , aes( y=mean_sqrtCPM, x=Day) )+geom_point()
ggplot( summar_all[gene=="PTBP2"] , aes( y=mean_sqrtCPM, x=Day) )+geom_point()

m1<- lm(mean_sqrtCPM~ as.factor(Day)*TreatmentLb2 ,    summar_all[gene=="ERBB4"])
m1<- lm(mean_sqrtCPM~ as.factor(Day) ,    summar_all[gene=="ERBB4"])
summary(m1)


require(lmerTest)













 

u_dat[Day %in% c(0,180)]
[1:10,1:10]

lm(gene~)

# 
# 
# ### Load Ligand Receptor database list of Ramilowski et al 2015
# #load( "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/LigandReceptor/Filtered_Human-2015-Ramilowski-LR-pairs.RData")
# #LRgenelist <- unique( c(LRpairsFiltered$HPMR.Receptor, LRpairsFiltered$HPMR.Ligand) )
# 
# 
# ### Load counts per million read data for a specific patient's samples and gen Ligand+Receptor genes
# CPMlocs <- "/Users/jason/Dropbox/FELINE Project/Data_analysis/FELINE_data_folder/scRNA_count_CPM/output/"
# CPMfiles <- list.files(CPMlocs)[ grep("CPM", list.files(CPMlocs) ) ]
# n10Xpats <- length(CPMfiles)
# 
# 
# # settings
# shouldscale <- FALSE ## should cpm be scaled
# 
# x_celldd<-allphenotypes[PhenoCelltype=="Cancer cells"]
# #saveloc <- "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/CommunicationOutputAllArms/"
# 
# ## Select a patient index to load cpm data
# fulldd<-rbindlist(lapply(1:length(CPMfiles) ,function(i){
#   ## Load gene expression of all macrophages
#   cat("patient index    ");cat(i);cat("    loading data    ") #i= 2 #patient index
#   cpm_i <- data.table( fread( paste0(CPMlocs, CPMfiles[i]) ) )   # load full gene expression
#   t_cpm_i <- cpm_i[, data.table(t(.SD), keep.rownames=TRUE), .SDcols=-"Gene.ID"]
#   colnames(t_cpm_i) <- c("Cell.ID",cpm_i$Gene.ID)
#   y<- merge(x_celldd,t_cpm_i,by="Cell.ID")
#   return(y)  
# }))  
# 
# 
# #fulldd[1:10,1:50]
# 
# umap_in  <- fulldd %>% dplyr::select(-c(colnames(x_celldd)))
# gene_summary <- data.table( GeneName= colnames(umap_in),Mean= colMeans(umap_in), Variance= colVars( as.matrix(umap_in) ),coverage= colMeans(umap_in>0), nCellsExpressed= colSums(umap_in>0) )
# gene_summary[Mean>0][coverage>0.1]$Mean %>% log()%>%hist()
# gene_summary[Mean>0][coverage>0.05]$Variance %>% log()%>%hist()
# rm(list="umap_in")
# 
# #save(gene_summary,file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Cancer gene expression wrangling/Cancer summary gene expression.RData")
# #save(fulldd,file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Cancer gene expression wrangling/Cancer gene expression.RData")
# 
# umap_in2 <-  fulldd[,data.table(gene_summary)[Mean>0][coverage>0.1]$GeneName ,with=FALSE ]
# rm(list="fulldd")
# 
# umap_mod2 <- umap::umap(log(1+umap_in2),n_components=4)
# #save(umap_mod2,file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Cancer gene expression wrangling/Cancer UMAP mod gene expression.RData")
# load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Cancer gene expression wrangling/Cancer gene expression.RData")
# u_dat <- cbind(fulldd,umap_mod2$layout)
# 
# cor(umap_mod2$layout)
# pairs(umap_mod2$layout)
# #u_dat <- cbind(fulldd,umap_mod2$layout)
# #u_dat[,MyloidStatus:= "Macrophage"]
# #u_dat[Celltype_subtype=="DC",MyloidStatus:= "Dendritic cell"]
# names(u_dat ) <- gsub("-","_", names(u_dat))
# 
# rm(list="fulldd")
# rm(list="umap_mod2")
# rm(list="umap_in2")
# 
# #save(u_dat , file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Gene expression landscapes/cancer gene expression landscape/cancer gene expression landscape.RData")
# 
# 
# ggplot( u_dat,aes(V3,V4,col= ARM ) )+  geom_point(size=0.5)+theme_classic()+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# ggplot( u_dat,aes(V2,V4,col= Patient.Study.ID ) )+  geom_point(size=0.5)+theme_classic()+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# 
# 
# which_compare <- which(abs(u_dat$V1)<5&abs(u_dat$V2)<5&abs(u_dat$V3)<5&abs(u_dat$V4)<5)
# corgenes <- data.table(t(
#   cor(u_dat[which_compare,]%>%dplyr::select(V1:V4), log(1+umap_in2[which_compare,]) ,method="spearman")
# ),keep.rownames = T)
# corgenes[order(-abs(V1))][1:20]
# corgenes[order(-abs(V2))][1:20]
# 
# corgenes[order(-abs(V3))][1:10]
# corgenes[order(-abs(V4))][1:10]
# 
# 
# ggplot( u_dat[which_compare,],aes(V1,V2,col= ARM ) )+  geom_point(size=0.5)+theme_classic()+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# ggplot( u_dat[which_compare,],aes(V1,V3,col= ARM ) )+  geom_point(size=0.5)+theme_classic()+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# 
# 
# ### new block
# # rm(list=ls())
# load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Cancer gene expression wrangling/Cancer gene expression.RData")
# load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Cancer gene expression wrangling/Cancer summary gene expression.RData")
# 
# load( "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/LigandReceptor/Filtered_Human-2015-Ramilowski-LR-pairs.RData")
# LRgenelist <- unique( c(LRpairsFiltered$HPMR.Receptor, LRpairsFiltered$HPMR.Ligand) )
# 
# umap_in2 <-  fulldd[,data.table(gene_summary[GeneName%in%LRgenelist])[Mean>0][coverage>0.05]$GeneName ,with=FALSE ]
# rm(list="fulldd")
# rm(list="LRpairs")
# rm(list="LRpairsFiltered")
# 
# umap_mod2 <- umap::umap(log(1+umap_in2),n_components=4)
# rm(list="umap_in2")
# load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Cancer gene expression wrangling/Cancer gene expression.RData")
# u_dat <- cbind(fulldd,umap_mod2$layout)
# rm(list="fulldd")
# 
# #save(u_dat,file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Cancer gene expression wrangling/Cancer LR gene expression umap.RData")
# #save(umap_mod2,file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Cancer gene expression wrangling/Cancer LR gene expression umap model.RData")
# rm(list="umap_mod2")
# 
# #rm(list=ls())
# #load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Cancer gene expression wrangling/Cancer LR gene expression umap.RData")
# #load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Cancer gene expression wrangling/Cancer summary gene expression.RData")
# #load( "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/LigandReceptor/Filtered_Human-2015-Ramilowski-LR-pairs.RData")
# #LRgenelist <- unique( c(LRpairsFiltered$HPMR.Receptor, LRpairsFiltered$HPMR.Ligand) )
# LRumapgenes <- data.table(gene_summary[GeneName%in%LRgenelist])[Mean>0][coverage>0.05]$GeneName
# #names(LRumapgenes ) <- gsub("-","_", names(LRumapgenes))
# #names(u_dat ) <- gsub("-","_", names(u_dat))
# #save(LRumapgenes,file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Cancer gene expression wrangling/Cancer LR gene list.RData")
# load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Cancer gene expression wrangling/Cancer LR gene list.RData")
# 
# TrimUmapdd<-function(){
#   load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Cancer gene expression wrangling/Cancer LR gene expression umap.RData")
#   u_dat%>%dplyr::select(c("Cell.ID","Celltype","Celltype_subtype","Sample","orig.ident","Day",
#                           "Patient.Study.ID","ARM","Treatment","dynamic_class","dynamic_class3",
#                           "key_","Disc_V1","Disc_V2","Disc_V3","Disc_V4","Disc_V5","PhenoCelltype" ,
#                           LRumapgenes,"V1","V2","V3","V4"))
# }
# u_datTrimmed <- TrimUmapdd()
# save(u_datTrimmed,file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Cancer gene expression wrangling/Cancer LR gene expression umap Trimmed.RData")
# names(u_datTrimmed ) <- gsub("-","_", names(u_datTrimmed))
# 
# 
# ggplot( u_datTrimmed[ARM!="A"][5*(1:1000),],aes(V1,V2,col= log(1+TGFBR1) ) )+  geom_point(size=0.5)+theme_classic()#+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# ggplot( u_datTrimmed[ARM!="A"][V1<5&V2<5& V1> -5&V2>-5],aes(V1,V2,col= dynamic_class3 ) )+  geom_point(size=0.5)+theme_classic()+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# 
# 
# u_dat <- u_datTrimmed
# rm(list="u_datTrimmed")
# 
# which_compare <- which(abs(u_dat$V1)<10&abs(u_dat$V2)<10) #&abs(u_dat$V3)<5&abs(u_dat$V4)<5)
# corgenes <- data.table(t(
#   cor(u_dat[which_compare,]%>%dplyr::select(V1:V4), log(1+u_dat[which_compare,]%>%dplyr::select(gsub("-","_", LRumapgenes) )) ,method="spearman")
# ),keep.rownames = T)
# corgenes[order(-abs(V1))][1:10]
# corgenes[order(-abs(V2))][1:10]
# 
# corgenes[order(-abs(V3))][1:10]
# corgenes[order(-abs(V4))][1:10]
# 
# 
# ggplot( u_dat[which_compare,][ARM!="A"],aes(x=dynamic_class3, V2,y= log(1+TGFB1) ) )+  geom_violin()
# ggplot( u_dat[which_compare,][ARM!="A"],aes(x=dynamic_class3, V2,y= log(1+CD82) ) )+  geom_violin()
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V1,V2,col= log(1+CD82) ) )+  geom_point(size=0.5)
# 
# 
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V1,V2,col= dynamic_class3 ) )+  geom_point(size=0.5)+theme_classic()+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V1,V2,col= log(1+EFNA5) ) )+  geom_point(size=0.5)+theme_classic()+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V1,V2,col= log(1+SEMA6D) ) )+  geom_point(size=0.5)
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V1,V2,col= log(1+TNFSF10) ) )+  geom_point(size=0.5)
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V1,V2,col= log(1+NTNG1) ) )+  geom_point(size=0.5)
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V1,V2,col= log(1+CD44) ) )+  geom_point(size=0.5)
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V1,V2,col= log(1+ERBB4) ) )+  geom_point(size=0.5)
# 
# 
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V1,V2,col= log(1+IL4R) ) )+  geom_point(size=0.5)
# 
# 
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V1,V2,col= Patient.Study.ID ) )+  geom_point(size=0.5)+theme_classic()+theme(aspect.ratio=1)+facet_wrap(~Patient.Study.ID)
# 
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V1,V2,col= log(1+EGF) ) )+  geom_point(size=0.5)
# 
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V1,V2,col= log(1+EFNB2) ) )+  geom_point(size=0.5)
# 
# 
# 
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V1,V2,col= log(1+BMPR1B) ) )+  geom_point(size=0.5)
# 
# 
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V1,V2,col= log(1+TNFSF10) ) )+  geom_point(size=0.5)+theme_classic()+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V1,V2,col= log(1+SEMA6D) ) )+  geom_point(size=0.5)+theme_classic()+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V3,V2,col= Patient.Study.ID ) )+  geom_point(size=0.5)+theme_classic()+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# 
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V1,V2,col= log(1+LRP2) ) )+  geom_point(size=0.5)+theme_classic()+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# 
# 
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V1,V2,col= log(1+TNFSF10) ) )+  geom_point(size=0.5)+theme_classic()+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# 
# 
# 
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V1,V2,col= log(1+LRP2) ) )+  geom_point(size=0.5)+theme_classic()#+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V1,V2,col= log(1+LRP2) ) )+  geom_point(size=0.5)+theme_classic()#+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# 
# 
# 
# ggplot( u_dat,aes(V3,V4,col= ARM ) )+  geom_point(size=0.5)+theme_classic()+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# ggplot( u_dat,aes(V2,V4,col= Patient.Study.ID ) )+  geom_point(size=0.5)+theme_classic()+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# 
# 
# 
# 
# cor(umap_mod2$layout)
# pairs(umap_mod2$layout)
# 
# 
# 
# 
# unique(LRpairsFiltered$HPMR.Receptor)
# 
# ggplot( u_dat[which_compare,][ARM!="A"][5*(1:1000),],aes(V1,V2,col= log(1+TGFBR1) ) )+  geom_point(size=0.5)+theme_classic()#+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# 
# 
# ggplot( u_dat[which_compare,][ARM!="A"][sample(1:10000,n=2000) ,],aes(x=dynamic_class3, y= log(1+EPHB1) ) )+  geom_violin()+theme_classic()#+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# 
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V1,V2,col= log(1+EGFR) ) )+  geom_point(size=0.5)+theme_classic()#+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V1,V2,col= log(1+TGFBR1) ) )+  geom_point(size=0.5)+theme_classic()#+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V1,V2,col= log(1+TGFBR3) ) )+  geom_point(size=0.5)+theme_classic()#+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# 
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V1,V2,col= log(1+COL1A2) ) )+  geom_point(size=0.5)+theme_classic()#+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V1,V2,col= log(1+COL1A1) ) )+  geom_point(size=0.5)+theme_classic()#+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V1,V2,col= log(1+COL11A1) ) )+  geom_point(size=0.5)+theme_classic()#+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V1,V2,col= log(1+FN1) ) )+  geom_point(size=0.5)+theme_classic()#+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V1,V2,col= log(1+SPARC) ) )+  geom_point(size=0.5)+theme_classic()#+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V1,V2,col= log(1+COL3A1) ) )+  geom_point(size=0.5)+theme_classic()#+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V1,V2,col= log(1+NFIA) ) )+  geom_point(size=0.5)+theme_classic()#+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V1,V2,col= log(1+LAMA2) ) )+  geom_point(size=0.5)+theme_classic()#+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V1,V2,col= log(1+COL5A1) ) )+  geom_point(size=0.5)+theme_classic()#+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# 
# #periostin binds to integrins on cancer cells, activating the Akt/PKB- and FAK-mediated signaling pathways.  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3205268/ 
# # Caf secreted https://www.nature.com/articles/s41419-018-1116-6
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V1,V2,col= log(1+POSTN) ) )+  geom_point(size=0.5)+theme_classic()#+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V1,V2,col= log(1+SULF1) ) )+  geom_point(size=0.5)+theme_classic()#+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V1,V2,col= log(1+KIF26B) ) )+  geom_point(size=0.5)+theme_classic()#+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V1,V2,col= log(1+NTS) ) )+  geom_point(size=0.5)+theme_classic()#+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V1,V2,col= log(1+ABI3BP) ) )+  geom_point(size=0.5)+theme_classic()#+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V1,V2,col= log(1+DCLK1) ) )+  geom_point(size=0.5)+theme_classic()#+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V1,V2,col= log(1+ITGA11) ) )+  geom_point(size=0.5)+theme_classic()#+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# 
# 
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V1,V3,col= log(1+ABI3BP) ) )+  geom_point(size=0.5)+theme_classic()#+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V1,V3,col= log(1+MIR99AHG) ) )+  geom_point(size=0.5)+theme_classic()#+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V1,V3,col= log(1+DCLK1) ) )+  geom_point(size=0.5)+theme_classic()#+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V1,V3,col= log(1+MMP11) ) )+  geom_point(size=0.5)+theme_classic()#+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V1,V3,col= log(1+ABCA10) ) )+  geom_point(size=0.5)+theme_classic()#+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V1,V3,col= log(1+ITGA11) ) )+  geom_point(size=0.5)+theme_classic()#+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V1,V3,col= log(1+KIF26B) ) )+  geom_point(size=0.5)+theme_classic()#+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# 
# 
# ggplot( u_dat[which_compare,][ARM!="A"],aes(V1,V3,col=  dynamic_class3)) +  geom_point(size=0.5)+scale_color_npg(name="Tumor response")+theme_classic()+ facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)+labs(y="umap B",x="umap A")
# ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Fibroblasts phenotypes AllArms/FibroblastsGeneUmap by Response and Day.png")
# #ggplot( u_dat[ARM!="A"],aes(V1,V3,col=  Celltype_subtype)) +  geom_point(size=0.5)+scale_color_npg(name="Cell type")+theme_classic()+theme(aspect.ratio=1)+labs(y="umap B",x="umap A")#+facet_wrap(~Celltype_subtype)
# #ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Fibroblasts phenotypes AllArms/FibroblastsGeneUmap by Celltype_subtype and Day.png")
# ggplot( u_dat[ARM!="A"],aes(V1,V3,col=  Patient.Study.ID)) +  geom_point(size=0.5)+theme_classic()+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)+labs(y="umap B",x="umap A")
# ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Fibroblasts phenotypes AllArms/FibroblastsGeneUmap by Patient.Study.ID and Day.png")
# 
# 
# #ggplot( u_dat[ARM!="A"],aes(V1,V3,col= log(1+S100A4) ) )+  geom_point(size=0.5)+theme_classic()#+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# 
# #CAF-S1 (== FAP high) fibroblasts are defined by extracellular matrix and inflammation signatures, 
# # cause increasing attraction, survival, and overall content of CD4+ CD25+ FOXP3+ regulatory T cells (Tregs) in the tumor microenvironment
# #CAF-S4 fibroblasts are characterized by a perivascular signature and do not contribute to immune suppression
# # Fibroblast activation protein (FAP) marker are associated with an immunosuppressive environment  https://cancerdiscovery.aacrjournals.org/content/10/9/1330
# ggplot( u_dat[ARM!="A"],aes(V1,V3,col= log(1+FAP) ) )+  geom_point(size=0.5)+theme_classic()#+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
# ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Fibroblasts phenotypes AllArms/FibroblastsGeneUmap overlay FAP.png")
# 
# # Differentiation marker : Myofibroblasts are differentiated fibroblasts that express the intracellular contractile protein alpha smooth muscle actin (α-SMA) https://www.sciencedirect.com/topics/neuroscience/myofibroblast ref 55.
# # Smooth muscle α-actin gene (SMA) = ACTA2 gene
# #allows fibroblasts to interact with the extracellular matrix through cell membrane integrins and influence matrix organization and contracture
# ggplot( u_dat[ARM!="A"],aes(V1,V3,col= log(1+ACTA2) ) )+  geom_point(size=0.5)+theme_classic()
# 
# #integrin β1 (CD29=ITGB1 gene)
# ggplot( u_dat[ARM!="A"],aes(V1,V3,col= log(1+ITGB1) ) )+  geom_point(size=1)+theme_classic()
# ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Fibroblasts phenotypes AllArms/FibroblastsGeneUmap overlay ITGB1.png")
# 
# 
# 
# 
# corgenes <- data.table(t(
#   cor(u_dat[V1> -1]%>%dplyr::select(V1:V3), log(1+umap_in2[which(u_dat$V1> -1),]) ,method="spearman")
# ),keep.rownames = T)
# corgenes[order(-abs(V1))][1:10]
# corgenes[order(-abs(V2))][1:10]
# corgenes[order(-abs(V3))][1:10]
# 
# ggplot( u_dat[ARM!="A"],aes(V1,V2,col= log(1+COL14A1) ) )+  geom_point(size=1)+theme_classic()
# 
# ggplot( u_dat[ARM!="A"],aes(V1,V3,col= log(1+NTS) ) )+  geom_point(size=1)+theme_classic()
# 
# ggplot( u_dat[ARM!="A"],aes(V1,V3,col= log(1+COL11A1) ) )+  geom_point(size=1)+theme_classic()
# ggplot( u_dat[ARM!="A"],aes(V1,V3,col= log(1+FN1) ) )+  geom_point(size=1)+theme_classic()
# ggplot( u_dat[ARM!="A"],aes(V1,V3,col= log(1+COL10A1) ) )+  geom_point(size=1)+theme_classic()
# 
# ggplot( u_dat[ARM!="A"],aes(V1,V3,col= log(1+ABI3BP) ) )+  geom_point(size=1)+theme_classic()
# 
# #Undifferentiated cell marker: In pancreatic cancer, the expression of HEY1 was found to be higher in the DCLK1HI T-IC subpopulation than in its more differentiated counterpart (Bailey et al., 2014).
# ggplot( u_dat[ARM!="A"],aes(V1,V3,col= log(1+DCLK1) ) )+  geom_point(size=1)+theme_classic()
# 
# 
# 
# 
# # ADAM12 = marker of CAF https://www.nature.com/articles/s41416-019-0509-3
# ggplot( u_dat[ARM!="A"],aes(V1,V3,col= log(1+ADAM12) ) )+  geom_point(size=1)+theme_classic()
# ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Fibroblasts phenotypes AllArms/FibroblastsGeneUmap overlay ADAM12.png")
# 
# # Complex multiple effects, but controls antigen expression and processing and inflammatory response (https://molecular-cancer.biomedcentral.com/articles/10.1186/s12943-020-01290-7). can block
# ggplot( u_dat[ARM!="A"],aes(V1,V3,col= log(1+SPON1) ) )+  geom_point(size=1)+theme_classic()
# ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Fibroblasts phenotypes AllArms/FibroblastsGeneUmap overlay SPON1.png")
# 
# # Wnt signal enhancers and stem cell growth factors http://genesdev.cshlp.org/content/28/4/305.full.html
# ggplot( u_dat[ARM!="A"],aes(V1,V3,col= log(1+LGR4) ) )+  geom_point(size=1)+theme_classic()
# filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Fibroblasts phenotypes AllArms/FibroblastsGeneUmap overlay LGR4.png"
# 
# # CAF-S1 marker gene https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7690906/
# ggplot( u_dat[ARM!="A"],aes(V1,V3,col= log(1+SEMA3C) ) )+  geom_point(size=1)+theme_classic()
# filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Fibroblasts phenotypes AllArms/FibroblastsGeneUmap overlay SEMA3C.png"
# 
# 
# ggplot( u_dat[ARM!="A"],aes(V1,V3,col= log(1+ROR1) ) )+  geom_point(size=1)+theme_classic()
# ggplot( u_dat[ARM!="A"],aes(V1,V3,col= log(1+FGFR2) ) )+  geom_point(size=1)+theme_classic()
# ggplot( u_dat[ARM!="A"],aes(V1,V3,col= log(1+EGFR) ) )+  geom_point(size=1)+theme_classic()
# ggplot( u_dat[ARM!="A"],aes(V1,V3,col= log(1+ERBB4) ) )+  geom_point(size=1)+theme_classic()
# ggplot( u_dat[ARM!="A"],aes(V1,V3,col= log(1+CXCL12) ) )+  geom_point(size=1)+theme_classic()
# 
# 
# 
# 
# 
# 
# #ggplot( u_dat[ARM!="A"],aes(V1,V3,col= log(1+MEG8) ) )+  geom_point(size=1)+theme_classic()
# ggplot( u_dat[ARM!="A"],aes(V1,V3,col= log(1+SEMA3C) ) )+  geom_point(size=1)+theme_classic()
# ggplot( u_dat[ARM!="A"],aes(V1,V3,col= log(1+ENAH) ) )+  geom_point(size=1)+theme_classic()
# ggplot( u_dat[ARM!="A"],aes(V1,V3,col= log(1+ROR1) ) )+  geom_point(size=1)+theme_classic()
# ggplot( u_dat[ARM!="A"],aes(V1,V3,col= log(1+ARHGAP28) ) )+  geom_point(size=1)+theme_classic()
# ggplot( u_dat[ARM!="A"],aes(V1,V3,col= log(1+OSBPL3) ) )+  geom_point(size=1)+theme_classic()
# ggplot( u_dat[ARM!="A"],aes(V1,V3,col= log(1+KCNQ1OT1) ) )+  geom_point(size=1)+theme_classic()
# 
# 
# 
# load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/LigandReceptor/Interesting Fibroblast Cancer LR pairs.RData")
# #listGF_LR
# 
# 
# unname(unlist( listGF_LR ))
# 
# 
# umap_inGF  <- fulldd %>% dplyr::select( unique(unname(unlist( listGF_LR$LigandGene )) ) )
# gene_summaryGF <- data.table( GeneName= colnames(umap_inGF),Mean= colMeans(umap_inGF), Variance= colVars( as.matrix(umap_inGF) ),coverage= colMeans(umap_inGF>0), nCellsExpressed= colSums(umap_inGF>0) )
# gene_summaryGF[Mean>0][coverage>0.05]$Mean %>% log()%>%hist()
# gene_summaryGF[Mean>0][coverage>0.05]$Variance %>% log()%>%hist()
# 
# umap_inGF2 <-  umap_inGF[,gene_summaryGF[Mean>0][coverage>0.001]$GeneName ,with=FALSE ]
# umap_modGF2 <- umap::umap(log(1+umap_inGF2),n_components=2)
# u_datGF <- cbind(fulldd,umap_modGF2$layout)
# #umap_mod <- umap::umap(umap_in,n_components=2)
# cor(umap_modGF2$layout)
# pairs(umap_modGF2$layout[1:100,])
# #u_dat <- cbind(fulldd,umap_mod2$layout)
# 
# names(u_datGF ) <- gsub("-","_", names(u_datGF))
# 
# gf_l_dist <- dist( t(scale(log(1 + umap_inGF2 ) , center=T)))
# #gg_gf_ldist <- gather(as.data.table(as.matrix(gf_l_dist),keep.rownames = T),var,val,-rn)
# #ggplot(gg_gf_ldist,aes(x=rn ,y= var ,fill= val))+geom_tile()
# pheatmap::pheatmap(as.matrix(gf_l_dist), cutree_rows=2, cutree_cols=2)
# 
# 
# ggplot( u_datGF[ARM!="A"][V1> -10 & V2<10 & V1< 10 & V2> -10],aes(V1,V2,col=  dynamic_class3)) +  geom_point(size=1)+scale_color_npg(name="Tumor response")+theme_classic()+facet_wrap(~Day)
# ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Fibroblast phenotypes AllArms/FibroblastGrowthFactorGeneUmap by Response and Day.png")
# ggplot( u_datGF[ARM!="A"],aes(V1,V2,col=  Celltype_subtype)) +  geom_point(size=3)+scale_color_npg(name="Cell type")+theme_classic()#+facet_wrap(~Celltype_subtype)
# ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Fibroblast phenotypes AllArms/FibroblastGrowthFactorGeneUmap by Celltype_subtype and Day.png")
# ggplot( u_datGF[ARM!="A"],aes(V1,V2,col=  Patient.Study.ID)) +  geom_point(size=3)+theme_classic()
# ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Fibroblast phenotypes AllArms/FibroblastGrowthFactorGeneUmap by Patient.Study.ID and Day.png")
# 
# corgenesGF <- data.table(t(
#   cor(u_datGF%>%dplyr::select(V1:V2), log(1+umap_inGF2) ,method="spearman")
# ),keep.rownames = T)
# corgenesGF[order(-abs(V1))][1:10]
# corgenesGF[order(-abs(V2))][1:10]
# 
# 
# #
# ggplot( u_dat[ARM!="A"],aes(V1,V3,col= dynamic_class3 ) )+  geom_point(size=1)+theme_classic()+facet_wrap(dynamic_class3~Day)
# 
# 
# 
# ggplot( u_dat[which_compare,][ARM!="A"],aes(y=V1,x=log(1+Day),group=interaction(Day,dynamic_class3),fill=  dynamic_class3,col=  dynamic_class3)) +  
#   geom_violin()+scale_color_npg(name="Tumor response")+scale_fill_npg(name="Tumor response")+theme_classic()+scale_x_continuous(breaks=log(c(0,14,180)), labels=c(0,14,180))+
#   labs(y="Umap 1",x="Day")+theme(aspect.ratio=1)
# 
# ggplot( u_dat[which_compare,][ARM!="A"],aes(y=V2,x=log(1+Day),group=interaction(Day,dynamic_class3),fill=  dynamic_class3,col=  dynamic_class3)) +  
#   geom_violin()+scale_color_npg(name="Tumor response")+scale_fill_npg(name="Tumor response")+theme_classic()+scale_x_continuous(breaks=log(c(0,14,180)), labels=c(0,14,180))+
#   labs(y="Umap 2",x="Day")
# 
# ggplot( u_dat[which_compare,][ARM!="A"],aes(y=V3,x=log(1+Day),group=interaction(Day,dynamic_class3),fill=  dynamic_class3,col=  dynamic_class3)) +  
#   geom_violin()+scale_color_npg(name="Tumor response")+scale_fill_npg(name="Tumor response")+theme_classic()+scale_x_continuous(breaks=log(c(0,14,180)), labels=c(0,14,180))+
#   labs(y="Umap 3",x="Day")
# 
# ggplot( u_dat[which_compare,][ARM!="A"],aes(y=V4,x=log(1+Day),group=interaction(Day,dynamic_class3),fill=  dynamic_class3,col=  dynamic_class3)) +  
#   geom_violin()+scale_color_npg(name="Tumor response")+scale_fill_npg(name="Tumor response")+theme_classic()+scale_x_continuous(breaks=log(c(0,14,180)), labels=c(0,14,180))+
#   labs(y="Umap 4",x="Day")
# 
# 
# 
# ggplot( u_dat[ARM!="A"],aes(V1,V4,col=  dynamic_class3)) +  geom_point(size=2)+scale_color_npg(name="Tumor response")+theme_classic()#+facet_wrap(~Day)
# ggplot( u_dat[ARM!="A"],aes(V1,V4,col=  Patient.Study.ID)) +  geom_point(size=1)+theme_classic()+facet_wrap(dynamic_class3~Day)
# 
# methodcor= "spearman"
# NN <- 4
# corVallSub <- cor(u_dat[ARM!="A"][Celltype_subtype!="DC"]%>%dplyr::select(paste0("V",1:NN)),u_dat[ARM!="A"][Celltype_subtype!="DC"]%>%dplyr::select(common_ssgsea),method=methodcor)#corV1<-cor(u_dat$V1,full_dd_transp)
# 
# 
# corrWithUmap <- data.table( t(corVall) , keep.rownames = T)
# corrWithUmap[order(-abs(V1))][1:30]
# corrWithUmap[order(-abs(V2))][1:10]
# corrWithUmap[order(-abs(V3))][1:10]
# corrWithUmap[order(-abs(V4))][1:10]
# 
# corrWithUmapSub <- data.table( t(corVallSub) , keep.rownames = T)
# corrWithUmapSub[order(-abs(V1))][1:30]
# corrWithUmapSub[order(-abs(V2))][1:10]
# corrWithUmapSub[order(-abs(V3))][1:10]
# corrWithUmapSub[order(-abs(V4))][1:10]
# 
# u_dat[,MyloidStatus:= "Macrophage"]
# u_dat[Celltype_subtype=="DC",MyloidStatus:= "Dendritic cell"]
# 
# ggplot( u_dat[ARM!="A"],aes(V1,V4,shape=  MyloidStatus, col=  MyloidStatus)) +  geom_point(size=2)+theme_classic() + theme(aspect.ratio=1)+labs(y="Umap B" ,x="Umap A") +scale_color_aaas(name="Cell type") + scale_shape_discrete(name="Cell type")
# ggplot( u_dat[ARM!="A"],aes(V1,V4,shape=  MyloidStatus, col=dynamic_class3)) +  geom_point(size=2)+theme_classic() + theme(aspect.ratio=1)+labs(y="Umap B" ,x="Umap A") +scale_color_npg(name="Tumor response", labels=c("Resistant","Sensitive")) + scale_shape_discrete(name="Cell type")
# ggplot( u_dat[ARM!="A"],aes(V1,V4,shape=  MyloidStatus, col=dynamic_class3)) +  geom_point(size=1)+theme_classic() + theme(aspect.ratio=1)+labs(y="Umap B" ,x="Umap A") +scale_color_npg(name="Tumor response", labels=c("Resistant","Sensitive")) + scale_shape_discrete(name="Cell type")+facet_wrap(~Day)
# 
# ggplot( u_dat[ARM!="A"],aes(V1,V4,shape=  MyloidStatus,col=dynamic_class3)) +  geom_point(size=2)+theme_classic()+facet_wrap(~Day)
# ggplot( u_dat[ARM!="A"],aes(V1,V4,shape=  MyloidStatus,col= Patient.Study.ID)) +  geom_point(size=2)+theme_classic()+facet_wrap(~Day)
# 
# ggplot( u_dat[ARM!="A"],aes(V1,V4,col=  ZHANG_RESPONSE_TO_IKK_INHIBITOR_AND_TNF_UP)) +  geom_point(size=2)+theme_classic()
# 
# ggplot( u_dat[ARM!="A"],aes(V1,V4,col=  FULCHER_INFLAMMATORY_RESPONSE_LECTIN_VS_LPS_DN)) +  geom_point(size=2)+theme_classic()+ theme(aspect.ratio=1)+labs(y="Umap B" ,x="Umap A")
# ggplot( u_dat[ARM!="A"],aes(V1,V4,col=  GAURNIER_PSMD4_TARGETS)) +  geom_point(size=2)+theme_classic()+ theme(aspect.ratio=1)+labs(y="Umap B" ,x="Umap A")
# ggplot( u_dat[ARM!="A"],aes(V1,V4,col=  LINDSTEDT_DENDRITIC_CELL_MATURATION_B)) +  geom_point(size=2)+theme_classic()+ theme(aspect.ratio=1)+labs(y="Umap B" ,x="Umap A")
# ggplot( u_dat[ARM!="A"],aes(V1,V4,col=  RUTELLA_RESPONSE_TO_HGF_VS_CSF2RB_AND_IL4_DN)) +  geom_point(size=2)+theme_classic()+ theme(aspect.ratio=1)+labs(y="Umap B" ,x="Umap A")
# ggplot( u_dat[ARM!="A"],aes(V1,V4,col=  DIRMEIER_LMP1_RESPONSE_EARLY)) +  geom_point(size=2)+theme_classic()+ theme(aspect.ratio=1)+labs(y="Umap B" ,x="Umap A")
# ggplot( u_dat[ARM!="A"],aes(V1,V4,col=  BASSO_CD40_SIGNALING_UP)) +  geom_point(size=2)+theme_classic()+ theme(aspect.ratio=1)+labs(y="Umap B" ,x="Umap A")
# ggplot( u_dat[ARM!="A"],aes(V1,V4,col=  TIAN_TNF_SIGNALING_VIA_NFKB)) +  geom_point(size=2)+theme_classic()+ theme(aspect.ratio=1)+labs(y="Umap B" ,x="Umap A")
# ggplot( u_dat[ARM!="A"],aes(V1,V4,col=  RASHI_NFKB1_TARGETS)) +  geom_point(size=2)+theme_classic()+ theme(aspect.ratio=1)+labs(y="Umap B" ,x="Umap A")
# ggplot( u_dat[ARM!="A"],aes(V1,V4,col=  REACTOME_DOWNSTREAM_TCR_SIGNALING)) +  geom_point(size=2)+theme_classic()+ theme(aspect.ratio=1)+labs(y="Umap B" ,x="Umap A")
# ggplot( u_dat[ARM!="A"],aes(V1,V4,col=  PHONG_TNF_TARGETS_UP)) +  geom_point(size=2)+theme_classic()+ theme(aspect.ratio=1)+labs(y="Umap B" ,x="Umap A")
# 
# ggplot( u_dat[ARM!="A"],aes(V1,V4,col=  KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION)) +  geom_point(size=2)+theme_classic()
# ggplot( u_dat[ARM!="A"],aes(V1,V4,col=  PHONG_TNF_TARGETS_UP)) +  geom_point(size=2)+theme_classic()
# ggplot( u_dat[ARM!="A"],aes(V1,V4,col=  BIOCARTA_CD40_PATHWAY)) +  geom_point(size=2)+theme_classic()
# ggplot( u_dat[ARM!="A"],aes(V1,V4,col=  BIOCARTA_TNFR2_PATHWAY)) +  geom_point(size=2)+theme_classic()
# 
# 
# ggplot( u_dat[ARM!="A"],aes(V1,V4,col=  REACTOME_PHOSPHORYLATION_OF_CD3_AND_TCR_ZETA_CHAINS)) +  geom_point(size=2)+theme_classic()
# ggplot( u_dat[ARM!="A"],aes(V1,V4,col=  BIOCARTA_TH1TH2_PATHWAY)) +  geom_point(size=2)+theme_classic()
# ggplot( u_dat[ARM!="A"],aes(V1,V4,col=  REACTOME_PD1_SIGNALING)) +  geom_point(size=2)+theme_classic()
# ggplot( u_dat[ARM!="A"],aes(V1,V4,col=  GAURNIER_PSMD4_TARGETS)) +  geom_point(size=2)+theme_classic()
# ggplot( u_dat[ARM!="A"],aes(V1,V4,col=  EINAV_INTERFERON_SIGNATURE_IN_CANCER)) +  geom_point(size=2)+theme_classic()
# ggplot( u_dat[ARM!="A"],aes(V1,V4,col=  ZHANG_INTERFERON_RESPONSE)) +  geom_point(size=2)+theme_classic()
# 
# 
# 
# ggplot( u_dat[ARM!="A"],aes(V1,V4,col=  FOSTER_TOLERANT_MACROPHAGE_DN)) +  geom_point(size=2)+theme_classic()
# ggplot( u_dat[ARM!="A"],aes(V1,V4,col=  REACTOME_INTERFERON_SIGNALING)) +  geom_point(size=2)+theme_classic()
# 
# 
# 
# ## Next step is to overlay genes
# 
# 
# 
# load("~/Dropbox/Cancer_pheno_evo/data/FELINE2/PhenotypesAllArms/T_cellsCD4+ T cells.RData")
# 
# 
# #u_dat[,MyloidStatus:= "Macrophage"]
# #u_dat[Celltype_subtype=="DC",MyloidStatus:= "Dendritic cell"]
# 
# # set.seed(123); umod2 <- umap::umap(u_dat%>% dplyr::select(V1:V5) , n_components=2,n_neighbors=40)
# # #umod2$layout%>%plot
# # u_dat2<- cbind(u_dat,VV=umod2$layout)
# # ggplot( u_dat2[ARM!="A"],aes(VV.V1,VV.V2,col=  Celltype_subtype,shape=Celltype_subtype )) +  geom_point(size=2)+scale_color_npg(name="Tumor response")+theme_classic()
# # ggplot( u_dat2[ARM!="A"],aes(VV.V1,VV.V2,col=  Celltype_subtype,shape=Celltype_subtype )) +  geom_point(size=2)+scale_color_npg(name="Tumor response")+theme_classic()
# # ggplot( u_dat2[ARM!="A"],aes(VV.V1,VV.V2,col=  dynamic_class3,shape=Celltype_subtype )) +  geom_point(size=2)+scale_color_npg(name="Tumor response")+theme_classic()+facet_wrap(~Day)
# # 
# # ggplot( u_dat2[ARM!="A"],aes(VV.V2,VV.V1,col=  dynamic_class3,shape=Celltype_subtype )) +  geom_point(size=2)+scale_color_npg(name="Tumor response")+theme_classic()+facet_wrap(dynamic_class3~Day)
# 
# ggplot( u_dat[ARM!="A"],aes(V3,V1,col=  dynamic_class3)) +  geom_point(size=2)+scale_color_npg(name="Tumor response")+theme_classic()+facet_wrap(~Day)
# ggplot( u_dat[ARM!="A"],aes(V1,V2,col=  V3)) +  geom_point(size=2)+theme_classic()+facet_wrap(dynamic_class3~Day)
# #ggplot( u_dat,aes(V1,V2,col=  Celltype_subtype)) +  geom_point(size=2)+scale_color_npg(name="Tumor response")+theme_classic()+facet_wrap(~Day)
# 
# 
# ggplot( u_dat[ARM!="A"],aes(y=V1,x=log(1+Day),group=interaction(Day,dynamic_class3),fill=  dynamic_class3,col=  dynamic_class3)) +  
#   geom_violin()+scale_color_npg(name="Tumor response")+scale_fill_npg(name="Tumor response")+theme_classic()+scale_x_continuous(breaks=log(c(0,14,180)), labels=c(0,14,180))+
#   labs(y="Umap 1",x="Day")+theme(aspect.ratio=1)
# 
# ggplot( u_dat[ARM!="A"]%>%group_by(Patient.Study.ID,Day,dynamic_class3) %>% dplyr::summarise(V2=mean(V2)), aes(y=V2,x=log(1+Day),group=interaction(Day,dynamic_class3),fill=  dynamic_class3,col=  dynamic_class3)) +  
#   geom_violin()+scale_color_npg(name="Tumor response")+scale_fill_npg(name="Tumor response")+theme_classic()+scale_x_continuous(breaks=log(c(0,14,180)), labels=c(0,14,180))+
#   labs(y="Umap 2",x="Day")
# 
# ggplot( u_dat[ARM!="A"],aes(y=V3,x=log(1+Day),group=interaction(Day,dynamic_class3),fill=  dynamic_class3,col=  dynamic_class3)) +  
#   geom_violin()+scale_color_npg(name="Tumor response")+scale_fill_npg(name="Tumor response")+theme_classic()+scale_x_continuous(breaks=log(c(0,14,180)), labels=c(0,14,180))+
#   labs(y="Umap 3",x="Day")
# 
# ggplot( u_dat[ARM!="A"],aes(y=V4,x=log(1+Day),group=interaction(Day,dynamic_class3),fill=  dynamic_class3,col=  dynamic_class3)) +  
#   geom_violin()+scale_color_npg(name="Tumor response")+scale_fill_npg(name="Tumor response")+theme_classic()+scale_x_continuous(breaks=log(c(0,14,180)), labels=c(0,14,180))+
#   labs(y="Umap 4",x="Day")
# 
# 
# ggplot( u_dat[ARM!="A"],aes(y=V5,x=log(1+Day),group=interaction(Day,dynamic_class3),fill=  dynamic_class3,col=  dynamic_class3)) +  
#   geom_violin()+scale_color_npg(name="Tumor response")+scale_fill_npg(name="Tumor response")+theme_classic()+scale_x_continuous(breaks=log(c(0,14,180)), labels=c(0,14,180))+
#   labs(y="Umap 5",x="Day")
# 
# 
# 
# corrWithUmap <- data.table( t(corVall) , keep.rownames = T)
# corrWithUmap[order(-abs(V1))][1:10]
# corrWithUmap[order(-abs(V2))][1:10]
# corrWithUmap[order(-abs(V3))][1:10]
# corrWithUmap[order(-abs(V4))][1:10]
# corrWithUmap[order(-abs(V4))][1:10]
# 
# 
# ggplot( u_dat,aes(V1,V2,col=  GSE41978_KLRG1_HIGH_VS_LOW_EFFECTOR_CD8_TCELL_UP)) +  geom_point(size=2)+theme_classic()#+facet_wrap(~Day)
# 
# ggplot( u_dat[ARM!="A"],aes(V4,V5,col=  dynamic_class3)) +  geom_point(size=2)+scale_color_npg(name="Tumor response")+theme_classic()#+facet_wrap(~Day)
# ggplot( u_dat[ARM!="A"],aes(V3,V4,col=  Patient.Study.ID)) +  geom_point(size=1)+theme_classic()+facet_wrap(dynamic_class3~Day)
# 
# 
# 
# 
# 
# # 
# #  require(ouija);require( viridis)
# #  example_gex <- u_dat[ARM!="A"] %>%dplyr::select(V1:V5)
# #  cols <- paste0("V",1:5)
# #  example_gex[, (cols) := lapply(.SD, scale), .SDcols=cols]
# #  example_gex[, V1 := log(exp(V1)+1)]
# #  example_gex[, V2 := log(exp(V2)+1)]
# #  example_gex[, V3 := log(exp(V3)+1)]
# #  example_gex[, V4 := log(exp(V4)+1)]
# #  example_gex[, V5 := log(exp(V5)+1)]
# # 
# #  oui <- ouija(as.matrix( example_gex) , response_type="switch" , iter = 1000)#[sample(seq_len(nrow(example_gex)), 200), ]
# # # plot_diagnostics(oui)
# # # plot_expression(oui)
# # tmap <- map_pseudotime(oui) 
# # u_dat2 <- u_dat[ARM!="A"]
# # u_dat2[,pseudotime:= tmap]
# # u_dat2[, (cols) := lapply(.SD, scale), .SDcols=cols]
# # #u_dat2[, V1 := log(exp(V1)+1)]
# # #u_dat2[, V2 := log(exp(V2)+1)]
# # ggplot( u_dat2 ,aes(V1,V2,col=pseudotime,shape=dynamic_class3)) + geom_point()+scale_color_viridis()+facet_wrap(~dynamic_class3)+theme_classic()+theme(aspect.ratio=1)
# # ggplot( u_dat2 ,aes(V1,V2,col=Celltype_subtype,shape=dynamic_class3)) + geom_point()+facet_wrap(~dynamic_class3)
# # ggplot( u_dat2 ,aes(x=pseudotime,col=dynamic_class3,fill=dynamic_class3)) + geom_density(alpha=0.3)+scale_color_npg()+theme_classic()+theme(aspect.ratio=1)#+facet_wrap(~dynamic_class3)
# # 
# # 
# # ggplot( u_dat2 ,aes(V1,V2,col=dynamic_class3,shape=dynamic_class3)) + geom_point()+scale_color_npg()+facet_wrap(~Day)+theme_classic()+theme(aspect.ratio=1)
# 
# # ggplot( u_dat2 ,aes(log(exp(V1)-1),log(exp(V2)-1),col=Patient.Study.ID,shape=dynamic_class3)) + geom_point()+facet_wrap(~dynamic_class3)+theme_classic()+theme(legend.position="none",aspect.ratio=1)
# # ggplot( u_dat2 ,aes(log(exp(V1)-1),log(exp(V2)-1),col=pseudotime,shape=dynamic_class3)) + geom_point()+scale_color_viridis()+facet_wrap(~dynamic_class3)+theme_classic()+theme(aspect.ratio=1)
# # ggplot( u_dat2 ,aes(log(exp(V1)-1),log(exp(V2)-1),col=GSE5589_IL6_KO_VS_IL10_KO_LPS_STIM_MACROPHAGE_180MIN_UP,shape=dynamic_class3)) + geom_point()+scale_color_viridis()+facet_wrap(~dynamic_class3)+theme_classic()+theme(aspect.ratio=1)
# # #good
# # ggplot( u_dat2[Platform!="ICELL8"] ,aes(V1,V2,col=GSE9988_ANTI_TREM1_VS_ANTI_TREM1_AND_LPS_MONOCYTE_UP)) + geom_point()+scale_color_viridis() #ggplot( u_dat[Platform!="ICELL8"] ,aes(PC1,PC2,col=Infercnv_CNA)) + geom_point()+scale_color_jco()
# # ggplot( u_dat2[Platform!="ICELL8"] ,aes(V1,V2,col=GSE10856_CTRL_VS_TNFRSF6B_IN_MACROPHAGE_DN)) + geom_point()+scale_color_viridis() #ggplot( u_dat[Platform!="ICELL8"] ,aes(PC1,PC2,col=Infercnv_CNA)) + geom_point()+scale_color_jco()
# # ## Macrophage polarization : (DcR3, TLR4) M2 - M1 
# # ggplot( u_dat2[Platform!="ICELL8"] ,aes(V1,V2,col=GSE9509_LPS_VS_LPS_AND_IL10_STIM_IL10_KO_MACROPHAGE_20MIN_UP)) + geom_point()+scale_color_viridis() #ggplot( u_dat[Platform!="ICELL8"] ,aes(PC1,PC2,col=Infercnv_CNA)) + geom_point()+scale_color_jco()
# # 
# # 
# # ggplot( u_dat2 ,aes(x=pseudotime,col=dynamic_class3,fill=dynamic_class3,group=Patient.Study.ID)) + geom_density(alpha=0.3)+scale_color_npg()+theme_classic()+theme(aspect.ratio=1)+facet_wrap(~dynamic_class3)
# # ggplot( u_dat2 ,aes(x=(pseudotime),col=dynamic_class3,fill=dynamic_class3,group=Patient.Study.ID)) + geom_density(alpha=0.3)+scale_color_npg()+theme_classic()+theme(aspect.ratio=1)+facet_wrap(~dynamic_class3)
# # 
# # lmer1 <- lmer((pseudotime)~dynamic_class3+(1|Patient.Study.ID),data=u_dat2)
# # summary(lmer1)
# # example_gex <- u_dat[Platform!="ICELL8"][Celltype=="T cells"] %>%dplyr::select(V1,V2)
# # cols <- paste0("V",1:2)
# # example_gex[, (cols) := lapply(.SD, scale), .SDcols=cols]
# # example_gex[, V1 := log(exp(V1)+1)]
# # example_gex[, V2 := log(exp(V2)+1)]
# # 
# 
# 
# 
# umod2 <- umap::umap(u_dat%>% dplyr::select(V1:V5))
# umod2$layout%>%plot
# ggplot( u_dat[ARM!="A"],aes(V1,V4,col=  NAGASHIMA_EGF_SIGNALING_UP)) +  geom_point(size=2)+theme_classic()
# 
# ggplot( u_dat[ARM!="A"],aes(V1,V2,col=  TIAN_TNF_SIGNALING_VIA_NFKB)) +  geom_point(size=2)+theme_classic()+facet_wrap(dynamic_class3~Day)
# 
# 
# 
# ggplot( u_dat[ARM!="A"],aes(y=DAUER_STAT3_TARGETS_DN,x=log(1+Day),group=interaction(Day,dynamic_class3),fill=  dynamic_class3,col=  dynamic_class3)) +  
#   geom_violin()+scale_color_npg(name="Tumor response")+scale_fill_npg(name="Tumor response")+theme_classic()+scale_x_continuous(breaks=log(c(0,14,180)), labels=c(0,14,180))+
#   labs(y="Umap 1",x="Day")
# ggplot( u_dat[ARM!="A"],aes(V1,V2,col=  GSE21546_WT_VS_SAP1A_KO_DP_THYMOCYTES_UP)) +  geom_point(size=2)+theme_classic()
# 
# ggplot( u_dat[ARM!="A"],aes(V1,V2,col=  GSE29618_MONOCYTE_VS_MDC_UP)) +  geom_point(size=2)+theme_classic()
# ggplot( u_dat[ARM!="A"],aes(V1,V2,col=  GSE29618_MONOCYTE_VS_MDC_DAY7_FLU_VACCINE_UP)) +  geom_point(size=2)+theme_classic()
# ggplot( u_dat[ARM!="A"],aes(V1,V2,col=  GSE44732_UNSTIM_VS_IL27_STIM_IMATURE_DC_DN)) +  geom_point(size=2)+theme_classic()
# ggplot( u_dat[ARM!="A"],aes(V1,V2,col=  GSE34156_NOD2_LIGAND_VS_TLR1_TLR2_LIGAND_6H_TREATED_MONOCYTE_DN)) +  geom_point(size=2)+theme_classic()
# 
# 
# ggplot( u_dat[ARM!="A"],aes(V1,V2,col=  HECKER_IFNB1_TARGETS)) +  geom_point(size=1)+theme_classic()
# ggplot( u_dat[ARM!="A"],aes(V1,V2,col=  NAGASHIMA_EGF_SIGNALING_UP)) +  geom_point(size=1)+theme_classic()
# ggplot( u_dat[ARM!="A"],aes(V1,V3,col=  HECKER_IFNB1_TARGETS)) +  geom_point(size=1)+theme_classic()
# 
# 
# 
# 
# 
# corResponse_Vs <- cor(log(u_dat[u_dat$V2<5][]$prop_change),u_dat[u_dat$V2<5][]%>%dplyr::select(V1:V2))
# methodcor= "pearson"
# methodcor= "spearman"
# #corV1 <-cor(u_dat[]$PC1,ss_dd[,],method=methodcor)#corV1<-cor(u_dat$V1,full_dd_transp)
# #corV1ord<-corV1[,order(-abs(corV1))]
# #corV2<-cor(u_dat[]$PC2,ss_dd[,],method=methodcor)#corV2<-cor(u_dat$V2,full_dd_transp[,])
# #corV2ord<-corV2[,order(-abs(corV2))]
# 
# corV1 <-cor(u_dat[u_dat$V2<5][]$V1,u_dat[u_dat$V2<5]%>%dplyr::select(colnames(ss_dd)),method=methodcor)#corV1<-cor(u_dat$V1,full_dd_transp)
# corV1ord<-corV1[,order(-abs(corV1))]
# corV2<-cor(u_dat[u_dat$V2<5][]$V2,u_dat[u_dat$V2<5]%>%dplyr::select(colnames(ss_dd)),method=methodcor)#corV2<-cor(u_dat$V2,full_dd_transp[,])
# corV2ord<-corV2[,order(-abs(corV2))]
# nn=20
# corV1ord[1:nn]
# corV2ord[1:nn]
# 
# 
# 
# ggplot( u_dat,aes(y=V1,x=log(1+Day),group=Day,fill=  dynamic_class3,col=  dynamic_class3)) +  geom_violin()+scale_color_npg(name="Tumor response")+scale_fill_npg(name="Tumor response")+theme_classic()+facet_wrap(ARM~dynamic_class3,ncol=2)
# ggplot( u_dat,aes(y=V2,x=log(1+Day),group=Day,fill=  dynamic_class3,col=  dynamic_class3)) +  geom_violin()+scale_color_npg(name="Tumor response")+scale_fill_npg(name="Tumor response")+theme_classic()+facet_wrap(ARM~dynamic_class3,ncol=2)
# ggplot( u_dat,aes(y=V3,x=log(1+Day),group=Day,fill=  dynamic_class3,col=  dynamic_class3)) +  geom_violin()+scale_color_npg(name="Tumor response")+scale_fill_npg(name="Tumor response")+theme_classic()+facet_wrap(ARM~dynamic_class3,ncol=2)
# ggplot( u_dat,aes(y=V4,x=log(1+Day),group=Day,fill=  dynamic_class3,col=  dynamic_class3)) +  geom_violin()+scale_color_npg(name="Tumor response")+scale_fill_npg(name="Tumor response")+theme_classic()+facet_wrap(ARM~dynamic_class3,ncol=2)
# 
# ggplot( u_dat,aes(V3,V2,col=  Patient.Study.ID)) +  geom_point(size=3)+theme_classic()+facet_wrap(ARM~dynamic_class3)
# ggplot( u_dat,aes(V1,V2,col=  dynamic_class3)) +  geom_point(size=3)+scale_color_npg(name="Tumor response")+theme_classic()
# 
# 
# 
# 
# 
# ggplot( u_dat[Celltype_subtype!="DC"],aes(-V2,-V1,col=  dynamic_class3)) +  geom_point(size=3)+scale_color_npg(name="Tumor response")+theme_classic()+facet_wrap(~dynamic_class3)
# ggplot( u_dat[Celltype_subtype!="DC"],aes(-V2,-V1,col=  dynamic_class3)) +  geom_point(size=3)+scale_color_npg(name="Tumor response")+theme_classic()+facet_wrap(dynamic_class3~Day)
# 
# ggplot( u_dat,aes(-V2,-V3,col=  Celltype_subtype)) +  geom_point(size=3)+scale_color_npg(name="Tumor response")+theme_classic()
# ggplot( u_dat,aes(-V1,-V3,col=  Celltype_subtype)) +  geom_point(size=3)+scale_color_npg(name="Tumor response")+theme_classic()
# ggplot( u_dat,aes(-V1,-V2,col=  Celltype_subtype)) +  geom_point(size=3)+scale_color_npg(name="Tumor response")+theme_classic()
# 
# ggplot( u_dat,aes(-V1,-V3,col=  Patient.Study.ID)) +  geom_point(size=0.9)+theme_classic()
# ggplot( u_dat,aes(-V2,-V3,col=  Patient.Study.ID)) +  geom_point(size=0.9)+theme_classic()
# ggplot( u_dat[Celltype_subtype!="DC"],aes(-V2,-V1,col=  Patient.Study.ID)) +  geom_point(size=0.9)+theme_classic()
# 
# 
# # V1=IFN activation
# # V2= TGFB_AND_IL4 activation 
# # V3=DC activation
# # V4=Immune suppression
# m1<-glm(yy~V4*V3*V2*V1,data=u_dat)
# u_dat[,yy:=0]
# u_dat[dynamic_class3=="Response",yy:=1]
# 
# require(GGally)
# 
# ggpairs(u_dat%>%dplyr::select(dynamic_class3,V1:V4), mapping = ggplot2::aes(color = dynamic_class3) , columns = 2:5)
# ggpairs(u_dat%>%dplyr::select(Celltype_subtype,V1:V4), mapping = ggplot2::aes(color = Celltype_subtype) , columns = 2:5)
# ggpairs(u_dat%>%dplyr::select(Patient.Study.ID,V1:V4), mapping = ggplot2::aes(color = Patient.Study.ID) , columns = 2:5)
# 
# ggpairs(u_dat[Day==14]%>%dplyr::select(Celltype_subtype,V1:V4), mapping = ggplot2::aes(color = Celltype_subtype) , columns = 2:5)
# ggpairs(u_dat[Day==14]%>%dplyr::select(dynamic_class3,V1:V4), mapping = ggplot2::aes(color = dynamic_class3) , columns = 2:5)
# 
# ggpairs(u_dat[Day==180]%>%dplyr::select(dynamic_class3,V1:V4), mapping = ggplot2::aes(color = dynamic_class3) , columns = 2:5)
# ggpairs(u_dat[Day==0]%>%dplyr::select(dynamic_class3,V1:V4), mapping = ggplot2::aes(color = dynamic_class3) , columns = 2:5)
# 
# ggplot( u_dat ,aes(x=V4,fill=  as.factor(Day),group=Day)) +  geom_density(alpha=0.5)+facet_wrap(~dynamic_class3)
# ggplot( u_dat ,aes(x=V3,fill=  as.factor(Day),group=Day)) +  geom_density(alpha=0.5)+facet_wrap(~dynamic_class3)
# ggplot( u_dat ,aes(x=V2,fill=  as.factor(Day),group=Day)) +  geom_density(alpha=0.5)+facet_wrap(~dynamic_class3)
# ggplot( u_dat ,aes(x=V1,fill=  as.factor(Day),group=Day)) +  geom_density(alpha=0.5)+facet_wrap(~dynamic_class3)
