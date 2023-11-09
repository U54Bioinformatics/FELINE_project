
rm(list=ls())   
##devtools::install_github("rikenbit/nnTensor") #install.packages("rTensor")
#install.packages("https://cran.r-project.org/src/contrib/Archive/rTensor/rTensor_1.4.tar.gz", repo=NULL, type="source")
require(rTensor); require(nnTensor); library(abind)
require(data.table); require(dplyr); require(ggplot2); require(tidyr)
require(igraph);require(ggsci)
require(lmerTest)

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
allphenotypes <- merge(allphenotypes %>% dplyr::select(-c( paste0("V", 1:5), "nCount_RNA", "nFeature_RNA", "Infercnv_CNA", "Platform","Timepoint","Sample_p_t","prop_change", "file_string", "day_fact" , "Ribo","TreatLab", "Burden_t0" ,"BurdenTracked" ,"Day0", "DayLastBurd", "Dose_lab","TreatCode","TreatCodeOrd","dynamic_class2","rgr_A","rgr_B")), tmp, by= c("Celltype", "Celltype_subtype") )

### Unique clusters of cells (ALL) and their umap discretization level
uu <- unique( allphenotypes %>% dplyr::select( c("Celltype","PhenoCelltype", "key_", paste0("Disc_V", 1:5) ) ) )  #"Celltype_subtype",

### Load Ligand Receptor database list of Ramilowski et al 2015
load( "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/LigandReceptor/Filtered_Human-2015-Ramilowski-LR-pairs.RData")
LRgenelist <- unique( c(LRpairsFiltered$HPMR.Receptor, LRpairsFiltered$HPMR.Ligand) )

### Load some signalling data for downstream analysis
savelocCCI <- "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/CommunicationOutputCohort1indivlevelAllArms/"

# Are we studying the per indiv or per cell type level of signalling?
perIndiv=TRUE

# Names of LR data to analyze
filenamesCCI <- list.files(savelocCCI)#%>%length

#summarytable <- read.csv("/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/TensorAnalysisOutput/Ligand Receptor list fro NTD model of population cellcell communication_AB.csv")

# extract cell cell interactions from communication data.table
CommCancMacro <- rbindlist(lapply(1:length(filenamesCCI), function(ii){
  load(file= paste0(savelocCCI, filenamesCCI[ii]))
  cat(ii);
  SumComm <- data.table( Communication[(LigandPhenoCelltype=="Cancer cells"&ReceptorPhenoCelltype=="Macrophages")|(LigandPhenoCelltype=="Macrophages"&ReceptorPhenoCelltype=="CD8+ T cells")] )
  return(SumComm)
}))
CommCancMacro[, RiboTreated:=TRUE]
CommCancMacro[ARM=="A", RiboTreated:=FALSE]
CommCancMacro[, Treat:="CombinationRibo"]
CommCancMacro[ARM=="A", Treat:="LetrozoleAlone"]

pnamelist <- unique(CommCancMacro$Pair.Name)

macrosubset <- CommCancMacro[ReceptorPhenoCelltype=="Macrophages"]
Tcellsubset <- CommCancMacro[LigandPhenoCelltype=="Macrophages"]
pnamelist <- data.table(macrosubset%>%group_by(Pair.Name)%>%summarise(n=length(unique(Disc_V2))))[n>2]$Pair.Name  #data.table(table(macrosubset$Pair.Name))


macrosubset[,Differentiation:= "M1"]
macrosubset[Disc_V2<=3, Differentiation:= "M2"]

#macrosubset[,Differentiation:= "M1"]
#macrosubset[Disc_V1>2, Differentiation:= "M2"]

M1_M2keys<-unique(macrosubset %>% select(Differentiation,key_))[order(Differentiation,key_)]
setnames(M1_M2keys,old="key_",new="key_signaller")
Tcellsubset <- merge(Tcellsubset, M1_M2keys, by="key_signaller")


#save( Tcellsubset , file= "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Macrophage T cell communication/Ribo and Letrozole M1 vs M2 individual Macrophage T cell communication summary Cohort1.RData")
load(file= "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Macrophage T cell communication/Ribo and Letrozole M1 vs M2 individual Macrophage T cell communication summary Cohort1.RData")

# Cancer to M2 vs M1
Summar_MC <- data.table(Tcellsubset%>%
                          group_by(Differentiation,Pair.Name,Treat,Day,dynamic_class3,Disc_V2,Disc_V1,Patient.Study.ID,ReceptorPhenoCelltype) %>%
                          summarise(Signal=median(Signal),
                                    TransductionMu=mean(Signal) ))
Summar_MC[ , muD0:= sum(TransductionMu*(Day==0) )/sum((Day==0)), by=c("Patient.Study.ID" ,"Treat","dynamic_class3", "Pair.Name") ]
Summar_MC[,c("L","R"):= tstrsplit(Pair.Name,"_", fixed=T)]
Summar_MC[,ln1Day:= log(1+Day)]
Summar_MC[,Fact_Patient.Study.ID:= as.factor(Patient.Study.ID)]

countPairdd <- data.table(Summar_MC %>% group_by(Treat,Pair.Name,Differentiation) %>% summarise(n=n()) %>% group_by(Treat,Pair.Name)%>%
                            summarise(n=min(n))
)
Summar_MC<-Summar_MC[Pair.Name %in% countPairdd[n>10]$Pair.Name]


TrendsResMfromCV2diffRibo <- rbindlist( lapply( unique( Summar_MC$Pair.Name )  , function(pp){
  cat(pp)
  out <- tryCatch({
    mm1 <- lmer( I(log(1 + TransductionMu) ) ~ Differentiation +(1|Disc_V1) +(1|Disc_V2) +  (1|Patient.Study.ID),  #(1|Disc_V1) +
                 data= Summar_MC[Pair.Name==pp][Treat=="CombinationRibo"] )
    if( length(coef(summary(mm1)))==0 ){ stop() }else{
      data.table(Treat="CombinationRibo",Pair.Name= pp, data.table( coef(summary(mm1)) , keep.rownames = T))
      
    }
  },
  error=function(x){
    data.table(Treat="CombinationRibo",'Pair.Name'=pp, data.table("rn"=NA, "Estimate"=NA ,'Std. Error'=NA, df =NA,'t value'=NA,'Pr(>|t|)'=NA) )
  })
  return(out)
}))
TrendsResMfromCV2diffLetro <- rbindlist( lapply( unique( Summar_MC$Pair.Name )  , function(pp){
  cat(pp)
  out <- tryCatch({
    mm1 <- lmer( I(log(1 + TransductionMu) ) ~ Differentiation +(1|Disc_V1) +(1|Disc_V2)  + (1|Patient.Study.ID), #+ (1|Disc_V1)
                 data= Summar_MC[Pair.Name==pp][Treat=="LetrozoleAlone"] )
    if( length(coef(summary(mm1)))==0 ){ stop() }else{
      data.table(Treat="LetrozoleAlone",Pair.Name= pp, data.table( coef(summary(mm1)) , keep.rownames = T))
    }
  },
  error=function(x){
    data.table(Treat="LetrozoleAlone",'Pair.Name'=pp, data.table("rn"=NA, "Estimate"=NA ,'Std. Error'=NA, df =NA,'t value'=NA,'Pr(>|t|)'=NA) )
  })
  return(out)
}))


setnames(TrendsResMfromCV2diffRibo,old=c("Std. Error","t value","Pr(>|t|)"),new=c("Std.Error","tval","pval"))
TrendsResMfromCV2diffRibo[,c("L","R"):= tstrsplit(Pair.Name,"_",fixed=T)]
TrendsResMfromCV2diffRibo$ adjpval<- p.adjust(TrendsResMfromCV2diffRibo$pval)
setnames(TrendsResMfromCV2diffLetro,old=c("Std. Error","t value","Pr(>|t|)"),new=c("Std.Error","tval","pval"))
TrendsResMfromCV2diffLetro[,c("L","R"):= tstrsplit(Pair.Name,"_",fixed=T)]
TrendsResMfromCV2diffLetro$ adjpval<- p.adjust(TrendsResMfromCV2diffLetro$pval)

# here is my list
TrendsResMfromCV2diffRibo[grep("Differentiation",rn)][Estimate<0][!grep("ITG",Pair.Name)][order(pval)][pval<0.05][1:20]
TrendsResMfromCV2diffRibo[grep("Differentiation",rn)][Estimate>0][!grep("ITG",Pair.Name)][order(pval)][pval<0.05][1:20]
TrendsResMfromCV2diffLetro[grep("Differentiation",rn)][Estimate<0][order(pval)][pval<0.05][1:20]
TrendsResMfromCV2diffLetro[grep("Differentiation",rn)][Estimate>0][order(pval)][pval<0.05][1:20]




#### CXCL10_CXCR3: CXCL10 promotes T cell activation
# CXCL9_CXCR3 CXCL9 promotes T cell activation
# INHBA_ACVR2A : INHBA prevents Treg induction (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6325588/)
# SORBS1_INSR SORBS1: binding increasea glucose consumption and accompanies T cell activation
#SPP1_CD44 : SPP1 == early T lymphocyte activation 1, induces T cell survival (https://www.nature.com/articles/ni0107-19)

# CXCL13_CXCR3 : CXCL13  is an agonist for the human CXCR3 receptor preventing T cell activation (https://pubmed.ncbi.nlm.nih.gov/11554781/)
# SEMA7A_ITGA1:  SEMA7A as a critical role in negative regulation of T cell activation and function (https://pubmed.ncbi.nlm.nih.gov/16713976/#:~:text=Semaphorin%207A%20(Sema7A)%20promotes%20axonal,and%20antigen%2Dinduced%20proliferative%20response.)
# THBS1_CD47: THBS1 is a potent inhibitor of T cell and dendritic cell activation and mediates clearance of apoptotic cells by phagocytes
#             CD47  https://www.jimmunol.org/content/180/12/8073#:~:text=Although%20CD47%20ligation%20appears%20to,functions%20as%20an%20inhibitory%20molecule.



listplotpairs <- c("CXCL9_CXCR3", "INHBA_ACVR2A","SORBS1_INSR",  "CXCL13_CXCR3", "SEMA7A_ITGA1","THBS1_CD47")
ddplot <- Summar_MC[Pair.Name%in% listplotpairs  ]#[Pair.Name != "VTN_ITGB6"]
ddplot[,scalelntransduction := scale( log(1 + TransductionMu) ) , by= Pair.Name ]
ddplot$Pair.Name <- factor( ddplot$Pair.Name, levels = unique(listplotpairs))
ddplot[,Treatment:= "Combination ribociclib"]
ddplot[Treat=="LetrozoleAlone",Treatment:= "Letrozole alone"]


ggplot(ddplot, aes(x=Pair.Name, y= scalelntransduction ,col= Differentiation,fill= Differentiation  )) + 
  theme_classic(base_size=18)+
  geom_point(position= position_dodge(width=0.9))+
  geom_boxplot(aes(group= interaction(Differentiation,Pair.Name)),position= position_dodge(width=0.9),scale="width",alpha=0.65)+
  labs(y="Myeloid to CD8+ T cell \n communication",x="Communication pathway"  )+
  scale_color_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
  scale_fill_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
  facet_wrap(Treatment~dynamic_class3)+theme(aspect.ratio=0.6,axis.text.x = element_text(angle=90,vjust=0.5),legend.position = "bottom")
ggplot(ddplot, aes(x=Pair.Name, y= scalelntransduction ,col= Differentiation,fill= Differentiation  )) + 
  theme_classic(base_size=18)+
  geom_point(position= position_dodge(width=0.9))+
  geom_boxplot(aes(group= interaction(Differentiation,Pair.Name)),position= position_dodge(width=0.9),scale="width",alpha=0.65)+
  labs(y="Myeloid to CD8+ T cell \n cell-cell communication",x="Communication pathway"  )+
  scale_color_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
  scale_fill_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
  facet_wrap(~Treatment,ncol=1)+theme(aspect.ratio=0.6,axis.text.x = element_text(angle=90,vjust=0.5),legend.position = "bottom")+
  geom_vline(xintercept=3.5,linetype="dashed")
ggsave(filename=paste0("/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cohort 2 results/M2 vs M1 Macrophages to T cell communication differences.png"),width=8,height=8)

ggplot(ddplot, aes(x=Pair.Name, y= scalelntransduction ,col= Differentiation,fill= Differentiation  )) + 
  theme_classic(base_size=18)+
  geom_point(position= position_dodge(width=0.9))+
  geom_boxplot(aes(group= interaction(Differentiation,Pair.Name)),position= position_dodge(width=0.9),scale="width",alpha=0.65)+
  labs(y="Myeloid to CD8+ T cell \n cell-cell communication",x="Communication pathway"  )+
  scale_color_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
  scale_fill_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
  theme(aspect.ratio=0.6,axis.text.x = element_text(angle=90,vjust=0.5),legend.position = "bottom")+
  geom_vline(xintercept=3.5,linetype="dashed")








TrendsRibo_M1diffResist<- rbindlist( lapply( unique( Summar_MC$Pair.Name )  , function(pp){
  cat(pp)
  out <- tryCatch({
    mm1 <- lmer( I(log(1 + TransductionMu) ) ~ Day*dynamic_class3 +(1|Disc_V1) +(1|Disc_V2) +  (1|Patient.Study.ID),  #(1|Disc_V1) +
                 data= Summar_MC[Pair.Name==pp][Treat=="CombinationRibo"][Differentiation=="M1"] )
    if( length(coef(summary(mm1)))==0 ){ stop() }else{
      data.table(Treat="CombinationRibo",differentiation="M1",Pair.Name= pp, data.table( coef(summary(mm1)) , keep.rownames = T))
      
    }
  },
  error=function(x){
    data.table(Treat="CombinationRibo",differentiation="M1", 'Pair.Name'=pp, data.table("rn"=NA, "Estimate"=NA ,'Std. Error'=NA, df =NA,'t value'=NA,'Pr(>|t|)'=NA) )
  })
  return(out)
}))
TrendsLetrozole_M1diffResist<- rbindlist( lapply( unique( Summar_MC$Pair.Name )  , function(pp){
  cat(pp)
  out <- tryCatch({
    mm1 <- lmer( I(log(1 + TransductionMu) ) ~ Day*dynamic_class3 +(1|Disc_V1) +(1|Disc_V2) +  (1|Patient.Study.ID),  #(1|Disc_V1) +
                 data= Summar_MC[Pair.Name==pp][Treat=="LetrozoleAlone"][Differentiation=="M1"] )
    if( length(coef(summary(mm1)))==0 ){ stop() }else{
      data.table(Treat="LetrozoleAlone",differentiation="M1",Pair.Name= pp, data.table( coef(summary(mm1)) , keep.rownames = T))
    }
  },
  error=function(x){
    data.table(Treat="LetrozoleAlone",differentiation="M1", 'Pair.Name'=pp, data.table("rn"=NA, "Estimate"=NA ,'Std. Error'=NA, df =NA,'t value'=NA,'Pr(>|t|)'=NA) )
  })
  return(out)
}))



setnames(TrendsRibo_M1diffResist,old=c("Std. Error","t value","Pr(>|t|)"),new=c("Std.Error","tval","pval"))
TrendsRibo_M1diffResist[,c("L","R"):= tstrsplit(Pair.Name,"_",fixed=T)]
TrendsRibo_M1diffResist$ adjpval<- p.adjust(TrendsRibo_M1diffResist$pval)
setnames(TrendsLetrozole_M1diffResist,old=c("Std. Error","t value","Pr(>|t|)"),new=c("Std.Error","tval","pval"))
TrendsLetrozole_M1diffResist[,c("L","R"):= tstrsplit(Pair.Name,"_",fixed=T)]
TrendsLetrozole_M1diffResist$ adjpval<- p.adjust(TrendsLetrozole_M1diffResist$pval)

# here is my list
TrendsRibo_M1diffResist[grep("Day",rn)][Estimate<0][order(pval)][pval<0.05][1:20]
TrendsRibo_M1diffResist[grep("Day",rn)][Estimate>0][order(pval)][pval<0.05][1:20]

TrendsRibo_M1diffResist[grep("Response",rn)][order(Estimate)][pval<0.05][Estimate>0][1:20]
TrendsRibo_M1diffResist[grep("Response",rn)][Estimate>0][!grep("ITG",Pair.Name)][order(pval)][pval<0.05][1:20]
TrendsRibo_M1diffResist[grep("Response",rn)][Estimate<0][!grep("ITG",Pair.Name)][order(pval)][pval<0.05][1:20]


ddplot <- Summar_MC[Pair.Name%in%c("CCL19_CXCR3","CXCL10_CXCR3","CCL19_CCR7","CCL21_CCR7") | grepl("IL15",Pair.Name)| grepl("IL18",Pair.Name) ]#[grepl("CCR5",Pair.Name)]
#ddplot <- Summar_MC[Pair.Name%in%c("CCL4_CCR5","CCL3_CCR5","CXCL9_CXCR3","CCL3_CCR1","CCL4_CCR1") ]
#Summar_MC[grepl("IL15",Pair.Name)| grepl("IL18",Pair.Name) ]#[Pair.Name != "VTN_ITGB6"]
ddplot[,scalelntransduction := scale( log(1 + Signal) ) , by= Pair.Name ]
ddplot[,Treatment:= "Combination ribociclib"]
ddplot[Treat=="LetrozoleAlone",Treatment:= "Letrozole alone"]

unique(ddplot$Pair.Name)
summary(lm(scalelntransduction~dynamic_class3, ddplot[Differentiation=="M1"][Treatment=="Combination ribociclib"][Pair.Name%in% "IL15_IL2RB"] ))
summary(lm(scalelntransduction~dynamic_class3, ddplot[Differentiation=="M1"][Treatment=="Combination ribociclib"][Pair.Name%in% "IL15_IL2RB"][Day==180] ))


ddplot[,DayLab:= paste0("Day ",Day)]
ddplot[,Pair.NameLab:= gsub("_","-",Pair.Name)]

TrendsRibo_M1diffResist[Pair.Name%in% unique(ddplot$Pair.Name)]
ggplot(ddplot[][Differentiation=="M1"]%>%group_by(Differentiation,Pair.NameLab)%>%mutate(scalelntransduction=scale(scalelntransduction)), 
       aes(x=dynamic_class3, y= scalelntransduction ,
           col= Pair.NameLab,
           fill= Pair.NameLab  )) + 
  theme_classic(base_size=18)+
  geom_point(position= position_dodge(width=0.9))+
  geom_boxplot(aes(group= interaction(dynamic_class3,Differentiation,Pair.NameLab)),position= position_dodge(width=0.9),alpha=0.65)+
  labs(y="M1 myeloid to CD8+ T cell \n cell-cell communication",x="Communication pathway"  )+
  scale_x_discrete(name="Tumor response", labels=c("Resistant","Sensitive") ) +
  scale_color_viridis(name="Communication pathway", discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
  scale_fill_viridis(name="Communication pathway", discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
  facet_wrap(Treatment~DayLab,nrow=2) + theme(aspect.ratio=1,
                                              #axis.text.x = element_text(angle=90,vjust=0.5),
                                              legend.position = "bottom")+
  geom_vline(xintercept=1.5,linetype="dashed")

ggsave(filename=paste0("/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cohort 2 results/M1 Macrophages to T cell communication differences between Resistant and sensitive tumors.png"),width=12.5,height=10)


ggplot(ddplot[][Differentiation=="M1"]%>%group_by(Differentiation,Pair.NameLab)%>%mutate(scalelntransduction=scale(scalelntransduction)), 
       aes(x=dynamic_class3, y= scalelntransduction ,
           col= Pair.NameLab,
           fill= Pair.NameLab  )) + 
  theme_classic(base_size=18)+
  geom_point(position= position_dodge(width=0.9))+
  geom_boxplot(aes(group= interaction(dynamic_class3,Differentiation,Pair.NameLab)),position= position_dodge(width=0.9),alpha=0.65)+
  labs(y="M1 myeloid to CD8+ T cell \n cell-cell communication",x="Communication pathway"  )+
  scale_x_discrete(name="Tumor response", labels=c("Resistant","Sensitive") ) +
  scale_color_viridis(name="Communication pathway", discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
  scale_fill_viridis(name="Communication pathway", discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
  facet_wrap(~Treatment,nrow=2) + theme(aspect.ratio=1
                                              #axis.text.x = element_text(angle=90,vjust=0.5),
                                              )+
  geom_vline(xintercept=1.5,linetype="dashed")

ggsave(filename=paste0("/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cohort 2 results/M1 Macrophages to T cell communication differences between Resistant and sensitive tumors gather times.png"),width=12.5,height=6)


paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Ribo and Letrozole M1 Macrophages to T cell communication differences between Resistant and sensitive tumors gather times.png"),height=10,width=10)

ggplot(ddplot[][Differentiation=="M1"]%>%group_by(Differentiation,Pair.NameLab)%>%mutate(scalelntransduction=scale(scalelntransduction)), 
       aes(x=dynamic_class3, y= scalelntransduction ,
           col= Pair.NameLab,
           fill= Pair.NameLab  )) + 
  theme_classic(base_size=18)+
  geom_point(position= position_dodge(width=0.9))+
  geom_boxplot(aes(group= interaction(dynamic_class3,Differentiation,Pair.NameLab)),position= position_dodge(width=0.9),alpha=0.65)+
  labs(y="M1 myeloid to CD8+ T cell \n cell-cell communication",x="Communication pathway"  )+
  scale_x_discrete(name="Tumor response", labels=c("Resistant","Sensitive") ) +
  scale_color_viridis(name="Communication pathway", discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
  scale_fill_viridis(name="Communication pathway", discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
  facet_wrap(~Treatment,nrow=2) + theme(aspect.ratio=1
                                        #axis.text.x = element_text(angle=90,vjust=0.5),
  )+
  geom_vline(xintercept=1.5,linetype="dashed")+
  theme(axis.title=element_blank(),  axis.text=element_blank(),strip.text = element_blank(),legend.title=element_blank(),legend.text=element_blank() )

ggsave(paste0(paperfile,"BLANK Ribo and Letrozole M1 Macrophages to T cell communication differences between Resistant and sensitive tumors gather times.png"),height=10,width=10)


ggplot(ddplot[][Differentiation=="M2"]%>%group_by(Differentiation,Pair.Name)%>%mutate(scalelntransduction=scale(scalelntransduction)), aes(x=dynamic_class3, y= scalelntransduction ,col= Pair.Name,fill= Pair.Name  )) + 
  theme_classic(base_size=18)+
  geom_point(position= position_dodge(width=0.9))+
  geom_boxplot(aes(group= interaction(dynamic_class3,Differentiation,Pair.Name)),position= position_dodge(width=0.9),alpha=0.65)+
  labs(y="Myeloid to CD8+ T cell \n communication",x="Communication pathway"  )+
  #scale_color_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
  #scale_fill_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
  facet_wrap(Treatment~Day,nrow=2)+theme(aspect.ratio=0.6,axis.text.x = element_text(angle=90,vjust=0.5),legend.position = "bottom")



# 
# 
# ggplot(ddplot[Day==180], aes(x=Pair.Name, y= scalelntransduction ,col= Differentiation,fill= Differentiation  )) + 
#   theme_classic(base_size=18)+
#   geom_point(position= position_dodge(width=0.9))+
#   geom_boxplot(aes(group= interaction(Differentiation,Pair.Name)),position= position_dodge(width=0.9),scale="width",alpha=0.65)+
#   labs(y="Myeloid to CD8+ T cell \n communication",x="Communication pathway"  )+
#   scale_color_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
#   scale_fill_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
#   facet_wrap(Treatment~dynamic_class3,nrow=2)+theme(aspect.ratio=0.6,axis.text.x = element_text(angle=90,vjust=0.5),legend.position = "bottom")
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ggplot(ddplot[Pair.Name=="THBS1_CD47"][Treatment=="Combination ribociclib"], aes(x=Disc_V2, y=Disc_V1 , fill= scalelntransduction ,col= Differentiation,fill= Differentiation  )) + 
#   theme_classic(base_size=18)+
#   geom_tile()+
#   #  geom_point(position= position_dodge(width=0.9))+
#   #  geom_boxplot(aes(group= interaction(Differentiation,Pair.Name)),position= position_dodge(width=0.9),scale="width",alpha=0.65)+
#   #  labs(y="Myeloid to CD8+ T cell \n communication",x="Communication pathway"  )+
#   #  scale_color_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
#   #  scale_fill_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
#   facet_wrap(Differentiation~Day,nrow=2)+theme(aspect.ratio=0.6,axis.text.x = element_text(angle=90,vjust=0.5),legend.position = "bottom")
# 
# ggplot(ddplot[Pair.Name=="CXCL9_CXCR3"][Treatment=="Combination ribociclib"], aes(x=Disc_V2, y=Disc_V1 , fill= scalelntransduction ,col= Differentiation,fill= Differentiation  )) + 
#   theme_classic(base_size=18)+
#   geom_tile()+
#   #  geom_point(position= position_dodge(width=0.9))+
#   #  geom_boxplot(aes(group= interaction(Differentiation,Pair.Name)),position= position_dodge(width=0.9),scale="width",alpha=0.65)+
#   #  labs(y="Myeloid to CD8+ T cell \n communication",x="Communication pathway"  )+
#   #  scale_color_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
#   #  scale_fill_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
#   facet_wrap(Differentiation~Day,nrow=2)+theme(aspect.ratio=0.6,axis.text.x = element_text(angle=90,vjust=0.5),legend.position = "bottom")
# 
# 
# ddplot <- Tcellsubset[grepl("IL2",Pair.Name)| grepl("IL15",Pair.Name)| grepl("IL18",Pair.Name) ]
# #Summar_MC[grepl("IL15",Pair.Name)| grepl("IL18",Pair.Name) ]#[Pair.Name != "VTN_ITGB6"]
# ddplot[,scalelntransduction := scale( log(1 + Signal) ) , by= Pair.Name ]
# ddplot[,Treatment:= "Combination ribociclib"]
# ddplot[Treat=="LetrozoleAlone",Treatment:= "Letrozole alone"]
# 
# 
# ggplot(ddplot[Day==180], aes(x=Pair.Name, y= scalelntransduction ,col= Differentiation,fill= Differentiation  )) + 
#   theme_classic(base_size=18)+
#   geom_point(position= position_dodge(width=0.9))+
#   geom_boxplot(aes(group= interaction(Differentiation,Pair.Name)),position= position_dodge(width=0.9),scale="width",alpha=0.65)+
#   labs(y="Myeloid to CD8+ T cell \n communication",x="Communication pathway"  )+
#   scale_color_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
#   scale_fill_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
#   facet_wrap(Treatment~dynamic_class3,nrow=2)+theme(aspect.ratio=0.6,axis.text.x = element_text(angle=90,vjust=0.5),legend.position = "bottom")
# 
# 
# ggplot(ddplot[][Differentiation=="M1"]%>%group_by(Differentiation,Pair.Name)%>%mutate(scalelntransduction=scale(scalelntransduction)), aes(x=dynamic_class3, y= scalelntransduction ,col= Pair.Name,fill= Pair.Name  )) + 
#   theme_classic(base_size=18)+
#   geom_point(position= position_dodge(width=0.9))+
#   geom_boxplot(aes(group= interaction(dynamic_class3,Differentiation,Pair.Name)),position= position_dodge(width=0.9),alpha=0.65)+
#   labs(y="Myeloid to CD8+ T cell \n communication",x="Communication pathway"  )+
#   #scale_color_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
#   #scale_fill_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
#   facet_wrap(Treatment~Day,nrow=2)+theme(aspect.ratio=0.6,axis.text.x = element_text(angle=90,vjust=0.5),legend.position = "bottom")
# 
# 
# 
# ggplot(ddplot[Pair.Name=="IL15_IL15RA"][Treatment=="Combination ribociclib"], aes(x=Disc_V2, y=Disc_V1 , fill= scalelntransduction ,col= Differentiation,fill= Differentiation  )) + 
#   theme_classic(base_size=18)+
#   geom_tile()+
# #  geom_point(position= position_dodge(width=0.9))+
# #  geom_boxplot(aes(group= interaction(Differentiation,Pair.Name)),position= position_dodge(width=0.9),scale="width",alpha=0.65)+
# #  labs(y="Myeloid to CD8+ T cell \n communication",x="Communication pathway"  )+
# #  scale_color_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
# #  scale_fill_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
#   facet_wrap(Differentiation~dynamic_class3,nrow=2)+theme(aspect.ratio=0.6,axis.text.x = element_text(angle=90,vjust=0.5),legend.position = "bottom")
# 
# 
# 
# 
# 
# # 
# # 
# # 
# # 
# # TrendsResMfromCV2diffLetro[][order(pval)][adjpval<0.05][rn=="DifferentiationM2"]
# # 
# # TrendsResMfromCV2diffRibo[grep("Differentiation",rn)][order(pval)][adjpval<0.05][rn=="DifferentiationM2"]
# # TrendsResMfromCV2diffRibo[grep("Differentiation",rn)][order(pval)][pval<0.05][1:20]
# # TrendsResMfromCV2diffRibo[grep("Differentiation",rn)][order(pval)][adjpval<0.05][rn=="DifferentiationM2"]
# # TrendsResMfromCV2diffRibo[grep("Differentiation",rn)][Estimate<0][order(pval)][pval<0.05][1:20]
# # TrendsResMfromCV2diffRibo[grep("Differentiation",rn)][Estimate<0][order(pval)][pval<0.05][1:20]
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # TrendsResMfromCV2diffRibo <- rbindlist( lapply( unique( Summar_MC$Pair.Name )  , function(pp){
# #   cat(pp)
# #   out <- tryCatch({
# #     mm1 <- lmer( I(log(1 + TransductionMu) ) ~ Differentiation +(1|Disc_V1) +(1|Disc_V2) +  (1|Patient.Study.ID),  #(1|Disc_V1) +
# #                  data= Summar_MC[Pair.Name==pp][Treat=="CombinationRibo"] )
# #     if( length(coef(summary(mm1)))==0 ){ stop() }else{
# #       data.table(Treat="CombinationRibo",Pair.Name= pp, data.table( coef(summary(mm1)) , keep.rownames = T))
# #       
# #     }
# #   },
# #   error=function(x){
# #     data.table(Treat="CombinationRibo",'Pair.Name'=pp, data.table("rn"=NA, "Estimate"=NA ,'Std. Error'=NA, df =NA,'t value'=NA,'Pr(>|t|)'=NA) )
# #   })
# #   return(out)
# # }))
# # 
# # TrendsResMfromCV2diffLetro <- rbindlist( lapply( unique( Summar_MC$Pair.Name )  , function(pp){
# #   cat(pp)
# #   out <- tryCatch({
# #     mm1 <- lmer( I(log(1 + TransductionMu) ) ~ Differentiation +(1|Disc_V1) +(1|Disc_V2)  + (1|Patient.Study.ID), #+ (1|Disc_V1)
# #                  data= Summar_MC[Pair.Name==pp][Treat=="LetrozoleAlone"] )
# #     if( length(coef(summary(mm1)))==0 ){ stop() }else{
# #       data.table(Treat="LetrozoleAlone",Pair.Name= pp, data.table( coef(summary(mm1)) , keep.rownames = T))
# #     }
# #   },
# #   error=function(x){
# #     data.table(Treat="LetrozoleAlone",'Pair.Name'=pp, data.table("rn"=NA, "Estimate"=NA ,'Std. Error'=NA, df =NA,'t value'=NA,'Pr(>|t|)'=NA) )
# #   })
# #   return(out)
# # }))
# # 
# # setnames(TrendsResMfromCV2diffRibo,old=c("Std. Error","t value","Pr(>|t|)"),new=c("Std.Error","tval","pval"))
# # TrendsResMfromCV2diffRibo[,c("L","R"):= tstrsplit(Pair.Name,"_",fixed=T)]
# # TrendsResMfromCV2diffRibo$ adjpval<- p.adjust(TrendsResMfromCV2diffRibo$pval)
# # TrendsResMfromCV2diffRibo[][order(pval)][adjpval<0.05][rn=="DifferentiationM2"]
# # TrendsResMfromCV2diffRibo[][order(pval)][pval<0.05][rn=="DifferentiationM2"][1:20]
# # TrendsResMfromCV2diffRibo[][order(pval)][adjpval<0.05][rn=="DifferentiationM2"]
# # setnames(TrendsResMfromCV2diffLetro,old=c("Std. Error","t value","Pr(>|t|)"),new=c("Std.Error","tval","pval"))
# # TrendsResMfromCV2diffLetro[,c("L","R"):= tstrsplit(Pair.Name,"_",fixed=T)]
# # TrendsResMfromCV2diffLetro$ adjpval<- p.adjust(TrendsResMfromCV2diffLetro$pval)
# # TrendsResMfromCV2diffLetro[][order(pval)][adjpval<0.05][rn=="DifferentiationM2"]
# # 
# # TrendsResMfromCV2diffRibo[Estimate<0][order(pval)][pval<0.05][rn=="DifferentiationM2"][1:20]
# # 
# # nn <- 5#25
# # 
# # short.table <- rbind(data.table(ChangeInM2 = "up", TrendsResMfromCV2diffRibo[Estimate>0][order(pval)][pval<0.05][grep("DifferentiationM2",rn)][1:nn]),
# #                      data.table(ChangeInM2 = "up",TrendsResMfromCV2diffLetro[Estimate>0][order(pval)][pval<0.05][grep("DifferentiationM2",rn)][1:nn]),
# #                      data.table(ChangeInM2 = "down",TrendsResMfromCV2diffRibo[Estimate<0][order(pval)][pval<0.05][grep("DifferentiationM2",rn)][1:nn]),
# #                      data.table(ChangeInM2 = "down",TrendsResMfromCV2diffLetro[Estimate<0][order(pval)][pval<0.05][grep("DifferentiationM2",rn)][1:nn]))[order(Estimate)]
# # 
# # listplotpairs <- unique( short.table$Pair.Name )
# # listplotpairs <- unique( c(
# #   TrendsResMfromCV2diffRibo[Estimate>0][order(pval)][pval<0.05][grep("DifferentiationM2",rn)][1:nn]$Pair.Name,
# #   TrendsResMfromCV2diffLetro[Estimate>0][order(pval)][pval<0.05][grep("DifferentiationM2",rn)][1:nn]$Pair.Name,
# #   TrendsResMfromCV2diffRibo[Estimate<0][order(pval)][pval<0.05][grep("DifferentiationM2",rn)][1:nn]$Pair.Name,
# #   TrendsResMfromCV2diffLetro[Estimate<0][order(pval)][pval<0.05][grep("DifferentiationM2",rn)][1:nn]$Pair.Name))
# # 
# # listplotpairs <- c("IL15_IL15RA","IL18_IL18R1")
# # ddplot <- Summar_MC[Pair.Name%in% listplotpairs  ]#[Pair.Name != "VTN_ITGB6"]
# # ddplot[,scalelntransduction := scale( log(1 + TransductionMu) ) , by= Pair.Name ]
# # ddplot$Pair.Name <- factor( ddplot$Pair.Name, levels = unique(listplotpairs))
# # ddplot[,Treatment:= "Combination ribociclib"]
# # ddplot[Treat=="LetrozoleAlone",Treatment:= "Letrozole alone"]
# # 
# # 
# # ggplot(ddplot, aes(x=Pair.Name, y= scalelntransduction ,col= Differentiation,fill= Differentiation  )) + 
# #   theme_classic(base_size=18)+
# #   geom_point(position=position_dodge(width=0.9))+
# #   geom_violin(aes(group=interaction(Differentiation,Pair.Name)),position=position_dodge(width=0.9),col=NA,scale="width",alpha=0.65)+
# #   labs(y="Cancer to Myeloid communication",x="Communication pathway"  )+
# #   scale_color_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
# #   scale_fill_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
# #   facet_wrap(~Treatment,scale="free_y",nrow=4)+theme(axis.text.x = element_text(angle=90),legend.position = "bottom")
# # 
# # 
# # ggplot(ddplot[Treat=="CombinationRibo"], aes(x=Disc_V1, y= scalelntransduction ,col= Differentiation,fill= Differentiation  )) + 
# #   theme_classic(base_size=18)+
# #   geom_point(position=position_dodge(width=0.9))+
# #   geom_violin(aes(group=interaction(Differentiation,Disc_V1)),position=position_dodge(width=0.9),col=NA,scale="width",alpha=0.65)+
# #   labs(y="Cancer to Myeloid communication",x="Communication pathway"  )+
# #   scale_color_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
# #   scale_fill_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
# #   facet_wrap(~Disc_V2, scale="free_y",nrow=4)+theme(axis.text.x = element_text(angle=90),legend.position = "bottom")
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # IL10_IL10RA
# # CXCL9_CXCR3
# # TGFB2_TGFBR2
# # HBEGF_CD44
# # CCL8_CCR2
# # SORBS1_INSR
# # WNT5A_ROR2
# # 
# # VCAM1_ITGA4
# # VCAM1_ITGB2
# # listplotpairs <- c("IL10_IL10RA","CXCL9_CXCR3","TGFB2_TGFBR2","HBEGF_CD44","CCL8_CCR2")
# # 
# # listplotpairs <- c("IL15_IL15RA","IL18_IL18R1","IL10_IL10RA","CXCL9_CXCR3","CCL8_CCR2")
# # 
# # listplotpairs <- c("CXCL9_CXCR3","INHBA_ACVR2A","SORBS1_INSR","SPP1_CD44","OSM_IL6ST","CCL21_CCR7","IL1R1",
# #                    "CXCL13_CXCR3" )
# # 
# # listplotpairs <- c("CXCL9_CXCR3","CXCL13_CXCR3", "INHBA_ACVR2A","NODAL_ACVR2A","SORBS1_INSR","SPP1_CD44","OSM_IL6ST","CCL21_CCR7","IL1R1"
# #                     )
# # ddplot <- Tcellsubset[Pair.Name%in% listplotpairs  ]#[Pair.Name != "VTN_ITGB6"]
# # ddplot <- Tcellsubset[grepl("IL15_",Pair.Name)|grepl("IL18_",Pair.Name) ]#[Pair.Name != "VTN_ITGB6"]
# # ddplot[,scalelntransduction := scale( log(1 + Signal) ) , by= Pair.Name ]
# # #ddplot$Pair.Name <- factor( ddplot$Pair.Name, levels = unique(listplotpairs))
# # ddplot[,Treatment:= "Combination ribociclib"]
# # ddplot[Treat=="LetrozoleAlone",Treatment:= "Letrozole alone"]
# # 
# # ggplot(ddplot, aes(x=Pair.Name, y= scalelntransduction ,col= Differentiation,fill= Differentiation  )) + 
# #   theme_classic(base_size=18)+
# #   geom_point(position=position_dodge(width=0.9))+
# #   geom_boxplot(aes(group=interaction(Differentiation,Pair.Name)),position=position_dodge(width=0.9),col=NA,scale="width",alpha=0.65)+
# #   labs(y="Cancer to Myeloid communication",x="Communication pathway"  )+
# #   scale_color_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
# #   scale_fill_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
# #   facet_wrap(~Treatment,scale="free_y",nrow=4)+theme(axis.text.x = element_text(angle=90),legend.position = "bottom")
# # 
# # ggplot(ddplot, aes(x=Pair.Name, y= scalelntransduction ,col= Differentiation,fill= Differentiation  )) + 
# #   theme_classic(base_size=18)+
# #   geom_point(position=position_dodge(width=0.9))+
# #   geom_boxplot(aes(group=interaction(Differentiation,Pair.Name)),position=position_dodge(width=0.9),col=NA,scale="width",alpha=0.65)+
# #   labs(y="Cancer to Myeloid communication",x="Communication pathway"  )+
# #   scale_color_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
# #   scale_fill_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
# #   facet_wrap(Treatment~dynamic_class3)+theme(axis.text.x = element_text(angle=90),legend.position = "bottom")
# # 
# # ggplot(ddplot[Treat=="CombinationRibo"], aes(x=Disc_V2, y= scalelntransduction ,col= Differentiation,fill= Differentiation  )) + 
# #   theme_classic(base_size=18)+
# #   geom_point(position=position_dodge(width=0.9))+
# #   geom_violin(aes(group=interaction(Differentiation,Disc_V2)),position=position_dodge(width=0.9),col=NA,scale="width",alpha=0.65)+
# #   labs(y="Cancer to Myeloid communication",x="Communication pathway"  )+
# #   scale_color_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
# #   scale_fill_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
# #   facet_wrap(~Disc_V1, scale="free_y",nrow=4)+theme(axis.text.x = element_text(angle=90),legend.position = "bottom")
# # 
# # 
# # # Cancer to M2 vs M1
# # Summar_MC <- data.table(Tcellsubset%>%
# #                           group_by(Differentiation,Pair.Name,Day,dynamic_class3,Disc_V2,Disc_V1,Patient.Study.ID,ReceptorPhenoCelltype) %>%
# #                           summarise(Signal=median(Signal),
# #                                     TransductionMu=mean(Signal) ))
# # Summar_MC[ , muD0:= sum(TransductionMu*(Day==0) )/sum((Day==0)), by=c("Patient.Study.ID" ,"dynamic_class3", "Pair.Name") ]
# # Summar_MC[,c("L","R"):= tstrsplit(Pair.Name,"_",fixed=T)]
# # Summar_MC[,ln1Day:= log(1+Day)]
# # Summar_MC[,Fact_Patient.Study.ID:= as.factor(Patient.Study.ID)]
# # 
# # countPairdd <- data.table(Summar_MC %>% group_by(Pair.Name,Differentiation) %>% summarise(n=n()) %>% group_by(Pair.Name)%>%
# #                             summarise(n=min(n))
# # )
# # Summar_MC<-Summar_MC[Pair.Name %in% countPairdd[n>10]$Pair.Name]
# # 
# # 
# # TrendsResMfromCV2diffRiboAndLetrozole <- rbindlist( lapply( unique( Summar_MC$Pair.Name )  , function(pp){
# #   cat(pp)
# #   out <- tryCatch({
# #     mm1 <- lmer( I(log(1 + TransductionMu) ) ~ Differentiation +(1|Disc_V1) + (1|Disc_V2) +  (1|Patient.Study.ID),  #(1|Disc_V1) +
# #                  data= Summar_MC[Pair.Name==pp] )
# #     if( length(coef(summary(mm1)))==0 ){ stop() }else{
# #       data.table(Treat="Both treatments",Pair.Name= pp, data.table( coef(summary(mm1)) , keep.rownames = T))
# #       
# #     }
# #   },
# #   error=function(x){
# #     data.table(Treat="Both treatments",'Pair.Name'=pp, data.table("rn"=NA, "Estimate"=NA ,'Std. Error'=NA, df =NA,'t value'=NA,'Pr(>|t|)'=NA) )
# #   })
# #   return(out)
# # }))
# # 
# # setnames(TrendsResMfromCV2diffRiboAndLetrozole,old=c("Std. Error","t value","Pr(>|t|)"),new=c("Std.Error","tval","pval"))
# # TrendsResMfromCV2diffRiboAndLetrozole[,c("L","R"):= tstrsplit(Pair.Name,"_",fixed=T)]
# # TrendsResMfromCV2diffRiboAndLetrozole$ adjpval<- p.adjust(TrendsResMfromCV2diffRiboAndLetrozole$pval)
# # TrendsResMfromCV2diffRiboAndLetrozole[rn=="DifferentiationM2"][][order(pval)][adjpval<0.05]
# # TrendsResMfromCV2diffRiboAndLetrozole[rn=="DifferentiationM2"][order(pval)][pval<0.05][1:20]
# # TrendsResMfromCV2diffRiboAndLetrozole[rn=="DifferentiationM2"][Estimate<0][order(pval)][pval<0.05]
# # TrendsResMfromCV2diffRiboAndLetrozole[rn=="DifferentiationM2"][Estimate>0][order(pval)][adjpval<0.05]
# # 
# # 
# listplotpairs <- c("INHBA_ACVR2A","INHBA_TGFBR3", "INHBA_ACVR2B","INHBA_ACVR1","SPP1_CD44","SPP1_ITGB1","CXCL9_CXCR3",
#                    "TGFBR2","TGFBR3","TGFBR1","BMPR2"
# )
# ddplot <- Tcellsubset[Pair.Name%in% listplotpairs  ]#[Pair.Name != "VTN_ITGB6"]
# ddplot[,scalelntransduction := scale( log(1 + Signal) ) , by= Pair.Name ]
# ddplot$Pair.Name <- factor( ddplot$Pair.Name, levels = unique(listplotpairs))
# ddplot[,Treatment:= "Combination ribociclib"]
# ddplot[Treat=="LetrozoleAlone",Treatment:= "Letrozole alone"]
# 
# ggplot(ddplot, aes(x=Pair.Name, y= scalelntransduction ,col= Differentiation,fill= Differentiation  )) + 
#   theme_classic(base_size=18)+
#   geom_point(position=position_dodge(width=0.9))+
#   geom_violin(aes(group=interaction(Differentiation,Pair.Name)),position=position_dodge(width=0.9),col=NA,scale="width",alpha=0.65)+
#   labs(y="Cancer to Myeloid communication",x="Communication pathway"  )+
#   scale_color_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
#   scale_fill_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
#   facet_wrap(~Treatment,scale="free_y",nrow=4)+theme(axis.text.x = element_text(angle=90),legend.position = "bottom")
# 
# 
# 
# 
# 
# 
# 
# 
# # ggsave(filename=paste0("/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cohort 2 results/Cancer communication with M2 vs M1 Macrophages.png"),width=8,height=8)
# 
# ggplot(ddplot[Pair.Name%in%c("ADAM10_AXL","CSF1_CSF1R","TGFB1_TGFBR1","ZP3_MERTK","C3_ITGAX","SPP1_CD44","FN1_SDC2")], aes(x=Pair.Name, y= scalelntransduction ,col= Differentiation,fill= Differentiation  )) + 
#   theme_classic(base_size=18)+
#   geom_point(position=position_dodge(width=0.9))+
#   geom_violin(aes(group=interaction(Differentiation,Pair.Name)),position=position_dodge(width=0.9),col=NA,scale="width",alpha=0.65)+
#   labs(y="Cancer to Myeloid communication",x="Communication pathway"  )+
#   scale_color_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
#   scale_fill_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
#   facet_wrap(~Treatment,scale="free_y",nrow=4)+theme(aspect.ratio=0.5,axis.text.x = element_text(angle=90),legend.position = "bottom")
# ggsave(filename=paste0("/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cohort 2 results/Trimmed Cancer communication with M2 vs M1 Macrophages.png"),width=8,height=8)
# 
# ggplot(ddplot, aes(x=Pair.Name, y= scalelntransduction ,col= Differentiation,fill= Differentiation  )) + 
#   theme_classic(base_size=16)+
#   geom_point(position=position_dodge(width=0.8))+
#   geom_violin(aes(group=interaction(Differentiation,Pair.Name)),position=position_dodge(width=0.8),color=NA,scale="width",alpha=0.7)+
#   labs(y="Cancer to Myeloid communications",x="Communication pathway"  )+
#   scale_color_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
#   scale_fill_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
#   theme(aspect.ratio=0.5,axis.text.x = element_text(angle=90),legend.position = "bottom")
# ggsave(filename=paste0("/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cohort 2 results/Cancer communication with M2 vs M1 Macrophages combined ribo and letrozole alone.png"),width=8,height=8)
# 
# 
# datatoorder <- data.table( ddplot%>%group_by(Pair.Name,Treat,Differentiation)%>%summarise(mu=median(scalelntransduction))%>%spread(Differentiation,mu))
# datatoorder[,diffM2_1:=M2 - M1]
# datatoorder_wide <- datatoorder%>%select(Pair.Name,Treat,diffM2_1)%>%spread(Treat,diffM2_1)
# ggplot( datatoorder_wide , aes(y= CombinationRibo, x=LetrozoleAlone)) + geom_point(size=2.5) + 
#   geom_smooth(method="lm") + theme_classic(base_size=18) + theme(aspect.ratio = 1) +
#   theme(aspect.ratio = 1) +
#   labs(y="Cancer to Myeloid communication \n in ribociclib treatment arm", 
#        x="Cancer to Myeloid communication \n in letrozole treatment arm" )
# ggsave(filename= paste0("/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cohort 2 results/Cancer communication with M2 vs M1 Macrophages comparing combined ribo and letrozole alone.png"), width= 6 , height= 6)
# 
# 
# 
# 
# unique(datatoorder[order(-diffM2_1)]$Pair.Name)
# 
# ddplot$Pair.Name <- factor(ddplot$Pair.Name, levels= unique(datatoorder[Treat=="CombinationRibo"][order(-diffM2_1)]$Pair.Name) )
# ddplot[,Treatment:="Combination ribociclib"]
# ddplot[Treat=="LetrozoleAlone",Treatment:="Letrozole alone"]
# 
# ggplot(ddplot, aes(x=Pair.Name, y= scalelntransduction ,col= Differentiation,fill= Differentiation  )) + 
#   theme_classic()+
#   geom_point(position=position_dodge(width=0.8))+
#   geom_violin(aes(group=interaction(Differentiation,Pair.Name)),col=NA,position=position_dodge(width=0.8),scale="width",alpha=0.6)+
#   labs(y="Cancer to Myeloid communications",x="Communication pathway"  )+
#   scale_color_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
#   scale_fill_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
#   facet_wrap(~Treatment,scale="free_y",nrow=4) +theme(axis.text.x = element_text(angle=90),legend.position = "bottom")
# 
# 
# 
# 
# ggplot(ddplot, aes(x=Pair.Name, y= scalelntransduction ,col= Differentiation,fill= Differentiation  )) + 
#   theme_classic()+
#   geom_point(position=position_dodge(width=0.8))+
#   geom_violin(aes(group=interaction(Differentiation,Pair.Name)),position=position_dodge(width=0.8),scale="width",alpha=0.7)+
#   labs(y="Cancer to Myeloid communications",x="Communication pathway"  )+
#   scale_color_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
#   scale_fill_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +theme(axis.text.x = element_text(angle=90))
# #ggsave(filename=paste0("/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cohort 2 results/Cancer communication with M2 vs M1 Macrophages.png"),width=8,height=8)
# 
# nn<-20#25
# 
# short.table <- rbind(data.table(ChangeInM2 = "up", TrendsResMfromCV2diffRibo[Estimate>0][order(pval)][pval<0.05][rn=="DifferentiationM2"][1:nn]),
#                      data.table(ChangeInM2 = "up",TrendsResMfromCV2diffLetro[Estimate>0][order(pval)][pval<0.05][rn=="DifferentiationM2"][1:nn]),
#                      data.table(ChangeInM2 = "down",TrendsResMfromCV2diffRibo[Estimate<0][order(pval)][pval<0.05][rn=="DifferentiationM2"][1:nn]),
#                      data.table(ChangeInM2 = "down",TrendsResMfromCV2diffLetro[Estimate<0][order(pval)][pval<0.05][rn=="DifferentiationM2"][1:nn]))[order(Estimate)]
# 
# listplotpairs<- unique( short.table$Pair.Name )
# listplotpairs<- unique( c(
#   TrendsResMfromCV2diffRibo[Estimate>0][order(pval)][pval<0.05][rn=="DifferentiationM2"][1:nn]$Pair.Name,
#   TrendsResMfromCV2diffLetro[Estimate>0][order(pval)][pval<0.05][rn=="DifferentiationM2"][1:nn]$Pair.Name,
#   TrendsResMfromCV2diffRibo[Estimate<0][order(pval)][pval<0.05][rn=="DifferentiationM2"][1:nn]$Pair.Name,
#   TrendsResMfromCV2diffLetro[Estimate<0][order(pval)][pval<0.05][rn=="DifferentiationM2"][1:nn]$Pair.Name))
# 
# ddplot <- Summar_MC[Pair.Name%in% listplotpairs  ][Pair.Name!="VTN_ITGB6"]
# ddplot[,scalelntransduction:= scale( log(1+TransductionMu) ) , by= Pair.Name ]
# ddplot$Pair.Name <- factor(ddplot$Pair.Name, levels= unique(listplotpairs))
# 
# datatoorder<-data.table( ddplot%>%group_by(Pair.Name,Treat,Differentiation)%>%summarise(mu=mean(scalelntransduction))%>%spread(Differentiation,mu))
# datatoorder[,diffM2_1:=M2 - M1]
# datatoorder_wide <- datatoorder%>%select(Pair.Name,Treat,diffM2_1)%>%spread(Treat,diffM2_1)
# ggplot( datatoorder_wide , aes(y= CombinationRibo, x=LetrozoleAlone)) + geom_point(size=2.5) + geom_smooth(method="lm") +
#   theme_classic(base_size=18) + theme(aspect.ratio = 1) +
#   labs(y="Cancer to Myeloid communication \n in ribociclib treatment arm",x="Cancer to Myeloid communication \n in letrozole treatment arm")
# 
# 
# 
# 
# 
# 
# 
# # ggplot(ddplot, aes(x=Pair.Name, y= scalelntransduction ,col= Differentiation,fill= Differentiation  )) + 
# #   theme_classic()+
# #   geom_point(position=position_dodge(width=0.8))+
# #   geom_violin(aes(group=interaction(Differentiation,Pair.Name)),position=position_dodge(width=0.8),scale="width",alpha=0.7)+
# #   labs(y="Cancer to Myeloid communications",x="Communication pathway"  )+
# #   scale_color_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
# #   scale_fill_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
# #   facet_wrap(~Treat,scale="free_y",nrow=4) #+ geom_smooth(method="lm",formula=y~-1+x,se=T)+scale_x_continuous(breaks=log(1+c(0,14,180)) , labels=c(0,14,180) )#, formula=y~s(x,k=3))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # 
# # 
# # # M2 vs M1 to T cells
# # Tcellsubset[,Differentiation:= "M1"]
# # Tcellsubset[ grep("Macrophages_Macrophages_1_1",   key_signaller) , Differentiation:= "M2"]
# # Tcellsubset[ grep("Macrophages_Macrophages_1_2",   key_signaller) , Differentiation:= "M2"]
# # Tcellsubset[ grep("Macrophages_Macrophages_2_1",   key_signaller) , Differentiation:= "M2"]
# # Tcellsubset[ grep("Macrophages_Macrophages_2_2",   key_signaller) , Differentiation:= "M2"]
# # Tcellsubset[ grep("Macrophages_Macrophages_3_1",   key_signaller) , Differentiation:= "M2"]
# # Tcellsubset[ grep("Macrophages_Macrophages_3_2",   key_signaller) , Differentiation:= "M2"]
# # Tcellsubset[ grep("Macrophages_Macrophages_4_1",   key_signaller) , Differentiation:= "M2"]
# # Tcellsubset[ grep("Macrophages_Macrophages_4_2",   key_signaller) , Differentiation:= "M2"]
# # Tcellsubset[ grep("Macrophages_Macrophages_5_1",   key_signaller) , Differentiation:= "M2"]
# # Tcellsubset[ grep("Macrophages_Macrophages_5_2",   key_signaller) , Differentiation:= "M2"]
# # 
# # Tcellsubset[ grep("Macrophages_Macrophages_1_3",   key_signaller) , Differentiation:= "M2"]
# # Tcellsubset[ grep("Macrophages_Macrophages_2_3",   key_signaller) , Differentiation:= "M2"]
# # Tcellsubset[ grep("Macrophages_Macrophages_3_3",   key_signaller) , Differentiation:= "M2"]
# # Tcellsubset[ grep("Macrophages_Macrophages_4_3",   key_signaller) , Differentiation:= "M2"]
# # Tcellsubset[ grep("Macrophages_Macrophages_5_3",   key_signaller) , Differentiation:= "M2"]
# # 
# # 
# # Summar_MT <- data.table(Tcellsubset %>%
# #                           group_by(Pair.Name,Day,dynamic_class3,Differentiation,Disc_V3, Disc_V2,Disc_V1,Patient.Study.ID) %>%
# #                           summarise(Signal= median(Signal),
# #                                     TransductionMu= mean(Signal) ))
# # Summar_MT[ , muD0:= sum(TransductionMu*(Day==0) )/sum((Day==0)), by= c("Patient.Study.ID", "dynamic_class3", "Pair.Name") ]
# # Summar_MT[,c("L","R"):= tstrsplit(Pair.Name,"_",fixed=T)]
# # Summar_MT[,ln1Day:= log(1 + Day)]
# # Summar_MT[,Fact_Patient.Study.ID:= as.factor(Patient.Study.ID)]
# # 
# # countPairdd <- data.table(Summar_MT %>% group_by(Pair.Name) %>% summarise(n=n()))
# # Summar_MT<-Summar_MT[Pair.Name %in% countPairdd[n>10]$Pair.Name]
# # 
# # ggplot(Tcellsubset[grep("IL2", Pair.Name)  ], aes(x=Differentiation, y=log(1+Signal) ,col=Differentiation,fill=Differentiation  )) + #,linetype=as.factor(6-Disc_V2)
# #   theme_classic()+# linetype=as.factor(6-Disc_V2)
# #   geom_point(position=position_dodge(width=0.5))+
# #   geom_violin(aes(group=interaction(Differentiation)),position=position_dodge(width=0.5),scale="width",alpha=0.7)+
# #   labs(y="Cancer to Myeloid communications \n stimulating M2 differentiation",x="Myeloid differentiation"  )+
# #   scale_color_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
# #   scale_fill_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
# #   #scale_color_npg(name="Tumor response" ,labels=c("Resistant","Sensitive"))+
# #   #scale_fill_npg(name="Tumor response" ,labels=c("Resistant","Sensitive"))+
# #   facet_wrap(~Pair.Name,scale="free_y",nrow=4) #+ geom_smooth(method="lm",formula=y~-1+x,se=T)+scale_x_continuous(breaks=log(1+c(0,14,180)) , labels=c(0,14,180) )#, formula=y~s(x,k=3))
# # 
# # ggplot(Tcellsubset[grep("IFN", Pair.Name)  ], aes(x=Differentiation, y=log(1+Signal) ,col=Differentiation,fill=Differentiation  )) + #,linetype=as.factor(6-Disc_V2)
# #   theme_classic()+# linetype=as.factor(6-Disc_V2)
# #   geom_point(position=position_dodge(width=0.5))+
# #   geom_violin(aes(group=interaction(Differentiation)),position=position_dodge(width=0.5),scale="width",alpha=0.7)+
# #   labs(y="Cancer to Myeloid communications \n stimulating M2 differentiation",x="Myeloid differentiation"  )+
# #   scale_color_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
# #   scale_fill_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
# #   #scale_color_npg(name="Tumor response" ,labels=c("Resistant","Sensitive"))+
# #   #scale_fill_npg(name="Tumor response" ,labels=c("Resistant","Sensitive"))+
# #   facet_wrap(~Pair.Name,scale="free_y",nrow=4) #+ geom_smooth(method="lm",formula=y~-1+x,se=T)+scale_x_continuous(breaks=log(1+c(0,14,180)) , labels=c(0,14,180) )#, formula=y~s(x,k=3))
# # 
# # 
# # ggplot(Tcellsubset[grep("IL18", Pair.Name)  ][Disc_V2=="1"][Disc_V3=="1"], aes(x=Differentiation, y=log(1+Signal) ,col=Differentiation,fill=Differentiation  )) + #,linetype=as.factor(6-Disc_V2)
# #   theme_classic()+# linetype=as.factor(6-Disc_V2)
# #   geom_point(position=position_dodge(width=0.5))+
# #   geom_violin(aes(group=interaction(Differentiation)),position=position_dodge(width=0.5),scale="width",alpha=0.7)+
# #   labs(y="Cancer to Myeloid communications \n stimulating M2 differentiation",x="Myeloid differentiation"  )+
# #   scale_color_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
# #   scale_fill_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
# #   #scale_color_npg(name="Tumor response" ,labels=c("Resistant","Sensitive"))+
# #   #scale_fill_npg(name="Tumor response" ,labels=c("Resistant","Sensitive"))+
# #   facet_wrap(~Pair.Name,scale="free_y",nrow=4) #+ geom_smooth(method="lm",formula=y~-1+x,se=T)+scale_x_continuous(breaks=log(1+c(0,14,180)) , labels=c(0,14,180) )#, formula=y~s(x,k=3))
# # 
# # 
# # ggplot(Summar_MT[grep("IL1B", Pair.Name)  ], aes(x=Differentiation, y=log(1+Signal) ,col=Differentiation,fill=Differentiation  )) + #,linetype=as.factor(6-Disc_V2)
# #   theme_classic()+# linetype=as.factor(6-Disc_V2)
# #   geom_point(position=position_dodge(width=0.5))+
# #   geom_violin(aes(group=interaction(Differentiation)),position=position_dodge(width=0.5),scale="width",alpha=0.7)+
# #   labs(y="Cancer to Myeloid communications \n stimulating M2 differentiation",x="Myeloid differentiation"  )+
# #   scale_color_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
# #   scale_fill_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
# #   #scale_color_npg(name="Tumor response" ,labels=c("Resistant","Sensitive"))+
# #   #scale_fill_npg(name="Tumor response" ,labels=c("Resistant","Sensitive"))+
# #   facet_wrap(~Pair.Name,scale="free_y",nrow=4) #+ geom_smooth(method="lm",formula=y~-1+x,se=T)+scale_x_continuous(breaks=log(1+c(0,14,180)) , labels=c(0,14,180) )#, formula=y~s(x,k=3))
# # 
# # ggplot(Summar_MT[Pair.Name%in% "CXCL10_CXCR3"
# #                  ], aes(x=Differentiation, y=log(1+TransductionMu) ,col=Differentiation,fill=Differentiation  )) + #,linetype=as.factor(6-Disc_V2)
# #   theme_classic()+# linetype=as.factor(6-Disc_V2)
# #   geom_point(position=position_dodge(width=0.5))+
# #   geom_violin(aes(group=interaction(Differentiation)),position=position_dodge(width=0.5),scale="width",alpha=0.7)+
# #   labs(y="Cancer to Myeloid communications \n stimulating M2 differentiation",x="Myeloid differentiation"  )+
# #   scale_color_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
# #   scale_fill_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
# #   #scale_color_npg(name="Tumor response" ,labels=c("Resistant","Sensitive"))+
# #   #scale_fill_npg(name="Tumor response" ,labels=c("Resistant","Sensitive"))+
# #   facet_wrap(~Pair.Name,scale="free_y",nrow=4) #+ geom_smooth(method="lm",formula=y~-1+x,se=T)+scale_x_continuous(breaks=log(1+c(0,14,180)) , labels=c(0,14,180) )#, formula=y~s(x,k=3))
# # 
# # TrendsResTfromMV2diff <- rbindlist( lapply( unique( Tcellsubset$Pair.Name )  , function(pp){
# #   cat(pp)
# #   out <- tryCatch({
# #     mm1 <- lmer( I(log(1 + Signal) ) ~ Differentiation + #(1|Disc_V1) + 
# #                    # mm1 <- lmer( I(log(1 + TransductionMu) ) ~ Differentiation + #(1|Disc_V1) + 
# #                    (1|Patient.Study.ID), 
# #                  data= Tcellsubset[Pair.Name==pp] )
# #     if( length(coef(summary(mm1)))==0 ){ stop() }else{
# #       data.table(Pair.Name= pp, data.table( coef(summary(mm1)) , keep.rownames = T))
# #     }
# #   },
# #   error=function(x){
# #     data.table('Pair.Name'=pp, data.table("rn"=NA, "Estimate"=NA ,'Std. Error'=NA, df =NA,'t value'=NA,'Pr(>|t|)'=NA) )
# #   })
# #   return(out)
# # }))
# # 
# # setnames(TrendsResTfromMV2diff,old=c("Std. Error","t value","Pr(>|t|)"),new=c("Std.Error","tval","pval"))
# # TrendsResTfromMV2diff[,c("L","R"):= tstrsplit(Pair.Name,"_",fixed=T)]
# # TrendsResTfromMV2diff$ adjpval<- p.adjust(TrendsResTfromMV2diff$pval)
# # 
# # TrendsResTfromMV2diff[][order(pval)][adjpval<0.05][rn=="DifferentiationM2"]
# # TrendsResTfromMV2diff[][order(pval)][pval<0.05][rn=="DifferentiationM2"][1:20]
# # TrendsResTfromMV2diff[Estimate<0][order(pval)][pval<0.05][rn=="DifferentiationM2"]
# # 
# # 
# # ggplot(Summar_MT[grepl("IL15",R)#Pair.Name%in% "ADAM17_ERBB4"
# #                  ], aes(x=Differentiation, y=log(1+TransductionMu) ,col=Differentiation,fill=Differentiation  )) + #,linetype=as.factor(6-Disc_V2)
# #   theme_classic()+# linetype=as.factor(6-Disc_V2)
# #   geom_point(position=position_dodge(width=0.5))+
# #   geom_violin(aes(group=interaction(Differentiation)),position=position_dodge(width=0.5),scale="width",alpha=0.7)+
# #   labs(y="Cancer to Myeloid communications \n stimulating M2 differentiation",x="Myeloid differentiation"  )+
# #   scale_color_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
# #   scale_fill_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
# #   #scale_color_npg(name="Tumor response" ,labels=c("Resistant","Sensitive"))+
# #   #scale_fill_npg(name="Tumor response" ,labels=c("Resistant","Sensitive"))+
# #   facet_wrap(~Pair.Name,scale="free_y",nrow=4) #+ geom_smooth(method="lm",formula=y~-1+x,se=T)+scale_x_continuous(breaks=log(1+c(0,14,180)) , labels=c(0,14,180) )#, formula=y~s(x,k=3))
# # 
# # 
# # ggplot(Summar_MT[Pair.Name%in% TrendsResTfromMV2diff[][order(pval)][pval<0.05][rn=="DifferentiationM2"][Estimate<0]$Pair.Name
# #                  ], aes(x=Differentiation, y=log(1+TransductionMu) ,col=Differentiation,fill=Differentiation  )) + #,linetype=as.factor(6-Disc_V2)
# #   theme_classic()+# linetype=as.factor(6-Disc_V2)
# #   theme(aspect.ratio=1)+
# #   geom_point(position=position_dodge(width=0.5))+
# #   geom_violin(aes(group=interaction(Differentiation)),position=position_dodge(width=0.5),scale="width",alpha=0.7)+
# #   labs(y="Cancer to Myeloid communications \n stimulating M2 differentiation",x="Myeloid differentiation"  )+
# #   scale_color_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
# #   scale_fill_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
# #   #scale_color_npg(name="Tumor response" ,labels=c("Resistant","Sensitive"))+
# #   #scale_fill_npg(name="Tumor response" ,labels=c("Resistant","Sensitive"))+
# #   facet_wrap(~Pair.Name,scale="free_y",nrow=2) #+ geom_smooth(method="lm",formula=y~-1+x,se=T)+scale_x_continuous(breaks=log(1+c(0,14,180)) , labels=c(0,14,180) )#, formula=y~s(x,k=3))
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # TrendsRes<-rbindlist( lapply( unique( Summar_MC$Pair.Name )  , function(pp){
# #   out <- tryCatch({
# #     mm1 <- lmer( I(log(1+Signal)-log(1+muD0))~ -1 + ln1Day + ln1Day:dynamic_class3 +(1|Disc_V1) +(1|Disc_V2), 
# #                  data= Summar_MC[Pair.Name==pp] )
# #     if( length(coef(summary(mm1)))==0 ){ stop() }else{
# #       data.table(Pair.Name=pp, data.table( coef(summary(mm1)) ,keep.rownames = T))
# #     }
# #   },
# #   error=function(x){
# #     data.table('Pair.Name'=pp, data.table("rn"=NA, "Estimate"=NA ,'Std. Error'=NA, df =NA,'t value'=NA,'Pr(>|t|)'=NA) )
# #   })
# #   return(out)
# # }))
# # # TrendsRes<- na.omit(TrendsRes)
# # # setnames(TrendsRes,old=c("Std. Error","t value","Pr(>|t|)"),new=c("Std.Error","tval","pval"))
# # # TrendsRes[,c("L","R"):= tstrsplit(Pair.Name,"_",fixed=T)]
# # # TrendsRes$ adjpval<- p.adjust(TrendsRes$pval)
# # # TrendsRes[][order(pval)][adjpval<0.05][rn=="ln1Day:dynamic_class3Non-response"][1:20]#$Pair.Name
# # # TrendsRes[grepl("TGFB",L)&grepl("TGFB",R)][order(pval)]
# # 
# # 
# # 
# # 
# # 
# # Summar_CfromM <- data.table(Tcellsubset%>%
# #                               group_by(Pair.Name,Day,dynamic_class3,Disc_V2,Disc_V1,Patient.Study.ID)%>%
# #                               summarise(TransductionMu=mean(Signal) ))
# # Summar_CfromM[ , muD0:=sum(TransductionMu*(Day==0) )/sum((Day==0)), by=c("Patient.Study.ID","dynamic_class3","Pair.Name") ]
# # Summar_CfromM[,c("L","R"):= tstrsplit(Pair.Name,"_",fixed=T)]
# # Summar_CfromM[,ln1Day:=log(1+Day)]
# # Summar_CfromM[,Fact_Patient.Study.ID:=as.factor(Patient.Study.ID)]
# # Summar_CfromM[,Differentiation:= "a non mycaf"]
# # Summar_CfromM[Disc_V2==1, Differentiation:= "mycaf differentiated"]
# # countPairdd<-data.table(Summar_CfromM%>%group_by(Pair.Name)%>%summarise(n=n()))
# # Summar_CfromM<-Summar_CfromM[Pair.Name%in%countPairdd[n>20]$Pair.Name]
# # 
# # 
# # 
# # 
# # library(viridis)
# # 
# # ggplot(Summar_CfromM[grepl("ERBB",R)#Pair.Name%in% "ADAM17_ERBB4"
# #                      ], aes(x=Differentiation, y=log(1+TransductionMu) ,col=Differentiation,fill=Differentiation  )) + #,linetype=as.factor(6-Disc_V2)
# #   theme_classic()+# linetype=as.factor(6-Disc_V2)
# #   geom_point(position=position_dodge(width=0.5))+
# #   geom_violin(aes(group=interaction(Differentiation)),position=position_dodge(width=0.5),scale="width",alpha=0.7)+
# #   labs(y="Cancer to fibroblast communications \n stimulating myCAF differentiation",x="myCAF differentiation"  )+
# #   scale_color_viridis(name="Fibroblast type",labels=c("inactive","myCAF"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
# #   scale_fill_viridis(name="Fibroblast type",labels=c("inactive","myCAF"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
# #   #scale_color_npg(name="Tumor response" ,labels=c("Resistant","Sensitive"))+
# #   #scale_fill_npg(name="Tumor response" ,labels=c("Resistant","Sensitive"))+
# #   facet_wrap(~Pair.Name,scale="free_y",nrow=4) #+ geom_smooth(method="lm",formula=y~-1+x,se=T)+scale_x_continuous(breaks=log(1+c(0,14,180)) , labels=c(0,14,180) )#, formula=y~s(x,k=3))
# # 
# # 
# # # TrendsRes<-rbindlist( lapply( unique( Summar_MC$Pair.Name )  , function(pp){
# # #   out <- tryCatch({
# # #     mm1 <- lmer( I(log(1+Signal)-log(1+muD0))~ -1 + ln1Day + ln1Day:dynamic_class3 +(1|Disc_V1) +(1|Disc_V2), 
# # #                  data= Summar_MC[Pair.Name==pp] )
# # #     if( length(coef(summary(mm1)))==0 ){ stop() }else{
# # #       data.table(Pair.Name=pp, data.table( coef(summary(mm1)) ,keep.rownames = T))
# # #     }
# # #   },
# # #   error=function(x){
# # #     data.table('Pair.Name'=pp, data.table("rn"=NA, "Estimate"=NA ,'Std. Error'=NA, df =NA,'t value'=NA,'Pr(>|t|)'=NA) )
# # #   })
# # #   return(out)
# # # }))
# # # TrendsRes<- na.omit(TrendsRes)
# # # setnames(TrendsRes,old=c("Std. Error","t value","Pr(>|t|)"),new=c("Std.Error","tval","pval"))
# # # TrendsRes[,c("L","R"):= tstrsplit(Pair.Name,"_",fixed=T)]
# # # TrendsRes$ adjpval<- p.adjust(TrendsRes$pval)
# # # TrendsRes[][order(pval)][adjpval<0.05][rn=="ln1Day:dynamic_class3Non-response"][1:20]#$Pair.Name
# # # TrendsRes[grepl("TGFB",L)&grepl("TGFB",R)][order(pval)]
# # # ggplot(Summar_MC[Pair.Name%in%TrendsRes[Pair.Name!="L1CAM_EGFR"][order(pval)][adjpval<0.05][rn=="ln1Day:dynamic_class3Non-response"][Estimate>0][1:20]$Pair.Name], aes(x=log(1+Day), y=log(1+Signal)-log(1+muD0), col=dynamic_class3   )) + #,linetype=as.factor(6-Disc_V2)
# # #   theme_classic()+# linetype=as.factor(6-Disc_V2)
# # #   geom_violin(aes(group=interaction(Day,dynamic_class3),fill=dynamic_class3),position=position_dodge(width=2.5),scale="width",alpha=0.7)+
# # #   geom_point(aes(shape=),position=position_dodge(width=2.5))+
# # #   labs(y="Cancer to fibroblast communications \n stimulating myCAF differentiation",x="Day"  )+
# # #   scale_color_npg(name="Tumor response" ,labels=c("Resistant","Sensitive"))+
# # #   scale_fill_npg(name="Tumor response" ,labels=c("Resistant","Sensitive"))+
# # #   facet_wrap(~Pair.Name,scale="free_y",nrow=4) + geom_smooth(method="lm",formula=y~-1+x,se=T)+scale_x_continuous(breaks=log(1+c(0,14,180)) , labels=c(0,14,180) )#, formula=y~s(x,k=3))
# # 
# # TrendsResCfromFV2diff <- rbindlist( lapply( unique( Summar_CfromM$Pair.Name )  , function(pp){
# #   out <- tryCatch({
# #     # mm1 <- lmer( I(log(1+TransductionMu)-log(1+muD0))~ -1 + ln1Day + ln1Day:dynamic_class3 + (1|Disc_V1) + (1|Disc_V2), 
# #     # mm1 <- lmer( I(log(1+TransductionMu) )~ -1 + ln1Day + ln1Day:dynamic_class3 + (1|Disc_V1) + (1|Disc_V2), 
# #     # mm1 <- lmer( I(log(1+TransductionMu) )~ Disc_V2 + (1|Disc_V1) +(1|Patient.Study.ID), 
# #     mm1 <- lmer( I(log(1 + TransductionMu) ) ~ Differentiation + (1|Disc_V1) + (1|Patient.Study.ID), 
# #                  data= Summar_CfromM[Pair.Name==pp] )
# #     if( length(coef(summary(mm1)))==0 ){ stop() }else{
# #       data.table(Pair.Name= pp, data.table( coef(summary(mm1)) , keep.rownames = T))
# #     }
# #   },
# #   error=function(x){
# #     data.table('Pair.Name'=pp, data.table("rn"=NA, "Estimate"=NA ,'Std. Error'=NA, df =NA,'t value'=NA,'Pr(>|t|)'=NA) )
# #   })
# #   return(out)
# # }))
# # TrendsResCfromFV3diff <- na.omit(TrendsResCfromFV3diff)
# # TrendsResCfromFV3diff <-TrendsResCfromFV3diff[rn!="(Intercept)"]
# # setnames(TrendsResCfromFV3diff,old=c("Std. Error","t value","Pr(>|t|)"),new=c("Std.Error","tval","pval"))
# # TrendsResCfromFV3diff[,c("L","R"):= tstrsplit(Pair.Name,"_",fixed=T)]
# # TrendsResCfromFV3diff$ adjpval<- p.adjust(TrendsResCfromFV3diff$pval)
# # 
# # TrendsResCfromFV3diff[Estimate>0][order(pval)][pval<0.05][1:20]
# # 
# # ggplot(Summar_CfromM[Pair.Name %in% countPairdd[n>=10]$Pair.Name ][ Pair.Name %in% TrendsResCfromFV3diff[Estimate>0][order(pval)][pval<0.05][1:20] $Pair.Name #c("EGF_ERBB4","BMP4_BMPR1B")
# #                                                                     ], aes(x=Differentiation, y=log(1+TransductionMu)   )) + #,linetype=as.factor(6-Disc_V2)
# #   theme_classic()+# linetype=as.factor(6-Disc_V2)
# #   geom_point(position=position_dodge(width=0.5))+
# #   geom_violin(aes(group=interaction(Differentiation)),position=position_dodge(width=0.5),scale="width",alpha=0.7)+
# #   labs(y="Cancer to fibroblast communications \n stimulating myCAF differentiation",x="myCAF differentiation"  )+
# #   #scale_color_npg(name="Tumor response" ,labels=c("Resistant","Sensitive"))+
# #   #scale_fill_npg(name="Tumor response" ,labels=c("Resistant","Sensitive"))+
# #   facet_wrap(~Pair.Name,scale="free_y",nrow=4) #+ geom_smooth(method="lm",formula=y~-1+x,se=T)+scale_x_continuous(breaks=log(1+c(0,14,180)) , labels=c(0,14,180) )#, formula=y~s(x,k=3))
# # 
# # ggplot(Summar_CfromM[Pair.Name %in% countPairdd[n>=10]$Pair.Name ][ Pair.Name %in% TrendsResCfromFV3diff[Estimate>0][order(pval)][pval<0.05][1:10] $Pair.Name #c("EGF_ERBB4","BMP4_BMPR1B")
# #                                                                     ], aes(x=Differentiation, y=log(1+TransductionMu) ,col=Differentiation,fill=Differentiation  )) + #,linetype=as.factor(6-Disc_V2)
# #   theme_classic()+# linetype=as.factor(6-Disc_V2)
# #   geom_point(position=position_dodge(width=0.5))+
# #   geom_violin(aes(group=interaction(Differentiation)),position=position_dodge(width=0.5),scale="width",alpha=0.7)+
# #   labs(y="Cancer to fibroblast communications \n stimulating myCAF differentiation",x="myCAF differentiation"  )+
# #   scale_color_viridis(name="Fibroblast type",labels=c("inactive","myCAF"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
# #   scale_fill_viridis(name="Fibroblast type",labels=c("inactive","myCAF"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
# #   #scale_color_npg(name="Tumor response" ,labels=c("Resistant","Sensitive"))+
# #   #scale_fill_npg(name="Tumor response" ,labels=c("Resistant","Sensitive"))+
# #   facet_wrap(~Pair.Name,scale="free_y",nrow=2) #+ geom_smooth(method="lm",formula=y~-1+x,se=T)+scale_x_continuous(breaks=log(1+c(0,14,180)) , labels=c(0,14,180) )#, formula=y~s(x,k=3))
# # 
# # 
# # plt0<- Summar_CfromM[Pair.Name %in% countPairdd[n>=10]$Pair.Name ][ Pair.Name %in% TrendsResCfromFV3diff[Estimate>0][order(pval)][pval<0.05][1:10] $Pair.Name ]
# # plt0[,scaledcomm:=scale(log(1+TransductionMu)) ,by=c("Pair.Name")]
# # plt0$Pair.Name <- factor(plt0$Pair.Name , 
# #                          #levels=c("ADAM17_ERBB4","EGF_ERBB4", "BMP4_BMPR1B" , "BMP7_BMPR1B", "NPNT_ITGA8","FN1_ITGA8",  "SEMA6D_KDR",  "ANGPT1_TIE1",      "HLA-A_LILRB2", "HLA-E_KLRC1")
# #                          levels=c("FN1_ITGA8", "NPNT_ITGA8", "BMP4_BMPR1B" , "BMP7_BMPR1B", "SEMA6D_KDR",  "ADAM17_ERBB4","EGF_ERBB4",  "ANGPT1_TIE1",      "HLA-A_LILRB2", "HLA-E_KLRC1")
# # )
# # ggplot( plt0, aes(x=Pair.Name, y=scaledcomm ,col=Differentiation,fill=Differentiation  )) + #,linetype=as.factor(6-Disc_V2)
# #   theme_classic()+# linetype=as.factor(6-Disc_V2)
# #   geom_point(position=position_dodge(width=0.95))+
# #   geom_violin(aes(group=interaction(Differentiation,Pair.Name)),position=position_dodge(width=0.95),scale="width",alpha=0.7)+
# #   labs(y="Fibroblast to cancer communications",x="Communication pathways"  )+
# #   scale_color_viridis(name="Fibroblast type",labels=c("inactive","myCAF"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
# #   scale_fill_viridis(name="Fibroblast type",labels=c("inactive","myCAF"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
# #   theme(axis.text.x=element_text(angle = 90))
# # #scale_color_npg(name="Tumor response" ,labels=c("Resistant","Sensitive"))+
# # #scale_fill_npg(name="Tumor response" ,labels=c("Resistant","Sensitive"))+
# # #+ geom_smooth(method="lm",formula=y~-1+x,se=T)+scale_x_continuous(breaks=log(1+c(0,14,180)) , labels=c(0,14,180) )#, formula=y~s(x,k=3))
# # 
# # save(plt0 ,Summar_CfromM,TrendsResCfromFV3diff ,cancersubset ,file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/mycaf vs inactive Fibroblast communication with cancer/mycaf vs inactive Fibroblast communication with cancer.RData")
# # load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/mycaf vs inactive Fibroblast communication with cancer/mycaf vs inactive Fibroblast communication with cancer.RData")
# # 
# # ggplot(Summar_CfromM[grepl("ERBB",R)
# #                      ], aes(x=Differentiation, y=log(1+TransductionMu) ,col=Differentiation,fill=Differentiation  )) + #,linetype=as.factor(6-Disc_V2)
# #   theme_classic()+# linetype=as.factor(6-Disc_V2)
# #   geom_point(position=position_dodge(width=0.5))+
# #   geom_violin(aes(group=interaction(Differentiation)),position=position_dodge(width=0.5),scale="width",alpha=0.7)+
# #   labs(y="Cancer to fibroblast communications \n stimulating myCAF differentiation",x="myCAF differentiation"  )+
# #   scale_color_viridis(name="Fibroblast type",labels=c("inactive","myCAF"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
# #   scale_fill_viridis(name="Fibroblast type",labels=c("inactive","myCAF"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
# #   #scale_color_npg(name="Tumor response" ,labels=c("Resistant","Sensitive"))+
# #   #scale_fill_npg(name="Tumor response" ,labels=c("Resistant","Sensitive"))+
# #   facet_wrap(~Pair.Name,scale="free_y",nrow=4) #+ geom_smooth(method="lm",formula=y~-1+x,se=T)+scale_x_continuous(breaks=log(1+c(0,14,180)) , labels=c(0,14,180) )#, formula=y~s(x,k=3))
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # ggplot(Summar_CfromM[ Pair.Name %in% TrendsResCfromFV3diff[Estimate>0][order(pval)][pval<0.05][1:20] $Pair.Name #c("EGF_ERBB4","BMP4_BMPR1B")
# #                       ], aes(x=Differentiation, y=log(1+TransductionMu) ,col=dynamic_class3 ,fill=dynamic_class3  )) + #,linetype=as.factor(6-Disc_V2)
# #   theme_classic()+# linetype=as.factor(6-Disc_V2)
# #   geom_point(position=position_dodge(width=0.5))+
# #   geom_violin(aes(group=interaction(dynamic_class3,Differentiation)),position=position_dodge(width=0.5),scale="width",alpha=0.7)+
# #   labs(y="Cancer to fibroblast communications \n stimulating myCAF differentiation",x="myCAF differentiation"  )+
# #   #scale_color_npg(name="Tumor response" ,labels=c("Resistant","Sensitive"))+
# #   #scale_fill_npg(name="Tumor response" ,labels=c("Resistant","Sensitive"))+
# #   facet_wrap(~Pair.Name,scale="free_y",nrow=4) #+ geom_smooth(method="lm",formula=y~-1+x,se=T)+scale_x_continuous(breaks=log(1+c(0,14,180)) , labels=c(0,14,180) )#, formula=y~s(x,k=3))
# # 
# # 
# # 
# # ggplot(Summar_CfromM[Pair.Name%in% "PTN_ALK"
# #                      ], aes(x=(6-Disc_V2), y=log(1+TransductionMu)   )) + #,linetype=as.factor(6-Disc_V2)
# #   theme_classic()+# linetype=as.factor(6-Disc_V2)
# #   geom_point(position=position_dodge(width=0.5))+
# #   geom_violin(aes(group=interaction(Disc_V2)),position=position_dodge(width=0.5),scale="width",alpha=0.7)+
# #   labs(y="Cancer to fibroblast communications \n stimulating myCAF differentiation",x="myCAF differentiation"  )+
# #   #scale_color_npg(name="Tumor response" ,labels=c("Resistant","Sensitive"))+
# #   #scale_fill_npg(name="Tumor response" ,labels=c("Resistant","Sensitive"))+
# #   facet_wrap(dynamic_class3~Day,scale="free_y",nrow=3) #+ geom_smooth(method="lm",formula=y~-1+x,se=T)+scale_x_continuous(breaks=log(1+c(0,14,180)) , labels=c(0,14,180) )#, formula=y~s(x,k=3))
# # 
# # ggplot(Summar_CfromM[Pair.Name%in% TrendsResCfromFV3diff[][order(pval)][pval<0.05][1:20]$Pair.Name
# #                      ], aes(x=(6-Disc_V2), y=log(1+TransductionMu)   )) + #,linetype=as.factor(6-Disc_V2)
# #   theme_classic()+# linetype=as.factor(6-Disc_V2)
# #   geom_point(,position=position_dodge(width=0.5))+
# #   geom_violin(aes(group=interaction(Disc_V2)),position=position_dodge(width=0.5),scale="width",alpha=0.7)+
# #   labs(y="Cancer to fibroblast communications \n stimulating myCAF differentiation",x="myCAF differentiation"  )+
# #   #scale_color_npg(name="Tumor response" ,labels=c("Resistant","Sensitive"))+
# #   #scale_fill_npg(name="Tumor response" ,labels=c("Resistant","Sensitive"))+
# #   facet_wrap(~Pair.Name,scale="free_y",nrow=3) #+ geom_smooth(method="lm",formula=y~-1+x,se=T)+scale_x_continuous(breaks=log(1+c(0,14,180)) , labels=c(0,14,180) )#, formula=y~s(x,k=3))
# # 
# # 
# # TrendsResCfromFV3diff[][order(pval)][adjpval<0.05][rn=="ln1Day:dynamic_class3Non-response"][1:20]#$Pair.Name
# # TrendsResCfromFV3diff[grepl("TGFB",L)&grepl("TGFB",R)][order(pval)]
# # 
# # ggplot(Summar_CfromM[#grepl("TGFB1",L)&grepl("TGFBR1",R)
# #   grepl("ERBB",Pair.Name)#Pair.Name=="EGF_EGFR"
# #   #Pair.Name=="ZP3_MERTK"
# #   ], aes(x=(6-Disc_V2), y=log(1+TransductionMu), col=dynamic_class3   )) + #,linetype=as.factor(6-Disc_V2)
# #   theme_classic()+# linetype=as.factor(6-Disc_V2)
# #   geom_point(,position=position_dodge(width=0.5))+
# #   geom_violin(aes(group=interaction(Disc_V2,dynamic_class3),fill=dynamic_class3),position=position_dodge(width=0.5),scale="width",alpha=0.7)+
# #   labs(y="Cancer to fibroblast communications \n stimulating myCAF differentiation",x="Day"  )+
# #   scale_color_npg(name="Tumor response" ,labels=c("Resistant","Sensitive"))+
# #   scale_fill_npg(name="Tumor response" ,labels=c("Resistant","Sensitive"))+
# #   facet_wrap(~Pair.Name,scale="free_y",nrow=3) #+ geom_smooth(method="lm",formula=y~-1+x,se=T)+scale_x_continuous(breaks=log(1+c(0,14,180)) , labels=c(0,14,180) )#, formula=y~s(x,k=3))
# # 
# # ggplot(Summar_CfromM[grepl("TGF",L)&grepl("TGF",R)], aes(x=Pair.Name, y=log(1+TransductionMu), col=dynamic_class3   )) + #,linetype=as.factor(6-Disc_V2)
# #   theme_classic()+# linetype=as.factor(6-Disc_V2)
# #   geom_violin(aes(group=interaction(Pair.Name,dynamic_class3),fill=dynamic_class3),position=position_dodge(width=2.5),scale="width",alpha=0.7)+
# #   geom_point(,position=position_dodge(width=2.5))+
# #   labs(y="Cancer to fibroblast communications \n stimulating myCAF differentiation",x="Day"  )+
# #   scale_color_npg(name="Tumor response" ,labels=c("Resistant","Sensitive"))+
# #   scale_fill_npg(name="Tumor response" ,labels=c("Resistant","Sensitive"))+
# #   facet_wrap(~Day,scale="free_y",nrow=3) #+ geom_smooth(method="lm",formula=y~-1+x,se=T)+scale_x_continuous(breaks=log(1+c(0,14,180)) , labels=c(0,14,180) )#, formula=y~s(x,k=3))
# # 
# # 
# # 
# # ggplot(Summar_CfromM[Pair.Name=="PTN_ALK"], aes(x=Pair.Name, y=log(1+TransductionMu), col=dynamic_class3   )) + #,linetype=as.factor(6-Disc_V2)
# #   theme_classic()+# linetype=as.factor(6-Disc_V2)
# #   geom_violin(aes(group=interaction(Pair.Name,dynamic_class3),fill=dynamic_class3), position=position_dodge(width=2.5),scale="width",alpha=0.7)+
# #   geom_point(,position=position_dodge(width=2.5))+
# #   labs(y="Cancer to fibroblast communications \n stimulating myCAF differentiation",x="Day"  )+
# #   scale_color_npg(name="Tumor response" ,labels=c("Resistant","Sensitive"))+
# #   scale_fill_npg(name="Tumor response" ,labels=c("Resistant","Sensitive"))+
# #   facet_wrap(~Day,scale="free_y",nrow=3) #+ geom_smooth(method="lm",formula=y~-1+x,se=T)+scale_x_continuous(breaks=log(1+c(0,14,180)) , labels=c(0,14,180) )#, formula=y~s(x,k=3))
# # 
# # # 
# # # 
# # # 
# # # ggplot(macrosubset[grepl("TGFBR",Pair.Name)][Disc_V1<3]%>%group_by(Pair.Name,Day,dynamic_class3,Disc_V2,Patient.Study.ID)%>%summarise(Signal=median(Signal)), aes(x=log(1+Day), y=log(1+Signal), col=dynamic_class3 ,linetype=as.factor(6-Disc_V2)  )) +theme_classic()+# linetype=as.factor(6-Disc_V2)
# # #   geom_point(position=position_dodge(width=1.5))+geom_violin(aes(group=interaction(Day,dynamic_class3)),position=position_dodge(width=1.5),scale="width")+
# # #   facet_wrap(~Pair.Name,scale="free_y",nrow=2) + geom_smooth(method="lm",se=F)#, formula=y~s(x,k=3))
# # # 
# # # ggplot(macrosubset[grepl("TGFBR",Pair.Name)][Disc_V1<3]%>%group_by(Pair.Name,Day,dynamic_class3,Disc_V2,Patient.Study.ID)%>%summarise(Signal=median(Signal)), aes(x=log(1+Day), y=log(1+Signal), col=dynamic_class3 ,linetype=as.factor(6-Disc_V2)  ))+# linetype=as.factor(6-Disc_V2)
# # #   geom_point(position=position_dodge(width=2.5))+
# # #   geom_violin(aes(group=interaction(Day,dynamic_class3),fill=dynamic_class3),alpha=0.9,position=position_dodge(width=2.5),scale="width")+
# # #   facet_wrap(~Pair.Name,scale="free_y",nrow=2) +theme_classic()
# # # 
# # # 
# # # ggplot(macrosubset[grepl("TGFBR",Pair.Name)][]%>%group_by(Day,Pair.Name,dynamic_class3,Disc_V2,Patient.Study.ID)%>%summarise(Signal=median(Signal)), aes(x=log(1+Day), y=log(1+Signal), col=dynamic_class3   ))+# linetype=as.factor(6-Disc_V2)
# # #   #geom_point()+
# # #   facet_wrap(~Pair.Name,scale="free_y",nrow=2) + geom_smooth(method="lm",se=F)#, formula=y~s(x,k=3))
# # # 
# # # ggplot(macrosubset[grepl("TGFBR",Pair.Name)], aes(y=6-Disc_V2,x=Disc_V1, fill=log(1+Signal) ))+  geom_tile()+
# # #   facet_wrap(~Pair.Name,scale="free_y")
# # # 
# # # 
# # # ggplot(macrosubset[Pair.Name%in%slopes[Estimate>0][order(pval)]$Pair.Name[1:10]], aes(x=6-Disc_V2, y=log(1+Signal)))+geom_point()+geom_violin(aes(group=Disc_V2),scale="width")+
# # #   facet_wrap(~Pair.Name,scale="free_y")+geom_smooth(method="lm")+#, formula=y~s(x,k=3))
# # #   labs(x="myCAF differentiation",y="Cancer to fibroblast \n communication")
# # # 
# # # 
# # # ggplot(macrosubset[Pair.Name%in%slopes[Estimate>0][order(pval)]$Pair.Name[1:10]], aes(x=6-Disc_V2, y=log(1+Signal),group=Disc_V1,col=Disc_V1))+geom_point()+
# # #   facet_wrap(~Pair.Name,scale="free_y")+geom_smooth(method="lm")+#, formula=y~s(x,k=3))
# # #   labs(x="myCAF differentiation",y="Cancer to fibroblast \n communication")
# # # ggplot(macrosubset[Pair.Name%in%slopes[Estimate>0][order(pval)]$Pair.Name[1:10]], aes(x=6-Disc_V2, y=log(1+Signal) ))+geom_point()+
# # #   facet_wrap(~Pair.Name,scale="free_y")+geom_smooth(method="lm")#, formula=y~s(x,k=3))
# # # 
# # # ggplot(macrosubset[Pair.Name%in%slopes[Estimate<0][order(pval)]$Pair.Name[1:20]], aes(x=6-Disc_V2, y=log(1+Signal),col=dynamic_class3 ,group=dynamic_class3  ))+geom_point()+
# # #   facet_wrap(~Pair.Name,scale="free_y")+geom_smooth(method="lm")#, formula=y~s(x,k=3))
# # # 
# # # # not varying by resistant vs sensitive class :,col=dynamic_class3 ,group=dynamic_class3 
# # # ggplot(macrosubset[grepl("EGFR",Pair.Name)], aes(x=Disc_V1, y=log(1+Signal)))+geom_point()+geom_violin(aes(group=Disc_V1),scale="width")+
# # #   facet_wrap(~Pair.Name,scale="free_y")+geom_smooth(method="lm",se=F)+
# # #   labs(x="EGFR activated proliferation",y="Cancer to fibroblast \n growth signal via EGFR") #, formula=y~s(x,k=3))
# # # 
# # # 
# # # 
# # # ggplot(macrosubset[Day==0][Pair.Name%in%slopes[Estimate<0][order(pval)]$Pair.Name[1:20]], aes(x=Disc_V1, y=log(1+Signal) ))+geom_point()+
# # #   facet_wrap(~Pair.Name,scale="free_y")+geom_smooth(method="lm")#, formula=y~s(x,k=3))
# # # 
# # # #CCI[,scaleTransduction:=scale(TransductionMu,center=F),by=c("Pair.Name")] #"ReceptorPhenoCelltype",
# # # #CCI[!is.finite(TransductionMu),scaleTransduction:=0]
# # # #CCI[,scalelnTransductionMu:=scale(log(1+TransductionMu)),by=Pair.Name]
# # # #CCI[,c("L","R"):= tstrsplit(Pair.Name,"_",fixed=T)]
