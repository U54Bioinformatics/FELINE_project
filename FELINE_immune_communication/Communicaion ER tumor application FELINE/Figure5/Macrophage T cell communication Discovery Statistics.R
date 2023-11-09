rm(list=ls())
require(data.table)
require(dplyr)
require(ggplot2)
require(tidyr)
require(ggsci)
require(viridis)
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
allphenotypes[PhenoCelltype%in%c("Macrophages","CD8+ T cells","CD4+ T cells")]
allphenotypes[PhenoCelltype%in%c("Macrophages")]
allphenotypes[PhenoCelltype%in%c("CD8+ T cells","CD4+ T cells")]
allphenotypes[PhenoCelltype%in%c("CD8+ T cells")]

### Unique clusters of cells (ALL) and their umap discretization level
uu <- unique( allphenotypes %>% dplyr::select( c("Celltype","PhenoCelltype", "key_", paste0("Disc_V", 1:5) ) ) )  #"Celltype_subtype",

### Load counts per million read data for a specific patient's samples and gen Ligand+Receptor genes
CPMlocs <- "/Users/jason/Dropbox/FELINE Project (1)/Data_analysis/FELINE_data_folder/scRNA_count_CPM/output/"
CPMfiles <- list.files(CPMlocs)[ grep("CPM", list.files(CPMlocs) ) ]
n10Xpats <- length(CPMfiles)


# settings
shouldscale <- FALSE ## should cpm be scaled


# Overlaying signal received on cancer umap
load(file= "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Tcell gene expression wrangling/TcellsSubclassCCI.RData")
#CCIperReceiverclass 
TcellSubtypeCommCCI <- CCIperReceiverclass%>%dplyr::select(Patient.Study.ID,Day,dynamic_class3,ARM,Pair.Name, key_,LigandPhenoCelltype,ReceptorPhenoCelltype,SumSignal)
rm(list="CCIperReceiverclass")

### Load Ligand Receptor database list of Ramilowski et al 2015
load( "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/LigandReceptor/Filtered_Human-2015-Ramilowski-LR-pairs.RData")
LRgenelist <- unique( c(LRpairsFiltered$HPMR.Receptor, LRpairsFiltered$HPMR.Ligand) )



load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/GeneOntology/cytokineReceptors.RData" )
#cytokineReceptors,dd1reduced
TcellSubtypeCommCCI[, scaleTransduction:= exp(scale(log(SumSignal), center=F)), by=c("Pair.Name")] 
TcellsummaryInteractionGF <- data.table(TcellSubtypeCommCCI[Pair.Name %in% LRpairsFiltered[
  HPMR.Receptor%in% cytokineReceptors
  ]$Pair.Name]%>%
    group_by(Patient.Study.ID,LigandPhenoCelltype,ReceptorPhenoCelltype,Day,dynamic_class3,ARM,key_) %>%
    dplyr::summarise(muln_scaleTransduction=mean(log(1+scaleTransduction))))
#rm(list="TcellSubtypeCommCCI")

TcellSubtypeCommCCI[,Treatmentlab:= "Combination ribociclib"]
TcellSubtypeCommCCI[ARM=="A",Treatmentlab:= "Letrozole alone"]


subsetInflam <- TcellSubtypeCommCCI[Pair.Name %in% LRpairsFiltered[
  HPMR.Receptor%in% cytokineReceptors
  ]$Pair.Name][LigandPhenoCelltype=="Macrophages"]
subsetInflam[,startval:=sum((Day==0)*log(1+scaleTransduction))/sum(Day==0), by=c("Treatmentlab","Pair.Name","Patient.Study.ID","ReceptorPhenoCelltype")]
subsetInflam[ ,startvalAv:=median(startval,na.rm=T), by=c("Treatmentlab","Pair.Name","dynamic_class3","ReceptorPhenoCelltype")]

pairstotest <- c("IL15_IL2RA","IL15_IL2RB","IL18_IL18R1","IL18_IL18RAP", "IL15_IL15RA" )

activeoutRiboDisc <- rbindlist(lapply( 1:length(pairstotest),function(x){
  print (x)
  res<-tryCatch({
    m1 <- lmer( log(1+scaleTransduction) ~ log(1+Day) * dynamic_class3 +(1|Patient.Study.ID), data=   subsetInflam[ARM!="A"][Pair.Name== pairstotest[x] ])
    data.table( Pair.Name=pairstotest[x],  data.table(summary(m1)$coef,keep.rownames=T) )
  },
  error=function(x){
    return(NULL)
  })
  setnames(res,old=c("t value","Pr(>|t|)"),new=c("tval","pval"))
  res
}))[rn=="log(1 + Day):dynamic_class3Response"]

activeoutLetrozoleDisc <- rbindlist(lapply( 1:length(pairstotest),function(x){
  print (x)
  res<-tryCatch({
    m1 <- lmer( log(1+scaleTransduction) ~ log(1+Day) * dynamic_class3 +(1|Patient.Study.ID), data=   subsetInflam[ARM=="A"][Pair.Name== pairstotest[x] ])
    data.table( Pair.Name=pairstotest[x],  data.table(summary(m1)$coef,keep.rownames=T) )
  },
  error=function(x){
    return(NULL)
  })
  setnames(res,old=c("t value","Pr(>|t|)"),new=c("tval","pval"))
  res
}))[rn=="log(1 + Day):dynamic_class3Response"]



pairstotest <- c("ADAM12_ITGB1" , "COL18A1_ITGB1", "F13A1_ITGA4","F13A1_ITGB1","ICAM2_ITGB2"  , "ICAM3_ITGB2","VCAM1_ITGA4","VCAM1_ITGB1","THBS2_ITGA4")
recruitoutRiboDisc <- rbindlist(lapply( 1:length(pairstotest),function(x){
  print (x)
  res<-tryCatch({
    m1 <- lmer( log(1+scaleTransduction) ~ log(1+Day) * dynamic_class3 +(1|Patient.Study.ID), data=   subsetInflam[ARM!="A"][Pair.Name== pairstotest[x] ])
    data.table( Pair.Name=pairstotest[x],  data.table(summary(m1)$coef,keep.rownames=T) )
  },
  error=function(x){
    return(NULL)
  })
  setnames(res,old=c("t value","Pr(>|t|)"),new=c("tval","pval"))
  res
}))[rn=="log(1 + Day):dynamic_class3Response"]

recruitoutLetrozoleDisc <- rbindlist(lapply( 1:length(pairstotest),function(x){
  print (x)
  res<-tryCatch({
    m1 <- lmer( log(1+scaleTransduction) ~ log(1+Day) * dynamic_class3 +(1|Patient.Study.ID), data=   subsetInflam[ARM=="A"][Pair.Name== pairstotest[x] ])
    data.table( Pair.Name=pairstotest[x],  data.table(summary(m1)$coef,keep.rownames=T) )
  },
  error=function(x){
    return(NULL)
  })
  setnames(res,old=c("t value","Pr(>|t|)"),new=c("tval","pval"))
  res
}))[rn=="log(1 + Day):dynamic_class3Response"]

recruitoutRiboDisc[pval<0.05]
recruitoutLetrozoleDisc[pval<0.05]
activeoutRiboDisc[pval<0.05]
activeoutLetrozoleDisc[pval<0.05]

totalstatsDiscovery <- rbind( data.table(Cohort="Discovery",Treatment="Combination ribociclib" ,Function="Recruitment" , recruitoutRiboDisc),
       data.table(Cohort="Discovery",Treatment="Letrozole alone" ,Function="Recruitment" , recruitoutLetrozoleDisc),
       data.table(Cohort="Discovery",Treatment="Combination ribociclib" ,Function="Activation" , activeoutRiboDisc),
       data.table(Cohort="Discovery",Treatment="Letrozole alone" ,Function="Activation" , activeoutLetrozoleDisc) )
write.csv(totalstatsDiscovery,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/temporary/DiscoveryStatsMyeloidTcellComm.csv")
