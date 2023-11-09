rm(list=ls())
require(data.table)
require(dplyr)
require(ggplot2)
require(tidyr)
require(ggsci)

require(lme4)
require(lmerTest)

### Phenotype data
load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/PhenotypesOfAllCellTypesAllArmsCohort2/PhenotypesOfAllCellTypesAllArmsCohort2.RData")
# Encode the subtypes that have been analysed using UMAP
tmp <- data.table(allphenotypes %>% group_by(key_) %>% slice(1)) %>% dplyr::select(Celltype , Celltype_subtype) %>% unique
tmp[ , PhenoCelltype:= Celltype]
tmp[Celltype_subtype %in% c("CD4+ T cells", "Tregs"), PhenoCelltype:= "CD4+ T cells"]
tmp[Celltype_subtype %in% c("CD8+ T cells", "NK cells"), PhenoCelltype:= "CD8+ T cells"]
tmp[Celltype_subtype %in% c("Macrophages", "DC", "Monocytes"), PhenoCelltype:= "Macrophages"]
tmp[Celltype_subtype %in% c("Vas-Endo", "Lym-Endo","Endothelial cells"), PhenoCelltype:= "Endothelial cells"]
tmp[Celltype_subtype %in% c("Cancer cells"), PhenoCelltype:= "Cancer cells"]
tmp[Celltype_subtype %in% c("Normal epithelial cells"), PhenoCelltype:= "Normal epithelial cells"]
allphenotypes <- merge(allphenotypes %>%dplyr::select(-c( paste0("V", 1:5), "nCount_RNA", "nFeature_RNA",  "Platform","Sample_p_t", "file_string", "day_fact")), tmp, by= c("Celltype", "Celltype_subtype") )
allphenotypes[Patient.Study.ID=="001-171_853686",Patient.Study.ID:="001-171"]
load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE/Feline clinical input.RData")
allphenotypes <- merge(allphenotypes, unique(Clin_resp_dd_classAdd%>%dplyr::select(Patient.Study.ID,prop_change,dynamic_class)), by="Patient.Study.ID")
allphenotypes[PhenoCelltype%in%c("Macrophages","CD8+ T cells","CD4+ T cells")]
allphenotypes[PhenoCelltype%in%c("Macrophages")]
allphenotypes[PhenoCelltype%in%c("CD8+ T cells","CD4+ T cells")]
allphenotypes[PhenoCelltype%in%c("CD8+ T cells")]

### Unique clusters of cells (ALL) and their umap discretization level
uu <- unique( allphenotypes %>% dplyr::select( c("Celltype","PhenoCelltype", "key_", paste0("Disc_V", 1:5) ) ) )  #"Celltype_subtype",

# Overlaying signal received on cancer umap
load(file= "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Tcell gene expression wrangling Cohort2/TcellsSubclassCCICohort2.RData")
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

subsetInflam <- TcellSubtypeCommCCI[][Pair.Name %in% LRpairsFiltered[
  HPMR.Receptor%in% cytokineReceptors
  ]$Pair.Name][LigandPhenoCelltype=="Macrophages"]

subsetInflam[,startval:=sum((Day==0)*log(1+scaleTransduction))/sum(Day==0), by=c("Pair.Name","Patient.Study.ID","ReceptorPhenoCelltype")]
subsetInflam[ ,startvalAv:=median(startval,na.rm=T), by=c("Pair.Name","dynamic_class3","ReceptorPhenoCelltype","ARM")]

subsetInflam[Patient.Study.ID=="001-171_853686",Patient.Study.ID:="001-171"]
pairstotest <- c("IL15_IL2RA","IL15_IL2RB","IL18_IL18R1","IL18_IL18RAP", "IL15_IL15RA" )
subsetInflam[][Pair.Name%in% pairstotest]$Patient.Study.ID%>%unique()%>%length()
subsetInflam[Pair.Name%in% pairstotest]%>%dplyr::select(Patient.Study.ID,Day)%>%unique()%>%nrow()

activeoutRiboValid <- rbindlist(lapply( 1:length(pairstotest),function(x){
  print (x)
  res<-tryCatch({
    m1 <- lmer( log(1+scaleTransduction) ~ log(1+Day) * dynamic_class3 +(1|Patient.Study.ID), data=   subsetInflam[Patient.Study.ID!="001-167"][ARM!="A"][Pair.Name== pairstotest[x] ])
    data.table( Pair.Name=pairstotest[x],  data.table(summary(m1)$coef,keep.rownames=T) )
  },
  error=function(x){
    return(NULL)
  })
  setnames(res,old=c("t value","Pr(>|t|)"),new=c("tval","pval"))
  res
}))[rn=="log(1 + Day):dynamic_class3Response"]

activeoutLetrozoleValid <- rbindlist(lapply( 1:length(pairstotest),function(x){
  print (x)
  res<-tryCatch({
    m1 <- lmer( log(1+scaleTransduction) ~ log(1+Day) * dynamic_class3 +(1|Patient.Study.ID), data=   subsetInflam[Patient.Study.ID!="001-167"][ARM=="A"][Pair.Name== pairstotest[x] ])
    data.table( Pair.Name=pairstotest[x],  data.table(summary(m1)$coef,keep.rownames=T) )
  },
  error=function(x){
    return(NULL)
  })
  setnames(res,old=c("t value","Pr(>|t|)"),new=c("tval","pval"))
  res
}))[rn=="log(1 + Day):dynamic_class3Response"]



pairstotest <- c("ADAM12_ITGB1" , "COL18A1_ITGB1", "F13A1_ITGA4","F13A1_ITGB1","ICAM2_ITGB2"  , "ICAM3_ITGB2","VCAM1_ITGA4","VCAM1_ITGB1","THBS2_ITGA4")
recruitoutRiboValid <- rbindlist(lapply( 1:length(pairstotest),function(x){
  print (x)
  res<-tryCatch({
    m1 <- lmer( log(1+scaleTransduction) ~ log(1+Day) * dynamic_class3 +(1|Patient.Study.ID), data=   subsetInflam[Patient.Study.ID!="001-167"][ARM!="A"][Pair.Name== pairstotest[x] ])
    data.table( Pair.Name=pairstotest[x],  data.table(summary(m1)$coef,keep.rownames=T) )
  },
  error=function(x){
    return(NULL)
  })
  setnames(res,old=c("t value","Pr(>|t|)"),new=c("tval","pval"))
  res
}))[rn=="log(1 + Day):dynamic_class3Response"]

recruitoutLetrozoleValid<- rbindlist(lapply( 1:length(pairstotest),function(x){
  print (x)
  res<-tryCatch({
    m1 <- lmer( log(1+scaleTransduction) ~ log(1+Day) * dynamic_class3 +(1|Patient.Study.ID), data=   subsetInflam[Patient.Study.ID!="001-167"][ARM=="A"][Pair.Name== pairstotest[x] ])
    data.table( Pair.Name=pairstotest[x],  data.table(summary(m1)$coef,keep.rownames=T) )
  },
  error=function(x){
    return(NULL)
  })
  setnames(res,old=c("t value","Pr(>|t|)"),new=c("tval","pval"))
  res
}))[rn=="log(1 + Day):dynamic_class3Response"]

recruitoutRiboValid[pval<0.05]
recruitoutLetrozoleValid[pval<0.05]
recruitoutLetrozoleValid[pval>0.05]
activeoutRiboValid[pval<0.05]
activeoutLetrozoleValid[pval<0.05]

totalstatsValidation <- rbind( data.table(Cohort="Validation",Treatment="Combination ribociclib" ,Function="Recruitment" , recruitoutRiboValid),
                              data.table(Cohort="Validation",Treatment="Letrozole alone" ,Function="Recruitment" , recruitoutLetrozoleValid),
                              data.table(Cohort="Validation",Treatment="Combination ribociclib" ,Function="Activation" , activeoutRiboValid),
                              data.table(Cohort="Validation",Treatment="Letrozole alone" ,Function="Activation" , activeoutLetrozoleValid) )
write.csv(totalstatsValidation,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/temporary/ValidationStatsMyeloidTcellComm.csv")

rm(list=ls())
totalstatsValidation<- data.table(read.csv(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/temporary/ValidationStatsMyeloidTcellComm.csv"))
totalstatsDiscovery<- data.table( read.csv(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/temporary/DiscoveryStatsMyeloidTcellComm.csv"))
allstats<-rbind(totalstatsDiscovery,totalstatsValidation)
allstats[,X:=NULL]
write.csv(allstats,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/Figure5/DiscoveryandValidationStatsMyeloidTcellComm.csv")
