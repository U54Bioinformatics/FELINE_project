rm(list=ls())   
require(rTensor); require(nnTensor); library(abind)
require(data.table); require(dplyr); require(ggplot2); require(tidyr)
require(igraph)
require(ggsci)
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

### Unique clusters of cells (ALL) and their umap discretization level
uu <- unique( allphenotypes %>% dplyr::select( c("Celltype","PhenoCelltype", "key_", paste0("Disc_V", 1:5) ) ) )  #"Celltype_subtype",

### Load Ligand Receptor database list of Ramilowski et al 2015
load( "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/LigandReceptor/Filtered_Human-2015-Ramilowski-LR-pairs.RData")
LRgenelist <- unique( c(LRpairsFiltered$HPMR.Receptor, LRpairsFiltered$HPMR.Ligand) )

### Load some signalling data for downstream analysis
savelocCCI <- "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/CommunicationOutputAllArmsCohort2/"

# Are we studying the per indiv or per cell type level of signalling?
perIndiv=FALSE
load(file ="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/cohort2 Communication output merged/PopulationCommunicationMerged ALL ARMS cohort2.RData" )#CCI,allphenotypes, uu,perIndiv,
CCI[,Treat:="CombinationRibo"]
CCI[ARM=="A",Treat:="LetrozoleAlone"]

CCI[,scaleTransduction:=scale(TransductionMu,center=F),by=c("Pair.Name")] 
CCI[!is.finite(TransductionMu),scaleTransduction:=0]
CCI[Patient.Study.ID=="001-171_853686",Patient.Study.ID:="001-171"]
CCI <-merge(CCI, unique(Clin_resp_dd_classAdd%>%dplyr::select(Patient.Study.ID,prop_change,dynamic_class)), by="Patient.Study.ID")
CCI[ prop_change < (2/3) ,dynamic_class3:="Response"]
CCI[,scalelnTransductionMu:=scale(log(1+TransductionMu)),by=Pair.Name]

CCIvalid <- CCI
rm(list=ls()[ls()!="CCIvalid"])



load(file ="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Communication output merged ALL ARMS/PopulationCommunicationMerged ALL ARMS.RData" )#CCI,allphenotypes, uu,perIndiv,
CCI[,Treat:="CombinationRibo"]
CCI[ARM=="A",Treat:="LetrozoleAlone"]
CCI[,scaleTransduction:=scale(TransductionMu,center=F),by=c("Pair.Name")] #"ReceptorPhenoCelltype",
CCI[!is.finite(TransductionMu),scaleTransduction:=0]
CCI[,scalelnTransductionMu:=scale(log(1+TransductionMu)),by=Pair.Name]

CCIdisc<-CCI

CCI <- rbind(data.table(Cohort="Discovery", CCIdisc%>%dplyr::select(Patient.Study.ID,LigandPhenoCelltype,ReceptorPhenoCelltype,Day,Pair.Name,dynamic_class3,ARM,Treat,scalelnTransductionMu)),
             data.table(Cohort="Validation", CCIvalid%>%dplyr::select(Patient.Study.ID,LigandPhenoCelltype,ReceptorPhenoCelltype,Day,Pair.Name,dynamic_class3,ARM,Treat,scalelnTransductionMu) ))

commsvect<-c("CLCF1_CRLF1" ,  "CLCF1_IL6ST"  , "CLCF1_LIFR" ,   "CSF1_CSF1R" ,   "FGF1_FGFR1" ,   "FGF17_FGFR1"  , "FGF18_FGFR1" ,  "FGF2_CD44"  ,   "FGF2_FGFR1" ,   "FGF2_NRP1" ,   
  "FGF2_SDC1" ,    "FGF2_SDC2"  ,   "FGF2_SDC3" ,    "FGF2_SDC4" ,    "FGF23_FGFR1" ,  "FGF7_FGFR1"  ,  "GDF9_ACVR2A"  , "GDF9_BMPR1A"  , "GDF9_BMPR1B"  , "GDF9_BMPR2" ,  
  "GDF9_TGFBR1" ,  "IL11_IL11RA"  , "IL11_IL6ST" ,   "IL12A_IL12RB1" , "IL12A_IL12RB2" , "IL5_CSF2RB",    "MADCAM1_CD44" , "MADCAM1_ITGA4", "PLAU_IGF2R",    "PLAU_LRP1"  ,  
  "PLAU_LRP2"   ,  "PLAU_PLAUR")

sampcodesdiscRibo<-c("001-1320"  ,    "001-1130","001-1020","2972-006-6010", "001-1190" ,   "001-1400","001-1430","001-1450","001-1180","001-1250"     ,
  "001-1420","2972-007-7030")
sampcodesdiscLetrozole<-c("001-1050","001-1150","001-1160","001-1370","2972-005-2020" ,"001-1030","001-1070","001-1210","001-1220","001-1310"  ,   
  "001-1340" )
sampcodesall<- c(sampcodesdiscRibo,sampcodesdiscLetrozole)

sampcodesvalidLetrozole<-c("001-1610","001-1650","001-1710","2972-002-3070", "2972-002-3110" ,"2972-009-9020", "2972-009-9040" ,"001-1590","001-1690","2972-002-3040",
   "2972-009-9070")

sampcodesvalidRibo<-c("001-1080","001-1270","001-1670","001-1720","2972-007-7090", "001-1510","001-1520","001-1530","001-1600","001-1640"     ,
 "001-1660","001-1680","2972-003-4080", "2972-007-7060", "2972-009-9010" ,"2972-009-9060")
sampcodesall2<- c(sampcodesvalidLetrozole,sampcodesvalidRibo)

CancertoMacrophageM2comm <- CCI[][Day!=180][LigandPhenoCelltype=="Cancer cells"][ReceptorPhenoCelltype=="Macrophages"][Pair.Name%in%commsvect][paste0(Patient.Study.ID,Day)%in%c(sampcodesall,sampcodesall2)]
#save(CancertoMacrophageM2comm,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/Figure4/Discovery and Validation CancertoMyeloid M2macrophageCommunications.RData")

CancertoMacrophageM2commF<- CCI[][Day!=180][LigandPhenoCelltype=="Cancer cells"][ReceptorPhenoCelltype=="Macrophages"][Pair.Name%in%commsvect][]

statsF <- rbindlist( lapply( commsvect , function(i){
  data.table(Pair.Name=i, coef(summary(lm(scalelnTransductionMu~dynamic_class3 ,CancertoMacrophageM2commF[Pair.Name==i][Treat=="CombinationRibo"]))), keep.rownames = T)
}) )
setnames(statsF, old=c("Std. Error","t value","Pr(>|t|)"), new=c("Std.Error","tvalue","pvalue"))
statsF[rn!="(Intercept)"][rn=="dynamic_class3Response"][pvalue<0.05]
RiboRes<-statsF[rn!="(Intercept)"][rn=="dynamic_class3Response"]
#write.csv(RiboRes,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/temporary/CancerMyeloidM2communication Combination ribociclib.csv")

statsFLetrozole <- rbindlist( lapply( commsvect , function(i){
  data.table(Pair.Name=i, coef(summary(lm(scalelnTransductionMu~dynamic_class3 ,CancertoMacrophageM2commF[Pair.Name==i][Treat!="CombinationRibo"]))), keep.rownames = T)
}) )
setnames(statsFLetrozole, old=c("Std. Error","t value","Pr(>|t|)"), new=c("Std.Error","tvalue","pvalue"))
statsFLetrozole[rn!="(Intercept)"][rn=="dynamic_class3Response"][pvalue<0.05]
LetrozoleRes<-statsFLetrozole[rn!="(Intercept)"][rn=="dynamic_class3Response"]
#write.csv(LetrozoleRes,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/temporary/CancerMyeloidM2communication Letrozole alone.csv")




