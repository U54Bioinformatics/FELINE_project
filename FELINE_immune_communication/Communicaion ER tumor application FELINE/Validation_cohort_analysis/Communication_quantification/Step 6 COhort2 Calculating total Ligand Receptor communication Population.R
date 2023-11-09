rm(list=ls())   
##devtools::install_github("rikenbit/nnTensor") #install.packages("rTensor")
#install.packages("https://cran.r-project.org/src/contrib/Archive/rTensor/rTensor_1.4.tar.gz", repo=NULL, type="source")
require(rTensor); require(nnTensor); library(abind)
require(data.table); require(dplyr); require(ggplot2); require(tidyr)
require(igraph);require(ggsci)
### Load clinical data
#load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE/Feline clinical input.RData") #save(Clin_resp_dd,Clin_resp_dd_class,Clin_resp_dd_classAdd,patientCode_LU,res_file_46,response_code_dd,response_dd,file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE/Feline clinical input.RData") #gsea_path<-"~/Dropbox/FELINE Project/Data_analysis/scRNA/05_ssGSEA_score/Signature_c2_hallmark/results/"   #"~/Dropbox/FELINE/Data_share/Modeling_Data/pathway_seperated_files/Data_gene_count_per_celltype_model_zinbwave_ssGSEA/FEL001043/"

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

### Unique clusters of cells (ALL) and their umap discretization level
uu <- unique( allphenotypes %>% dplyr::select( c("Celltype","PhenoCelltype", "key_", paste0("Disc_V", 1:5) ) ) )  #"Celltype_subtype",

### Load Ligand Receptor database list of Ramilowski et al 2015
load( "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/LigandReceptor/Filtered_Human-2015-Ramilowski-LR-pairs.RData")
LRgenelist <- unique( c(LRpairsFiltered$HPMR.Receptor, LRpairsFiltered$HPMR.Ligand) )

### Load some signalling data for downstream analysis
savelocCCI <- "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/CommunicationOutputAllArmsCohort2/"

# Are we studying the per indiv or per cell type level of signalling?
perIndiv=FALSE

# Names of LR data to analyze
filenamesCCI <- list.files(savelocCCI)#%>%length

#summarytable <- read.csv("/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/TensorAnalysisOutput/Ligand Receptor list fro NTD model of population cellcell communication_AB.csv")


##If not done, merge data across patients
# extract total cell-cell signal transduction (communication interaction) from any cell to a focal cell type. Integrate over the detailed communication data.table
TotCCI <- rbindlist(lapply(1:length(filenamesCCI), function(ii){
  load(file= paste0(savelocCCI, filenamesCCI[ii]))
  cat(ii);
  # for a given sample, select signal data for each cell type and composition and response data
  SumComm <- unique(data.table(Communication %>% dplyr::select(key_,Patient.Study.ID,ReceptorPhenoCelltype,    Day ,Pair.Name,dynamic_class3,ARM, FracSample,Transduction,Receptor) ))
  return(SumComm)#[!( is.na(ReceptorPhenoCelltype) ) ]
}))
# merge phenotype data for each cell type
TotCCImrg <- merge(TotCCI,uu,by=c("key_")) #[!( is.na(ReceptorPhenoCelltype) ) ]
#TotCCImrg[,ReceptorPhenoCelltype:=PhenoCelltype]

# summarise amount of signal received by the average cell of each cell type , weighting for each cell subtype's abundance
SumTotComm0 <- data.table(TotCCImrg%>%group_by(Patient.Study.ID ,Day, PhenoCelltype, Pair.Name, dynamic_class3,ARM ) %>%dplyr::summarise(
  Transduction= weighted.mean(Transduction,FracSample),
  FracSample= sum(FracSample)
))
#save(TotCCImrg, SumTotComm0,file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/cohort2 TotalCommunicationAllArms/cohort2 TotalCommunicationAllArms.RData")

load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/cohort2 TotalCommunicationAllArms/cohort2 TotalCommunicationAllArms.RData")


#######Continute from here 





# scale each communication pathway 
SumTotComm0[,scaleTransduction:=scale(Transduction,center=F),by=c("PhenoCelltype","Pair.Name")]#
SumTotComm0[!is.finite(scaleTransduction),scaleTransduction:=0]
SumTotComm0[,Treat:="CombinationRibo"]
SumTotComm0[ARM=="A",Treat:="LetrozoleAlone"]

# look up of the LR pairs and signal-receivers pairs to consider 
lu_lm <- unique(SumTotComm0%>%dplyr::select(Pair.Name,PhenoCelltype,Treat))


## now for the more intricate question time*response
assessment <- rbindlist(lapply(1:nrow(lu_lm), function(i){
  cat(i) ; cat("  ")
  yi <-  SumTotComm0[Pair.Name==lu_lm[i]$Pair.Name][PhenoCelltype==lu_lm[i]$PhenoCelltype][Treat==lu_lm[i]$Treat]
  if(length(unique(yi[Day==max(Day)]$dynamic_class3))>1  & length(unique(yi[Day!=max(Day)]$dynamic_class3))>1  &
     length(unique(yi[dynamic_class3=="Non-response"]$Day))>1 & length(unique(yi[dynamic_class3=="Response"]$Day))>1 &
     min( (yi%>%group_by(dynamic_class3)%>%dplyr::summarise(np=length(unique(Patient.Study.ID)))) $np ) >2
  ){
    m1 <- lm(log(1+scaleTransduction) ~ log(1+Day)*dynamic_class3 , data=yi)
    m2 <- lm(log(1+scaleTransduction) ~ log(1+Day) , data=yi[dynamic_class3=="Response"])      #ggplot(yi,aes(y=log(1+TransductionMu) , x= log(1+Day) , col=dynamic_class3 ))+geom_point()
    m3 <- lm(log(1+scaleTransduction) ~ dynamic_class3 , data=yi[Day==180])      #ggplot(yi,aes(y=log(1+TransductionMu) , x= log(1+Day) , col=dynamic_class3 ))+geom_point()
    m4 <- lm(log(1+scaleTransduction) ~ dynamic_class3 , data=yi[Day==0])        #ggplot(yi,aes(y=log(1+TransductionMu) , x= log(1+Day) , col=dynamic_class3 ))+geom_point()
    
    res <- data.table(lu_lm[i], as.data.table(summary(m1)$ coefficients,keep.rownames = T) , adj.r.squared=summary(m1)$ adj.r.squared,i)
    res[,Model:="two way anova"]
    res$paramNames <- c( "Non-response_int" ,"Non-response_slope", "NtoR_change_int", "NtoR_change_slope")
    res2 <- data.table(lu_lm[i], as.data.table(summary(m2)$ coefficients,keep.rownames = T) , adj.r.squared=summary(m2)$ adj.r.squared,i)
    res2[,Model:="one way regression_Response"]
    res2$paramNames <- c( "Response_int" ,"Response_slope")
    
    res3 <- data.table(lu_lm[i], as.data.table(summary(m3)$ coefficients,keep.rownames = T) , adj.r.squared=summary(m3)$ adj.r.squared,i)
    res3[,Model:="one way anova_DayEnd_Response"]
    res3$paramNames <- c( "Non-response_int_final" ,"NtoR_change_int_final")

    res4 <- data.table(lu_lm[i], as.data.table(summary(m4)$ coefficients,keep.rownames = T) , adj.r.squared=summary(m4)$ adj.r.squared,i)
    res4[,Model:="one way anova_DayStart_Response"]
    res4$paramNames <- c( "Non-response_int_start" ,"NtoR_change_int_start")
    
    resAll <- rbind(res,res2, res3, res4)
    return(resAll)
  }else{
    return(NULL)
  }
}))
setnames(assessment, old=c("rn","Pr(>|t|)"), new=c("coefficient","pval"))
# FDR correction
assessment <- data.table(assessment%>%group_by(Treat,Model)%>%dplyr::mutate(adjust.p.val= p.adjust(pval, method="fdr")))
#assessment[Treat=="CombinationRibo"][Pair.Name=="A2M_LRP1"][PhenoCelltype=="Adipocytes"]
# Subset signaificant terms
responseRelated <- assessment[adjust.p.val<0.05][ coefficient!="(Intercept)"]


# Count minimum number of samples that are contribuing to statistical effect at different timepoints 
repdd <- SumTotComm0%>%group_by(Day,Pair.Name,PhenoCelltype,dynamic_class3,Treat)%>%dplyr::summarise(replic=length(Transduction))%>%group_by(Pair.Name,PhenoCelltype,dynamic_class3,Treat)%>%
   dplyr::summarise(replic=min(replic))
# Get the significant effects that are identified when there are more than three samples at each timepoint
responseRelated2_0 <- merge(responseRelated[order(adjust.p.val)] ,repdd, by=c("Pair.Name","PhenoCelltype","Treat"),all.x=T)[replic>3][order(adjust.p.val)]
responseRelated2 <- responseRelated2_0[grepl("NtoR",paramNames)|(grepl("Non-response",paramNames)&dynamic_class3=="Non-response") |(grepl("Response",paramNames)&dynamic_class3=="Response")]
responseRelated3 <- unique(data.table(responseRelated2[!coefficient=="(Intercept)" ]%>%select(-i) ))   #|paramNames=="Non-response_int_final"
responseRelated3[grepl("NtoR",paramNames),dynamic_class3:="NRcomparison"] 
responseRelated4 <- unique(responseRelated3)

# unique(responseRelated4$PhenoCelltype)
# data.table( responseRelated4[] %>% group_by(PhenoCelltype) %>% dplyr::summarise(n_signif=length(unique(Pair.Name))) ) [order(-n_signif)]
# 
print.this <- unique(data.table( responseRelated4[order(Treat,PhenoCelltype,-adjust.p.val)] %>% select("Treat","PhenoCelltype","Pair.Name","paramNames","dynamic_class3","Model",
                                                                                                  "Estimate", "Std. Error", "t value", "pval", "adj.r.squared", "adjust.p.val") ) )
print.this[PhenoCelltype=="Fibroblasts"]
print.this[PhenoCelltype=="Cancer cells"][1:5]
print.this[paramNames=="NtoR_change_int"]
print.this[paramNames=="NtoR_change_int_start"]
print.this[paramNames=="NtoR_change_int_final"]

#save(assessment,print.this,responseRelated,responseRelated4,file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Cohort2 When cells communicate_communication changes AllArms/Cohort2 Total signal differences and trends by response mods.RData")
#write.csv(print.this,file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Trends in communication AllArms/Trends in TOTAL communication by response AllArms.csv")

load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/When cells communicate_communication changes AllArms/Total signal differences and trends by response mods.RData")


# 
# #CC<-"Cancer cells"
Treatplot<-"CombinationRibo" #Treatplot<-"LetrozoleAlone"
for(CC in unique(responseRelated4$PhenoCelltype) ){
  plotdd <- SumTotComm0[Pair.Name%in% unique(responseRelated4[Treat%in%Treatplot][PhenoCelltype==CC]$Pair.Name) ] [PhenoCelltype==CC]
  plotdd$Pair.Name <- factor(plotdd$Pair.Name, levels=unique(responseRelated4[PhenoCelltype==CC]$Pair.Name))
  if(nrow(plotdd[Treat==Treatplot])>0){
    ggplot(plotdd[Treat==Treatplot],#[grepl("IL",Pair.Name)],
           aes(y=log(1+scaleTransduction),x=log(1+Day),col=dynamic_class3))+geom_point()+geom_smooth(method="lm")+
      facet_wrap(~Pair.Name,scales="free") +theme_classic()+theme(aspect.ratio=1,axis.text.x = element_text(angle = 90))+
      labs(y=paste0(CC, " total signal received \n (standardized)"),x="Day" )+scale_color_npg(name="Tumor response")+scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))
    
      ggsave(filename=paste0("/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/LR TOTAL trends AllArms/RiboCombination Increasing and decreasing total communication",CC,".png"),width=10,height=10)
  }
}

Treatplot<-"LetrozoleAlone"
for(CC in unique(responseRelated4$PhenoCelltype) ){
  plotdd <- SumTotComm0[Pair.Name%in% unique(responseRelated4[Treat%in%Treatplot][PhenoCelltype==CC]$Pair.Name) ] [PhenoCelltype==CC]
  plotdd$Pair.Name <- factor(plotdd$Pair.Name, levels=unique(responseRelated4[PhenoCelltype==CC]$Pair.Name))
  if(nrow(plotdd[Treat==Treatplot])>0){
    ggplot(plotdd[Treat==Treatplot],#[grepl("IL",Pair.Name)],
           aes(y=log(1+scaleTransduction),x=log(1+Day),col=dynamic_class3))+geom_point()+geom_smooth(method="lm")+
      facet_wrap(~Pair.Name) +theme_classic()+theme(aspect.ratio=1,axis.text.x = element_text(angle = 90))+
      labs(y=paste0(CC, " total signal received \n (standardized)"),x="Day" )+scale_color_npg(name="Tumor response")+scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))
    
    ggsave(filename=paste0("/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/LR TOTAL trends AllArms/LetrozoleAlone Increasing and decreasing total communication",CC,".png"),width=8,height=8)
  }
}

print.this[paramNames=="NtoR_change_int_start"]

responseRelated4[Model=="one way anova_DayStart_Response"][paramNames=="NtoR_change_int_start"][Treat=="CombinationRibo"]
responseRelated4[Model=="two way anova"][paramNames=="NtoR_change_int"][Treat=="CombinationRibo"]

InitDiffsRibo <- na.omit( assessment[(Model=="one way anova_DayStart_Response"&paramNames=="NtoR_change_int_start")|
           (Model=="two way anova"&paramNames=="NtoR_change_int")][Treat=="CombinationRibo"][adjust.p.val<0.05][order(pval)][1:20])
InitDiffsLetr <- na.omit( assessment[(Model=="one way anova_DayStart_Response"&paramNames=="NtoR_change_int_start")|
                                     (Model=="two way anova"&paramNames=="NtoR_change_int")][Treat!="CombinationRibo"][adjust.p.val<0.05][order(pval)][1:20])
#assessment[Model=="two way anova"&paramNames=="NtoR_change_int"][Treat="CombinationRibo"][adjust.p.val<0.05][order(pval)][1:20]

# CC <- "Cancer cells"
# plotdd <- SumTotComm0[Pair.Name%in% unique(InitDiffsRibo[PhenoCelltype=="Macrophages"]$Pair.Name) ] [PhenoCelltype==CC]
# plotdd$Pair.Name <- factor(plotdd$Pair.Name, levels=unique(InitDiffsRibo[PhenoCelltype==CC]$Pair.Name))
# if(nrow(plotdd)>0){
#   ggplot(plotdd,#[grepl("IL",Pair.Name)],
#          aes(y=log(1+scaleTransduction),x=log(1+Day),col=dynamic_class3))+geom_point()+geom_smooth(method="lm")+
#     facet_wrap(~Pair.Name) +theme_classic()+theme(aspect.ratio=1,axis.text.x = element_text(angle = 90))+
#     labs(y=paste0(CC, " total signal received \n (standardized)"),x="Day" )+scale_color_npg(name="Tumor response")+scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))
#   
  
  
  
for(CC in unique(InitDiffsRibo$PhenoCelltype) ){
  plotdd <- SumTotComm0[Pair.Name%in% unique(InitDiffsRibo[PhenoCelltype==CC]$Pair.Name) ] [PhenoCelltype==CC]
  plotdd$Pair.Name <- factor(plotdd$Pair.Name, levels=unique(InitDiffsRibo[PhenoCelltype==CC]$Pair.Name))
  if(nrow(plotdd)>0){
    ggplot(plotdd,#[grepl("IL",Pair.Name)],
           aes(y=log(1+scaleTransduction),x=log(1+Day),col=dynamic_class3))+geom_point()+geom_smooth(method="lm")+
      facet_wrap(~Pair.Name) +theme_classic()+theme(aspect.ratio=1,axis.text.x = element_text(angle = 90))+
      labs(y=paste0(CC, " total signal received \n (standardized)"),x="Day" )+scale_color_npg(name="Tumor response")+scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))
    
    ggsave(filename=paste0("/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/LR TOTAL initial condition diffs AllArms/RibociclibCombination Initial differences in total communication",CC,".png"),width=8,height=8)
  }
}

for(CC in unique(InitDiffsLetr$PhenoCelltype) ){
  plotdd <- SumTotComm0[Pair.Name%in% unique(InitDiffsLetr[PhenoCelltype==CC]$Pair.Name) ] [PhenoCelltype==CC]
  plotdd$Pair.Name <- factor(plotdd$Pair.Name, levels=unique(InitDiffsLetr[PhenoCelltype==CC]$Pair.Name))
  if(nrow(plotdd)>0){
    ggplot(plotdd,#[grepl("IL",Pair.Name)],
           aes(y=log(1+scaleTransduction),x=log(1+Day),col=dynamic_class3))+geom_point()+geom_smooth(method="lm")+
      facet_wrap(~Pair.Name) +theme_classic()+theme(aspect.ratio=1,axis.text.x = element_text(angle = 90))+
      labs(y=paste0(CC, " total signal received \n (standardized)"),x="Day" )+scale_color_npg(name="Tumor response")+scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))
    
    ggsave(filename=paste0("/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/LR TOTAL initial condition diffs AllArms/Letrozolealone Initial differences in total communication",CC,".png"),width=8,height=8)
  }
}


plotdd <- merge(SumTotComm0,InitDiffsRibo%>%dplyr::select(Pair.Name,PhenoCelltype,Treat),by=c("Pair.Name","PhenoCelltype","Treat"))
plotdd$Pair.Name <- factor(plotdd$Pair.Name, levels=unique(InitDiffsRibo$Pair.Name))

ggplot(plotdd,#[grepl("IL",Pair.Name)],
         aes(y=log(1+scaleTransduction),x=log(1+Day),col=dynamic_class3))+geom_point()+geom_smooth(method="lm")+
    facet_grid(PhenoCelltype~Pair.Name) +theme_classic()+theme(aspect.ratio=1,axis.text.x = element_text(angle = 90))+
    labs(y="Total signal received \n (standardized)",x="Day" )+
  scale_color_npg(name="Tumor response")+
  scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))
  
ggsave(filename=paste0("/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/LR TOTAL initial condition diffs AllArms/RibociclibCombination Initial differences in total communication.png"),width=16,height=8)
  


plotdd <- merge(SumTotComm0,InitDiffsLetr%>%dplyr::select(Pair.Name,PhenoCelltype,Treat),by=c("Pair.Name","PhenoCelltype","Treat"))
plotdd$Pair.Name <- factor(plotdd$Pair.Name, levels=unique(InitDiffsLetr$Pair.Name))

ggplot(plotdd,#[grepl("IL",Pair.Name)],
       aes(y=log(1+scaleTransduction),x=log(1+Day),col=dynamic_class3))+geom_point()+geom_smooth(method="lm")+
  facet_grid(PhenoCelltype~Pair.Name) +theme_classic()+theme(aspect.ratio=1,axis.text.x = element_text(angle = 90))+
  labs(y="Total signal received \n (standardized)",x="Day" )+
  scale_color_npg(name="Tumor response")+
  scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))

ggsave(filename=paste0("/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/LR TOTAL initial condition diffs AllArms/Letrozolealone  Initial differences in total communication.png"),width=16,height=8)




