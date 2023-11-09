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

# Names of LR data to analyze
filenamesCCI <- list.files(savelocCCI)#%>%length

# First let's do this for the population level
#summarytable <- read.csv("/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/TensorAnalysisOutput/Ligand Receptor list fro NTD model of population cellcell communication_AB.csv")
load(file ="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/cohort2 Communication output merged/PopulationCommunicationMerged ALL ARMS cohort2.RData" )#CCI,allphenotypes, uu,perIndiv,
CCI[,Treat:="CombinationRibo"]
CCI[ARM=="A",Treat:="LetrozoleAlone"]

CCI[,scaleTransduction:=scale(TransductionMu,center=F),by=c("Pair.Name")] #"ReceptorPhenoCelltype",
CCI[!is.finite(TransductionMu),scaleTransduction:=0]
CCI[Patient.Study.ID=="001-171_853686",Patient.Study.ID:="001-171"]
CCI2<- CCI <-merge(CCI, unique(Clin_resp_dd_classAdd%>%dplyr::select(Patient.Study.ID,prop_change,dynamic_class)), by="Patient.Study.ID")
CCI[ prop_change < (2/3) ,dynamic_class3:="Response"]
CCI2[ prop_change < (2/3) ,dynamic_class3:="Response"]

# look up of the LR pairs and signal-receivers pairs to consider 
lu_lm <- unique(CCI%>%dplyr::select(Pair.Name,LigandPhenoCelltype,ReceptorPhenoCelltype,Treat))
CCI$Pair.Name%>%unique%>%length

plotdat <- CCI[Treat=="CombinationRibo"][LigandPhenoCelltype%in%c("Cancer cells") ][ReceptorPhenoCelltype%in%c("Macrophages")][Pair.Name%in%c("CCL19_CCR7","CCL5_CCR5","GNAI2_CCR5","TNFSF14_LTBR","MDK_SDC1","SORBS1_INSR")]
plotdat[,scalelnTransductionMu:=scale(log(1+TransductionMu)),by=Pair.Name]
plotdat[,startval:=sum((Day==0)*scalelnTransductionMu)/sum(Day==0), by=c("Pair.Name","Patient.Study.ID","ReceptorPhenoCelltype")]
plotdat[ ,startvalAv:=median(startval,na.rm=T), by=c("Pair.Name","dynamic_class3","ReceptorPhenoCelltype")]
require(ggsci)
ggplot(plotdat[Pair.Name!="MDK_SDC1"], 
       aes(y=scalelnTransductionMu,x=log(1+Day),col=dynamic_class3,group=dynamic_class3,fill=dynamic_class3,shape=Pair.Name))+ theme_classic()+
  geom_point(aes(shape=Pair.Name),size=4,position = position_dodge(width=1.9)) + geom_violin(aes(group=interaction(dynamic_class3,Day)),alpha=0.5, scale="width",position=position_dodge(width=1.9))+
  #facet_grid(LigandPhenoCelltype~ Pair.Name,scales="free")+
  labs(y= "Immune activating communication"   , x= "Day")+
  scale_x_continuous(breaks= log(1+c(0,14,180)), labels= c(0,14,180) )+
  theme(aspect.ratio= 1)+
  scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  scale_shape_manual(values=c(15:19,8),name="Communication \n pathway" ,labels= gsub("_","-",sort(unique(plotdat[Pair.Name!="MDK_SDC1"]$Pair.Name))))+
  geom_smooth(method="lm",se=F)
#geom_smooth(data=plotdat[Day<180],method="lm")
   #geom_smooth(method="gam", formula=y~s(x,k=3))+facet_wrap(~Pair.Name)
ggsave(file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Macrophage phenotypes AllArms/Cohort2 CancertoMacrophage immune activating signals over time by Response.png")

plotdat <- CCI[Treat=="CombinationRibo"][LigandPhenoCelltype%in%c("Cancer cells") ][ReceptorPhenoCelltype%in%c("Macrophages")][grep( "CCR",Pair.Name)]
plotdat[,scalelnTransductionMu:=scale(log(1+TransductionMu)),by=Pair.Name]
plotdat[,startval:=sum((Day==0)*scalelnTransductionMu)/sum(Day==0), by=c("Pair.Name","Patient.Study.ID","ReceptorPhenoCelltype")]
plotdat[ ,startvalAv:=median(startval,na.rm=T), by=c("Pair.Name","dynamic_class3","ReceptorPhenoCelltype")]

ggplot(plotdat, 
       aes(y=scalelnTransductionMu+startvalAv-startval,x=log(1+Day),col=dynamic_class3,group=dynamic_class3,fill=dynamic_class3))+ theme_classic()+
  geom_point(size=4,position = position_dodge(width=1.9)) + geom_violin(aes(group=interaction(dynamic_class3,Day)),alpha=0.5, scale="width",position=position_dodge(width=1.9))+
  #facet_grid(LigandPhenoCelltype~ Pair.Name,scales="free")+
  labs(y= "Immune activating communication"   , x= "Day")+
  scale_x_continuous(breaks= log(1+c(0,14,180)), labels= c(0,14,180) )+
  theme(aspect.ratio= 1)+
  scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  geom_smooth(method="lm",se=F)


## now for the more intricate question time*response
assessment <- rbindlist(lapply(1:nrow(lu_lm), function(i){
  cat(i) ; cat("  ")
  yi <-  CCI[Pair.Name==lu_lm[i]$Pair.Name][LigandPhenoCelltype==lu_lm[i]$LigandPhenoCelltype][ReceptorPhenoCelltype==lu_lm[i]$ReceptorPhenoCelltype][Treat==lu_lm[i]$Treat]
  if(length(unique(yi[Day==max(Day)]$dynamic_class3))>1  & length(unique(yi[Day!=max(Day)]$dynamic_class3))>1  &
     length(unique(yi[dynamic_class3=="Non-response"]$Day))>1 & length(unique(yi[dynamic_class3=="Response"]$Day))>1 &
     min( (yi%>%group_by(dynamic_class3)%>%dplyr::summarise(np=length(unique(Patient.Study.ID)))) $np ) >2 &
     length(unique(yi[Day==0]$dynamic_class3))>1 &length(unique(yi[Day==180]$dynamic_class3))>1 
     #& length(unique(yi[Day==14]$dynamic_class3))>1
  ){
     m1 <- lm(log(1+TransductionMu) ~ log(1+Day)*dynamic_class3 , data=yi)
     m2 <- lm(log(1+TransductionMu) ~ log(1+Day) , data=yi[dynamic_class3=="Response"])
     m3 <- lm(log(1+TransductionMu) ~ dynamic_class3 , data=yi[Day==180])      #ggplot(yi,aes(y=log(1+TransductionMu) , x= log(1+Day) , col=dynamic_class3 ))+geom_point()
     m4 <- lm(log(1+TransductionMu) ~ dynamic_class3 , data=yi[Day==0])        #ggplot(yi,aes(y=log(1+TransductionMu) , x= log(1+Day) , col=dynamic_class3 ))+geom_point()
     #m5 <- lm(log(1+TransductionMu) ~ dynamic_class3 , data=yi[Day==14])        #ggplot(yi,aes(y=log(1+TransductionMu) , x= log(1+Day) , col=dynamic_class3 ))+geom_point()
     
   
    #ggplot(yi,aes(y=log(1+TransductionMu) , x= log(1+Day) , col=dynamic_class3 ))+geom_point()+geom_smooth(method="lm")
    #spread(yi%>%select(Patient.Study.ID,TransductionMu,dynamic_class3,Day),Day,TransductionMu)
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
    
    #res5 <- data.table(lu_lm[i], as.data.table(summary(m5)$ coefficients,keep.rownames = T) , adj.r.squared=summary(m5)$ adj.r.squared,i)
    #res5[,Model:="one way anova_DayMiddle_Response"]
    #res5$paramNames <- c( "Non-response_int_middle" ,"NtoR_change_int_middle")
    
    resAll <- rbind(res,res2, res3, res4)#, res5)
    return(resAll)
  }else{
    return(NULL)
  }
}))
setnames(assessment,old=c("rn","Pr(>|t|)"),new=c("coefficient","pval"))
# FDR correction
assessment <- data.table(assessment%>%group_by(Treat,Model)%>%dplyr::mutate(adjust.p.val= p.adjust(pval, method="fdr")))
#assessment[,adjust.p.val :=p.adjust(pval,method="fdr") ] #,method="fdr"


# Subset signaificant terms
responseRelated <- assessment[adjust.p.val<0.05][ coefficient!="(Intercept)"]

# Count minimum number ofsamples that are contribuing to statistical effect at different timepoints 
repdd <- CCI%>%group_by(Day,Pair.Name,LigandPhenoCelltype,ReceptorPhenoCelltype,dynamic_class3,Treat)%>%dplyr::summarise(replic=length(TransductionMu))%>%group_by(Pair.Name,LigandPhenoCelltype,ReceptorPhenoCelltype,dynamic_class3,Treat)%>%
  dplyr::summarise(replic=min(replic))
# Get the significant effects that are identified when there are more than three samples at each timepoint
responseRelated2_0 <- merge(responseRelated[order(adjust.p.val)] ,repdd, by=c("Pair.Name","LigandPhenoCelltype","ReceptorPhenoCelltype","Treat"),all.x=T)[replic>4][order(adjust.p.val)]
responseRelated2_0[Model=="one way anova_DayStart_Response"][Treat=="CombinationRibo"]$LigandPhenoCelltype%>%table()
responseRelated2_0[Model=="one way anova_DayStart_Response"][Treat=="CombinationRibo"]$ReceptorPhenoCelltype%>%table()
responseRelated2_0[Model=="one way anova_DayStart_Response"][Treat=="CombinationRibo"] %>% select(LigandPhenoCelltype,ReceptorPhenoCelltype)%>%table()

responseRelated2_0[Model=="one way anova_DayStart_Response"][Treat=="CombinationRibo"][LigandPhenoCelltype=="Cancer cells"][ReceptorPhenoCelltype=="Macrophages"]
responseRelated[Model=="one way anova_DayStart_Response"][Treat=="CombinationRibo"][LigandPhenoCelltype=="Cancer cells"][ReceptorPhenoCelltype=="Macrophages"]
responseRelated2_0[Model=="one way anova_DayStart_Response"][Treat=="CombinationRibo"]%>%select(LigandPhenoCelltype,ReceptorPhenoCelltype)%>%table()
responseRelated2_0[Model=="one way anova_DayStart_Response"][Treat=="CombinationRibo"][LigandPhenoCelltype=="Cancer cells"][ReceptorPhenoCelltype=="Macrophages"]
responseRelated2_0[Model=="one way anova_DayStart_Response"][Treat=="CombinationRibo"][LigandPhenoCelltype=="Cancer cells"][ReceptorPhenoCelltype=="Cancer cells"]

#assessment[pval<0.05][ coefficient!="(Intercept)"]
responseRelated2_0[Model=="one way anova_DayStart_Response"][Treat=="CombinationRibo"][LigandPhenoCelltype=="Cancer cells"][ReceptorPhenoCelltype=="Macrophages"][order(pval)]%>%select()

#assessment[pval<0.05][ coefficient!="(Intercept)"][Model=="one way anova_DayStart_Response"][Treat=="CombinationRibo"][LigandPhenoCelltype=="Cancer cells"][ReceptorPhenoCelltype=="Macrophages"][order(pval)]

assessment[ coefficient!="(Intercept)"][Model=="one way anova_DayStart_Response"][Treat=="CombinationRibo"][LigandPhenoCelltype=="Cancer cells"][ReceptorPhenoCelltype=="Macrophages"][Pair.Name%in%c(
  "SERPINC1_LRP1",  "SERPINE2_LRP1","SEMA3B_NRP1", "SEMA3B_NRP2",#"SEMA4F_NRP2",
  #"FGF2_SDC2","FGF2_SDC3","COL5A3_SDC3","FGF2_FGFR1","FGF2_CD44",
  "MADCAM1_CD44" )|pval<0.05 ][order(pval)]


#save(responseRelated2_0,responseRelated,assessment,file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Cohort2 When cells communicate_communication changes AllArms/Cohort2 trendsby response AllArms.RData")
load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Cohort2 When cells communicate_communication changes AllArms/Cohort2 trendsby response AllArms.RData")
extendelistC2<-assessment[ coefficient!="(Intercept)"][Treat=="CombinationRibo"][LigandPhenoCelltype=="Cancer cells"][ReceptorPhenoCelltype=="Macrophages"][Model=="one way anova_DayStart_Response"][pval<0.05]

assessment[pval<0.05][ coefficient!="(Intercept)"][Model=="one way anova_DayStart_Response"][Treat=="CombinationRibo"][LigandPhenoCelltype=="Cancer cells"][ReceptorPhenoCelltype=="Macrophages"][order(pval)]%>%
  dplyr::select(c('Pair.Name', 'LigandPhenoCelltype', 'ReceptorPhenoCelltype','Treat' , 'coefficient','Estimate', 'Std. Error',    't value','pval','Model',         'paramNames'))
responseRelated2_0[Model=="one way anova_DayStart_Response"][Treat=="CombinationRibo"][LigandPhenoCelltype=="Cancer cells"][ReceptorPhenoCelltype=="Macrophages"][order(pval)]%>%
  dplyr::select(c('Pair.Name', 'LigandPhenoCelltype', 'ReceptorPhenoCelltype','Treat' , 'coefficient','Estimate', 'Std. Error', 't value','pval','Model', 'paramNames'))

CCI2[,scalelnTransductionMu:=scale(log(1+TransductionMu)),by=Pair.Name]

ggplot(CCI2[][Treat=="CombinationRibo"][LigandPhenoCelltype=="Cancer cells"][ReceptorPhenoCelltype=="Fibroblasts"][
  grepl("_EGF",Pair.Name) ][Day==0],
       aes(y=scalelnTransductionMu, x=Pair.Name,col=dynamic_class3,fill=dynamic_class3,group=interaction(dynamic_class3,Pair.Name)))+ theme_classic()+
  geom_point(size=4,position = position_dodge(width=0.9)) +
  geom_violin(alpha=0.5, position=position_dodge())+
  #facet_grid(LigandPhenoCelltype~ Pair.Name,scales="free")+
  labs(y="Pre-treatment cancer communication with fibroblasts")+
  theme(aspect.ratio=0.51,axis.text.x=element_text(, angle=90, vjust=0.3))+
  scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant", "Sensitive")) 

ggplot(CCI2[][Treat=="CombinationRibo"][LigandPhenoCelltype=="Cancer cells"][ReceptorPhenoCelltype=="Fibroblasts"][
  grepl("TNF",Pair.Name) ][Day==0],
  
       aes(y=scalelnTransductionMu, x=prop_change,col=dynamic_class3,group=Pair.Name))+ theme_classic()+
  geom_point(size=4,position = position_dodge(width=0.9)) +
  #geom_violin(alpha=0.5, position=position_dodge())+
  #facet_grid(LigandPhenoCelltype~ Pair.Name,scales="free")+
  labs(y="Pre-treatment cancer communication with fibroblasts")+
  theme(aspect.ratio=1,axis.text.x=element_text(, angle=90, vjust=0.3))+
  scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant", "Sensitive")) +facet_wrap(~Pair.Name)+
  geom_smooth(method="lm",se=F )+
  xlab("Fraction of tumor remaining after treatment")



ggplot(CCI2[][Treat=="CombinationRibo"][LigandPhenoCelltype=="Cancer cells"][ReceptorPhenoCelltype=="Macrophages"][Pair.Name%in%c("SERPINC1_LRP1",  "SERPINE2_LRP1","SEMA3B_NRP1", "SEMA3B_NRP2","SEMA4F_NRP2",  
                                                                                                                                  "FGF2_SDC2","FGF2_SDC3","COL5A3_SDC3","FGF2_FGFR1","FGF2_CD44","MADCAM1_CD44" ) ][Day==0],
       aes(y=scalelnTransductionMu, x=Pair.Name,col=dynamic_class3,fill=dynamic_class3,group=interaction(dynamic_class3,Pair.Name)))+ theme_classic()+
  geom_point(size=4,position = position_dodge(width=0.9)) +
  geom_violin(alpha=0.5, position=position_dodge())+
  #facet_grid(LigandPhenoCelltype~ Pair.Name,scales="free")+
  labs(y="Pre-treatment cancer communication with macrophages")+
  theme(aspect.ratio=1,axis.text.x=element_text(, angle=90, vjust=0.3))+
  scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant", "Sensitive")) 



ggplot(CCI2[][Treat=="CombinationRibo"][LigandPhenoCelltype=="Cancer cells"][ReceptorPhenoCelltype=="Macrophages"][Pair.Name%in%c("SERPINC1_LRP1",  "SERPINE2_LRP1","SEMA3B_NRP1", "SEMA3B_NRP2","SEMA4F_NRP2",  
                                                                                                                                  "FGF2_SDC2","FGF2_SDC3","COL5A3_SDC3","FGF2_FGFR1","FGF2_CD44","MADCAM1_CD44" ) ][Day==0],
       aes(y=scalelnTransductionMu, x=prop_change,col=dynamic_class3,group=Pair.Name))+ theme_classic()+
  geom_point(size=4,position = position_dodge(width=0.9)) +
  #geom_violin(alpha=0.5, position=position_dodge())+
  #facet_grid(LigandPhenoCelltype~ Pair.Name,scales="free")+
  labs(y="Pre-treatment cancer communication with macrophages")+
  theme(aspect.ratio=1,axis.text.x=element_text(, angle=90, vjust=0.3))+
  scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant", "Sensitive")) +facet_wrap(~Pair.Name)+
  geom_smooth(method="lm",se=F )+
  xlab("Fraction of tumor remaining after treatment")
  
ggplot(CCI2[][Treat=="CombinationRibo"][LigandPhenoCelltype=="Cancer cells"][ReceptorPhenoCelltype=="Macrophages"][Pair.Name%in%c(
  "SERPINC1_LRP1",  "SERPINE2_LRP1","SEMA3B_NRP1", "SEMA3B_NRP2",#"SEMA4F_NRP2",
   #"FGF2_SDC2","FGF2_SDC3","COL5A3_SDC3","FGF2_FGFR1","FGF2_CD44",
  "MADCAM1_CD44" ) ][Day==0],
       aes(y=scalelnTransductionMu, x=prop_change,group=Pair.Name))+ theme_classic()+
  geom_point(aes(col=dynamic_class3),size=4,position = position_dodge(width=0.9)) +
  #geom_violin(alpha=0.5, position=position_dodge())+
  #facet_grid(LigandPhenoCelltype~ Pair.Name,scales="free")+
  labs(y="Pre-treatment cancer communication with macrophages")+
  theme(aspect.ratio=1,axis.text.x=element_text(, angle=90, vjust=0.3))+
  scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant", "Sensitive")) +
  geom_smooth(method="lm",se=F )+facet_wrap(~Pair.Name)




extended.signiftable<- data.table(read.csv( file="~/Dropbox/Cancer_pheno_evo/data/FELINE2/Final communication files/Communication summary files by arm timepoint and resistance1/Ribo Extended Communication differences between resistant and sensitive tumors at day0.csv"))
extendelistC1<-as.character(extended.signiftable$Pair.Name)
  
overlaplist <- extendelistC1[ extendelistC1%in%extendelistC2$Pair.Name ]
#save(overlaplist, file="~/Dropbox/Cancer_pheno_evo/data/FELINE2/Final communication files/C1 and C2 Ribo Extended Communication differences between resistant and sensitive tumors at day0.RData")
CCI2[,scalelnTransductionMu:=scale(log(1+TransductionMu)),by=Pair.Name]
ggplot(CCI2[][Treat=="CombinationRibo"][LigandPhenoCelltype=="Cancer cells"][ReceptorPhenoCelltype=="Macrophages"][Pair.Name%in%overlaplist ][Day==0],
       aes(y=scalelnTransductionMu, x=Pair.Name,col=dynamic_class3,fill=dynamic_class3,group=interaction(dynamic_class3,Pair.Name)))+ theme_classic(base_size=18)+
  geom_point(size=4,position = position_dodge(width=0.9)) +
  geom_violin(alpha=0.5, position=position_dodge())+
  #facet_grid(LigandPhenoCelltype~ Pair.Name,scales="free")+
  labs(y="Pre-treatment cancer to myeloid communication")+
  theme(aspect.ratio=1,axis.text.x=element_text(, angle=90, vjust=0.3))+
  scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant", "Sensitive")) 
#pathfig <- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cohort 2 results/"
ggsave(filename=paste0(pathfig,"Cohort 2 validated Pre-treatment cancer communication with macrophages pathways initially differing in resistant and sensitive tumors under Ribo.png"), height = 8,width=8)





ggplot(CCI2[][Treat=="CombinationRibo"][LigandPhenoCelltype=="Cancer cells"][ReceptorPhenoCelltype=="Macrophages"][Pair.Name%in%extended.signiftable$Pair.Name ][Day==0],
       aes(y=scalelnTransductionMu, x=Pair.Name,col=dynamic_class3,fill=dynamic_class3,group=interaction(dynamic_class3,Pair.Name)))+ theme_classic()+
  geom_point(size=4,position = position_dodge(width=0.9)) +
  geom_violin(alpha=0.5, position=position_dodge())+
  #facet_grid(LigandPhenoCelltype~ Pair.Name,scales="free")+
  labs(y="Pre-treatment cancer communication with macrophages")+
  theme(aspect.ratio=1,axis.text.x=element_text(, angle=90, vjust=0.3))+
  scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant", "Sensitive")) 


ggplot(CCI2[][Treat=="CombinationRibo"][LigandPhenoCelltype=="Cancer cells"][ReceptorPhenoCelltype=="Macrophages"][Pair.Name%in%extended.signiftable$Pair.Name ][Day==0],
       aes(y=scalelnTransductionMu, x=prop_change,col=dynamic_class3,group=Pair.Name))+ theme_classic()+
  geom_point(size=4,position = position_dodge(width=0.9)) +
  #geom_violin(alpha=0.5, position=position_dodge())+
  #facet_grid(LigandPhenoCelltype~ Pair.Name,scales="free")+
  labs(y="Pre-treatment cancer communication with macrophages")+
  theme(aspect.ratio=1,axis.text.x=element_text(, angle=90, vjust=0.3))+
  scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant", "Sensitive")) +facet_wrap(~Pair.Name)+
  geom_smooth(method="lm",se=F )+
  xlab("Fraction of tumor remaining after treatment")






ggplot(CCI2[
  #Patient.Study.ID!="001-152"
  ][Treat=="CombinationRibo"][LigandPhenoCelltype=="Cancer cells"][ReceptorPhenoCelltype=="Macrophages"][Pair.Name%in%c(
    "SERPINC1_LRP1",  "SERPINE2_LRP1","SEMA3B_NRP1", "SEMA3B_NRP2","SEMA4F_NRP2",
    #"FGF2_SDC2","FGF2_SDC3","COL5A3_SDC3","FGF2_FGFR1","FGF2_CD44",
    "MADCAM1_CD44" ) ][Day==0]%>%group_by(Patient.Study.ID,Day,prop_change,dynamic_class3)%>%dplyr::summarise(scalelnTransductionMu1=median(scalelnTransductionMu),ucl=median(scalelnTransductionMu)+sd(scalelnTransductionMu)/sqrt(n()),lcl=median(scalelnTransductionMu)-sd(scalelnTransductionMu)/sqrt(n())),
  aes(y=scalelnTransductionMu1, x=prop_change,ymax=ucl,ymin=lcl))+ theme_classic()+
  geom_point(aes(col=dynamic_class3),size=4) +
  geom_errorbar(aes(ymax=ucl,ymin=lcl))+
  labs(y="Pre-treatment cancer communication with macrophages")+
  theme(aspect.ratio=1,axis.text.x=element_text(, angle=90, vjust=0.3))+
  scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant", "Sensitive")) +
  geom_smooth(method="gam",se=F )
#ggsave(filename=paste0(pathfig," Pre-treatment cancer communication with macrophages"," pathways initially differing in resistant and sensitive tumors under Ribo.png"), height = 12,width=15)



initdiffs <- unique(responseRelated2_0[pval<0.05][ coefficient!="(Intercept)"][Model=="one way anova_DayStart_Response"][Treat=="CombinationRibo"][] )
#initdiffs <- unique(assessment[pval<0.05][ coefficient!="(Intercept)"][Model=="one way anova_DayStart_Response"][Treat=="CombinationRibo"][] )
plotdat<-merge(CCI[Day==0], initdiffs[ReceptorPhenoCelltype=="Macrophages"][LigandPhenoCelltype%in%"Cancer cells"]%>%dplyr::select(Pair.Name,LigandPhenoCelltype, ReceptorPhenoCelltype,Treat,adjust.p.val, ),by=c("Pair.Name","LigandPhenoCelltype", "ReceptorPhenoCelltype","Treat" )) 
plotdat[,scalelnTransductionMu:=scale(log(1+TransductionMu)),by=Pair.Name]

plotdat$Pair.Name <-factor(plotdat$Pair.Name, 
        levels= data.table(plotdat%>%group_by(Pair.Name)%>%dplyr::summarise(n=max(scalelnTransductionMu)))[order(-n)]$Pair.Name )
plotdat <- plotdat%>%group_by(Pair.Name,dynamic_class3)%>%mutate(n=n())
plotdat <- data.table(plotdat%>%group_by(Pair.Name)%>%mutate(nmin= min (n)))

#plotdat$Pair.Name <-factor(plotdat$Pair.Name, 
#                           levels= c("SERPINC1_LRP1",  "SERPINE2_LRP1","SEMA3B_NRP1", "SEMA3B_NRP2","SEMA4F_NRP2",
 #                                    "FGF2_SDC2","FGF2_SDC3","COL5A3_SDC3","FGF2_FGFR1","FGF2_CD44","MADCAM1_CD44" ,"LAMB1_ITGA6","ADIPOQ_ADIPOR2", "CALR_SCARF1" ,"IL12A_IL12RB1" ,"PLAU_PLAUR") )

ggplot(plotdat[nmin>1], aes(y=scalelnTransductionMu,x=Pair.Name,col=dynamic_class3,fill=dynamic_class3))+ theme_classic()+
  geom_point(size=4,position = position_dodge(width=0.9)) + geom_violin(alpha=0.5, position=position_dodge())+
  #facet_grid(LigandPhenoCelltype~ Pair.Name,scales="free")+
  labs(y="Pre-treatment cancer communication with macrophages"   , x= "Ligand-receptor pair")+
  theme(aspect.ratio=1,axis.text.x=element_text(, angle=90, vjust=0.3))+
  scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant", "Sensitive")) 
#ggsave(filename=paste0(pathfig," C2 Pre-treatment cancer communication with macrophages"," pathways initially differing in resistant and sensitive tumors under Ribo.png"), height = 12,width=15)

ggplot(plotdat[nmin>1], aes(y=scalelnTransductionMu,x=prop_change))+ theme_classic()+
  geom_point(aes(col=dynamic_class3,fill=dynamic_class3), size=4,position = position_dodge(width=0.9)) +
  geom_smooth(method="lm",se=F)  +
  labs(y="Pre-treatment cancer communication with macrophages"   , x= "Fraction of tumor remaining")+
  theme(aspect.ratio=1,axis.text.x=element_text(, angle=90, vjust=0.3))+
  scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant", "Sensitive")) +facet_wrap(~Pair.Name)







initdiffs <- unique(responseRelated2_0[pval<0.05][ coefficient!="(Intercept)"][Model=="one way anova_DayStart_Response"][Treat=="CombinationRibo"][] )
#initdiffs <- unique(assessment[pval<0.05][ coefficient!="(Intercept)"][Model=="one way anova_DayStart_Response"][Treat=="CombinationRibo"][] )
plotdat <- merge(CCI[Day==0], allcomms,by=c("Pair.Name","LigandPhenoCelltype", "ReceptorPhenoCelltype","Treat" )) 
plotdat[,scalelnTransductionMu:=scale(log(1+TransductionMu)),by=Pair.Name]

plotdat$Pair.Name <-factor(plotdat$Pair.Name, 
                           levels= data.table(plotdat%>%group_by(Pair.Name)%>%dplyr::summarise(n=max(scalelnTransductionMu)))[order(-n)]$Pair.Name )
plotdat <- plotdat%>%group_by(Pair.Name,dynamic_class3)%>%mutate(n=n())
plotdat <- unique(data.table(plotdat%>%group_by(Pair.Name)%>%mutate(nmin= min (n))%>%dplyr::select(-adjust.p.val)))

ggplot(plotdat[nmin>1], aes(y=scalelnTransductionMu,x=prop_change))+ theme_classic()+
  geom_point(aes(col=dynamic_class3,fill=dynamic_class3), size=4,position = position_dodge(width=0.9)) +
  geom_smooth(method="lm",se=F)  +
  labs(y="Pre-treatment cancer communication with macrophages"   , x= "Fraction of tumor remaining")+
  theme(aspect.ratio=1,axis.text.x=element_text(, angle=90, vjust=0.3))+
  scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant", "Sensitive")) +facet_wrap(~Pair.Name)



ggplot(plotdat[nmin>1], aes(y=scalelnTransductionMu,x=prop_change))+ theme_classic()+
  geom_point(aes(col=dynamic_class3,fill=dynamic_class3), size=4) +
  geom_smooth(method="lm",se=F)  +
  labs(y="Pre-treatment cancer communication with macrophages"   , x= "Fraction of tumor remaining")+
  theme(aspect.ratio=1,axis.text.x=element_text( vjust=0.3))+
  scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant", "Sensitive")) 


ggplot(plotdat[nmin>1], aes(y=scalelnTransductionMu,x=prop_change))+ theme_classic()+
  geom_point(aes(col=Patient.Study.ID), size=4) +
  geom_smooth(method="gam",se=F,formula=y~s(x,k=3))  +
  labs(y="Pre-treatment cancer communication with macrophages"   , x= "Fraction of tumor remaining")+
  theme(aspect.ratio=1,legend.position="none")
  
  
ggplot(plotdat[nmin>1] %>% group_by(Patient.Study.ID,Day,prop_change,dynamic_class3) %>% dplyr::summarise(scalelnTransductionMu1=mean(scalelnTransductionMu),ucl=mean(scalelnTransductionMu)+sd(scalelnTransductionMu)/sqrt(n()),lcl=mean(scalelnTransductionMu)-sd(scalelnTransductionMu)/sqrt(n())),
  aes(y=scalelnTransductionMu1, x=prop_change,ymax=ucl,ymin=lcl))+ theme_classic()+
  geom_point(aes(col=dynamic_class3),size=4) +
  geom_errorbar(aes(ymax=ucl,ymin=lcl))+
  labs(y="Pre-treatment cancer communication with macrophages")+
  theme(aspect.ratio=1,axis.text.x=element_text(, angle=90, vjust=0.3))+
  scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant", "Sensitive")) +
  geom_smooth(method="gam",se=F )
########


CtoMpathsC1C2 <- unique(c( # Cancer to macrophage signals detected from C1
  as.character(extended.signiftable$Pair.Name),
#  "FGF9_FGFR1","FGF9_FGFR2","DKK2_LRP6","RARRES2_CCRL2", "IL11_IL11RA",   "MADCAM1_ITGB7", "GNB3_TGFBR1",   "RSPO1_LGR4","CLCF1_IL6ST","IL11_IL6ST","MADCAM1_CD44","SEMA3B_NRP1","SEMA3B_NRP2" , 
#   "SERPINC1_LRP1", "SERPINE2_LRP1","FGF2_SDC2","FGF2_SDC3","COL5A3_SDC3","FGF2_FGFR1","FGF2_CD44",
  # Cancer to macrophage signals detected from C2
   "FGF9_FGFR1","FGF9_FGFR2","DKK2_LRP6","RARRES2_CCRL2","IL11_IL11RA","MADCAM1_ITGB7", "GNB3_TGFBR1","RSPO1_LGR4","CLCF1_IL6ST","IL11_IL6ST" ))


C1commextract<-function(){
  load(file ="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Communication output merged ALL ARMS/PopulationCommunicationMerged ALL ARMS.RData" )#CCI,allphenotypes, uu,perIndiv,
   return(CCI)
}
CCI1 <- merge(C1commextract(), unique(Clin_resp_dd_classAdd%>%dplyr::select(Patient.Study.ID,prop_change,dynamic_class)), by="Patient.Study.ID")
#unique(CCI1 %>% dplyr::select(dynamic_class3,prop_change,Patient.Study.ID))[order(prop_change)]

CCI1[,scalelnTransductionMu:=scale(log(1+TransductionMu)),by=Pair.Name]
CCI1[ prop_change < (2/3) ,dynamic_class3:="Response"]
CCI1[,Cohort:="C1"]
CCI2[,Cohort:="C2"]

CCIAll_CtoMpathsC1C2 <- data.table(rbind(
CCI1[Pair.Name%in%CtoMpathsC1C2][ARM!="A"][LigandPhenoCelltype=="Cancer cells"][ReceptorPhenoCelltype=="Macrophages"][Day==0] %>% dplyr::select("Patient.Study.ID", "LigandPhenoCelltype", "ReceptorPhenoCelltype", "Day",    "Pair.Name","prop_change", "dynamic_class","dynamic_class3","ARM","Cohort","scalelnTransductionMu","TransductionMu") ,
CCI2[Pair.Name%in%CtoMpathsC1C2][ARM!="A"][LigandPhenoCelltype=="Cancer cells"][ReceptorPhenoCelltype=="Macrophages"][Day==0]  %>% dplyr::select("Patient.Study.ID", "LigandPhenoCelltype", "ReceptorPhenoCelltype", "Day",    "Pair.Name","prop_change", "dynamic_class","dynamic_class3","ARM","Cohort","scalelnTransductionMu","TransductionMu")
))


#CCIAll_CtoMpathsC1C2[ Pair.Name%in%c("MADCAM1_CD44","SEMA3B_NRP1","SEMA3B_NRP2" ,   "SERPINC1_LRP1", "SERPINE2_LRP1","FGF2_SDC2","FGF2_SDC3","COL5A3_SDC3","FGF2_FGFR1","FGF2_CD44"),   CommPredictor:="C1 predictors"   ]
CCIAll_CtoMpathsC1C2[ Pair.Name%in%as.character(extended.signiftable$Pair.Name),   CommPredictor:="C1 predictors"   ]
CCIAll_CtoMpathsC1C2[ Pair.Name%in%c("FGF9_FGFR1","FGF9_FGFR2","DKK2_LRP6","RARRES2_CCRL2", "IL11_IL11RA",   "MADCAM1_ITGB7", "GNB3_TGFBR1",   "RSPO1_LGR4","CLCF1_IL6ST","IL11_IL6ST"),   CommPredictor:="C2 predictors"   ]

CCIAll_CtoMpathsC1C2[  , validated:=F] 
CCIAll_CtoMpathsC1C2[ Pair.Name%in% overlaplist , validated:=T] 






summarCommVsGrowth<- data.table( CCIAll_CtoMpathsC1C2[! Pair.Name%in%c("FGF2_SDC2","FGF2_SDC3","COL5A3_SDC3","FGF2_FGFR1","FGF2_CD44")] %>% 
                                   group_by(Cohort,Patient.Study.ID,Day,prop_change,dynamic_class3) %>% #%>% group_by(Cohort,Patient.Study.ID,Day,prop_change,dynamic_class3) %>%
  dplyr::summarise(scalelnTransductionMu1=mean(scalelnTransductionMu),ucl=mean(scalelnTransductionMu)+sd(scalelnTransductionMu)/sqrt(n()),lcl=mean(scalelnTransductionMu)-sd(scalelnTransductionMu)/sqrt(n())))
summarCommVsGrowthC1C2<- data.table( CCIAll_CtoMpathsC1C2 [
 # ! Pair.Name%in%c("FGF2_SDC2","FGF2_SDC3","COL5A3_SDC3","FGF2_FGFR1","FGF2_CD44")
  ] %>%  
                                       group_by(Cohort,CommPredictor,Patient.Study.ID,Day,prop_change,dynamic_class3) %>% 
                                   dplyr::summarise(scalelnTransductionMu1=mean(scalelnTransductionMu),ucl=mean(scalelnTransductionMu)+sd(scalelnTransductionMu)/sqrt(n()),lcl=mean(scalelnTransductionMu)-sd(scalelnTransductionMu)/sqrt(n())))
#dplyr::summarise(scalelnTransductionMu1=mean(scalelnTransductionMu),ucl=mean(scalelnTransductionMu)+sd(scalelnTransductionMu)*1.96,lcl=mean(scalelnTransductionMu)-sd(scalelnTransductionMu)*1.96),
summarCommVsGrowth[order(scalelnTransductionMu1)]
#unique(summarCommVsGrowth %>% dplyr::select(dynamic_class3,prop_change,Patient.Study.ID))



summarCommVsGrowth<- data.table( CCIAll_CtoMpathsC1C2[validated==T] %>% 
                                   group_by(Cohort,Patient.Study.ID,Day,prop_change,dynamic_class3) %>% #%>% group_by(Cohort,Patient.Study.ID,Day,prop_change,dynamic_class3) %>%
                                   dplyr::summarise(scalelnTransductionMu1=mean(scalelnTransductionMu),ucl=mean(scalelnTransductionMu)+sd(scalelnTransductionMu)/sqrt(n()),lcl=mean(scalelnTransductionMu)-sd(scalelnTransductionMu)/sqrt(n())))
summarCommVsGrowthC1C2<- data.table( CCIAll_CtoMpathsC1C2 [validated==T
  # ! Pair.Name%in%c("FGF2_SDC2","FGF2_SDC3","COL5A3_SDC3","FGF2_FGFR1","FGF2_CD44")
  ] %>%  
    group_by(Cohort,CommPredictor,Patient.Study.ID,Day,prop_change,dynamic_class3) %>% 
    dplyr::summarise(scalelnTransductionMu1=mean(scalelnTransductionMu),ucl=mean(scalelnTransductionMu)+sd(scalelnTransductionMu)/sqrt(n()),lcl=mean(scalelnTransductionMu)-sd(scalelnTransductionMu)/sqrt(n())))
#dplyr::summarise(scalelnTransductionMu1=mean(scalelnTransductionMu),ucl=mean(scalelnTransductionMu)+sd(scalelnTransductionMu)*1.96,lcl=mean(scalelnTransductionMu)-sd(scalelnTransductionMu)*1.96),
summarCommVsGrowth[order(scalelnTransductionMu1)]
#unique(summarCommVsGrowth %>% dplyr::select(dynamic_class3,prop_change,Patient.Study.ID))

ggplot(summarCommVsGrowth%>%group_by(Cohort)%>%mutate(Mu1= mean(scalelnTransductionMu1)),
       aes(y=scalelnTransductionMu1-Mu1, x=prop_change,ymax=ucl-Mu1,ymin=lcl-Mu1))+ theme_classic()+
  geom_point(aes(col=dynamic_class3,shape=Cohort),size=4) +
  geom_errorbar()+
  labs(y="Pre-treatment cancer to myeloid communication",x="Post treatment fraction of tumor remaining")+
  theme(aspect.ratio=1,axis.text.x=element_text(, angle=90, vjust=0.3))+
  scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant", "Sensitive")) +
  geom_smooth(method="gam",se=T)
pathfig<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cohort 2 results/"
#ggsave(filename=paste0(pathfig,"Cancer to macrophage communications predict tumor response during treatment.png"))


ggplot(na.omit(summarCommVsGrowth)%>%group_by(Cohort)%>%mutate(Mu1= mean(scalelnTransductionMu1)),
       aes(x=scalelnTransductionMu1-Mu1, y=prop_change))+ theme_classic()+
  geom_point(aes(col=dynamic_class3,shape=Cohort),size=4) +
  geom_errorbarh(aes(xmax=ucl-Mu1,xmin=lcl-Mu1))+
  labs(x="Pre-treatment cancer to myeloid communication",y="Post treatment fraction of tumor remaining")+
  theme(aspect.ratio=1,axis.text.x=element_text(, angle=90, vjust=0.3))+
  scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive")) +
  scale_fill_npg(name="Tumor response", labels=c("Resistant", "Sensitive")) +
  geom_smooth(method="gam",se=T)
pathfig<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cohort 2 results/"
#ggsave(filename=paste0(pathfig,"Cancer to macrophage communications predict tumor response during treatment 2.png"))








ggplot(summarCommVsGrowth%>%group_by(Cohort)%>%mutate(Mu1= mean(scalelnTransductionMu1)),
       aes(x=scalelnTransductionMu1-Mu1, y=prop_change))+ theme_classic()+
  geom_point(aes(col=dynamic_class3,shape=Cohort),size=4) +
  labs(x="Pre-treatment cancer communication with macrophages",y="Fraction of tumor remaining")+
  theme(aspect.ratio=1,axis.text.x=element_text(, angle=90, vjust=0.3))+
  scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant", "Sensitive")) +
  geom_smooth(method="gam",se=T )



ggplot(CCIAll_CtoMpathsC1C2[! Pair.Name%in%c("FGF2_SDC2","FGF2_SDC3","COL5A3_SDC3","FGF2_FGFR1","FGF2_CD44")] ,
       aes(y=scalelnTransductionMu, x=prop_change))+ theme_classic()+
  geom_point(aes(col=Patient.Study.ID),size=4) +
  labs(y="Pre-treatment cancer communication with macrophages")+
  theme(aspect.ratio=1)+
  geom_smooth(method="gam",se=F )+facet_wrap(~Cohort)




###############


ggplot(CCI2[Patient.Study.ID!="001-152"][Treat=="CombinationRibo"][LigandPhenoCelltype=="Cancer cells"][ReceptorPhenoCelltype=="Macrophages"][Pair.Name%in%c("SERPINC1_LRP1",  "SERPINE2_LRP1","SEMA3B_NRP1", "SEMA3B_NRP2","SEMA4F_NRP2",
                                                                                                                                                             "FGF2_SDC2","FGF2_SDC3","COL5A3_SDC3","FGF2_FGFR1","FGF2_CD44","MADCAM1_CD44" ) ][Day==0],
       aes(y=scalelnTransductionMu, x=dynamic_class,group=Pair.Name))+ theme_classic()+
  geom_point(size=4,position = position_dodge(width=0.9)) +
  #geom_violin(alpha=0.5, position=position_dodge())+
  #facet_grid(LigandPhenoCelltype~ Pair.Name,scales="free")+
  labs(y="Pre-treatment cancer communication with macrophages")+
  theme(aspect.ratio=1,axis.text.x=element_text(, angle=90, vjust=0.3))+
  scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant", "Sensitive")) +facet_wrap(~Pair.Name)+
  geom_smooth(method="lm",se=F )

CCI2[Treat=="CombinationRibo"][LigandPhenoCelltype=="Cancer cells"][ReceptorPhenoCelltype=="Macrophages"][Pair.Name%in%c("FGF2_SDC3" ) ][Day==0][scalelnTransductionMu==min(scalelnTransductionMu)]
CCI2[dynamic_class=="Rebound disease"][Day==0][Pair.Name%in%c("FGF2_SDC3" ) ][Treat=="CombinationRibo"][LigandPhenoCelltype=="Cancer cells"][ReceptorPhenoCelltype=="Macrophages"]


responseRelated2_0[LigandPhenoCelltype=="Cancer cells"&ReceptorPhenoCelltype=="Fibroblasts"][Treat=="CombinationRibo"]$Pair.Name%>%unique()
responseRelated2_0[ReceptorPhenoCelltype== "Cancer cells" & LigandPhenoCelltype== "Fibroblasts"]$Pair.Name %>% unique()
sort( table( responseRelated2_0[ReceptorPhenoCelltype=="Endothelial cells"][LigandPhenoCelltype=="Cancer cells"]$Pair.Name ) )

# extract intitial differences in ribo arm
initdiffs <- unique(responseRelated2_0[Model=="one way anova_DayStart_Response"][Treat=="CombinationRibo"][] %>% dplyr::select(-dynamic_class3))
initdiffsCountLong <- data.table( reshape2::melt( table( initdiffs %>% dplyr::select(LigandPhenoCelltype ,  ReceptorPhenoCelltype) ) ) )

ggplot(initdiffsCountLong[value>1], aes(x = LigandPhenoCelltype, y = ReceptorPhenoCelltype)) + 
   geom_raster(aes(fill=value)) + theme_classic()+
   scale_fill_viridis_c("Communication pathways \n initially differing in \n resistant and sensitive tumors", option="B") +
   labs(x="Signaling cell type", y="Reciever cell type") +
  theme(axis.text.x=element_text(size=9, angle=90, vjust=0.3),
                      axis.text.y=element_text(size=9),
                      aspect.ratio=1)

initdiffs$LigandPhenoCelltype <-factor(initdiffs$LigandPhenoCelltype, 
                                       levels= data.table(initdiffsCountLong%>%group_by(LigandPhenoCelltype)%>%dplyr::summarise(n=sum(value)))[order(-n)]$LigandPhenoCelltype )
initdiffs$ReceptorPhenoCelltype <-factor(initdiffs$ReceptorPhenoCelltype, 
                                         levels= data.table(initdiffsCountLong%>%group_by(ReceptorPhenoCelltype)%>%dplyr::summarise(n=sum(value)))[order(-n)]$ReceptorPhenoCelltype )
pathfig<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/LR initial differences/"

ggplot(initdiffs, aes(fill = LigandPhenoCelltype, x = ReceptorPhenoCelltype)) + 
  geom_bar(aes()) + theme_classic()+
  scale_fill_viridis_d("Signal sending cell type", option="B",direction = -1) +
  labs(x="Signal receiving cell type", y="Number of communication pathways \n initially differing in resistant and sensitive tumors") +
  theme(axis.text.x=element_text(size=9, angle=90, vjust=0.3),
        axis.text.y=element_text(size=9),
        aspect.ratio=1)

#ggsave(filename=paste0(pathfig,"Number of communication pathways initially differing in resistant and sensitive tumors by cell type under Ribo.png"))
initdiffs[,labelLigandPhenoCelltype:=LigandPhenoCelltype]
initdiffs[LigandPhenoCelltype%in% c("Cancer cells","Normal epithelial cells"),labelLigandPhenoCelltype:="Cancer/Epithelial cells"]
initdiffs[LigandPhenoCelltype%in% c("CD4+ T cells","CD8+ T cells","B cells"),labelLigandPhenoCelltype:="Lymphocytes"]
initdiffs$labelLigandPhenoCelltype <- factor(initdiffs$labelLigandPhenoCelltype, 
                                       levels= data.table(initdiffs%>%group_by(labelLigandPhenoCelltype)%>%dplyr::summarise(n=n()))[order(-n)]$labelLigandPhenoCelltype )

ggplot(initdiffs, aes(fill = labelLigandPhenoCelltype, x = ReceptorPhenoCelltype)) + 
  geom_bar(aes()) + theme_classic(base_size=18)+
  scale_fill_viridis_d("Signal sending cell type", option="B",direction = -1) +
  labs(x="Signal receiving cell type", y="Number of communication pathways \n initially differing in resistant and sensitive tumors") +
  theme(axis.text.x=element_text(, angle=90, vjust=0.3),
        axis.text.y=element_text(),
        aspect.ratio=1)

#ggsave(filename=paste0(pathfig,"Number of communication pathways initially differing in resistant and sensitive tumors by cell type under Ribo reclass.png"))


ggplot(initdiffs[ReceptorPhenoCelltype!="B cells"], aes(fill = labelLigandPhenoCelltype, x = ReceptorPhenoCelltype)) + 
  geom_bar(aes()) + theme_classic()+
  scale_fill_viridis_d("Signal sending cell type", option="B",direction = -1) +
  labs(x="Signal receiving cell type", y="Number of communication pathways \n initially differing in resistant and sensitive tumors") +
  theme(aspect.ratio=1,
        axis.text=element_blank(),axis.title=element_blank(),legend.text=element_blank(),legend.title=element_blank(),
          legend.key.size = unit(1.5, 'cm')
        )
#ggsave(filename=paste0(pathfig,"BLANK Number of communication pathways initially differing in resistant and sensitive tumors by cell type under Ribo reclass.png"))


rr <- names(rev(sort(initdiffs$ReceptorPhenoCelltype%>%table))[1])
wch<-names(rev(sort(initdiffs$LigandPhenoCelltype%>%table))[1:3])

plotdat<-merge(CCI[Day==0], initdiffs[ReceptorPhenoCelltype=="Macrophages"][LigandPhenoCelltype%in%"Cancer cells"]%>%dplyr::select(Pair.Name,LigandPhenoCelltype, ReceptorPhenoCelltype,Treat,adjust.p.val ),by=c("Pair.Name","LigandPhenoCelltype", "ReceptorPhenoCelltype","Treat" )) 
plotdat[,scalelnTransductionMu:=scale(log(1+TransductionMu)),by=Pair.Name]

#plotdat$Pair.Name <-factor(plotdat$Pair.Name, 
 #        levels= data.table(plotdat%>%group_by(Pair.Name)%>%dplyr::summarise(n=max(scalelnTransductionMu)))[order(-n)]$Pair.Name )

plotdat$Pair.Name <-factor(plotdat$Pair.Name, 
                           levels= c("SERPINC1_LRP1",  "SERPINE2_LRP1","SEMA3B_NRP1", "SEMA3B_NRP2","SEMA4F_NRP2",
                                     "FGF2_SDC2","FGF2_SDC3","COL5A3_SDC3","FGF2_FGFR1","FGF2_CD44","MADCAM1_CD44" ,"LAMB1_ITGA6","ADIPOQ_ADIPOR2", "CALR_SCARF1" ,"IL12A_IL12RB1" ,"PLAU_PLAUR") )


ggplot(plotdat[!Pair.Name%in%c("ADIPOQ_ADIPOR2", "CALR_SCARF1" ,"IL12A_IL12RB1" ,"PLAU_PLAUR","LAMB1_ITGA6") ], aes(y=scalelnTransductionMu,x=Pair.Name,col=dynamic_class3,fill=dynamic_class3))+ theme_classic()+
  geom_point(size=4,position = position_dodge(width=0.9)) + geom_violin(alpha=0.5, position=position_dodge())+
  #facet_grid(LigandPhenoCelltype~ Pair.Name,scales="free")+
  labs(y="Pre-treatment cancer communication with macrophages"   , x= "Ligand-receptor pair")+
  theme(aspect.ratio=1,axis.text.x=element_text(, angle=90, vjust=0.3))+
  scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant", "Sensitive")) 
#ggsave(filename=paste0(pathfig," Pre-treatment cancer communication with macrophages"," pathways initially differing in resistant and sensitive tumors under Ribo.png"), height = 12,width=15)

ggplot(plotdat[!Pair.Name%in%c("ADIPOQ_ADIPOR2", "CALR_SCARF1" ,"IL12A_IL12RB1" ,"PLAU_PLAUR","LAMB1_ITGA6") ], aes(y=scalelnTransductionMu,x=Pair.Name,col=dynamic_class3,fill=dynamic_class3))+ theme_classic()+
  geom_point(size=4,position = position_dodge(width=0.9)) + geom_violin(alpha=0.5, position=position_dodge())+
  #facet_grid(LigandPhenoCelltype~ Pair.Name,scales="free")+
  labs(y="Pre-treatment cancer communication with macrophages"   , x= "Ligand-receptor pair")+
  theme(aspect.ratio=1,axis.text.x=element_text(, angle=90, vjust=0.3))+
  scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))+coord_flip() 
ggsave(filename=paste0(pathfig," Flip Pre-treatment cancer communication with macrophages"," pathways initially differing in resistant and sensitive tumors under Ribo.png"), height = 12,width=15)


ggplot(plotdat[!Pair.Name%in%c("ADIPOQ_ADIPOR2", "CALR_SCARF1" ,"IL12A_IL12RB1" ,"PLAU_PLAUR","LAMB1_ITGA6") ], aes(y=scalelnTransductionMu,x=Pair.Name,col=dynamic_class3,fill=dynamic_class3))+ theme_classic()+
  geom_point(size=4,position = position_dodge(width=0.9)) + geom_violin(alpha=0.5, position=position_dodge())+
  #facet_grid(LigandPhenoCelltype~ Pair.Name,scales="free")+
  labs(y="Pre-treatment cancer communication with macrophages"   , x= "Ligand-receptor pair")+
  scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant", "Sensitive")) + 
  theme(aspect.ratio=1,legend.position="none",
      axis.text=element_blank(),axis.title=element_blank(),legend.text=element_blank(),legend.title=element_blank(),
      legend.key.size = unit(1, 'cm')
)
#ggsave(filename=paste0(pathfig,"BLANK Number of communication pathways initially differing in resistant and sensitive tumors by cell type under Ribo reclass.png"))
ggsave(filename=paste0(pathfig,"BLANK Pre-treatment cancer communication with macrophages"," pathways initially differing in resistant and sensitive tumors under Ribo2.png"), height = 12,width=15)


ggplot(plotdat[!Pair.Name%in%c("ADIPOQ_ADIPOR2", "CALR_SCARF1" ,"IL12A_IL12RB1" ,"PLAU_PLAUR","LAMB1_ITGA6") ], aes(y=scalelnTransductionMu,x=Pair.Name,col=dynamic_class3,fill=dynamic_class3))+ theme_classic()+
  geom_point(size=4,position = position_dodge(width=0.9)) + geom_violin(alpha=0.5, position=position_dodge())+
  #facet_grid(LigandPhenoCelltype~ Pair.Name,scales="free")+
  labs(y="Pre-treatment cancer communication with macrophages"   , x= "Ligand-receptor pair")+
  scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant", "Sensitive")) + 
  theme(aspect.ratio=1,legend.position="none",
        axis.text=element_blank(),axis.title=element_blank(),legend.text=element_blank(),legend.title=element_blank(),
        legend.key.size = unit(1, 'cm')
  )+coord_flip()
#ggsave(filename=paste0(pathfig,"BLANK Number of communication pathways initially differing in resistant and sensitive tumors by cell type under Ribo reclass.png"))
ggsave(filename=paste0(pathfig,"BLANK Flip Pre-treatment cancer communication with macrophages"," pathways initially differing in resistant and sensitive tumors under Ribo.png"), height = 12,width=15)



# end data
initdiffs <- unique(responseRelated[Model=="one way anova_DayEnd_Response"][Treat=="CombinationRibo"][])
initdiffsCountLong <- data.table( reshape2::melt( table( initdiffs%>% dplyr::select(LigandPhenoCelltype ,  ReceptorPhenoCelltype) ) ) )

initdiffs$LigandPhenoCelltype <-factor(initdiffs$LigandPhenoCelltype, 
                                       levels= data.table(initdiffsCountLong%>%group_by(LigandPhenoCelltype)%>%dplyr::summarise(n=sum(value)))[order(-n)]$LigandPhenoCelltype )
initdiffs$ReceptorPhenoCelltype <-factor(initdiffs$ReceptorPhenoCelltype, 
                                         levels= data.table(initdiffsCountLong%>%group_by(ReceptorPhenoCelltype)%>%dplyr::summarise(n=sum(value)))[order(-n)]$ReceptorPhenoCelltype )
initdiffs[,labelLigandPhenoCelltype:=LigandPhenoCelltype]
initdiffs[LigandPhenoCelltype%in% c("Cancer cells","Normal epithelial cells"),labelLigandPhenoCelltype:="Cancer/Epithelial cells"]
initdiffs[LigandPhenoCelltype%in% c("CD4+ T cells","CD8+ T cells","B cells"),labelLigandPhenoCelltype:="Lymphocytes"]
initdiffs$labelLigandPhenoCelltype <-factor(initdiffs$labelLigandPhenoCelltype, 
                                            levels= data.table(initdiffs%>%group_by(labelLigandPhenoCelltype)%>%dplyr::summarise(n=n()))[order(-n)]$labelLigandPhenoCelltype )

# ggplot(initdiffsCountLong[value>1], aes(x = LigandPhenoCelltype, y = ReceptorPhenoCelltype)) + 
#   geom_raster(aes(fill=value)) + theme_classic()+
#   scale_fill_viridis_c("Communication pathways \n initially differing in \n resistant and sensitive tumors", option="B") +
#   labs(x="Signaling cell type", y="Reciever cell type") +
#  theme(axis.text.x=element_text(size=9, angle=90, vjust=0.3),
#                      axis.text.y=element_text(size=9),
#                      aspect.ratio=1)

ggplot(initdiffs, aes(fill = labelLigandPhenoCelltype, x = ReceptorPhenoCelltype)) + 
  geom_bar(aes()) + theme_classic(base_size=18)+
  scale_fill_viridis_d("Signal sending cell type", option="B",direction = -1) +
  labs(x="Signal receiving cell type", y="Number of communication pathways \n initially differing in resistant and sensitive tumors") +
  theme(axis.text.x=element_text(, angle=90, vjust=0.3),
        axis.text.y=element_text(),
        aspect.ratio=1)
#ggsave(filename=paste0(pathfig,"Number of communication pathways finally differing in resistant and sensitive tumors by cell type under Ribo reclass.png"))

ggplot(initdiffs[ReceptorPhenoCelltype!="B cells"], aes(fill = labelLigandPhenoCelltype, x = ReceptorPhenoCelltype)) + 
  geom_bar(aes()) + theme_classic()+
  scale_fill_viridis_d("Signal sending cell type", option="B",direction = -1) +
  labs(x="Signal receiving cell type", y="Number of communication pathways \n initially differing in resistant and sensitive tumors") +
  theme(aspect.ratio=1,
        axis.text=element_blank(),axis.title=element_blank(),legend.text=element_blank(),legend.title=element_blank(),
        legend.key.size = unit(1.5, 'cm')
  )

#ggsave(filename=paste0(pathfig,"BLANK Number of communication pathways finally differing in resistant and sensitive tumors by cell type under Ribo reclass.png"))

initdiffs[ReceptorPhenoCelltype=="Macrophages"][LigandPhenoCelltype=="Fibroblasts"]
initdiffs[ReceptorPhenoCelltype=="Macrophages"][LigandPhenoCelltype=="Cancer cells"]
plotdat<-merge(CCI[Day==180], initdiffs[ReceptorPhenoCelltype=="Macrophages"][LigandPhenoCelltype%in%"Fibroblasts"]%>%dplyr::select(Pair.Name,LigandPhenoCelltype, ReceptorPhenoCelltype,Treat,adjust.p.val ),by=c("Pair.Name","LigandPhenoCelltype", "ReceptorPhenoCelltype","Treat" )) 
plotdat<-merge(CCI[], initdiffs[ReceptorPhenoCelltype=="Endothelial cells"][LigandPhenoCelltype%in%"Cancer cells"]%>%dplyr::select(Pair.Name,LigandPhenoCelltype, ReceptorPhenoCelltype,Treat,adjust.p.val ),by=c("Pair.Name","LigandPhenoCelltype", "ReceptorPhenoCelltype","Treat" )) 

plotdat<-merge(CCI[Day==180], initdiffs[ReceptorPhenoCelltype=="Macrophages"][LigandPhenoCelltype%in%"Cancer cells"]%>%dplyr::select(Pair.Name,LigandPhenoCelltype, ReceptorPhenoCelltype,Treat,adjust.p.val ),by=c("Pair.Name","LigandPhenoCelltype", "ReceptorPhenoCelltype","Treat" )) 
plotdat[,scalelnTransductionMu:=scale(log(1+TransductionMu)),by=Pair.Name]

#plotdat$Pair.Name <-factor(plotdat$Pair.Name, 
#        levels= data.table(plotdat%>%group_by(Pair.Name)%>%dplyr::summarise(n=max(scalelnTransductionMu)))[order(-n)]$Pair.Name )

plotdat$Pair.Name <-factor(plotdat$Pair.Name, 
                           levels= c("MST1_MST1R",  "CCL5_CCR5","CCL19_CCR7","GNAI2_CCR5", "TNFSF14_LTBR",
                                     "EDN1_EDNRB","NXPH1_NRXN2","MDK_SDC1","SORBS1_INSR") )


ggplot(plotdat[!Pair.Name%in%c("MST1_MST1R","NXPH1_NRXN2","EDN1_EDNRB") ], aes(y=scalelnTransductionMu,x=Pair.Name,col=dynamic_class3,fill=dynamic_class3))+ theme_classic()+
  geom_point(size=4,position = position_dodge(width=0.9)) + geom_violin(alpha=0.5, position=position_dodge())+
  #facet_grid(LigandPhenoCelltype~ Pair.Name,scales="free")+
  labs(y="Pre-treatment cancer communication with macrophages"   , x= "Ligand-receptor pair")+
  theme(aspect.ratio=1,axis.text.x=element_text(, angle=90, vjust=0.3))+
  scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant", "Sensitive")) 
ggsave(filename=paste0(pathfig," Post-treatment cancer communication with macrophages"," pathways finally differing in resistant and sensitive tumors under Ribo.png"), height = 12,width=15)

ggplot(plotdat[!Pair.Name%in%c("MST1_MST1R","NXPH1_NRXN2","EDN1_EDNRB") ], aes(y=scalelnTransductionMu,x=Pair.Name,col=dynamic_class3,fill=dynamic_class3))+ theme_classic()+
  geom_point(size=4,position = position_dodge(width=0.9)) + geom_violin(alpha=0.5, position=position_dodge())+
  #facet_grid(LigandPhenoCelltype~ Pair.Name,scales="free")+
  labs(y="Pre-treatment cancer communication with macrophages"   , x= "Ligand-receptor pair")+
  theme(aspect.ratio=1,axis.text.x=element_text(, angle=90, vjust=0.3))+
  scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))+coord_flip() 
ggsave(filename=paste0(pathfig," Flip Post-treatment cancer communication with macrophages"," pathways finally differing in resistant and sensitive tumors under Ribo.png"), height = 12,width=15)



ggplot(plotdat[!Pair.Name%in%c("MST1_MST1R","NXPH1_NRXN2","EDN1_EDNRB") ], aes(y=scalelnTransductionMu,x=Pair.Name,col=dynamic_class3,fill=dynamic_class3))+ theme_classic()+
  geom_point(size=4,position = position_dodge(width=0.9)) + geom_violin(alpha=0.5, position=position_dodge())+
  #facet_grid(LigandPhenoCelltype~ Pair.Name,scales="free")+
  labs(y="Pre-treatment cancer communication with macrophages"   , x= "Ligand-receptor pair")+
  scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant", "Sensitive")) + 
  theme(aspect.ratio=1,legend.position="none",
        axis.text=element_blank(),axis.title=element_blank(),legend.text=element_blank(),legend.title=element_blank(),
        legend.key.size = unit(1, 'cm')
  )
#ggsave(filename=paste0(pathfig,"BLANK Number of communication pathways initially differing in resistant and sensitive tumors by cell type under Ribo reclass.png"))
ggsave(filename=paste0(pathfig,"BLANK Post-treatment cancer communication with macrophages"," pathways Finally differing in resistant and sensitive tumors under Ribo2.png"), height = 12,width=15)





tt<- wch[2]
for(tt in wch){
  ggplot(merge(CCI[], initdiffs[ReceptorPhenoCelltype=="Macrophages"][LigandPhenoCelltype%in%tt]%>%dplyr::select(Pair.Name,LigandPhenoCelltype, ReceptorPhenoCelltype,Treat ),by=c("Pair.Name","LigandPhenoCelltype", "ReceptorPhenoCelltype","Treat" )) ,
         aes(y=log(1+TransductionMu),x=log(1+Day),col=dynamic_class3))+ theme_classic()+
    geom_point()+
    #facet_grid(LigandPhenoCelltype~ Pair.Name,scales="free")+
    facet_wrap(~ Pair.Name,scales="free")+
    geom_smooth(method = "lm") + theme(aspect.ratio=1) +
    labs(y=paste0("Macrophage communication from ",tt)   , x= "Day")+
    scale_x_continuous(breaks=log(1+c(0,14,180)), labels=c(0,14,180) )+
    scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))
  
  ggsave(filename=paste0(pathfig,"Macrophage communication from ",tt," pathways initially differing in resistant and sensitive tumors under Ribo.png"), height = 12,width=15)
}





# extract intitial differences in letrozole arm
initdiffs<-unique(responseRelated2_0[Model=="one way anova_DayStart_Response"][Treat!="CombinationRibo"][] %>% dplyr::select(-dynamic_class3))
initdiffsCountLong <- data.table( reshape2::melt( table( initdiffs%>% dplyr::select(LigandPhenoCelltype ,  ReceptorPhenoCelltype) ) ) )
initdiffs$LigandPhenoCelltype <-factor(initdiffs$LigandPhenoCelltype, 
                                       levels= data.table(initdiffsCountLong%>%group_by(LigandPhenoCelltype)%>%dplyr::summarise(n=sum(value)))[order(-n)]$LigandPhenoCelltype )
initdiffs$ReceptorPhenoCelltype <-factor(initdiffs$ReceptorPhenoCelltype, 
                                         levels= data.table(initdiffsCountLong%>%group_by(ReceptorPhenoCelltype)%>%dplyr::summarise(n=sum(value)))[order(-n)]$ReceptorPhenoCelltype )
pathfig<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/LR initial differences/"

ggplot(initdiffs, aes(fill = LigandPhenoCelltype, x = ReceptorPhenoCelltype)) + 
  geom_bar(aes()) + theme_classic()+
  scale_fill_viridis_d("Signal sending cell type", option="B",direction = -1) +
  labs(x="Signal receiving cell type", y="Number of communication pathways \n initially differing in resistant and sensitive tumors") +
  theme(axis.text.x=element_text(size=9, angle=90, vjust=0.3),
        axis.text.y=element_text(size=9),
        aspect.ratio=1)

#ggsave(filename=paste0(pathfig,"Number of communication pathways initially differing in resistant and sensitive tumors by cell type under Letrozole.png"))


rr<-names(rev(sort(initdiffs$ReceptorPhenoCelltype%>%table))[1])
wch<-names(rev(sort(initdiffs$LigandPhenoCelltype%>%table))[1:3])

tt<- wch[2]
for(tt in wch){
  ggplot(merge(CCI[], initdiffs[ReceptorPhenoCelltype=="Fibroblasts"][LigandPhenoCelltype%in%tt]%>%dplyr::select(Pair.Name,LigandPhenoCelltype, ReceptorPhenoCelltype,Treat ),by=c("Pair.Name","LigandPhenoCelltype", "ReceptorPhenoCelltype","Treat" )) ,
         aes(y=log(1+TransductionMu),x=log(1+Day),col=dynamic_class3))+ theme_classic()+
    geom_point()+
    #facet_grid(LigandPhenoCelltype~ Pair.Name,scales="free")+
    facet_wrap(~ Pair.Name,scales="free")+
    geom_smooth(method = "lm") + theme(aspect.ratio=1) +
    labs(y=paste0("Fibroblast communication from ",tt)   , x= "Day")+
    scale_x_continuous(breaks=log(1+c(0,14,180)), labels=c(0,14,180) )+
    scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))
  
  ggsave(filename=paste0(pathfig,"Fibroblast communication from ",tt," pathways initially differing in resistant and sensitive tumors under Letrozole.png.png"), height = 12,width=15)
}


Macrophage_frac <- data.table( allphenotypes%>%group_by(Patient.Study.ID,Day,dynamic_class3,ARM)%>%summarise(frac=sum(Celltype=="Macrophages")/n()) )
ggplot( Macrophage_frac[ARM!="A"] ,aes(y=(frac), x=log(1+Day), col=dynamic_class3) ) + geom_jitter(size=2,width=0.1)+ scale_x_continuous(breaks=log(c(0,14,180)), labels=c(0,14,180) )+theme_classic()+
  labs(y="Fraction of macrophages",x="Day")

#####
#CCI[,scaleTransduction:=exp(scale(log(TransductionMu),center=F))]#,by=c("Pair.Name")] # CCI[, scaleTransduction := exp((log(TransductionMu)-centVal)/scalsd) ]
CCI[,scaleTransduction:=exp(scale(log(TransductionMu),center=F)),by=c("Pair.Name")] # CCI[, scaleTransduction := exp((log(TransductionMu)-centVal)/scalsd) ]

# summaryInteraction<-data.table(CCI[Pair.Name%in%summarytable$key_]%>%
#                                  group_by(LigandPhenoCelltype,ReceptorPhenoCelltype,Day,dynamic_class3,Treat) %>%
#                                  dplyr::summarise(muln_scaleTransduction=median(scaleTransduction)))
summaryInteraction<-data.table(CCI%>%
                                 group_by(Patient.Study.ID,LigandPhenoCelltype,ReceptorPhenoCelltype,Day,dynamic_class3,ARM,Treat) %>%
                                 dplyr::summarise(muln_scaleTransduction=mean(log(1+scaleTransduction))))

average_muln_scaleTransduction <- data.table(summaryInteraction %>% 
                                               group_by(LigandPhenoCelltype,ReceptorPhenoCelltype,Day,dynamic_class3,Treat) %>%
                                               dplyr::summarise(average_muln_scaleTransduction=mean(muln_scaleTransduction)))
average_muln_scaleTransduction$LigandPhenoCelltype <-factor(average_muln_scaleTransduction$LigandPhenoCelltype, 
                                                            levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Normal epithelial cells","Cancer cells","B cells","Macrophages","CD8+ T cells","CD4+ T cells"))
average_muln_scaleTransduction$ReceptorPhenoCelltype <-factor(average_muln_scaleTransduction$ReceptorPhenoCelltype, 
                                                              levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Normal epithelial cells","Cancer cells","B cells","Macrophages","CD8+ T cells","CD4+ T cells") )

initState <- data.table(average_muln_scaleTransduction[Day==0]%>%
                          group_by(LigandPhenoCelltype,ReceptorPhenoCelltype)%>% #Treat, dynamic_class3
                          dplyr::summarise(intiState=mean(average_muln_scaleTransduction)) )

average_muln_scaleTransduction2 <- merge(average_muln_scaleTransduction,initState,by=c("LigandPhenoCelltype", "ReceptorPhenoCelltype")) #,"Treat","dynamic_class3"
average_muln_scaleTransduction2[ , lnfoldchange:= (average_muln_scaleTransduction-intiState) ] 



ggplot( average_muln_scaleTransduction2[][Treat=="CombinationRibo"][LigandPhenoCelltype!="B cells"][ReceptorPhenoCelltype!="B cells"] ,aes(LigandPhenoCelltype,ReceptorPhenoCelltype,fill=lnfoldchange  )) + geom_tile()+
  facet_wrap(dynamic_class3~Day)+scale_fill_viridis_c(name="Log fold change \n in communication \n relative to baseline",option="B") + theme_classic()+
  labs(y="Signal receiver",x="Signal sender")+theme(aspect.ratio=1,axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
#ggsave(filename = "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Strength communication cell types AllArms/signaling strength between cell types under combination ribo therapy.png")

ggplot( average_muln_scaleTransduction2[][Treat=="CombinationRibo"][LigandPhenoCelltype!="B cells"][ReceptorPhenoCelltype!="B cells"] ,aes(LigandPhenoCelltype,ReceptorPhenoCelltype,fill=lnfoldchange  )) + geom_tile()+
  facet_wrap(dynamic_class3~Day)+scale_fill_viridis_c(name="Log fold change \n in communication \n relative to baseline",option="B") + theme_classic()+
  labs(y="Signal receiver",x="Signal sender")+theme(aspect.ratio=1)+
  theme(axis.title=element_blank(),  axis.text=element_blank() ,
        legend.title = element_blank(),legend.text = element_blank(),
        strip.background = element_blank(),strip.text.x = element_blank() ,
  legend.key.size = unit(1, 'cm'))
#ggsave(filename = "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Strength communication cell types AllArms/BLANK2_signaling strength between cell types under combination ribo therapy.png")



ggplot( average_muln_scaleTransduction2[][Treat=="LetrozoleAlone"][LigandPhenoCelltype!="B cells"][ReceptorPhenoCelltype!="B cells"] ,aes(LigandPhenoCelltype,ReceptorPhenoCelltype,fill=lnfoldchange  ))+geom_tile()+
  facet_wrap(dynamic_class3~Day)+scale_fill_viridis_c(name="Log fold change \n in communication \n relative to baseline",option="B") + theme_classic()+
  labs(y="Signal receiver",x="Signal sender")+theme(aspect.ratio=1,axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#ggsave(filename = "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Strength communication cell types AllArms/signaling strength between cell types under letrozole alone therapy.png")


# 
# 
# 
# 
# ggplot( average_muln_scaleTransduction2[Day!=0][Treat=="CombinationRibo"] ,aes(LigandPhenoCelltype,ReceptorPhenoCelltype,fill=log(average_muln_scaleTransduction/intiState)  ))+geom_tile()+
#   facet_wrap(dynamic_class3~Day)+scale_fill_viridis_c(option="A")
# 
# 
# ggplot( average_muln_scaleTransduction[Treat=="CombinationRibo"] ,aes(LigandPhenoCelltype,ReceptorPhenoCelltype,fill=log(average_muln_scaleTransduction)  ))+geom_tile()+
#   facet_wrap(dynamic_class3~Day)+scale_fill_viridis_c()
# 
# #ggplot( CCI[Treat=="CombinationRibo"][Day==180][dynamic_class3=="Non-response"], aes(y=log(scaleTransduction), x=ReceptorPhenoCelltype)) + geom_violin()
# 
# ggplot( CCI[Treat=="CombinationRibo"][Day==180][dynamic_class3=="Non-response"][Pair.Name=="EFNB2_EPHA4"],aes(y=log(TransductionMu), x=ReceptorPhenoCelltype))+geom_violin()
# 
# 
# # Scale each ligand-receptor pair across samples and cell types
# CCI[,scaleTransduction:=exp(scale(log(TransductionMu),center=F))]#,by=c("Pair.Name")] # CCI[, scaleTransduction := exp((log(TransductionMu)-centVal)/scalsd) ]
# summaryInteraction<-data.table(CCI[Pair.Name%in%summarytable$key_]%>%
#                                  group_by(LigandPhenoCelltype,ReceptorPhenoCelltype,Day,dynamic_class3,Treat) %>%
#                                  dplyr::summarise(muln_scaleTransduction=median(scaleTransduction)))
# average_muln_scaleTransduction <- data.table(summaryInteraction %>% 
#                                                group_by(LigandPhenoCelltype,ReceptorPhenoCelltype,Day,dynamic_class3,Treat) %>%
#                                                dplyr::summarise(average_muln_scaleTransduction=mean(muln_scaleTransduction)))
# average_muln_scaleTransduction$LigandPhenoCelltype <-factor(average_muln_scaleTransduction$LigandPhenoCelltype, 
#                                                             levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Normal epithelial cells","Cancer cells","Macrophages","B cells","CD4+ T cells","CD8+ T cells"))
# average_muln_scaleTransduction$ReceptorPhenoCelltype <-factor(average_muln_scaleTransduction$ReceptorPhenoCelltype, 
#                                                               levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Normal epithelial cells","Cancer cells","Macrophages","B cells","CD4+ T cells","CD8+ T cells") )
# 
# ggplot( average_muln_scaleTransduction[Day==180][Treat=="CombinationRibo"] ,aes(LigandPhenoCelltype,ReceptorPhenoCelltype,fill=log(average_muln_scaleTransduction)  ))+geom_tile()+
#   facet_wrap(dynamic_class3~Day)+scale_fill_viridis_c()
# 
# 
# 
# ggplot(merge(CCI[], initdiffs[ReceptorPhenoCelltype=="Macrophages"][LigandPhenoCelltype%in%tt]%>%dplyr::select(Pair.Name,LigandPhenoCelltype, ReceptorPhenoCelltype,Treat ),by=c("Pair.Name","LigandPhenoCelltype", "ReceptorPhenoCelltype","Treat" )) ,
#        aes(y=log(1+scaleTransduction),x=log(1+Day),col=dynamic_class3))+ theme_classic()+
#   geom_point()+
#   #facet_grid(LigandPhenoCelltype~ Pair.Name,scales="free")+
#   facet_wrap(~ Pair.Name,scales="free")+
#   geom_smooth(method = "lm") + theme(aspect.ratio=1) +
#   labs(y=paste0("Macrophage communication from ",tt)   , x= "Day")+
#   scale_x_continuous(breaks=log(1+c(0,14,180)), labels=c(0,14,180) )+
#   scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))
# 
# 
# 
# 
# initdiffs<-unique(responseRelated2_0[Model=="one way anova_DayEnd_Response"][Treat=="CombinationRibo"] %>% dplyr::select(-dynamic_class3))
# initdiffsCountLong <- data.table( reshape2::melt( table( initdiffs%>% dplyr::select(LigandPhenoCelltype ,  ReceptorPhenoCelltype) ) ) )
# 
# ggplot(initdiffsCountLong[value>1], aes(x = LigandPhenoCelltype, y = ReceptorPhenoCelltype)) + 
#   geom_raster(aes(fill=value)) + theme_classic()+
#   scale_fill_viridis_c("Communication pathways \n differing post-treatment in \n resistant and sensitive tumors", option="B") +
#   labs(x="Signaling cell type", y="Reciever cell type") +
#   theme(axis.text.x=element_text(size=9, angle=90, vjust=0.3),
#         axis.text.y=element_text(size=9),
#         aspect.ratio=1)
# 
# 
# 
# 
responseRelated2 <- responseRelated2_0[grepl("NtoR",paramNames)|(grepl("Non-response",paramNames)&dynamic_class3=="Non-response") |(grepl("Response",paramNames)&dynamic_class3=="Response")]
# 
# # extract significant trends and create a excel readable format to write
print.this <- responseRelated2[][grep("slope",paramNames)]%>%dplyr::select(-c(coefficient, replic))
# #print.this[]$paramNames%>%unique
print.this[paramNames=="Response_slope",paramNames:="signif trend for responder tumors"]
print.this[paramNames=="Non-response_slope",paramNames:="signif trend for non responder tumors"]
print.this[paramNames=="NtoR_change_slope",paramNames:="signif difference in trend between responder and non responder tumors"]
# #print.this[]$paramNames%>%unique
# 
print.this <- print.this[order(Treat,ReceptorPhenoCelltype,adjust.p.val)]#[ReceptorPhenoCelltype=="Cancer cells"]
# write.csv(print.this,file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Trends in communication/Trends in communication by response ALLARMS.csv")
# 
# 
print.thisI <- responseRelated2[][grep("NtoR_change_int_start",paramNames)]%>%dplyr::select(-c(coefficient,Model, replic))
#responseRelated2[]$paramNames%>%unique
print.thisI[paramNames=="NtoR_change_int_start",paramNames:="signif initial difference in communication between responder and non responder tumors"]
# 
print.thisI <- unique(print.thisI[order(Treat,ReceptorPhenoCelltype,adjust.p.val)]%>%dplyr::select(
   c("ReceptorPhenoCelltype","LigandPhenoCelltype","Pair.Name","paramNames",  "Estimate","Std. Error","t value","pval","adjust.p.val")))#[ReceptorPhenoCelltype=="Cancer cells"]
#write.csv(print.thisI ,file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Initial differences in communication/Initial differences in communication by responseALLARMS.csv")
# print.thisI[order(ReceptorPhenoCelltype,adjust.p.val)][ReceptorPhenoCelltype=="Cancer cells"]
# 
# 
# #tab_signal_changes <- responseRelated2[!grep("Intercept",coefficient)] [!grepl("change", paramNames)]%>%dplyr::select(LigandPhenoCelltype,ReceptorPhenoCelltype) %>%table()
tab_signal_changesRN <- responseRelated2 [grepl("change", paramNames)]%>%dplyr::select(LigandPhenoCelltype,ReceptorPhenoCelltype,Treat) %>%table()
ggplot( data.table(tab_signal_changesRN)[Treat=="CombinationRibo"] ,aes( y=LigandPhenoCelltype,x=ReceptorPhenoCelltype ))+geom_tile(col="white",fill="white")+theme_classic()+
   geom_tile( data=data.table(tab_signal_changesRN)[N>3][Treat=="CombinationRibo"] ,aes( x=LigandPhenoCelltype,y=ReceptorPhenoCelltype , col=N, fill=(N))) + 
   theme(aspect.ratio=1,axis.text.x = element_text(angle = 90))+
   labs(y="Receiver cell type",x="Signaling cell type") +
   scale_fill_continuous(name="Number of \n signalling trends \n differing by \n tumor response")  +
   scale_color_continuous(name="Number of \n signalling trends \n differing by \n tumor response")#+
 # facet_wrap(~Treat)
# 
ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/LR trends AllArms/Number of signal trend differences between response groups AllARMS.png",
        width=6,height=6)
# 
# 
# #tab_signal_changes <- responseRelated2[!grep("Intercept",coefficient)] [!grepl("change", paramNames)]%>%dplyr::select(LigandPhenoCelltype,ReceptorPhenoCelltype) %>%table()
tab_signal_changes <- data.table(responseRelated2[!grep("Intercept",coefficient)] [!grepl("change", paramNames)]%>%group_by(LigandPhenoCelltype,ReceptorPhenoCelltype,Treat)%>%dplyr::summarise(N=length(unique(Pair.Name)),meanEffectSize=mean(abs(Estimate ))) )
ggplot( data.table(tab_signal_changes)[Treat=="CombinationRibo"] ,aes( y=LigandPhenoCelltype,x=ReceptorPhenoCelltype ))+geom_tile(col="white",fill="white")+theme_classic()+
   geom_tile( data=data.table(tab_signal_changes)[N>3][Treat=="CombinationRibo"] ,aes( x=LigandPhenoCelltype,y=ReceptorPhenoCelltype , col=N, fill=N)) + 
   theme(aspect.ratio=1,axis.text.x = element_text(angle = 90))+
   labs(y="Receiver cell type",x="Signaling cell type") +
   scale_fill_continuous(name="Number of \n signalling \n trends")  +
   scale_color_continuous(name="Number of \n signalling \n trends") 
ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/LR trends AllArms/Number of signal trends AllARMS.png",
        width=6,height=6)
# 
ggplot( data.table(tab_signal_changes)[Treat=="CombinationRibo"]%>%group_by(LigandPhenoCelltype)%>%summarise(N=sum(N)) ,aes( y=N ,x=LigandPhenoCelltype ))+geom_point()+theme_classic()+
   theme(aspect.ratio=1,axis.text.x = element_text(angle = 90))+
   labs(x="Signaling cell type",y="Number of \n signalling \n changes") +
   scale_fill_continuous(name="Number of \n signalling \n changes")  +
   scale_color_continuous(name="Number of \n signalling \n changes") 
ggplot( data.table(tab_signal_changes)[Treat=="CombinationRibo"]%>%group_by(ReceptorPhenoCelltype)%>%summarise(N=sum(N)) ,aes( y=N ,x=ReceptorPhenoCelltype ))+geom_point()+theme_classic()+
  theme(aspect.ratio=1,axis.text.x = element_text(angle = 90))+
  labs(x="Receiver cell type",y="Number of \n signalling \n changes") +
  scale_fill_continuous(name="Number of \n signalling \n changes")  +
  scale_color_continuous(name="Number of \n signalling \n changes") 

ggplot( data.table(tab_signal_changesRN)[Treat=="CombinationRibo"]%>%group_by(LigandPhenoCelltype)%>%summarise(N=sum(N)) ,aes( y=N ,x=LigandPhenoCelltype ))+geom_point()+theme_classic()+
  theme(aspect.ratio=1,axis.text.x = element_text(angle = 90))+
  labs(x="Signaling cell type",y="Number of \n signalling \n changes") +
  scale_fill_continuous(name="Number of \n signalling \n changes")  +
  scale_color_continuous(name="Number of \n signalling \n changes") 
ggplot( data.table(tab_signal_changesRN)[Treat=="CombinationRibo"]%>%group_by(ReceptorPhenoCelltype)%>%summarise(N=sum(N)) ,aes( y=N ,x=ReceptorPhenoCelltype ))+geom_point()+theme_classic()+
  theme(aspect.ratio=1,axis.text.x = element_text(angle = 90))+
  labs(x="Receiver cell type",y="Number of \n signalling \n changes") +
  scale_fill_continuous(name="Number of \n signalling \n changes")  +
  scale_color_continuous(name="Number of \n signalling \n changes") 
#responseRelated2 [grepl("change", paramNames)][ReceptorPhenoCelltype=="Fibroblasts"]
 

# #assessment[Pair.Name==unique_combs[1]$Pair.Name][LigandPhenoCelltype==unique_combs[1]$LigandPhenoCelltype  & ReceptorPhenoCelltype==unique_combs[1]$ReceptorPhenoCelltype ]
# #CCI[Pair.Name==unique_combs[1]$Pair.Name][LigandPhenoCelltype==unique_combs[1]$LigandPhenoCelltype  & ReceptorPhenoCelltype==unique_combs[1]$ReceptorPhenoCelltype ]
# notresponseRelated2 <- merge(assessment[paramNames%in%c("Non-response_slope","Response_slope")][adjust.p.val>=0.05][order(adjust.p.val)] ,repdd, by=c("Pair.Name","LigandPhenoCelltype","ReceptorPhenoCelltype"),all.x=T)[replic>3][order(adjust.p.val)]
# 
# require(boot)
# ggplot(responseRelated2[paramNames%in%c("Non-response_slope","Response_slope")],aes(Estimate,log(1e-20+adjust.p.val),col=LigandPhenoCelltype))+geom_point()#+facet_grid(dynamic_class3~ReceptorPhenoCelltype)
# ggplot(notresponseRelated2,aes(Estimate,logit(adjust.p.val)))+geom_point(alpha=1,size=0.7)+theme_classic()+#+facet_grid(dynamic_class3~ReceptorPhenoCelltype)
#   geom_point(data=responseRelated2[paramNames%in%c("Non-response_slope","Response_slope")],
#              aes(Estimate,logit(adjust.p.val),col=ReceptorPhenoCelltype),alpha=1,size=0.7)+
#   scale_color_discrete(name="Reciever cell type")+theme(aspect.ratio=1)+
#   scale_y_continuous(breaks=logit(c(0.001,0.01,0.1,0.5,0.9,0.99,0.999)) ,labels=c(0.001,0.01,0.1,0.5,0.9,0.99,0.999))+ labs(y="Adjusted p-value", x="Communication change effect size")#+facet_grid(dynamic_class3~ReceptorPhenoCelltype)
# 
# ggplot(notresponseRelated2,aes(Estimate,logit(adjust.p.val)))+geom_point(alpha=1,size=0.7)+theme_classic()+#+facet_grid(dynamic_class3~ReceptorPhenoCelltype)
#   geom_point(data=responseRelated2[paramNames%in%c("Non-response_slope","Response_slope")],
#              aes(Estimate,logit(adjust.p.val),col=LigandPhenoCelltype),alpha=1,size=0.7)+
#   scale_color_discrete(name="Sender cell type")+theme(aspect.ratio=1)+
#   scale_y_continuous(breaks=logit(c(0.001,0.01,0.1,0.5,0.9,0.99,0.999)) ,labels=c(0.001,0.01,0.1,0.5,0.9,0.99,0.999))+ labs(y="Adjusted p-value", x="Communication change effect size")#+facet_grid(dynamic_class3~ReceptorPhenoCelltype)
# 
# 
# unique( responseRelated2[paramNames%in%c("Non-response_slope","Response_slope")][abs(Estimate)>0.5]%>%dplyr::select(Pair.Name,LigandPhenoCelltype,ReceptorPhenoCelltype) )
# ### Look at trends
unique_combs <- unique( responseRelated2[grep("slope",paramNames)][!grep("Intercept",coefficient)][paramNames=="NtoR_change_slope" ][!grep("anova_Day",Model)] %>% dplyr::select( Pair.Name, LigandPhenoCelltype, ReceptorPhenoCelltype,Treat ) )
 
# join trend to CI data so that positive and negative trends can be extracted
signalDat_pvals_byR <- rbindlist(lapply( 1:nrow(unique_combs), function(j){
   Signal <- CCI[Treat==unique_combs[j]$Treat][Pair.Name==unique_combs[j]$Pair.Name][LigandPhenoCelltype==unique_combs[j]$LigandPhenoCelltype][ReceptorPhenoCelltype==unique_combs[j]$ReceptorPhenoCelltype]
   #slopediffres 
   asess_i<- assessment[Treat==unique_combs[j]$Treat][Pair.Name==unique_combs[j]$Pair.Name][LigandPhenoCelltype==unique_combs[j]$LigandPhenoCelltype][ReceptorPhenoCelltype==unique_combs[j]$ReceptorPhenoCelltype]  
   est<-asess_i$Estimate;names(est)<- asess_i$paramNames              #interceptdiffres <- assessment[Pair.Name==unique_combs[j]$Pair.Name][LigandPhenoCelltype==unique_combs[j]$LigandPhenoCelltype][ReceptorPhenoCelltype==unique_combs[j]$ReceptorPhenoCelltype] [grep("dynamic_class3",coefficient)][!grep(":",coefficient)] 
   Signal[, Rank := j ]
   Signal <-data.table( Signal,t(est))
   return(Signal)
 }))
setnames(signalDat_pvals_byR,old=c("Non-response_int","Non-response_slope","Non-response_int_final","Non-response_int_start"),new=c("Non_response_int","Non_response_slope","Non_response_int_final","Non_response_int_start"))
signalDat_pvals_byR[,trendNR:="positive"]
signalDat_pvals_byR[Non_response_slope<0, trendNR:="negative"]
signalDat_pvals_byR[,trendR:="positive"]
signalDat_pvals_byR[Response_slope<0, trendR:="negative"]

Celltype_Receptor<- "Cancer cells" #"Fibroblasts"#
TRT <- "CombinationRibo"
signalDat_pvals_byR[Treat==TRT][ReceptorPhenoCelltype==Celltype_Receptor][trendR=="positive"]$Pair.Name%>%unique()
# 
# 
# 
for(TRT in unique(signalDat_pvals_byR$Treat)){
for(Celltype_Receptor in unique(signalDat_pvals_byR[Treat==TRT]$ReceptorPhenoCelltype)){
  cat(TRT)
  cat(Celltype_Receptor)
  # find 20 most trending gene pairs
  nplot<-10
  top10increase_byR <- (signalDat_pvals_byR[Treat==TRT][ReceptorPhenoCelltype==Celltype_Receptor][trendR=="positive"]$Pair.Name%>%unique)[1:nplot]
  top10decrease_byR <- (signalDat_pvals_byR[Treat==TRT][ReceptorPhenoCelltype==Celltype_Receptor][trendR=="negative"]$Pair.Name%>%unique)[1:nplot]
  top10increase_byNR <- (signalDat_pvals_byR[Treat==TRT][ReceptorPhenoCelltype==Celltype_Receptor][trendNR=="positive"]$Pair.Name%>%unique)[1:nplot]
  top10decrease_byNR <- (signalDat_pvals_byR[Treat==TRT][ReceptorPhenoCelltype==Celltype_Receptor][trendNR=="negative"]$Pair.Name%>%unique)[1:nplot]
  #top10increase_byR <- (responseRelated2[paramNames%in%c("Non-response_slope","Response_slope","NtoR_change_slope")][ReceptorPhenoCelltype==Celltype_Receptor][Estimate>0]$Pair.Name%>%unique)[1:nplot]
  #top10decrease_byR <- (responseRelated2[paramNames%in%c("Non-response_slope","Response_slope","NtoR_change_slope")][ReceptorPhenoCelltype==Celltype_Receptor][Estimate<0]$Pair.Name%>%unique)[1:nplot]
  #top10increase_byNR <- (responseRelated2[paramNames%in%c("Non-response_slope","Response_slope","NtoR_change_slope")][ReceptorPhenoCelltype==Celltype_Receptor][Estimate>0]$Pair.Name%>%unique)[1:nplot]
  #top10decrease_byNR <- (responseRelated2[paramNames%in%c("Non-response_slope","Response_slope","NtoR_change_slope")][ReceptorPhenoCelltype==Celltype_Receptor][Estimate<0]$Pair.Name%>%unique)[1:nplot]
  #responseRelated2$paramNames%>%unique

  # extract data for visualisation
  require(ggsci)
  plotdd <- signalDat_pvals_byR[Treat==TRT][trendNR=="positive"][Pair.Name%in%c(top10increase_byNR)][ReceptorPhenoCelltype==Celltype_Receptor]
  plotdd$Pair.Name <- factor(plotdd$Pair.Name, levels=top10increase_byNR)
  if(nrow(plotdd)>0){
    pU<-ggplot(plotdd,
             aes(y=log(1+TransductionMu),x=log(1+Day) ,col=dynamic_class3  , group=interaction(Rank,dynamic_class3) ) ) + geom_point()+theme_classic()+
    geom_smooth(method="lm",se=F)+#,alpha=0.3)+
    facet_grid(LigandPhenoCelltype~Pair.Name)+theme(aspect.ratio=1)+
    labs(y=paste0(Celltype_Receptor, " signal receipt"),x="Day" )+scale_color_npg(name="Tumor response")+scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))
  }else{
    pU<- NULL
  }
  plotdd<- signalDat_pvals_byR[Treat==TRT][trendNR=="negative"][Pair.Name%in%c(top10decrease_byNR)][ReceptorPhenoCelltype==Celltype_Receptor]
  plotdd$Pair.Name <- factor(plotdd$Pair.Name, levels=top10decrease_byNR)
  if(nrow(plotdd)>0){
    pD<-ggplot(plotdd,
             aes(y=log(1+TransductionMu),x=log(1+Day) ,col=dynamic_class3  , group=interaction(Rank,dynamic_class3) ) ) + geom_point()+theme_classic()+
    geom_smooth(method="lm",se=F)+
    facet_grid(LigandPhenoCelltype~Pair.Name)+theme(aspect.ratio=1)+
    labs(y=paste0(Celltype_Receptor, " signal receipt"),x="Day" )+scale_color_npg(name="Tumor response")+scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))
  }else{
    pU<- NULL
  }
  if(!is.null(pU)&!is.null(pD)){
  cowplot::plot_grid(pU,pD, labels = c('A) Increasing communication', 'B) Decreasing communication'), label_size = 12,ncol=1)
  }else{
    if(!is.null(pU) & is.null(pD) ){
      pU
    }
    if(!is.null(pU) & is.null(pD) ){
      pD
    }
  }
   ggsave(filename=paste0("/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/LR trends AllArms/Increasing and decreasing communication",Celltype_Receptor,TRT," AllArms.png"),width=15,height=15)
   
}}




unique_combs <- unique( responseRelated2[grep("slope",paramNames)][!grep("Intercept",coefficient)][!grep("anova_Day",Model)] %>% dplyr::select( Pair.Name, LigandPhenoCelltype, ReceptorPhenoCelltype,Treat ) )

# join trend to CI data so that positive and negative trends can be extracted
signalDat_pvals_byR <- rbindlist(lapply( 1:nrow(unique_combs), function(j){
  Signal <- CCI[Treat==unique_combs[j]$Treat][Pair.Name==unique_combs[j]$Pair.Name][LigandPhenoCelltype==unique_combs[j]$LigandPhenoCelltype][ReceptorPhenoCelltype==unique_combs[j]$ReceptorPhenoCelltype]
  #slopediffres 
  asess_i<- assessment[Treat==unique_combs[j]$Treat][Pair.Name==unique_combs[j]$Pair.Name][LigandPhenoCelltype==unique_combs[j]$LigandPhenoCelltype][ReceptorPhenoCelltype==unique_combs[j]$ReceptorPhenoCelltype]  
  est<-asess_i$Estimate;names(est)<- asess_i$paramNames              #interceptdiffres <- assessment[Pair.Name==unique_combs[j]$Pair.Name][LigandPhenoCelltype==unique_combs[j]$LigandPhenoCelltype][ReceptorPhenoCelltype==unique_combs[j]$ReceptorPhenoCelltype] [grep("dynamic_class3",coefficient)][!grep(":",coefficient)] 
  Signal[, Rank := j ]
  Signal <-data.table( Signal,t(est))
  return(Signal)
}))
setnames(signalDat_pvals_byR,old=c("Non-response_int","Non-response_slope","Non-response_int_final","Non-response_int_start"),new=c("Non_response_int","Non_response_slope","Non_response_int_final","Non_response_int_start"))
signalDat_pvals_byR[,trendNR:="positive"]
signalDat_pvals_byR[Non_response_slope<0, trendNR:="negative"]
signalDat_pvals_byR[,trendR:="positive"]
signalDat_pvals_byR[Response_slope<0, trendR:="negative"]

Celltype_Receptor<- "Cancer cells" #"Fibroblasts"#
TRT <- "CombinationRibo"
signalDat_pvals_byR[Treat==TRT][ReceptorPhenoCelltype==Celltype_Receptor][trendR=="positive"]$Pair.Name%>%unique()
# 
# 
# 
for(TRT in unique(signalDat_pvals_byR$Treat)){
  for(Celltype_Receptor in unique(signalDat_pvals_byR[Treat==TRT]$ReceptorPhenoCelltype)){
    cat(TRT)
    cat(Celltype_Receptor)
    # find 20 most trending gene pairs
    nplot<-20
    top10increase_byR <- (signalDat_pvals_byR[Treat==TRT][ReceptorPhenoCelltype==Celltype_Receptor][trendR=="positive"]$Pair.Name%>%unique)[1:nplot]
    top10decrease_byR <- (signalDat_pvals_byR[Treat==TRT][ReceptorPhenoCelltype==Celltype_Receptor][trendR=="negative"]$Pair.Name%>%unique)[1:nplot]
    top10increase_byNR <- (signalDat_pvals_byR[Treat==TRT][ReceptorPhenoCelltype==Celltype_Receptor][trendNR=="positive"]$Pair.Name%>%unique)[1:nplot]
    top10decrease_byNR <- (signalDat_pvals_byR[Treat==TRT][ReceptorPhenoCelltype==Celltype_Receptor][trendNR=="negative"]$Pair.Name%>%unique)[1:nplot]
    #top10increase_byR <- (responseRelated2[paramNames%in%c("Non-response_slope","Response_slope","NtoR_change_slope")][ReceptorPhenoCelltype==Celltype_Receptor][Estimate>0]$Pair.Name%>%unique)[1:nplot]
    #top10decrease_byR <- (responseRelated2[paramNames%in%c("Non-response_slope","Response_slope","NtoR_change_slope")][ReceptorPhenoCelltype==Celltype_Receptor][Estimate<0]$Pair.Name%>%unique)[1:nplot]
    #top10increase_byNR <- (responseRelated2[paramNames%in%c("Non-response_slope","Response_slope","NtoR_change_slope")][ReceptorPhenoCelltype==Celltype_Receptor][Estimate>0]$Pair.Name%>%unique)[1:nplot]
    #top10decrease_byNR <- (responseRelated2[paramNames%in%c("Non-response_slope","Response_slope","NtoR_change_slope")][ReceptorPhenoCelltype==Celltype_Receptor][Estimate<0]$Pair.Name%>%unique)[1:nplot]
    #responseRelated2$paramNames%>%unique
    
    # extract data for visualisation
    require(ggsci)
    plotdd <- signalDat_pvals_byR[Treat==TRT][trendNR=="positive"][Pair.Name%in%c(top10increase_byNR)][ReceptorPhenoCelltype==Celltype_Receptor]
    plotdd$Pair.Name <- factor(plotdd$Pair.Name, levels=top10increase_byNR)
    if(nrow(plotdd)>0){
      pU<-ggplot(plotdd,
                 aes(y=log(1+TransductionMu),x=log(1+Day) ,col=dynamic_class3  , group=interaction(Rank,dynamic_class3) ) ) + geom_point()+theme_classic()+
        geom_smooth(method="lm",se=F)+#,alpha=0.3)+
        facet_grid(LigandPhenoCelltype~Pair.Name)+theme(aspect.ratio=1)+
        labs(y=paste0(Celltype_Receptor, " signal receipt"),x="Day" )+scale_color_npg(name="Tumor response")+scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))
    }else{
      pU<- NULL
    }
    plotdd<- signalDat_pvals_byR[Treat==TRT][trendNR=="negative"][Pair.Name%in%c(top10decrease_byNR)][ReceptorPhenoCelltype==Celltype_Receptor]
    plotdd$Pair.Name <- factor(plotdd$Pair.Name, levels=top10decrease_byNR)
    if(nrow(plotdd)>0){
      pD<-ggplot(plotdd,
                 aes(y=log(1+TransductionMu),x=log(1+Day) ,col=dynamic_class3  , group=interaction(Rank,dynamic_class3) ) ) + geom_point()+theme_classic()+
        geom_smooth(method="lm",se=F)+
        facet_grid(LigandPhenoCelltype~Pair.Name)+theme(aspect.ratio=1)+
        labs(y=paste0(Celltype_Receptor, " signal receipt"),x="Day" )+scale_color_npg(name="Tumor response")+scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))
    }else{
      pU<- NULL
    }
    if(!is.null(pU)&!is.null(pD)){
      cowplot::plot_grid(pU,pD, labels = c('A) Increasing communication', 'B) Decreasing communication'), label_size = 12,ncol=1)
    }else{
      if(!is.null(pU) & is.null(pD) ){
        pU
      }
      if(!is.null(pU) & is.null(pD) ){
        pD
      }
    }
    ggsave(filename=paste0("/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/LR trends AllArms/Overall Increasing and decreasing communication",Celltype_Receptor,TRT," AllArms.png"),width=35,height=15)
    
  }}


# ggplot(signalDat_pvals_byR[Pair.Name%in%c(top10decrease_byNR[1:10])][ReceptorPhenoCelltype=="Fibroblasts"],
#        aes(y=log(1+TransductionMu),x=as.factor(Day) ,col=dynamic_class3  , group=interaction(Rank,dynamic_class3) ) ) + geom_point()+
#   geom_smooth(method="lm",se=F)+
#   facet_grid(LigandPhenoCelltype~Pair.Name)+theme(aspect.ratio=1)
# 
# 
# ggplot(signalDat_pvals_byR[Pair.Name%in%c(top10increase_byNR[1:10])][ReceptorPhenoCelltype=="Fibroblasts"],
#        aes(y=log(1+TransductionMu),x=as.factor(Day) ,col=dynamic_class3  , group=interaction(Rank,dynamic_class3) ) ) + geom_point()+
#   geom_smooth(method="lm",se=F)+
#   facet_grid(LigandPhenoCelltype~Pair.Name)+theme(aspect.ratio=1)
# 
# 
# ggplot(signalDat_pvals_byR[Pair.Name%in%c(top10increase_byR[1:3])],
#        aes(y=log(1+TransductionMu),x=as.factor(Day) ,col=LigandPhenoCelltype ,group=Rank) ) + geom_point()+
#   geom_smooth(method="lm",se=F)+
#   facet_grid(dynamic_class3~ReceptorPhenoCelltype)
# 
# ### Plot those increasing LR pairs
# ggplot(signalDat_pvals_byR[ReceptorPhenoCelltype=="Cancer cells"][Pair.Name%in% top10increase_byR],
#        aes(y=log(1+TransductionMu), x=log(1+Day) , col=LigandPhenoCelltype , group=Rank) ) + geom_jitter(size=2,height =0,width=0.2)+theme_classic(base_size=14)+
#   geom_smooth(method="lm",se=F)+ theme(aspect.ratio=1)+scale_color_discrete(name="Signalling cell type")+
#   facet_wrap(~Pair.Name,nrow=4)+labs(y="Communication with cancer cells", x="Day")+scale_x_continuous(breaks=log(1+c(0,14,180) ) , labels=c(0,14,180)) +scale_y_continuous(breaks=log(c(1,10,100,1000,1000,10000)), labels=c(1,10,100,1000,1000,10000) )
# 
# ### Plot those decreasing LR pairs
# ggplot(signalDat_pvals_byR[ReceptorPhenoCelltype=="Cancer cells"][Pair.Name%in% top10decrease_byR],
#        aes(y=log(1+TransductionMu),x=log(1+Day) ,col=LigandPhenoCelltype ,group=Rank) ) + geom_jitter(size=2,height =0,width=0.2)+theme_classic(base_size=14)+
#   geom_smooth(method="lm",se=F)+ theme(aspect.ratio=1)+scale_color_discrete(name="Signalling cell type")+
#   facet_wrap(~Pair.Name,nrow=4)+labs(y="Communication with cancer cells", x="Day")+scale_x_continuous(breaks=log(1+c(0,14,180) ) , labels=c(0,14,180)) +scale_y_continuous(breaks=log(c(1,10,100,1000,1000,10000)), labels=c(1,10,100,1000,1000,10000) )
# 
# 
# save(signalDat_pvals,assessmentSimple,file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/When cells communicate_communication changes/overalltrends.RData")
# 
# 
# 
# 
# 
# 
# 
# 
# responseRelated2$LigandPhenoCelltype%>%table()
# 
plt<-responseRelated2[order(adjust.p.val)]
# plt[1:30]
# j <- 13
# 
# 
# 
# 
# ggplot(CCI[Pair.Name%in%plt[j]$Pair.Name][LigandPhenoCelltype==plt[j]$LigandPhenoCelltype][ReceptorPhenoCelltype==plt[j]$ReceptorPhenoCelltype],
#        aes(y=log(1+TransductionMu),x=log(1+Day) ,col=dynamic_class3 ) ) + geom_point()+geom_smooth(method="lm",se=F)
# 
# image(responseRelatedSimple2[order(adjust.p.val)]%>%dplyr::select(LigandPhenoCelltype,ReceptorPhenoCelltype)%>%table())
# 
# unique( signalDat_pvals%>%dplyr::select(Pair.Name,trend,LigandPhenoCelltype,ReceptorPhenoCelltype)) %>%dplyr::select(LigandPhenoCelltype,ReceptorPhenoCelltype)%>%table()
# 
# ggplot( as.data.table(responseRelatedSimple2[order(adjust.p.val)]%>%dplyr::select(LigandPhenoCelltype,ReceptorPhenoCelltype)%>%table()),
#         aes(y=LigandPhenoCelltype,x=  ReceptorPhenoCelltype ,fill= N))+geom_tile()
# 
# 
# 
# grep("dynamic_class3Response",coefficient)
centVal<-mean(log(CCI$TransductionMu))
scalsd<-sd(log(CCI$TransductionMu))

load("/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/GeneOntology/growthFactorReceptors.RData")
growthFactorReceptors <- c( #"EGFR" ,"ERBB2","ERBB3","ERBB4",
                            #"FGFR1","FGFR2","FGFR3","FGFR4","FGFRL1",#"IGF1R","IGF2R","INSR",
                            #"NRP1","NRP2",
                            "TGFBR1","TGFBR2" ,"TGFBR3"
                            ) 

RRR<-"Response"
RRR<-"Non-response"
Day_var=180
TTT<- "CombinationRibo"
# #
CCI_plot <- CCI[dynamic_class3==RRR][Day==Day_var][Pair.Name%in%
   LRpairsFiltered[ HPMR.Receptor %in% growthFactorReceptors]$Pair.Name
   ][
     LigandPhenoCelltype%in%c("Cancer cells","Fibroblasts")
   ][
     ReceptorPhenoCelltype%in%c("Cancer cells","Fibroblasts")
     ][order(LigandPhenoCelltype,ReceptorPhenoCelltype)][Treat==TTT]

CCI_plot[,signalCol:="blue"]
CCI_plot[LigandPhenoCelltype=="Cancer cells",signalCol:="red"]

gfCCI <- CCI[Pair.Name%in%
      LRpairsFiltered[ HPMR.Receptor %in% growthFactorReceptors]$Pair.Name
    ][Treat==TTT][LigandPhenoCelltype== "Cancer cells"& ReceptorPhenoCelltype== "Fibroblasts"]
gfCCIsummar <- data.table( gfCCI %>% group_by(Pair.Name, Day, dynamic_class3) %>% summarise(scaleTransduction=mean(log(1+scaleTransduction)) ))
gfCCIsummarB <- gfCCIsummar %>% spread(dynamic_class3, scaleTransduction)
setnames( gfCCIsummarB , old=c("Non-response", "Response"), new=c("NonResponse", "Response"))
gfCCIsummarBD0 <- na.omit(gfCCIsummarB[Day==0])
gfCCIsummarBD0[,NminusR:= NonResponse - Response]
gfCCIsummarB[,NminusR:= NonResponse - Response]
gfCCIsummarBD0$Pair.Name <- factor( gfCCIsummarBD0$Pair.Name, levels= gfCCIsummarBD0[order(NminusR)]$Pair.Name )

ggplot( gfCCI[Pair.Name=="TGFB3_TGFBR3"] , aes(y=log(1+scaleTransduction), x=Day, col=dynamic_class3,group=interaction(dynamic_class3,Pair.Name)))+geom_point()+geom_smooth(method="lm",se=F)
ggplot( gfCCI[Pair.Name=="TGFB1_TGFBR3"] , aes(y=log(1+scaleTransduction), x=Day, col=dynamic_class3,group=interaction(dynamic_class3,Pair.Name)))+geom_point()+geom_smooth(method="lm",se=F)

ddplot<-gfCCI[Pair.Name %in% as.character(gfCCIsummarBD0$Pair.Name)][
  Pair.Name %in%  assessment[pval<0.05][ coefficient!="(Intercept)"][LigandPhenoCelltype=="Cancer cells" & ReceptorPhenoCelltype=="Fibroblasts"][Treat=="CombinationRibo"]$Pair.Name]

ddplot2<-data.table( ddplot %>% group_by(Pair.Name) %>% mutate(plotscaleTransduction=scale(log(1+scaleTransduction))) )
ddplot2[,lnDay2:=log(1+Day)-0.1]
#ddplot2[dynamic_class3=="Response",lnDay2:=log(1+Day)+0.1]
ddplot2[Pair.Name=="GDF9_TGFBR1",lnDay2:=log(1+Day)+0.1]
ddplot2[Pair.Name=="INHBA_TGFBR3",lnDay2:=log(1+Day)+0.2]
ddplot2[Pair.Name=="TGFB1_TGFBR2",lnDay2:=log(1+Day)]
ddplot2[Pair.Name=="TGFB1_TGFBR3",lnDay2:=log(1+Day)-0.2]
ddplot2[Pair.Name=="TGFB3_TGFBR2",lnDay2:=log(1+Day)+0.3]




GFsignals <- ggplot( ddplot2 ,
                 aes(y= plotscaleTransduction, x= lnDay2,  fill= Pair.Name ,
                     col= Pair.Name, group= interaction(dynamic_class3, Pair.Name))) +
  geom_smooth(method="lm", se=F, alpha=0.3, size=2) + 
  geom_point(size=3.5)+#, position = position_dodge(width=0.05) ) +
  theme_classic() + facet_wrap( ~ dynamic_class3) + 
  geom_violin(aes(group= interaction(Pair.Name, Day) ), alpha=0.5) +
  scale_color_viridis_d(name= "Communication pathway" ) +
  scale_fill_viridis_d(name= "Communication pathway" ) +
  labs(x= "Day", y= "Cancer signaling to fibroblasts" )+theme(aspect.ratio=1)
ggsave(GFsignals,
       filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cross cell type communication/GF signals cancer to fibroblasts by pathway and response.png")

BlankGFsignals <- ggplot( ddplot2 ,
                     aes(y= plotscaleTransduction, x= lnDay2,  fill= Pair.Name ,
                         col= Pair.Name, group= interaction(dynamic_class3, Pair.Name))) +
  geom_smooth(method="lm", se=F, alpha=0.3, size=2) + 
  geom_point(size=4.5)+#, position = position_dodge(width=0.05) ) +
  theme_classic() + facet_wrap( ~ dynamic_class3) + 
  geom_violin(aes(group= interaction(Pair.Name, Day) ), alpha=0.5) +
  scale_color_viridis_d(name= "Communication pathway" ) +
  scale_fill_viridis_d(name= "Communication pathway" ) +
  labs(x= "Day", y= "Cancer to fibroblast communication" )+
  theme(axis.text=element_blank(),
      axis.title=element_blank(),
      aspect.ratio=1,
      legend.position="none",
      strip.background = element_blank(),
      strip.text=element_blank())
ggsave(BlankGFsignals,
       filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cross cell type communication/BLANK GF signals cancer to fibroblasts by pathway and response.png")



GFsignals <- ggplot( ddplot2[Day==0] ,
                     aes(y= plotscaleTransduction, x= Pair.Name,  fill= Pair.Name ,
                         col= Pair.Name, group= interaction(dynamic_class3, Pair.Name))) +
  geom_point(size=3.5)+#, position = position_dodge(width=0.05) ) +
  theme_classic() + facet_wrap( ~ dynamic_class3) + 
  geom_violin(aes(group= interaction(Pair.Name, Day) ), alpha=0.5) +
  scale_color_viridis_d(name= "Communication pathway" ) +
  scale_fill_viridis_d(name= "Communication pathway" ) +
  labs(x= "Communication pathway", y= "Cancer signaling to fibroblasts" )+theme(aspect.ratio=1, axis.text.x=element_text(angle=90) )
ggsave(GFsignals,
       filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cross cell type communication/GF signals day0 cancer to fibroblasts by pathway and response.png")

BlankGFsignals <- ggplot( ddplot2[Day==0] ,
                          aes(y= plotscaleTransduction, x= Pair.Name,  fill= Pair.Name ,
                              col= Pair.Name, group= interaction(dynamic_class3, Pair.Name))) +
  geom_point(size=4.5)+#, position = position_dodge(width=0.05) ) +
  theme_classic() + facet_wrap( ~ dynamic_class3) + 
  geom_violin(aes(group= interaction(Pair.Name, Day) ), alpha=0.5) +
  scale_color_viridis_d(name= "Communication pathway" ) +
  scale_fill_viridis_d(name= "Communication pathway" ) +
  labs(x= "Communication pathway", y= "Cancer signaling to fibroblasts" )+theme(aspect.ratio=1 )+
  theme(axis.text=element_blank(),
        axis.title=element_blank(),
        aspect.ratio=1,
        legend.position="none",
        strip.background = element_blank(),
        strip.text=element_blank())
ggsave(BlankGFsignals,
       filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cross cell type communication/BLANK GF signals day0 cancer to fibroblasts by pathway and response.png")


pheatdd<- ( ddplot2[Day==0]%>%
     dplyr::select(-c(scaleTransduction ,lnDay2,TransductionMu,SumSignal, Receptor,Ligandtot, Ligand_Ntot))%>%
     spread(Pair.Name,plotscaleTransduction,fill= 0)%>%
     dplyr::select(Patient.Study.ID,GDF9_TGFBR1:TGFB3_TGFBR3))

pheatdd1 <- t(as.matrix(pheatdd[,-1]));colnames(pheatdd1)<-pheatdd$Patient.Study.ID
anotdd <- unique( ddplot2[Day==0]%>%
                    dplyr::select(-c(scaleTransduction ,lnDay2,TransductionMu,SumSignal, Receptor,Ligandtot, Ligand_Ntot))%>%
                    spread(Pair.Name,plotscaleTransduction,fill= 0) %>%
                    dplyr::select(Patient.Study.ID,dynamic_class3 ))
rownames(anotdd) <- anotdd$Patient.Study.ID;anotdd$Patient.Study.ID <- NULL
colannotdd <- anotdd
colannotdd[dynamic_class3=="Non-response",dynamic_class3:="red"]
colannotdd[dynamic_class3=="response",dynamic_class3:="blue"]

pheatmap( pheatdd1,
          annotation_col=anotdd ,
          annotation_colors=list(colannotdd$dynamic_class3 ) )


Blank <- ggplot( ddplot2 ,
        aes(y=plotscaleTransduction, x=lnDay2,  fill=Pair.Name ,
            col=Pair.Name, group=interaction(dynamic_class3, Pair.Name)))+
  geom_smooth(method="lm", se=F, alpha=0.3, size=2)+ 
  geom_point(size=2)+#, position = position_dodge(width=0.05) ) +
  theme_classic() + facet_wrap( ~ dynamic_class3)+ 
  geom_violin(aes(group=interaction(Pair.Name, Day) ), alpha=0.5)




ggplot( ddplot2 ,
        aes(y=plotscaleTransduction, x=lnDay2,  fill=Pair.Name ,
            col=dynamic_class3,group=interaction(dynamic_class3,Pair.Name)))+
  geom_smooth(method="lm",se=F,alpha=0.3,size=2)+ #,aes(linetype=Pair.Name)
  geom_point(size=2,position = position_dodge(width=0.05)) +
  theme_classic()+facet_wrap(~dynamic_class3)+ 
  geom_violin(aes(group=interaction(dynamic_class3,Day) ),alpha=0.5)
  


ggplot( gfCCI[Pair.Name=="TGFB2_TGFBR2"] , aes(y=log(1+scaleTransduction), x=Day, col=dynamic_class3,group=interaction(dynamic_class3,Pair.Name)))+geom_point()+geom_smooth(method="lm",se=F)
gfCCI$Pair.Name <- factor( gfCCI$Pair.Name, levels= gfCCIsummarBD0[order(NminusR)]$Pair.Name )

gfCCIwide<-spread( gfCCI[Day==0]%>%dplyr::select(-c(TransductionMu,SumSignal, Receptor,Ligandtot,Ligand_Ntot)) ,Pair.Name,scaleTransduction, fill=0)

u1<-umap::umap(gfCCIwide%>%select(GNB3_TGFBR1:GDF9_TGFBR1))
gfCCIwide1 <- data.table(gfCCIwide , u1$layout)

ggplot(gfCCIwide1[],#[Day==180],
       aes(y=V1, x=V2, col=dynamic_class3))+geom_point()

ggplot(gfCCI[Pair.Name%in% as.character(gfCCIsummarBD0$Pair.Name)],#[Day==180],
       aes(y=log(1+scaleTransduction), x=Pair.Name, col=dynamic_class3))+
  #geom_point() +geom_violin()+
  geom_point(size=4,position = position_dodge(width=0.9)) + geom_violin(alpha=0.5, position=position_dodge())+
  theme(axis.text.x=element_text(angle=90))# + facet_wrap(~Day) 

ggplot( gfCCI[Pair.Name=="TGFB2_TGFBR2"] , aes(y=log(1+scaleTransduction), x=Day, col=dynamic_class3,group=interaction(dynamic_class3,Pair.Name)))+geom_point()+geom_smooth(method="lm",se=F)

ggplot(gfCCIsummarBD0, aes(y=NonResponse-Response, x=Response))+
  geom_point()

ggplot(gfCCIsummarBD0, aes(y=NminusR, x=Pair.Name))+
  geom_point() +
  theme(axis.text.x=element_text(angle=90))


ggplot(CCI[Pair.Name%in%
             LRpairsFiltered[ HPMR.Receptor %in% growthFactorReceptors]$Pair.Name
           ][Treat==TTT][LigandPhenoCelltype== "Cancer cells"& ReceptorPhenoCelltype== "Fibroblasts"],
       aes(x=Pair.Name, y=log(1+scaleTransduction), col=dynamic_class3))+
  geom_point() + facet_wrap(~dynamic_class3)

g <- graph.data.frame(CCI_plot%>%dplyr::select(-Patient.Study.ID), directed=TRUE)
E(g)$weight <- exp((log(CCI_plot$TransductionMu)-centVal)/scalsd)
is.weighted((g))
# 
# # specifcy color of nodes
resp_cols<- ggsci::pal_npg("nrc")(2)
if(RRR=="Non-response"){V(g)$color <- resp_cols[1] }else{if(RRR=="Response"){V(g)$color <- resp_cols[2]}}
# 
# 
# # circle layout.
n= length(unique(CCI_plot$ReceptorPhenoCelltype)) -1
pts.circle <- t(sapply(1:n,function(r)c(cos(2*r*pi/n),sin(2*r*pi/n))))
# 
# #NodeList <- data.table(sort( unique(SignalReceiverRes$ReceptorPhenoCelltype) ), c(4,9.5,5,8,9,4,2,9,1,1.5 ) /5-1 ,c(9,5,5,8,6,1,7,2,4,1.5)  /5-1 )
NodeList <- data.table(c("Cancer cells", "CD4+ T cells" ,"Adipocytes","Fibroblasts","Normal epithelial cells","Pericytes"  , "Endothelial cells" ,"Macrophages", "B cells","CD8+ T cells" )  , c(0, pts.circle[,1] )  ,  c(0, pts.circle[,2] ) )
# 
presloc<-NodeList[na.omit((match(names(V(g) ), NodeList$V1  ))), ]
# 
# #present_i <- unique( as.vector(unlist(as.matrix(  unique(plotddi[order(from,to)]%>%dplyr::select(from,to)) ) )) )
# # plot.igraph(g,rescale=FALSE,#vertex.label=V(g)$name,
# #             layout = as.matrix(presloc%>%select(-V1))   , xlim =c(-1,1), ylim =c(-1,1),#layout=layout.fruchterman.reingold, 
# #             vertices=NodeList[V1%in%presloc$V1],
# #             edge.label.color = adjustcolor("black", .5),
# #             edge.color="black",edge.width=0.03*E(g)$weight)
# # 
# # plot.igraph(g,rescale=FALSE,#vertex.label=V(g)$name,
# #             layout = as.matrix(presloc%>%select(-V1))   , xlim =c(-1,1), ylim =c(-1,1),#layout=layout.fruchterman.reingold, 
# #             vertices=NodeList[V1%in%presloc$V1],
# #             edge.label.color =adjustcolor("black", .5),
# #             edge.color=adjustcolor("black", .05),edge.width=0.01*E(g)$weight)
# # 
# 
# 
# 
# 
plot.igraph(g,rescale=FALSE,#vertex.label=V(g)$name,
            layout = as.matrix(presloc%>%select(-V1))   , xlim =c(-1,1), ylim =c(-1,1),#layout=layout.fruchterman.reingold, 
            vertices=NodeList[V1%in%presloc$V1],        
            edge.label.color =adjustcolor("black", .5),
            edge.color= adjustcolor(CCI_plot$signalCol, .5),
           # edge.color=adjustcolor("black", .05),
            edge.width=0.1*E(g)$weight)
 
# 
# 
# 
# # plot.igraph(g,rescale=FALSE,#vertex.label=V(g)$name,
# #             layout = as.matrix(presloc%>%select(-V1))   , xlim =c(-1,1), ylim =c(-1,1),#layout=layout.fruchterman.reingold, 
# #             vertices=NodeList[V1%in%presloc$V1],
# #             edge.label.color =adjustcolor("black", .5),
# #             edge.color=adjustcolor("black", .05),edge.width=0.2)
# 
# 
# 
# 
# ### per indiv plots
# 
# perIndiv=TRUE
# 
# CCI <- rbindlist(lapply(1:length(filenamesCCI), function(ii){
#   load(file= paste0(savelocCCI, filenamesCCI[ii]))
#   cat(ii);#abc<-(Communication[Patient.Study.ID=="001-101" &Day==0][Pair.Name=="EGF_EGFR"]%>%dplyr::select(key_,key_signaller, Signal)%>%spread(key_signaller, Signal))#sum( abc[1,-1],na.rm=T)
#   SumComm <- data.table( Communication %>%
#                            group_by(Patient.Study.ID,LigandPhenoCelltype,ReceptorPhenoCelltype,    Day ,Pair.Name,dynamic_class3,ARM)%>%
#                            dplyr::summarise(TransductionMu=mean(TransductionperSignaller),  Receptor=mean(Receptor),Ligandtot=sum(Ligand),Ligand_Ntot=sum(Ligand_N))  )[!( is.na(LigandPhenoCelltype)|is.na(ReceptorPhenoCelltype) ) ]
# }))
# 
# #save(CCI,allphenotypes, uu,perIndiv,file ="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Communication output merged/IndividualCommunicationMerged.RData" )
# 
# centVal<-mean(log(CCI$TransductionMu))
# scalsd<-sd(log(CCI$TransductionMu))
# RRR<-"Response"
# RRR<-"Non-response"
# Day_var=14
# #CCI_plot <- CCI[dynamic_class3==RRR][Day==Day_var][Pair.Name%in%summarytable$key_][order(LigandPhenoCelltype,ReceptorPhenoCelltype)]
# CCI_plot <- CCI[dynamic_class3==RRR][Day==Day_var][][order(LigandPhenoCelltype,ReceptorPhenoCelltype)]
# g <- graph.data.frame(CCI_plot%>%dplyr::select(-Patient.Study.ID), directed=TRUE)
# E(g)$weight <- exp((log(CCI_plot$TransductionMu)-centVal)/scalsd)
# is.weighted((g))
# 
# # specifcy color of nodes
# resp_cols<- ggsci::pal_npg("nrc")(2)
# if(RRR=="Non-response"){V(g)$color <- resp_cols[1] }else{if(RRR=="Response"){V(g)$color <- resp_cols[2]}}
# 
# 
# # circle layout.
# n= length(unique(CCI_plot$ReceptorPhenoCelltype)) -1
# pts.circle <- t(sapply(1:n,function(r)c(cos(2*r*pi/n),sin(2*r*pi/n))))
# 
# #NodeList <- data.table(sort( unique(SignalReceiverRes$ReceptorPhenoCelltype) ), c(4,9.5,5,8,9,4,2,9,1,1.5 ) /5-1 ,c(9,5,5,8,6,1,7,2,4,1.5)  /5-1 )
# NodeList <- data.table(c("Cancer cells", "CD4+ T cells" ,"Adipocytes","Fibroblasts","Normal epithelial cells","Pericytes"  , "Endothelial cells" ,"Macrophages", "B cells","CD8+ T cells" )  , c(0, pts.circle[,1] )  ,  c(0, pts.circle[,2] ) )
# 
# presloc<-NodeList[na.omit((match(names(V(g) ), NodeList$V1  ))), ]
# 
# #present_i <- unique( as.vector(unlist(as.matrix(  unique(plotddi[order(from,to)]%>%dplyr::select(from,to)) ) )) )
# # plot.igraph(g,rescale=FALSE,#vertex.label=V(g)$name,
# #             layout = as.matrix(presloc%>%select(-V1))   , xlim =c(-1,1), ylim =c(-1,1),#layout=layout.fruchterman.reingold, 
# #             vertices=NodeList[V1%in%presloc$V1],
# #             edge.label.color = adjustcolor("black", .5),
# #             edge.color="black",edge.width=0.03*E(g)$weight)
# # 
# # plot.igraph(g,rescale=FALSE,#vertex.label=V(g)$name,
# #             layout = as.matrix(presloc%>%select(-V1))   , xlim =c(-1,1), ylim =c(-1,1),#layout=layout.fruchterman.reingold, 
# #             vertices=NodeList[V1%in%presloc$V1],
# #             edge.label.color =adjustcolor("black", .5),
# #             edge.color=adjustcolor("black", .05),edge.width=0.01*E(g)$weight)
# # 
# 
# 
# 
# 
# plot.igraph(g,rescale=FALSE,#vertex.label=V(g)$name,
#             layout = as.matrix(presloc%>%select(-V1))   , xlim =c(-1,1), ylim =c(-1,1),#layout=layout.fruchterman.reingold, 
#             vertices=NodeList[V1%in%presloc$V1],
#             edge.label.color =adjustcolor("black", .5),
#             edge.color=adjustcolor("black", .005),edge.width=0.000000001*E(g)$weight)
# 
# 
