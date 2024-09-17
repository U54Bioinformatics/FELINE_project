rm(list=ls())   
require(rTensor); require(nnTensor); library(abind)
require(data.table); require(dplyr); require(ggplot2); require(tidyr)
require(igraph)
require(ggraph)
require(tidygraph)

### Load Ligand Receptor database list of Ramilowski et al 2015
LRpairsFiltered <- data.table(read.csv( file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/SI data/FigureS4/SourceData_Figure S4 LigandReceptorPairsRamilowski2015.csv"))
LRgenelist <- unique( c(LRpairsFiltered$HPMR.Receptor, LRpairsFiltered$HPMR.Ligand) )

### Load GO database list of GF receptors
#growthFactorReceptors2 <- read.csv(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/SI data/FigureS4/Figure S4 GeneOntologyGrowthFactorReceptors.csv")$x

### Load some signalling data for downstream analysis
savloc<-"/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/Figure3/"
CCI_all <-  read.csv(file=paste0(savloc, "SourceData_Figure3_CellCummunicationIndexTumorWideSingnaling.csv"))
CCI_discovery <- data.table(CCI_all)[Cohort=="Discovery"]
CCI<-CCI_discovery

# look up of the LR pairs and signal-receivers pairs to consider 
lu_lm <- unique(CCI%>%dplyr::select(Pair.Name,LigandPhenoCelltype,ReceptorPhenoCelltype,Treat))
CCI$Pair.Name%>%unique%>%length


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
responseRelated2_0 <- merge(responseRelated[order(adjust.p.val)] ,repdd, by=c("Pair.Name","LigandPhenoCelltype","ReceptorPhenoCelltype","Treat"),all.x=T)[replic>3][order(adjust.p.val)]
responseRelated2_0[Model=="one way anova_DayStart_Response"][Treat=="CombinationRibo"]$LigandPhenoCelltype%>%table()
responseRelated2_0[Model=="one way anova_DayStart_Response"][Treat=="CombinationRibo"]$ReceptorPhenoCelltype%>%table()

# save(responseRelated2_0,responseRelated,assessment,file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/When cells communicate_communication changes AllArms/trendsby response AllArms.RData")
# load(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/When cells communicate_communication changes AllArms/trendsby response AllArms.RData")

# extract intitial differences in ribo arm
initdiffs <- unique(responseRelated2_0[Model=="one way anova_DayStart_Response"][][] %>% dplyr::select(-dynamic_class3))

initdiffs$LigandPhenoCelltype <-factor(initdiffs$LigandPhenoCelltype, 
                                       levels= data.table(initdiffsCountLong%>%group_by(LigandPhenoCelltype)%>%dplyr::summarise(n=sum(value)))[order(-n)]$LigandPhenoCelltype )
initdiffs$ReceptorPhenoCelltype <-factor(initdiffs$ReceptorPhenoCelltype, 
                                         levels= data.table(initdiffsCountLong%>%group_by(ReceptorPhenoCelltype)%>%dplyr::summarise(n=sum(value)))[order(-n)]$ReceptorPhenoCelltype )

ggplot(initdiffs, aes(fill = LigandPhenoCelltype, x = ReceptorPhenoCelltype)) + 
  geom_bar(aes()) + theme_classic()+
  scale_fill_viridis_d("Signal sending cell type", option="B",direction = -1) +
  labs(x="Signal receiving cell type", y="Number of communication pathways \n initially differing in resistant and sensitive tumors") +
  theme(axis.text.x=element_text(size=9, angle=90, vjust=0.3),
        axis.text.y=element_text(size=9),
        aspect.ratio=1)



initdiffs[,labelLigandPhenoCelltype:=LigandPhenoCelltype]
initdiffs[LigandPhenoCelltype%in% c("Cancer cells","Normal epithelial cells"),labelLigandPhenoCelltype:="Cancer/Epithelial cells"]
initdiffs[LigandPhenoCelltype%in% c("CD4+ T cells","CD8+ T cells","B cells"),labelLigandPhenoCelltype:="Lymphocytes"]
initdiffs[LigandPhenoCelltype%in% c("Macrophages"),labelLigandPhenoCelltype:="Myeloid cells"]
initdiffs[ReceptorPhenoCelltype%in% c("Macrophages"),ReceptorPhenoCelltype:="Myeloid cells"]
initdiffs[ReceptorPhenoCelltype%in% c("Normal epithelial cells"),ReceptorPhenoCelltype:="Diploid epithelial cells"]
initdiffs$ReceptorPhenoCelltype <- factor(initdiffs$ReceptorPhenoCelltype, 
                                          levels= c("Myeloid cells", "Diploid epithelial cells", "Adipocytes","Pericytes", "Fibroblasts", "Cancer cells", "Endothelial cells", "CD4+ T cells", "CD8+ T cells" ,"B cells") )
initdiffs$labelLigandPhenoCelltype <- factor(initdiffs$labelLigandPhenoCelltype, 
                                             levels= data.table(initdiffs%>%group_by(labelLigandPhenoCelltype)%>%dplyr::summarise(n=n()))[order(-n)]$labelLigandPhenoCelltype )

initdiffs[,Treatlab:="Combination ribociclib"]
initdiffs[Treat == "LetrozoleAlone",Treatlab:="Letrozole alone"]

ggplot(initdiffs, aes(fill = labelLigandPhenoCelltype, x = ReceptorPhenoCelltype)) + 
  geom_bar(aes()) + theme_classic(base_size=18)+
  scale_fill_viridis_d("Signal sending cell type", option="B",direction = -1) +
  labs(x="Signal receiving cell type", y="Number of communication pathways \n initially differing in resistant and sensitive tumors") +
  theme(axis.text.x=element_text(, angle=90, vjust=0.3),
        axis.text.y=element_text(),
        aspect.ratio=1) +facet_wrap(~Treatlab)


savloc <- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/SI data/FigureS7/"
write.csv(initdiffs, file=paste0(savloc,"SourceData_NumberCellTypeCommunicationsInitiallyDifferingRvsS.csv" ))



#####
# Figure S7 heatmap of cell type communication
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

TemporalCellTypeCommunication<-average_muln_scaleTransduction2[][Treat=="CombinationRibo"][LigandPhenoCelltype!="B cells"][ReceptorPhenoCelltype!="B cells"] 
ggplot( TemporalCellTypeCommunication,aes(LigandPhenoCelltype,ReceptorPhenoCelltype,fill=lnfoldchange  )) + geom_tile()+
  facet_wrap(dynamic_class3~Day)+scale_fill_viridis_c(name="Log fold change \n in communication \n relative to baseline",option="B") + theme_classic()+
  labs(y="Signal receiver",x="Signal sender")+theme(aspect.ratio=1,axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


TemporalCellTypeCommunication

savloc <- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/SI data/FigureS7/"
write.csv(TemporalCellTypeCommunication, file=paste0(savloc,"SourceData_TemporalCellTypeCommunicationR vsSTumors.csv" ))
