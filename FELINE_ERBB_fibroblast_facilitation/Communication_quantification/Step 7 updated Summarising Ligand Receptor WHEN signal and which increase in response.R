rm(list=ls())   
require(data.table); require(dplyr); require(ggplot2); require(tidyr);require(igraph)

### Phenotype data
load(file="~/Dropbox/Cancer_pheno_evo/data/FELINE2/PhenotypesOfAllCellTypesAllArms/PhenotypesOfAllCellTypesAllArms.RData") #allphenotypes,UMAPlocs ,UMAPfiles,umapDImRedloc,umapDImRedfiles,nCellTypes,
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
uu <- unique( allphenotypes %>% dplyr::select( c("Celltype","PhenoCelltype", "key_", paste0("Disc_V", 1:5) ) ) )  

### Load Ligand Receptor database list of Ramilowski et al 2015
load( "~/Dropbox/Cancer_pheno_evo/data/FELINE2/LigandReceptor/Filtered_Human-2015-Ramilowski-LR-pairs.RData")
LRgenelist <- unique( c(LRpairsFiltered$HPMR.Receptor, LRpairsFiltered$HPMR.Ligand) )

### Load some signalling data for downstream analysis
savelocCCI <- "~/Dropbox/Cancer_pheno_evo/data/FELINE2/CommunicationOutputAllArms/"

# Are we studying the per indiv or per cell type level of signalling?
perIndiv=FALSE

# Names of LR data to analyze
filenamesCCI <- list.files(savelocCCI)#%>%length

# First let's do this for the population level
load(file ="~/Dropbox/Cancer_pheno_evo/data/FELINE2/Communication output merged ALL ARMS/PopulationCommunicationMerged ALL ARMS.RData" )#CCI,allphenotypes, uu,perIndiv,
CCI[,Treat:="CombinationRibo"]
CCI[ARM=="A",Treat:="LetrozoleAlone"]
CCI[,scaleTransduction:=scale(TransductionMu,center=F),by=c("Pair.Name")] 
CCI[!is.finite(TransductionMu),scaleTransduction:=0]

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
     m3 <- lm(log(1+TransductionMu) ~ dynamic_class3 , data=yi[Day==180])      
     m4 <- lm(log(1+TransductionMu) ~ dynamic_class3 , data=yi[Day==0])     
     
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
    
       resAll <- rbind(res,res2, res3, res4)#, res5)
    return(resAll)
  }else{
    return(NULL)
  }
}))
setnames(assessment,old=c("rn","Pr(>|t|)"),new=c("coefficient","pval"))
# FDR correction
assessment <- data.table(assessment%>%group_by(Treat,Model)%>%dplyr::mutate(adjust.p.val= p.adjust(pval, method="fdr")))


# Subset signaificant terms
responseRelated <- assessment[adjust.p.val<0.05][ coefficient!="(Intercept)"]

# Count minimum number ofsamples that are contribuing to statistical effect at different timepoints 
repdd <- CCI%>%group_by(Day,Pair.Name,LigandPhenoCelltype,ReceptorPhenoCelltype,dynamic_class3,Treat)%>%dplyr::summarise(replic=length(TransductionMu))%>%group_by(Pair.Name,LigandPhenoCelltype,ReceptorPhenoCelltype,dynamic_class3,Treat)%>%
  dplyr::summarise(replic=min(replic))
# Get the significant effects that are identified when there are more than three samples at each timepoint
responseRelated2_0 <- merge(responseRelated[order(adjust.p.val)] ,repdd, by=c("Pair.Name","LigandPhenoCelltype","ReceptorPhenoCelltype","Treat"),all.x=T)[replic>3][order(adjust.p.val)]
responseRelated2_0[Model=="one way anova_DayStart_Response"][Treat=="CombinationRibo"]$LigandPhenoCelltype%>%table()
responseRelated2_0[Model=="one way anova_DayStart_Response"][Treat=="CombinationRibo"]$ReceptorPhenoCelltype%>%table()
# save(responseRelated2_0,responseRelated,assessment,file="~/Dropbox/Cancer_pheno_evo/data/FELINE2/When cells communicate_communication changes AllArms/trendsby response AllArms.RData")
#
load(file="~/Dropbox/Cancer_pheno_evo/data/FELINE2/When cells communicate_communication changes AllArms/trendsby response AllArms.RData")

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

initdiffs <- unique(responseRelated2_0[Model=="one way anova_DayStart_Response"][][] %>% dplyr::select(-dynamic_class3))

initdiffs$LigandPhenoCelltype <-factor(initdiffs$LigandPhenoCelltype, 
                                       levels= data.table(initdiffsCountLong%>%group_by(LigandPhenoCelltype)%>%dplyr::summarise(n=sum(value)))[order(-n)]$LigandPhenoCelltype )
initdiffs$ReceptorPhenoCelltype <-factor(initdiffs$ReceptorPhenoCelltype, 
                                         levels= data.table(initdiffsCountLong%>%group_by(ReceptorPhenoCelltype)%>%dplyr::summarise(n=sum(value)))[order(-n)]$ReceptorPhenoCelltype )
pathfig<- "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/LR initial differences/"

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

pathfig <- "~/Dropbox/FELINE Project (1)/Manuscript  Feline immune communication/Figures Communication Project/Paper figures Cancer immune communication/"
#ggsave(filename=paste0(pathfig,"Ribo and Letrozole SI number of differences in communication.png"))


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

