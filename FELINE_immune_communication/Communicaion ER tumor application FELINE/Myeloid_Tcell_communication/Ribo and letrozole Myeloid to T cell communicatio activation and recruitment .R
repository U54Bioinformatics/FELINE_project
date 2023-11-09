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

resultoutFUN <-function(){
  load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Macrophage T cell communication/Macrophage T cell communication.RData")
  return(resultout)#subsetInflam,plotComm,sigdiffs ,resultout,
}
resultout<-resultoutFUN()
sigdiffsFUN <-function(){
  load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Macrophage T cell communication/Macrophage T cell communication.RData")
  return(sigdiffs)#subsetInflam,plotComm,sigdiffs ,resultout,
}
sigdiffs<-sigdiffsFUN()

Signifresults <- resultout[grep("Response",rn)][adj.pval<0.05] [order(Pair.Name,rn)]

subsetInflam[,"Tumorresponse":= "Resistant"]
subsetInflam[dynamic_class3=="Response","Tumorresponse":= "Sensitive"]



ggplot( subsetInflam[ grep("IL15_" , Pair.Name)] , aes(y= log(1 + scaleTransduction), log(1 + Day) )) + theme_classic() +
  geom_smooth(method="gam",formula=y~s(x,k=3),se=F) +theme(aspect.ratio=1)+
  geom_jitter(aes(group=Day,col=Pair.Name ),height=0,width=0.4) + facet_wrap(Treatmentlab~Tumorresponse) + 
  scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))+
  scale_color_discrete(name="Communication \n pathway")+
  labs(y="Macrophage to T cell communication",x="Day")


ggplot( subsetInflam[ grep("IL18_" , Pair.Name)][] , aes(y= log(1 + scaleTransduction), log(1 + Day) )) + theme_classic() +
  geom_smooth(method="gam",formula=y~s(x,k=3),se=F) +theme(aspect.ratio=1)+
  geom_jitter(aes(group=Day,col=Pair.Name ),height=0,width=0.4) + facet_wrap(Treatmentlab~Tumorresponse) +  
  scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))+
  scale_color_discrete(name="Communication \n pathway")+
  labs(y="Macrophage to T cell communication",x="Day")


plotthis<- subsetInflam[ grepl("IL18_" , Pair.Name)| grepl("IL15_" , Pair.Name)]
plotthis[, Pair.Name2 :=gsub("_","-",Pair.Name) ]
ggplot( plotthis , aes(y= log(1 + scaleTransduction), log(1 + Day),fill= Tumorresponse)) +
  theme_classic(base_size=18) +
  geom_smooth(alpha=1,col= "black",method="gam",formula=y~s(x,k=3),se=T) +theme(aspect.ratio=1)+
  geom_jitter(aes(group=Day,col=Pair.Name2 ),height=0,width=0.4) + facet_wrap(Treatmentlab~Tumorresponse) +  
  scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))+
  scale_color_discrete(name="Communication \n pathway")+
  scale_fill_npg(name="Tumor response")+
  labs(y="Myeloid to T cell \n activation communication",x="Day")
ggsave(file="~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Tcell communication phenotypes ALLArms/Ribo and letrozole Myeloid to T cell IL15 and 18 activation.png",width=8)


paths <- sigdiffs[grep( "ITG", Pair.Name) ][order(adj.pval)][1:10]$Pair.Name
plotthis<- subsetInflam[ Pair.Name%in%paths]
plotthis[, Pair.Name2 :=gsub("_","-",Pair.Name) ]
ggplot( plotthis , aes(y= log(1 + scaleTransduction), log(1 + Day),fill= Tumorresponse)) +
  theme_classic(base_size=18) +
  geom_smooth(alpha=1,col= "black",method="gam",formula=y~s(x,k=3),se=T) +theme(aspect.ratio=1)+
  geom_jitter(aes(group=Day,col=Pair.Name2 ),height=0,width=0.4) + facet_wrap(Treatmentlab~Tumorresponse) +  
  scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))+
  scale_color_discrete(name="Communication \n pathway")+
  scale_fill_npg(name="Tumor response")+
  labs(y="Myeloid to T cell \n recruitment communication",x="Day")
ggsave(file="~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Tcell communication phenotypes ALLArms/Ribo and letrozole Myeloid to T cell recruitment.png",width=8)



paths <- sigdiffs[grep( "ITG", Pair.Name) ][order(adj.pval)][1:10]$Pair.Name
plotthis<- subsetInflam[ Pair.Name%in%paths]
plotthis[, Pair.Name2 :=gsub("_","-",Pair.Name) ]

plotthis<- subsetInflam[Pair.Name%in%sigdiffs[grep( "ITG", Pair.Name) ][order(adj.pval)][1:10]$Pair.Name |
                          ( (grepl("IL18_" , Pair.Name)| grepl("IL15_" , Pair.Name))& Pair.Name%in% sigdiffs[grep( "IL", Pair.Name) ]$Pair.Name )
                        ]
plotthis[, Pair.Name2 :=gsub("_","-",Pair.Name) ]
plotthis[, Function := "Activation" ]
plotthis[Pair.Name%in%sigdiffs[grep( "ITG", Pair.Name) ][order(adj.pval)][1:10]$Pair.Name, Function := "Recruitment" ]

plotthis$Pair.Name2 <- factor(plotthis$Pair.Name2 , levels= gsub("_","-", c(sort(unique(subsetInflam[ grepl("IL18_" , Pair.Name)| grepl("IL15_" , Pair.Name)]$Pair.Name)) ,
                                                              sort(paths) )) )
ggplot( plotthis[Day==180]%>%group_by(Pair.Name2)%>%mutate(y=scale( log(1 + scaleTransduction))) , aes(y= y, Tumorresponse,group=interaction(Tumorresponse,Pair.Name2),col= Pair.Name2,fill= Pair.Name2)) +
  theme_classic(base_size=18) +
  geom_boxplot(alpha = 0.6,position=position_dodge(width=1))+
  #geom_smooth(alpha=1,col= "black",method="gam",formula=y~s(x,k=3),se=T) +
  theme(aspect.ratio=1)+
  geom_point(position=position_dodge(width=1))+ facet_wrap(Treatmentlab~Function)+
  geom_vline(xintercept = 1.5,linetype="dashed") +
  scale_color_discrete(name="Communication \n pathway")+
  scale_fill_discrete(name="Communication \n pathway")+
  labs(y="Myeloid to T cell \n communication post treatment",x="Tumor response")
ggsave(file="~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Tcell communication phenotypes ALLArms/Ribo and letrozole Myeloid to T cell post treatment activation and recruitment.png",width=8)
paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Ribo and letrozole Myeloid to T cell post treatment activation and recruitment.png"),height=10,width=10)

ggplot( plotthis[Day==180]%>%group_by(Pair.Name2)%>%mutate(y=scale( log(1 + scaleTransduction))) , aes(y= y, Tumorresponse,group=interaction(Tumorresponse,Pair.Name2),col= Pair.Name2,fill= Pair.Name2)) +
  theme_classic(base_size=18) +
  geom_boxplot(alpha = 0.6,position=position_dodge(width=1))+
  #geom_smooth(alpha=1,col= "black",method="gam",formula=y~s(x,k=3),se=T) +
  theme(aspect.ratio=1)+
  geom_point(position=position_dodge(width=1))+ facet_wrap(Treatmentlab~Function)+
  geom_vline(xintercept = 1.5,linetype="dashed") +
  scale_color_discrete(name="Communication \n pathway")+
  scale_fill_discrete(name="Communication \n pathway")+
  labs(y="Myeloid to T cell \n communication post treatment",x="Tumor response") +
  theme(axis.title=element_blank(),  axis.text=element_blank(),strip.text = element_blank(),legend.text=element_blank(),legend.title=element_blank())

paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"BLANK Ribo and letrozole Myeloid to T cell post treatment activation and recruitment.png"),height=10,width=10)

plotthis[,Daylab:=paste0("Day ",Day)]

ggplot( plotthis[Function=="Activation"][Day==180]%>%group_by(Pair.Name2)%>%mutate(y=scale( log(1 + scaleTransduction))) , aes(y= y, Tumorresponse,group=interaction(Tumorresponse,Pair.Name2),col= Pair.Name2,fill= Pair.Name2)) +
  theme_classic(base_size=26) +
  stat_boxplot(geom ='errorbar',position=position_dodge(width=1),col="black") + 
  geom_boxplot(alpha = 0.6,col="black",position=position_dodge(width=1))+
  #geom_smooth(alpha=1,col= "black",method="gam",formula=y~s(x,k=3),se=T) +
  theme(aspect.ratio=1)+
  geom_point(position=position_dodge(width=1))+ facet_wrap(Treatmentlab~Daylab, ncol=1)+
  geom_vline(xintercept = 1.5,linetype="dashed") +
  #scale_color_discrete(name="Communication \n pathway")+
  scale_colour_brewer(name="Communication \n pathway", palette = "YlGn") +
  scale_fill_brewer(name="Communication \n pathway", palette = "YlGn") +
  labs(y="Myeloid to T cell \n activation communications post treatment",x="Tumor response")
#ggsave(file="~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Tcell communication phenotypes ALLArms/Ribo and letrozole Myeloid to T cell post treatment activation communications",width=8)
paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Ribo and letrozole Myeloid to T cell post treatment activation communications.png"),height=10,width=10)


ggplot( plotthis[Function=="Activation"][Day==180]%>%group_by(Pair.Name2)%>%mutate(y=scale( log(1 + scaleTransduction))) , aes(y= y, Tumorresponse,group=interaction(Tumorresponse,Pair.Name2),col= Pair.Name2,fill= Pair.Name2)) +
  theme_classic(base_size=26) +
  stat_boxplot(geom ='errorbar',position=position_dodge(width=1),col="black") + 
  geom_boxplot(alpha = 0.6,col="black",position=position_dodge(width=1))+
  #geom_smooth(alpha=1,col= "black",method="gam",formula=y~s(x,k=3),se=T) +
  theme(aspect.ratio=1)+
  geom_point(position=position_dodge(width=1))+ facet_wrap(Treatmentlab~Daylab, ncol=1)+
  geom_vline(xintercept = 1.5,linetype="dashed") +
  #scale_color_discrete(name="Communication \n pathway")+
  scale_colour_brewer(name="Communication \n pathway", palette = "YlGn") +
  scale_fill_brewer(name="Communication \n pathway", palette = "YlGn") +
  labs(y="Myeloid to T cell \n communication post treatment",x="Tumor response") +
  theme(axis.title=element_blank(),  axis.text=element_blank(),strip.text = element_blank(),legend.text=element_blank(),legend.title=element_blank())
paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"BLANK Ribo and letrozole Myeloid to T cell post treatment activation communications.png"),height=10,width=10)




ggplot( plotthis[Function=="Recruitment"][Pair.Name!="VCAM1_ITGB2"][Day==180]%>%group_by(Pair.Name2)%>%mutate(y=scale( log(1 + scaleTransduction))) , aes(y= y, Tumorresponse,group=interaction(Tumorresponse,Pair.Name2),col= Pair.Name2,fill= Pair.Name2)) +
  theme_classic(base_size=26) +
  stat_boxplot(geom ='errorbar',position=position_dodge(width=1),col="black") + 
  geom_boxplot(alpha = 0.6,col="black",position=position_dodge(width=1))+
  #geom_smooth(alpha=1,col= "black",method="gam",formula=y~s(x,k=3),se=T) +
  theme(aspect.ratio=1)+
  geom_point(position=position_dodge(width=1))+ facet_wrap(Treatmentlab~Daylab, ncol=1)+
  geom_vline(xintercept = 1.5,linetype="dashed") +
  #scale_color_discrete(name="Communication \n pathway")+
  scale_colour_brewer(name="Communication \n pathway", palette = "YlGn") +
  scale_fill_brewer(name="Communication \n pathway", palette = "YlGn") +
  labs(y="Myeloid to T cell \n recruitment communications post treatment",x="Tumor response")
paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Ribo and letrozole Myeloid to T cell post treatment recruitment communications.png"),height=10,width=10)


ggplot( plotthis[Function=="Recruitment"][Pair.Name!="VCAM1_ITGB2"][Day==180]%>%group_by(Pair.Name2)%>%mutate(y=scale( log(1 + scaleTransduction))) , aes(y= y, Tumorresponse,group=interaction(Tumorresponse,Pair.Name2),col= Pair.Name2,fill= Pair.Name2)) +
  theme_classic(base_size=28) +
  stat_boxplot(geom ='errorbar',position=position_dodge(width=1),col="black") + 
  geom_boxplot(alpha = 0.6,col="black",position=position_dodge(width=1))+
  #geom_smooth(alpha=1,col= "black",method="gam",formula=y~s(x,k=3),se=T) +
  theme(aspect.ratio=1)+
  geom_point(position=position_dodge(width=1))+ facet_wrap(~Treatmentlab, ncol=1)+
  geom_vline(xintercept = 1.5,linetype="dashed") +
  #scale_color_discrete(name="Communication \n pathway")+
  scale_colour_brewer(name="Communication \n pathway", palette = "YlGn") +
  scale_fill_brewer(name="Communication \n pathway", palette = "YlGn") +
  labs(y="Myeloid to T cell \n recruitment communication post treatment",x="Tumor response") +
  theme(axis.title=element_blank(),  axis.text=element_blank(),strip.text = element_blank(),legend.text=element_blank(),legend.title=element_blank())
paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"BLANK Ribo and letrozole Myeloid to T cell post treatment recruitment communications.png"),height=10,width=10)

#save( plotthis, file= "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Macrophage T cell communication/Ribo and Letrozole Macrophage T cell activation and recruitment communication Cohort1.RData" )
load(  file= "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Macrophage T cell communication/Ribo and Letrozole Macrophage T cell activation and recruitment communication Cohort1.RData" )

summary(lm( log(1 + scaleTransduction) ~dynamic_class3, plotthis[Day==180][Treatmentlab=="Combination ribociclib"][Pair.Name%in% "IL15_IL15RA"] ))
summary(lm( log(1 + scaleTransduction) ~dynamic_class3, plotthis[Day==180][Treatmentlab=="Combination ribociclib"][Pair.Name%in% "IL15_IL2RA"] ))
#summary(lm( log(1 + scaleTransduction) ~dynamic_class3, plotthis[Day==180][Treatmentlab=="Combination ribociclib"][Pair.Name%in% "IL15_IL2RB"] ))

summary(lm( log(1 + scaleTransduction) ~dynamic_class3, plotthis[Day==180][Treatmentlab=="Combination ribociclib"][Pair.Name%in% "IL18_IL18R1"] ))
#summary(lm( log(1 + scaleTransduction) ~dynamic_class3, plotthis[Day==180][Treatmentlab=="Combination ribociclib"][Pair.Name%in% "IL18_IL18RAP"] ))
summary(lm( log(1 + scaleTransduction) ~dynamic_class3, plotthis[Day==180][Treatmentlab=="Combination ribociclib"][Pair.Name%in% "ADAM12_ITGB1"] ))
summary(lm( log(1 + scaleTransduction) ~dynamic_class3, plotthis[Day==180][Treatmentlab=="Combination ribociclib"][Pair.Name%in% "VCAM1_ITGB2"] ))
summary(lm( log(1 + scaleTransduction) ~dynamic_class3, plotthis[Day==180][Treatmentlab=="Combination ribociclib"][Pair.Name%in% "F13A1_ITGB1"] ))
summary(lm( log(1 + scaleTransduction) ~dynamic_class3, plotthis[Day==180][Treatmentlab=="Combination ribociclib"][Pair.Name%in% "F13A1_ITGB1"] ))
summary(lm( log(1 + scaleTransduction) ~dynamic_class3, plotthis[Day==180][Treatmentlab=="Combination ribociclib"][Pair.Name%in% "VCAM1_ITGA4"] ))



# reshape
TcellsummaryInteractionGF2 <- data.table( TcellsummaryInteractionGF %>%
                                            dplyr::select(Patient.Study.ID, Day, dynamic_class3, ARM,  key_, LigandPhenoCelltype, ReceptorPhenoCelltype, muln_scaleTransduction )%>%
                                            spread(LigandPhenoCelltype,muln_scaleTransduction) )
TcellsummaryInteractionGF2[,Treatmentlab:= "Combination ribociclib"]
TcellsummaryInteractionGF2[ARM=="A",Treatmentlab:= "Letrozole alone"]

macrocelldd <- allphenotypes[Celltype_subtype%in% c("CD8+ T cells", "NK cells")]
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
macrophageTcelcomm <- merge( fulldd, TcellsummaryInteractionGF2 , by=c("key_", "Patient.Study.ID" , "Day", "ARM", "dynamic_class3" )) 


ggplot( macrophageTcelcomm , aes(y=FromMacrophages, x=log(1+Day), fill=dynamic_class3, group=interaction(dynamic_class3, Day) ))+
  theme_classic(base_size=18)+
  geom_violin(position=position_dodge() )+#geom_smooth(method="lm")+
  geom_point(position =position_jitterdodge(dodge.width=2.3,jitter.width=0.0)  )+
  theme(aspect.ratio=1)+scale_x_continuous(breaks=log(1+c(0,14,180)) ,labels=c(0,14,180))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant","Sensitive"))+
  labs(y="Inflammatory cytokine signaling \n from macrophages to T cells", x="Day")+facet_wrap(~Treatmentlab)


ggplot( macrophageTcelcomm , aes(y=FromMacrophages, x=log(1+Day), fill=dynamic_class3, group=interaction(dynamic_class3, Day) ))+
  theme_classic(base_size=18)+
  geom_boxplot(outlier.colour=NA, position= position_dodge() )+
  stat_boxplot(geom="errorbar",position=position_dodge(1.75),width=0.5)+#geom_smooth(method="lm")+
  geom_point(position =position_jitterdodge(dodge.width=2.3,jitter.width=0.0)  )+
  theme(aspect.ratio=1)+scale_x_continuous(breaks=log(1+c(0,14,180)) ,labels=c(0,14,180))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant","Sensitive"))+
  labs(y="Macrophages to CD8 T cell \n inflammatory cytokine \n communications", x="Day")+facet_wrap(~Treatmentlab, nrow=2)

ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Tcell communication phenotypes ALLArms/Letrozole and Ribo Inflammatory signals from macrophages to T cells over time.png",width=6)


paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Ribo and Letrozole Inflammatory signals from macrophages to T cells over time.png"),height=10,width=10)

ggplot( macrophageTcelcomm , aes(y=FromMacrophages, x=log(1+Day), fill=dynamic_class3, group=interaction(dynamic_class3, Day) ))+
  theme_classic(base_size=18)+
  geom_boxplot(outlier.colour=NA, position= position_dodge() )+
  stat_boxplot(geom="errorbar",position=position_dodge(1.75),width=0.5)+#geom_smooth(method="lm")+
  geom_point(position =position_jitterdodge(dodge.width=2.3,jitter.width=0.0)  )+
  theme(aspect.ratio=1)+scale_x_continuous(breaks=log(1+c(0,14,180)) ,labels=c(0,14,180))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant","Sensitive"))+
  labs(y="Macrophages to CD8 T cell \n inflammatory cytokine \n communications", x="Day")+facet_wrap(~Treatmentlab, nrow=2)+
  theme(axis.title=element_blank(),  axis.text=element_blank(),strip.text = element_blank(),legend.title=element_blank(),legend.text=element_blank(),legend.position="none" )

paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
ggsave(paste0(paperfile,"BLANK Ribo and Letrozole Inflammatory signals from macrophages to T cells over time.png"),height=10,width=10)






#save( macrophageTcelcomm, file= "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Macrophage T cell communication/Ribo and Letrozole Macrophage T cell communication summary Cohort1.RData" )
load(  file= "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Macrophage T cell communication/Ribo and Letrozole Macrophage T cell communication summary Cohort1.RData" )

require(lme4)
require(lmerTest)
MR<- lmer(FromMacrophages~log(1+Day)*dynamic_class3+(1|Patient.Study.ID) , data=macrophageTcelcomm[ARM!="A"])
ML <- lmer(FromMacrophages~log(1+Day)*dynamic_class3+(1|Patient.Study.ID) , data=macrophageTcelcomm[ARM=="A"])
summary(MR)
summary(ML)


M1 <- lmer(FromMacrophages~dynamic_class3+(1|Patient.Study.ID) , data=macrophageTcelcomm[Day==14][ARM!="A"])
M1 <- lmer(FromMacrophages~dynamic_class3+(1|Patient.Study.ID) , data=macrophageTcelcomm[Day==14][ARM=="A"])
summary(M1)


# load T cell phenotype data
load( file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/T cell ssgsea phenotype C1/T cell ssgsea phenotype C1.RData")
#TcellPhenotytpedata
macrophageTcelcomm2<-merge( macrophageTcelcomm,  TcellPhenotytpedata%>%dplyr::select(Cell.ID,Tcellactivation), by="Cell.ID" )
macrophageTcelcomm2[,MacroCytokineCommunication:="Hi"]
macrophageTcelcomm2[FromMacrophages<mean(FromMacrophages),MacroCytokineCommunication:="Low"]

macrophageTcelcomm2[,TcellStatus:= "Effector"]
macrophageTcelcomm2[Tcellactivation<mean(Tcellactivation),TcellStatus:= "Naive"]
macrophageTcelcomm2[,roundFromMacrophages:=round(17*FromMacrophages)/17]
#macrophageTcelcomm2[,roundFromMacrophages:=round(30*FromMacrophages)/30]
#macrophageTcelcomm2[,roundFromMacrophages:=(100*FromMacrophages)/100]

#macrophageTcelcomm2[,roundFromMacrophages:=(round(40*sqrt(FromMacrophages))/40)^2]


ggplot(macrophageTcelcomm2[],aes(x= (FromMacrophages), y=Tcellactivation, col= roundFromMacrophages,fill= roundFromMacrophages))+
  theme_classic(base_size = 22)+
  #geom_boxplot(alpha=0.5, aes(group= interaction(roundFromMacrophages)),position=position_dodge())+
  #geom_dotplot(dotsize=0.2,stackratio=0.1,binaxis="y", stackdir = "center",aes(group= interaction(roundFromMacrophages)))+
  geom_point(size=3,alpha=1)+
  scale_color_viridis(option= "A" ,end=0.9)+
  scale_fill_viridis(option= "A" ,end=0.9)+
  geom_smooth(method= "lm") +
  theme(aspect.ratio= 1, legend.position="none") + 
  labs(y="CD8 T cell \n cytotoxicity phenotype", x= "Macrophage to CD8 T cell \n inflammatory communication" )
ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Tcell communication phenotypes ALLArms/Letrozole and Ribo Together Inflammatory signals from macrophages to T cells vs T cell phenotype.png",width=8)



ggplot(macrophageTcelcomm2[],aes(x= (roundFromMacrophages), y=Tcellactivation, col= roundFromMacrophages,fill= roundFromMacrophages))+
  theme_classic(base_size = 22)+
  geom_boxplot(alpha=0.5, aes(group= interaction(roundFromMacrophages)),position=position_dodge())+
  #geom_dotplot(dotsize=0.2,stackratio=0.1,binaxis="y", stackdir = "center",aes(group= interaction(roundFromMacrophages)))+
  geom_point(alpha=1)+
  scale_color_viridis(option= "A" ,end=0.9)+
  scale_fill_viridis(option= "A" ,end=0.9)+
  geom_smooth(method= "lm") + theme(aspect.ratio= 1, legend.position="none") + 
  labs(y="CD8+ T cell \n cytotoxicity phenotype", x= "Myeloid to CD8+ T cell \n inflammatory communication" )
#ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Tcell communication phenotypes ALLArms/Letrozole and Ribo TogetherB Inflammatory signals from macrophages to T cells vs T cell phenotype.png",width=8)
paperfile <- "/Users/jason/Dropbox/FELINE Project (1)/Manuscript  Feline immune communication/Figures Communication Project/Paper figures Cancer immune communication/"
ggsave(paste0(paperfile,"Ribo and Letrozole SI myeloid inflammatory communication to CD8 phenotype.png"),height=10,width=10)

summary(lm(Tcellactivation~roundFromMacrophages,data=macrophageTcelcomm2))
summary(lm(Tcellactivation~FromMacrophages,data=macrophageTcelcomm2))

ggplot(macrophageTcelcomm2[],aes(x= (roundFromMacrophages), y=Tcellactivation, col= roundFromMacrophages, fill= roundFromMacrophages))+
  theme_classic(base_size = 22)+
  geom_boxplot(alpha=0.5, aes(group= roundFromMacrophages))+
  geom_point(alpha=1)+
  scale_color_viridis(option= "A" )+
  scale_fill_viridis(option= "A" )+
  geom_smooth(method= "lm") + theme(aspect.ratio= 1, legend.position="none") + facet_wrap(~Treatmentlab)+
  labs(y="CD8 T cell \n cytotoxicity phenotype", x= "Myeloid to CD8 T cell \n inflammatory communication" )
ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Tcell communication phenotypes ALLArms/Letrozole and Ribo Inflammatory signals from macrophages to T cells vs T cell phenotype.png",width=8)

# 
# summtcelcompheno<-data.table(macrophageTcelcomm2%>%group_by(Day,dynamic_class3,Treatmentlab,TcellStatus,MacroCytokineCommunication)%>%summarise(n=n()))
# summtcelcompheno
# ggplot(summtcelcompheno[Day==180],aes(x=MacroCytokineCommunication,y=TcellStatus,fill=n ))+theme_classic()+
#   geom_tile()+facet_wrap(Treatmentlab~Day)+theme(aspect.ratio=1)
# geom_point(aes(col=TcellStatus))+
#   geom_smooth(method="lm")+theme(aspect.ratio=1)+
#   #geom_smooth(method="gam",formula=y~s(x,k=3))+theme(aspect.ratio=1)+
#   facet_wrap(Treatmentlab~Day)
# ggplot(macrophageTcelcomm2[Day==180],aes(x=(FromMacrophages),y=Tcellactivation ))+theme_classic()+
#   geom_point(aes(col=TcellStatus))+
#   geom_smooth(method="lm")+theme(aspect.ratio=1)+
#   #geom_smooth(method="gam",formula=y~s(x,k=3))+theme(aspect.ratio=1)+
#   facet_wrap(Treatmentlab~Day)
# ggplot(macrophageTcelcomm2[Day==180],aes(x=(FromMacrophages),y=Tcellactivation ))+theme_classic()+
#   geom_point(aes(col=MacroCytokineCommunication))+
#   geom_smooth(method="lm")+theme(aspect.ratio=1)+
#   #geom_smooth(method="gam",formula=y~s(x,k=3))+theme(aspect.ratio=1)+
#   facet_wrap(Treatmentlab~Day)
# 
# ggplot(macrophageTcelcomm2[],aes(x=(FromMacrophages),y=Tcellactivation ))+theme_classic()+
#   geom_point(aes(col=dynamic_class3))+
#   geom_smooth(method="lm")+theme(aspect.ratio=1)+
#   #geom_smooth(method="gam",formula=y~s(x,k=3))+theme(aspect.ratio=1)+
#   facet_wrap(Treatmentlab~Day)
# 
