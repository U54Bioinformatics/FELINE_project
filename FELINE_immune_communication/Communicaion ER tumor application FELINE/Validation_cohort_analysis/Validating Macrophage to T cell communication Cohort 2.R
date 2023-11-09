rm(list=ls())
require(data.table)
require(dplyr)
require(ggplot2)
require(tidyr)
require(ggsci)

require(lme4)
require(lmerTest)
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
allphenotypes[Patient.Study.ID=="001-171_853686",Patient.Study.ID:="001-171"]
load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE/Feline clinical input.RData")
allphenotypes <- merge(allphenotypes, unique(Clin_resp_dd_classAdd%>%dplyr::select(Patient.Study.ID,prop_change,dynamic_class)), by="Patient.Study.ID")

### Unique clusters of cells (ALL) and their umap discretization level
uu <- unique( allphenotypes %>% dplyr::select( c("Celltype","PhenoCelltype", "key_", paste0("Disc_V", 1:5) ) ) )  #"Celltype_subtype",

#load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/PhenotypesAllArmsCohort2/geneUMAP_T_cellsCD8+ T cells.RData" )
#u_datssgsea,

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
#rm(list="TcellSubtypeCommCCI")

ggplot( TcellsummaryInteractionGF , aes(y=muln_scaleTransduction, log(1+Day) ))+
  geom_violin(aes(group=Day))+facet_wrap(~dynamic_class3)+geom_smooth(method="lm")


subsetInflam <- TcellSubtypeCommCCI[ARM!="A"][Pair.Name %in% LRpairsFiltered[
  HPMR.Receptor%in% cytokineReceptors
  ]$Pair.Name][LigandPhenoCelltype=="Macrophages"]

subsetInflam[,startval:=sum((Day==0)*log(1+scaleTransduction))/sum(Day==0), by=c("Pair.Name","Patient.Study.ID","ReceptorPhenoCelltype")]
subsetInflam[ ,startvalAv:=median(startval,na.rm=T), by=c("Pair.Name","dynamic_class3","ReceptorPhenoCelltype")]


pairstotest <- unique(subsetInflam$Pair.Name)

resultout <- rbindlist(lapply( 1:length(pairstotest),function(x){
  print (x)
  res<-tryCatch({
    m1 <- lmer( log(1+scaleTransduction) ~ log(1+Day) * dynamic_class3 +(1|Patient.Study.ID), data=   subsetInflam[ARM!="A"][Pair.Name== pairstotest[x] ])
    data.table( Pair.Name=pairstotest[x],  data.table(summary(m1)$coef,keep.rownames=T) )
  },
  error=function(x){
    return(NULL)
  })
  res
}))
setnames(resultout,old=c("t value","Pr(>|t|)"),new=c("tval","pval"))
resultout$adj.pval <- p.adjust(p=resultout$pval,method="fdr")#,n=length(unique(resultout$Pair.Name)))

#resultout[grep("IL18", Pair.Name)][grep("Resp", rn)]
#resultout[grep("IL15", Pair.Name)][grep("Resp", rn)]
#resultout[grep("CCR", Pair.Name)][grep("Resp", rn)][pval<0.05]



resultout[rn!="(Intercept)"][pval<0.05][order(pval)][1:30]
resultout[Pair.Name=="IL1B_IL1R1"]
resultout[rn!="(Intercept)"][pval<0.05][order(pval)][grepl("Resp",rn)]$Pair.Name%>%unique()
resultout[rn=="log(1 + Day):dynamic_class3Response"][adj.pval<0.05][order(pval)][][Estimate>0][!grepl("ITG", Pair.Name)][!grepl("CD44", Pair.Name)]
resultout[rn=="log(1 + Day):dynamic_class3Response"][adj.pval<0.05][order(pval)][][Estimate<0][!grepl("ITG", Pair.Name)][!grepl("CD44", Pair.Name)]
resultout[rn!="(Intercept)"][pval<0.05][order(pval)][grepl("Resp",rn)][Estimate<0]

Signifresults <- resultout[grep("Response",rn)][adj.pval<0.05] [order(Pair.Name,rn)]
write.csv(Signifresults, "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Final communication files/Macrophage to T cell communication/Macrophage to T cell communication Cohort 2.csv")

sigdiffs <- resultout[rn=="log(1 + Day):dynamic_class3Response"][adj.pval<0.05][order(pval)]
resultout[rn=="log(1 + Day):dynamic_class3Response"][adj.pval<0.05][order(Pair.Name)]
C1comparison <- resultout[order(Pair.Name)][Pair.Name%in%c("ADAM12_ITGB1","COL18A1_ITGB1","F13A1_ITGA4","F13A1_ITGB1", "ICAM2_ITGB2", "ICAM3_ITGB2","THBS2_ITGA4","VCAM1_ITGA4","VCAM1_ITGB1","VCAM1_ITGB2","IL15_IL15RA","IL15_IL2RA","IL15_IL2RB","IL18_IL18R1","IL18_IL18RAP","IL18_IL1R1","OSM_IL16ST","CCL19_CCR7","CCL2_CCR5","CCL21_CCR7","CCL3_CCR5","CCL8_CCR2","CCL8_CCR5")]
write.csv(C1comparison, "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Final communication files/Macrophage to T cell communication/Macrophage to T cell Validation communication from Cohort 1 in Cohort 2.csv")

sigdiffs[order(Pair.Name)]
plotComm<-function(R,n=NULL,exclude=NULL){
  if(is.null(n)){
    paths <- sigdiffs[grep( R , Pair.Name)][order(adj.pval)]$Pair.Name
    paths <- paths[!paths%in%exclude]
  }else{
    paths <- sigdiffs[grep( R , Pair.Name)][order(adj.pval)]$Pair.Name
    paths <- paths[!paths%in%exclude]
    if(n<length(paths)){ paths <- paths[1:n]}
  }
  tmp<- subsetInflam[ Pair.Name%in%paths][] 
  tmp[,"Tumorresponse":= "Resistant"]
  tmp[dynamic_class3=="Response","Tumorresponse":= "Sensitive"]
  ggplot( tmp , aes(y= log(1 + scaleTransduction)+startvalAv-startval, log(1 + Day) )) + theme_classic() +
    geom_jitter(aes(group=Day,col=Pair.Name ),height=0,width=0.4) + facet_wrap(~Tumorresponse) + 
    #geom_smooth(method="gam",formula=y~s(x,k=3),se=F)+ 
    geom_smooth(method="lm",se=F) + #aes(group=Pair.Name,col=Pair.Name),
    theme(aspect.ratio=1) +
    scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))+
    scale_color_discrete(name="Communication \n pathway")+
    labs(y="Macrophage to T cell communication",x="Day")
}

plotComm("ITG",n=10)
plotComm("CCR")
plotComm("CXC")
plotComm("TNFR")
plotComm("TGFBR")
plotComm("IL",exclude=c("HLA-C_LILRB1","IL10_IL10RA","ICAM1_IL2RA"))

# change    save(subsetInflam,plotComm,sigdiffs ,resultout,file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Macrophage T cell communication Cohort 2/Macrophage T cell communication Cohort2.RData")
load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Macrophage T cell communication Cohort 2/Macrophage T cell communication Cohort2.RData")
#subsetInflam,plotComm,sigdiffs ,resultout,
subsetInflam[,MacroEnriched :=F]
subsetInflam[Patient.Study.ID=="001-167",MacroEnriched :=T]
subsetInflam[,"Tumorresponse":= "Resistant"]
subsetInflam[dynamic_class3=="Response","Tumorresponse":= "Sensitive"]



plotthis <- TcellSubtypeCommCCI[][Pair.Name %in%c("IL15_IL15RA","IL15_IL2RA","IL15_IL2RB","IL18_IL18R1" ,"IL18_IL18RAP")    ][LigandPhenoCelltype=="Macrophages"]
plotthis[,"Tumorresponse":= "Resistant"]
plotthis[dynamic_class3=="Response","Tumorresponse":= "Sensitive"]

plotthis[, Pair.Name2 :=gsub("_","-",Pair.Name) ]
plotthis[, Function := "Activation" ]
plotthis$Pair.Name2 <- factor(plotthis$Pair.Name2 , levels= c("IL15-IL15RA","IL15-IL2RA","IL15-IL2RB","IL18-IL18R1" ,"IL18-IL18RAP") )

plotthis[,Treatmentlab:= "Combination ribociclib"]
plotthis[ARM=="A",Treatmentlab:= "Letrozole alone"]

plotthis[, y:=  scale( log(1 + scaleTransduction)),by=Pair.Name2]
plotthis[,y0:=sum(y*(Day==0))/sum(Day==0) , by=c("Patient.Study.ID","Pair.Name2")]
plotthis[,deltay:=y-y0 , by=c("Patient.Study.ID","Pair.Name2")]
# 

ggplot( plotthis[Function=="Activation"][Day==180][!(Patient.Study.ID=="001-167")],#%>%group_by(Pair.Name2)%>%mutate(y=scale( log(1 + scaleTransduction))) , 
        aes(y= deltay#y
            , Tumorresponse,group=interaction(Tumorresponse,Pair.Name2),col= Pair.Name2,fill= Pair.Name2)) +
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
  labs(y="Myeloid to CD8+ T cell \n activation communications post treatment",x="Tumor response")

paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
paperfile<- "/Users/jason/Dropbox/FELINE Project (1)/Manuscript  Feline immune communication/Figures Communication Project/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Validation IL communicaiton Ribo and letrozole Myeloid to T cell post treatment activation communications Validation C2.png"),height=10,width=10)



ggplot( plotthis[Function=="Activation"][Day==180][!(Patient.Study.ID=="001-167")]%>%group_by(Pair.Name2)%>%mutate(y=scale( log(1 + scaleTransduction))) , 
        aes(y= y, Tumorresponse,group=interaction(Tumorresponse,Pair.Name2),col= Pair.Name2,fill= Pair.Name2)) +
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
  labs(y="Myeloid to T cell \n activation communications post treatment",x="Tumor response")


#ggsave(file="~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Tcell communication phenotypes ALLArms/Ribo and letrozole Myeloid to T cell post treatment activation communications",width=8)
paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Ribo and letrozole Myeloid to T cell post treatment activation communications.png"),height=10,width=10)





plotthis <- TcellSubtypeCommCCI[][Pair.Name %in%c("IL15_IL15RA","IL15_IL2RA","IL15_IL2RB","IL18_IL18R1" ,"IL18_IL18RAP")    ][LigandPhenoCelltype=="Macrophages"]
plotthis[,"Tumorresponse":= "Resistant"]
plotthis[dynamic_class3=="Response","Tumorresponse":= "Sensitive"]

plotthis[, Pair.Name2 :=gsub("_","-",Pair.Name) ]
plotthis[, Function := "Activation" ]
plotthis$Pair.Name2 <- factor(plotthis$Pair.Name2 , levels= c("IL15-IL15RA","IL15-IL2RA","IL15-IL2RB","IL18-IL18R1" ,"IL18-IL18RAP") )

plotthis[,Treatmentlab:= "Combination ribociclib"]
plotthis[ARM=="A",Treatmentlab:= "Letrozole alone"]

plotthis[, y:=  scale( log(1 + scaleTransduction)),by=Pair.Name2]
plotthis[,y0:=sum(y*(Day==0))/sum(Day==0) , by=c("Patient.Study.ID","Pair.Name2")]
plotthis[,deltay:=y-y0 , by=c("Patient.Study.ID","Pair.Name2")]
# 

ggplot( plotthis[Function=="Activation"][Day==180][!(Patient.Study.ID=="001-167")],#%>%group_by(Pair.Name2)%>%mutate(y=scale( log(1 + scaleTransduction))) , 
        aes(y= deltay#y
            , Tumorresponse,group=interaction(Tumorresponse,Pair.Name2),col= Pair.Name2,fill= Pair.Name2)) +
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
  labs(y="Myeloid to CD8+ T cell \n activation communications post treatment",x="Tumor response")

paperfile<- "/Users/jason/Dropbox/FELINE Project (1)/Manuscript  Feline immune communication/Figures Communication Project/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Validation IL communicaiton Ribo and letrozole Myeloid to T cell post treatment activation communications Validation C2.png"),height=10,width=10)








plotthis <- TcellSubtypeCommCCI[][Pair.Name %in%c("ADAM12_ITGB1","COL18A1_ITGB1","F13A1_ITGA4","F13A1_ITGB1","ICAM2_ITGB2","ICAM3_ITGB2","THBS2_ITGA4","VCAM1_ITGA4","VCAM1_ITGB1")    ][LigandPhenoCelltype=="Macrophages"]
plotthis[,"Tumorresponse":= "Resistant"]
plotthis[dynamic_class3=="Response","Tumorresponse":= "Sensitive"]

plotthis[, Pair.Name2 :=gsub("_","-",Pair.Name) ]
plotthis[, Function := "Recruitment" ]
plotthis$Pair.Name2 <- factor(plotthis$Pair.Name2 , levels= c("ADAM12-ITGB1","COL18A1-ITGB1","F13A1-ITGA4","F13A1-ITGB1","ICAM2-ITGB2","ICAM3-ITGB2","THBS2-ITGA4","VCAM1-ITGA4","VCAM1-ITGB1") )

plotthis[,Treatmentlab:= "Combination ribociclib"]
plotthis[ARM=="A",Treatmentlab:= "Letrozole alone"]

plotthis[, y:=  scale( log(1 + scaleTransduction)),by=Pair.Name2]
plotthis[,y0:=sum(y*(Day==0))/sum(Day==0) , by=c("Patient.Study.ID","Pair.Name2")]
plotthis[,deltay:=y-y0 , by=c("Patient.Study.ID","Pair.Name2")]
# 

ggplot( plotthis[][Day==180][!(Patient.Study.ID=="001-167")],#%>%group_by(Pair.Name2)%>%mutate(y=scale( log(1 + scaleTransduction))) , 
        aes(y= deltay#y
            , Tumorresponse,group=interaction(Tumorresponse,Pair.Name2),col= Pair.Name2,fill= Pair.Name2)) +
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
  labs(y="Myeloid to CD8+ T cell \n activation communications post treatment",x="Tumor response")


#ggsave(file="~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Tcell communication phenotypes ALLArms/Ribo and letrozole Myeloid to T cell post treatment activation communications",width=8)
#paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Validation ITG communicaiton Ribo and letrozole Myeloid to T cell post treatment recruitment communications Validation C2.png"),height=10,width=10)












ggplot( subsetInflam[ grep("IL15_" , Pair.Name)][!(Patient.Study.ID=="001-167")], aes(y= log(1 + scaleTransduction), log(1 + Day) )) + theme_classic() +
  geom_smooth(method="gam",formula=y~s(x,k=3),se=F) +theme(aspect.ratio=1)+
  geom_jitter(aes(group=Day,col=Pair.Name ),height=0,width=0.4) + facet_wrap(~Tumorresponse) + 
  scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))+
  scale_color_discrete(name="Communication \n pathway")+
  labs(y="Macrophage to T cell communication",x="Day")


ggplot( subsetInflam[ grep("IL18_" , Pair.Name)][!(Patient.Study.ID=="001-167")], aes(y= log(1 + scaleTransduction), log(1 + Day) )) + theme_classic() +
  geom_smooth(method="gam",formula=y~s(x,k=3),se=F) +theme(aspect.ratio=1)+
  geom_jitter(aes(group=Day,col=Pair.Name ),height=0,width=0.4) + facet_wrap(~Tumorresponse) + 
  scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))+
  scale_color_discrete(name="Communication \n pathway")+
  labs(y="Macrophage to T cell communication",x="Day")

subsetInflam[, Pair.Name2 :=gsub("_","-",Pair.Name) ]




plotdatIL<- TcellSubtypeCommCCI[ grepl("IL18_" , Pair.Name)|grepl("IL15_" , Pair.Name)][LigandPhenoCelltype=="Macrophages"]#[!(Patient.Study.ID=="001-167")]
plotdatIL[,Treatmentlab:= "Combination ribociclib"]
plotdatIL[ARM=="A",Treatmentlab:= "Letrozole alone"]
plotdatIL[,"Tumorresponse":= "Resistant"]
plotdatIL[dynamic_class3=="Response","Tumorresponse":= "Sensitive"]
plotdatIL[, Pair.Name2 :=gsub("_","-",Pair.Name) ]

ggplot(plotdatIL , aes(y= log(1 + scaleTransduction), log(1 + Day) )) + 
  theme_classic(base_size=18) +
  geom_smooth(method="gam",formula=y~s(x,k=3),se=F) +theme(aspect.ratio=1)+
  geom_jitter(aes(group=Day,col=Pair.Name2 ),height=0,width=0.4) + facet_wrap(Treatmentlab~Tumorresponse) + 
  scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))+
  labs(y="Macrophage to T cell communication",x="Day") +
  scale_colour_brewer(name="Communication \n pathway", palette = "YlGn") +
  scale_fill_brewer(name="Communication \n pathway", palette = "YlGn")
paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Ribo and letrozole Myeloid to T cell validated post treatment activation communications.png"),height=10,width=10)


ggplot( subsetInflam[ grepl("IL18_" , Pair.Name)|grepl("IL15_" , Pair.Name)][Day==180], aes(y= log(1 + scaleTransduction), log(1 + Day) )) + theme_classic() +
  theme(aspect.ratio=1)+
  geom_vi
  geom_jitter(aes(group=Day,col=Pair.Name ),height=0,width=0.4) + facet_wrap(~Tumorresponse) + 
  scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))+
  scale_color_discrete(name="Communication \n pathway")+
  labs(y="Macrophage to T cell communication",x="Day")

paths <- unique(resultout[rn!="(Intercept)"][adj.pval<0.05][order(pval)][rn=="log(1 + Day):dynamic_class3Response"][grep("ITG" , Pair.Name)][]$Pair.Name)
#unique(subsetInflam[grep("ITG" , Pair.Name)]$Pair.Name)
#paths <- c("ADAM15_ITGAV","ADAM15_ITGB1","ADAM17_ITGB1","ADAM28_ITGA4","CHAD_ITGB1","COL4A1_ITGAV","COL4A1_ITGB1","F13A1_ITGA4","F13A1_ITGB1","FN1_ITGAV",
#           "ICAM4_ITGA4", "ICAM4_ITGB2",  "LAMA2_ITGB1","LAMB1_ITGB1","RELN_ITGB1","SEMA7A_ITGB1","THBS2_ITGA4","THBS2_ITGB1")
paths <- c("ADAM12_ITGB1","COL18A1_ITGB1","F13A1_ITGA4","F13A1_ITGB1", "ICAM2_ITGB2", "ICAM3_ITGB2","THBS2_ITGA4","VCAM1_ITGA4","VCAM1_ITGB1","VCAM1_ITGB2")
ggplot( subsetInflam[ Pair.Name%in%paths[-5]][!(Patient.Study.ID=="001-167")] 
        , aes(y= log(1 + scaleTransduction)+startvalAv-startval, log(1 + Day)  )) + theme_classic() +
  geom_jitter(aes(group=Day,col=Pair.Name ),height = 0,width=0.2) + facet_wrap(~Tumorresponse) +
  #geom_smooth(method="gam",formula=y~s(x,k=3),se=F) +
  geom_smooth(method="lm",se=F) +
  theme(aspect.ratio=1) +    labs(y="Macrophage to T cell communication",x="Day")+
  scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))+
  scale_color_discrete(name="Communication \n pathway")
resultout[rn!="(Intercept)"][ Pair.Name%in%paths[]][rn=="log(1 + Day):dynamic_class3Response"]




paths <- subsetInflam[grep("CD44" , Pair.Name)]$Pair.Name
ggplot( subsetInflam[ Pair.Name%in%paths][!(Patient.Study.ID=="001-167")] 
        , aes(y= log(1 + scaleTransduction)+startvalAv-startval, log(1 + Day)  )) + theme_classic() +
  geom_jitter(aes(group=Day,col=Pair.Name ),height = 0,width=0.2) + facet_wrap(~Tumorresponse) +
  #geom_smooth(method="gam",formula=y~s(x,k=3),se=F) +
  geom_smooth(method="lm",se=F) +
  theme(aspect.ratio=1) +    labs(y="Macrophage to T cell communication",x="Day")+
  scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))+
  scale_color_discrete(name="Communication \n pathway")



paths <- subsetInflam[grep("IL18" , Pair.Name)]$Pair.Name
ggplot( subsetInflam[ Pair.Name%in%paths][!(Patient.Study.ID=="001-167")] 
        , aes(y= log(1 + scaleTransduction)+startvalAv-startval, log(1 + Day)  )) + theme_classic() +
  geom_jitter(aes(group=Day,col=Pair.Name ),height = 0,width=0.2) + facet_wrap(~Tumorresponse) +
  #geom_smooth(method="gam",formula=y~s(x,k=3),se=F) +
  geom_smooth(method="lm",se=F) +
  theme(aspect.ratio=1) +    labs(y="Macrophage to T cell communication",x="Day")+
  scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))+
  scale_color_discrete(name="Communication \n pathway")

paths <- subsetInflam[grep("IL15" , Pair.Name)]$Pair.Name
paths <- c("IL15_IL15RA","IL15_IL2RA","IL15_IL2RB","IL18_IL18R1","IL18_IL18RAP","IL18_IL1R1","OSM_IL16ST")
ggplot( subsetInflam[ Pair.Name%in%paths][!(Patient.Study.ID=="001-167")] 
        , aes(y= log(1 + scaleTransduction)+startvalAv-startval, log(1 + Day)  )) + theme_classic() +
  geom_jitter(aes(group=Day,col=Pair.Name ),height = 0,width=0.2) + facet_wrap(~Tumorresponse) +
  #geom_smooth(method="gam",formula=y~s(x,k=3),se=F) +
  geom_smooth(method="lm",se=F) +
  theme(aspect.ratio=1) +    labs(y="Macrophage to T cell communication",x="Day")+
  scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))+
  scale_color_discrete(name="Communication \n pathway")

paths <- subsetInflam[grepl("CCR5" , Pair.Name)|grepl("CCR7" , Pair.Name)|grepl("CCR2" , Pair.Name)]$Pair.Name
paths <- c("CCL19_CCR7","CCL2_CCR5","CCL21_CCR7","CCL3_CCR5","CCL8_CCR2","CCL8_CCR5")
ggplot( subsetInflam[ Pair.Name%in%paths][!(Patient.Study.ID=="001-167")] 
        , aes(y= log(1 + scaleTransduction)+startvalAv-startval, log(1 + Day)  )) + theme_classic() +
  geom_jitter(aes(group=Day,col=Pair.Name ),height = 0,width=0.2) + facet_wrap(~Tumorresponse) +
  #geom_smooth(method="gam",formula=y~s(x,k=3),se=F) +
  geom_smooth(method="lm",se=F) +
  theme(aspect.ratio=1) +    labs(y="Macrophage to T cell communication",x="Day")+
  scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))+
  scale_color_discrete(name="Communication \n pathway")


ggplot( subsetInflam[ Pair.Name%in%paths][Patient.Study.ID=="001-167"] , aes(y= log(1 + scaleTransduction)+startvalAv-startval, log(1 + Day) )) + theme_classic() +
  geom_jitter(aes(group=Day,col=Pair.Name ),height = 0,width=0.2) + facet_wrap(ReceptorPhenoCelltype~dynamic_class3) +
  #geom_smooth(method="gam",formula=y~s(x,k=3),se=F) +
  geom_smooth(method="lm",se=F) +
  theme(aspect.ratio=1)


paths <- subsetInflam[grep("IL18" , Pair.Name)]$Pair.Name
ggplot( subsetInflam[ Pair.Name%in%paths][] , aes(y= log(1 + scaleTransduction)+startvalAv-startval, log(1 + Day) )) + theme_classic() +
  geom_jitter(aes(group=Day,col=Pair.Name ),height = 0,width=0.2) + facet_wrap(~dynamic_class3) +
  #geom_smooth(method="gam",formula=y~s(x,k=3),se=F) +
  geom_smooth(method="lm",se=F) +
  theme(aspect.ratio=1)


paths <- unique( subsetInflam[grep("ITGB" , Pair.Name)]$Pair.Name )
ggplot( subsetInflam[ Pair.Name%in%paths][] , aes(y= log(1 + scaleTransduction)+startvalAv-startval, log(1 + Day) )) + theme_classic() +
  geom_jitter(aes(group=Day,col=Pair.Name ),height = 0,width=0.2) + facet_wrap(~dynamic_class3) +
  #geom_smooth(method="gam",formula=y~s(x,k=3),se=F) +
  geom_smooth(method="lm",se=F) +
  theme(aspect.ratio=1)

paths <- unique(subsetInflam[grep("TNFRSF" , Pair.Name)]$Pair.Name)
ggplot( subsetInflam[ Pair.Name%in%paths][] , aes(y= log(1 + scaleTransduction)+startvalAv-startval, log(1 + Day) )) + theme_classic() +
  geom_jitter(aes(group=Day,col=Pair.Name ),height = 0,width=0.2) + facet_wrap(~dynamic_class3) +
  #geom_smooth(method="gam",formula=y~s(x,k=3),se=F) +
  geom_smooth(method="lm",se=F) +
  theme(aspect.ratio=1)




#paths <- subsetInflam[grep("CCR" , Pair.Name)]$Pair.Name
ggplot( subsetInflam[ Pair.Name%in%paths][Day!=180] , aes(y= log(1 + scaleTransduction)+startvalAv-startval, log(1 + Day) )) + theme_classic() +
  geom_jitter(aes(group=Day,col=Pair.Name )) +
  facet_wrap(~dynamic_class3) + 
  #geom_smooth(method="gam",formula=y~s(x,k=3),se=F) +
  geom_smooth(method="lm",se=F) +
  theme(aspect.ratio=1)



# reshape
TcellsummaryInteractionGF2 <- data.table( TcellsummaryInteractionGF %>%
                                            dplyr::select(Patient.Study.ID, Day, dynamic_class3, ARM,  key_, LigandPhenoCelltype, ReceptorPhenoCelltype, muln_scaleTransduction )%>%
                                            spread(LigandPhenoCelltype,muln_scaleTransduction) )
setnames(TcellsummaryInteractionGF2,old=c("Adipocytes", "B cells", "Cancer cells", "CD4+ T cells", "CD8+ T cells", "Endothelial cells", "Fibroblasts", "Macrophages", "Normal epithelial cells", "Pericytes"),
         new=c("FromAdipocytes", "FromBcells", "FromCancercells", "FromCD4Tcells", "FromCD8Tcells", "FromEndothelialcells", "FromFibroblasts", "FromMacrophages", "FromNormalepithelialcells", "FromPericytes"))
udataComm <- merge( allphenotypes, TcellsummaryInteractionGF2 , by=c("key_", "Patient.Study.ID" , "Day", "ARM", "dynamic_class3" )) 


#load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/PhenotypesAllArms/geneUMAP_T_cellsCD8+ T cells.RData" )
load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/PhenotypesAllArmsCohort2/geneUMAP_T_cellsCD8+ T cells.RData" )
# u_datssgsea
#macrophageTcelcomm <- merge( merge( u_datssgsea,
#                                    allphenotypes%>%select(Cell.ID,key_),by="Cell.ID")
                                   # , TcellsummaryInteractionGF2 , by=c("key_", "Patient.Study.ID" , "Day", "ARM", "dynamic_class3" )) 
macrophageTcelcomm <- udataComm[ARM!="A"]
macrophageTcelcomm[,startval:=sum((Day==0)*FromMacrophages)/sum(Day==0), by=c("Patient.Study.ID","ReceptorPhenoCelltype")]
macrophageTcelcomm[ ,startvalAv:=median(startval,na.rm=T), by=c("dynamic_class3","ReceptorPhenoCelltype")]

#save( macrophageTcelcomm,file= "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Macrophage T cell communication Cohort 2/Macrophage T cell communication summary Cohort2.RData" )
load(file= "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Macrophage T cell communication Cohort 2/Macrophage T cell communication summary Cohort2.RData" )

require(lme4)
require(lmerTest)
M1 <- lmer(FromMacrophages~log(1+Day)*dynamic_class3+(1+Day|Patient.Study.ID) , data=macrophageTcelcomm)
summary(M1)
#commSummaryPatlev<-data.table(macrophageTcelcomm[ARM!="A"]%>%group_by(Day,Patient.Study.ID)%>%summarise(FromMacrophages=median(FromMacrophages)) )
ggplot( macrophageTcelcomm[ARM!="A"][Patient.Study.ID!="001-167"] , aes(y=FromMacrophages+startvalAv-startval, log(1+Day) ))+geom_violin(aes(group=Day))+facet_wrap(ReceptorPhenoCelltype~dynamic_class3)+#geom_smooth(method="lm")+
  geom_jitter( aes(y=FromMacrophages+startvalAv-startval, log(1+Day)))

ggplot( macrophageTcelcomm[ARM!="A"][ReceptorPhenoCelltype=="CD8+ T cells"] [], aes(y=FromMacrophages+startvalAv-startval, x=log(1+Day), fill=dynamic_class3, group=interaction(dynamic_class3) ))+
  theme_classic(base_size=18)+#facet_wrap(~ReceptorPhenoCelltype)+
  geom_violin(alpha=0.4,aes(group=interaction(dynamic_class3, Day) ),scale="width" ,position=position_dodge() )+#geom_smooth(method="lm")+
  geom_point(position =position_jitterdodge(dodge.width=2.3,jitter.width=0.0)  )+
  theme(aspect.ratio=1)+scale_x_continuous(breaks=log(1+c(0,14,180)) ,labels=c(0,14,180))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant","Sensitive"))+
  scale_color_npg(name="Tumor response", labels=c("Resistant","Sensitive"))+
  labs(y="Inflammatory cytokine signaling \n from macrophages to T cells", x="Day")+
  geom_smooth(method="gam",formula=y~s(x,k=3),aes( col=dynamic_class3))
ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Tcell communication phenotypes ALLArms/Inflammatory signals from macrophages to T cells over time Cohort 2.png")

ggplot( macrophageTcelcomm[ARM!="A"][ReceptorPhenoCelltype=="CD8+ T cells"][]  , aes(y=FromMacrophages+startvalAv-startval, x=log(1+Day), fill=dynamic_class3, group=interaction(dynamic_class3) ))+
  theme_classic(base_size=18)+
  geom_violin(alpha=0.4,aes(group=interaction(dynamic_class3, Day)) ,scale="width" ,position=position_dodge() )+#geom_smooth(method="lm")+
  geom_point(position =position_jitterdodge(dodge.width=2.3,jitter.width=0.0)  )+
  theme(aspect.ratio=1)+scale_x_continuous(breaks=log(1+c(0,14,180)) ,labels=c(0,14,180))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant","Sensitive"))+
  scale_color_npg(name="Tumor response", labels=c("Resistant","Sensitive"))+
  labs(y="Inflammatory cytokine signaling \n from macrophages to cancer cells", x="Day")+
  theme(axis.title = element_blank(),axis.text = element_blank(), legend.position = "none")+
    geom_smooth(method="gam",formula=y~s(x,k=3),aes( col=dynamic_class3))
  
ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Tcell communication phenotypes ALLArms/BLANK Inflammatory signals from macrophages to T cells over time Cohort 2.png")




ggplot( macrophageTcelcomm[ARM!="A"][ReceptorPhenoCelltype=="CD8+ T cells"] [Patient.Study.ID!="001-167"], aes(y=FromMacrophages+startvalAv-startval, x=log(1+Day), fill=dynamic_class3, group=interaction(dynamic_class3) ))+
  theme_classic(base_size=18)+#facet_wrap(~ReceptorPhenoCelltype)+
  geom_violin(alpha=0.4,aes(group=interaction(dynamic_class3, Day) ),scale="width" ,position=position_dodge() )+#geom_smooth(method="lm")+
  geom_point(position =position_jitterdodge(dodge.width=2.3,jitter.width=0.0)  )+
  theme(aspect.ratio=1)+scale_x_continuous(breaks=log(1+c(0,14,180)) ,labels=c(0,14,180))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant","Sensitive"))+
  scale_color_npg(name="Tumor response", labels=c("Resistant","Sensitive"))+
  labs(y="Inflammatory cytokine signaling \n from macrophages to cancer cells", x="Day")+
  geom_smooth(method="gam",formula=y~s(x,k=3),aes( col=dynamic_class3))
ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Tcell communication phenotypes ALLArms/Inflammatory signals from macrophages to T cells over time Cohort 2 remove Menrighed samp.png")

ggplot( macrophageTcelcomm[ARM!="A"][ReceptorPhenoCelltype=="CD8+ T cells"][Patient.Study.ID!="001-167"]  , aes(y=FromMacrophages+startvalAv-startval, x=log(1+Day), fill=dynamic_class3, group=interaction(dynamic_class3) ))+
  theme_classic(base_size=18)+
  geom_violin(alpha=0.4,aes(group=interaction(dynamic_class3, Day)) ,scale="width" ,position=position_dodge() )+#geom_smooth(method="lm")+
  geom_point(position =position_jitterdodge(dodge.width=2.3,jitter.width=0.0)  )+
  theme(aspect.ratio=1)+scale_x_continuous(breaks=log(1+c(0,14,180)) ,labels=c(0,14,180))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant","Sensitive"))+
  scale_color_npg(name="Tumor response", labels=c("Resistant","Sensitive"))+
  labs(y="Inflammatory cytokine signaling \n from macrophages to cancer cells", x="Day")+
  theme(axis.title = element_blank(),axis.text = element_blank(), legend.position = "none")+
  geom_smooth(method="gam",formula=y~s(x,k=3),aes( col=dynamic_class3))

ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Tcell communication phenotypes ALLArms/BLANK Inflammatory signals from macrophages to T cells over time Cohort 2 remove Menrighed samp.png")
