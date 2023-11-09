rm(list=ls())
require(data.table);require(dplyr);require(tidyr);require(ggplot2);require(ggsci)
require(boot);require("compositions");require(vegan)
require(mgcv);require(lme4);require(lmerTest);require(parallel)
library(effects);require(umap);require(Rfast);require(ider);library("dendextend");library(ggdendro);require(ggsci);require(viridis);require("Rdimtools")
library(caret);library(pROC)


# macrophage annotations
load( file=paste0("~/Dropbox/Cancer_pheno_evo/data/FELINE2/PhenotypesAllArmsCohort2/CPMPhenotpyeLandscape_C2ProjectC2RevisednewMacrophagesDC.RData"))#C1umap,u_dat,DAY,cell_types_all,ARMS,Subtype,
load( file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Cohort2Metadata/Cohort2Metadata.RData") #metadd,cohort1metadd,annotation.file,compdataLU

dd1 <- u_dat %>% dplyr::select(Cell.ID:file_string,V1,V2)

Cohort1Extract <- function(){
  load(file= "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/UMAP genes macrophage only ALLARMS/UMAP genes macrophage only ALLARMS.RData")
  dd1<-u_dat%>%dplyr::select(Cell.ID:PhenoCelltype,V1,V2)
  return(dd1)
}
ddCohort1 <- Cohort1Extract()

ddCohort2 <- dd1
ddCohort2[,Treatment:="letrozole + ribo"]
ddCohort2[ARM=="A",Treatment:="letrozole"]


joineddd <- rbind(
  data.table(Cohort="Discovery" ,    ddCohort1%>%dplyr::select(Treatment,Cell.ID,Sample,Patient.Study.ID,Day,ARM,dynamic_class3,Celltype,Celltype_subtype,V1,V2) ),
  data.table(Cohort="Validation" ,  ddCohort2%>%dplyr::select(Treatment,Cell.ID,Sample,Patient.Study.ID,Day,ARM,dynamic_class3,Celltype,Celltype_subtype,V1,V2))
)

# all cell type annotations
dd0 <-rbind(
  data.table(Cohort="Validation", metadd %>%dplyr::select(Cell.ID,Celltype,Celltype_subtype,dynamic_class3,Patient.Study.ID,Day,ARM)) , 
  data.table(Cohort="Discovery",cohort1metadd%>%dplyr::select(Cell.ID,Celltype,Celltype_subtype,dynamic_class3,Patient.Study.ID,Day,ARM) ))

dd0[,Treatmentlab:= "Combination ribociclib"]
dd0[ARM=="A",Treatmentlab:= "Letrozole alone"]
dd0[Cell.ID%in% joineddd[Celltype=="Macrophages"][Celltype_subtype=="Macrophages"][V2>0]$Cell.ID, Celltype_subtype:= "M2 macrophages" ]
dd0[Cell.ID%in% joineddd[Celltype=="Macrophages"][Celltype_subtype=="Macrophages"][V2<0]$Cell.ID, Celltype_subtype:= "M1 macrophages" ]
dd0[Celltype=="Macrophages"]$Celltype_subtype%>%table()
dd0[Celltype_subtype=="Endothelial cells",Celltype_subtype:="Vas-Endo"]

dd0$CellAnnot1<-dd0$Celltype
dd0$CellAnnot2<-dd0$Celltype_subtype
dd0[Celltype%in%c("Cancer cells","Normal epithelial cells"),CellAnnot1:="Epithelial cells"]


# Total cell count per sample= sampling effort
dd0[,totalCount:= length(Cell.ID) ,by=c("Patient.Study.ID","Day","Cohort","dynamic_class3","Treatmentlab")]

# count of each cell type's abundance and frequency
countTable0 <- data.table( dd0 %>% group_by(Patient.Study.ID,Day,Cohort,Treatmentlab,dynamic_class3,
                                            CellAnnot2,totalCount)%>% summarise(count=n() , frac=n()/unique(totalCount) ) ) 
countTableFeline1 <- countTable0[Cohort=="Discovery"]
#write.csv( countTableFeline1  , file="~/Dropbox/Cancer_pheno_evo/data/FELINE2/Composition/Feline_Cohort1_ERpositiveBC_Tumorcelltypecomposition_classified.csv" )

countTable0[CellAnnot2%in%c("Tregs","CD4+ T_cells", "CD8+ T cells","NK cells")]$count%>%sum()
abundDD <- countTable0 %>% dplyr::select(-count)%>%spread(CellAnnot2,frac,fill=0)
names(abundDD) <- gsub("\\+", "_", names(abundDD))
names(abundDD) <- gsub("\\-", "_", names(abundDD))
names(abundDD) <- gsub(" ", "_", names(abundDD))

eps <-1e-5
ggplot(abundDD[], aes(y=logit(eps+Tregs+CD4__T_cells+CD8__T_cells+NK_cells), x=log(1+Day),group=interaction(Day,dynamic_class3), fill=dynamic_class3 ,shape= Cohort) )+  theme_classic(base_size=26) +  #,fill=as.factor(archetype)
  geom_boxplot(col="black") + geom_point(size=3,position = position_dodge(width=2)) + facet_wrap(~Treatmentlab, ncol=1) +
  scale_fill_npg(name="Tumor response", labels=c("Resistant","Sensitive"))  +
  theme(aspect.ratio=1)+
  labs(y="T cell frequency", x= "Day")+
  scale_y_continuous(breaks=logit(eps+c(0.001,0.01,0.1) ), labels=c(0.001,0.01,0.1) ) +
  scale_x_continuous(breaks=log(1+c(0,14,180) ), labels=c(0,14,180) )
#ggsave(paste0(paperfile,"Ribo and Letrozole T cell Composition by response over time version2 Discovery and Validation.png"),height=10,width=10)

#save(abundDD , file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/Figure5/Discovery and Validation T cell frequency.RData")

abundDD
abundDD[,success:= round((Tregs+CD4__T_cells+CD8__T_cells+NK_cells)*totalCount)]
abundDD[,fail:= totalCount-success]


summary( glm(cbind(success,fail)   ~0+dynamic_class3*as.factor(Day) + Cohort , family="binomial", abundDD[Treatmentlab=="Combination ribociclib"]) )
summary( glm(cbind(success,fail)   ~0+dynamic_class3*as.factor(Day) + Cohort, family="binomial", abundDD[Treatmentlab=="Letrozole alone"]) )

summary( glm(cbind(success,fail)   ~dynamic_class3*as.factor(Day) + Cohort , family="binomial", abundDD[Treatmentlab=="Combination ribociclib"]) )

summary( glm(cbind(success,fail)   ~dynamic_class3*as.factor(Day) + Cohort, family="binomial", abundDD[Treatmentlab=="Letrozole alone"]) )


