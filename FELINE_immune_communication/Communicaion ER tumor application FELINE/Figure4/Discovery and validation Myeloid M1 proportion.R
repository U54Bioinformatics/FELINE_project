rm(list=ls())
require(data.table);require(dplyr);require(tidyr);require(ggplot2)
require(boot);require("compositions");require(vegan);require(ggsci)
require(mgcv);require(lme4);require(lmerTest);require(parallel)
library(effects);require(umap);require(Rfast);require(ider);library("dendextend");library(ggdendro);require(ggsci);require(viridis);require("Rdimtools")
library(caret)
library(pROC)

load( file=paste0("~/Dropbox/Cancer_pheno_evo/data/FELINE2/PhenotypesAllArmsCohort2/CPMPhenotpyeLandscape_C2ProjectC2RevisednewMacrophagesDC.RData"))#C1umap,u_dat,DAY,cell_types_all,ARMS,Subtype,
load( file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Cohort2Metadata/Cohort2Metadata.RData") #metadd,cohort1metadd,annotation.file,compdataLU

dd1 <- u_dat %>% dplyr::select(Cell.ID:file_string,V1,V2)

Cohort1Extract <-function(){
  load(file= "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/UMAP genes macrophage only ALLARMS/UMAP genes macrophage only ALLARMS.RData")
  dd1<-u_dat%>%dplyr::select(Cell.ID:PhenoCelltype,V1,V2)
  return(dd1)
}
ddCohort1 <- Cohort1Extract()
ddCohort2 <- dd1
ddCohort2[,Treatment:="letrozole + ribo"]
ddCohort2[ARM=="A",Treatment:="letrozole"]

joineddd <- rbind(
  data.table(Cohort="Discovery" ,    ddCohort1%>%select(Treatment,Cell.ID,Sample,Patient.Study.ID,Day,ARM,dynamic_class3,Celltype,Celltype_subtype,V1,V2) ),
  data.table(Cohort="Validation" ,  ddCohort2%>%select(Treatment,Cell.ID,Sample,Patient.Study.ID,Day,ARM,dynamic_class3,Celltype,Celltype_subtype,V1,V2))
)

dd0 <-rbind(
  data.table(Cohort="Validation", metadd %>%select(Cell.ID,Celltype,Celltype_subtype,dynamic_class3,Patient.Study.ID,Day,ARM)) , 
  data.table(Cohort="Discovery",cohort1metadd%>%select(Cell.ID,Celltype,Celltype_subtype,dynamic_class3,Patient.Study.ID,Day,ARM) ))

dd0[,Treatmentlab:= "Combination ribociclib"]
dd0[ARM=="A",Treatmentlab:= "Letrozole alone"]
dd0[Cell.ID%in% joineddd[Celltype=="Macrophages"][Celltype_subtype=="Macrophages"][V2>0]$Cell.ID, Celltype_subtype:= "M2 macrophages" ]
dd0[Cell.ID%in% joineddd[Celltype=="Macrophages"][Celltype_subtype=="Macrophages"][V2<0]$Cell.ID, Celltype_subtype:= "M1 macrophages" ]
dd0[Celltype=="Macrophages"]$Celltype_subtype%>%table()
dd0[Celltype_subtype=="Endothelial cells",Celltype_subtype:="Vas-Endo"]

dd0$CellAnnot1<-dd0$Celltype
dd0$CellAnnot2<-dd0$Celltype_subtype
dd0[Celltype%in%c("Cancer cells","Normal epithelial cells"),CellAnnot1:="Epithelial cells"]
dd0[Celltype%in%c("Cancer cells","Normal epithelial cells"),CellAnnot2:="Epithelial cells"]

# Total cell count per sample= sampling effort
dd0[,totalCount:= length(Cell.ID) ,by=c("Patient.Study.ID","Day","Cohort","dynamic_class3","Treatmentlab")]

# count of each cell type's abundance and frequency
countTable0 <- data.table( dd0 %>% group_by(Patient.Study.ID, Day, Cohort, Treatmentlab, dynamic_class3,
                                            CellAnnot2,totalCount) %>% summarise(count=n() , frac=n()/unique(totalCount) ) ) 

countTable0wideM1M2<-countTable0[CellAnnot2 %in% unique(dd0[CellAnnot1=="Macrophages"]$CellAnnot2) ]%>% select(-count)%>%spread(CellAnnot2,frac,fill=0) 
names(countTable0wideM1M2) <- gsub(" ","_",names(countTable0wideM1M2))
countTable0wideM1M2[,M1M2ratio:= M1_macrophages/(  M1_macrophages + M2_macrophages ) ]
countTable0wideM1M2[,Myeloidtot := ( DC + M1_macrophages + M2_macrophages + Monocytes)  ]

countTable0wideM1M2[,Success:=1]
countTable0wideM1M2[dynamic_class3=="Non-response",Success:=0]
countTable0wideM1M2<-countTable0wideM1M2[is.finite(M1M2ratio)]
countTable0wideM1M2[,x:=M1M2ratio ]
countTable0wideM1M2[,y:=log(Myeloidtot) ]


ggplot(countTable0wideM1M2[(Myeloidtot*totalCount)>20], 
       aes(M1M2ratio, x=log(1+Day), fill=dynamic_class3,col=dynamic_class3)) + theme_classic(base_size=16)+
  geom_violin(aes(group=interaction(Day,dynamic_class3)),alpha=0.6,position=position_dodge(width=0.9))+
  geom_point(aes(shape=Cohort,group=interaction(Day,dynamic_class3)),size=2,position=position_dodge(width=0.9)) + 
  facet_grid(~Treatmentlab) + theme(aspect.ratio=1,axis.text.x=element_text(angle=90)) +
  geom_smooth(method="gam", formula=y~s(x,k=3),se=F,size=2)+
  scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"),guide="none")+
  scale_fill_npg(name="Tumor response", labels=c("Resistant", "Sensitive"),guide="none")+
  scale_x_continuous(name="Day", breaks= log(1+c(0,14,180)),labels=c(0,14,180))

ggplot(countTable0wideM1M2[Day!=180][(Myeloidtot*totalCount)>20], 
       aes(y= (M1M2ratio), x= dynamic_class3, fill= dynamic_class3 )) + theme_classic(base_size= 26)+
  geom_boxplot()+
  geom_boxplot(outlier.colour=NA, position= position_dodge() ,col="black")+
  stat_boxplot(geom="errorbar",position=position_dodge(1.75),width=0.5)+#geom_smooth(method="lm")+
  facet_wrap(Cohort~ Treatmentlab) + theme(aspect.ratio= 1) +
  scale_fill_npg(name= "Tumor response", labels= c("Resistant", "Sensitive"))+
  labs(y="Proportion myeloid cells in M1 state",x="Tumor response") + 
  geom_point(aes(shape=Cohort),size=2.5, col="black")+theme(legend.position="none")+
  scale_x_discrete(labels=c("Resistant", "Sensitive"))

ggplot(countTable0wideM1M2[Day!=180][(Myeloidtot*totalCount)>20], 
       aes(y= (M1M2ratio), x= dynamic_class3, fill= dynamic_class3 )) + theme_classic(base_size= 28)+
  geom_boxplot()+
  geom_boxplot(outlier.colour=NA, position= position_dodge() ,col="black")+
  stat_boxplot(geom="errorbar",position=position_dodge(1.75),width=0.5)+
  facet_wrap(~ Treatmentlab,ncol=1) + theme(aspect.ratio= 1) +
  scale_fill_npg(name= "Tumor response", labels= c("Resistant", "Sensitive"))+
  labs(y="Proportion myeloid cells in M1 state",x="Tumor response") + 
  geom_point(aes(shape=Cohort),size=2.5, col="black")+theme(legend.position="none")+
  scale_x_discrete(labels=c("Resistant", "Sensitive"))

paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"

ggsave(paste0(paperfile,"Discovery and Validation Ribociclib and Letrozole Early treatment M1 abundance high in response to ribo not letrozole boxplot.png"),height=10,width=10, dpi=320)


ggplot(countTable0wideM1M2[Day!=180], 
       aes(y= (M1M2ratio), x= dynamic_class3, fill= dynamic_class3 )) + theme_classic(base_size= 28)+
  geom_boxplot()+
  geom_boxplot(outlier.colour=NA, position= position_dodge() ,col="black")+
  stat_boxplot(geom="errorbar",position=position_dodge(1.75),width=0.5)+
  facet_wrap(~ Treatmentlab,ncol=1) + theme(aspect.ratio= 1) +
  scale_fill_npg(name= "Tumor response", labels= c("Resistant", "Sensitive"))+
  labs(y="Proportion myeloid cells in M1 state",x="Tumor response") + 
  geom_point(aes(shape=Cohort),size=2.5, col="black")+theme(legend.position="none")+
  scale_x_discrete(labels=c("Resistant", "Sensitive"))

paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"

ggsave(paste0(paperfile,"Discovery and Validation Ribociclib and Letrozole Early treatment M1 abundance high in response to ribo not letrozole boxplotB.png"),height=10,width=10, dpi=320)


glmmodC<-lm( M1M2ratio~ dynamic_class3*Day*Cohort, data=countTable0wideM1M2[Day!=180][Treatmentlab=="Combination ribociclib"])
glmmodL <-lm( M1M2ratio ~ dynamic_class3*Day*Cohort, data=countTable0wideM1M2[Day!=180][Treatmentlab=="Letrozole alone"] )
summary(glmmodC)
summary(glmmodL)

resdat<- countTable0wideM1M2

#save(resdat, file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/Figure4/Discovery and Validation Myeloid M1 proportion early treatment.RData")
load(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/Figure4/Discovery and Validation Myeloid M1 proportion early treatment.RData")
#resdat




