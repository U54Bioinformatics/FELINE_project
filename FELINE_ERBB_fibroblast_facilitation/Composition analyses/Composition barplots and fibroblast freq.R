rm(list=ls())
require(data.table);require(dplyr);require(tidyr);require(ggplot2);require(boot);#require("compositions")
require(vegan);require(ggsci)#library(microbiome)#library(phyloseq)
require(mgcv);require(lme4);require(lmerTest);require(parallel);library(effects);require(umap);require(Rfast);require(ider);library("dendextend");library(ggdendro);require(viridis);require("Rdimtools")
library(caret);library(pROC)

load("/Users/jason/Jason Griffiths Dropbox/CancerMacrophageCPM_FELINE1/scRNAseq_FelineCohort1and2_CellMetaData.RData")
dd0<-CellMetaData

## Start by getting the M1-M2 phenotypes of myeloid cells
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
  data.table(Cohort="Discovery" ,    ddCohort1%>%select(Treatment,Cell.ID,Sample,Patient.Study.ID,Day,ARM,dynamic_class3,Celltype,Celltype_subtype,V1,V2) ),
  data.table(Cohort="Validation" ,  ddCohort2%>%select(Treatment,Cell.ID,Sample,Patient.Study.ID,Day,ARM,dynamic_class3,Celltype,Celltype_subtype,V1,V2))
)

# all cell type annotations
dd0 <-rbind(
  data.table(Cohort="Validation", metadd %>%select(Cell.ID,Celltype,Celltype_subtype,dynamic_class3,Patient.Study.ID,Day,ARM)) , 
  data.table(Cohort="Discovery",cohort1metadd%>%select(Cell.ID,Celltype,Celltype_subtype,dynamic_class3,Patient.Study.ID,Day,ARM) ))

dd0[,Treatmentlab:= "Combination ribociclib"]
dd0[ARM=="A",Treatmentlab:= "Letrozole alone"]
dd0[Cell.ID%in% joineddd[Celltype=="Macrophages"][Celltype_subtype=="Macrophages"][V2>0]$Cell.ID, Celltype_subtype:= "M2 macrophages" ]
dd0[Cell.ID%in% joineddd[Celltype=="Macrophages"][Celltype_subtype=="Macrophages"][V2<0]$Cell.ID, Celltype_subtype:= "M1 macrophages" ]
dd0[Celltype_subtype=="Endothelial cells",Celltype_subtype:="Vas-Endo"]

dd0$CellAnnot1<-dd0$Celltype
dd0$CellAnnot2<-dd0$Celltype_subtype
dd0[Celltype%in%c("Cancer cells","Normal epithelial cells"),CellAnnot1:="Epithelial cells"]

# Total cell count per sample= sampling effort
dd0[,totalCount:= length(Cell.ID) ,by=c("Patient.Study.ID","Day","Cohort","dynamic_class3","Treatmentlab")]

# count of each cell type's abundance and frequency
countTable0 <- data.table( dd0 %>% group_by(Patient.Study.ID,Day,Cohort,Treatmentlab,dynamic_class3,
                                            CellAnnot2,totalCount)%>% summarise(count=n() , frac=n()/unique(totalCount) ) ) 

#countTable0[,CellAnnot2Broader:= "Immune"]
#countTable0[CellAnnot2%in%c("Fibroblasts","Lym-Endo","Pericytes", "Vas-Endo" ,"Adipocytes") ,CellAnnot2Broader:= "Stromal"]
#countTable0[CellAnnot2%in%c("Cancer cells","Normal epithelial cells"),CellAnnot2Broader:= "Epithelial"]

countTable0[,CellAnnot2Broader:= CellAnnot2]
countTable0[CellAnnot2%in%c("CD4+ T cells","CD8+ T cells","M1 macrophages","M2 macrophages","Monocytes","NK cells","Plasma cells","Tregs","B cells","DC" ) ,CellAnnot2Broader:= "Immune"]
countTable0[CellAnnot2%in%c("Vas-Endo" , "Lym-Endo","Pericytes") ,CellAnnot2Broader:= "Endothelial/Pericytes"]
#countTable0[CellAnnot2%in%c("Adipocytes" , "Normal epithelial cells") ,CellAnnot2Broader:= "Other"]

GroupedcountTable0 <- data.table(countTable0%>%group_by(Patient.Study.ID,Day,Treatmentlab,CellAnnot2Broader,Cohort)%>%dplyr::summarise(frac=sum(frac)) )
GroupedcountTable0[, Sample:= paste(Patient.Study.ID,Day,sep="_")]
GroupedcountTable0[, FibroFrac:= sum(frac*(CellAnnot2Broader=="Fibroblasts")) ,by="Sample"]
GroupedcountTable0[, CancerFrac:= sum(frac*(CellAnnot2Broader=="Cancer cells")) ,by="Sample"]
#GroupedcountTable0$CellAnnot2Broader <- factor(GroupedcountTable0$CellAnnot2Broader, levels = c("Cancer cells","Fibroblasts" ,"Normal epithelial cells","Endothelial", "Adipocytes" ,"Pericytes","Immune" ))
GroupedcountTable0$CellAnnot2Broader <- factor(GroupedcountTable0$CellAnnot2Broader, levels = c("Fibroblasts" ,"Cancer cells","Normal epithelial cells","Endothelial/Pericytes" ,"Adipocytes","Immune" ))
#GroupedcountTable0$CellAnnot2Broader <- factor(GroupedcountTable0$CellAnnot2Broader, levels = c("Fibroblasts" ,"Cancer cells","Endothelial/Pericytes" ,"Immune","Other" ))
# spread data to construct dendrogram
#GroupedcountTable0wide<- data.table(spread( GroupedcountTable0[Day%in%c(0)] ,CellAnnot2Broader,frac, fill=0)  )
#dend <- as.dendrogram(hclust(dist(GroupedcountTable0wide[,-c(1:3)] )))
#GroupedcountTable0$Patient.Study.ID <- factor(GroupedcountTable0$Patient.Study.ID ,levels= unique(GroupedcountTable0wide[labels(dend),]$Patient.Study.ID) )
GroupedcountTable0$Sample <- factor(GroupedcountTable0$Sample ,levels= unique(GroupedcountTable0[order(Day, FibroFrac)]$Sample))

ggplot(GroupedcountTable0[][Day%in%c(0,14,180)],aes(y=frac,x=Sample,fill=CellAnnot2Broader) )+ geom_bar(stat="identity") +
  theme_classic(base_size=26)+
  facet_wrap(~paste0("Day ",Day),scales="free_x")+ scale_fill_aaas(name="") +
  theme(aspect.ratio=1,axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  labs(y="Tumor fraction", x="Tumor sample")

ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Fibroblast focused composition plot.pdf", height=6, width = 12)

GroupedcountTable0[, Timepoint:="Pre treatment"]
GroupedcountTable0[Day!=0, Timepoint:="Post treatment"]
GroupedcountTable0$Timepoint<- factor(GroupedcountTable0$Timepoint , levels=c("Pre treatment","Post treatment"))
ggplot(GroupedcountTable0[][Day%in%c(0, 180)],aes(y=100*frac,x=Sample,fill=CellAnnot2Broader) )+ geom_bar(stat="identity") +
  theme_classic(base_size=26)+
  facet_wrap(~Timepoint,scales="free_x")+ scale_fill_aaas(name="") +
  theme(aspect.ratio=1,axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  labs(y="Composition (%)", x="Patient tumor biopsy")
ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Fibroblast focused composition plotB.pdf", height=6, width = 12)


write.csv( countTable0 ,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/BothCohortTumorCellComposition.csv" )
write.csv( GroupedcountTable0 ,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/BothCohortTumorBroadCellComposition.csv" )

countTable0 <- data.table( read.csv( file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/BothCohortTumorCellComposition.csv" ))
countTable0[CellAnnot2Broader=="Normal epithelial cells",CellAnnot2Broader:="Diploid epithelial cells"]
countTable0[CellAnnot2=="Normal epithelial cells",CellAnnot2:="Diploid epithelial cells"]
countTable0[,CellAnnot2Medium:= CellAnnot2]
countTable0[CellAnnot2%in%c("CD4+ T cells","CD8+ T cells","NK cells","Tregs"),CellAnnot2Medium:= "T cells"]
countTable0[CellAnnot2%in%c("M1 macrophages","M2 macrophages","Monocytes","DC"),CellAnnot2Medium:= "Myeloid cells"]
countTable0[CellAnnot2%in%c("Plasma cells","B cells"),CellAnnot2Medium:= "B cells"]
countTable0[CellAnnot2%in%c("Vas-Endo" , "Lym-Endo"),CellAnnot2Medium:= "Endothelial cells"]
countTable0[CellAnnot2%in%c("Pericytes"),CellAnnot2Medium:= "Pericytes"]

countTable0[ ,BroaderN:= sum(count) , by="CellAnnot2Broader"]
countTable0[ ,MediumN:= sum(count) , by="CellAnnot2Medium"]
countTable0[ ,SpecificN:= sum(count) , by="CellAnnot2"]

ordFreq<-unique(countTable0%>%dplyr::select(CellAnnot2Medium, MediumN))[order(-MediumN)] #ordFreq$CellAnnot2Medium <- factor(ordFreq$CellAnnot2Medium  , levels=ordFreq$CellAnnot2Medium)
ordFreq$CellAnnot2Medium <- factor(ordFreq$CellAnnot2Medium  , levels=c("Cancer cells","Diploid epithelial cells",
                                                                        "Fibroblasts","Endothelial cells","Adipocytes","Pericytes",
                                                                         "Myeloid cells","T cells","B cells"))
ordFreq2<-unique(countTable0%>%dplyr::select(CellAnnot2Broader, BroaderN))[order(-BroaderN)] #ordFreq$CellAnnot2Medium <- factor(ordFreq$CellAnnot2Medium  , levels=ordFreq$CellAnnot2Medium)
ordFreq2$CellAnnot2Broader <- factor(ordFreq2$CellAnnot2Broader  , levels=c("Cancer cells","Diploid epithelial cells",
                                                                        "Fibroblasts","Endothelial/Pericytes","Adipocytes",
                                                                        "Immune"))

countsummary <- ordFreq[order(CellAnnot2Medium)] #unique(countTable0%>%dplyr::select(CellAnnot2,CellAnnot2Medium, CellAnnot2Broader,SpecificN, MediumN,BroaderN))
countsummary[,MediumFrac:=MediumN/sum(MediumN)]
# Compute the cumulative percentages (top of each rectangle)
countsummary$ymax = cumsum(countsummary$MediumFrac)
# Compute the bottom of each rectangle
countsummary$ymin = c(0, head(countsummary$ymax, n=-1))
# Compute a good label
countsummary$labelPosition <- (countsummary$ymax + countsummary$ymin) / 2
#countsummary$label <- paste0(countsummary$CellAnnot2Medium, "\n n= ", countsummary$MediumN)
countsummary$label <- paste0(countsummary$CellAnnot2Medium, " \n (", round(countsummary$MediumFrac*100,2), "%)")


countsummary2 <- ordFreq2[order(CellAnnot2Broader)] #unique(countTable0%>%dplyr::select(CellAnnot2,CellAnnot2Medium, CellAnnot2Broader,SpecificN, MediumN,BroaderN))
countsummary2[,BroaderFrac:=BroaderN/sum(BroaderN)]
# Compute the cumulative percentages (top of each rectangle)
countsummary2$ymax = cumsum(countsummary2$BroaderFrac)
# Compute the bottom of each rectangle
countsummary2$ymin = c(0, head(countsummary2$ymax, n=-1))
# Compute a good label
countsummary2$labelPosition <- (countsummary2$ymax + countsummary2$ymin) / 2
#countsummary$label <- paste0(countsummary$CellAnnot2Medium, "\n n= ", countsummary$MediumN)
countsummary2$label <- paste0(countsummary2$CellAnnot2Broader, " \n (", round(countsummary2$BroaderFrac*100,2), "%)")

#  plot
p1<-ggplot(countsummary, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=CellAnnot2Medium)) +
  geom_rect() +
  geom_label( x= #c(4,4,4,4,4,4, 3.66,3.33,3)
              c(4,3,4,3,4,3, 4,3,4)
                #2+2*(nrow(countsummary):1)/nrow(countsummary)
                , aes(y=labelPosition, label=label), size=4) +
 # scale_fill_brewer(palette=4) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")

p1+  geom_rect(data=countsummary2, aes(ymax=ymax, ymin=ymin, xmax=2, xmin=2.5, fill=CellAnnot2Broader)) 
 


ggplot(countsummary, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=CellAnnot2Medium)) +
  geom_rect() +
  geom_label( x= c(4,4,4,4,4,4, 3.66,3.33,3)
              #2+2*(nrow(countsummary):1)/nrow(countsummary)
              , aes(y=labelPosition, label=label), size=4) +
  # scale_fill_brewer(palette=4) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")

+
  geom_rect()







countsummary$SpecificN%>%sum()
FibroFracdd <- data.table( GroupedcountTable0%>%group_by(Patient.Study.ID,Day, Sample,Cohort)%>%dplyr::summarise(CancerFrac = mean(CancerFrac),
                                                                                                            FibroFrac=mean(FibroFrac)) )
ggplot(FibroFracdd[CancerFrac !=0 ][ FibroFrac!=0][Day %in% c(0, 14, 180) ] , aes(y= log(FibroFrac), x= log(1+Day), group= Day, fill= log(1+Day) , shape=Cohort ) ) + geom_violin() +
  geom_point(size=1.5) +  theme_classic(base_size=26)+
  # geom_path(aes(group=Patient.Study.ID),size=0.5) +
  facet_wrap(~Cohort) +
  scale_x_continuous(breaks=log(1+c(0,14,180)), labels=c(0,14,180))+
  scale_fill_continuous(name="Day",breaks=log(1+c(0,14,180)), labels=c(0,14,180))+
  scale_y_continuous(breaks=log(c(0.001,0.01,0.1,1)), labels=c(0.001,0.01,0.1,1))+
  labs(y="Fibroblast fraction", x="Day")+theme(aspect.ratio=1)+guides(fill = "none")
ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Fibroblast frequency over time.pdf", height=6, width = 12)


ggplot(FibroFracdd[CancerFrac !=0 ][ FibroFrac!=0][Day %in% c(0, 14, 180) ] , aes(y= (FibroFrac), x= log(1+Day), group= Day, fill= Day , shape=Cohort ) ) + geom_violin() +
  geom_point(size=1.5) +  theme_classic(base_size=26)+
  # geom_path(aes(group=Patient.Study.ID),size=0.5) +
  facet_wrap(~Cohort) +
  scale_x_continuous(breaks=log(1+c(0,14,180)), labels=c(0,14,180))+
  #scale_y_continuous(breaks=(c(0.001,0.01,0.1,1)), labels=c(0.001,0.01,0.1,1))+
  labs(y="Fibroblast fraction", x="Day")+theme(aspect.ratio=1)


plotfib<-FibroFracdd[CancerFrac !=0 ][ FibroFrac!=0]
plotfib[ ,initialFibroFrac:= sum( (Day==0)*FibroFrac) ,by=Patient.Study.ID]


ggplot(plotfib, aes(y= log(FibroFrac/initialFibroFrac), x= log(1+Day), group= Day, fill= Day , shape=Cohort ) ) + geom_violin() +
  geom_point(size=1.5) +
  geom_path(aes(group=Patient.Study.ID),size=0.5)

wideFibrofrac<- data.table( spread(FibroFracdd[Day %in% c(0,  180) ][CancerFrac !=0 & FibroFrac!=0]%>%dplyr::select(-c(CancerFrac,Sample)),Day,FibroFrac ) , fill=0 )
setnames(wideFibrofrac , old=c("0","180"), new=c("F0","F180"))


ggplot(wideFibrofrac , aes(y= log(F180/F0), 1  ) ) + geom_violin() +
  geom_point(aes( shape= Cohort) ,size=1.5)+geom_hline(yintercept=0) 


summary( lmer( log((FibroFrac)/CancerFrac) ~ Day  +(1|Patient.Study.ID) ,  data= FibroFracdd[Cohort=="Validation"][CancerFrac !=0& FibroFrac!=0 ][Day %in% c(0,180)] ) )
summary( lmer( log((FibroFrac)/CancerFrac) ~ Day  +(1|Patient.Study.ID) ,  data= FibroFracdd[Cohort!="Validation"][CancerFrac !=0& FibroFrac!=0 ][Day %in% c(0,180)] ) )
summary( lmer( log((FibroFrac)/CancerFrac) ~ Day *Cohort +(1|Patient.Study.ID) ,  data= FibroFracdd[CancerFrac !=0& FibroFrac!=0 ][Day %in% c(0,180)] ) )

summary( lm( log((1e-5+FibroFrac)/CancerFrac) ~ as.factor(Day) , family= "binomial", data= FibroFracdd[CancerFrac !=0 ][Day %in% c(0,180)] ) )
summary( lm( log((1e-5+FibroFrac)/CancerFrac) ~ as.factor(Day) , family= "binomial", data= FibroFracdd[CancerFrac !=0&FibroFrac>0 ][Day %in% c(0,180)] ) )
summary( lm( log((1e-5+FibroFrac)/CancerFrac) ~ as.factor(Day) , family= "binomial", data= FibroFracdd[CancerFrac !=0&FibroFrac!=0 ][Day %in% c(0,180)] ) )

summary( glm( (FibroFrac/CancerFrac) ~as.factor(Day) ,family="binomial", data=FibroFracdd[Day%in%c(0,180)] ) )
  +
  facet_wrap(~Day,scales="free")


ggplot(GroupedcountTable0[][Day%in%c(0,180)],aes(y=frac,x=Sample,fill=CellAnnot2Broader) )+geom_bar(stat="identity") +
  facet_wrap(Cohort~CellAnnot2Broader,scales="free")+scale_fill_aaas()

GroupedcountTable0%>%group_by(CellAnnot2Broader)%>%dplyr::summarise(m=median(frac))
