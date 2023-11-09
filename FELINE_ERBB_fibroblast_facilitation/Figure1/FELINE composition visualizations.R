rm(list=ls())
require(data.table);require(dplyr);require(tidyr);require(ggplot2);require(boot);#require("compositions")
require(vegan);require(ggsci)#library(microbiome)#library(phyloseq)
require(mgcv);require(lme4);require(lmerTest);require(parallel);library(effects);require(umap);require(Rfast);require(ider);library("dendextend");library(ggdendro);require(viridis);require("Rdimtools")
library(caret);library(pROC)

load("/Users/jason/Jason Griffiths Dropbox/CancerMacrophageCPM_FELINE1/scRNAseq_FelineCohort1and2_CellMetaData.RData")
dd0 <- CellMetaData
dd0[Celltype=="Normal epithelial cells",Celltype:="Diploid epithelial cells"]
dd0[Celltype_subtype=="Normal epithelial cells",Celltype_subtype:="Diploid epithelial cells"]

# Make new catagorical labels
dd0$CellAnnot1<-dd0$Celltype
dd0$CellAnnot2<-dd0$Celltype_subtype
dd0[Celltype%in%c("Cancer cells","Diploid epithelial cells"),CellAnnot1:="Epithelial cells"]


# Total cell count per sample= sampling effort
dd0[,totalCount:= length(Cell.ID) ,by=c("Patient.Study.ID","Day","Cohort","dynamic_class3","Treatmentlab")]

# count of each cell type's abundance and frequency
countTable0 <- data.table( dd0 %>% group_by(Patient.Study.ID,Day,Cohort,Treatmentlab,dynamic_class3,
                                            CellAnnot2,totalCount)%>% summarise(count=n() , frac=n()/unique(totalCount) ) ) 

countTable0[,CellAnnot2Broader:= CellAnnot2]
countTable0[CellAnnot2%in%c("CD4+ T cells","CD8+ T cells","M1 macrophages","M2 macrophages","Monocytes","NK cells","Plasma cells","Tregs","B cells","DC" ) ,CellAnnot2Broader:= "Immune"]
countTable0[CellAnnot2%in%c("Vas-Endo" , "Lym-Endo","Pericytes") ,CellAnnot2Broader:= "Endothelial/Pericytes"]

GroupedcountTable0 <- data.table(countTable0%>%group_by(Patient.Study.ID,Day,Treatmentlab,CellAnnot2Broader,Cohort)%>%dplyr::summarise(frac=sum(frac)) )
GroupedcountTable0[, Sample:= paste(Patient.Study.ID,Day,sep="_")]
GroupedcountTable0[, FibroFrac:= sum(frac*(CellAnnot2Broader=="Fibroblasts")) ,by="Sample"]
GroupedcountTable0[, CancerFrac:= sum(frac*(CellAnnot2Broader=="Cancer cells")) ,by="Sample"]
GroupedcountTable0$CellAnnot2Broader <- factor(GroupedcountTable0$CellAnnot2Broader, levels = c("Fibroblasts" ,"Cancer cells","Diploid epithelial cells","Endothelial/Pericytes" ,"Adipocytes","Immune" ))
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


colPal<-pal_aaas()(6)
colPal[6]<-pal_jco()(2)[2]
ggplot(GroupedcountTable0[][Day%in%c(0, 180)],aes(y=100*frac,x=Sample,fill=CellAnnot2Broader) )+ geom_bar(stat="identity") +
  theme_classic(base_size=26)+
  facet_wrap(~Timepoint,scales="free_x")+ 
  scale_fill_manual(name="", values=colPal) +
  theme(aspect.ratio=1,axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  labs(y="Composition (%)", x="Patient tumor biopsy")
ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Fibroblast focused composition plotC.pdf", height=6, width = 12)




write.csv( countTable0 ,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/BothCohortTumorCellComposition.csv" )
write.csv( GroupedcountTable0 ,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/BothCohortTumorBroadCellComposition.csv" )

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

countsummary <- ordFreq[order(CellAnnot2Medium)] 
countsummary[,MediumFrac:=MediumN/sum(MediumN)]
# Compute the cumulative percentages (top of each rectangle)
countsummary$ymax = cumsum(countsummary$MediumFrac)
# Compute the bottom of each rectangle
countsummary$ymin = c(0, head(countsummary$ymax, n=-1))
# Compute a good label
countsummary$labelPosition <- (countsummary$ymax + countsummary$ymin) / 2
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

p1

ggplot(countsummary2, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=CellAnnot2Broader)) +
  geom_rect() +
  geom_label( x= c(3.85,3.85,3.85,3.85,3.85,3.85)            , aes(y=labelPosition, label=label), size=4) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void(base_size=32) +
  theme(legend.position = "none",aspect.ratio=1)+
  scale_fill_manual(values=pal_aaas()(6)[c(2,3,1,4,5,6)])


ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Pie chart composition plot.pdf", height=5.5, width = 5.5)

colPal<-pal_aaas()(6)
colPal[6]<-pal_jco()(2)[2]

ggplot(countsummary2, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=CellAnnot2Broader)) +
  geom_rect() +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void(base_size=32) +
  theme(legend.position = "none",aspect.ratio=1)+
  scale_fill_manual(values=colPal[c(2,3,1,4,5,6)])


ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/BLANK Pie chart composition plot.pdf", height=5.5, width = 5.5)



ggplot(countsummary2, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=CellAnnot2Broader)) +
  geom_rect() +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void(base_size=32) +
  theme(legend.position = "none",aspect.ratio=1)+
  scale_fill_manual(values=pal_aaas()(6)[c(2,3,1,4,5,6)])
ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Pie chart composition plot no labels.pdf", height=5.5, width = 5.5)





