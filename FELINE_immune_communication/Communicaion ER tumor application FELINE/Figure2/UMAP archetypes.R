rm(list=ls())
require(data.table);require(dplyr);require(tidyr);require(ggplot2);require(ggsci)
require(boot);require("compositions");require(vegan)
require(mclust);require(mgcv);require(lme4);require(lmerTest);require(parallel)
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

abundDD <- countTable0 %>% dplyr::select(-count)%>%spread(CellAnnot2,frac,fill=0)
names(abundDD) <- gsub("\\+", "_", names(abundDD))
names(abundDD) <- gsub("\\-", "_", names(abundDD))
names(abundDD) <- gsub(" ", "_", names(abundDD))

### Full cohort all samples
# fraction table spread with cell types across columns #[Timepoint==0]
fracTable1 <- spread( countTable0[][]%>%dplyr::select(-count), CellAnnot2,frac,fill=0)
fracTable1$rowID <- 1:nrow(fracTable1)

# fraction matrix without metadata. name rows by ID to retain match to metadata
fracTable2 <- as.matrix( fracTable1%>%dplyr::select(-c("Patient.Study.ID", "Day", "Cohort","Treatmentlab" ,"dynamic_class3", "totalCount","rowID")) )
rownames(fracTable2)<- fracTable1$rowID

rowannot <- data.frame(fracTable1 %>% dplyr::select(rowID,Cohort, Day,dynamic_class3, Treatmentlab,totalCount))
rownames(rowannot) <- rowannot$rowID; rowannot$rowID <- NULL
pheatmap::pheatmap( logit( 1e-10+ fracTable2*0.9999),annotation_row=rowannot)

distmat <- vegdist(fracTable2, method = "bray")
pheatmap::pheatmap(as.matrix(distmat),annotation_col=rowannot)

######

fracTable1 <- spread( countTable0[
  Cohort=="Discovery"
  ]%>%dplyr::select(-count), CellAnnot2,frac,fill=0)
fracTable1$rowID <- 1:nrow(fracTable1)

# fraction matrix without metadata. name rows by ID to retain match to metadata
fracTable2 <- as.matrix( fracTable1 %>% dplyr::select(-c("Patient.Study.ID", "Day", "Cohort","Treatmentlab" ,"dynamic_class3", "totalCount","rowID")) )
rownames(fracTable2)<- fracTable1$rowID

# adjust range
eps <- 1e-5
adjfracTable2 <- eps+ fracTable2 *(1 - 2*eps)

set.seed(123) 
Duinput <- logit(adjfracTable2[,!colnames(adjfracTable2)%in%c("Cancer cells","Normal_epithelial_cells")]) 
u1 <- umap::umap(logit(adjfracTable2[,!colnames(adjfracTable2)%in%c("Cancer cells","Normal_epithelial_cells")]) ,n_components=2, n_neighbors=8)#, n_neighbors=20)

shannonDiv <- vegan::diversity(adjfracTable2,index = "shannon")
#Pielou’s evenness
Evenness <- shannonDiv/log(specnumber(fracTable2))
gmm1 <- Mclust(u1$layout,G=3,modelNames="EEV")
# gather output of NMDS
umapout <- data.table( cbind( fracTable1,  u1$layout,shannon=shannonDiv ,Evenness), archetype= gmm1$classification ) #,princomp(scores(comm.bc.mds))$scores
names(umapout) <- gsub("\\+", "_", names(umapout))
names(umapout) <- gsub("\\-", "_", names(umapout))
names(umapout) <- gsub(" ", "_", names(umapout))
names(umapout)
cors<-data.table( t(cor(umapout%>%dplyr::select(V1,V2), umapout%>%dplyr::select(Adipocytes:Vas_Endo))) , keep.rownames =T)
cors[order(-abs(V2))]
cors[order(-abs(V1))]

table( umapout %>% dplyr::select(archetype,dynamic_class3))%>%chisq.test()

###
Vdd0 <-rbind(
  data.table(Cohort="Validation", metadd[Quality=="HQ"] %>%dplyr::select(Cell.ID,Celltype,Celltype_subtype,dynamic_class3,Patient.Study.ID,Day,ARM)) , 
  data.table(Cohort="Discovery",cohort1metadd%>%dplyr::select(Cell.ID,Celltype,Celltype_subtype,dynamic_class3,Patient.Study.ID,Day,ARM) ))

Vdd0[,Treatmentlab:= "Combination ribociclib"]
Vdd0[ARM=="A",Treatmentlab:= "Letrozole alone"]
Vdd0[Cell.ID%in% joineddd[Celltype=="Macrophages"][Celltype_subtype=="Macrophages"][V2>0]$Cell.ID, Celltype_subtype:= "M2 macrophages" ]
Vdd0[Cell.ID%in% joineddd[Celltype=="Macrophages"][Celltype_subtype=="Macrophages"][V2<0]$Cell.ID, Celltype_subtype:= "M1 macrophages" ]
Vdd0[Celltype_subtype=="Endothelial cells",Celltype_subtype:="Vas-Endo"]
Vdd0$CellAnnot1<-Vdd0$Celltype
Vdd0$CellAnnot2<-Vdd0$Celltype_subtype
Vdd0[Celltype%in%c("Cancer cells","Normal epithelial cells"),CellAnnot1:="Epithelial cells"]
# Total cell count per sample= sampling effort
Vdd0[,totalCount:= length(Cell.ID) ,by=c("Patient.Study.ID","Day","Cohort","dynamic_class3","Treatmentlab")]



VcountTable0 <- data.table( Vdd0 %>% group_by(Patient.Study.ID,Day,Cohort,Treatmentlab,dynamic_class3,
                                            CellAnnot2,totalCount)%>% summarise(count=n() , frac=n()/unique(totalCount) ) ) 

VfracTable1 <- spread( VcountTable0[      Cohort=="Validation" ]%>%dplyr::select(-count), CellAnnot2,frac,fill=0)
VfracTable1$rowID <- 1:nrow(VfracTable1)
VfracTable1$Tregs <- 0

# fraction matrix without metadata. name rows by ID to retain match to metadata
VfracTable2 <- as.matrix( VfracTable1 %>% dplyr::select(-c("Patient.Study.ID", "Day", "Cohort","Treatmentlab" ,"dynamic_class3", "totalCount","rowID")) )
rownames(VfracTable2)<- VfracTable1$rowID
VadjfracTable2 <- eps + VfracTable2 *(1 - 2*eps)

Vinput <-  logit( VadjfracTable2[,!colnames(VadjfracTable2)%in%c("Cancer cells","Normal_epithelial_cells")] )[,u1$data%>%colnames()]
Vu1 <- predict(u1,data=Vinput)

VshannonDiv <- vegan::diversity(VadjfracTable2,index = "shannon")
#Pielou’s evenness
VEvenness <- VshannonDiv/log(specnumber(VfracTable2))
# gather output of NMDS
Vumapout <- data.table( cbind( VfracTable1,  Vu1,shannon=VshannonDiv ,Evenness=VEvenness), archetype= predict(gmm1,Vu1)$classification ) #,princomp(scores(comm.bc.mds))$scores
names(Vumapout) <- gsub("\\+", "_", names(Vumapout))
names(Vumapout) <- gsub("\\-", "_", names(Vumapout))
names(Vumapout) <- gsub(" ", "_", names(Vumapout))

####


allumapout <- rbind( umapout,Vumapout)
#save( allumapout , file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/Figure2/Discovery and Validation UMAP composition archetypes.RData")

sz<-5
ggplot(allumapout[], aes(x=V1, y=V2, col=(dynamic_class3) , shape=as.factor(Day) ) ) + geom_point(size=sz) + theme_classic(base_size=26) +
  scale_shape_manual(name="Day",values=c(8,18,19))+
  scale_color_npg(name="Tumor response", labels=c("Resistant","Sensitive")) +
  labs(y="TME composition axis 2",x="TME composition axis 1" ) +theme(aspect.ratio=1) +
  facet_wrap(~Cohort,scales="free")

ggplot(allumapout[Day==180], aes(x=V1, y=V2,col=(dynamic_class3) , shape=as.factor(Day) )) + geom_point(size=sz) + theme_classic(base_size=26) +
  scale_shape_manual(name="Day",values=c(8,18,19))+
  scale_color_npg(name="Tumor response", labels=c("Resistant","Sensitive")) +
  labs(y="TME composition axis 2",x="TME composition axis 1" ) +theme(aspect.ratio=1) +
  facet_wrap(Treatmentlab~Cohort)


ggplot(allumapout[Day%in%c(0,14,180)], aes(x=V1, y=V2,col=as.factor(archetype) , shape=Cohort )) + geom_point(size=sz) + theme_classic(base_size=26) +
  labs(y="TME composition axis 2",x="TME composition axis 1" ) +theme(aspect.ratio=1) 

ggplot(allumapout[Day%in%c(0,14,180)], aes(x=V1, y=V2,col=(dynamic_class3) , shape=as.factor(Day) )) + geom_point(size=sz) + theme_classic(base_size=26) +
  scale_shape_manual(name="Day",values=c(8,18,19))+
  scale_color_npg(name="Tumor response", labels=c("Resistant","Sensitive")) +
  labs(y="TME composition axis 2",x="TME composition axis 1" ) +theme(aspect.ratio=1) +
  facet_wrap(~Cohort)



paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/Discovery and Validation "
sz<-5
# Composition archetypes
ggplot(allumapout, aes(x=V1, y=V2,col=as.factor(archetype),shape=Cohort )) + geom_point(size=sz) + theme_classic(base_size=18)+theme(aspect.ratio=1)+
  scale_color_uchicago(name="Composition \n archetype", labels=c("A","B","C","D")) + #jama# jco
  labs(y="TME composition axis 2",x="TME composition axis 1" )
ggsave(paste0(paperfile,"Ribo and Letrozole Composition archetypes with archetypes annotation.png"),height=10,width=10, dpi=320)

ggplot(allumapout, aes(x=V1, y=V2,col=as.factor(archetype),shape=Cohort )) + geom_point(size=1.5*sz) + theme_classic(base_size=18)+theme(aspect.ratio=1)+
  scale_color_uchicago(name="Composition \n archetype", labels=c("A","B","C","D")) + #jama# jco
  labs(y="TME composition axis 2",x="TME composition axis 1" )+
  theme(axis.title=element_blank(),  axis.text=element_blank(),strip.text = element_blank(),legend.title=element_blank(),legend.text=element_blank() ,legend.position = "none")
ggsave(paste0(paperfile,"BLANK Ribo and Letrozole Composition archetypes with archetypes annotation.png"),height=5,width=5, dpi=320)





# Composition corelates to outcome
ggplot(allumapout[], aes(x=V1, y=V2,col=(dynamic_class3) , shape=as.factor(Day) )) + geom_point(size=sz) + theme_classic(base_size=26) +
  scale_shape_manual(name="Day",values=c(8,18,19))+
  scale_color_npg(name="Tumor response", labels=c("Resistant","Sensitive")) +
  labs(y="TME composition axis 2",x="TME composition axis 1" ) +theme(aspect.ratio=1) #+

ggsave(paste0(paperfile,"Ribo and Letrozole updated Composition archetypes with Tumor response.png"),height=10,width=10,dpi=320)


ggplot(allumapout, aes(x=V1, y=V2,col=(dynamic_class3) , shape=as.factor(Day) )) + geom_point(size=1.5*sz) + theme_classic(base_size=18) +
  scale_shape_manual(name="Day",values=c(8,18,19))+
  scale_color_npg(name="Tumor response", labels=c("Resistant","Sensitive")) +
  labs(y="TME composition axis 2",x="TME composition axis 1" ) +theme(aspect.ratio=1)+
  theme(axis.title=element_blank(),  axis.text=element_blank(),strip.text = element_blank(),legend.title=element_blank(),legend.text=element_blank() ,legend.position = "none")

ggsave(paste0(paperfile,"BLANK Ribo and Letrozole updated BOTH Composition archetypes with Tumor response.png"),height=5,width=5,dpi=320)


ggplot(allumapout[Treatmentlab=="Combination ribociclib"], aes(x=V1, y=V2,col=(dynamic_class3) , shape=as.factor(Day) )) + geom_point(size=1.5*sz) + theme_classic(base_size=18) +
  scale_shape_manual(name="Day",values=c(8,18,19))+
  scale_color_npg(name="Tumor response", labels=c("Resistant","Sensitive")) +
  labs(y="TME composition axis 2",x="TME composition axis 1" ) +theme(aspect.ratio=1)+ facet_wrap(~Treatmentlab, ncol=1)+
  theme(axis.title=element_blank(),  axis.text=element_blank(),strip.text = element_blank(),legend.title=element_blank(),legend.text=element_blank() ,legend.position = "none")

ggsave(paste0(paperfile,"BLANK Ribo and Letrozole updated RIBO Composition archetypes with Tumor response.png"),height=5,width=5,dpi=320)


ggplot(allumapout[Treatmentlab!="Combination ribociclib"], aes(x=V1, y=V2,col=(dynamic_class3) , shape=as.factor(Day) )) + geom_point(size=1.5*sz) + theme_classic(base_size=18) +
  scale_shape_manual(name="Day",values=c(8,18,19))+
  scale_color_npg(name="Tumor response", labels=c("Resistant","Sensitive")) +
  labs(y="TME composition axis 2",x="TME composition axis 1" ) +theme(aspect.ratio=1)+ facet_wrap(~Treatmentlab, ncol=1)+
  theme(axis.title=element_blank(),  axis.text=element_blank(),strip.text = element_blank(),legend.title=element_blank(),legend.text=element_blank() ,legend.position = "none")

ggsave(paste0(paperfile,"BLANK Ribo and Letrozole updated LETROZOLE Composition archetypes with Tumor response.png"),height=5,width=5,dpi=320)



ggplot(allumapout, aes(x=V1, y=V2,col=(dynamic_class3) , shape=as.factor(Day) )) + geom_point(size=sz) + theme_classic(base_size=28) +
  scale_shape_manual(name="Day",values=c(8,18,19))+
  scale_color_npg(name="Tumor response", labels=c("Resistant","Sensitive")) +
  labs(y="TME composition axis 2",x="TME composition axis 1" ) +theme(aspect.ratio=1)+facet_wrap(~Treatmentlab, ncol=1)
ggsave(paste0(paperfile,"Ribo and Letrozole updated LETROZOLE Composition archetypes with Tumor response by treatment.png"),height=5,width=5,dpi=320)


# Composition archetypes drivers
ggplot(allumapout, aes(x=V1, y=V2,col=logit(eps+ Cancer_cells*(1 - (100)*eps)) , shape=as.factor(archetype)) )+ geom_point(size=sz) + theme_classic(base_size=18)+
  theme(aspect.ratio=1)+labs(y="TME composition axis 2",x="TME composition axis 1" )+
  scale_color_viridis_c(name="Cancer cell \n frequency",option='C',breaks=logit(eps+c(0.005,0.05,0.5)*(1 - 100*eps)), labels=c(0.005,0.05,0.5)) +
  scale_shape_discrete(name="Composition \n archetype", labels= c( "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched", "Cancer \n dominated")) 
ggsave(paste0(paperfile,"Ribo and Letrozole Composition archetypes with Cancer_cells frequency.png"),height=10,width=10)

ggplot(allumapout, aes(x=V1, y=V2,col=logit(eps+ Cancer_cells*(1 - (100)*eps)) , shape=as.factor(archetype)) )+ geom_point(size=1.5*sz) + theme_classic(base_size=18)+
  theme(aspect.ratio=1)+labs(y="TME composition axis 2",x="TME composition axis 1" )+
  scale_color_viridis_c(name="Cancer cell \n frequency",option='C',breaks=logit(eps+c(0.005,0.05,0.5)*(1 - 100*eps)), labels=c(0.005,0.05,0.5)) +
  scale_shape_discrete(name="Composition \n archetype", labels= c( "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched", "Cancer \n dominated"),guide=F) +
  theme(axis.title=element_blank(),  axis.text=element_blank(),strip.text = element_blank(),legend.title=element_blank(),legend.text=element_blank() )
ggsave(paste0(paperfile,"BLANK Ribo and Letrozole Composition archetypes with Cancer_cells frequency.png"),height=10,width=10)



ggplot(allumapout, aes(x=V1, y=V2, col=logit(eps+(B_cells+Plasma_cells+CD4__T_cells+CD8__T_cells+NK_cells+Tregs+Monocytes+DC+M1_macrophages+M2_macrophages)*(1 - eps)) , shape=as.factor(archetype)) )+ geom_point(size=sz) + theme_classic(base_size=18)+
  theme(aspect.ratio=1)+labs(y="TME composition axis 2",x="TME composition axis 1" ) +
  scale_color_viridis_c(name="Immune cell \n frequency",option='B',breaks=logit(eps+c(0.001,0.01,0.1)*(1 - eps)), labels=c(0.001,0.01,0.1) ) +
  scale_shape_discrete(name="Composition \n archetype", labels= c( "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched", "Cancer \n dominated")) 
ggsave(paste0(paperfile,"Ribo and Letrozole Composition archetypes with Immune_cells frequency.png"),height=10,width=10)

ggplot(allumapout, aes(x=V1, y=V2, col=logit(eps+(B_cells+Plasma_cells+CD4__T_cells+CD8__T_cells+NK_cells+Tregs+Monocytes+DC+M1_macrophages+M2_macrophages)*(1 - eps)) , shape=as.factor(archetype)) )+
  geom_point(size=1.5*sz) + theme_classic(base_size=18)+
  theme(aspect.ratio=1)+labs(y="TME composition axis 2",x="TME composition axis 1" ) +
  scale_color_viridis_c(name="Immune cell \n frequency",option='B',breaks=logit(eps+c(0.001,0.01,0.1)*(1 - eps)), labels=c(0.001,0.01,0.1) ) +
  scale_shape_discrete(name="Composition \n archetype", labels= c( "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched", "Cancer \n dominated"),guide=F)  +
  theme(axis.title=element_blank(),  axis.text=element_blank(),strip.text = element_blank(),legend.title=element_blank(),legend.text=element_blank() )+theme(legend.position="none")
ggsave(paste0(paperfile,"BLANK Ribo and Letrozole Composition archetypes with Immune_cells frequency.png"),height=5,width=5)




ggplot(allumapout, aes(x=V1, y=V2, col=logit(eps+(CD4__T_cells + CD8__T_cells + NK_cells +Tregs)*(1 - eps)) , shape=as.factor(archetype)) )+ geom_point(size=sz) + theme_classic(base_size=18)+
  theme(aspect.ratio=1)+labs(y="TME composition axis 2",x="TME composition axis 1" ) +
  scale_color_viridis_c(name="T cell \n frequency",option='B',breaks=logit(eps+c(0.001,0.01,0.1)*(1 - eps)), labels=c(0.001,0.01,0.1)) +
  scale_shape_discrete(name="Composition \n archetype", labels= c( "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched", "Cancer \n dominated")) 
ggsave(paste0(paperfile,"Ribo and Letrozole Composition archetypes with T_cells frequency.png"),height=10,width=10)


ggplot(allumapout, aes(x=V1, y=V2, col=logit(eps+(CD4__T_cells + CD8__T_cells + NK_cells +Tregs)*(1 - eps)) , shape=as.factor(archetype)) )+ geom_point(size=1.5*sz) + theme_classic(base_size=18)+
  theme(aspect.ratio=1)+labs(y="TME composition axis 2",x="TME composition axis 1" ) +
  scale_color_viridis_c(name="T cell \n frequency",option='B',breaks=logit(eps+c(0.001,0.01,0.1)*(1 - eps)), labels=c(0.001,0.01,0.1)) +
  scale_shape_discrete(name="Composition \n archetype", labels= c( "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched", "Cancer \n dominated"),guide=F) +
  theme(axis.title=element_blank(),  axis.text=element_blank(),strip.text = element_blank(),legend.title=element_blank(),legend.text=element_blank() )+theme(legend.position="none")
ggsave(paste0(paperfile,"BLANK Ribo and Letrozole Composition archetypes with T_cells frequency.png"),height=5,width=5)



ggplot(allumapout, aes(x=V1, y=V2, col=logit(eps+(Adipocytes+Fibroblasts+Lym_Endo+Pericytes+Vas_Endo)*(1 - eps)), shape=as.factor(archetype) ) )+ geom_point(size=sz) + theme_classic(base_size=18)+
  theme(aspect.ratio=1)+labs(y="TME composition axis 2",x="TME composition axis 1" ) +
  scale_color_viridis_c(name="Stromal cell \n frequency",option='D',breaks=logit(eps+c(0.001,0.01,0.1)*(1 - eps)), labels=c(0.001,0.01,0.1)) +
  scale_shape_discrete(name="Composition \n archetype", labels= c( "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched", "Cancer \n dominated")) 
ggsave(paste0(paperfile,"Ribo and Letrozole Composition archetypes with Stromal_cells frequency.png"),height=10,width=10)

ggplot(allumapout, aes(x=V1, y=V2, col=logit(eps+(Adipocytes+Fibroblasts+Lym_Endo+Pericytes+Vas_Endo)*(1 - eps)), shape=as.factor(archetype) ) )+ geom_point(size=1.5*sz) + theme_classic(base_size=18)+
  theme(aspect.ratio=1)+labs(y="TME composition axis 2",x="TME composition axis 1" ) +
  scale_color_viridis_c(name="Stromal cell \n frequency",option='D',breaks=logit(eps+c(0.001,0.01,0.1)*(1 - eps)), labels=c(0.001,0.01,0.1)) +
  scale_shape_discrete(name="Composition \n archetype", labels= c( "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched", "Cancer \n dominated"),guide=F) +
  theme(axis.title=element_blank(),  axis.text=element_blank(),strip.text = element_blank(),legend.title=element_blank(),legend.text=element_blank() )+theme(legend.position="none")
ggsave(paste0(paperfile,"BLANK Ribo and Letrozole Composition archetypes with Stromal_cells frequency.png"),height=10,width=10)




ggplot(allumapout, aes(x=V1, y=V2, shape=as.factor(archetype), col=shannon ) )+ geom_point(size=sz)+ theme_classic(base_size=18) +  #,fill=as.factor(archetype)
  theme(aspect.ratio=1)+labs(y="TME composition axis 2",x="TME composition axis 1" ) +
  scale_color_viridis_c(name="Shannon \n Diversity")+
  scale_shape_discrete(name="Composition \n archetype", labels= c( "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched", "Cancer \n dominated")) +
  labs(y="TME composition axis 2",x="TME composition axis 1")+ theme(aspect.ratio=1)#+
ggsave(paste0(paperfile,"Ribo and Letrozole Composition archetypes with Shannon Diversity.png"),height=10,width=10, dpi=320)


ggplot(allumapout, aes(x=V1, y=V2, shape=as.factor(archetype), col=shannon ) )+ geom_point(size=1.5*sz)+ theme_classic(base_size=18) +  #,fill=as.factor(archetype)
  theme(aspect.ratio=1)+labs(y="TME composition axis 2",x="TME composition axis 1" ) +
  scale_color_viridis_c(name="Shannon \n Diversity")+
  scale_shape_discrete(name="Composition \n archetype", labels= c( "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched", "Cancer \n dominated"),guide=F) +
  labs(y="TME composition axis 2",x="TME composition axis 1")+ theme(aspect.ratio=1)+
  theme(axis.title=element_blank(),  axis.text=element_blank(),strip.text = element_blank(),legend.title=element_blank(),legend.text=element_blank() )+theme(legend.position="none")
ggsave(paste0(paperfile,"BLANK Ribo and Letrozole Composition archetypes with Shannon Diversity.png"),height=5,width=5, dpi=320)

allumapout[,isarchetype1:=F]
allumapout[archetype==1,isarchetype1:=T]
allumapout[,isarchetype3:=F]
allumapout[archetype==3,isarchetype3:=T]

summary( glm(I(dynamic_class3=="Response")~ as.factor(isarchetype1), data=allumapout[Day==0][Treatmentlab=="Combination ribociclib"], family="binomial") )
summary( glm(I(dynamic_class3=="Response")~ as.factor(isarchetype1), data=allumapout[Day==14][Treatmentlab=="Combination ribociclib"], family="binomial") )
summary( glm(I(dynamic_class3=="Response")~ as.factor(isarchetype1), data=allumapout[Day==180][Treatmentlab=="Combination ribociclib"], family="binomial") )

summary( glm(I(dynamic_class3=="Response")~ as.factor(isarchetype1), data=allumapout[Day==0][Treatmentlab=="Letrozole alone"], family="binomial") )
summary( glm(I(dynamic_class3=="Response")~ as.factor(isarchetype1), data=allumapout[Day==14][Treatmentlab=="Letrozole alone"], family="binomial") )
summary( glm(I(dynamic_class3=="Response")~ as.factor(isarchetype1), data=allumapout[Day==180][Treatmentlab=="Letrozole alone"], family="binomial") )


summary( glm(cbind( I(totalCount*( B_cells+Plasma_cells+CD4__T_cells+CD8__T_cells+NK_cells+Tregs+Monocytes+DC+M1_macrophages+M2_macrophages ) ),totalCount )~ as.factor(isarchetype1), family="binomial", data=allumapout[][])  )
summary( glm(cbind( I(totalCount*( CD4__T_cells+CD8__T_cells+NK_cells+Tregs ) ),totalCount )~ as.factor(isarchetype1), family="binomial", data=allumapout[][])  )
summary( lm(log(shannon)~ as.factor(isarchetype1), data=allumapout[][]) )
summary( lm(log(shannon)~ as.factor(isarchetype3), data=allumapout[][]) )


m1<- glm(I(dynamic_class3=="Response") ~ 0 + as.factor(archetype), data=allumapout[], family="binomial")
summary(m1)
unique(predict(m1, type="response"))
# save data to make plot 
save(umapout,file="~/Dropbox/Cancer_pheno_evo/data/FELINE2/Composition/Composition archetype data.RData")

m180<- glm(I(dynamic_class3=="Response")~0+as.factor(archetype), data=allumapout[Day==180], family="binomial")
summary(m180)
predict(m180,newdata=data.table(archetype=as.factor(1:3)), type="response")%>%unique()

m0<- glm(I(dynamic_class3=="Response")~0+as.factor(archetype), data=allumapout[Day==0], family="binomial")
summary(m0)
predict(m0,newdata=data.table(archetype=as.factor(1:3)), type="response")%>%unique()

transitionsInput<-na.omit(allumapout[Day!=14]%>% 
                            mutate(DayClass = paste0("Day",Day),
                                   archetypeClass = paste0("Archetype",archetype)) %>%
                            dplyr::select(Treatmentlab,dynamic_class3,Patient.Study.ID,DayClass,archetypeClass)%>% spread(DayClass,archetypeClass))

transitionsInput[dynamic_class3=="Response"]$Day180%>%table()
transitionsInput[dynamic_class3!="Response"]$Day180%>%table()
transitionsprobs <- data.table( data.table( table( transitionsInput%>%dplyr::select(Day0,Day180,Treatmentlab,dynamic_class3 )) )  %>%group_by(Treatmentlab,dynamic_class3,Day0) %>%mutate(prop=N/sum(N)) )


umapout[Treatmentlab=="Letrozole alone"]%>%dplyr::select(dynamic_class3,archetype)%>%table()
transitionsprobs[,From:="Immune hot \n and diverse"]
transitionsprobs[Day0=="Archetype2",From:="Fibroblast/endothlial \n enriched"]
transitionsprobs[Day0=="Archetype3",From:="Cancer \n dominated"]
transitionsprobs[,To:="Immune hot \n and diverse"]
transitionsprobs[Day180=="Archetype2",To:="Fibroblast/endothlial \n enriched"]
transitionsprobs[Day180=="Archetype3",To:="Cancer \n dominated"]
transitionsprobs[prop==0,prop:=NA]
ggplot(transitionsprobs,aes(x=From,y=To, fill=(prop) ))+geom_tile() + facet_wrap(Treatmentlab~dynamic_class3)+
  scale_fill_viridis_c(option="C",na.value="white")

ggplot(transitionsprobs,aes(x=From,y=To, alpha=(prop) ))+geom_tile() + theme_classic(base_size=12)+ facet_wrap(Treatmentlab~dynamic_class3) + 
  scale_alpha_continuous(name="Transition probability",na.value=0)+theme(aspect.ratio=1,axis.text.x=element_text(angle=90,vjust = 0.5))
ggsave(paste0(paperfile,"Ribo and Letrozole T cell Composition State transition probability.png"),height=10,width=10)

ggplot(transitionsprobs,aes(x=From,y=To, alpha=(prop) ))+geom_tile() + theme_classic(base_size=12)+ facet_wrap(Treatmentlab~dynamic_class3) + 
  scale_alpha_continuous(name="Transition probability",na.value=0)+theme(aspect.ratio=1,axis.text.x=element_text(angle=90,vjust = 0.5))+
  theme(axis.title=element_blank(),  axis.text=element_blank(),strip.text = element_blank(),legend.title=element_blank(),legend.text=element_blank() )
#ggsave(paste0(paperfile,"BLANK Ribo and Letrozole Composition State transition probability.png"),height=10,width=10)



fracTable1 <- spread( rbind(countTable0[Cohort=="Discovery"],VcountTable0[Cohort!="Discovery"]) %>%dplyr::select(-count), CellAnnot2,frac,fill=0)
fracTable1$rowID <- 1:nrow(fracTable1)
fracTable1<- merge(fracTable1, allumapout%>%dplyr::select(Patient.Study.ID,Day,archetype,M1_macrophages,shannon) , by=c("Patient.Study.ID","Day"))[order(archetype,-M1_macrophages,-shannon)]
fracTable1$archetype <- as.factor(fracTable1$archetype)

# fraction matrix without metadata. name rows by ID to retain match to metadata
fracTable2 <- as.matrix( fracTable1%>%dplyr::select(-c("Patient.Study.ID", "Day", "Cohort","Treatmentlab" ,"dynamic_class3", "totalCount","rowID","archetype","shannon","M1_macrophages")) )
rownames(fracTable2)<- fracTable1$rowID
# adjust range
eps <- 1e-5
adjfracTable2<- eps+ fracTable2 *(1 - 2*eps)
rowannot <- data.frame(fracTable1 %>% dplyr::select(rowID, dynamic_class3,Day, Treatmentlab,archetype,Cohort))
colnames(rowannot)[colnames(rowannot)=="archetype"]<-"Archetype"
rownames(rowannot) <- rowannot$rowID; rowannot$rowID <- NULL
rowannot$Day <- as.factor(rowannot$Day)
rowannot$dynamic_class3[rowannot$dynamic_class3=="Response"]<- "Sensitive"
rowannot$dynamic_class3[rowannot$dynamic_class3=="Non-response"]<- "Resistant"

rowannot$Archetype<- as.character(rowannot$Archetype)
rowannot$Archetype[rowannot$Archetype=="1"]<- "Immune hot and diverse"
rowannot$Archetype[rowannot$Archetype=="2"]<- "Fibroblast/endothlial enriched"
rowannot$Archetype[rowannot$Archetype=="3"]<- "Cancer dominated"

rowannot$Archetype<- factor(rowannot$Archetype, levels=c("Immune hot and diverse","Fibroblast/endothlial enriched","Cancer dominated"))

setnames(rowannot,old=c("dynamic_class3","Treatmentlab"), new=c("Tumor response","Treatment"))

annotation_col = list( "Tumor response" = rev(pal_npg("nrc")(2) )   ,
                       Treatment=rev(pal_jco("default")(2))   ,
                       Day=c("khaki4","purple", "grey" ),
                       Cohort=c("black","grey"),
                       Archetype=c(pal_uchicago("default")(3)))

names(annotation_col[[1]] ) <- c("Sensitive","Resistant")
names(annotation_col[[2]] ) <- c("Combination ribociclib","Letrozole alone")
names(annotation_col[[3]] ) <- c("0","14","180")
names(annotation_col[[4]] ) <- c("Discovery", "Validation")
names(annotation_col[[5]] ) <- c("Immune hot and diverse","Fibroblast/endothlial enriched","Cancer dominated")


cellord <-c("Cancer cells","DC","M1 macrophages","Monocytes","B cells","Plasma cells","M2 macrophages","Tregs","NK cells","CD4+ T cells","CD8+ T cells",
            "Normal epithelial cells","Fibroblasts"  ,"Pericytes","Vas-Endo","Adipocytes" ,"Lym-Endo" )

y<- adjfracTable2[,cellord]

colnames(y)[colnames(y)=="Normal epithelial cells"]<-"Diploid epithelial cells"
colnames(y)[colnames(y)=="DC"]<-"Dendritic cells"
colnames(y)[colnames(y)=="Vas-Endo"]<-"Vascular Endodothelial cells"
colnames(y)[colnames(y)=="Lym-Endo"]<-"Lymphatic Endodothelial cells"

cellord <-c("Cancer cells","Dendritic cells","M1 macrophages","Monocytes","B cells","Plasma cells","M2 macrophages","Tregs","NK cells","CD4+ T cells","CD8+ T cells",
            "Diploid epithelial cells","Fibroblasts"  ,"Pericytes","Vascular Endodothelial cells","Adipocytes" ,"Lymphatic Endodothelial cells" )
p1<-pheatmap::pheatmap( scale(logit(y  )) ,annotation_row=rowannot,cutree_rows = 4,cutree_cols = 3,
                        annotation_colors = annotation_col,
                        border_color=NA,
                        treeheight_row=10,treeheight_col=6,
                        cluster_rows = F,cluster_cols = F,
                        labels_row=rep("",200),
                        cellheight=3,cellwidth=12)

#dims 9x9: filename:Discovery and Validation composition heatmap Ribo and Letrozole alone.pdf

p1<-pheatmap::pheatmap( scale(logit( adjfracTable2)),annotation_row=rowannot[,1:3],cutree_rows = 4,cutree_cols = 3,
                        annotation_colors = annotation_col,
                        border_color=NA,
                        treeheight_row=10,treeheight_col=6,
                        #cellwidth=10,cellheight=4,
                        #annotation_legend=F,
                        labels_row=rep("",100),
                        cellheight=3,cellwidth=12)




fracTable1 <- spread( countTable0[
  Cohort=="Discovery"
  ]%>%dplyr::select(-count), CellAnnot2,frac,fill=0)
fracTable1$rowID <- 1:nrow(fracTable1)
fracTable1<- merge(fracTable1, umapout%>%dplyr::select(Patient.Study.ID,Day,archetype) , by=c("Patient.Study.ID","Day"))
fracTable1$archetype <- as.factor(fracTable1$archetype)

# fraction matrix without metadata. name rows by ID to retain match to metadata
fracTable2 <- as.matrix( fracTable1%>%dplyr::select(-c("Patient.Study.ID", "Day", "Cohort","Treatmentlab" ,"dynamic_class3", "totalCount","rowID","archetype")) )
rownames(fracTable2)<- fracTable1$rowID
# adjust range
eps <- 1e-5
adjfracTable2<- eps+ fracTable2 *(1 - 2*eps)
rowannot <- data.frame(fracTable1 %>% dplyr::select(rowID, dynamic_class3,Day, Treatmentlab,archetype))
rownames(rowannot) <- rowannot$rowID; rowannot$rowID <- NULL
rowannot$Day <- as.factor(rowannot$Day)
rowannot$dynamic_class3[rowannot$dynamic_class3=="Response"]<- "Sensitive"
rowannot$dynamic_class3[rowannot$dynamic_class3=="Non-response"]<- "Resistant"
setnames(rowannot,old=c("dynamic_class3","Treatmentlab"), new=c("Tumor response","Treatment"))

annotation_col = list( "Tumor response" = rev(pal_npg("nrc")(2) )   ,Treatment=rev(pal_jco("default")(2))   ,Day=c("khaki4","purple", "grey" ))
names(annotation_col[[1]] ) <- c("Sensitive","Resistant")
names(annotation_col[[2]] ) <- c("Combination ribociclib","Letrozole alone")
names(annotation_col[[3]] ) <- c("0","14","180")

#pheatmap::pheatmap((pp),cluster_rows=F,annotation_row = annotation_row, annotation_colors = annotation_col,scale="none",show_rownames=F,annotation_legend=F,border_color=NA)


p1<-pheatmap::pheatmap( scale(logit( adjfracTable2)),annotation_row=rowannot,cutree_rows = 4,cutree_cols = 3,
                        annotation_colors = annotation_col,
                        border_color=NA,
                        treeheight_row=10,treeheight_col=6,
                        labels_row=rep("",100),
                        cellheight=3,cellwidth=12)


