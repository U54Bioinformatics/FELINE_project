rm(list=ls())
require(data.table);require(dplyr);require(tidyr);require(ggplot2);require(boot);require("compositions");require(vegan);require(ggsci);require(mclust)



## Start by getting the M1-M2 phenotypes of myeloid cells
require(mgcv);require(data.table);require(dplyr);require(ggplot2);require(tidyr);require(lme4);require(lmerTest);require(parallel)
library(effects);require(umap);require(Rfast);require(ider);library("dendextend");library(ggdendro);require(ggsci);require(viridis);require("Rdimtools")
library(caret)
library(pROC)


# macrophage annotations
load( file=paste0("~/Dropbox/Cancer_pheno_evo/data/FELINE2/PhenotypesAllArmsCohort2/CPMPhenotpyeLandscape_C2ProjectC2RevisednewMacrophagesDC.RData"))#C1umap,u_dat,DAY,cell_types_all,ARMS,Subtype,
load( file="~/Dropbox/Cancer_pheno_evo/data/FELINE2/Cohort2Metadata/Cohort2Metadata.RData") #metadd,cohort1metadd,annotation.file,compdataLU

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

# all cell type annotations
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


# Total cell count per sample= sampling effort
dd0[,totalCount:= length(Cell.ID) ,by=c("Patient.Study.ID","Day","Cohort","dynamic_class3","Treatmentlab")]

# count of each cell type's abundance and frequency
countTable0 <- data.table( dd0 %>% group_by(Patient.Study.ID,Day,Cohort,Treatmentlab,dynamic_class3,
                                            CellAnnot2,totalCount)%>% summarise(count=n() , frac=n()/unique(totalCount) ) ) 

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

fracTable1 <- spread( countTable0[
 Cohort=="Discovery"
  ]%>%dplyr::select(-count), CellAnnot2,frac,fill=0)
fracTable1$rowID <- 1:nrow(fracTable1)




# fraction matrix without metadata. name rows by ID to retain match to metadata
fracTable2 <- as.matrix( fracTable1%>%dplyr::select(-c("Patient.Study.ID", "Day", "Cohort","Treatmentlab" ,"dynamic_class3", "totalCount","rowID")) )
rownames(fracTable2)<- fracTable1$rowID
# adjust range
eps <- 1e-5
adjfracTable2<- eps+ fracTable2 *(1 - 2*eps)
rowannot <- data.frame(fracTable1 %>% dplyr::select(rowID, dynamic_class3,Day, Treatmentlab))
rownames(rowannot) <- rowannot$rowID; rowannot$rowID <- NULL
rowannot$Day <- as.factor(rowannot$Day)
setnames(rowannot,old=c("dynamic_class3","Treatmentlab"), new=c("Tumor response","Treatment"))


annotation_col = list( "Tumor response" = rev(pal_npg("nrc")(2) )   ,Treatment=rev(pal_jco("default")(2))   ,Day=c("khaki4","purple", "grey" ))
names(annotation_col[[1]] ) <- c("Response","Non-response")
names(annotation_col[[2]] ) <- c("Combination ribociclib","Letrozole alone")
names(annotation_col[[3]] ) <- c("0","14","180")

p1<-pheatmap::pheatmap( scale(logit( adjfracTable2)),annotation_row=rowannot,cutree_rows = 4,cutree_cols = 3,
                    annotation_colors = annotation_col,
                    border_color=NA,
                    treeheight_row=10,treeheight_col=6,
                    #cellwidth=10,cellheight=4,
                    #annotation_legend=F,
                    labels_row=rep("",100),
                    cellheight=3,cellwidth=12)
# save 6x6 : ~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Composition/heatmap of composition across treatments outcomes and days.pdf
save_pheatmap <- function(x, filename, width=480, height=480) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  png(filename,width = width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}


paperfile<- "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
save_pheatmap(p1,filename= paste0(paperfile,"Ribo and Letrozole heatmap of composition across treatments outcomes and days.png"),height=480,width=480)

blankrowannot<- rowannot
names( blankrowannot)<- c(" ","  ","   ")

annotation_col = list( "Tumor response" = rev(pal_npg("nrc")(2) )   ,Day=c("khaki4","purple", "grey" ) ,Treatment=rev(pal_jco("default")(2))   )
names( annotation_col)<- names( blankrowannot)
names(annotation_col[[1]] ) <- c("Response","Non-response")
names(annotation_col[[3]] ) <- c("Combination ribociclib","Letrozole alone")
names(annotation_col[[2]] ) <- c("0","14","180")

p1<-pheatmap::pheatmap( scale(logit( adjfracTable2)),annotation_row=blankrowannot,cutree_rows = 4,cutree_cols = 3,
                        annotation_colors = annotation_col,
                        border_color=NA,
                        treeheight_row=10,treeheight_col=6,
                        #cellwidth=10,cellheight=4,
                        labels_col=rep("",100),
                        annotation_legend=F,
                        legend_breaks = -4:2,
                        legend_labels = rep("",7),
                        labels_row=rep("",100),
                        cellheight=3,cellwidth=12)
paperfile<- "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
save_pheatmap(p1,filename= paste0(paperfile,"BLANK Ribo and Letrozole heatmap of composition across treatments outcomes and days.png"),height=480,width=480)


distmat <- vegdist(fracTable2, method = "bray")
pheatmap::pheatmap(as.matrix(distmat),annotation_col=rowannot)


set.seed(123) 
u1 <- umap::umap(logit(adjfracTable2[,!colnames(adjfracTable2)%in%c("Cancer cells","Normal_epithelial_cells")]) ,n_components=2, n_neighbors=8)#, n_neighbors=20)
shannonDiv <- diversity(adjfracTable2,index = "shannon")
#Pielou’s evenness
Evenness <- shannonDiv/log(specnumber(fracTable2))
gmm1 <- Mclust(u1$layout,G=3,modelNames="EEV")
# gathr output of NMDS
umapout <- data.table( cbind( fracTable1,  u1$layout,shannon=shannonDiv ,Evenness), archetype= gmm1$classification ) #,princomp(scores(comm.bc.mds))$scores
names(umapout) <- gsub("\\+", "_", names(umapout))
names(umapout) <- gsub("\\-", "_", names(umapout))
names(umapout) <- gsub(" ", "_", names(umapout))
names(umapout)
cors<-data.table( t(cor(umapout%>%select(V1,V2), umapout%>%select(Adipocytes:Vas_Endo))) , keep.rownames =T)
cors[order(-abs(V2))]
cors[order(-abs(V1))]

table( umapout%>%select(archetype,dynamic_class3))%>%chisq.test()
sz<-6

# Composition corelates to outcome
ggplot(umapout[], aes(x=V1, y=V2,col=(dynamic_class3) , shape=as.factor(Day) )) + geom_point(size=sz) + theme_classic(base_size=18) +
  scale_shape_manual(name="Day",values=c(8,18,19))+
  scale_color_npg(name="Tumor response", labels=c("Resistant","Sensitive")) +
  labs(y="TME composition axis 2",x="TME composition axis 1" ) +theme(aspect.ratio=1)

ggsave("~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Composition/Composition archetypes with Tumor response.png")
paperfile<- "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
ggsave(paste0(paperfile,"Ribo and Letrozole Composition archetypes with Tumor response.png"),height=10,width=10)

ggplot(umapout[], aes(x=V1, y=V2,col=(dynamic_class3) , shape=as.factor(Day) )) + geom_point(size=1.5*sz) + theme_classic(base_size=18) +
  scale_shape_manual(name="Day",values=c(8,18,19))+
  scale_color_npg(name="Tumor response", labels=c("Resistant","Sensitive")) +
  labs(y="TME composition axis 2",x="TME composition axis 1" ) +theme(aspect.ratio=1)+
  theme(axis.title=element_blank(),  axis.text=element_blank(),strip.text = element_blank(),legend.title=element_blank(),legend.text=element_blank() ,legend.position = "none")

ggsave(paste0(paperfile,"BLANK Ribo and Letrozole Composition archetypes with Tumor response.png"),height=10,width=10)




ggplot(umapout, aes(x=V1, y=V2,col=(dynamic_class3) , shape=as.factor(Day) )) + geom_point(size=sz) + theme_classic(base_size=18) +
  scale_shape_manual(name="Day",values=c(8,18,19))+
  scale_color_npg(name="Tumor response", labels=c("Resistant","Sensitive")) +
  labs(y="TME composition axis 2",x="TME composition axis 1" ) +theme(aspect.ratio=1)+facet_wrap(~Treatmentlab)

ggsave("~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Composition/Composition archetypes with Tumor response by treatment.png")

# Composition archetypes
ggplot(umapout, aes(x=V1, y=V2,col=as.factor(archetype) )) + geom_point(size=sz) + theme_classic(base_size=18)+theme(aspect.ratio=1)+
  scale_color_uchicago(name="Composition \n archetype", labels=c("A","B","C","D")) + #jama# jco
  labs(y="TME composition axis 2",x="TME composition axis 1" )
ggsave("~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Composition/Composition archetypes with archetypes annotation.png")
ggsave(paste0(paperfile,"Ribo and Letrozole Composition archetypes with archetypes annotation.png"),height=10,width=10)

ggplot(umapout, aes(x=V1, y=V2,col=as.factor(archetype) )) + geom_point(size=1.5*sz) + theme_classic(base_size=18)+theme(aspect.ratio=1)+
  scale_color_uchicago(name="Composition \n archetype", labels=c("A","B","C","D")) + #jama# jco
  labs(y="TME composition axis 2",x="TME composition axis 1" )+
  theme(axis.title=element_blank(),  axis.text=element_blank(),strip.text = element_blank(),legend.title=element_blank(),legend.text=element_blank() ,legend.position = "none")
ggsave(paste0(paperfile,"BLANK Ribo and Letrozole Composition archetypes with archetypes annotation.png"),height=10,width=10)


# Composition archetypes drivers
ggplot(umapout, aes(x=V1, y=V2,col=logit(eps+ Cancer_cells*(1 - (100)*eps)) , shape=as.factor(archetype)) )+ geom_point(size=sz) + theme_classic(base_size=18)+
  theme(aspect.ratio=1)+labs(y="TME composition axis 2",x="TME composition axis 1" )+
  scale_color_viridis_c(name="Cancer cell \n frequency",option='C',breaks=logit(eps+c(0.005,0.05,0.5)*(1 - 100*eps)), labels=c(0.005,0.05,0.5)) +
  scale_shape_discrete(name="Composition \n archetype", labels= c( "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched", "Cancer \n dominated")) 
ggsave("~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Composition/Composition archetypes with Cancer_cells frequency.png")
ggsave(paste0(paperfile,"Ribo and Letrozole Composition archetypes with Cancer_cells frequency.png"),height=10,width=10)

ggplot(umapout, aes(x=V1, y=V2,col=logit(eps+ Cancer_cells*(1 - (100)*eps)) , shape=as.factor(archetype)) )+ geom_point(size=1.5*sz) + theme_classic(base_size=18)+
  theme(aspect.ratio=1)+labs(y="TME composition axis 2",x="TME composition axis 1" )+
  scale_color_viridis_c(name="Cancer cell \n frequency",option='C',breaks=logit(eps+c(0.005,0.05,0.5)*(1 - 100*eps)), labels=c(0.005,0.05,0.5)) +
  scale_shape_discrete(name="Composition \n archetype", labels= c( "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched", "Cancer \n dominated"),guide=F) +
  theme(axis.title=element_blank(),  axis.text=element_blank(),strip.text = element_blank(),legend.title=element_blank(),legend.text=element_blank() )
ggsave(paste0(paperfile,"BLANK Ribo and Letrozole Composition archetypes with Cancer_cells frequency.png"),height=10,width=10)



ggplot(umapout, aes(x=V1, y=V2, col=logit(eps+(B_cells+Plasma_cells+CD4__T_cells+CD8__T_cells+NK_cells+Tregs+Monocytes+DC+M1_macrophages+M2_macrophages)*(1 - eps)) , shape=as.factor(archetype)) )+ geom_point(size=sz) + theme_classic(base_size=18)+
  theme(aspect.ratio=1)+labs(y="TME composition axis 2",x="TME composition axis 1" ) +
  scale_color_viridis_c(name="Immune cell \n frequency",option='B',breaks=logit(eps+c(0.001,0.01,0.1)*(1 - eps)), labels=c(0.001,0.01,0.1) ) +
  scale_shape_discrete(name="Composition \n archetype", labels= c( "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched", "Cancer \n dominated")) 
ggsave("~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Composition/Composition archetypes with Immune_cells frequency.png")
ggsave(paste0(paperfile,"Ribo and Letrozole Composition archetypes with Immune_cells frequency.png"),height=10,width=10)

ggplot(umapout, aes(x=V1, y=V2, col=logit(eps+(B_cells+Plasma_cells+CD4__T_cells+CD8__T_cells+NK_cells+Tregs+Monocytes+DC+M1_macrophages+M2_macrophages)*(1 - eps)) , shape=as.factor(archetype)) )+
  geom_point(size=1.5*sz) + theme_classic(base_size=18)+
  theme(aspect.ratio=1)+labs(y="TME composition axis 2",x="TME composition axis 1" ) +
  scale_color_viridis_c(name="Immune cell \n frequency",option='B',breaks=logit(eps+c(0.001,0.01,0.1)*(1 - eps)), labels=c(0.001,0.01,0.1) ) +
  scale_shape_discrete(name="Composition \n archetype", labels= c( "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched", "Cancer \n dominated"),guide=F)  +
  theme(axis.title=element_blank(),  axis.text=element_blank(),strip.text = element_blank(),legend.title=element_blank(),legend.text=element_blank() )
ggsave(paste0(paperfile,"BLANK Ribo and Letrozole Composition archetypes with Immune_cells frequency.png"),height=10,width=10)




ggplot(umapout, aes(x=V1, y=V2, col=logit(eps+(CD4__T_cells + CD8__T_cells + NK_cells +Tregs)*(1 - eps)) , shape=as.factor(archetype)) )+ geom_point(size=sz) + theme_classic(base_size=18)+
  theme(aspect.ratio=1)+labs(y="TME composition axis 2",x="TME composition axis 1" ) +
  scale_color_viridis_c(name="T cell \n frequency",option='B',breaks=logit(eps+c(0.001,0.01,0.1)*(1 - eps)), labels=c(0.001,0.01,0.1)) +
  scale_shape_discrete(name="Composition \n archetype", labels= c( "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched", "Cancer \n dominated")) 
ggsave("~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Composition/Composition archetypes with T_cells frequency.png")
ggsave(paste0(paperfile,"Ribo and Letrozole Composition archetypes with T_cells frequency.png"),height=10,width=10)


ggplot(umapout, aes(x=V1, y=V2, col=logit(eps+(CD4__T_cells + CD8__T_cells + NK_cells +Tregs)*(1 - eps)) , shape=as.factor(archetype)) )+ geom_point(size=1.5*sz) + theme_classic(base_size=18)+
  theme(aspect.ratio=1)+labs(y="TME composition axis 2",x="TME composition axis 1" ) +
  scale_color_viridis_c(name="T cell \n frequency",option='B',breaks=logit(eps+c(0.001,0.01,0.1)*(1 - eps)), labels=c(0.001,0.01,0.1)) +
  scale_shape_discrete(name="Composition \n archetype", labels= c( "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched", "Cancer \n dominated"),guide=F) +
  theme(axis.title=element_blank(),  axis.text=element_blank(),strip.text = element_blank(),legend.title=element_blank(),legend.text=element_blank() )
ggsave(paste0(paperfile,"BLANK Ribo and Letrozole Composition archetypes with T_cells frequency.png"),height=10,width=10)




ggplot(umapout, aes(x=V1, y=V2, col=logit(eps+(Adipocytes+Fibroblasts+Lym_Endo+Pericytes+Vas_Endo)*(1 - eps)), shape=as.factor(archetype) ) )+ geom_point(size=sz) + theme_classic(base_size=18)+
  theme(aspect.ratio=1)+labs(y="TME composition axis 2",x="TME composition axis 1" ) +
  scale_color_viridis_c(name="Stromal cell \n frequency",option='D',breaks=logit(eps+c(0.001,0.01,0.1)*(1 - eps)), labels=c(0.001,0.01,0.1)) +
  scale_shape_discrete(name="Composition \n archetype", labels= c( "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched", "Cancer \n dominated")) 
ggsave("~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Composition/Composition archetypes with Stromal_cells frequency.png")
ggsave(paste0(paperfile,"Ribo and Letrozole Composition archetypes with Stromal_cells frequency.png"),height=10,width=10)

ggplot(umapout, aes(x=V1, y=V2, col=logit(eps+(Adipocytes+Fibroblasts+Lym_Endo+Pericytes+Vas_Endo)*(1 - eps)), shape=as.factor(archetype) ) )+ geom_point(size=1.5*sz) + theme_classic(base_size=18)+
  theme(aspect.ratio=1)+labs(y="TME composition axis 2",x="TME composition axis 1" ) +
  scale_color_viridis_c(name="Stromal cell \n frequency",option='D',breaks=logit(eps+c(0.001,0.01,0.1)*(1 - eps)), labels=c(0.001,0.01,0.1)) +
  scale_shape_discrete(name="Composition \n archetype", labels= c( "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched", "Cancer \n dominated"),guide=F) +
  theme(axis.title=element_blank(),  axis.text=element_blank(),strip.text = element_blank(),legend.title=element_blank(),legend.text=element_blank() )
ggsave(paste0(paperfile,"BLANK Ribo and Letrozole Composition archetypes with Stromal_cells frequency.png"),height=10,width=10)





ggplot(umapout, aes(x=V1, y=V2, shape=as.factor(archetype), col=shannonDiv ) )+ geom_point(size=sz)+ theme_classic(base_size=18) +  #,fill=as.factor(archetype)
  theme(aspect.ratio=1)+labs(y="TME composition axis 2",x="TME composition axis 1" ) +
  scale_color_viridis_c(name="Shannon \n Diversity")+
  scale_shape_discrete(name="Composition \n archetype", labels= c( "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched", "Cancer \n dominated")) +
  labs(y="TME composition axis 2",x="TME composition axis 1")+ theme(aspect.ratio=1)#+
ggsave("~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Composition/Composition archetypes with Shannon Diversity.png")
ggsave(paste0(paperfile,"Ribo and Letrozole Composition archetypes with Shannon Diversity.png"),height=10,width=10)


ggplot(umapout, aes(x=V1, y=V2, shape=as.factor(archetype), col=shannonDiv ) )+ geom_point(size=1.5*sz)+ theme_classic(base_size=18) +  #,fill=as.factor(archetype)
  theme(aspect.ratio=1)+labs(y="TME composition axis 2",x="TME composition axis 1" ) +
  scale_color_viridis_c(name="Shannon \n Diversity")+
  scale_shape_discrete(name="Composition \n archetype", labels= c( "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched", "Cancer \n dominated"),guide=F) +
  labs(y="TME composition axis 2",x="TME composition axis 1")+ theme(aspect.ratio=1)+
  theme(axis.title=element_blank(),  axis.text=element_blank(),strip.text = element_blank(),legend.title=element_blank(),legend.text=element_blank() )
ggsave(paste0(paperfile,"BLANK Ribo and Letrozole Composition archetypes with Shannon Diversity.png"),height=10,width=10)

summary(lm(shannonDiv~as.factor(archetype), data=umapout))

# save data to make plot 
save(umapout,file="~/Dropbox/Cancer_pheno_evo/data/FELINE2/Composition/Composition archetype data.RData")



load(file="~/Dropbox/Cancer_pheno_evo/data/FELINE2/Composition/Composition archetype data.RData")
sz<-6

umapout[,"TumorResponse":="Sensitive"]
umapout[dynamic_class3=="Non-response","TumorResponse":="Resistant"]


longdata<- data.table(gather( umapout,celltype,frac,Adipocytes:Vas_Endo) %>%group_by(celltype)%>%summarise(frac=mean(frac)) )[order(-frac)]
longdata[celltype=="DC",celltype:="Dendritic cells"]
longdata[,celltype:= gsub("_"," ",longdata$celltype)]
longdata[,celltype:= gsub("CD4","CD4+",longdata$celltype)]
longdata[,celltype:= gsub("CD8","CD8+",longdata$celltype)]
longdata[,celltype:= gsub("Normal","Diploid",longdata$celltype)]
longdata$celltype <- factor(longdata$celltype, levels=longdata$celltype)

longdata <- data.table( dd0 %>% group_by(Celltype)%>% summarise(count=n() , frac=n() /nrow(dd0) ) ) [order(-frac)]
longdata[,Celltype:= gsub("Normal","Diploid",longdata$Celltype)]
longdata$Celltype <- factor(longdata$Celltype, levels=longdata$Celltype)
ggplot(   longdata , aes(x="",y=frac,fill=Celltype)) + geom_bar(stat="identity",col="white") + coord_polar("y",start=0)+
  theme_void(base_size=22)

paperfile<- "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
ggsave(paste0(paperfile,"Ribo and Letrozole Cell type composition.png"),height=10,width=10)

ggplot(   longdata , aes(x="",y=frac,fill=Celltype)) + geom_bar(stat="identity",col="white") + coord_polar("y",start=0)+
  theme_void(base_size=22)+ theme(legend.position = "none")
ggsave(paste0(paperfile,"BLANK Ribo and Letrozole Cell type composition.png"),height=10,width=10)


ggplot(umapout[Day!=180], aes(y=logit(eps+((Tregs)/(Tregs + CD4__T_cells+CD8__T_cells+NK_cells))), x=TumorResponse, fill=dynamic_class3  ) )+  theme_classic(base_size=18) +  #,fill=as.factor(archetype)
  geom_boxplot(col="black") + geom_point(size=3) + facet_wrap(~Treatmentlab) +
  scale_fill_npg(name="Tumor response", labels=c("Resistant","Sensitive"))  +
  theme(aspect.ratio=1, legend.position = "none")+
  labs(y="T cell frequency", x= "Tumor response")+
  scale_y_continuous(breaks=logit(eps+c(0.001,0.01,0.1) ), labels=c(0.001,0.01,0.1) ) 



ggplot(umapout[Day==180], aes(y=logit(eps+Tregs+CD4__T_cells+CD8__T_cells+NK_cells), x=TumorResponse, fill=dynamic_class3  ) )+  theme_classic(base_size=18) +  #,fill=as.factor(archetype)
 geom_boxplot(col="black") + geom_point(size=3) + facet_wrap(~Treatmentlab) +
  scale_fill_npg(name="Tumor response", labels=c("Resistant","Sensitive"))  +
  theme(aspect.ratio=1, legend.position = "none")+
  labs(y="T cell frequency", x= "Tumor response")+
  scale_y_continuous(breaks=logit(eps+c(0.001,0.01,0.1) ), labels=c(0.001,0.01,0.1) ) 
  
umapout[, DayLab:= paste0("Day ",Day)]
ggplot(umapout[], aes(y=logit(eps+Tregs+CD4__T_cells+CD8__T_cells+NK_cells), x=TumorResponse, fill=dynamic_class3  ) )+  theme_classic(base_size=18) +  #,fill=as.factor(archetype)
  geom_boxplot(col="black") + geom_point(size=3) + facet_wrap(Treatmentlab~DayLab) +
  scale_fill_npg(name="Tumor response", labels=c("Resistant","Sensitive"))  +
  theme(aspect.ratio=1, legend.position = "none")+
  labs(y="T cell frequency", x= "Tumor response")+
  scale_y_continuous(breaks=logit(eps+c(0.001,0.01,0.1) ), labels=c(0.001,0.01,0.1) ) 
ggsave("~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Composition/T cell Composition by response over time.png",width=8)
paperfile<- "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Ribo and Letrozole T cell Composition by response over time.png"),height=10,width=10)


ggplot(umapout[], aes(y=logit(eps+Tregs+CD4__T_cells+CD8__T_cells+NK_cells), x=TumorResponse, fill=dynamic_class3  ) )+  theme_classic(base_size=18) +  #,fill=as.factor(archetype)
  geom_boxplot(col="black") + geom_point(size=3) + facet_wrap(Treatmentlab~DayLab) +
  scale_fill_npg(name="Tumor response", labels=c("Resistant","Sensitive"))  +
  theme(aspect.ratio=1, legend.position = "none")+
  labs(y="T cell frequency", x= "Tumor response")+
  scale_y_continuous(breaks=logit(eps+c(0.001,0.01,0.1) ), labels=c(0.001,0.01,0.1) ) +
  theme(axis.title=element_blank(),  axis.text=element_blank(),strip.text = element_blank(),legend.title=element_blank(),legend.text=element_blank() )
paperfile<- "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"BLANK Ribo and Letrozole T cell Composition by response over time.png"),height=10,width=10)



ggplot(umapout[], aes(y=logit(eps+Tregs+CD4__T_cells+CD8__T_cells+NK_cells), x=log(1+Day),group=interaction(Day,TumorResponse), fill=dynamic_class3  ) )+  theme_classic(base_size=26) +  
  geom_boxplot(col="black") + geom_point(size=3,position = position_dodge(width=2)) + facet_wrap(~Treatmentlab, ncol=1) +
  scale_fill_npg(name="Tumor response", labels=c("Resistant","Sensitive"))  +
  theme(aspect.ratio=1)+#, legend.position = "none")+
  labs(y="T cell frequency", x= "Day")+
  scale_y_continuous(breaks=logit(eps+c(0.001,0.01,0.1) ), labels=c(0.001,0.01,0.1) ) +
  scale_x_continuous(breaks=log(1+c(0,14,180) ), labels=c(0,14,180) )
paperfile<- "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Ribo and Letrozole T cell Composition by response over time version2.png"),height=10,width=10)


ggplot(umapout[], aes(y=logit(eps+Tregs+CD4__T_cells+CD8__T_cells+NK_cells), x=log(1+Day),group=interaction(Day,TumorResponse), fill=dynamic_class3  ) )+  theme_classic(base_size=18) +  
  geom_boxplot(col="black") + geom_point(size=3,position = position_dodge(width=2)) + facet_wrap(~Treatmentlab, ncol=1) +
  scale_fill_npg(name="Tumor response", labels=c("Resistant","Sensitive"))  +
  theme(aspect.ratio=1, legend.position = "none")+
  labs(y="T cell frequency", x= "Tumor response")+
  scale_y_continuous(breaks=logit(eps+c(0.001,0.01,0.1) ), labels=c(0.001,0.01,0.1) ) +
  scale_x_continuous(breaks=log(1+c(0,14,180) ), labels=c(0,14,180) )+
  theme(axis.title=element_blank(),  axis.text=element_blank(),strip.text = element_blank(),legend.title=element_blank(),legend.text=element_blank() )

paperfile<- "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"BLANK Ribo and Letrozole T cell Composition by response over time version2.png"),height=10,width=10)




umapout[,success:= round((Tregs+CD4__T_cells+CD8__T_cells+NK_cells)*totalCount)]
umapout[,fail:= totalCount-success]

summary( lm(logit(eps+Tregs+CD4__T_cells+CD8__T_cells+NK_cells)~Treatmentlab*dynamic_class3 , umapout[Day==180]) )
summary( glm(cbind(success,fail)   ~Treatmentlab*dynamic_class3 , family="binomial", umapout[Day==180]) )

summary( glm(cbind(success,fail)   ~dynamic_class3 , family="binomial", umapout[Treatmentlab=="Letrozole alone"]) )


# Supplementary Treg assessment
require(pscl)
zmod<- zeroinfl(round(1e6*Tregs)~ I(eps+M2_macrophages/(eps+M1_macrophages+M2_macrophages)),dist="negbin",data= umapout )
umapout$Tregspred <- predict(zmod, type="response")/1e6
ggplot(umapout, aes(y=logit(eps+Tregs), x= logit(eps+Tregspred)  ) )  + geom_point(size=sz) + theme_classic(base_size=18)+  
  geom_smooth(method="gam",se=F)#,fill=as.factor(archetype)

options(scripen=10000)

ggplot(umapout, aes(y= logit(eps + Tregs), x= M2_macrophages/(M1_macrophages + M2_macrophages)) )  + geom_point(size= sz) + 
  theme_classic(base_size=22)+  
  geom_smooth( data= umapout, aes(y= logit(eps + Tregspred) ) , se= F  ) +
  labs(y= "Regulatory T cell abundance \n (Treg fraction)", x= "Proportion of myeloid cells in M2 state" )+
  theme(aspect.ratio= 1) + scale_y_continuous(breaks=logit(eps + c(0.0005, 0.005, 0.05, 0.5 )) ,labels=  format( c(0.0005, 0.005, 0.05, 0.5 ) , scientific=F)   )
ggsave("~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Composition/Composition Treg vs M1M2 proportion.png")

zmodCD8 <- zeroinfl(round(1e8*(CD8__T_cells + NK_cells)) ~ I(eps + DC + M1_macrophages + M2_macrophages), dist="negbin", data= umapout )
zmodMy <- zeroinfl(round(1e8*(DC + M1_macrophages + M2_macrophages)) ~ I(eps + CD8__T_cells + NK_cells), dist="negbin", data= umapout )
umapout$CD8pred <- predict(zmodCD8, type="response")/1e8
umapout$Mypred <- predict(zmodMy, type="response")/1e8


ggplot(umapout, aes(y=logit(eps + CD8__T_cells + NK_cells), x= logit(eps+CD8pred)  ) )  + geom_point(size=sz) + theme_classic(base_size=18)+  
  geom_smooth(method="gam",se=F)#,fill=as.factor(archetype)

ggplot(umapout, aes(y= logit(eps +CD8__T_cells + NK_cells), x= logit( DC + M1_macrophages + M2_macrophages)) )  + 
  geom_point(size= sz) + theme_classic(base_size=22) +  
  geom_smooth( data= umapout, aes(y= logit(eps + CD8pred) ) , se= F  ) +
  labs(y= "CD8 T cell abundance", x= "Myeloid cell abundance" ) +
  theme(aspect.ratio= 1) + 
  scale_y_continuous(breaks= logit(eps + c(0.0005, 0.005, 0.05, 0.5 )) , labels=  format( c(0.0005, 0.005, 0.05, 0.5 ) , scientific=F)   ) +
  scale_x_continuous(breaks= logit(eps + c(0.0005, 0.005, 0.05, 0.5 )) , labels=  format( c(0.0005, 0.005, 0.05, 0.5 ) , scientific=F)   )
ggsave("~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Composition/Composition CD8 Tcell vs myeloid abundance.png")


transitionsInput<-na.omit(umapout[Day!=14]%>% 
                            mutate(DayClass = paste0("Day",Day),
                                   archetypeClass = paste0("Archetype",archetype)) %>%
                            select(Treatmentlab,dynamic_class3,Patient.Study.ID,DayClass,archetypeClass)%>% spread(DayClass,archetypeClass))

transitionsInput[dynamic_class3=="Response"]$Day180%>%table()
transitionsInput[dynamic_class3!="Response"]$Day180%>%table()
transitionsprobs <- data.table( data.table( table( transitionsInput%>%select(Day0,Day180,Treatmentlab,dynamic_class3 )) )  %>%group_by(Treatmentlab,dynamic_class3,Day0) %>%mutate(prop=N/sum(N)) )

ggplot(umapout, aes(x=V1, y=V2,col=as.factor(archetype) )) + geom_point(size=sz) + theme_classic(base_size=18)+theme(aspect.ratio=1)+
  scale_color_uchicago()

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
ggsave("~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Composition/State transition probability.png",width=10)
paperfile<- "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Ribo and Letrozole T cell Composition State transition probability.png"),height=10,width=10)

ggplot(transitionsprobs,aes(x=From,y=To, alpha=(prop) ))+geom_tile() + theme_classic(base_size=12)+ facet_wrap(Treatmentlab~dynamic_class3) + 
  scale_alpha_continuous(name="Transition probability",na.value=0)+theme(aspect.ratio=1,axis.text.x=element_text(angle=90,vjust = 0.5))+
  theme(axis.title=element_blank(),  axis.text=element_blank(),strip.text = element_blank(),legend.title=element_blank(),legend.text=element_blank() )
paperfile<- "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"BLANK Ribo and Letrozole T cell Composition State transition probability.png"),height=10,width=10)



ggplot(transitionsprobs[Day180=="Archetype1"],aes(y=From,x=dynamic_class3, alpha=(prop) ))+geom_tile() + theme_classic(base_size=12)+ facet_wrap(~Treatmentlab) + 
  scale_alpha_continuous(name="Transition probability",na.value=0)+theme(aspect.ratio=1)
ggsave("~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Composition/State transition probability.png",width=10)


