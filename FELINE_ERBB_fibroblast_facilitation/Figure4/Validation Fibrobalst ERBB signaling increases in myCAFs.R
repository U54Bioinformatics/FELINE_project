rm(list=ls())
require(mgcv);require(data.table);require(dplyr);require(ggplot2);require(tidyr);require(lme4);require(lmerTest);require(parallel)
library(effects);require(umap)
require(Rfast);require(ider)
library("dendextend");library(ggdendro);require(ggsci);require(viridis)
require("Rdimtools")

# Load clinical data
load(file= "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE/Feline clinical input.RData" )

# Load metadata for a specific cell type
gsea_path <- "~/Dropbox/FELINE Project (1)/Data_analysis/scRNA/05_ssGSEA_score/Signature_c2_hallmark/results/"   #"~/Dropbox/FELINE/Data_share/Modeling_Data/pathway_seperated_files/Data_gene_count_per_celltype_model_zinbwave_ssGSEA/FEL001043/"

cell_types_all <-c("Fibroblasts")
DAY=c(0,14,180)
ARMS <-c("A","B","C")
metadd <- rbindlist( lapply(cell_types_all , function (cell_type_i){
  annotation.file <- paste0( "~/Dropbox/FELINE Project (1)/Data_analysis/scRNA/01_metadata/results/FEL001046_meta_for_each_celltype/FEL001046_", cell_type_i, "_scRNA.metadata.txt")
  cell_type_meta_dd <- data.table(  fread(annotation.file))[Celltype_subtype!="Low-quality cells"]
  setnames( cell_type_meta_dd,old="Sample",new="Sample_p_t")
  cell_type_meta_dd[, c("Sample", "Timepoint") := tstrsplit(Sample_p_t, "_", fixed=TRUE)]
  cell_type_meta_dd[,Day:=0]  ; cell_type_meta_dd[grepl("_M",Sample_p_t),Day:=14]   ; cell_type_meta_dd[grepl("_E",Sample_p_t),Day:=180]   
  cell_type_meta_dd[,day_fact:=as.factor(Day)]
  cell_type_meta_dd[,file_string:=cell_type_i]
  # merge response scores
  FULL_cell_type_meta_dd <- merge(cell_type_meta_dd,response_code_dd,by="Sample")      
}) )[Platform!="ICELL8"]

load("~/Dropbox/Cancer_pheno_evo/data/FELINE2/PhenotypesAllArmsCohort2/UpdatedRevisednewFibroblastsFibroblasts.RData")
u_dat[, Timepoint:="Pre treatment"]
u_dat[Day==180, Timepoint:="Post treatment"]

# V2 correlates strongly with EMT in the cohort 2 fibroblast umap projection
ggplot( u_dat,aes(V1,V2,col=  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION)) +  geom_point(size=1)+theme_classic()+theme(aspect.ratio=1)+facet_wrap(~Day)

# extract  just the core genes to test
get_gene_umap <- function(){
  geneList<-c("NRG1",  "NRG2",  "NRG3",  "EGF", "HBEGF", "TGFA","ANXA1")
  ### Load counts per million read data for a specific patient's samples and gen Ligand+Receptor genes
  CPMlocs1 <- "/Users/jason/Dropbox/FELINE Project (1)/FELINE Cohort 2/DNANexusCopy/CPM/HQ/"
  CPMlocs2 <- "/Users/jason/Dropbox/FELINE Project (1)/FELINE Cohort 2/DNANexusCopy/CPM/HQLQ/"

  CPMfiles1 <- list.files(CPMlocs1)[ grep("CPM", list.files(CPMlocs1) ) ]
  CPMfiles2 <- list.files(CPMlocs2)[ grep("CPM", list.files(CPMlocs2) ) ]
  n10Xpats <- length(CPMfiles1)
  
  cpm_i1 <- data.table( fread( paste0(CPMlocs1, CPMfiles1[grep("Fibroblast",CPMfiles1)]) ) )[V1%in%geneList]   # load full gene expression
  names(cpm_i1)[1]<- "Gene.ID"
  transposedHQ <- as.data.table( t(cpm_i1[,-1])  , keep.rownames = T)
  colnames(transposedHQ) <- c("Cell.ID", cpm_i1$Gene.ID)
  rm(list="cpm_i1")
  
  cpm_i2<- data.table( fread( paste0(CPMlocs2, CPMfiles2[grep("Fibroblast",CPMfiles2)]) ) )[V1%in%geneList]   # load full gene expression
  names(cpm_i2)[1]<- "Gene.ID"
  transposedLQ <- as.data.table( t(cpm_i2[,-1])  , keep.rownames = T)
  colnames(transposedLQ) <- c("Cell.ID", cpm_i2$Gene.ID)
  rm(list="cpm_i2")
  
  srtGene.ID0<-sort(intersect(names(transposedHQ)[-1],names(transposedLQ)[-1]))
  srtGene.ID <- sort(intersect(srtGene.ID0,geneList))
  
  transposedHQ<- transposedHQ[,c("Cell.ID",srtGene.ID),with=F]
  transposedLQ<- transposedLQ[,c("Cell.ID",srtGene.ID),with=F]
  transposedHQ[1:10,]
  transposedLQ[1:10,]
  transposedall <- rbind(transposedHQ, transposedLQ)
  rm(list="transposedHQ")
  rm(list="transposedLQ")
  return(transposedall)
}

gene_u_dat <- get_gene_umap()

u_dat$Timepoint<- factor(u_dat$Timepoint, levels=c("Pre treatment",  "Post treatment"))

ggplot( merge(u_dat[Day!=14],gene_u_dat , by="Cell.ID" ),
        aes(V1,V2,col=  log(1+ NRG1 + NRG2 + NRG3 + EGF + HBEGF + TGFA+ ANXA1     # 
        )))+theme_classic(base_size=22) +  geom_point(size=1.5)+
  facet_wrap(~Timepoint,nrow=2)+theme(aspect.ratio=1) + scale_color_viridis(option="D",name="Fibroblast GF signaling \n to ERBB family receptors", labels=c(1,10,100,1000), breaks=log(c(1,10,100,1000) )) +
  labs(y="",x="")+ theme(legend.position="top")+guides(colour=guide_colourbar(barwidth=20))
#ggsave("~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Validation Fibroblast signaling UMAP start vs end.png",height=10,width=10)

plotthis<- merge(u_dat[Day!=14][]%>%dplyr::select( Cell.ID,Quality,Timepoint,V1,V2) ,gene_u_dat , by="Cell.ID" )
plotthis <- data.table(plotthis%>%gather(var,val,ANXA1:TGFA))

plotthis[,emtlevel:= round(V2*1.87)/1.87]
length(unique(plotthis$emtlevel))
plotthis[,egfproliflevel:= round(V1*1.925)/1.925]
length(unique(plotthis$egfproliflevel))
plotthis[,scaleval:=scale( log(1+val) ) ,by=var]

ggplot(plotthis[][var=="NRG1"],aes(y= scaleval ,x=emtlevel))+geom_boxplot(aes(group=emtlevel)) +geom_smooth(method="gam", formula= y~s(x,k=3))+facet_wrap(~Timepoint,scales="free_y")
ggplot(plotthis[val>0][var=="HBEGF"],aes(y= scaleval ,x=emtlevel))+geom_point() +geom_smooth(method="gam", formula= y~s(x,k=4))+facet_wrap(~var,scales="free_y")

plotthisCell <- data.table( plotthis[]%>%group_by(Cell.ID,Timepoint)%>%dplyr::summarise(emtlevel=emtlevel, egfproliflevel= egfproliflevel, lnvalmu= mean( scaleval )  )  )
plotthisCell[,egfproliflevelB:="Quiescent"]
plotthisCell[egfproliflevel==T,egfproliflevelB:="EGF proliferative"]
plotthisCell$egfproliflevelB <- factor(plotthisCell$egfproliflevelB , levels= c("Quiescent","EGF proliferative"  ) )

validStatsData <- plotthisCell 
#write.csv(validStatsData,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/ValidationCohort_Stats_FibroblastERBBligandexpressionincreasesmyCAFs.csv")
validStatsData <- data.table(read.csv(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/Figure4/ValidationCohort_Stats_FibroblastERBBligandexpressionincreasesmyCAFs.csv"))


ggplot(validStatsData[],aes(y= lnvalmu ,x=emtlevel,group=emtlevel,fill=emtlevel))+
  theme_classic(base_size=26)+theme(aspect.ratio=1)+
  labs(x="Fibroblast myCAF phenotype \n (TGFb stimulated)",y="ERBB ligand signaling \n (scaled)")+
  geom_boxplot() +stat_boxplot(geom = "errorbar",
                               width = 0.5) +geom_point(size=0.4) +
  scale_fill_gradient(low =pal_aaas()(6)[3] ,high =pal_aaas()(6)[1])+theme(legend.position = "none")
#ggsave("~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Validation Fibroblast ERBBligandexpression increases in myCAFs.pdf",height=8,width=8,dpi = 320)

lm(lnvalmu~as.factor(emtlevel),data= validStatsData[])%>%summary()

lm(lnvalmu~(emtlevel),data= validStatsData[])%>%summary()

m1<- aov(lnvalmu~as.factor(round(10*emtlevel)),data= discStatsData)
plot(TukeyHSD( m1, conf.level=.95), las = 2)


