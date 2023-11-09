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

load("~/Dropbox/Cancer_pheno_evo/data/FELINE2/PhenotypesAllArms/FibroblastsFibroblasts.RData")

corgenes <- data.table(t(
  corVall
),keep.rownames = T)
corgenes[order(-abs(V1))][1:20]
corgenes[order(-abs(V2))][1:20]

corgenes[order(-abs(V3))][1:20]
corgenes[order(-abs(V4))][1:10]

loc1<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Fibroblasts phenotypes AllArms/"
setwd(loc1)
#ggplot( u_dat,aes(V1,V2,col=  Patient.Study.ID))+theme_classic() +  geom_point(alpha=0.4,size=0.4)+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
#ggsave(filename = "Fibroblast ssgsea phenotype landscape per patient by response and day.png")

u_dat[, Timepoint:="Pre treatment"]
u_dat[Day==180, Timepoint:="Post treatment"]

cordf<-data.table( t(corVall) , keep.rownames = T)
cordf[order(-abs(V1))][1:20]
cordf[order(-abs(V2))][1:20]
cordf[order(-abs(V3))][1:20]

ggplot( u_dat[ARM!="A"],aes(V1,V2,col=  dynamic_class3))+theme_classic() +  geom_point(size=1)+facet_wrap(dynamic_class3~Day)+theme(aspect.ratio=1)
ggsave(filename = "Fibroblast ssgsea phenotype landscape by response and day.png")

ggplot( u_dat,aes(V1,-V3,col=  AMIT_SERUM_RESPONSE_40_MCF10A)) +  geom_point(size=1)+facet_wrap(dynamic_class3~Day)+theme_classic()+theme(aspect.ratio=1)
ggsave(filename = "Fibroblast ssgsea phenotype landscape overlay AMIT_SERUM_RESPONSE_40_MCF10A.png")
ggplot( u_dat,aes(V1,-V3,col=  NAGASHIMA_EGF_SIGNALING_UP)) +  geom_point(size=1)+facet_wrap(dynamic_class3~Day)+theme_classic()+theme(aspect.ratio=1)
ggsave(filename = "Fibroblast ssgsea phenotype landscape overlay NAGASHIMA_EGF_SIGNALING_UP.png")
ggplot( u_dat,aes(V1,-V3,col=  BIOCARTA_EGF_PATHWAY)) +  geom_point(size=1)+facet_wrap(dynamic_class3~Day)+theme_classic()+theme(aspect.ratio=1)
ggsave(filename = "Fibroblast ssgsea phenotype landscape overlay BIOCARTA_EGF_PATHWAY.png")
ggplot( u_dat[ARM!="A"],aes(V1,V2,col=  NAGASHIMA_NRG1_SIGNALING_UP)) +  geom_point(size=1)+facet_wrap(dynamic_class3~Day)+theme_classic()+theme(aspect.ratio=1)
ggsave(filename = "Fibroblast ssgsea phenotype landscape overlay NAGASHIMA_NRG1_SIGNALING_UP.png")


ggplot( u_dat[ARM!="A"],aes(V1,V2,col=  PHONG_TNF_TARGETS_UP)) +  geom_point(size=1)+facet_wrap(dynamic_class3~Day)+theme_classic()+theme(aspect.ratio=1)
ggsave(filename = "Fibroblast ssgsea phenotype landscape overlay PHONG_TNF_TARGETS_UP.png")
ggplot( u_dat[ARM!="A"],aes(V1,V2,col=  HALLMARK_TNFA_SIGNALING_VIA_NFKB)) +  geom_point(size=1)+facet_wrap(dynamic_class3~Day)+theme_classic()+theme(aspect.ratio=1)
ggsave(filename = "Fibroblast ssgsea phenotype landscape overlay HALLMARK_TNFA_SIGNALING_VIA_NFKB.png")
ggplot( u_dat[ARM!="A"],aes(V1,V2,col=  PID_IL2_1PATHWAY)) +  geom_point(size=1)+facet_wrap(dynamic_class3~Day)+theme_classic()+theme(aspect.ratio=1)
ggsave(filename = "Fibroblast ssgsea phenotype landscape overlay PID_IL2_1PATHWAY.png")
ggplot( u_dat[ARM!="A"],aes(V1,V2,col=  BIOCARTA_IL6_PATHWAY)) +  geom_point(size=1)+facet_wrap(dynamic_class3~Day)+theme_classic()+theme(aspect.ratio=1)
ggsave(filename = "Fibroblast ssgsea phenotype landscape overlay BIOCARTA_IL6_PATHWAY.png")
ggplot( u_dat[ARM!="A"],aes(V1,V2,col=  PID_IL6_7_PATHWAY)) +  geom_point(size=1)+facet_wrap(dynamic_class3~Day)+theme_classic()+theme(aspect.ratio=1)
ggsave(filename = "Fibroblast ssgsea phenotype landscape overlay PID_IL6_7_PATHWAY.png")
ggplot( u_dat[ARM!="A"],aes(V1,V2,col=  REACTOME_IL_2_SIGNALING)) +  geom_point(size=1)+facet_wrap(dynamic_class3~Day)+theme_classic()+theme(aspect.ratio=1)
ggsave(filename = "Fibroblast ssgsea phenotype landscape overlay REACTOME_IL_2_SIGNALING.png")
ggplot( u_dat[ARM!="A"],aes(V1,V2,col=  PLASARI_TGFB1_TARGETS_10HR_UP)) +  geom_point(size=1)+facet_wrap(dynamic_class3~Day)+theme_classic()+theme(aspect.ratio=1)
ggsave(filename = "Fibroblast ssgsea phenotype landscape overlay PLASARI_TGFB1_TARGETS_10HR_UP.png")
ggplot( u_dat[ARM!="A"],aes(V1,V2,col=  KEGG_CHEMOKINE_SIGNALING_PATHWAY)) +  geom_point(size=1)+facet_wrap(dynamic_class3~Day)+theme_classic()+theme(aspect.ratio=1)
ggsave(filename = "Fibroblast ssgsea phenotype landscape overlay KEGG_CHEMOKINE_SIGNALING_PATHWAY.png")


ggplot( u_dat[ARM!="A"],aes(dynamic_class3,V3,fill=  dynamic_class3)) +  geom_violin()+facet_wrap(~Day)+theme_classic()+theme(aspect.ratio=1)


get_gene_umap <- function(){
  load( file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Gene expression landscapes/fibroblast gene expression landscape/fibroblast gene expression landscape.RData")
  return(u_dat)
}
gene_u_dat <- get_gene_umap()

v3cortest<-cor(gene_u_dat$V3, log(1+ gene_u_dat[,19:21297] ))
v3cortest2<-data.table( t(v3cortest) ,keep.rownames=T)
v3cortest2[order(-abs(V1))][1:30]
v1cortest<-cor(gene_u_dat[abs(V1)<5][abs(V3)<5]$V1, log(1+ gene_u_dat[abs(V1)<5][abs(V3)<5][,19:21297] ))
v1cortest2<-data.table( t(v1cortest) ,keep.rownames=T)
v1cortest2[order(-abs(V1))][1:30]

u_dat$Timepoint<- factor(u_dat$Timepoint, levels=c("Pre treatment",  "Post treatment"))
ggplot( merge(u_dat[Day!=14],gene_u_dat%>%dplyr::select(c("Cell.ID","FAP")) , by="Cell.ID" ),
        aes(V1,-V3,col=  REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION))+theme_classic() +  geom_point(size=1)+
  facet_wrap(~Timepoint)+theme(aspect.ratio=1) + scale_color_viridis(option="B")+
  theme(axis.text = element_blank(),axis.title = element_blank(),legend.text = element_blank(),legend.title = element_blank(),strip.background = element_blank(),strip.text = element_blank())

ggplot( merge(u_dat[Day!=14],gene_u_dat%>%dplyr::select(c("Cell.ID","FAP","ACTA2")) , by="Cell.ID" ),
        aes(V1,-V3,col=  log(1+FAP+ACTA2)))+theme_classic() +  geom_point(size=1)+
  facet_wrap(~Timepoint)+theme(aspect.ratio=1) + scale_color_viridis(option="B")+
  theme(axis.text = element_blank(),axis.title = element_blank(),legend.text = element_blank(),legend.title = element_blank(),strip.background = element_blank(),strip.text = element_blank())


ggplot( merge(u_dat[Day!=14],gene_u_dat%>%dplyr::select(c("Cell.ID","NRG1")) , by="Cell.ID" ),
        aes(V1,-V3,col=  log(1+NRG1)))+theme_classic() +  geom_point(size=1)+
  facet_wrap(~Timepoint)+theme(aspect.ratio=1) + scale_color_viridis(option="B")+
  theme(axis.text = element_blank(),axis.title = element_blank(),legend.text = element_blank(),legend.title = element_blank(),strip.background = element_blank(),strip.text = element_blank())






ggplot( merge(u_dat[Day!=14],gene_u_dat%>%dplyr::select(c("Cell.ID","EGFR","TGFBR1","TGFBR2","TGFBR3")) , by="Cell.ID" ),
        aes(V1,-V3,col=  log(1+EGFR+TGFBR1+TGFBR2+TGFBR3     
        )))+theme_classic(base_size=22) +  geom_point(size=1.5)+
  facet_wrap(~Timepoint)+theme(aspect.ratio=1) + scale_color_viridis(option="D",name="Fibroblast growth factor signaling \n to ERBB family receptors", labels=c(1,10,100,1000), breaks=log(c(1,10,100,1000) )) +
  labs(y="",x="")+ theme(legend.position="top")+guides(colour=guide_colourbar(barwidth=20))


ggplot( merge(u_dat[Day!=14],gene_u_dat%>%dplyr::select(c("Cell.ID","TGFBR1","TGFBR2","TGFBR3")) , by="Cell.ID" ),
        aes(V1,-V3,col=   log(1+TGFBR1+TGFBR2+TGFBR3  )  ))+theme_classic(base_size=16) +  geom_point(size=1.5)+
  theme(aspect.ratio=1) + #scale_color_viridis(option="B", name="Hallmark EMT \n ssGSEA score") +
  labs(y="",x="")+ theme(legend.position="top")+
  scale_color_viridis(option="A",name="Fibroblast TGFBR1-3 \n receptor expression", labels=c(1,10,100,1000), breaks=log(c(1,10,100,1000) ))+
  guides(colour=guide_colourbar(barwidth=10))#
#ggsave("~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Fibroblast TGFBR UMAP expression.png",height=10,width=5)


ggplot( merge(u_dat[Day!=14],gene_u_dat%>%dplyr::select(c("Cell.ID","EGFR","TGFBR1","TGFBR2","TGFBR3")) , by="Cell.ID" ),
        aes(V1,-V3,col=  log(1+TGFBR1+TGFBR2+TGFBR3     
        )))+theme_classic(base_size=22) +  geom_point(size=1.5)+
  theme(aspect.ratio=1) + scale_color_viridis(option="D",name="Fibroblast \n TGFBR1-3 GF\n receptor expression", labels=c(1,10,100,1000), breaks=log(c(1,10,100,1000) )) +
  labs(y="",x="")#+ #theme(legend.position="top")+
#ggsave("~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Fibroblast TGFR UMAP start vs end.png",height=10,width=5)




ggplot( merge(u_dat[Day!=14],gene_u_dat%>%dplyr::select(c("Cell.ID","NRG1","NRG2","NRG3","EGF","HBEGF","AREG","TGFA","SEMA4D","HLA_A","EFNB1","ICAM1","GNAI2","CDH1","ANXA1","ADAM17")) , by="Cell.ID" ),
        aes(V1,-V3,col=  log(1+ NRG1 + NRG2 + NRG3 + EGF + HBEGF + TGFA+ ANXA1     
        )))+theme_classic(base_size=22) +  geom_point(size=1.5)+
  facet_wrap(~Timepoint,nrow=2)+theme(aspect.ratio=1) + scale_color_viridis(option="D",name="Fibroblast GF signaling \n to ERBB family receptors", labels=c(1,10,100,1000), breaks=log(c(1,10,100,1000) )) +
  labs(y="",x="")+ theme(legend.position="top")+guides(colour=guide_colourbar(barwidth=20))
#ggsave("~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Fibroblast signaling UMAP start vs end.png",height=10,width=10)

szval<-30
ggplot( merge(u_dat[Day!=14],gene_u_dat%>%dplyr::select(c("Cell.ID","NRG1","NRG2","NRG3","EGF","HBEGF","AREG","TGFA","SEMA4D","HLA_A","EFNB1","ICAM1","GNAI2","CDH1","ANXA1","ADAM17")) , by="Cell.ID" ),
        aes(V1,-V3,col=  log(1+NRG1+NRG2+NRG3+EGF+HBEGF +TGFA+ANXA1     
        )))+theme_classic(base_size=szval) +  geom_point(size=1.5)+
  facet_wrap(~Timepoint)+theme(aspect.ratio=1) + scale_color_viridis(option="D",name="Fibroblast signaling of GF \n to ERBB family receptors", labels=c(1,10,100,1000), breaks=log(c(1,10,100,1000) )) +
  labs(y="",x="")+ theme(legend.position="top")+guides(colour=guide_colourbar(barwidth=20))+
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        #legend.text = element_blank(),
        legend.title = element_blank(),
        #strip.background = element_blank(),strip.text = element_blank()
        strip.text=element_text(size=szval*1.2)
  )
ggsave("~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/BLANK Fibroblast signaling UMAP start vs end.png",height=10,width=10)

ggplot( merge(u_dat[Day!=14],gene_u_dat%>%dplyr::select(c("Cell.ID","NRG1","NRG2","NRG3","EGF","HBEGF","AREG","TGFA","SEMA4D","HLA_A","EFNB1","ICAM1","GNAI2","CDH1","ANXA1","ADAM17")) , by="Cell.ID" ),
        aes(V1,-V3,col=  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION  ))+theme_classic(base_size=szval) +  geom_point(size=1.5)+
  theme(aspect.ratio=1) + scale_color_viridis(option="B", name="Hallmark EMT \n pathway") +
  labs(y="",x="")+ theme(legend.position="top")+guides(colour=guide_colourbar(barwidth=20))+
  theme(axis.text = element_blank(),#,axis.title = element_blank(),legend.text = element_blank(),
        legend.title = element_blank(),strip.background = element_blank(),strip.text = element_blank())
ggsave("~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/BLANK Fibroblast HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION  UMAP.png",height=10,width=6)

ggplot( merge(u_dat[Day!=14],gene_u_dat%>%dplyr::select(c("Cell.ID","NRG1","NRG2","NRG3","EGF","HBEGF","AREG","TGFA","SEMA4D","HLA_A","EFNB1","ICAM1","GNAI2","CDH1","ANXA1","ADAM17")) , by="Cell.ID" )[(NRG1+NRG2+NRG3+EGF+HBEGF +TGFA+ANXA1)>0],
        aes(x= V1,y=  log(NRG1+NRG2+NRG3+EGF+HBEGF +TGFA+ANXA1 )  ))+theme_classic(base_size=szval) +  geom_point(size=1.5)+
  theme(aspect.ratio=1) + scale_color_viridis(option="B", name="Hallmark EMT \n pathway") +
  labs(y="",x="")+ theme(legend.position="top")+guides(colour=guide_colourbar(barwidth=20))+
  geom_smooth(method="lm")+facet_wrap((-V3>0)~Day)




ggplot( merge(u_dat[Day!=14],gene_u_dat%>%dplyr::select(c("Cell.ID","NRG1","NRG2","NRG3","EGF","HBEGF","AREG","TGFA","SEMA4D","HLA_A","EFNB1","ICAM1","GNAI2","CDH1","ANXA1","ADAM17")) , by="Cell.ID" ),
        aes(V1,-V3,col=  BIOCARTA_EGF_PATHWAY ))+theme_classic(base_size=szval) +  geom_point(size=1.5)+
  theme(aspect.ratio=1) + scale_color_viridis(option="C", name="Biocarta EGF \n pathway") +
  labs(y="",x="")+ theme(legend.position="top")+guides(colour=guide_colourbar(barwidth=20))+
  theme(axis.text = element_blank(),#axis.title = element_blank(),legend.text = element_blank(),
        legend.title = element_blank(),strip.background = element_blank(),strip.text = element_blank())
ggsave("~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/BLANK Fibroblast BIOCARTA_EGF_PATHWAY UMAP.png",height=10,width=6)





ggplot( merge(u_dat[Day!=14],gene_u_dat%>%dplyr::select(c("Cell.ID","NRG1","NRG2","NRG3","EGF","HBEGF","AREG","TGFA","SEMA4D","HLA_A","EFNB1","ICAM1","GNAI2","CDH1","ANXA1","ADAM17")) , by="Cell.ID" ),
        aes(V1,-V3,col=  BIOCARTA_ERK_PATHWAY ))+theme_classic(base_size=16) +  geom_point(size=1.5)+
  theme(aspect.ratio=1) + scale_color_viridis(option="C", name="Biocarta ERK pathway \n ssGSEA score") +
  labs(y="",x="")+ theme(legend.position="top")#+
#ggsave("~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/BLANK Fibroblast BIOCARTA_ERK_PATHWAY UMAP.png",height=10,width=5)

u_dat[ ,Timepoint:="Post treatment"]
u_dat[Day==0 ,Timepoint:="Pre treatment"]

plotthis<- merge(u_dat[Day!=14]%>%dplyr::select( Cell.ID,Timepoint,V1,V3) ,gene_u_dat%>%dplyr::select(c("Cell.ID","NRG1","NRG2","NRG3","EGF","HBEGF","AREG","TGFA","SEMA4D","HLA_A","EFNB1","ICAM1","GNAI2","CDH1","ANXA1","ADAM17")) , by="Cell.ID" )
plotthis <- data.table(plotthis%>%gather(var,val,NRG1:ADAM17))
plotthis[,emtlevel:= round(-V3*1.15)/1.15]
length(unique(plotthis$emtlevel))
plotthis[,egfproliflevel:= round(V1*1.925)/1.925]
length(unique(plotthis$egfproliflevel))
plotthis[,scaleval:=scale( log(1+val) ) ,by=var]

ggplot(plotthis[][var=="NRG1"],aes(y= scaleval ,x=emtlevel))+geom_boxplot(aes(group=emtlevel)) +geom_smooth(method="gam", formula= y~s(x,k=3))+facet_wrap(~Timepoint,scales="free_y")
ggplot(plotthis[val>0][var=="HBEGF"],aes(y= scaleval ,x=emtlevel))+geom_point() +geom_smooth(method="gam", formula= y~s(x,k=4))+facet_wrap(~var,scales="free_y")

plotthisCell <- data.table( plotthis%>%group_by(Cell.ID,Timepoint)%>%dplyr::summarise(emtlevel=emtlevel, egfproliflevel= egfproliflevel, lnvalmu= mean( scaleval )  )  )
plotthisCell[,egfproliflevelB:="Quiescent"]
plotthisCell[egfproliflevel==T,egfproliflevelB:="EGF proliferative"]
plotthisCell$egfproliflevelB <- factor(plotthisCell$egfproliflevelB , levels= c("Quiescent","EGF proliferative"  ) )


discStatsData <- plotthisCell 
#write.csv(discStatsData,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/Figure4/DiscoveryCohort_Stats_FibroblastERBBligandexpressionincreasesmyCAFs.csv")
discStatsData <- data.table(read.csv(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/Figure4/DiscoveryCohort_Stats_FibroblastERBBligandexpressionincreasesmyCAFs.csv"))

ggplot(discStatsData[],aes(y= lnvalmu ,x=emtlevel,group=emtlevel,fill=emtlevel))+
  theme_classic(base_size=26)+theme(aspect.ratio=1)+
  labs(x="Fibroblast myCAF phenotype \n (TGFb stimulated)",y="ERBB ligand signaling \n (scaled)")+
  geom_boxplot() +stat_boxplot(geom = "errorbar",
                               width = 0.5) +geom_point(size=0.4) +
  scale_fill_gradient(low =pal_aaas()(6)[3] ,high =pal_aaas()(6)[1])+theme(legend.position = "none")
#ggsave("~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Fibroblast ERBBligandexpression increases in myCAFs.pdf",height=8,width=8,dpi = 320)

lm(lnvalmu~as.factor(emtlevel),data= discStatsData[])%>%summary()

lm(lnvalmu~(emtlevel),data= plotthisCell[])%>%summary()

m1<- aov(lnvalmu~as.factor(round(10*emtlevel)),data= discStatsData)
plot(TukeyHSD( m1, conf.level=.95), las = 2)


