rm(list=ls())
library(glmnet); require(dplyr); require(data.table); require(tidyr)
require(ggplot2); require(pheatmap); require(ggsci); require(lme4);require(lmerTest)
require(scales)
require(viridis)

# Select one pathway to analyze
wch_signalTrans <-c("HALLMARK_ESTROGEN_RESPONSE_EARLY","HALLMARK_ESTROGEN_RESPONSE_LATE","BIOCARTA_EGF_PATHWAY","KEGG_TGF_BETA_SIGNALING_PATHWAY","HALLMARK_TGF_BETA_SIGNALING","KEGG_ERBB_SIGNALING_PATHWAY","BIOCARTA_HER2_PATHWAY")

## Load general datasets
# Load subclone information from Nature Cancer 2021 paper
subconeInfo <- data.table( read.csv("~/Dropbox/FELINE Project (1)/Manuscript/Sumbission folder/Resubmission version/Acceptance revision/Third publication revisions/Source data files/Bild_SourceData_Fig5.csv") %>%
                             dplyr::select( Cell.ID,	Patient.Study.ID,	ARM,	Day,	Response,	Subclone, CDK6_mu,mapk_axisV1,mapk_axisV2, ST_JNK_MAPK_PATHWAY.ssGSEA))


# Select just the intracellular signal transduction signatures of interest.
ssgseapathways <- data.table( read.csv("~/Dropbox/Cancer_pheno_evo/data/FELINE2/SSGSEA list/Cancer SSGSEA pathway list_ab.csv"))[(!is.na(extra))|(pathway%in%wch_signalTrans),]#[!is.na(fib_ECM)| !is.na(immune)| !is.na(extra)]
# Load ssgsea data and umap landscape of cancer and normal epithelial cells
load(  file= "~/Dropbox/Cancer_pheno_evo/data/FELINE2/Cancer gene expression wrangling/cancer and normal Communication on ssgsea landscape.RData")# udataComm

PatIDunique <- unique(udataComm[Celltype== "Cancer cells"]$Patient.Study.ID)

# Load A. Nath's ER response signature
ERscoresAN_raw <- data.table(readRDS("/Users/jason/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Data_analysis/scRNA/05_ssGSEA_score/Signature_Aritro_AN_GeneSets/results/FEL001046_Cancer_cells_scRNA.zinbwave.normalized.ssGSEA_scores.Aritro_signature.RDS"))
ERscoresAN <- data.table(t(ERscoresAN_raw[,-1]),keep.rownames = T)
names(ERscoresAN) <- c("Cell.ID" , ERscoresAN_raw[,1]$'Gene Set')


ssGSEAdat <- udataComm[Celltype== "Cancer cells"] %>%
  dplyr::select(Cell.ID,orig.ident,dynamic_class3,Patient.Study.ID,Day,ARM.x,Treatment,orig.ident,key_,Celltype,one_of(c(
    names(udataComm)[grep("ERBB",names(udataComm))],names(udataComm)[grep("_EGF",names(udataComm))],
    wch_signalTrans,as.character(ssgseapathways$pathway)) ) )
PatIDunique <- unique(ssGSEAdat$Patient.Study.ID)
rm(list="udataComm")

ssGSEAdat <- merge(ssGSEAdat, ERscoresAN, by="Cell.ID")


ERBBsetlist<-c(names(ssGSEAdat)[grep("ERBB",names(ssGSEAdat))],names(ssGSEAdat)[grep("_EGF",names(ssGSEAdat) )] ) 
ERBBsetlist <- ERBBsetlist[!grepl("_DN",ERBBsetlist)]
rmthese<-c("RAY_TUMORIGENESIS_BY_ERBB2_CDC25A_UP" , "BORLAK_LIVER_CANCER_EGF_UP","REACTOME_EGFR_DOWNREGULATION","ZWANG_DOWN_BY_2ND_EGF_PULSE")
ERBBsetlist <- ERBBsetlist[ERBBsetlist!=rmthese]
corrplot::corrplot(cor( ssGSEAdat %>% dplyr::select(ERBBsetlist  ))  ,
                   hclust.method ="ward",diag=F,order="FPC",
                   type="upper",
                   tl.cex = 0.5)
# PCA
res.pca <- prcomp(ssGSEAdat %>% dplyr::select(ERBBsetlist  ), center=T,scale = TRUE)
ssGSEAdat$PCA1 <- res.pca$x[,1]
ssGSEAdat$PCA2 <- res.pca$x[,2]


ssGSEAdat[,responseLab:="Sensitive"]
ssGSEAdat[dynamic_class3=="Non-response",responseLab:="Resistant"]

ssGSEAdat[,HasDay0:=(sum(Day==0))>0, by=Patient.Study.ID]
ssGSEAdat[,HasDay14:=(sum(Day==14))>0, by=Patient.Study.ID]
ssGSEAdat[,HasDay180:=(sum(Day==180))>0, by=Patient.Study.ID]
ssGSEAdat[,nDay0:=(sum(Day==0)), by=Patient.Study.ID]
ssGSEAdat[,nDay14:=(sum(Day==14)), by=Patient.Study.ID]
ssGSEAdat[,nDay180:=(sum(Day==180)), by=Patient.Study.ID]


ggplot(ssGSEAdat,aes(y=PCA2,x=NAGASHIMA_EGF_SIGNALING_UP, col=as.factor(Day)))+geom_point(size=0.01)
ggplot(ssGSEAdat,aes(y=PCA1,x=PCA2, col=as.factor(Day)))+geom_point(size=0.01)

cor( ssGSEAdat$PCA1, ssGSEAdat %>% dplyr::select(ERBBsetlist  ))[, order( -cor( ssGSEAdat$PCA1, ssGSEAdat %>% dplyr::select(ERBBsetlist  )))] [1:20]
out1<-cor( ssGSEAdat$PCA2, ssGSEAdat %>% dplyr::select(ERBBsetlist  ))[, order( -cor( ssGSEAdat$PCA2, ssGSEAdat %>% dplyr::select(ERBBsetlist  )))] #[1:20]
#write.csv(as.data.table(out1,keep.rownames = T),file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/correlatedERBBssGSEAscores.csv")

colnames(cor( ssGSEAdat$PCA1, ssGSEAdat %>% dplyr::select(ERBBsetlist  ))) [which.max(cor( ssGSEAdat$PCA1, ssGSEAdat %>% dplyr::select(ERBBsetlist  )) )]
colnames(cor( ssGSEAdat$PCA2, ssGSEAdat %>% dplyr::select(ERBBsetlist  ))) [which.max(cor( ssGSEAdat$PCA2, ssGSEAdat %>% dplyr::select(ERBBsetlist  )) )]

# PCA 2 is associated with response to ERBB whereas PCA1 incodes signaling
#ggplot(ssGSEAdat,aes(y=-Bad_ER,x=PCA2, col=as.factor(Day)))+geom_point(size=0.01)
#ggplot(ssGSEAdat,aes(y=-Bad_ER,x=PCA2, col=dynamic_class3))+theme_classic()+geom_point(size=0.01)+facet_grid(dynamic_class3~paste0("Day ",Day))+theme(aspect.ratio=1)


ggplot(ssGSEAdat,aes(y=Bad_ER, x=PCA2))+# ,col=HALLMARK_ESTROGEN_RESPONSE_EARLY))+
  theme_classic(base_size=26)+
  geom_point(size=0.1,alpha=0.5)+
  theme(aspect.ratio=1)+
  labs(y="Endorse empirical signature \n (predicts poor survival on endocrine therapy)", x="ERBB family pathway activation \n (Composite ERBB response signature)")+
  facet_wrap(~paste0("Day ",Day))
ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Endorse vs ERBB byDay.pdf",width=12, height=8.5)

ggplot(ssGSEAdat,aes(y=Bad_ER, x=HALLMARK_ESTROGEN_RESPONSE_EARLY))+# ,col=HALLMARK_ESTROGEN_RESPONSE_EARLY))+
  theme_classic(base_size=26)+
  geom_point(size=0.1,alpha=0.5)+
  theme(aspect.ratio=1)+
  labs(y="Endorse empirical signature \n (predicts poor survival on endocrine therapy)", x="Hallmark estrogen \n response early")+
  facet_wrap(~paste0("Day ",Day))
ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Endorse vs Hallmark estrogen early byDay.pdf",width=12, height=8.5)

ggplot(ssGSEAdat,aes( x=HALLMARK_ESTROGEN_RESPONSE_EARLY,group=Day,fill=as.factor(Day)))+
  theme_classic(base_size=26)+
  geom_density(size=0.1, alpha=0.6,col=NA)+
  theme(aspect.ratio=1)+
  labs(y="Frequency density", x="Hallmark estrogen \n response early")+
  scale_color_manual(name="Day",values=(pal_aaas("default")(3))[c(3,1,2)])  +
  scale_fill_manual(name="Day",values=(pal_aaas("default")(3))[c(3,1,2)])  

ggplot(ssGSEAdat,aes( y=HALLMARK_ESTROGEN_RESPONSE_EARLY,x=log(1+Day),fill=as.factor(Day)))+
  theme_classic(base_size=26)+
  geom_violin(size=0.1, alpha=0.6,col=NA)+
  theme(aspect.ratio=1, legend.position="none")+
  labs(x="Day", y="Estrogen pathway activation \n (Hallmark estrogen response early)")+
  scale_color_manual(name="Day",values=(pal_aaas("default")(3))[c(3,1,2)])  +
  scale_fill_manual(name="Day",values=(pal_aaas("default")(3))[c(3,1,2)])  +
  scale_x_continuous(breaks=log(1+c(0,14,180)),labels=(c(0,14,180)))
ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Hallmark estrogen early by Day violin.pdf",width=12, height=8.5)


ggplot(ssGSEAdat[HasDay0==T& HasDay180==T],aes( x=HALLMARK_ESTROGEN_RESPONSE_EARLY,group=interaction(Patient.Study.ID,Day),fill=interaction(Patient.Study.ID,Day)))+
  theme_classic(base_size=26)+
  geom_density(bw=0.01,size=0.1,alpha=0.5)+
  theme(aspect.ratio=1,legend.position="none")+
  labs(y="", x="Hallmark estrogen \n response early")+facet_wrap(~Patient.Study.ID,ncol=5)


### Pre-post plots
ggplot(ssGSEAdat[Day!=14],aes(y=Bad_ER, x=PCA2, col=(Day==180) ))+theme_classic(base_size=26)+
  geom_point(size=0.1,alpha=0.5)+
  theme(aspect.ratio=1)+
  labs(y="Endorse empirical signature \n (predicts poor survival on endocrine therapy)", x="ERBB family pathway activation \n (Composite ERBB response signature)")+
  scale_color_manual(name="Time point", values=c("black",(pal_aaas("default")(3))[2]) ,labels=c("pre-treatment", "post-treatment"))+
  geom_point(data=ssGSEAdat[Day==180],size=0.1,alpha=0.5)
ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Endorse vs ERBB response prevspost.pdf",width=12, height=8.5)
ggplot(ssGSEAdat[Day!=14],aes(y=HALLMARK_ESTROGEN_RESPONSE_EARLY, x=PCA2, col=(Day==180) ))+theme_classic(base_size=26)+
  geom_point(size=0.1,alpha=0.5)+
  theme(aspect.ratio=1)+
  labs(y="Estrogen pathway activation \n (Hallmark estrogen response early)", x="ERBB family pathway activation \n (Composite ERBB response signature)")+
  scale_color_manual(name="Time point", values=c("black",(pal_aaas("default")(3))[2]) ,labels=c("pre-treatment", "post-treatment"))+
  geom_point(data=ssGSEAdat[Day==180],size=0.1,alpha=0.5)
ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Estrogen vs ERBB activity prevspost.pdf",width=12, height=8.5)


# All timepoints&response
ggplot(ssGSEAdat,aes(y=Bad_ER,x=PCA2, col=dynamic_class3))+ theme_classic(base_size=26)+ geom_point(size=0.01)+ facet_grid(responseLab~paste0("Day ",Day))+
  theme(legend.position="none",aspect.ratio=1)+
  labs(y="Endorse empirical signature \n (predicts poor survival on endocrine therapy)", x="ERBB family pathway activation \n (Composite ERBB response signature)")+
  scale_color_npg()
ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Endorse vs ERBB response over time and by tumor response.pdf",width=12, height=8.5)

ggplot(ssGSEAdat,aes(y=HALLMARK_ESTROGEN_RESPONSE_EARLY,x=PCA2, col=dynamic_class3))+ theme_classic(base_size=26)+ geom_point(size=0.01)+ facet_grid(responseLab~paste0("Day ",Day))+
  theme(legend.position="none",aspect.ratio=1)+
  labs(y="Estrogen pathway activation \n (Hallmark estrogen response early)", x="ERBB family pathway activation \n (Composite ERBB response signature)")+
  scale_color_npg()
ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Estrogen vs ERBB pathway activity over time and by tumor response.pdf",width=12, height=8.5)


# by patient
ggplot(ssGSEAdat[HasDay0==T& HasDay180==T],aes(y=Bad_ER, x=PCA2,col=as.factor(Day) ))+theme_classic(base_size=26)+
  geom_point(size=0.2,alpha=0.5)+
  facet_wrap(~Patient.Study.ID,ncol=5)+theme(aspect.ratio=1)+scale_color_manual(name="Day",values=(pal_aaas("default")(3))[c(3,1,2)])+
  geom_point(data=ssGSEAdat[HasDay0==T& HasDay180==T][Day==180],size=0.2,alpha=0.5)+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  labs(y="Endorse empirical signature \n (predicts poor survival on endocrine therapy)", x="ERBB family pathway activation \n (Composite ERBB response signature)")
ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Endorse vs ERBB response over time by tumor.pdf",width=20, height=20)

ggplot(ssGSEAdat[HasDay0==T& HasDay180==T],aes(y=HALLMARK_ESTROGEN_RESPONSE_EARLY, x=PCA2,col=as.factor(Day) ))+theme_classic(base_size=26)+
  geom_point(size=0.2,alpha=0.5)+
  facet_wrap(~Patient.Study.ID,ncol=5)+theme(aspect.ratio=1)+scale_color_manual(name="Day",values=(pal_aaas("default")(3))[c(3,1,2)])+
  geom_point(data=ssGSEAdat[HasDay0==T& HasDay180==T][Day==180],size=0.2,alpha=0.5)+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  labs(y="Estrogen pathway activation \n (Hallmark estrogen response early)", x="ERBB family pathway activation \n (Composite ERBB response signature)")
ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Estrogen vs ERBB pathway activity over time by tumor.pdf",width=20, height=20)


# All timepoints
ggplot(ssGSEAdat[],aes(y=Bad_ER, x=PCA2,col=as.factor(Day) ))+theme_classic(base_size=26)+
  geom_point(size=0.2,alpha=0.5)+
  theme(aspect.ratio=1)+scale_color_manual(name="Day",values=(pal_aaas("default")(3))[c(3,1,2)])+
  geom_point(data=ssGSEAdat[][Day==180],size=0.2,alpha=0.5)+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  labs(y="Endorse empirical signature \n (predicts poor survival on endocrine therapy)", x="ERBB family pathway activation \n (Composite ERBB response signature)")
ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Endorse vs ERBB response over time.pdf",width=20, height=20)

ggplot(ssGSEAdat,aes(y=HALLMARK_ESTROGEN_RESPONSE_EARLY, x=PCA2,col=as.factor(Day) ))+theme_classic(base_size=26)+
  geom_point(size=0.2,alpha=0.5)+
  theme(aspect.ratio=1)+scale_color_manual(name="Day",values=(pal_aaas("default")(3))[c(3,1,2)])+
  geom_point(data=ssGSEAdat[Day==180],size=0.2,alpha=0.5)+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  labs(y="Estrogen pathway activation \n (Hallmark estrogen response early)", x="ERBB family pathway activation \n (Composite ERBB response signature)")
ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Estrogen vs ERBB pathway activity over time.pdf",width=20, height=20)

ggplot(ssGSEAdat[Day==0],aes(y=HALLMARK_ESTROGEN_RESPONSE_EARLY, x=PCA2,col=as.factor(Day) ))+theme_classic(base_size=26)+
  geom_point(size=0.2,alpha=0.5)+
  theme(legend.position="none",aspect.ratio=1)+scale_color_manual(name="Day",values=(pal_aaas("default")(3))[c(3,1,2)])+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  facet_wrap(~paste0("Day ",Day), ncol=1)+lims(y=c(-0.11,0.3), x=c(-10,15))+
  labs(y="Estrogen pathway activation \n (Hallmark estrogen response early)", x="ERBB family pathway activation \n (Composite ERBB response signature)")
ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Estrogen vs ERBB pathway activity Day0.pdf",width=8, height=8)

ggplot(ssGSEAdat[Day==14],aes(y=HALLMARK_ESTROGEN_RESPONSE_EARLY, x=PCA2,col=as.factor(Day) ))+theme_classic(base_size=26)+
  geom_point(size=0.2,alpha=0.5)+
  theme(legend.position="none",aspect.ratio=1)+scale_color_manual(name="Day",values=(pal_aaas("default")(3))[c(1,2)])+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  facet_wrap(~paste0("Day ",Day), ncol=1)+lims(y=c(-0.11,0.3), x=c(-10,15))+
  labs(y="Estrogen pathway activation \n (Hallmark estrogen response early)", x="ERBB family pathway activation \n (Composite ERBB response signature)")
ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Estrogen vs ERBB pathway activity Day14.pdf",width=8, height=8)


ggplot(ssGSEAdat[Day%in%c(0,14)],aes(y=HALLMARK_ESTROGEN_RESPONSE_EARLY, x=PCA2,col=as.factor(Day) ))+theme_classic(base_size=26)+
  geom_point(size=0.2,alpha=0.5)+
  geom_point(data=ssGSEAdat[Day==0],size=0.12,alpha=0.05)+
  theme(aspect.ratio=1)+scale_color_manual(name="Day",values=(pal_aaas("default")(3))[c(3,1,2)])+
  guides(colour = guide_legend(override.aes = list(size=5,alpha=0.5)))+
  lims(y=c(-0.11,0.3), x=c(-10,15))+
  labs(y="Estrogen pathway activation \n (Hallmark estrogen response early)", x="ERBB family pathway activation \n (Composite ERBB response signature)")
ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Estrogen vs ERBB pathway activity Day0and14.pdf",width=8, height=8)


ggplot(ssGSEAdat[Day%in%c(0,14,180)][],aes(y=HALLMARK_ESTROGEN_RESPONSE_EARLY, x=PCA2
                                     #,col=as.factor(Day)
                                     ))+theme_classic(base_size=26)+
  theme(aspect.ratio=1)+
  facet_wrap(~paste0("Day ",Day), ncol=3)+
  lims(y=c(-0.11,0.3), x=c(-10,15))+
  stat_bin2d(aes(fill = log(1+..count..),col=log(1+..count..)), bins = 50) + 
  scale_fill_gradient(name="Cell count",low = "lightblue", high = "red", breaks=log(1+c(1,2,4,8,16,32,64,128,256,512)),labels=c(1,2,4,8,16,32,64,128,256,512))+
  scale_color_gradient(name="Cell count",low = "lightblue", high = "red", breaks=log(1+c(1,2,4,8,16,32,64,128,256,512)),labels=c(1,2,4,8,16,32,64,128,256,512))+
  labs(y="Estrogen pathway activation \n (Hallmark estrogen response early)", x="ERBB family pathway activation \n (Composite ERBB response signature)")+
  guides(colour=guide_colourbar(barheight = 15),fill=guide_colourbar(barheight = 15))
ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Estrogen vs ERBB pathway activity Day0to180 heatmap.pdf",width=16, height=8)


ggplot(ssGSEAdat[Day%in%c(0,14,180)][Day%in%c(0,14)][],aes(y=HALLMARK_ESTROGEN_RESPONSE_EARLY, x=PCA2
                                           #,col=as.factor(Day)
))+theme_classic(base_size=26)+
  theme(aspect.ratio=1)+
  facet_grid(paste0("Day ",Day)~TumorResponse)+
  lims(y=c(-0.11,0.3), x=c(-10,15))+
  stat_bin2d(aes(fill = log(1+..count..),col=log(1+..count..)), bins = 100) + 
  scale_fill_gradient(name="Cell count",low = "lightblue", high = "red", breaks=log(1+c(1,2,4,8,16,32,64,128,256,512)),labels=c(1,2,4,8,16,32,64,128,256,512))+
  scale_color_gradient(name="Cell count",low = "lightblue", high = "red", breaks=log(1+c(1,2,4,8,16,32,64,128,256,512)),labels=c(1,2,4,8,16,32,64,128,256,512))+
  labs(y="Estrogen pathway activation \n (Hallmark estrogen response early)", x="ERBB family pathway activation \n (Composite ERBB response signature)")+
  guides(colour=guide_colourbar(barheight = 15),fill=guide_colourbar(barheight = 15))
#ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Estrogen vs ERBB pathway activity Day0to14 heatmap by response.pdf",width=16, height=8)

ggplot(ssGSEAdat[Day%in%c(0,14,180)][Day%in%c(14)][],aes(y=HALLMARK_ESTROGEN_RESPONSE_EARLY, x=PCA2
                                                           #,col=as.factor(Day)
))+theme_classic(base_size=26)+
  theme(aspect.ratio=1)+
  facet_grid(paste0("Day ",Day)~TumorResponse)+
  lims(y=c(-0.11,0.3), x=c(-10,15))+
  stat_bin2d(aes(fill = log(1+..count..),col=log(1+..count..)), bins = 100) + 
  scale_fill_gradient(name="Cell count",low = "lightblue", high = "red", breaks=log(1+c(1,2,4,8,16,32,64,128,256,512)),labels=c(1,2,4,8,16,32,64,128,256,512))+
  scale_color_gradient(name="Cell count",low = "lightblue", high = "red", breaks=log(1+c(1,2,4,8,16,32,64,128,256,512)),labels=c(1,2,4,8,16,32,64,128,256,512))+
  labs(y="Estrogen pathway activation \n (Hallmark estrogen response early)", x="ERBB family pathway activation \n (Composite ERBB response signature)")+
  guides(colour=guide_colourbar(barheight = 15),fill=guide_colourbar(barheight = 15))
#ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Estrogen vs ERBB pathway activity Day14 heatmap by response.pdf",width=16, height=8)

# just be response 
ggplot(ssGSEAdat[responseLab!="Sensitive"],aes(y=Bad_ER, x=PCA2,col=responseLab ))+theme_classic(base_size=26)+
  geom_point(size=0.2,alpha=0.5)+
  theme(aspect.ratio=1)+scale_color_npg(name="Tumor response")+
  geom_point(data=ssGSEAdat[responseLab=="Sensitive"][],size=0.2,alpha=0.5)+
  facet_wrap(~paste0("Day ",Day))+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  labs(y="Endorse empirical signature \n (predicts poor survival on endocrine therapy)", x="ERBB family pathway activation \n (Composite ERBB response signature)")
#ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Endorse vs ERBB response by tumor response over time.pdf",width=20, height=20)

ggplot(ssGSEAdat[responseLab!="Sensitive"],aes(y=HALLMARK_ESTROGEN_RESPONSE_EARLY, x=PCA2,col=responseLab ))+theme_classic(base_size=26)+
  geom_point(size=0.2,alpha=0.5)+
  theme(aspect.ratio=1)+scale_color_npg(name="Tumor response")+
  geom_point(data=ssGSEAdat[responseLab=="Sensitive"][],size=0.2,alpha=0.5)+
  facet_wrap(~paste0("Day ",Day))+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  labs(y="Estrogen pathway activation \n (Hallmark estrogen response early)", x="ERBB family pathway activation \n (Composite ERBB response signature)")
#ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Estrogen vs ERBB response by tumor response over time.pdf",width=20, height=20)




medianEstrogen <- ssGSEAdat$HALLMARK_ESTROGEN_RESPONSE_EARLY%>%median()
medianERBB <- ssGSEAdat$PCA2%>%median()
ssGSEAdat[,EST_high:= HALLMARK_ESTROGEN_RESPONSE_EARLY > medianEstrogen]
ssGSEAdat[,ERB_high:= PCA2 > medianERBB]
ssGSEAdat[,TumorResponse:="Resistant"]
ssGSEAdat[dynamic_class3=="Response",TumorResponse:="Sensitive"]
ssGSEAdat[,Treat:="Combination ribociclib"]
ssGSEAdat[ARM.x=="A",Treat:="Letrozole alone"]

ggplot(ssGSEAdat[],aes(y=PCA2, x=EST_high,col=TumorResponse,fill=TumorResponse, group=interaction(EST_high,Day,TumorResponse) ))+theme_classic(base_size=26)+
  geom_violin(size=0.2,alpha=0.5)+
  geom_jitter(position=position_dodge(width=0.85), size=0.2,alpha=0.5)+
  facet_grid(Treat~Day)+
  theme(aspect.ratio=1)+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  scale_x_discrete(labels=c("Low","High"))+
  labs(y="ERBB family pathway activation \n (Composite ERBB response signature)",x="Estrogen pathway activity")+
  scale_color_npg()+scale_fill_npg()




prop_EST_highdat<-data.table(ssGSEAdat%>%group_by(Patient.Study.ID,dynamic_class3,ARM.x,Day)%>%
                               summarise(prop_EST_high= sum(EST_high)/ length(EST_high),
                                         prop_ERB_high= sum(ERB_high)/ length(ERB_high),
                                         mean_ERB_high= mean(ERB_high)))
prop_EST_highdat[,prop_EST_low:=1-prop_EST_high]
prop_EST_highdat[,prop_ERB_low:=1-prop_ERB_high]
prop_EST_highdat[, initERB:= ((Day==0)*mean_ERB_high), by=Patient.Study.ID ]
prop_EST_highdat[, change:=mean_ERB_high -initERB]

erbphenodat<-ssGSEAdat[HasDay0==T&HasDay14==T&HasDay180==T]#[Day%in%c(0,180)]
prop_EST_highdat1<-data.table(erbphenodat%>%group_by(Patient.Study.ID,dynamic_class3,ARM.x,Day)%>%
                               summarise(mean_ERB_high= median(PCA2),
                                         N=n() ))
prop_EST_highdat1[, initERB:= sum((Day==0)*mean_ERB_high), by=Patient.Study.ID ]
prop_EST_highdat1[, change:=mean_ERB_high -initERB]
erbphenodat$Patient.Study.ID <- factor(erbphenodat$Patient.Study.ID , levels= prop_EST_highdat1[Day==180][order(-change)]$Patient.Study.ID)
prop_EST_highdat1$Patient.Study.ID <- factor(prop_EST_highdat1$Patient.Study.ID , levels= prop_EST_highdat1[Day==180][order(-change)]$Patient.Study.ID)
prop_EST_highdat1[order(Patient.Study.ID,-change)]
ggplot(erbphenodat,
       aes(y=PCA2,x=Patient.Study.ID , fill=as.factor(Day),group=interaction(Patient.Study.ID,Day) ))+
  theme_classic(base_size=26)+ geom_boxplot()+stat_boxplot(geom="errorbar")+
  geom_point(size=0.1,alpha=0.1,position=position_dodge(width=0.8))+
  theme(aspect.ratio=.5)+
  scale_color_manual(name="Timepoint",values=(pal_aaas("default")(3))[c(3,1,2)], labels=c("pre-treatment","followup","post-treatment"))  +
  scale_fill_manual(name="Timepoint",values=(pal_aaas("default")(3))[c(3,1,2)],labels=c("pre-treatment","followup","post-treatment"))  +
  labs(y="ERBB family pathway activation \n (Composite ERBB response signature)",
       x="Patient")+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/ERBB phenotype shift pre to post treatment1.pdf",width=14, height=10)
ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/ERBB phenotype shift pre early and post treatment1.pdf",width=14, height=10)

statsdt<-data.table( coef( summary( lm(PCA2~Patient.Study.ID+as.factor(Day):Patient.Study.ID, erbphenodat) ) ) , keep.rownames = T)[grep("Day)180",rn)][order(-Estimate)]
write.csv( statsdt , file="~/Dropbox/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/ERBB phenotype shift pre to post treatment stats.csv")

erbphenodat <- ssGSEAdat[nDay0>20&nDay180>20][Day%in%c(0,180)]
prop_EST_highdat1 <- data.table(erbphenodat%>%group_by(Patient.Study.ID,dynamic_class3,ARM.x,Day)%>%
                                summarise(mean_ERB_high= median(PCA2),
                                          medan_NAGASHERB_high= median(NAGASHIMA_EGF_SIGNALING_UP),
                                          N=n() ))
prop_EST_highdat1[, initERB:= sum((Day==0)*mean_ERB_high), by=Patient.Study.ID ]
prop_EST_highdat1[, change:=mean_ERB_high -initERB]
erbphenodat$Patient.Study.ID <- factor(erbphenodat$Patient.Study.ID , levels= prop_EST_highdat1[Day==180][order(-change)]$Patient.Study.ID)
prop_EST_highdat1$Patient.Study.ID <- factor(prop_EST_highdat1$Patient.Study.ID , levels= prop_EST_highdat1[Day==180][order(-change)]$Patient.Study.ID)
prop_EST_highdat1[order(Patient.Study.ID,-change)]

ERBB_mu_save <- data.table(erbphenodat%>%group_by(Patient.Study.ID,dynamic_class3,ARM.x,Day)%>%
                                  summarise(median_ERB_high= median(PCA2),
                                            median_NAGASHERB_high= median(NAGASHIMA_EGF_SIGNALING_UP),
                                            N=n() ))

#save( ERBB_mu_save , res.pca,      ERBBsetlist,file="~/Dropbox/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/ERBB median phenotype shift pre to post treatment.RData")

discPlot <- ggplot(erbphenodat,
       aes(y=PCA2,x=Patient.Study.ID , fill=as.factor(Day),group=interaction(Patient.Study.ID,Day) ))+
  theme_classic(base_size=26)+ geom_boxplot()+stat_boxplot(geom="errorbar")+
  geom_jitter(size=0.4,alpha=0.6,position=position_dodge(width=0.8))+
  theme(aspect.ratio=.5)+
  scale_color_manual(name="Timepoint",values=c("yellow","red"), labels=c("pre-treatment","post-treatment"))  +
  scale_fill_manual(name="Timepoint",values=c("yellow","red"),labels=c("pre-treatment","post-treatment"))  +
  labs(y="ERBB family pathway activation \n (Composite ERBB response signature)",
       x="Patient")+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/ERBB phenotype shift pre to post treatment1.pdf",width=14, height=10)
ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/ERBB phenotype shift pre and post treatment Rena.pdf",width=14, height=10)

discData <- data.table( erbphenodat%>%dplyr::select(Cell.ID,Patient.Study.ID,dynamic_class3,Day,ARM.x,Treatment,Celltype,PCA2,
                            responseLab) )
setnames( discData , old=c("PCA2","ARM.x"), new=c("ERBBaxis","ARM"))

#write.csv(discData,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/DiscoveryCohortERBBPrePost.csv")
#summary( lm(ERBBaxis~ -1+Patient.Study.ID+ Day : Patient.Study.ID ,data=fread(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/Figure2/DiscoveryCohortERBBPrePost.csv")))
#summary( lmer(ERBBaxis~ Day +(1+Day|Patient.Study.ID) ,data=fread(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/Figure2/DiscoveryCohortERBBPrePost.csv")))
summary( lmer(ERBBaxis~ Day +(1+Day|Patient.Study.ID) ,data=discData))
Mod1 <- lmer(ssGSEA ~ -1+dynamic_class3*day_fact*ARM + (1+day_fact|Patient.Study.ID), REML=FALSE,data= u_dat) 
summary(Mod1)
summary( lm(ERBBaxis~ -1+Patient.Study.ID+ Day : Patient.Study.ID,data=discData))
#summary( lmer(ERBBaxis~ (Day+dynamic_class3) +(1+Day|Patient.Study.ID) ,data=fread(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/Figure2/ValidationCohortERBBPrePost.csv")))

ggplot(erbphenodat,
       aes(y=PCA2,x=Patient.Study.ID , fill=as.factor(Day),group=interaction(Patient.Study.ID,Day) ))+
  theme_classic(base_size=26)+ geom_boxplot()+stat_boxplot(geom="errorbar")+
  geom_jitter(size=0.4,alpha=0.6,position=position_dodge(width=0.8))+
  theme(aspect.ratio=.5)+
  scale_color_manual(name="Timepoint",values=c("yellow","red"), labels=c("pre-treatment","post-treatment"))  +
  scale_fill_manual(name="Timepoint",values=c("yellow","red"),labels=c("pre-treatment","post-treatment"))  +
  labs(y="ERBB family pathway activation \n (Composite ERBB response signature)",
       x="Patient")+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#



ggplot(erbphenodat,
              aes(y=HALLMARK_ESTROGEN_RESPONSE_EARLY,x=Patient.Study.ID , fill=as.factor(Day),group=interaction(Patient.Study.ID,Day) ))+
       theme_classic(base_size=26)+ geom_boxplot()+stat_boxplot(geom="errorbar")+
       geom_point(size=0.1,alpha=0.1,position=position_dodge(width=0.8))+
       theme(aspect.ratio=.5)+
       scale_color_manual(name="Timepoint",values=(pal_aaas("default")(3))[c(1,2)], labels=c("pre-treatment","post-treatment"))  +
       scale_fill_manual(name="Timepoint",values=(pal_aaas("default")(3))[c(1,2)],labels=c("pre-treatment","post-treatment"))  +
       labs(y="Estrogen pathway activation \n (Hallmark estrogen response early)",               x="Patient")+ 
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Estrogen phenotype shift pre to post treatment1.pdf",width=14, height=10)




ggplot( prop_EST_highdat, aes(y=log(prop_EST_high),x=log(1+Day) ))+geom_point()+#geom_violin(aes( group=Day))+
  geom_point()+theme_classic(base_size=26)+theme(aspect.ratio=1)+
  geom_smooth(size=2)+coord_trans(y="exp")+
  scale_y_continuous(name= "Proportion cells with \n high estrogen pathway activity",breaks=log(c(0:5)/5), labels=c(0:5)/5)+
  scale_x_continuous(name="Day",breaks=log(1+c(0,14,180)), labels=c(0,14,180))
ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Estrogen percent cells high activity over time.pdf",width=8, height=8)


ggplot( prop_EST_highdat, aes(y=log(prop_ERB_high), x=log(1+Day) ,col=prop_EST_high))+geom_point()+#geom_violin(aes( group=Day))+
  geom_jitter(width=0.25)+theme_classic(base_size=26)+theme(aspect.ratio=1)+
  geom_smooth(size=2)+coord_trans(y="exp")+
  scale_y_continuous(name= "Proportion cells with \n high ERBB pathway activity",breaks=log(c(0:5)/5), labels=c(0:5)/5)+
  scale_x_continuous(name="Day",breaks=log(1+c(0,14,180)), labels=c(0,14,180))+facet_wrap(Treat~TumorResponse)
ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/ERBB percent cells high activity over time by treat and outcome.pdf",width=8, height=8)

ggplot( prop_EST_highdat, aes(y=(prop_EST_high),x=log(1+Day) ,col=prop_ERB_high))+geom_point()+#geom_violin(aes( group=Day))+
  geom_jitter(width=0.25)+theme_classic(base_size=26)+theme(aspect.ratio=1)+
  geom_smooth(method="glm",col="black",family="binomial",size=2)+coord_trans(y="exp")+
  scale_color_viridis_c(name= "Proportion cells with \n high ERBB \npathway activity")+
  coord_trans(y="exp")+
  scale_y_continuous(name= "Proportion cells with \n high estrogen pathway activity",breaks=(c(0:5)/5), labels=c(0:5)/5)+
  scale_x_continuous(name="Day",breaks=log(1+c(0,14,180)), labels=c(0,14,180))+facet_wrap(Treat~TumorResponse)
ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/ESTROGEN and ERBB percent cells high activity over time by treat and outcome.pdf",width=12, height=12)

prop_EST_highdat[,Treat:="Combination ribociclib"]
prop_EST_highdat[ARM.x=="A",Treat:="Letrozole alone"]
prop_EST_highdat[,TumorResponse:="Resistant"]
prop_EST_highdat[dynamic_class3=="Response",TumorResponse:="Sensitive"]

ggplot( prop_EST_highdat[Treat=="Letrozole alone"], aes(y=log(prop_EST_high),x=log(1+Day)  , col=TumorResponse,fill=TumorResponse ))+geom_point()+#geom_violin(aes( group=Day))+
  geom_point()+theme_classic(base_size=26)+theme(aspect.ratio=1)+
  geom_smooth(size=2,se=F)+coord_trans(y="exp")+
  scale_y_continuous(name= "Proportion cells with \n high estrogen pathway activity",breaks=log(c(0:5)/5), labels=c(0:5)/5)+
  scale_x_continuous(name="Day",breaks=log(1+c(0,14,180)), labels=c(0,14,180))+
  scale_color_npg(name="Tumor response")+
  scale_fill_npg(name="Tumor response")
ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Estrogen percent cells high activity over time by response.pdf",width=8, height=8)

