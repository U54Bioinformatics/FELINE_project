rm(list=ls())
require(mgcv);require(data.table);require(dplyr);require(ggplot2);require(tidyr);require(lme4);require(lmerTest);require(parallel)
library(effects);require(umap);require(Rfast);require(ider);library("dendextend");library(ggdendro);require(ggsci);require(viridis);require("Rdimtools")
require(merTools)

# Define data location
SourceDataLoc<-"/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript Feline ERBB facilitation/Molecular Systems Biology/Revision and submission folder/Source Data/"
Intermediateloc<-paste0(SourceDataLoc,"Figure 2/")

# Load discovery data and ERBB list
Discovery_inputdata <- data.table(read.csv(file=paste0(Intermediateloc,"SourceData_Figure_2A_DiscoveryCancerERBBssGSEA.csv") ))

# Load validation data
Validation_inputdata <- data.table(read.csv(file=paste0(Intermediateloc,"SourceData_Figure_2A_ValidationCancerERBBssGSEA.csv") ))

# List erbb pathways
ERBBsetlist <- c("KEGG_ERBB_SIGNALING_PATHWAY","LANDIS_ERBB2_BREAST_PRENEOPLASTIC_UP",
                 "LANDIS_ERBB2_BREAST_TUMORS_324_UP","LANDIS_ERBB2_BREAST_TUMORS_65_UP" ,               
                 "PEDERSEN_METASTASIS_BY_ERBB2_ISOFORM_1","PEDERSEN_METASTASIS_BY_ERBB2_ISOFORM_4" ,         
                 "PEDERSEN_METASTASIS_BY_ERBB2_ISOFORM_6","PEDERSEN_METASTASIS_BY_ERBB2_ISOFORM_7",          
                 "PEDERSEN_TARGETS_OF_611CTF_ISOFORM_OF_ERBB2","PID_ERBB1_DOWNSTREAM_PATHWAY",                    
                 "PID_ERBB1_INTERNALIZATION_PATHWAY","PID_ERBB1_RECEPTOR_PROXIMAL_PATHWAY" ,            
                 "PID_ERBB2_ERBB3_PATHWAY","PID_ERBB4_PATHWAY"  ,                            
                 "RAY_TUMORIGENESIS_BY_ERBB2_CDC25A_UP","REACTOME_GRB2_EVENTS_IN_ERBB2_SIGNALING" ,        
                 "REACTOME_NUCLEAR_SIGNALING_BY_ERBB4","REACTOME_PI3K_EVENTS_IN_ERBB2_SIGNALING",         
                 "REACTOME_PI3K_EVENTS_IN_ERBB4_SIGNALING","REACTOME_SHC1_EVENTS_IN_ERBB4_SIGNALING" ,       
                 "REACTOME_SIGNALING_BY_ERBB2","REACTOME_SIGNALING_BY_ERBB4",                     
                 "SMID_BREAST_CANCER_ERBB2_UP","AMIT_EGF_RESPONSE_120_HELA",                      
                 "AMIT_EGF_RESPONSE_120_MCF10A","AMIT_EGF_RESPONSE_240_HELA",                      
                 "AMIT_EGF_RESPONSE_240_MCF10A","AMIT_EGF_RESPONSE_40_HELA" ,                      
                 "AMIT_EGF_RESPONSE_40_MCF10A","AMIT_EGF_RESPONSE_480_HELA",                      
                 "AMIT_EGF_RESPONSE_480_MCF10A","AMIT_EGF_RESPONSE_60_HELA",                       
                 "AMIT_EGF_RESPONSE_60_MCF10A","BIOCARTA_EGF_PATHWAY" ,                           
                 "BORLAK_LIVER_CANCER_EGF_UP","KOBAYASHI_EGFR_SIGNALING_24HR_UP",                
                 "NAGASHIMA_EGF_SIGNALING_UP","REACTOME_EGFR_DOWNREGULATION" ,                   
                 "REACTOME_SHC1_EVENTS_IN_EGFR_SIGNALING","REACTOME_SIGNALING_BY_CONSTITUTIVELY_ACTIVE_EGFR",
                 "REACTOME_SIGNALING_BY_EGFR_IN_CANCER","ZWANG_CLASS_1_TRANSIENTLY_INDUCED_BY_EGF",        
                 "ZWANG_CLASS_2_TRANSIENTLY_INDUCED_BY_EGF","ZWANG_CLASS_3_TRANSIENTLY_INDUCED_BY_EGF" ,       
                 "ZWANG_DOWN_BY_2ND_EGF_PULSE","ZWANG_EGF_INTERVAL_UP"  ,                         
                 "ZWANG_EGF_PERSISTENTLY_UP","ZWANG_TRANSIENTLY_UP_BY_2ND_EGF_PULSE_ONLY")

# explore correaltions
corrplot::corrplot(cor( Discovery_inputdata %>% dplyr::select(ERBBsetlist  ))  ,
                   hclust.method ="ward",diag=F,order="FPC",
                   type="upper",
                   tl.cex = 0.5)

# PCA
res.pca <- prcomp(Discovery_inputdata %>% dplyr::select(ERBBsetlist  ), center=T,scale = TRUE)
Discovery_inputdata$PCA1 <- res.pca$x[,1]
Discovery_inputdata$PCA2 <- res.pca$x[,2]

# label outcomes and indicate if a patient sample has cancer cells at a given timepoint 
Discovery_inputdata[,responseLab:="Sensitive"]
Discovery_inputdata[dynamic_class3=="Non-response",responseLab:="Resistant"]
Discovery_inputdata[,HasDay0:=(sum(Day==0))>0, by=Patient.Study.ID]
Discovery_inputdata[,HasDay14:=(sum(Day==14))>0, by=Patient.Study.ID]
Discovery_inputdata[,HasDay180:=(sum(Day==180))>0, by=Patient.Study.ID]
Discovery_inputdata[,nDay0:=(sum(Day==0)), by=Patient.Study.ID]
Discovery_inputdata[,nDay14:=(sum(Day==14)), by=Patient.Study.ID]
Discovery_inputdata[,nDay180:=(sum(Day==180)), by=Patient.Study.ID]

# UMAP and scatterplots showing correlation of pathways to principle components
#set.seed(123)
#umapEstGF1 <- umap::umap( Discovery_inputdata %>% dplyr::select(ERBBsetlist  ),n_components=1 )
#plot(umapEstGF1$layout)
#Discovery_inputdata$UMAP1 <- umapEstGF1$layout[,1]
#Discovery_inputdata$UMAP2 <- umapEstGF1$layout[,2]
#ggplot(Discovery_inputdata,aes(y=PCA1,x=UMAP2, col=as.factor(Day)))+geom_point(size=0.01)
#ggplot(Discovery_inputdata,aes(y=PCA2,x=UMAP1, col=as.factor(Day)))+geom_point(size=0.01)
ggplot(Discovery_inputdata,aes(y=PCA2,x=NAGASHIMA_EGF_SIGNALING_UP, col=as.factor(Day)))+geom_point(size=0.01)
ggplot(Discovery_inputdata,aes(y=PCA2,x=PEDERSEN_METASTASIS_BY_ERBB2_ISOFORM_4, col=as.factor(Day)))+geom_point(size=0.01)
cor( Discovery_inputdata$PCA2, Discovery_inputdata %>% dplyr::select(ERBBsetlist  ))[, order( -cor( Discovery_inputdata$PCA2, Discovery_inputdata %>% dplyr::select(ERBBsetlist  )))] [1:20]

# Extract data with paried samples with sufficient cells
erbphenodat <- Discovery_inputdata[nDay0>20&nDay180>20][Day%in%c(0,180)]

# calculate median ERBB activity to order patients on plot
ERBB_mu_save <- data.table(erbphenodat%>%group_by(Patient.Study.ID,dynamic_class3,ARM.x,Day)%>%
                             summarise(median_ERB_high= median(PCA2),
                                       median_NAGASHERB_high= median(NAGASHIMA_EGF_SIGNALING_UP),
                                       N=n() ))

ERBB_mu_save[, initERB:= sum((Day==0)*median_ERB_high), by=Patient.Study.ID ]
ERBB_mu_save[, change:=median_ERB_high -initERB]
erbphenodat$Patient.Study.ID <- factor(erbphenodat$Patient.Study.ID , levels= ERBB_mu_save[Day==180][order(-change)]$Patient.Study.ID)
ERBB_mu_save$Patient.Study.ID <- factor(ERBB_mu_save$Patient.Study.ID , levels= ERBB_mu_save[Day==180][order(-change)]$Patient.Study.ID)
ERBB_mu_save[order(Patient.Study.ID,-change)]
#save( ERBB_mu_save , res.pca,      ERBBsetlist,file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/ERBB median phenotype shift pre to post treatment.RData")

# plot results
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
#ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/ERBB phenotype shift pre and post treatment Rena.pdf",width=14, height=10)

# gather data
discData <- data.table( erbphenodat%>%dplyr::select(Cell.ID,Patient.Study.ID,dynamic_class3,Day,ARM.x,Treatment,Celltype,PCA2,
                                                    responseLab) )
setnames( discData , old=c("PCA2","ARM.x"), new=c("ERBBaxis","ARM"))

#write.csv(discData,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/DiscoveryCohortERBBPrePost.csv")


# Perform statistical analyses
#summary( lm(ERBBaxis~ -1+Patient.Study.ID+ Day : Patient.Study.ID ,data=fread(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/Figure2/DiscoveryCohortERBBPrePost.csv")))
#summary( lmer(ERBBaxis~ Day +(1+Day|Patient.Study.ID) ,data=fread(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/Figure2/DiscoveryCohortERBBPrePost.csv")))
summary( lmer(ERBBaxis~ Day +(1+Day|Patient.Study.ID) ,data=discData))
Mod1 <- lmer(ssGSEA ~ -1+dynamic_class3*day_fact*ARM + (1+day_fact|Patient.Study.ID), REML=FALSE,data= u_dat) 
summary(Mod1)
summary( lm(ERBBaxis~ -1+Patient.Study.ID+ Day : Patient.Study.ID,data=discData))
#summary( lmer(ERBBaxis~ (Day+dynamic_class3) +(1+Day|Patient.Study.ID) ,data=fread(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/Figure2/ValidationCohortERBBPrePost.csv")))


# project validation samples into same PCA space
Validation_inputdata$PCA2 <- predict( res.pca,newdata=data.frame(Validation_inputdata %>% dplyr::select(ERBBsetlist  ) ))[,2]

# extract patient data with paired pre-/post samples
erbphenodat <- Validation_inputdata[nDay0>20&nDay180>20][Day%in%c(0,180)]

# order for plotting as in discovery cohort 
ERBB_mu_saveC2 <- data.table(erbphenodat%>%group_by(Patient.Study.ID,dynamic_class3,ARM,Day)%>%
                               dplyr::summarise(median_ERB_high= median(PCA2),
                                         median_NAGASHERB_high= median(NAGASHIMA_EGF_SIGNALING_UP),
                                         N=n() ))
ERBB_mu_saveC2[, initERB:= sum((Day==0)*median_ERB_high), by=Patient.Study.ID ]
ERBB_mu_saveC2[, change:=median_ERB_high -initERB]
erbphenodat$Patient.Study.ID <- factor(erbphenodat$Patient.Study.ID , levels= ERBB_mu_saveC2[Day==180][order(-change)]$Patient.Study.ID)
ERBB_mu_saveC2$Patient.Study.ID <- factor(ERBB_mu_saveC2$Patient.Study.ID , levels= ERBB_mu_saveC2[Day==180][order(-change)]$Patient.Study.ID)
ERBB_mu_saveC2[order(Patient.Study.ID,-change)]

ERBB_mu_saveC2 <- data.table(erbphenodat%>%group_by(Patient.Study.ID,dynamic_class3,ARM,Day)%>%
                               summarise(median_ERB_high= median(PCA2),
                                         median_NAGASHERB_high= median(NAGASHIMA_EGF_SIGNALING_UP),
                                         N=n() ))

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
#ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/ERBB phenotype shift pre to post treatment1.pdf",width=14, height=10)
ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/ERBB phenotype shift pre and post treatment Validation Rena.pdf",width=14, height=10)


