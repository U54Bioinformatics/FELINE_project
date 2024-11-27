rm(list=ls())
require(dplyr);require(data.table);require(tidyr); require(ggplot2); require(lmerTest);require(ggsci)


# Define data location
SourceDataLoc<-"/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript Feline ERBB facilitation/Molecular Systems Biology/Revision and submission folder/Source Data/"
Intermediateloc<-paste0(SourceDataLoc,"Figure 2/")

# Load discovery communication data
#load(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/TotalCommunicationAllArms/TotalCommunicationAllArms.RData")
#DiscoveryTumotCancerComms <- SumTotComm0[ PhenoCelltype=="Cancer cells"],
#write.csv(SumTotComm0[ PhenoCelltype=="Cancer cells"], file= paste0(Intermediateloc,"SourceData_Figure_2D_DiscoveryTumortoCancerComms.csv"))
DiscoveryTumotCancerComms <- data.table(read.csv( file= paste0(Intermediateloc,"SourceData_Figure_2D_DiscoveryTumortoCancerComms.csv")))

# define erbb Ligand-receptro pairs :
LR_ERBB <- c("ADAM17_ERBB4", "EGF_ERBB4" , "NRG1_ERBB4", "NRG2_ERBB4", "NRG3_ERBB4", "TGFA_ERBB4", "HBEGF_ERBB4", 
  "ANXA1_EGFR", "AREG_EGFR",  "CDH1_EGFR","GNAI2_EGFR", "ICAM1_EGFR", "TGFA_EGFR", "HBEGF_EGFR", "EGF_EGFR", 
  "EFNB1_ERBB2", "HLA-A_ERBB2", "NRG1_ERBB2", "SEMA4D_ERBB2", "TGFA_ERBB2", "HBEGF_ERBB2", "EGF_ERBB2",   
  "AREG_ERBB3",  "EGF_ERBB3", "NRG1_ERBB3", "NRG2_ERBB3"  ) 

# Extract ERBB GF communications to cancer and compare pre-/post-treatment signaling
GFTumorCancerComm2 <- DiscoveryTumotCancerComms[Day%in%c(0,180)][Pair.Name%in% LR_ERBB] 

# Normalize for patient level variation
GFTumorCancerComm2[ ,scalelnTransductionFC:= scale(log(1+Transduction)),by=c("Pair.Name")]
# normalise communication relative to baseline
GFTumorCancerComm2[ ,d0PatMean:=sum( (Day==0)*scalelnTransductionFC )/sum(Day==0)   ,by=Patient.Study.ID ]
GFTumorCancerComm2[!is.finite(d0PatMean), d0PatMean:=NA]
GFTumorCancerComm2[,d0PatMean:=mean(d0PatMean, na.rm=T)]
GFTumorCancerComm2[, TimePoint:= "Pre treatment"]
GFTumorCancerComm2[Day==180, TimePoint:= "Post treatment"]
GFTumorCancerComm2$TimePoint<- factor(GFTumorCancerComm2$TimePoint ,levels=c("Pre treatment","Post treatment"))

# duplicate dataset to annotate for plotting
ERBBsubset<- GFTumorCancerComm2[grepl("_ERBB",Pair.Name)|grepl("_EGFR",Pair.Name)]
GRPPLOT <- GFTumorCancerComm2
ERBBsubset$Pair.Name <- gsub("_","-" ,ERBBsubset$Pair.Name)
ERBBsubset$Pair.Name<- factor(ERBBsubset$Pair.Name ,levels=
                                c("ADAM17-ERBB4", "EGF-ERBB4" , "NRG1-ERBB4", "NRG2-ERBB4", "NRG3-ERBB4", "TGFA-ERBB4", "HBEGF-ERBB4", 
                                  "ANXA1-EGFR", "AREG-EGFR", "CDH1-EGFR","GNAI2-EGFR", "ICAM1-EGFR", "TGFA-EGFR", "HBEGF-EGFR", "EGF-EGFR", 
                                   "EFNB1-ERBB2", "HLA-A-ERBB2", "NRG1-ERBB2", "SEMA4D-ERBB2", "TGFA-ERBB2", "HBEGF-ERBB2", "EGF-ERBB2",   
                                  "AREG-ERBB3",  "EGF-ERBB3", "NRG1-ERBB3", "NRG2-ERBB3"  ) )
ggplot(ERBBsubset, 
       aes(y=Pair.Name, x=Patient.Study.ID ,fill= scalelnTransductionFC-d0PatMean))+geom_tile()+theme_classic(base_size=15)+
  scale_fill_viridis_c(name="Fibroblast to Cancer \n Growth factor \n communication \n (scaled)",option="D") +
  labs(y="Communication pathway", x="Patient tumor")+
  facet_wrap(~TimePoint,scales="free_x", ncol=2) +  theme(aspect.ratio=1.5 ,axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Aggregate for boxplotting
GRPPLOT[,c("Ligand","Receptor"):= tstrsplit(Pair.Name,"_",fixed=T)]
GRPPLOT1 <- data.table(GRPPLOT %>%group_by(Patient.Study.ID,ARM,TimePoint,Receptor)%>% dplyr::summarise(
  y0=mean(Transduction , na.rm=T),
  y=mean(scalelnTransductionFC-d0PatMean , na.rm=T)))
GRPPLOT1[,grpx:=paste(Patient.Study.ID,TimePoint,sep="AT")]
GRPPLOT2<- data.table( GRPPLOT1 %>% select(-c(Patient.Study.ID,TimePoint,y0,ARM)) %>% spread(grpx,y,fill=0) %>%gather(grpx,y,-Receptor))
GRPPLOT2[,c("Patient.Study.ID","TimePoint"):= tstrsplit(grpx,"AT",fixed=T)]
GRPPLOT2$TimePoint <- factor(GRPPLOT2$TimePoint , levels= rev( unique(GRPPLOT2$TimePoint) ) )

ggplot(GRPPLOT2, 
       aes(x=Receptor ,y= y , fill=TimePoint,group=interaction(Receptor,TimePoint) ))+
  geom_boxplot(outlier.color=NA)+
  geom_point(position=position_dodge(width=0.9))+
  theme_classic(base_size=26)+
  scale_fill_manual(name="Timepoint",values= c("yellow","red")) +
  labs(x="Cancer growth factor receptor", y="Growth factor signalling \n to cancer cells")+ theme(aspect.ratio=1,legend.position="bottom")

# collect output data and perform statistical analysis (lmer)
discData <- data.table( GRPPLOT1 %>% 
                          dplyr::select(c(Receptor,Patient.Study.ID ,ARM,y0,y,TimePoint ) ))
discData[,treat:="Letrozole alone"]
discData[ARM!="A",treat:="Combination ribociclib"]
setnames( discData , old=c("Receptor","y0","y"), new=c("CommunicationReceptor","CommunicationScore","ScaledlnCommunicationScore"))
discData<- discData%>%dplyr::select(CommunicationReceptor, Patient.Study.ID, ARM,treat, CommunicationScore, ScaledlnCommunicationScore,TimePoint)
#write.csv(discData,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/DiscoveryCohortTMEwideERBBcommunicationwithtoCancerPrePostNEW.csv")

discStatsData <- GRPPLOT%>%dplyr::select( Patient.Study.ID, ARM,Receptor,Pair.Name ,scalelnTransductionFC,d0PatMean,TimePoint)
discStatsData$TimePoint <- factor(discStatsData$TimePoint, levels= c("Pre treatment","Post treatment"))
m0d <- lmer( I(scalelnTransductionFC-d0PatMean ) ~ 0+ Receptor+ TimePoint:Receptor + (1|Pair.Name:Patient.Study.ID)+ (1|Patient.Study.ID) , data= discStatsData[] )
summary(m0)
#write.csv(discStatsData,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/DiscoveryCohort_Stats_TMEwideERBBcommunicationwithtoCancerPrePost.csv")



# Load discovery communication data
#load(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/TotalCommunicationAllArms/TotalCommunicationAllArms.RData")



### Replicate analysis for validartion cohort
#rm(list=ls())
load(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/cohort2 TotalCommunicationAllArms/cohort2 TotalCommunicationAllArms.RData")

#ValidationTumorCancerComms <- SumTotComm0[ PhenoCelltype=="Cancer cells"]
#write.csv(ValidationTumotCancerComms, file= paste0(Intermediateloc,"SourceData_Figure_2D_ValidationTumortoCancerComms.csv"))
ValidationTumorCancerComms <- data.table(read.csv( file= paste0(Intermediateloc,"SourceData_Figure_2D_ValidationTumortoCancerComms.csv")))



GFFibroCancerComm2 <-ValidationTumorCancerComms[Day%in%c(0,180)][
  Pair.Name%in% c("ADAM17_ERBB4", "EGF_ERBB4" , "NRG1_ERBB4", "NRG2_ERBB4", "NRG3_ERBB4", "TGFA_ERBB4", "HBEGF_ERBB4", 
                  "ANXA1_EGFR", "AREG_EGFR",  "CDH1_EGFR","GNAI2_EGFR", "ICAM1_EGFR", "TGFA_EGFR", "HBEGF_EGFR", "EGF_EGFR", 
                   "EFNB1_ERBB2", "HLA-A_ERBB2", "NRG1_ERBB2", "SEMA4D_ERBB2", "TGFA_ERBB2", "HBEGF_ERBB2", "EGF_ERBB2",   
                  "AREG_ERBB3",  "EGF_ERBB3", "NRG1_ERBB3", "NRG2_ERBB3"  ) 
] 

GFFibroCancerComm2[ ,scalelnTransductionFC:= scale(log(1+Transduction)),by=c("Pair.Name")]
# normalise coloration relative to baseline
GFFibroCancerComm2[ ,d0PatMean:=sum( (Day==0)*scalelnTransductionFC )/sum(Day==0)   ,by=Patient.Study.ID ]
GFFibroCancerComm2[!is.finite(d0PatMean), d0PatMean:=NA]
GFFibroCancerComm2[,d0PatMean:=mean(d0PatMean, na.rm=T)]
GFFibroCancerComm2[, TimePoint:= "Pre treatment"]
GFFibroCancerComm2[Day==180, TimePoint:= "Post treatment"]
GFFibroCancerComm2$TimePoint<- factor(GFFibroCancerComm2$TimePoint ,levels=c("Pre treatment","Post treatment"))

ERBBsubset<- GFFibroCancerComm2[grepl("_ERBB",Pair.Name)|grepl("_EGFR",Pair.Name)]
GRPPLOT <- GFFibroCancerComm2
ERBBsubset$Pair.Name <- gsub("_","-" ,ERBBsubset$Pair.Name)
ERBBsubset$Pair.Name<- factor(ERBBsubset$Pair.Name ,levels=
                                c("ADAM17-ERBB4","EGF-ERBB4" , "NRG1-ERBB4", "NRG2-ERBB4", "NRG3-ERBB4", "TGFA-ERBB4", "HBEGF-ERBB4", 
                                  "ANXA1-EGFR", "AREG-EGFR", "CDH1-EGFR","GNAI2-EGFR", "ICAM1-EGFR", "TGFA-EGFR", "HBEGF-EGFR", "EGF-EGFR", 
                                   "EFNB1-ERBB2", "HLA-A-ERBB2", "NRG1-ERBB2", "SEMA4D-ERBB2", "TGFA-ERBB2", "HBEGF-ERBB2", "EGF-ERBB2",   
                                  "AREG-ERBB3", "EGF-ERBB3", "NRG1-ERBB3", "NRG2-ERBB3"  ) )
ggplot(ERBBsubset, 
       aes(y=Pair.Name, x=Patient.Study.ID ,fill= scalelnTransductionFC-d0PatMean))+geom_tile()+theme_classic(base_size=15)+
  scale_fill_viridis_c(name="Fibroblast to Cancer \n Growth factor \n communication \n (scaled)",option="D") +
  labs(y="Communication pathway", x="Patient tumor")+
  facet_wrap(~TimePoint,scales="free_x", ncol=2) +  theme(aspect.ratio=1.5 ,axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

GRPPLOT[,c("Ligand","Receptor"):= tstrsplit(Pair.Name,"_",fixed=T)]

GRPPLOT1 <- data.table(GRPPLOT %>%group_by(Patient.Study.ID,ARM,TimePoint,Receptor)%>% dplyr::summarise(
  y0=mean(Transduction , na.rm=T),
  y=mean(scalelnTransductionFC-d0PatMean , na.rm=T)))
GRPPLOT1[,grpx:=paste(Patient.Study.ID,TimePoint,sep="AT")]
GRPPLOT2<- data.table( GRPPLOT1 %>% select(-c(Patient.Study.ID,TimePoint,y0,ARM)) %>% spread(grpx,y,fill=0) %>%gather(grpx,y,-Receptor))
GRPPLOT2[,c("Patient.Study.ID","TimePoint"):= tstrsplit(grpx,"AT",fixed=T)]
GRPPLOT2$TimePoint <- factor(GRPPLOT2$TimePoint , levels= rev( unique(GRPPLOT2$TimePoint) ) )

ggplot(GRPPLOT2, 
       aes(x=Receptor ,y= y , fill=TimePoint,group=interaction(Receptor,TimePoint) ))+
  geom_boxplot(outlier.color=NA)+
  geom_point(position=position_dodge(width=0.9))+
  theme_classic(base_size=26)+
  scale_fill_manual(name="Timepoint",values= c("yellow","red")) +
  labs(x="Cancer growth factor receptor", y="Growth factor signalling \n to cancer cells")+ theme(aspect.ratio=1,legend.position="bottom")

#ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Validation cohort2 Growth factor signalling to cancer cells increased relative to baseline ERBB onlyYR.pdf", height=10, width = 12)
#ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Validation cohort2 Growth factor signalling to cancer cells increased relative to baseline ERBB onlyYRsquare.pdf", height=8, width = 8,dpi=320)


validData <- data.table( GRPPLOT1   %>%
                           dplyr::select(c(Receptor,Patient.Study.ID ,ARM,y0,y,TimePoint ) ))
validData[,treat:="Letrozole alone"]
validData[ARM!="A",treat:="Combination ribociclib"]

setnames( validData , old=c("Receptor","y0","y"), new=c("CommunicationReceptor","CommunicationScore","ScaledlnCommunicationScore"))
validData<- validData%>%dplyr::select(CommunicationReceptor, Patient.Study.ID, ARM,treat, CommunicationScore, ScaledlnCommunicationScore,TimePoint)
#write.csv(validData,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/ValidationCohortTMEwideERBBcommunicationwithtoCancerPrePost.csv")
#read.csv(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/DiscoveryCohortTMEwideERBBcommunicationwithtoCancerPrePost.csv")[1:3,]
#read.csv(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/ValidationCohortTMEwideERBBcommunicationwithtoCancerPrePost.csv")[1:3,]


#ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Validation cohort2 Growth factor signalling to cancer cells increased relative to baseline ERBB onlyB.png", height=10, width = 12)
#ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Validation cohort2 Growth factor signalling to cancer cells increased relative to baseline ERBB onlyB.pdf", height=10, width = 12)
validStatsData <- GRPPLOT%>%dplyr::select( Patient.Study.ID, ARM,Receptor,Pair.Name ,scalelnTransductionFC,d0PatMean,TimePoint)
validStatsData$TimePoint <- factor(validStatsData$TimePoint, levels= c("Pre treatment","Post treatment"))
#m0 <- lmer( I(scalelnTransductionFC-d0PatMean ) ~ Receptor + TimePoint:Receptor + (1+TimePoint|Patient.Study.ID/Pair.Name) , data= validStatsData[] )
m0 <- lmer( I(scalelnTransductionFC-d0PatMean ) ~ 0+Receptor + TimePoint:Receptor + (1|Pair.Name:Patient.Study.ID)+ (1|Patient.Study.ID) , data= validStatsData[] )
summary(m0)
#write.csv(validStatsData,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/ValidationCohort_Stats_TMEwideERBBcommunicationwithtoCancerPrePost.csv")

