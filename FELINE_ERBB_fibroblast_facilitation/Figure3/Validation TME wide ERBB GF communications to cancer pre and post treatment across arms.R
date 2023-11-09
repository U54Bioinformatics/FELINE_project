rm(list=ls())
library(glmnet);require(dplyr);require(data.table);require(tidyr); require(ggplot2); require(lmerTest)
load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/cohort2 TotalCommunicationAllArms/cohort2 TotalCommunicationAllArms.RData")

GFFibroCancerComm2 <-SumTotComm0[PhenoCelltype=="Cancer cells"][Day%in%c(0,180)][
  Pair.Name%in% c("ADAM17_ERBB4","BTC_ERBB4" , "EGF_ERBB4" , "NRG1_ERBB4", "NRG2_ERBB4", "NRG3_ERBB4", "TGFA_ERBB4", "HBEGF_ERBB4", 
                 "ANXA1_EGFR", "AREG_EGFR", "BTC_EGFR" , "CDH1_EGFR","GNAI2_EGFR", "ICAM1_EGFR", "TGFA_EGFR", "HBEGF_EGFR", "EGF_EGFR", 
                 "BTC_ERBB2", "EFNB1_ERBB2", "HLA-A_ERBB2", "NRG1_ERBB2", "SEMA4D_ERBB2", "TGFA_ERBB2", "HBEGF_ERBB2", "EGF_ERBB2",   
                 "AREG_ERBB3", "BTC_ERBB3", "EGF_ERBB3", "NRG1_ERBB3", "NRG2_ERBB3"  ) 
] [!grepl("BTC",Pair.Name)]

GFFibroCancerComm2[ ,scalelnTransductionFC:= scale(log(1+Transduction)),by=c("Pair.Name")]
# normalise coloration relative to baseline
GFFibroCancerComm2[ ,d0PatMean:=sum( (Day==0)*scalelnTransductionFC )/sum(Day==0)   ,by=Patient.Study.ID ]
GFFibroCancerComm2[!is.finite(d0PatMean), d0PatMean:=NA]
GFFibroCancerComm2[,d0PatMean:=mean(d0PatMean, na.rm=T)]
GFFibroCancerComm2[, TimePoint:= "Pre treatment"]
GFFibroCancerComm2[Day==180, TimePoint:= "Post treatment"]
GFFibroCancerComm2$TimePoint<- factor(GFFibroCancerComm2$TimePoint ,levels=c("Pre treatment","Post treatment"))
ggplot(GFFibroCancerComm2,
       aes(y=Pair.Name, x=Patient.Study.ID ,fill= scalelnTransductionFC-d0PatMean))+geom_tile()+theme_classic()+
  scale_fill_viridis_c(name="Fibroblast to Cancer \n Growth factor \n communication",option="D") +
  facet_wrap(~TimePoint,scales="free_x", ncol=2) +  theme(aspect.ratio=2 ,axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ERBBsubset<- GFFibroCancerComm2[grepl("_ERBB",Pair.Name)|grepl("_EGFR",Pair.Name)]
GRPPLOT <- GFFibroCancerComm2[!grepl("BTC",Pair.Name)]
ERBBsubset$Pair.Name <- gsub("_","-" ,ERBBsubset$Pair.Name)
ERBBsubset$Pair.Name<- factor(ERBBsubset$Pair.Name ,levels=
                                c("ADAM17-ERBB4","BTC-ERBB4" , "EGF-ERBB4" , "NRG1-ERBB4", "NRG2-ERBB4", "NRG3-ERBB4", "TGFA-ERBB4", "HBEGF-ERBB4", 
                                  "ANXA1-EGFR", "AREG-EGFR", "BTC-EGFR" , "CDH1-EGFR","GNAI2-EGFR", "ICAM1-EGFR", "TGFA-EGFR", "HBEGF-EGFR", "EGF-EGFR", 
                                  "BTC-ERBB2", "EFNB1-ERBB2", "HLA-A-ERBB2", "NRG1-ERBB2", "SEMA4D-ERBB2", "TGFA-ERBB2", "HBEGF-ERBB2", "EGF-ERBB2",   
                                  "AREG-ERBB3", "BTC-ERBB3", "EGF-ERBB3", "NRG1-ERBB3", "NRG2-ERBB3"  ) )
ggplot(ERBBsubset[!grepl("BTC",Pair.Name)], 
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
  scale_fill_manual(name="Timepoint",values=rev(pal_npg()(2))) +
  labs(x="Cancer growth factor receptor", y="Growth factor signalling \n to cancer cells")+ theme(aspect.ratio=1,legend.position="bottom")
#ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Validation cohort2 Growth factor signalling to cancer cells increased relative to baseline ERBB only.pdf", height=10, width = 12)
ggplot(GRPPLOT2, 
       aes(x=Receptor ,y= y , fill=TimePoint,group=interaction(Receptor,TimePoint) ))+
  geom_boxplot(outlier.color=NA)+
  geom_point(position=position_dodge(width=0.9))+
  theme_classic(base_size=26)+
  scale_fill_manual(name="Timepoint",values= c("yellow","red")) +
  labs(x="Cancer growth factor receptor", y="Growth factor signalling \n to cancer cells")+ theme(aspect.ratio=1,legend.position="bottom")
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

ggplot(GRPPLOT2, 
       aes(x=Receptor ,y= y , fill=TimePoint,group=interaction(Receptor,TimePoint) ))+
  geom_boxplot(outlier.color=NA)+
  geom_point(position=position_dodge(width=0.9))+
  theme_classic(base_size=32)+
  scale_fill_manual(name="Timepoint",values=rev(pal_npg()(2))) +
  labs(x="Cancer growth factor receptor", y="Growth factor signalling \n to cancer cells")+ theme(aspect.ratio=0.81,legend.position="bottom")
#ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Validation cohort2 Growth factor signalling to cancer cells increased relative to baseline ERBB onlyB.pdf", height=10, width = 12)
validStatsData <- GRPPLOT%>%dplyr::select( Patient.Study.ID, ARM,Receptor,Pair.Name ,scalelnTransductionFC,d0PatMean,TimePoint)
validStatsData$TimePoint <- factor(validStatsData$TimePoint, levels= c("Pre treatment","Post treatment"))
m0 <- lmer( I(scalelnTransductionFC-d0PatMean ) ~ 0+Receptor + TimePoint:Receptor + (1|Pair.Name:Patient.Study.ID)+ (1|Patient.Study.ID) , data= validStatsData[] )
summary(m0)
#write.csv(validStatsData,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/ValidationCohort_Stats_TMEwideERBBcommunicationwithtoCancerPrePost.csv")




m1 <- lmer( I(scalelnTransductionFC-d0PatMean ) ~TimePoint+ (1|Pair.Name:Patient.Study.ID) , data=GRPPLOT[Receptor=="EGFR"] )
m2 <- lmer( I(scalelnTransductionFC-d0PatMean ) ~TimePoint+ (1|Pair.Name:Patient.Study.ID) , data=GRPPLOT[Receptor=="ERBB2"] )
m3 <- lmer( I(scalelnTransductionFC-d0PatMean ) ~TimePoint+ (1|Pair.Name:Patient.Study.ID) , data=GRPPLOT[Receptor=="ERBB3"] )
m4 <- lmer( I(scalelnTransductionFC-d0PatMean ) ~TimePoint+ (1|Pair.Name:Patient.Study.ID) , data=GRPPLOT[Receptor=="ERBB4"] )
summary(m1)
summary(m2)
summary(m3)
summary(m4)

summary( lm(y~0+TimePoint ,data=GRPPLOT2[Receptor=="EGFR"]))
summary( lm(y~0+TimePoint ,data=GRPPLOT2[Receptor=="ERBB2"]))
summary( lm(y~0+TimePoint ,data=GRPPLOT2[Receptor=="ERBB3"]))
summary( lm(y~0+TimePoint ,data=GRPPLOT2[Receptor=="ERBB4"]))

