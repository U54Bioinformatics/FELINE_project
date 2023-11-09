rm(list=ls())
library(glmnet);require(dplyr);require(data.table);require(tidyr); require(ggplot2); require(lmerTest)
load( "~/Dropbox/Cancer_pheno_evo/data/FELINE2/LigandReceptor/Filtered_Human-2015-Ramilowski-LR-pairs.RData" )
load(file ="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/cohort2 Communication output merged/PopulationCommunicationMerged ALL ARMS cohort2.RData" )#CCI,allphenotypes, uu,perIndiv,
LRgenelist <- unique( data.table(merge(CCI[],LRpairsFiltered,by="Pair.Name"))$HPMR.Receptor) 
rm(list= c("LRpairs", "LRpairsFiltered"))

CCI[,Treat:="CombinationRibo"]
CCI[ARM=="A",Treat:="LetrozoleAlone"]
CCI[,scaleTransduction:=scale(TransductionMu,center=F),by=c("Pair.Name")] #"ReceptorPhenoCelltype",
CCI[!is.finite(TransductionMu),scaleTransduction:=0]
CCI[,scalelnTransductionMu:=scale(log(1+TransductionMu)),by=Pair.Name]
CCI[,scaleTransduction2:=exp(scale(log(TransductionMu),center=F)),by=c("Pair.Name")] 

load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/GeneOntology/growthFactorReceptors.RData" )

wchComms<- unique( c( "EGFR","ERBB2","ERBB3","ERBB4"))
isApath<-sapply( 1:length(wchComms),function(ii){grepl(paste0("_",wchComms[ii]),unique(CCI$Pair.Name) ) } ) 
l1 <- unique(CCI$Pair.Name)[rowSums(isApath)>0 ]

# select start and end timepoint data for fibroblast cancer interactions
GFFibroCancerComm <- CCI[Pair.Name%in%l1][
    Day!=14
    ][ReceptorPhenoCelltype=="Cancer cells"][LigandPhenoCelltype=="Fibroblasts"]#|LigandPhenoCelltype=="Cancer cells"]
GFFibroCancerComm[,npatSamp:=length(unique(Patient.Study.ID)),by=Pair.Name]
GFFibroCancerComm[,propPatSamp:=npatSamp/length(unique(GFFibroCancerComm$Patient.Study.ID))]
GFFibroCancerComm[,npairSamp:=length(unique(Pair.Name)),by=Patient.Study.ID]
GFFibroCancerComm[,propPairSamp:=npairSamp/length(unique(GFFibroCancerComm$Pair.Name))]

# filter out rarely measured communication pathways
GFFibroCancerComm2 <- GFFibroCancerComm[propPatSamp>0.5][propPairSamp>0.5]
GFFibroCancerComm2[ ,scalelnTransductionFC:= scale(log(TransductionMu)),by=c("Pair.Name","LigandPhenoCelltype")]
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

#Select ERBB family memeber communications
ERBBsubset<- GFFibroCancerComm2[grepl("_ERBB",Pair.Name)|grepl("_EGFR",Pair.Name)]
GRPPLOT <- ERBBsubset[!grepl("BTC",Pair.Name)]
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
ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Fibroblast to Cancer Growth Factor communication heatmap increase relative to baseline.png", height=10, width = 10)

m1 <- lmer( I(scalelnTransductionFC-d0PatMean ) ~TimePoint+ (1|Patient.Study.ID) , data=ERBBsubset[!grepl("BTC",Pair.Name)] )
summary(m1)

ggplot(ERBBsubset[!grepl("BTC",Pair.Name)], 
       aes(y=Pair.Name, x=Patient.Study.ID ,fill= scalelnTransductionFC-d0PatMean))+geom_tile()+theme_classic(base_size=15)+
  scale_fill_viridis_c(name="Fibroblast to Cancer \n Growth factor \n communication \n (scaled)",option="D") +
  labs(y="Communication pathway", x="Patient tumor")+
  facet_wrap(~TimePoint,scales="free_x", ncol=2) +  theme(aspect.ratio=1.5 ,axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Fibroblast to Cancer Growth Factor communication heatmap increase relative to baseline.png", height=10, width = 10)


GRPPLOT[,c("Ligand","Receptor"):= tstrsplit(Pair.Name,"_",fixed=T)]

GRPPLOT1<- data.table(GRPPLOT %>%group_by(Patient.Study.ID,Treat,ARM,Day,TimePoint,Receptor)%>% dplyr::summarise(y=mean(scalelnTransductionFC-d0PatMean , na.rm=T)))
GRPPLOT1[,grpx:=paste(Patient.Study.ID,TimePoint,sep="AT")]
GRPPLOT2<- data.table( GRPPLOT1 %>% select(-c(Patient.Study.ID,TimePoint,Treat,ARM,Day)) %>% spread(grpx,y,fill=0) %>%gather(grpx,y,-Receptor))
GRPPLOT2[,c("Patient.Study.ID","TimePoint"):= tstrsplit(grpx,"AT",fixed=T)]
GRPPLOT2$TimePoint <- factor(GRPPLOT2$TimePoint , levels= rev( unique(GRPPLOT2$TimePoint) ) )

ggplot(GRPPLOT2, 
       aes(x=Receptor ,y= y , fill=TimePoint,group=interaction(Receptor,TimePoint) ))+
  geom_boxplot(outlier.color=NA)+
  geom_point(position=position_dodge(width=0.9))+
  theme_classic(base_size=26)+
  scale_fill_manual(name="Timepoint",values=rev(pal_npg()(2))) +
  labs(x="Cancer growth factor receptor", y="Fibroblast-cancer cell communication \n (relative to baseline average)")+ theme(aspect.ratio=1,legend.position="bottom")

ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Validation Fibroblast to Cancer Growth ERBB Factor communication receiver increase relative to baseline.png", height=10, width = 12)
ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Validation Fibroblast to Cancer Growth ERBB Factor communication receiver increase relative to baseline.pdf", height=10, width = 12)

ggplot(GRPPLOT2, 
       aes(x=Receptor ,y= y , fill=TimePoint,group=interaction(Receptor,TimePoint) ))+
  geom_boxplot(outlier.color=NA)+
  geom_point(position=position_dodge(width=0.9))+
  theme_classic(base_size=26)+
  scale_fill_manual(name="Timepoint",values=c("yellow","red")) +
  labs(x="Cancer growth factor receptor", y="Fibroblast-cancer cell communication \n (relative to baseline average)")+ theme(aspect.ratio=1,legend.position="bottom")
#ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Validation Fibroblast to Cancer Growth ERBB Factor communication receiver increase relative to baselineYRsquare.pdf", height=8, width = 8,dpi=320)


validData <- data.table( merge(GRPPLOT2 ,unique(GRPPLOT1%>%dplyr::select(Patient.Study.ID,Treat,ARM)),by="Patient.Study.ID", all.x=T) %>%
                           dplyr::select(c(Receptor,Patient.Study.ID ,ARM,Treat,y,TimePoint ) ))
setnames( validData , old=c("Receptor","y","Treat"), new=c("CommunicationReceptor","logFibroblasttoCancerCommunicationScore","treat"))
validData$TimePoint <- factor(validData$TimePoint, levels=c("Pre treatment","Post treatment"))
#write.csv(validData,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/ValidationCohortFibroblastERBBcommunicationwithtoCancerPrePost.csv")

validStatsData <- GRPPLOT%>%dplyr::select( Patient.Study.ID, ARM,Receptor,Pair.Name ,scalelnTransductionFC,d0PatMean,TimePoint)
validStatsData$TimePoint <- factor(validStatsData$TimePoint, levels= c("Pre treatment","Post treatment"))

m0 <- lmer( I(scalelnTransductionFC-d0PatMean ) ~ 0+Receptor + TimePoint:Receptor + (1|Pair.Name:Patient.Study.ID)+ (1|Patient.Study.ID) , data= validStatsData[] )
summary(m0)
#write.csv(validStatsData,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/ValidationCohort_Stats_FibroblastERBBcommunicationwithtoCancerPrePost.csv")


