rm(list=ls())
load( file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/temporary/whocommValidation.RData")
Valid <- data.table(Cohort="Validation", plt1)

load( file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/temporary/whocommDiscovery.RData")
Discov <- data.table(Cohort="Discovery", plt1)

plt1<- rbind(Discov%>%dplyr::select(-c(Ligandiscancer, Receptoriscancer)),Valid)

paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"


ggplot(plt1, aes(y= log(exp(meanSig)-1), x= LigandPhenoCelltype, col= LigandPhenoCelltype,fill=LigandPhenoCelltype,group=interaction(LigandPhenoCelltype)))+
  geom_boxplot(alpha=0.5,outlier.color =NA)+ 
  stat_boxplot(geom="errorbar")+
  geom_jitter(aes(shape=Cohort),size=1.5,height=0,width=0.2)+
  theme_classic(base_size=26)+theme(axis.text.x=element_text(angle=90,vjust=0.5))+
  theme(aspect.ratio=1,legend.position="none") +labs(y="Contribution to \n tumor wide communication \n (mean)",x= "Cell type") +
  facet_wrap(~Treatmentlab,nrow=1)+ 
  scale_y_continuous(breaks=log(c(0.001,0.01,0.1,1)),labels=c(0.001,0.01,0.1,1))
#ggsave(paste0(paperfile,"Discovery and Validation Ribo and Letrozole Cell type strength of communicaiton received.png"),height=10,width=10, dpi=320)

ggplot(plt1, aes(y= log(exp(meanSig)-1), x= ReceptorPhenoCelltype, col= ReceptorPhenoCelltype,fill=ReceptorPhenoCelltype,group=interaction(ReceptorPhenoCelltype)))+
  geom_boxplot(alpha=0.5,outlier.color =NA)+
  stat_boxplot(geom="errorbar")+
  geom_jitter(aes(shape=Cohort),size=1.5,height=0,width=0.2)+
  theme_classic(base_size=26)+
  theme(axis.text.x=element_text(angle=90,vjust=0.5)) +
  theme(aspect.ratio=1,legend.position="none") +labs(y="Strength of communicaiton \n received (mean)",x= "Cell type") +
  facet_wrap(~Treatmentlab)+ 
  scale_y_continuous(breaks=log(c(0.001,0.01,0.1,1)),labels=c(0.001,0.01,0.1,1))
#ggsave(paste0(paperfile,"Discovery and Validation  Ribo and Letrozole Cell type strength of communicaiton received.png"),height=10,width=10, dpi=320)



plt1[, Ligandiscancer:=F]
plt1[LigandPhenoCelltype=="Cancer cells", Ligandiscancer:=T]
summary(lm( log(exp(meanSig)-1)~  Ligandiscancer*Cohort , plt1[Treatmentlab=="Combination ribociclib"]))
summary(lm( log(exp(meanSig)-1)~  Ligandiscancer*Cohort , plt1[Treatmentlab=="Letrozole alone"]))

summary(lm( log(exp(meanSig)-1)~  Ligandiscancer , plt1[Treatmentlab=="Combination ribociclib"]))
summary(lm( log(exp(meanSig)-1)~  Ligandiscancer , plt1[Treatmentlab=="Letrozole alone"]))

plt1[, Receptoriscancer:=F]
plt1[ReceptorPhenoCelltype=="Cancer cells", Receptoriscancer:=T]
summary(lm( log(exp(meanSig)-1)~ Receptoriscancer*Cohort , plt1[Treatmentlab=="Combination ribociclib"]))
summary(lm( log(exp(meanSig)-1)~ Receptoriscancer*Cohort , plt1[Treatmentlab=="Letrozole alone"]))


summary(lm( log(exp(meanSig)-1)~  Receptoriscancer , plt1[Treatmentlab=="Combination ribociclib"]))
summary(lm( log(exp(meanSig)-1)~  Receptoriscancer , plt1[Treatmentlab=="Letrozole alone"]))


#save(plt1, file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/Figure2/Discovery and Validation cell type signal contribution.RData")
