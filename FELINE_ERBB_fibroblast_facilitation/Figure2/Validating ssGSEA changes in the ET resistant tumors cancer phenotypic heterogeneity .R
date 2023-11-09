rm(list=ls())
require(mgcv);require(data.table);require(dplyr);require(ggplot2);require(tidyr);require(lme4);require(lmerTest);require(parallel)
library(effects);require(umap);require(Rfast);require(ider);library("dendextend");library(ggdendro);require(ggsci);require(viridis);require("Rdimtools")
require(merTools)

load(file=paste0("~/Dropbox/Cancer_pheno_evo/data/FELINE2/PhenotypesAllArmsCohort2/UpdatedRevisednewCancer_cellsCancer cells.RData"))

u_dat[,DayLab:="Day 0"]
u_dat[Day==14, DayLab:="Day 14"]
u_dat[Day==180, DayLab:="Day 180"]
u_dat$ssGSEA <- unname(unlist(u_dat[,"NAGASHIMA_EGF_SIGNALING_UP"]))
u_dat[, mu_ssGSEA_t0:= sum(ssGSEA *(Day==0))/sum(Day==0) , by=Patient.Study.ID ]
u_dat[is.na(mu_ssGSEA_t0),mu_ssGSEA_t0:=mean(u_dat[Day==0]$ssGSEA)]
u_dat[,standardised_ssGSEA:=ssGSEA-mu_ssGSEA_t0]

load( file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Cohort2Metadata/Cohort2Metadata.RData") #metadd,cohort1metadd,annotation.file,compdataLU

u_dat <- merge(u_dat , unique(metadd%>%dplyr::select(Patient.Study.ID,dynamic_class3)) , by="Patient.Study.ID")
u_dat[,dynamic_class3:=dynamic_class3.y]
u_dat[,dynamic_class3.x:=NULL]
table(unique(u_dat%>%dplyr::select(Patient.Study.ID,dynamic_class3)))
u_dat$Patient.Study.ID <-as.factor(u_dat$Patient.Study.ID)






NEWXDAT2 <-rbind( data.table( unique(u_dat%>%dplyr::select(Patient.Study.ID,ARM,dynamic_class3)),day_fact="0"),
                  data.table( unique(u_dat%>%dplyr::select(Patient.Study.ID,ARM,dynamic_class3)),day_fact="14"),
                  data.table( unique(u_dat%>%dplyr::select(Patient.Study.ID,ARM,dynamic_class3)),day_fact="180"))
NEWXDAT2$Patient.Study.ID <-as.factor(NEWXDAT2$Patient.Study.ID)

Mod1 <- lmer(ssGSEA ~ -1+dynamic_class3*day_fact*ARM + (1+day_fact|Patient.Study.ID), REML=FALSE,data= u_dat) 
summary(Mod1)
PRED.lme4 <-  data.table(y = predict(Mod1, newdata=NEWXDAT2,re.form = NA) ,NEWXDAT2)
PRED.lme4 <-  data.table(y = predict(Mod1, newdata=NEWXDAT2,re.form = NA) ,
                         yi = predict(Mod1, newdata=NEWXDAT2,re.form =~(1+day_fact|Patient.Study.ID)) ,NEWXDAT2)

cinfreg<-data.table( predictInterval(Mod1,newdata=PRED.lme4, which ="fixed", level = 0.95) )

PRED.lme4<-data.table( data.table(cbind(PRED.lme4,cinfreg)) %>% group_by(ARM,day_fact,dynamic_class3)%>%dplyr::mutate(lwr=min(fit),upr=max(fit),fit=mean(fit) ))
PRED.lme4<-data.table(PRED.lme4 %>% group_by(ARM,day_fact,dynamic_class3)%>%dplyr::mutate(lwr=min(fit),upr=max(fit),fit=mean(fit) ))
PRED.lme4[,Intermediate_Clusters_ML:=u_dat[1]$Celltype]
PRED.lme4[,is_day0:=0];PRED.lme4[day_fact=="0",is_day0:=1]
PRED.lme4[, y_init:=sum(is_day0*yi), by=c("Patient.Study.ID","ARM","dynamic_class3")]
PRED.lme4[, mu_init:=sum(is_day0*y), by=c("Patient.Study.ID","ARM","dynamic_class3")]

PRED.lme4[, upr1:=quantile(yi-y_init,probs =0.75), by=c("day_fact","ARM","dynamic_class3")]
PRED.lme4[, lwr1:=quantile(yi-y_init,probs =0.25), by=c("day_fact","ARM","dynamic_class3")]
PRED.lme4[,uprminlwrhalf:=(upr1-lwr1)/2]
PRED.lme4[,upr1:=fit-mu_init+uprminlwrhalf]
PRED.lme4[,lwr1:=fit-mu_init-uprminlwrhalf]

PRED.lme4[,dum_day:=0]
PRED.lme4[day_fact==14,dum_day:=1]
PRED.lme4[day_fact==180,dum_day:=2]

PRED.lme4[,Treat:="Letrozole alone"]
PRED.lme4[ARM=="B",Treat:="Intermittent high \n dose ribociclib"]
PRED.lme4[ARM=="C",Treat:="Continuous low \n dose ribociclib"]
PRED.lme4$Treat <- factor(PRED.lme4$Treat, levels = c("Letrozole alone", "Intermittent high \n dose ribociclib", "Continuous low \n dose ribociclib"))
PRED.lme4[, mu_i_init:=sum(is_day0*yi), by=c("Patient.Study.ID","ARM","dynamic_class3")]

ggplot(PRED.lme4,
       aes(y=y-mu_init,x=dum_day,group=dynamic_class3,col=dynamic_class3))+
  geom_line(aes(group=interaction(dynamic_class3,ARM)),size=2,alpha=1)+#geom_line(aes(y=yi-y_init,group=interaction(pathwayscore,dynamic_class3,Patient.Study.ID)),size=1,alpha=1)+
  geom_ribbon(aes(ymax=upr1,ymin=lwr1,fill=dynamic_class3),alpha=0.4,col=NA)+
  facet_grid(.~Treat,scales="free") + 
  theme_classic(base_size=26)+
  scale_color_npg(name="Response",labels=c("Resistant","Sensitive"))+
  scale_fill_npg(name="Response",labels=c("Resistant","Sensitive"))+
  scale_x_continuous(name="Day",breaks=0:2,labels=c(0,14,180))+
  ylab("Hallmark estrogen response early") + 
  theme(aspect.ratio=1)+
  geom_line(aes(y=yi-mu_i_init,x=dum_day,group=Patient.Study.ID),linetype=2)
ggsave("/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/Lab Facilitation/Validation ERBB signaling upregulated in ET resistant tumors.pdf")

ggplot(PRED.lme4,
       aes(y=y-mu_init,x=dum_day,group=dynamic_class3,col=dynamic_class3))+
  geom_line(aes(group=interaction(dynamic_class3,ARM)),size=2,alpha=1)+#geom_line(aes(y=yi-y_init,group=interaction(pathwayscore,dynamic_class3,Patient.Study.ID)),size=1,alpha=1)+
  geom_ribbon(aes(ymax=upr1,ymin=lwr1,fill=dynamic_class3),alpha=0.4,col=NA)+
  facet_grid(.~Treat,scales="free") + 
  theme_classic(base_size=26)+
  scale_color_npg(name="Response",labels=c("Resistant","Sensitive"))+
  scale_fill_npg(name="Response",labels=c("Resistant","Sensitive"))+
  scale_x_continuous(name="Day",breaks=0:2,labels=c(0,14,180))+
  ylab("Hallmark estrogen response early") + 
  theme(aspect.ratio=1)#+

ggplot(PRED.lme4,
       aes(y=yi-mu_i_init,x=dum_day,group=Patient.Study.ID,col=dynamic_class3))+
  geom_line(aes(group=interaction(dynamic_class3,ARM,Patient.Study.ID)),size=2,alpha=1)+#geom_line(aes(y=yi-y_init,group=interaction(pathwayscore,dynamic_class3,Patient.Study.ID)),size=1,alpha=1)+
  theme_classic(base_size=26)+
  scale_color_npg(name="Response",labels=c("Resistant","Sensitive"))+
  scale_fill_npg(name="Response",labels=c("Resistant","Sensitive"))+
  scale_x_continuous(name="Day",breaks=0:2,labels=c(0,14,180))+
  ylab("EGF pathway activation") + 
  theme(aspect.ratio=1)


PRED.lme4[, TreatLAB:= "Letrozole alone"]
PRED.lme4[ARM!="A", TreatLAB:= "Combination ribociclib"]
PRED.lme4$TreatLAB <- factor(PRED.lme4$TreatLAB , levels=c( "Letrozole alone","Combination ribociclib" ) )
ggplot(PRED.lme4,
       aes(y=yi-mu_i_init,x=dum_day,group=Patient.Study.ID,col=TreatLAB))+
  geom_line(aes(group=interaction(TreatLAB,ARM,Patient.Study.ID)),size=2,alpha=1)+#geom_line(aes(y=yi-y_init,group=interaction(pathwayscore,dynamic_class3,Patient.Study.ID)),size=1,alpha=1)+
  theme_classic(base_size=26)+
  scale_color_jco(name="Treatment")+#,labels=c("Resistant","Sensitive"))+
  scale_fill_jco(name="Treatment")+#name="Treatment",labels=c("Resistant","Sensitive"))+
  scale_x_continuous(name="Day",breaks=0:2,labels=c(0,14,180))+
  ylab("ERBB pathway activation") + 
  theme(aspect.ratio=1)
ggsave("/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/Lab Facilitation/Validation ERBB signaling upregulated in ET resistant tumors.pdf")


u_dat[,responseLab:="Sensitive"]
u_dat[dynamic_class3=="Non-response",responseLab:="Resistant"]
u_dat[,HasDay0:=(sum(Day==0))>30, by=Patient.Study.ID]
u_dat[,HasDay14:=(sum(Day==14))>30, by=Patient.Study.ID]
u_dat[,HasDay180:=(sum(Day==180))>30, by=Patient.Study.ID]

u_dat[,nDay0:=(sum(Day==0)), by=Patient.Study.ID]
u_dat[,nDay14:=(sum(Day==14)), by=Patient.Study.ID]
u_dat[,nDay180:=(sum(Day==180)), by=Patient.Study.ID]

erbphenodat<-u_dat[HasDay0==T&HasDay14==T&HasDay180==T]#[Day%in%c(0,180)]

erbphenodat[, TreatLAB:= "Letrozole alone"]
erbphenodat[ARM!="A", TreatLAB:= "Combination ribociclib"]
erbphenodat$TreatLAB <- factor(erbphenodat$TreatLAB , levels=c( "Letrozole alone","Combination ribociclib" ) )
prop_EST_highdat1<-data.table(erbphenodat%>%group_by(Patient.Study.ID,dynamic_class3,TreatLAB,Day)%>%
                                summarise(mean_ERB_high= median(NAGASHIMA_EGF_SIGNALING_UP ),
                                          N=n() ))
prop_EST_highdat1[, initERB:= sum((Day==0)*mean_ERB_high), by=Patient.Study.ID ]
prop_EST_highdat1[, change:=mean_ERB_high -initERB]
erbphenodat$Patient.Study.ID <- factor(erbphenodat$Patient.Study.ID , levels= prop_EST_highdat1[Day==180][order(-change)]$Patient.Study.ID)
prop_EST_highdat1$Patient.Study.ID <- factor(prop_EST_highdat1$Patient.Study.ID , levels= prop_EST_highdat1[Day==180][order(-change)]$Patient.Study.ID)
prop_EST_highdat1[order(Patient.Study.ID,-change)]
ggplot(erbphenodat,
       aes(y=NAGASHIMA_EGF_SIGNALING_UP,x=Patient.Study.ID , fill=as.factor(Day),group=interaction(Patient.Study.ID,Day) ))+
  theme_classic(base_size=26)+ geom_boxplot()+stat_boxplot(geom="errorbar")+
  geom_point(size=0.1,alpha=0.1,position=position_dodge(width=0.8))+
  theme(aspect.ratio=.5)+
  scale_color_manual(name="Timepoint",values=(pal_aaas("default")(3))[c(3,1,2)], labels=c("pre-treatment","followup","post-treatment"))  +
  scale_fill_manual(name="Timepoint",values=(pal_aaas("default")(3))[c(3,1,2)],labels=c("pre-treatment","followup","post-treatment"))  +
  labs(y="ERBB family pathway activation \n (Composite ERBB response signature)",
       x="Patient")+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/ERBB phenotype shift pre to post treatment1.pdf",width=14, height=10)
ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Validation ERBB phenotype shift pre early and post treatment.pdf",width=14, height=10)



#######
load(       file="~/Dropbox/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/ERBB median phenotype shift pre to post treatment.RData")
#ERBB_mu_save , res.pca,ERBBsetlist
u_dat$PCA2 <- predict( res.pca,newdata=data.frame(u_dat %>% dplyr::select(ERBBsetlist  ) ))[,2]


erbphenodat <- u_dat[nDay0>20&nDay180>20][Day%in%c(0,180)]
prop_EST_highdat1 <- data.table(erbphenodat%>%group_by(Patient.Study.ID,dynamic_class3,ARM,Day)%>%
                                  summarise(mean_ERB_high= median(PCA2),
                                            medan_NAGASHERB_high= median(NAGASHIMA_EGF_SIGNALING_UP),
                                            N=n() ))
prop_EST_highdat1[, initERB:= sum((Day==0)*mean_ERB_high), by=Patient.Study.ID ]
prop_EST_highdat1[, change:=mean_ERB_high -initERB]
erbphenodat$Patient.Study.ID <- factor(erbphenodat$Patient.Study.ID , levels= prop_EST_highdat1[Day==180][order(-change)]$Patient.Study.ID)
prop_EST_highdat1$Patient.Study.ID <- factor(prop_EST_highdat1$Patient.Study.ID , levels= prop_EST_highdat1[Day==180][order(-change)]$Patient.Study.ID)
prop_EST_highdat1[order(Patient.Study.ID,-change)]

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

#fread(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/DiscoveryCohortERBBPrePost.csv")$Treatment%>%unique()
erbphenodat[,Treatment:="letrozole + ribo"]
erbphenodat[ARM=="A",Treatment:="letrozole"]
discData <- data.table( erbphenodat%>%dplyr::select(Cell.ID,Patient.Study.ID,dynamic_class3,Day,ARM,Treatment,Celltype,PCA2,        responseLab) )
setnames( discData , old=c("PCA2"), new=c("ERBBaxis"))
#write.csv(discData,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/ValidationCohortERBBPrePost.csv")

summary( lm(ERBBaxis~ -1+Patient.Study.ID+ Day : Patient.Study.ID ,data=discData[]))
summary( lm(ERBBaxis~ -1+Patient.Study.ID+ Day : Patient.Study.ID ,data=fread(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/DiscoveryCohortERBBPrePost.csv")))
summary( lm(ERBBaxis~ -1+Patient.Study.ID+ Day : Patient.Study.ID ,data=fread(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/ValidationCohortERBBPrePost.csv")))


setnames(ERBB_mu_save , old="ARM.x",new="ARM")
ERBB_mu_save[,Cohort:="Discovery"]
ERBB_mu_saveC2[,Cohort:="Validation"]
allAverageERBB <- rbind(ERBB_mu_save , ERBB_mu_saveC2)
allAverageERBB[, initERB:= sum((Day==0)*median_ERB_high), by=Patient.Study.ID ]
allAverageERBB[, initERBNAG:= sum((Day==0)*median_NAGASHERB_high), by=Patient.Study.ID ]
allAverageERBB[, change:=median_ERB_high -initERB]
allAverageERBB[, changeNAG:=median_NAGASHERB_high -initERBNAG]
allAverageERBB$Patient.Study.ID <- factor(allAverageERBB$Patient.Study.ID , levels= allAverageERBB[Day==180][order(-change)]$Patient.Study.ID)
allAverageERBB[order(Patient.Study.ID,-change)]

allAverageERBB[,muinitERB:=mean(initERB) , by=Cohort ]
allAverageERBB[,muinitERBNAG:=mean(initERBNAG) , by=Cohort ]

allAverageERBB[,TreatLAB:="Combination ribociclib"]
allAverageERBB[ARM=="A",TreatLAB:="Letrozole alone"]

allAverageERBB[, nTot:= sum(N), by=Patient.Study.ID]
allAverageERBB[, n0:= nTot-N]
ggplot(allAverageERBB[Day!=0][], aes(y= change, x= Cohort,group= interaction( Cohort,as.factor(Day))  ))+
  theme_classic(base_size=26)+ geom_boxplot()+
  geom_hline(linetype="dashed", yintercept=0)+
  geom_boxplot(fill="red",outlier.colour=NA)+
  geom_point(aes(size=N))+
  labs(y="ERBB pathway activation \n (post treatment median  \n relative to baseline)",x="Cohort")+
  theme(aspect.ratio=1)+
  scale_size_continuous(name="Cancer cells \n sampled")  

ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/ERBB phenotype shift pre and post treatment DiscoveryandValidation Rena.pdf",width=10, height=10)



