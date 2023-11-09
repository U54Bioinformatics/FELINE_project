rm(list=ls())
library(glmnet); require(dplyr); require(data.table); require(tidyr)
require(ggplot2); require(pheatmap); require(ggsci); require(lme4);require(lmerTest)
library("stringr")

discdat<- fread(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/Figure2/DiscoveryCohortERBBPrePost.csv")
disc_EST<-data.table(discdat%>%group_by(Patient.Study.ID,dynamic_class3,ARM ,Day)%>%
                                summarise(mean_ERB_high= median(ERBBaxis), N=n() ))
disc_EST[, initERB:= sum((Day==0)*mean_ERB_high), by=Patient.Study.ID ]
disc_EST[, change:=mean_ERB_high -initERB]
discdat$Patient.Study.ID <- factor(discdat$Patient.Study.ID , levels= disc_EST[Day==180][order(-change)]$Patient.Study.ID)





validdat<- fread(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/Figure2/ValidationCohortERBBPrePost.csv")
valid_EST<-data.table(validdat%>%group_by(Patient.Study.ID,dynamic_class3,ARM ,Day)%>%
                       summarise(mean_ERB_high= median(ERBBaxis), N=n() ))
valid_EST[, initERB:= sum((Day==0)*mean_ERB_high), by=Patient.Study.ID ]
valid_EST[, change:=mean_ERB_high -initERB]
validdat$Patient.Study.ID <- factor(validdat$Patient.Study.ID , levels= valid_EST[Day==180][order(-change)]$Patient.Study.ID)

D_codes<-formatC(1:length(discdat$Patient.Study.ID %>%unique()), width = 2, format = "d", flag = "0")
V_codes<-formatC(1:length(validdat$Patient.Study.ID %>%unique()), width = 2, format = "d", flag = "0")

D_ids <- data.table( Patient.Study.ID= levels(discdat$Patient.Study.ID), PlotID=paste0("D",D_codes) )
V_ids <- data.table( Patient.Study.ID= levels(validdat$Patient.Study.ID), PlotID=paste0("V",V_codes) )
erbphenodat <- rbind(data.table(Cohort="Discovery cohort", merge(D_ids,discdat, by="Patient.Study.ID")),
                     data.table(Cohort="Validation cohort",merge(V_ids,validdat, by="Patient.Study.ID")))

ggplot(erbphenodat,
       aes(y=ERBBaxis,x=PlotID , fill=as.factor(Day),group=interaction(PlotID,Day) ))+
  theme_classic(base_size=26)+ geom_boxplot()+stat_boxplot(geom="errorbar")+
  geom_jitter(size=0.4,alpha=0.6,position=position_dodge(width=0.8))+
  theme(aspect.ratio=0.5, legend.position = "none")+
  scale_color_manual(name="Timepoint",values=c("yellow","red"), labels=c("pre-treatment","post-treatment"))  +
  scale_fill_manual(name="Timepoint",values=c("yellow","red"),labels=c("pre-treatment","post-treatment"))  +
  labs(y="Cancer cell ERBB pathway activation",
       x="Patients")+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~Cohort,scales="free",ncol=1)+
  theme(strip.text = element_text(face="bold"))
#ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Discovery and Validation combined ERBB phenotype shift pre early and post treatment1.png",width=10, height=12, dpi=320)


erbphenodat[, c("FEL_ID","TimePoint","cellcode") := data.table(str_split_fixed(Cell.ID,"_",3))]
PatCodes <- unique( erbphenodat%>%dplyr::select(Cohort, Patient.Study.ID, FEL_ID,PlotID) )

#write.csv(PatCodes,  file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/Figure2/CrossreferencePairedCancerPatientTumorIDcodes.csv")
