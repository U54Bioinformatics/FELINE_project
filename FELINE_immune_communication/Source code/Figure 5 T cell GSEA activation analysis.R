rm(list=ls())
require(mgcv);require(data.table);require(dplyr);require(ggplot2);require(tidyr);require(lme4);require(lmerTest);require(parallel)

# load T cell ssGSEA pathway activation signature
savloc <- "/Users/jason/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/Figure5/Outputs/"
TcellPhenotytpedata <- data.table(read.csv(file=paste0(savloc,"SourceData_Figure5_TcellssGSEAActivationPhenotpe_Output")))

TcellPhenotytpedata$Patient.Study.ID%>%unique()%>%length() #54
TcellPhenotytpedata%>%dplyr::select(Patient.Study.ID,Day)%>%unique()%>%nrow() #116
TcellPhenotytpedata%>%dplyr::select(Patient.Study.ID,Day,Cohort)%>%group_by(Cohort)%>%summarise(n=n())#unique()%>%nrow() #116
TcellPhenotytpedata%>%dplyr::select(Patient.Study.ID,Day,Cohort)%>%group_by(Cohort)%>%summarise(n=length(unique(paste0(Patient.Study.ID,Day))))
TcellPhenotytpedata%>%dplyr::select(Patient.Study.ID,Day,Cohort)%>%group_by(Cohort)%>%summarise(n=length(unique(Patient.Study.ID)))


ggplot( TcellPhenotytpedata[] , aes(y=Tcellactivation, x=log(1+Day), col=dynamic_class3, fill=dynamic_class3, group=interaction(dynamic_class3) ))+
  theme_classic(base_size=18)+facet_grid(Treatmentlab~.)+
  geom_boxplot(alpha=0.6,aes(group=interaction(dynamic_class3, Day)),position=position_dodge() )+#geom_smooth(method="lm")+
  geom_point(aes(shape=Cohort),position =position_jitterdodge(dodge.width=2.3,jitter.width=0.4)  )+
  theme(aspect.ratio=1)+
  scale_x_continuous(breaks=log(1+c(0,14,180)) ,labels=c(0,14,180))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant","Sensitive"))+
  scale_color_npg(name="Tumor response", labels=c("Resistant","Sensitive"))+
  labs(y="CD8 T cell  cytotoxicity phenotype \n (GSE22886 NK vs naive ssGSEA signature)", x="Day") +
  geom_smooth(method="gam", formula=y~s(x,k=3),se=F)

ggplot( TcellPhenotytpedata[] , aes(y=Tcellactivation, x=log(1+Day), col=dynamic_class3, fill=dynamic_class3, group=interaction(dynamic_class3) ))+
  theme_classic(base_size=26)+facet_wrap(~Treatmentlab,nrow=2)+
  stat_boxplot(geom ='errorbar',position=position_dodge(),aes(group=interaction(dynamic_class3, Day)) ) + 
  geom_boxplot(alpha=0.6,outlier.shape=NA,aes(group=interaction(dynamic_class3, Day)),position=position_dodge() )+#geom_smooth(method="lm")+
  geom_point(aes(shape=Cohort), position =position_jitterdodge(dodge.width=2.3,jitter.width=0.4)  )+
  theme(aspect.ratio=1)+scale_x_continuous(breaks=log(1+c(0,14,180)) ,labels=c(0,14,180))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant","Sensitive"))+
  scale_color_npg(name="Tumor response", labels=c("Resistant","Sensitive"))+
  labs(y="CD8+ T cell  cytotoxicity phenotype \n (GSE22886 NK vs naive ssGSEA signature)", x="Day") 


ribomod<-lmer(Tcellactivation~(dynamic_class3*Day*Cohort+ (1|Patient.Study.ID)) , data= TcellPhenotytpedata[Treatmentlab=="Combination ribociclib"][ ] )
plot(ribomod)
qqnorm(residuals(ribomod))
coef(summary(ribomod))
plot(predict(ribomod),TcellPhenotytpedata[Treatmentlab=="Combination ribociclib"]$Tcellactivation)
letromod<-lmer(Tcellactivation~(dynamic_class3*Day*Cohort+ (1|Patient.Study.ID)) , data= TcellPhenotytpedata[Treatmentlab=="Letrozole alone"][ ] )
plot(letromod)
qqnorm(residuals(letromod))
coef(summary(letromod))
plot(predict(letromod),TcellPhenotytpedata[Treatmentlab=="Letrozole alone"]$Tcellactivation)

lmer(Tcellactivation~dynamic_class3*Day*Cohort+ (1|Patient.Study.ID) , data= TcellPhenotytpedata[][][Treatmentlab=="Combination ribociclib"] )%>%summary()
lmer(Tcellactivation~dynamic_class3*Day*Cohort+ (1|Patient.Study.ID) , data= TcellPhenotytpedata[][][Treatmentlab=="Letrozole alone"] )%>%summary()

#save(TcellPhenotytpedata, file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/Figure5/Discovery and Validation T cell cytotoxicity phenotype.RData")
