rm(list=ls())
require(mgcv);require(data.table);require(dplyr);require(ggplot2);require(tidyr);require(lme4);require(lmerTest);require(parallel)
load(file= "/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/temporary/Validation T cell cytotoxicity phenotype.RData")
VTcellPhenotytpedata <- TcellPhenotytpedata
  
load(file= "/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/temporary/Discovery T cell cytotoxicity phenotype.RData")
DTcellPhenotytpedata <- TcellPhenotytpedata %>% select(-Tcellactivation2)
setnames(VTcellPhenotytpedata, old="orig.ident", new="Sample_p_t")
TcellPhenotytpedata <- rbind( data.table(Cohort="Discovery",DTcellPhenotytpedata %>% dplyr::select(names(DTcellPhenotytpedata))),
                              data.table(Cohort="Validation",VTcellPhenotytpedata %>% dplyr::select(names(DTcellPhenotytpedata))))
                                                                                         
                                                                                                
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
paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Discovery and Validation Ribo and Letrozole Activated vs Naive T cell cytotoxicity phenotype over time boxplot.pdf"),height=10,width=10, dpi=320)

save(TcellPhenotytpedata, file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/Figure5/Discovery and Validation T cell cytotoxicity phenotype.RData")

lmer(Tcellactivation~dynamic_class3*Day*Cohort+ (1|Patient.Study.ID) , data= TcellPhenotytpedata[][][Treatmentlab=="Combination ribociclib"] )%>%summary()
lmer(Tcellactivation~dynamic_class3*Day*Cohort+ (1|Patient.Study.ID) , data= TcellPhenotytpedata[][][Treatmentlab=="Letrozole alone"] )%>%summary()


