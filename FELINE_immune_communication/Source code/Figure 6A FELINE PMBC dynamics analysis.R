rm(list=ls())
require(data.table)
require(dplyr)
require(tidyr)
require(ggplot2)
require(boot)
require("compositions")
require(vegan)
require(ggsci)

## Start by getting the M1-M2 phenotypes of myeloid cells
require(mgcv);require(lme4);require(lmerTest);require(parallel)

savloc<- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/Figure6/"
WBCdata <- data.table( read.csv(file=paste0(savloc,"SourceData_Figure6_FELINE PBMC.csv")) )

g1<- gam(I(White.Blood.Cells..WBC..Blood.K.UL)*1000~s(Day,k=4,by=as.factor(Treatmentlab)), data=WBCdata)
summary(g1)
plot(g1)
g2<- gam(I(log(White.Blood.Cells..WBC..Blood.K.UL*1000))~s(sqrt(Day),k=4,by=as.factor(Treatmentlab)), data=WBCdata)
summary(g2)
plot(g2)

m1<- lm(I(log(White.Blood.Cells..WBC..Blood.K.UL*1000))~Treatmentlab, data=WBCdata[Day==0])
summary(m1)

resout<-WBCdata[! is.na(dynamic_class3)][order(Day)][]%>%dplyr::select(Patient.Study.ID,Treatmentlab,Day,White.Blood.Cells..WBC..Blood.K.UL)
write.csv(resout,file=paste0(savloc,"Outputs/","SourceData_Figure6_FELINE PBMC_Out.csv"))

ggplot(resout, aes(y=(White.Blood.Cells..WBC..Blood.K.UL)*1000, x=Day, col=Treatmentlab))+
  #geom_path(aes(group=Patient.Study.ID)) +
  geom_point(position=position_dodge(width=10),size=2.5) +
  geom_smooth(method = "gam",  se = T, formula=y~s(x, k=3),
              method.args = list(family = "poisson"))+
  coord_trans(y="log")+
  theme_classic(base_size=20)+theme(aspect.ratio=1)+
  labs(y="Peripheral white blood cells (per ul)")+
  scale_color_manual(name="Treatment", values=rev(pal_jco("default")(2)))


ggplot(WBCdata[! is.na(dynamic_class3)][order(Day)][], aes(y=(White.Blood.Cells..WBC..Blood.K.UL), x=sqrt(Day), fill=Treatmentlab,col=Treatmentlab))+
  #geom_path(aes(group=Patient.Study.ID)) +
  geom_smooth(method = "gam",  se = T, alpha=0.6, formula=y~s(x, k=3))+
  geom_point()+#position=position_dodge(width=10),size=2.5) +
  coord_trans(y="log",x = scales::trans_new("square", function(x) x^2, "sqrt"))+
  #coord_trans(x="log")+
  theme_classic(base_size=26)+theme(aspect.ratio=1)+
  labs(y="Peripheral white blood cells (K/ul)",x="Day")+
  scale_color_manual(name="Treatment", values=rev(pal_jco("default")(2)))+
  scale_fill_manual(name="Treatment", values=rev(pal_jco("default")(2))) +
  scale_y_continuous(breaks=c(2,4,8,16,32),labels=c(2,4,8,16,32))+
  scale_x_continuous(breaks=sqrt(c(0,90,180)),labels=c(0,90,180))
paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Ribo and Letrozole Peripheral white blood cell count by treatment over time.png"),height=10,width=10)

ggplot(WBCdata[! is.na(dynamic_class3)][order(Day)][], aes(y=(White.Blood.Cells..WBC..Blood.K.UL), x=sqrt(Day), fill=Treatmentlab,col=Treatmentlab))+
  #geom_path(aes(group=Patient.Study.ID)) +
  geom_smooth(method = "gam",  se = T, alpha=0.6, formula=y~s(x, k=3))+
  geom_point(size=3.5)+#position=position_dodge(width=10),size=2.5) +
  coord_trans(y="log",x = scales::trans_new("square", function(x) x^2, "sqrt"))+
  #coord_trans(x="log")+
  theme_classic(base_size=26)+theme(aspect.ratio=1)+
  labs(y="Peripheral white blood cells (K/ul)",x="Day")+
  scale_color_manual(name="Treatment", values=rev(pal_jco("default")(2)))+
  scale_fill_manual(name="Treatment", values=rev(pal_jco("default")(2))) +
  scale_y_continuous(breaks=c(2,4,8,16,32),labels=c(2,4,8,16,32))+
  scale_x_continuous(breaks=sqrt(c(0,90,180)),labels=c(0,90,180))+
  theme(axis.title=element_blank(),
        axis.text=element_blank(), legend.position="none")
paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
ggsave(paste0(paperfile,"BLANK Ribo and Letrozole Peripheral white blood cell count by treatment over time.png"),height=10,width=10)


ggplot(WBCdata[! is.na(dynamic_class3)][order(Day)][], aes(y=(White.Blood.Cells..WBC..Blood.K.UL), x=sqrt(Day), fill=dynamic_class3,col=dynamic_class3))+
  #geom_path(aes(group=Patient.Study.ID)) +
  geom_smooth(method = "gam",  se = T, alpha=0.6, formula=y~s(x, k=3))+
  geom_point()+#position=position_dodge(width=10),size=2.5) +
  coord_trans(y="log",x = scales::trans_new("square", function(x) x^2, "sqrt"))+
  #coord_trans(x="log")+
  theme_classic(base_size=26)+theme(aspect.ratio=1)+
  labs(y="Peripheral white blood cells (K/ul)",x="Day")+
  scale_color_npg(name="Response", labels=c("Resistant","Sensitive"))+
  scale_fill_npg(name="Response", labels=c("Resistant","Sensitive")) +
  scale_y_continuous(breaks=c(2,4,8,16,32),labels=c(2,4,8,16,32))+
  scale_x_continuous(breaks=sqrt(c(0,90,180)),labels=c(0,90,180))+
  facet_wrap(~Treatmentlab, ncol=1)
#ggsave(paste0(paperfile,"Ribo and Letrozole Peripheral white blood cell count by response and treatment over time.png"),height=10,width=10)




ggplot(WBCdata[order(Day)][Monocytes..Absolute.Monocyte.Count..Blood..K.UL<5], aes(y=(Monocytes..Absolute.Monocyte.Count..Blood..K.UL), x=sqrt(Day), fill=dynamic_class3,col=dynamic_class3))+
  #geom_path(aes(group=Patient.Study.ID)) +
  #geom_smooth(method = "gam",  se = T, alpha=0.6, formula=y~s(x, k=3))+
  geom_point()+#position=position_dodge(width=10),size=2.5) +
  coord_trans(y="log",x = scales::trans_new("square", function(x) x^2, "sqrt"))+
  #coord_trans(x="log")+
  theme_classic(base_size=26)+theme(aspect.ratio=1)+
  labs(y="Peripheral white blood cells (K/ul)",x="Day")+
  scale_color_npg(name="Response", labels=c("Resistant","Sensitive"))+
  scale_fill_npg(name="Response", labels=c("Resistant","Sensitive")) +
  scale_y_continuous(breaks=c(2,4,8,16,32),labels=c(2,4,8,16,32))+
  scale_x_continuous(breaks=sqrt(c(0,90,180)),labels=c(0,90,180))+
  facet_wrap(~Treatmentlab, ncol=1)
#ggsave(paste0(paperfile,"Ribo and Letrozole Peripheral white blood cell count by response and treatment over time.png"),height=10,width=10)



ggplot(WBCdata[order(Day)][], aes(y=(Lymphocytes..Absolute.Lymph.Count..Blood.K.UL), x=sqrt(Day), fill=Treatmentlab,col=Treatmentlab))+
  #geom_path(aes(group=Patient.Study.ID)) +
  geom_smooth(method = "gam",  se = T, alpha=0.6, formula=y~s(x, k=3))+
  geom_point()+#position=position_dodge(width=10),size=2.5) +
  coord_trans(y="log",x = scales::trans_new("square", function(x) x^2, "sqrt"))+
  #coord_trans(x="log")+
  theme_classic(base_size=26)+theme(aspect.ratio=1)+
  labs(y="Peripheral lymphocytes (K/ul)",x="Day")+
  scale_color_manual(name="Treatment", values=rev(pal_jco("default")(2)))+
  scale_fill_manual(name="Treatment", values=rev(pal_jco("default")(2))) +
  scale_y_continuous(breaks=c(2,4,8,16,32),labels=c(2,4,8,16,32))+
  scale_x_continuous(breaks=sqrt(c(0,90,180)),labels=c(0,90,180))
paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Ribo and Letrozole Peripheral lymphocyte count by treatment over time.png"),height=10,width=10)

ggplot(WBCdata[order(Day)][Lymphocytes..Absolute.Lymph.Count..Blood.K.UL<25], aes(y=(Lymphocytes..Absolute.Lymph.Count..Blood.K.UL), x=sqrt(Day), fill=Treatmentlab,col=Treatmentlab))+
  #geom_path(aes(group=Patient.Study.ID)) +
  geom_smooth(method = "gam",  se = T, alpha=0.6, formula=y~s(x, k=3))+
  geom_point()+#position=position_dodge(width=10),size=2.5) +
  coord_trans(y="log",x = scales::trans_new("square", function(x) x^2, "sqrt"))+
  #coord_trans(x="log")+
  theme_classic(base_size=26)+theme(aspect.ratio=1)+
  labs(y="Peripheral lymphocytes (K/ul)",x="Day")+
  scale_color_manual(name="Treatment", values=rev(pal_jco("default")(2)))+
  scale_fill_manual(name="Treatment", values=rev(pal_jco("default")(2))) +
  scale_y_continuous(breaks=c(2,4,8,16,32),labels=c(2,4,8,16,32))+
  scale_x_continuous(breaks=sqrt(c(0,90,180)),labels=c(0,90,180))
paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Ribo and Letrozole Peripheral lymphocyte count without outliers by treatment over time.png"),height=10,width=10)



ggplot(WBCdata[order(Day)][], aes(y=(Monocytes..Absolute.Monocyte.Count..Blood..K.UL), x=sqrt(Day), fill=Treatmentlab,col=Treatmentlab))+
  #geom_path(aes(group=Patient.Study.ID)) +
  geom_smooth(method = "gam",  se = T, alpha=0.6, formula=y~s(x, k=3))+
  geom_point()+#position=position_dodge(width=10),size=2.5) +
  coord_trans(y="log",x = scales::trans_new("square", function(x) x^2, "sqrt"))+
  #coord_trans(x="log")+
  theme_classic(base_size=26)+theme(aspect.ratio=1)+
  labs(y="Peripheral monocytes (K/ul)",x="Day")+
  scale_color_manual(name="Treatment", values=rev(pal_jco("default")(2)))+
  scale_fill_manual(name="Treatment", values=rev(pal_jco("default")(2))) +
  scale_y_continuous(breaks=c(2,4,8,16,32),labels=c(2,4,8,16,32))+
  scale_x_continuous(breaks=sqrt(c(0,90,180)),labels=c(0,90,180))




