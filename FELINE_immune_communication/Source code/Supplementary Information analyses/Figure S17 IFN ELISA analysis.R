rm(list=ls())
require(data.table)
require(ggplot2)
require(tidyr)
require(dplyr)
require(ggsci)
savloc<-"/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/SI data/FigureS17/"
dd<- data.table(read.csv(paste0(savloc,"SourceData_Figure S17 IFN gamma ELISA data.csv")))


#dd<- data.table(read.csv("C:/Users/pcosgrove/OneDrive - City of Hope National Medical Center/Desktop/Griffiths_V3_ImmunePaper Response/20240909_Exp61_IFNg ELISA/20240913_Exp61_IFNg_Data_R.csv"))
dd$Composition<-factor(dd$Composition, levels=c("CAMA-1 Sensitive V2", "T Cell", "CAMA-1 Sensitive V2 + T Cell"))
dd$Treatment <- gsub("1uM", "1μM", dd$Treatment)
dd$TreatmentLab<-dd$Treatment
dd[Treatment=="1μM Ribociclib + 5ng/mL IL-15",TreatmentLab:="1μM Ribociclib \n + 5ng/mL IL-15" ]
dd$TreatmentLab<-factor(dd$TreatmentLab, levels=c("DMSO"    , "1μM Ribociclib"  , "5ng/mL IL-15", "1μM Ribociclib \n + 5ng/mL IL-15"))

dd[,Treat:="DMSO"]
dd[Treatment=="1μM Ribociclib" ,Treat:="ribociclib"]
dd[Treatment=="1μM Ribociclib + 5ng/mL IL-15" ,Treat:="combination"]
dd[Treatment=="5ng/mL IL-15",Treat:="IL15"]
dd$Treat<-factor(dd$Treat, levels=c("DMSO"    , "ribociclib"  , "IL15", "combination"))

dd[,Compos:="CT"]
dd[Composition=="CAMA-1 Sensitive V2" ,Compos:="C"]
dd[Composition=="T Cell" ,Compos:="T"]
dd$Compos<-factor(dd$Compos, levels=c("C"    , "T"  , "CT"))


lm1<-lm(log2(Conc_pgmL)~Compos*Treat ,data=dd)
plot(lm1)
summary(lm1)
dd$use=2^(predict(lm1, se=T)$fit+predict(lm1, se=T)$se.fit)
dd$lse=2^(predict(lm1, se=T)$fit-predict(lm1, se=T)$se.fit)

ggplot( dd , aes(y=Conc_pgmL, x=TreatmentLab, group=interaction(Treatment,Composition), fill=Composition))+
  geom_bar(data=dd%>% group_by(Treatment,TreatmentLab,Composition)%>%summarise(Conc_pgmL=2^(mean(log2(Conc_pgmL)))),position=position_dodge(),stat="identity")+
  theme_classic()+theme(aspect.ratio=1)+
  labs(y="Human IFN-Gamma (pg/mL)",x="Treatment")+
  #stat_summary(position=position_dodge(), 
   #                   geom = "errorbar") + 
  geom_errorbar(position=position_dodge(),aes(ymax=use, ymin=lse))+
  geom_point(position=position_dodge(width=1),pch=21)+
  scale_fill_npg(labels=c("Cancer","T cells","Cancer+ \n T cells"))
ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/Griffiths_V3_ImmunePaper Response/Additional Files/IFNElisabyCellCompositionAndTreatmentpgml.pdf", height=5.5,width=5.5)
ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/Griffiths_V3_ImmunePaper Response/Additional Files/IFNElisabyCellCompositionAndTreatmentpgml.png", dpi=320,height=5.5,width=5.5)
