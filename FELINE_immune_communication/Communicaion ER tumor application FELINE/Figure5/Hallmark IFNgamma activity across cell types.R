rm(list=ls())
require(data.table);require(dplyr);require(tidyr)
require(ggplot2);require(ggsci)

require(lme4)
require(lmerTest)
require(merTools)
load(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/temporary/Discovery Hallmark IFNgamma TME.RData")
DiscoveryplotThis<-plotThis
load(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/temporary/Validation Hallmark IFNgamma TME.RData")
ValidationplotThis<-plotThis

plotThis <- rbind(data.table(Cohort="Discovery",DiscoveryplotThis%>%dplyr::select(Cell.ID,Patient.Study.ID,Day,Daylab,value,valuemuDay0,NormSSGSEA,TumorResponse,dynamic_class3,Treatmentlab,ARM,Celltype_subtype2,Celltype,Celltype_subtype) ),
                  data.table(Cohort="Validation",ValidationplotThis%>%dplyr::select(Cell.ID,Patient.Study.ID,Day,Daylab,value,valuemuDay0,NormSSGSEA,TumorResponse,Treatmentlab,dynamic_class3,ARM,Celltype_subtype2,Celltype,Celltype_subtype)  ))


plotThis[,valuemuCelltypeSenstumor:= sum(value*(dynamic_class3=="Response"))/sum(dynamic_class3=="Response"),by=c("Celltype_subtype2","Cohort")]
p_dd[is.na(valuemuDay0),valuemuDay0:=0]
p_dd[,NormSSGSEA:=value-valuemuDay0]


plotThis$Celltype_subtype2 <-  factor(plotThis$Celltype_subtype2, levels=rev(c("Cancer cells", "Diploid epithelial cells",
                                                                               "Myeloid cells", "CD8+ T cells", "CD4+ T cells",
                                                                               "Adipocytes" , "Fibroblasts", "Endothelial cells","Pericytes"  #"NK cells",
                                                                                  ) ))

plotThis[,Celltype_subtype3:=Celltype_subtype2]
plotThis[Celltype_subtype2%in%c("Adipocytes" , "Fibroblasts", "Endothelial cells","Pericytes"  ),Celltype_subtype3:="Stromal cells"]
plotThis$Celltype_subtype3 <-  factor(plotThis$Celltype_subtype3, levels=rev(c("Cancer cells", "Diploid epithelial cells",
                                                                               "Myeloid cells", "CD8+ T cells", "CD4+ T cells",
                                                                               "Stromal cells" 
) ))


HallmarkMean <- data.table( plotThis%>%group_by(Cohort,Celltype_subtype3,dynamic_class3,Patient.Study.ID,Day,Daylab,TumorResponse,Treatmentlab,ARM)%>%dplyr::summarise(value=mean(value),valuemuCelltypeSenstumor=mean(valuemuCelltypeSenstumor)))
ggplot( plotThis%>%group_by(Cohort,Celltype_subtype3,dynamic_class3,Patient.Study.ID,Day,Daylab,TumorResponse,Treatmentlab,ARM)%>%dplyr::summarise(value=mean(value),valuemuCelltypeSenstumor=mean(valuemuCelltypeSenstumor))
        , aes(y=value-valuemuCelltypeSenstumor,x=Celltype_subtype3, group=interaction(Celltype_subtype3,dynamic_class3),fill=dynamic_class3,col=dynamic_class3))+theme_classic(base_size=24)+
  geom_boxplot(col="black",outlier.colour=NA, position= position_dodge() )+
  stat_boxplot(geom="errorbar",position=position_dodge())+#geom_smooth(method="lm")+
  scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  geom_point(aes(shape=Cohort,group=interaction(Celltype_subtype3,dynamic_class3)),col="black",position =position_dodge(width=1)  ,size=3)+
  theme(aspect.ratio=1)+
  theme(axis.text.x=element_text(angle=90,vjust=0.5))+
  facet_wrap(~Treatmentlab, ncol=1)+
  labs(x= "Cell type",y="Interferon gamma response \n (Hallmark ssGSEA score)")+ coord_flip()
paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Discovery and Validation Ribo and Letrozole Hallmark IFN gamma response across cell types summarised.png"),height=10,width=10, dpi=320)


ggplot( plotThis  , aes(y=value-valuemuCelltypeSenstumor,x=Celltype_subtype3, group=interaction(Celltype_subtype3,dynamic_class3),fill=dynamic_class3,col=dynamic_class3))+theme_classic(base_size=24)+
  geom_boxplot(outlier.colour=NA, position= position_dodge(), alpha=0.6 )+
  stat_boxplot(geom="errorbar",position=position_dodge())+#geom_smooth(method="lm")+
  scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  geom_point(aes(shape=Cohort,group=interaction(Celltype_subtype3,dynamic_class3)),position =position_dodge(width=1)  ,size=1.75)+
  theme(aspect.ratio=1)+
  theme(axis.text.x=element_text(angle=90,vjust=0.5))+
  facet_wrap(~Treatmentlab, ncol=1)+
  labs(x= "Cell type",y="Interferon gamma response \n (Hallmark ssGSEA score)")+ coord_flip()
paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Discovery and Validation Ribo and Letrozole Hallmark IFN gamma response across cell types.png"),height=10,width=10, dpi=320)

save(plotThis,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/Figure5/Discover and Validation Hallmark IFNgamma TME.RData")


MC <- lmer(value-valuemuCelltypeSenstumor~ -1+ Celltype_subtype3+Celltype_subtype3 : dynamic_class3 + Cohort +(1|Patient.Study.ID) , data= plotThis[Treatmentlab=="Combination ribociclib"])
summary( MC )
ML <- lmer(value-valuemuCelltypeSenstumor~ -1+ Celltype_subtype3+Celltype_subtype3 : dynamic_class3 + Cohort +(1|Patient.Study.ID) , data= plotThis[Treatmentlab=="Letrozole alone"])
summary( ML )

MC <- lmer(value-valuemuCelltypeSenstumor~ dynamic_class3 + Cohort +(1|Patient.Study.ID) , data= plotThis[Celltype_subtype3=="Stromal cells"][Treatmentlab=="Combination ribociclib"])
summary( MC )

Combores<-rbindlist(lapply( unique(plotThis$Celltype_subtype3) , function(i){
  MC <- lmer(value-valuemuCelltypeSenstumor~ dynamic_class3 + Cohort +(1|Patient.Study.ID) , data= plotThis[Celltype_subtype3==i][Treatmentlab=="Combination ribociclib"])
  out<-data.table(Celltype_subtype3=i,coef(summary( MC )), keep.rownames=T)
  setnames(out, old=c( "Std. Error" , "t value","Pr(>|t|)" ), new=c( "Std.Error" , "tvalue","pvalue" ) )
  return(out)
}))
Letrozoleres<-rbindlist(lapply( unique(plotThis$Celltype_subtype3) , function(i){
  MC <- lmer(value-valuemuCelltypeSenstumor~ dynamic_class3 + Cohort +(1|Patient.Study.ID) , data= plotThis[Celltype_subtype3==i][Treatmentlab!="Combination ribociclib"])
  out<-data.table(Celltype_subtype3=i,coef(summary( MC )), keep.rownames=T)
  setnames(out, old=c( "Std. Error" , "t value","Pr(>|t|)" ), new=c( "Std.Error" , "tvalue","pvalue" ) )
  return(out)
}))

outall<-rbind( data.table(Treatmentlab="Combination ribociclib" ,Combores[rn=="dynamic_class3Response"]),
data.table(Treatmentlab="Letrozole alone" ,Letrozoleres[rn=="dynamic_class3Response"]))

write.csv(outall,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/temporary/Discovery and Validation stats Hallmark IFN activationTME.csv")



MCv <- lmer(value-valuemuCelltypeSenstumor~  dynamic_class3  +(1|Patient.Study.ID) , data= plotThis[Celltype_subtype2=="Cancer cells"][Cohort=="Validation"][Treatmentlab=="Combination ribociclib"])
summary( MCv )



MLv <- lmer(value-valuemuCelltypeSenstumor~  dynamic_class3  +(1|Patient.Study.ID) , data= plotThis[Celltype_subtype2=="Cancer cells"][Cohort=="Validation"][Treatmentlab!="Combination ribociclib"])
summary( MLv )
