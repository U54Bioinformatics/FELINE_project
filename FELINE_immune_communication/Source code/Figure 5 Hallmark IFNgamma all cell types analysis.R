rm(list=ls())
require(data.table);require(dplyr);require(tidyr)
require(ggplot2);require(ggsci)

require(lme4)
require(lmerTest)
require(merTools)

# read Hallmark IFNgamma ssGSEA scores of each cell of the discovery and validation cohort across TME cell types
savloc <- "/Users/jason/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/Figure5/Outputs/"
Hallmarkdd <- data.table(read.csv(file=paste0(savloc,"SourceData_Figure5_HallmarkIFNgammaTME_Output")))

# For each cell type and cohort, normalize to the average pathway activity in sensive tumors 
Hallmarkdd[,valuemuCelltypeSenstumor:= sum(value*(dynamic_class3=="Response"))/sum(dynamic_class3=="Response"),by=c("Celltype_subtype2","Cohort")]

# Order for plotting by effect size
Hallmarkdd$Celltype_subtype2 <-  factor(Hallmarkdd$Celltype_subtype2, levels=rev(c("Cancer cells", "Diploid epithelial cells",
                                                                               "Myeloid cells", "CD8+ T cells", "CD4+ T cells",
                                                                               "Adipocytes" , "Fibroblasts", "Endothelial cells","Pericytes"  #"NK cells",
) ))

Hallmarkdd[,Celltype_subtype3:=Celltype_subtype2]
#Hallmarkdd[Celltype_subtype2%in%c("CD8+ T cells","CD4+ T cells"  ),Celltype_subtype3:="T cells"]
Hallmarkdd[Celltype_subtype2%in%c("Adipocytes" , "Fibroblasts", "Endothelial cells","Pericytes"  ),Celltype_subtype3:="Stromal cells"]
Hallmarkdd$Celltype_subtype3 <-  factor(Hallmarkdd$Celltype_subtype3, levels=rev(c("Cancer cells", "Diploid epithelial cells",
                                                                               "Myeloid cells", "CD8+ T cells", "CD4+ T cells",
                                                                               "Stromal cells" 
) ))

#
HallmarkMean <- data.table( Hallmarkdd%>%group_by(Cohort,Celltype_subtype3,dynamic_class3,Patient.Study.ID,Day,Daylab,TumorResponse,Treatmentlab,ARM)%>%dplyr::summarise(value=mean(value),valuemuCelltypeSenstumor=mean(valuemuCelltypeSenstumor)))
ggplot( Hallmarkdd%>%group_by(Cohort,Celltype_subtype3,dynamic_class3,Patient.Study.ID,Day,Daylab,TumorResponse,Treatmentlab,ARM)%>%dplyr::summarise(value=mean(value),valuemuCelltypeSenstumor=mean(valuemuCelltypeSenstumor))
        , aes(y=value-valuemuCelltypeSenstumor,x=Celltype_subtype3, group=interaction(Celltype_subtype3,dynamic_class3),fill=dynamic_class3,col=dynamic_class3))+theme_classic(base_size=24)+
  # geom_violin(scale="width",bw=0.01,aes(group=interaction(Celltype,dynamic_class3))) +
  geom_boxplot(col="black",outlier.colour=NA, position= position_dodge() )+
  stat_boxplot(geom="errorbar",position=position_dodge())+#geom_smooth(method="lm")+
  scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  geom_point(aes(shape=Cohort,group=interaction(Celltype_subtype3,dynamic_class3)),col="black",position =position_dodge(width=1)  ,size=3)+
  theme(aspect.ratio=1)+
  theme(axis.text.x=element_text(angle=90,vjust=0.5))+
  facet_wrap(~Treatmentlab, ncol=1)+
  labs(x= "Cell type",y="Interferon gamma response \n (Hallmark ssGSEA score)")+ coord_flip()


ggplot( Hallmarkdd  , aes(y=value-valuemuCelltypeSenstumor,x=Celltype_subtype3, group=interaction(Celltype_subtype3,dynamic_class3),fill=dynamic_class3,col=dynamic_class3))+theme_classic(base_size=24)+
  # geom_violin(scale="width",bw=0.01,aes(group=interaction(Celltype,dynamic_class3))) +
  geom_boxplot(outlier.colour=NA, position= position_dodge(), alpha=0.6 )+
  stat_boxplot(geom="errorbar",position=position_dodge())+#geom_smooth(method="lm")+
  scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  geom_point(aes(shape=Cohort,group=interaction(Celltype_subtype3,dynamic_class3)),position =position_dodge(width=1)  ,size=1.75)+
  theme(aspect.ratio=1)+
  theme(axis.text.x=element_text(angle=90,vjust=0.5))+
  facet_wrap(~Treatmentlab, ncol=1)+
  labs(x= "Cell type",y="Interferon gamma response \n (Hallmark ssGSEA score)")+ coord_flip()
#save(Hallmarkdd,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/Figure5/Discover and Validation Hallmark IFNgamma TME.RData")


MC <- lmer(value-valuemuCelltypeSenstumor~ -1+ Celltype_subtype3+Celltype_subtype3 : dynamic_class3 + Cohort +(1|Patient.Study.ID) , data= Hallmarkdd[Treatmentlab=="Combination ribociclib"])
summary( MC )
ML <- lmer(value-valuemuCelltypeSenstumor~ -1+ Celltype_subtype3+Celltype_subtype3 : dynamic_class3 + Cohort +(1|Patient.Study.ID) , data= Hallmarkdd[Treatmentlab=="Letrozole alone"])
summary( ML )

MC <- lmer(value-valuemuCelltypeSenstumor~ dynamic_class3 + Cohort +(1|Patient.Study.ID) , data= Hallmarkdd[Celltype_subtype3=="Stromal cells"][Treatmentlab=="Combination ribociclib"])
summary( MC )

Combores<-rbindlist(lapply( unique(Hallmarkdd$Celltype_subtype3) , function(i){
  MC <- lmer(value-valuemuCelltypeSenstumor~ dynamic_class3 + Cohort +(1|Patient.Study.ID) , data= Hallmarkdd[Celltype_subtype3==i][Treatmentlab=="Combination ribociclib"])
  out<-data.table(Celltype_subtype3=i,coef(summary( MC )), keep.rownames=T)
  setnames(out, old=c( "Std. Error" , "t value","Pr(>|t|)" ), new=c( "Std.Error" , "tvalue","pvalue" ) )
  return(out)
}))
Letrozoleres<-rbindlist(lapply( unique(Hallmarkdd$Celltype_subtype3) , function(i){
  MC <- lmer(value-valuemuCelltypeSenstumor~ dynamic_class3 + Cohort +(1|Patient.Study.ID) , data= Hallmarkdd[Celltype_subtype3==i][Treatmentlab!="Combination ribociclib"])
  out<-data.table(Celltype_subtype3=i,coef(summary( MC )), keep.rownames=T)
  setnames(out, old=c( "Std. Error" , "t value","Pr(>|t|)" ), new=c( "Std.Error" , "tvalue","pvalue" ) )
  return(out)
}))

outall<-rbind( data.table(Treatmentlab="Combination ribociclib" ,Combores[rn=="dynamic_class3Response"]),
               data.table(Treatmentlab="Letrozole alone" ,Letrozoleres[rn=="dynamic_class3Response"]))

write.csv(outall,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/temporary/Discovery and Validation stats Hallmark IFN activationTME.csv")



MCv <- lmer(value-valuemuCelltypeSenstumor~  dynamic_class3  +(1|Patient.Study.ID) , data= Hallmarkdd[Celltype_subtype2=="Cancer cells"][Cohort=="Validation"][Treatmentlab=="Combination ribociclib"])
summary( MCv )



MLv <- lmer(value-valuemuCelltypeSenstumor~  dynamic_class3  +(1|Patient.Study.ID) , data= Hallmarkdd[Celltype_subtype2=="Cancer cells"][Cohort=="Validation"][Treatmentlab!="Combination ribociclib"])
summary( MLv )
