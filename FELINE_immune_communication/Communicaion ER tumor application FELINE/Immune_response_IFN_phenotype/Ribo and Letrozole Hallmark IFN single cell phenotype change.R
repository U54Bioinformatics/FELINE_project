rm(list=ls())
require(data.table);require(dplyr);require(tidyr)
require(ggplot2);require(ggsci)

require(lme4)
require(lmerTest)
require(merTools)


load("/Users/jason/Dropbox/FELINE Project (1)/Data_analysis/scRNA/01_metadata/results/AllCellTypeMetaData.RData")
hallmarksLoc<-"/Users/jason/Dropbox/FELINE Project (1)/Data_analysis/scRNA/05_ssGSEA_score/Signature_hallmarkAlone/results/"
hallmarknms <- list.files(hallmarksLoc)
dd0 <- rbindlist(lapply(1:length(hallmarknms),function(i){
  load(paste0(hallmarksLoc,list.files(hallmarksLoc)[i]))
  nms <- hallmarkdd$`Gene Set`
  out <- as.data.table( t(hallmarkdd[,-c("Gene Set","ssgseafile")]),
                        keep.rownames = T)
  colnames(out)<- c( "Cell.ID",nms)
  return(out)
}), fill = TRUE)

pathwaynms <- names(dd0)[-1]
dd1 <- merge(metadd[Platform=="10x"],dd0,by="Cell.ID")
rm(list="dd0")
dd1long <-data.table( gather(dd1,pathway,value,all_of(pathwaynms)) )

iii=27
p <- pathwaynms[iii]
cat(paste0("Pathway ", iii, ")     " , p ,"    "))
p_dd <- dd1long[pathway==p]
# Normalize SSGSEA relative to day 0 by celltype
p_dd[,valuemuDay0:= sum(value*(Day==0))/sum(Day==0),by=c("Patient.Study.ID","Celltype")]
p_dd[is.na(valuemuDay0),valuemuDay0:=0]
p_dd[,NormSSGSEA:=value-valuemuDay0]

ggplot(
  p_dd[][Day==180],
  aes(y=value,x=log(1+Day), fill=dynamic_class3,col=dynamic_class3))+
  geom_violin(aes(group=interaction(Day,dynamic_class3))) + facet_wrap(~Celltype) +geom_smooth(method="lm",se=F)+ facet_wrap(~ Treatment)

ggplot( p_dd[Treatment=="letrozole + ribo"][Day==180], aes(y=value,x=Celltype, fill=dynamic_class3,col=dynamic_class3))+theme_classic()+
 # geom_violin(scale="width",bw=0.01,aes(group=interaction(Celltype,dynamic_class3))) +
  geom_boxplot(col="black",aes(group=interaction(Celltype,dynamic_class3)))+
  geom_point()+
  theme(aspect.ratio=1)+
  ylab("HALLMARK_INTERFERON_GAMMA_RESPONSE")+theme(axis.text.x=element_text(angle=90))


plotThis<- p_dd[Day==180][Celltype!="B cells" ]
plotThis[,TumorResponse:="Resistant"]
plotThis[dynamic_class3=="Response",TumorResponse:="Sensitive"]

plotThis[,Treatmentlab:= "Combination ribociclib"]
plotThis[Treatment=="letrozole",Treatmentlab:= "Letrozole alone"]
plotThis$Celltype_subtype2 <- plotThis$Celltype_subtype 
plotThis[Celltype_subtype%in%c("Vas-Endo","Lym-Endo" ),Celltype_subtype2:= "Endothelial cells"]
plotThis[Celltype_subtype%in%c("Normal epithelial cells" ),Celltype_subtype2:= "Diploid epithelial cells"]
plotThis[Celltype_subtype%in%c("CD4+ T cells" ,"Tregs"),Celltype_subtype2:= "CD4+ T cells"]
plotThis[Celltype_subtype%in%c("CD8+ T cells" ,"NK cells"),Celltype_subtype2:= "CD8+ T cells"]
plotThis[Celltype_subtype%in%c("Monocytes" ,"DC","Macrophages"),Celltype_subtype2:= "Myeloid cells"]

plotThis[,Daylab:=paste0("Day ", Day)]
data.table(plotThis[dynamic_class3=="Response"][ARM!="A"]%>%group_by(Celltype_subtype2)%>%summarise(m=median(value)))[order(-m)]

plotThis$Celltype_subtype2 <-  factor(plotThis$Celltype_subtype2, levels=c("Cancer cells", "Adipocytes" ,"Diploid epithelial cells", "Fibroblasts",
                                                                           "CD8+ T cells", "Endothelial cells","Pericytes" , #"NK cells",
                                                                               "Myeloid cells", "CD4+ T cells"   ) )
#plotThis$Celltype_subtype2 <-  factor(plotThis$Celltype_subtype2, levels=c("Adipocytes" ,"Pericytes", "Fibroblasts","Endothelial cells","Diploid epithelial cells","Cancer cells" ,"Macrophages", "CD8+ T cells", "CD4+ T cells"  ) )
ggplot( plotThis, aes(y=value,x=Celltype_subtype2, fill=dynamic_class3,col=dynamic_class3))+theme_classic(base_size=22)+
 # geom_violin(scale="width",bw=0.01,aes(group=interaction(Celltype,dynamic_class3))) +
  geom_boxplot(outlier.colour=NA, position= position_dodge(), alpha=0.6 )+
  stat_boxplot(geom="errorbar",position=position_dodge(),width=0.5)+#geom_smooth(method="lm")+
  scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  geom_point(position =position_dodge(width=1)  ,size=0.5)+
  theme(aspect.ratio=1)+
  theme(axis.text.x=element_text(angle=90,vjust=0.5))+
  facet_wrap(~Treatmentlab, ncol=1)+
  labs(x= "Cell type",y="Interferon gamma response \n (Hallmark ssGSEA score)")+
  scale_x_discrete(labels=c("Cancer", "Adipocyte" ,"Diploid epithelial ", "Fibroblast",
                            "CD8+ T cell", "Endothelial","Pericyte" , #"NK cells",
                            "Myeloid", "CD4+ T cell") )
paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
paperfile<- "/Users/jason/Dropbox/FELINE Project (1)/Manuscript  Feline immune communication/Figures Communication Project/Paper figures Cancer immune communication/"

#ggsave(paste0(paperfile,"Ribo and Letrozole IFN gamma response all cell types.png"),height=10,width=10)


plotThis$Celltype_subtype2 <-  factor(plotThis$Celltype_subtype2, levels=rev(c("Cancer cells", "Adipocytes" ,"Diploid epithelial cells", "Fibroblasts",
                                                                           "CD8+ T cells", "Endothelial cells","Pericytes" , #"NK cells",
                                                                           "Myeloid cells", "CD4+ T cells"   ) ))
ggplot( plotThis, aes(y=value,x=Celltype_subtype2, fill=dynamic_class3,col=dynamic_class3))+theme_classic(base_size=24)+
  # geom_violin(scale="width",bw=0.01,aes(group=interaction(Celltype,dynamic_class3))) +
  geom_boxplot(outlier.colour=NA, position= position_dodge(), alpha=0.6 )+
  stat_boxplot(geom="errorbar",position=position_dodge(),width=0.5)+#geom_smooth(method="lm")+
  scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  geom_point(position =position_dodge(width=1)  ,size=0.5)+
  theme(aspect.ratio=1)+
  theme(axis.text.x=element_text(angle=90,vjust=0.5))+
  facet_wrap(~Treatmentlab, ncol=1)+
  labs(x= "Cell type",y="Interferon gamma response \n (Hallmark ssGSEA score)")+
  scale_x_discrete(labels=rev(c("Cancer", "Adipocyte" ,"Diploid epithelial ", "Fibroblast",
                            "CD8+ T cell", "Endothelial","Pericyte" , #"NK cells",
                            "Myeloid", "CD4+ T cell")))+ coord_flip()
ggsave(paste0(paperfile,"Ribo and Letrozole IFN gamma response all cell types FLIP.png"),height=10,width=10)



ggplot( plotThis, aes(y=value,x=Celltype_subtype2, fill=dynamic_class3,col=dynamic_class3))+theme_classic(base_size=21)+
  # geom_violin(scale="width",bw=0.01,aes(group=interaction(Celltype,dynamic_class3))) +
  geom_boxplot(outlier.colour=NA, position= position_dodge(), alpha=0.6 )+
  stat_boxplot(geom="errorbar",position=position_dodge(),width=0.5)+#geom_smooth(method="lm")+
  scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"),guide="none")+
  scale_fill_npg(name="Tumor response", labels=c("Resistant", "Sensitive"),guide="none")+
  #scale_x_discrete(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  #geom_boxplot(col="black",aes(group=interaction(Celltype,dynamic_class3)))+
  #geom_point(size=2.5)+
  geom_point(position =position_dodge(width=1)  ,size=0.5)+
  theme(aspect.ratio=1)+
  theme(axis.text.x=element_text(angle=90,vjust=0.5))+
  facet_wrap(Treatmentlab~Daylab, ncol=1)+
  labs(x= "Cell type",y="Interferon gamma response \n (Hallmark ssGSEA score)")+
  theme(axis.title=element_blank(),  axis.text=element_blank(),axis.text.x=element_blank(),strip.text = element_blank(), legend.position="none")

paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"BLANK Ribo and Letrozole IFN gamma response all cell types.png"),height=10,width=10)

MA <- lmer(value~ dynamic_class3* Treatment+ (1| Celltype_subtype2) + (1|Patient.Study.ID) , data= plotThis)
MB <- lmer(value~ dynamic_class3* Treatment + (1| Patient.Study.ID:Celltype_subtype2) + (1|Patient.Study.ID) , data= plotThis)
#MBii <- lmer(value~ TumorResponse* Treatment + (1| Patient.Study.ID:Celltype_subtype2) + (1|Patient.Study.ID) , data= plotThis)
summary( MA )
summary( MB )
#summary( MBii )

MC <- lmer(value~-1+ (Celltype_subtype2 + dynamic_class3+ Treatment )^3+(1|Patient.Study.ID) , data= plotThis)
summary( MC )


MC <- lmer(value~ -1+ Celltype_subtype2+(Celltype_subtype2 * dynamic_class3 )+(1|Patient.Study.ID) , data= plotThis[Treatment=="letrozole + ribo"])
summary( MC )


MD <- lmer(value~ 0+Celltype_subtype2 + (Celltype_subtype2 :dynamic_class3 )+(1|Patient.Study.ID) , data= plotThis[Treatment!="letrozole + ribo"])
summary( MD )
dd1longsumm <- data.table( dd1long%>%group_by(Celltype,Day,Patient.Study.ID,pathway,Treatment,dynamic_class3)%>%summarise(muval=mean(value),n=n(),sdval=sd(value)))

dd1widesumm_mean <- spread(dd1longsumm[pathway%in%c("HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_TNFA_SIGNALING_VIA_NFKB")]%>%dplyr::select(-c(n,sdval)) , pathway,muval )

# ggplot(
#   dd1longsumm[Treatment=="letrozole + ribo"][Day==180][pathway%in%c("HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_INFLAMMATORY_RESPONSE")],
#   aes(y=muval,x=Celltype, fill=dynamic_class3,col=dynamic_class3))+theme_classic()+
#   geom_violin(scale="width",aes(group=interaction(Celltype,dynamic_class3))) +geom_smooth(method="lm",se=F)  +theme(aspect.ratio=1)+facet_wrap(~pathway,scale="free_y")

ggplot(
  dd1widesumm_mean[Treatment=="letrozole + ribo"][Day==180],
  aes(y=HALLMARK_INTERFERON_GAMMA_RESPONSE,x=Celltype, fill=dynamic_class3,col=dynamic_class3))+theme_classic()+
  geom_violin(trim=F,scale="width",bw=0.0061,aes(group=interaction(Celltype,dynamic_class3))) +geom_smooth(method="lm",se=F)  +theme(aspect.ratio=1)

ggplot( dd1widesumm_mean[Treatment=="letrozole + ribo"][Day==180],
  aes(y=HALLMARK_TNFA_SIGNALING_VIA_NFKB,x=Celltype, fill=dynamic_class3,col=dynamic_class3))+theme_classic()+
  geom_violin(trim=F,scale="width",bw=0.05,aes(group=interaction(Celltype,dynamic_class3))) +geom_smooth(method="lm",se=F)  +theme(aspect.ratio=1)


#iii=27
for(iii in 1:length(pathwaynms)){
p <- pathwaynms[iii]
cat(paste0("Pathway ", iii, ")     " , p ,"    "))
p_dd <- dd1long[pathway==p]
# Normalize SSGSEA relative to day 0 by celltype
p_dd[,valuemuDay0:= sum(value*(Day==0))/sum(Day==0),by=c("Patient.Study.ID","Celltype")]
p_dd[is.na(valuemuDay0),valuemuDay0:=0]
p_dd[,NormSSGSEA:=value-valuemuDay0]

outCC <- lapply( unique(p_dd$Celltype),function(CC){   #CC<-"Endothelial cells"
  cat(CC)
  dd_i <- p_dd[Celltype==CC]
  if (!all(is.na(p_dd[Celltype==CC]$value))){#p_dd[Celltype==CC][!is.na(p_dd[Celltype==CC]$value)]
  m1  <- lmer(value ~ day_fact*dynamic_class3 + (1 + day_fact | Patient.Study.ID) , data=dd_i)
  #m1  <- lmer(value ~ Day*dynamic_class3 + (1 + day_fact | Patient.Study.ID) , data=dd_i)
  dd_i$preds <- predict(m1)
  dd_i$y <- predict(m1,re.form = NA)
  dd_i$y_i <- predict(m1,re.form =~(1+day_fact|Patient.Study.ID))
  newX <- data.table( unique(dd_i%>%dplyr::select(Day,day_fact,Celltype)) %>% 
    tidyr::crossing( dd_i%>%dplyr::select(Patient.Study.ID,ARM,dynamic_class3))  )
  newX$Patient.Study.ID <-as.factor(newX$Patient.Study.ID)
  newX$pred <-   predict(m1, newdata=newX,re.form = NA) 
  cinfreg <- data.table( predictInterval(m1,newX, which ="fixed", level = 0.95) )
  PRED.lme4 <-data.table( cbind(newX,cinfreg) %>% group_by(ARM,day_fact,dynamic_class3)%>%
                           dplyr::mutate(lwr=min(fit),upr=max(fit),fit=mean(fit) ))
  #ggplot(dd_i,aes(x=preds,y=value))+geom_point()+geom_abline()  #ggplot(dd_i,aes(x=preds,y=value-preds1))+geom_point()+geom_hline(yintercept=0)

  output <- data.frame(summary(m1)$ coefficients);  output$Effect <- rownames(output)
  output$Celltype <-dd_i[1]$Celltype   ;  output$pathway <- dd_i[1]$pathway
  output <- gather(output, statistic, value,"Estimate":"Pr...t..")
  output$metric <- "coefs"
  #anova info
  output2 <-   as.data.table(anova(m1)%>%dplyr::select(-c("NumDF","DenDF","Sum Sq","Mean Sq") ));  output2$Effect <- rownames(anova(m1))
  output2 <- na.omit(output2)
  output2$Celltype <- dd_i[1]$Celltype ;    output2$pathway <- dd_i[1]$pathway
  output2 <- gather(output2, statistic, value,"F value":"Pr(>F)")
  output2$metric <- "signifcov"
  # join info
  outputFull <- data.table(rbind(output,output2))
  return(list(info=outputFull,data=dd_i,uncert=PRED.lme4))
  }else{
    return(list(info=NULL,data=NULL,uncert=NULL))
  }
  
})
outlocFile <- "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE/FELINE2/LMER_SingleCell_SignifCelltypeTrendsByResponse"
save(outCC, file = paste0( outlocFile,p,"_LMER_SingleCell_", "SignifCelltypeTrendsByResponse.png"))

infer <- rbindlist(lapply(1:length(outCC),function(i){outCC[[i]][[1]]}))
infercoefs <- data.table(infer[metric=="coefs"] %>% spread(statistic,value))
infercoefs[, adjp:= p.adjust(Pr...t.., method = "fdr") ]
infercoefs[,SigClass:="NS"]; infercoefs[adjp<0.05,SigClass:="Signif"]

plot0 <- ggplot( infercoefs[grep("day",Effect,ignore.case=TRUE)],  aes(y=adjp,x=Estimate,col=SigClass))+geom_point()+
  labs(y="P-value (Holm's adjusted)", "Effect size")+theme_classic(base_size = 14)#+ geom_text(aes(label=Celltype),angle = 30,nudge_y = 0.05)

SS <- infercoefs[grep("day",Effect,ignore.case=TRUE)][adjp<0.05]

#unique(SS$Celltype)

fittedDD <- rbindlist(lapply(1:length(outCC),function(i){ outCC[[i]][[2]]}))
uncertDD <- rbindlist(lapply(1:length(outCC),function(i){ outCC[[i]][[3]]}))

pa <- ggplot(fittedDD[Celltype%in%unique(SS$Celltype)],aes_string(y="value",x="log(1+Day)",group="interaction(Day,dynamic_class3)"))+
  geom_violin(aes(fill=dynamic_class3),alpha=0.8)+#geom_jitter(aes(col=dynamic_class3),size=2,width=0.5,height=0)+
  theme_classic(base_size=20)+theme(aspect.ratio = 1)+
  facet_wrap(~Celltype,ncol=2)+
  scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))+
  labs(x="Day",y=gsub("_"," ",p ) )+
  theme(legend.position = "none",aspect.ratio=1)+
  geom_line(aes(y=preds,x=log(1+Day),col=dynamic_class3,group=Patient.Study.ID))


pb <- ggplot(fittedDD,
       aes_string(y="value",x="log(1+Day)",group="interaction(Day,dynamic_class3)"))+
  geom_violin(aes(fill=dynamic_class3),alpha=0.8)+#geom_jitter(aes(col=dynamic_class3),size=2,width=0.5,height=0)+
  theme_classic(base_size=20)+theme(aspect.ratio = 1)+
  facet_wrap(~Celltype)+
  scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))+
  labs(x="Day",y=gsub("_"," ",p ) )+
  theme(legend.position = "none",aspect.ratio=1)+
  geom_line(aes(y=preds,x=log(1+Day),col=dynamic_class3,group=Patient.Study.ID))

outloc<-"/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/Feline/SingleCellTrendsByResponseCelltype/"
ggsave(plot0, filename = paste0( outloc,"PvalPlot/",p,"_SignifEffectsizes_", "SingleCellTrends.png"))
if(nrow(SS)>0){ 
  ggsave(pa, filename = paste0( outloc,"Significant/",p,"_SignifCelltype_", "SingleCellTrendsByResponse.png")) 
}
ggsave(pb, filename = paste0( outloc,"All/",p,"_AllCelltype_", "SingleCellTrendsByResponse.png"),height=12,width=12)

}


# 
# 
# 
# 
# 
# 
# ggplot(fittedDD[Celltype%in%unique(SS$Celltype)],
#        aes_string(y="value",x="log(1+Day)",group="interaction(Day,dynamic_class3)"))+
#   geom_violin(aes(fill=dynamic_class3),alpha=0.8)+#geom_jitter(aes(col=dynamic_class3),size=2,width=0.5,height=0)+
#   theme_classic(base_size=20)+theme(aspect.ratio = 1)+
#   facet_wrap(~Celltype)+
#   scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))+
#   labs(x="Day",y=gsub("_"," ",p ) )+
#   theme(legend.position = "none",aspect.ratio=1)+
#   geom_line(aes(y=preds,x=log(1+Day),col=dynamic_class3,group=Patient.Study.ID))
# 
# ggplot(uncertDD[Celltype%in%unique(SS$Celltype)],
#        aes_string(y="fit",x="log(1+Day)",col="dynamic_class3"))+
#   theme_classic(base_size=20)+theme(aspect.ratio = 1)+
#   facet_wrap(~Celltype)+
#   scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))+
#   labs(x="Day",y=gsub("_"," ",p ) )+
#   theme(legend.position = "none",aspect.ratio=1)+
#   geom_line()+
#   geom_ribbon(aes(ymax=upr,ymin=lwr,x=log(1+Day),fill=dynamic_class3),alpha=0.6,col=NA)
# 
# ggplot(uncertDD[Celltype%in%unique(SS$Celltype)],
#        aes_string(y="fit",x="log(1+Day)",col="dynamic_class3"))+
#   theme_classic(base_size=20)+theme(aspect.ratio = 1)+
#   facet_wrap(~Celltype)+
#   scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))+
#   labs(x="Day",y=gsub("_"," ",p ) )+
#   theme(legend.position = "none",aspect.ratio=1)+
#   geom_line(data=fittedDD[Celltype%in%unique(SS$Celltype)],aes(y=preds,x=log(1+Day),col=dynamic_class3,group=Patient.Study.ID),size=0.5,linetype="dashed")+
#   geom_line()+
#   geom_ribbon(aes(ymax=upr,ymin=lwr,x=log(1+Day),fill=dynamic_class3),alpha=0.6,col=NA)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ggplot(p_dd,aes_string(y="value-valuemuDay0",x="log(1+Day)",group="Day"))+
#   geom_violin()+
#   #geom_jitter(aes(col=Celltype),size=2,width=0.15,height=0)+
#   theme_classic(base_size=20)+theme(aspect.ratio = 1)+
#   facet_wrap(~Celltype,ncol=3)+
#   scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))+
#   labs(x="Day",y=gsub("_"," ",p ) )+
#   theme(legend.position = "none",aspect.ratio=1)
# 
# ggplot(p_dd,aes_string(y="value-valuemuDay0",x="log(1+Day)",group="interaction(Day,dynamic_class3)"))+
#   geom_violin(aes(fill=dynamic_class3))+#geom_jitter(aes(col=dynamic_class3),size=2,width=0.5,height=0)+
#   theme_classic(base_size=20)+theme(aspect.ratio = 1)+
#   facet_wrap(~Celltype,ncol=3)+
#   scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))+
#   labs(x="Day",y=gsub("_"," ",p ) )+
#   theme(legend.position = "none",aspect.ratio=1)
# 
# 
# ggplot(p_dd,aes_string(y="value",x="log(1+Day)",group="interaction(Day,dynamic_class3)"))+
#   geom_violin(aes(fill=dynamic_class3))+#geom_jitter(aes(col=dynamic_class3),size=2,width=0.5,height=0)+
#   theme_classic(base_size=20)+theme(aspect.ratio = 1)+
#   facet_wrap(~Celltype,ncol=3)+
#   scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))+
#   labs(x="Day",y=gsub("_"," ",p ) )+
#   theme(legend.position = "none",aspect.ratio=1)



# 
# 
# m0 <- lmer(NormSSGSEA~ -1 + Day + Day:dynamic_class3 + (0+Day|Patient.Study.ID) , data=p_dd[Celltype=="Endothelial cells"])
# summary(m0)
# ranef(m0)
# 
# preds<- p_dd[Celltype=="Endothelial cells"]
# preds$preds<-predict(m0)
# preds$preds1<-predict(m1)
# 
# ggplot(preds,aes_string(y="value",x="log(1+Day)",group="interaction(Day,dynamic_class3)"))+
#   geom_violin(aes(fill=dynamic_class3),alpha=0.8)+#geom_jitter(aes(col=dynamic_class3),size=2,width=0.5,height=0)+
#   theme_classic(base_size=20)+theme(aspect.ratio = 1)+
#   facet_wrap(~Celltype,ncol=3)+
#   scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))+
#   labs(x="Day",y=gsub("_"," ",p ) )+
#   theme(legend.position = "none",aspect.ratio=1)+
#   geom_line(aes(y=preds1,x=log(1+Day),col=dynamic_class3,group=Patient.Study.ID))
# 
# 
# 
# #####
# 
# dd1longsumm <- data.table( dd1long%>%group_by(Celltype,Day,Patient.Study.ID,pathway,ARM,dynamic_class3,Celltype)%>%summarise(muval=mean(value),n=n(),sdval=sd(value)))
# 
# dd1widesumm_mean <- spread(dd1longsumm%>%dplyr::select(-c(n,sdval)) , pathway,muval )
# 
# x <- 27
# x <- 15
# 
# for(x in 1:length(pathwaynms)){
#   path_x <- pathwaynms[x]
#   pathxdd <- dd1widesumm_mean%>%dplyr::select(c(c("Celltype","Day","Patient.Study.ID","ARM","dynamic_class3","Celltype"),all_of(path_x) ))
#   pathxdd$pathname <- path_x
#   setnames(pathxdd, old=path_x, new="pathwayScore")
#   
#   normcelltypesnms <- unique(pathxdd[Celltype!="Cancer cells"]$Celltype)
#   pathxddwide <- data.table( gather( spread( pathxdd , Celltype,pathwayScore),
#                                      NormalCelltype, pathwayScore,all_of(normcelltypesnms)) )
#   setnames(pathxddwide, old="Cancer cells", new="CancerpathwayActivity")
#   
#   #ggplot(pathxdd,aes(y=HALLMARK_INTERFERON_GAMMA_RESPONSE,x=Patient.Study.ID))+
#   #  geom_point()
#   
#   p1<- ggplot(pathxddwide,aes(x=CancerpathwayActivity,y=pathwayScore,col=NormalCelltype,fill=NormalCelltype))+
#     geom_smooth(method="lm")+geom_point(size=1.5)+ theme_classic(base_size = 21)+facet_wrap(~NormalCelltype,scales="free",ncol=4)+
#     labs(x="Cancer pathway Activity",y="Focal cell type \n pathway Activity") +
#     theme(legend.position = "none",aspect.ratio = 1)
#   
#   outpath <- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/Feline/Correlated hallmark activity across cell types/"
#   ggsave( p1 , filename=paste0(outpath,path_x,"correlation_cancer_vs_normalcells.png"), height=7,width=14)
# }
# 
# # results in zip file : CorrelationOfHallmarkExpressionAcrossCellTypes
# 
# traits <- data.table(t(na.omit(t(dd1widesumm_mean[,pathwaynms,with=F] ))))
# umod <- umap::umap( traits,n_components=2 )
# 
# y2 <- cbind( dd1widesumm_mean,umod$layout)
# ggplot(y2, aes(V1,V2, col= Celltype))+geom_point() 
# ggplot(y2, aes(V1,V2, col= dynamic_class3))+geom_point() 
# ggplot(y2, aes(V1,V2, col= HALLMARK_GLYCOLYSIS))+geom_point() 
# 
# traits <- data.table((na.omit(t(dd1widesumm_mean[,pathwaynms,with=F] ))))
# umod <- umap::umap( traits,n_components=2 )
# plot(umod$layout)
# 
# 
# ggplot(dd1widesumm_mean,aes(y=HALLMARK_TNFA_SIGNALING_VIA_NFKB,x=dynamic_class3))+geom_jitter()+ facet_wrap(~Day)
# ggplot(dd1widesumm_mean,aes(y=HALLMARK_IL6_JAK_STAT3_SIGNALING,x=dynamic_class3))+geom_jitter()+ facet_wrap(~Day)
# ggplot(dd1widesumm_mean,aes(y=HALLMARK_INFLAMMATORY_RESPONSE,x=dynamic_class3))+geom_jitter()+ facet_wrap(~Day)
# ggplot(dd1widesumm_mean,aes(y=HALLMARK_IL2_STAT5_SIGNALING,x=dynamic_class3))+geom_jitter()+ facet_wrap(~Day)
# ggplot(dd1widesumm_mean,aes(y=HALLMARK_HYPOXIA,x=dynamic_class3))+geom_violin()+ facet_wrap(~Day)+theme_classic()+theme(aspect.ratio = 1)
# 
# 
# 
# ggplot(dd1widesumm_mean,aes(y=HALLMARK_HYPOXIA,
#                             x=log(1+Day),col=Celltype,group=interaction(Celltype,Patient.Study.ID)))+geom_line()+ 
#   facet_wrap(~dynamic_class3)
# 
# ggplot(dd1widesumm_mean,aes(y=HALLMARK_HYPOXIA,x=log(1+Day),group=Day))+geom_violin()+theme_classic()+theme(aspect.ratio = 1)
# ggplot(dd1widesumm_mean,aes(y=HALLMARK_HYPOXIA,x=log(1+Day),group=Day))+geom_violin()+geom_jitter(aes(col=Celltype),size=2,width=0.5,height=0)+
#   theme_classic()+theme(aspect.ratio = 1)+facet_wrap(~Celltype,ncol=3)
# 
# y<-pathwaynms[1]
# 
# f <- function( p ) {
#   
#   ggplot(dd1widesumm_mean,aes_string(y=p,x="log(1+Day)",group="Day"))+geom_violin()+geom_jitter(aes(col=Celltype),size=2,width=0.5,height=0)+
#     theme_classic(base_size=20)+theme(aspect.ratio = 1)+
#     facet_wrap(~Celltype,ncol=3)+
#     scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))+
#     labs(x="Day",y=gsub("_"," ",p ) )+
#     theme(legend.position = "none",aspect.ratio=1)
# }  
# for(i in 1: length(pathwaynms)){
#   pp1<- f (pathwaynms[i]) 
#   ggsave(pp1,filename = paste0("/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/Feline/Hallmarks vs time by cell type/",
#                                pathwaynms[i],"overtime_percelltype.png"),height=10,width=10)       
#   
# }
# 
# f1 <- function( p ) {
#   
#   ggplot(dd1widesumm_mean,aes_string(y=p,x="log(1+Day)",group="interaction(Day,dynamic_class3)"))+geom_violin(aes(fill=dynamic_class3))+
#     geom_point(aes(fill=dynamic_class3),colour="black",pch=21,position=position_jitterdodge(dodge.width=2,jitter.width=0.4) ,size=2,width=0.5,height=0)+
#     
#     #geom_jitter(aes(col=dynamic_class3),size=2,width=0.5,height=0)+
#     theme_classic(base_size=20)+theme(aspect.ratio = 1)+
#     facet_wrap(~Celltype,ncol=3)+
#     scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))+
#     labs(x="Day",y=gsub("_"," ",p ) )+
#     theme(legend.position = "none",aspect.ratio=1)
# }  
# 
# f2 <- function( p ) {
#   
#   ggplot(dd1widesumm_mean,aes_string(y=p,x="log(1+Day)",group="interaction(Day,dynamic_class3)"))+geom_violin(aes(fill=dynamic_class3))+geom_jitter(aes(col=dynamic_class3),size=2,width=0.5,height=0)+
#     theme_classic(base_size=20)+theme(aspect.ratio = 1)+
#     #facet_wrap(~Celltype,ncol=3)+
#     scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))+
#     labs(x="Day",y=gsub("_"," ",p ) )+
#     theme(legend.position = "none",aspect.ratio=1)
# } 
# 
# 
# for(i in 1: length(pathwaynms)){
#   pp1<- f1 (pathwaynms[i]) 
#   ggsave(pp1,filename = paste0("/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/Feline/Hallmark vs time by response/",
#                                pathwaynms[i],"overtime_byresponse.png"),height=10,width=10)       
#   cat(i)
# }
# 
# 
# 
# f1 (pathwaynms[24]) 
# f1 (pathwaynms[44]) 
# f1 (pathwaynms[45]) 
# 
# f1 (pathwaynms[27]) 
# 
# f2 (pathwaynms[24]) 
# f2 (pathwaynms[45]) 
# 
# f2 (pathwaynms[25]) 
# f2 (pathwaynms[22]) 
# f2 (pathwaynms[27]) 
# f2 (pathwaynms[46]) 
# f2 (pathwaynms[50]) 
# 
# 
# i<- 44 #45 #27 #
# pp1<-ggplot(dd1widesumm_mean[Celltype%in%c("B cells","Macrophages", "T cells")],
#             aes_string(y=pathwaynms[i],x="log(1+Day)",group="interaction(Day,dynamic_class3)"))+
#   geom_violin(aes(fill=dynamic_class3))+
#   geom_point(aes(fill=dynamic_class3),colour="black",pch=21,position=position_jitterdodge(dodge.width=2,jitter.width=0.4) ,size=2,width=0.5,height=0)+
#   theme_classic(base_size=20)+theme(aspect.ratio = 1)+
#   facet_wrap(~Celltype,ncol=3)+
#   scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))+
#   labs(x="Day",y=gsub("_"," ",pathwaynms[i] ) )+
#   theme(aspect.ratio=1)+scale_fill_npg(name="Tumor response")
# ggsave(pp1,filename = paste0("/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/Feline/Hallmark vs time by response/",
#                              pathwaynms[i],"Immune_overtime_byresponse.png"),height=4,width=10)       
# 
# 
# i<- 45 #27 #
# pp1<-ggplot(dd1widesumm_mean[Celltype%in%c("B cells","Macrophages", "T cells")],
#             aes_string(y=pathwaynms[i],x="log(1+Day)",group="interaction(Day,dynamic_class3)"))+
#   geom_violin(aes(fill=dynamic_class3))+
#   geom_point(aes(fill=dynamic_class3),colour="black",pch=21,position=position_jitterdodge(dodge.width=2,jitter.width=0.4) ,size=2,width=0.5,height=0)+
#   theme_classic(base_size=20)+theme(aspect.ratio = 1)+
#   facet_wrap(~Celltype,ncol=3)+
#   scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))+
#   labs(x="Day",y=gsub("_"," ",pathwaynms[i] ) )+
#   theme(legend.position = "none",aspect.ratio=1)
# ggsave(pp1,filename = paste0("/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/Feline/Hallmark vs time by response/",
#                              pathwaynms[i],"Immune_overtime_byresponse.png"),height=4,width=10)       
# 
# 
# i<- 27
# pp1<-ggplot(dd1widesumm_mean[Celltype%in%c("B cells","Macrophages", "T cells")],
#             aes_string(y=pathwaynms[i],x="log(1+Day)",group="interaction(Day,dynamic_class3)"))+
#   geom_violin(aes(fill=dynamic_class3))+
#   geom_point(aes(fill=dynamic_class3),colour="black",pch=21,position=position_jitterdodge(dodge.width=2,jitter.width=0.4) ,size=2,width=0.5,height=0)+
#   theme_classic(base_size=20)+theme(aspect.ratio = 1)+
#   facet_wrap(~Celltype,ncol=3)+
#   scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))+
#   labs(x="Day",y=gsub("_"," ",pathwaynms[i] ) )+
#   theme(legend.position = "none",aspect.ratio=1)
# ggsave(pp1,filename = paste0("/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/Feline/Hallmark vs time by response/",
#                              pathwaynms[i],"Immune_overtime_byresponse.png"),height=4,width=10)       
# 
# 
# i<- 22# 23
# pp1<- ggplot(dd1widesumm_mean[Celltype%in%c("B cells","Macrophages", "T cells")],
#              aes_string(y=pathwaynms[i],x="log(1+Day)",group="interaction(Day,dynamic_class3)"))+
#   geom_violin(aes(fill=dynamic_class3))+
#   geom_point(aes(fill=dynamic_class3),colour="black",pch=21,position=position_jitterdodge(dodge.width=2,jitter.width=0.4) ,size=2,width=0.5,height=0)+
#   theme_classic(base_size=20)+theme(aspect.ratio = 1)+
#   facet_wrap(~Celltype,ncol=3)+
#   scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))+
#   labs(x="Day",y=gsub("_"," ",pathwaynms[i] ) )+
#   theme(legend.position = "none",aspect.ratio=1)
# ggsave(pp1,filename = paste0("/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/Feline/Hallmark vs time by response/",
#                              pathwaynms[i],"Immune_overtime_byresponse.png"),height=4,width=10)       
# 
# 
# 
# 
# 
# ggplot(dd1widesumm_mean,aes_string(y=p,x="log(1+Day)",group=Day))+geom_violin()+geom_jitter(aes(col=Celltype),size=2,width=0.5,height=0)+
#   theme_classic()+theme(aspect.ratio = 1)+facet_wrap(~Celltype,ncol=3)+
#   scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))+labs(x="Day")+
#   theme(legend.position = "none")
# 
# ggplot( rates.by.groups, aes_string(x="name", y="rate", colour= column,
#                                     group=column ) )
# }
# ggplot(dd1widesumm_mean,aes(y=HALLMARK_UV_RESPONSE_UP,x=log(1+Day),group=Day))+geom_violin()+geom_jitter(aes(col=Celltype),size=2,width=0.5,height=0)+
#   theme_classic()+theme(aspect.ratio = 1)+facet_wrap(~Celltype,ncol=3)+
#   scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))+labs(x="Day")+
#   theme(legend.position = "none")
# 
# 
# 
# 
# 
# ggplot(dd1widesumm_mean,aes(y=HALLMARK_INFLAMMATORY_RESPONSE,x=log(1+Day),group=Day))+geom_violin()+geom_jitter(aes(col=Celltype),size=2,width=0.5,height=0)+
#   theme_classic()+theme(aspect.ratio = 1)+facet_wrap(~Celltype,ncol=3)+
#   scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))+labs(x="Day")
# 
# 
# ggplot(dd1widesumm_mean,aes(y=HALLMARK_TNFA_SIGNALING_VIA_NFKB,x=log(1+Day),group=Day))+geom_violin()+geom_jitter(aes(col=Celltype),size=2,width=0.5,height=0)+
#   theme_classic()+theme(aspect.ratio = 1)+facet_wrap(~Celltype,ncol=3)+
#   scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))+labs(x="Day")
# 
# 
# ggplot(dd1widesumm_mean,aes(y=HALLMARK_ESTROGEN_RESPONSE_LATE,x=log(1+Day),group=Day))+geom_violin()+geom_jitter(aes(col=Celltype),size=2,width=0.5,height=0)+
#   theme_classic()+theme(aspect.ratio = 1)+facet_wrap(~Celltype,ncol=3)+
#   scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))+labs(x="Day")
# 
# 
# ggplot(dd1widesumm_mean,aes(y=HALLMARK_UV_RESPONSE_UP,x=log(1+Day),group=Day))+geom_violin()+geom_jitter(aes(col=Celltype),size=2,width=0.5,height=0)+
#   theme_classic()+theme(aspect.ratio = 1)+facet_wrap(~Celltype,ncol=3)+
#   scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))+labs(x="Day")+
#   theme(legend.position = "none")
# 
# 
# ggplot(dd1widesumm_mean,aes(y=HALLMARK_ESTROGEN_RESPONSE_LATE,x=log(1+Day),group=Day))+geom_violin()+geom_jitter(aes(col=Celltype),size=2,width=0.5,height=0)+
#   theme_classic()+theme(aspect.ratio = 1)+facet_wrap(~Celltype,ncol=3)+
#   scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))+labs(x="Day")
# 
# 
# 
# 
# 
# ggplot(dd1widesumm_mean,aes(y=HALLMARK_ESTROGEN_RESPONSE_LATE,x=dynamic_class3))+geom_jitter()+ facet_wrap(~Day)
# ggplot(dd1widesumm_mean,aes(y=HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,x=dynamic_class3))+geom_jitter()+ facet_wrap(~Day)+
#   ggplot(dd1widesumm_mean,aes(y=HALLMARK_G2M_CHECKPOINT,x=dynamic_class3))+geom_jitter()+ facet_wrap(~Day)
# ggplot(dd1widesumm_mean,aes(y=HALLMARK_GLYCOLYSIS,x=dynamic_class3))+geom_jitter()+ facet_wrap(~Day)
# 


