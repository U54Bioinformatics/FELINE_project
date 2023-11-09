rm(list=ls())
require(mgcv);require(data.table);require(dplyr);require(ggplot2);require(tidyr);require(lme4);require(lmerTest);require(parallel)
library(effects);require(umap);require(Rfast);require(ider);library("dendextend");library(ggdendro);require(ggsci);require(viridis);require("Rdimtools")
library(caret)
library(pROC)

load( file=paste0("~/Dropbox/Cancer_pheno_evo/data/FELINE2/PhenotypesAllArmsCohort2/CPMPhenotpyeLandscape_C2ProjectC2RevisednewMacrophagesDC.RData"))#C1umap,u_dat,DAY,cell_types_all,ARMS,Subtype,
load( file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Cohort2Metadata/Cohort2Metadata.RData") #metadd,cohort1metadd,annotation.file,compdataLU

dd1 <- u_dat %>%dplyr::select(Cell.ID:file_string,V1,V2)

Cohort1Extract <-function(){
  load(file= "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/UMAP genes macrophage only ALLARMS/UMAP genes macrophage only ALLARMS.RData")
  dd1<-u_dat%>%dplyr::select(Cell.ID:PhenoCelltype,V1,V2)
  return(dd1)
}
ddCohort1 <- Cohort1Extract()

ddCohort2 <- dd1
ddCohort2[,Treatment:="letrozole + ribo"]
ddCohort2[ARM=="A",Treatment:="letrozole"]


joineddd <- rbind(
  data.table(Cohort="Discovery" ,    ddCohort1%>%select(Treatment,Cell.ID,Sample,Patient.Study.ID,Day,ARM,dynamic_class3,Celltype,V1,V2) ),
  data.table(Cohort="Validation" ,  ddCohort2%>%select(Treatment,Cell.ID,Sample,Patient.Study.ID,Day,ARM,dynamic_class3,Celltype,V1,V2))
)

summarizedCohorts<- data.table(joineddd[Day!=180] %>% group_by(Treatment,#Day,
                                                   Patient.Study.ID,dynamic_class3,Cohort,ARM,Celltype) %>%
  summarise(n= n(),V1= mean(V1), V2= mean(V2) ,n= n() ) %>%
    group_by(Cohort) %>% mutate(cohortmu= mean(V2))%>%
    group_by(Treatment,Cohort)%>%mutate(cohortTreatmu= median(V2)))

summarizedCohorts[,Success:=1]
summarizedCohorts[dynamic_class3=="Non-response",Success:=0]

ggplot(summarizedCohorts     , aes(y=V2-cohortmu, x=Cohort,col=dynamic_class3, group=interaction(Cohort,dynamic_class3) )) +
  geom_violin(scale="width")+
  geom_point(aes(size=n,shape=Treatment),positio=position_dodge(width=0.9))+
  theme_classic()+theme(aspect.ratio = 1)+facet_wrap(~Treatment)


ggplot(summarizedCohorts     , aes(y=V2-cohortmu, x=Cohort,col=dynamic_class3, group=interaction(Cohort,dynamic_class3) )) +
  geom_boxplot(scale="width")+
  geom_point(aes(size=n,shape=Treatment),positio=position_dodge(width=0.9))+
  theme_classic()+theme(aspect.ratio = 1)+facet_wrap(~Treatment)

summarizedCohorts[,TumorDifferScore:=V2-cohortmu ]

ggplot(summarizedCohorts, 
       aes(y=V2-cohortmu, x=dynamic_class3,col=Cohort, group=interaction(Cohort,dynamic_class3) )) +
  geom_violin(scale="width")+
  geom_point(aes(shape=Treatment),positio=position_dodge(width=0.9))+
  theme_classic()+theme(aspect.ratio = 1)+facet_wrap(~Treatment)

summarizedCohorts[,Treatmentlab:= "Combination ribociclib"]
summarizedCohorts[Treatment=="letrozole",Treatmentlab:= "Letrozole alone"]

ggplot(summarizedCohorts,  aes(y=TumorDifferScore, x=dynamic_class3,col=dynamic_class3, fill=dynamic_class3, shape=Cohort, group=interaction(dynamic_class3) )) +
  geom_violin(scale="width", col="NA")+
  geom_point(col="black",size=2.5)+
  scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"),guide="none")+
  scale_fill_npg(name="Tumor response", labels=c("Resistant", "Sensitive"),guide="none")+
  scale_x_discrete(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  theme_classic(base_size=16)+theme(aspect.ratio = 1)+facet_wrap(~Treatmentlab,ncol=1)+
  labs(y="Myeloid M2 differentiation \n (tumor mean)")+theme(legend.position="bottom")
ggsave(filename=paste0("/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cohort 2 results/Myeloid differentiation predicts response to ribo not letrozole.png"),width=5,height=5)

ggplot(summarizedCohorts,  aes(y=TumorDifferScore, x=dynamic_class3,col=dynamic_class3, fill=dynamic_class3, shape=Cohort, group=interaction(dynamic_class3) )) +
 # geom_boxplot(scale="width")+
 # stat_boxplot(scale="errorbar",linetype=1,width=2,col="black")+
  
  geom_boxplot(outlier.colour=NA, position= position_dodge() ,col="black")+
  stat_boxplot(geom="errorbar",position=position_dodge(1.75),width=0.5)+#geom_smooth(method="lm")+
  
  geom_point(col="black",size=2.5)+
  scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"),guide="none")+
  scale_fill_npg(name="Tumor response", labels=c("Resistant", "Sensitive"),guide="none")+
  scale_x_discrete(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  theme_classic(base_size=26)+theme(aspect.ratio = 1)+facet_wrap(~Treatmentlab,ncol=1)+
  labs(y="Myeloid M2 differentiation \n (tumor mean)")+theme(legend.position="bottom")
ggsave(filename=paste0("/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cohort 2 results/Boxplot Myeloid differentiation predicts response to ribo not letrozole.png"),width=5,height=5)
paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Ribo and Letrozole Boxplot Myeloid differentiation predicts response to ribo not letrozole.png"),height=10,width=10)


summary(lm(TumorDifferScore~dynamic_class3*Treatmentlab,summarizedCohorts))

ggplot(summarizedCohorts,  aes(y=TumorDifferScore, x=dynamic_class3,col=dynamic_class3, fill=dynamic_class3, shape=Cohort, group=interaction(dynamic_class3) )) +
  geom_boxplot(scale="width",alpha=0.6)+
  
  geom_boxplot(outlier.colour=NA, position= position_dodge() ,col="black")+
  stat_boxplot(geom="errorbar",position=position_dodge(1.75),width=0.5)+#geom_smooth(method="lm")+
  
  geom_point(size=2.5,col="black")+
  scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"),guide="none")+
  scale_fill_npg(name="Tumor response", labels=c("Resistant", "Sensitive"),guide="none")+
  scale_x_discrete(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  theme_classic(base_size=16) + theme(aspect.ratio = 1)+facet_wrap(~Treatmentlab,ncol=1)+
  labs(y="Myeloid M2 differentiation \n (tumor mean)")+theme(legend.position="bottom") +
  theme(axis.title=element_blank(),  axis.text=element_blank(),strip.text = element_blank(),legend.title=element_blank(),  legend.text=element_blank())
#ggsave(paste0(paperfile,"BLANK Ribo and Letrozole Boxplot Myeloid differentiation predicts response to ribo not letrozole.png"),height=10,width=10)



lm(TumorDifferScore~dynamic_class3*Treatmentlab, summarizedCohorts)%>%summary()
lm(TumorDifferScore~dynamic_class3, summarizedCohorts[Treatmentlab=="Letrozole alone"])%>%summary()
lm(TumorDifferScore~dynamic_class3, summarizedCohorts[Treatmentlab!="Letrozole alone"])%>%summary()

#RJ_glm <- glm(Success~V2*Treatment+Cohort*Treatment, summarizedCohorts,family="binomial")
RJ_glm <- glm(Success~I(V2-cohortmu)*Treatment , summarizedCohorts,family="binomial")
RJ_glm <- mgcv::gam(Success ~ s(TumorDifferScore, k=5, by=as.factor(Treatment)) , summarizedCohorts,family="binomial")
summary(RJ_glm)
#RJ_glm <- glm(Success~I(V2-cohortTreatmu)*Treatment*Cohort , summarizedCohorts,family="binomial")
summarizedCohorts$glmPred <- predict(RJ_glm,type="response")

ggplot(summarizedCohorts,aes(Success, x= TumorDifferScore,col=Treatment)) + geom_point() +
  geom_line(aes(y= glmPred, group=Cohort))+facet_wrap(~Treatment)


ggplot(summarizedCohorts,aes(y=Success, x= TumorDifferScore,col=Treatment)) + geom_point() +
  geom_smooth(method="gam", method.args=list(family="binomial"))+
  geom_line(aes(y= glmPred, group=Cohort))+facet_wrap(~Treatment)



summarizedCohorts[,MyeloidStatus:= "M1 polarized"]
summarizedCohorts[(V2)>0,MyeloidStatus:= "M2 polarized"]
#summarizedCohorts[,Pred:= "Response"]
#summarizedCohorts[(V2)>0,Pred:= "Non-response"]
summarizedCohorts[,Pred:= "Response"]
thresh <- max(summarizedCohorts[Treatment!="letrozole"][glmPred>0.5]$TumorDifferScore)
summarizedCohorts[(glmPred)<=0.5,Pred:= "Non-response"]


tableL1 <- data.frame(confusionMatrix(as.factor(summarizedCohorts[Treatment=="letrozole"&Cohort=="Discovery"]$Pred), as.factor(summarizedCohorts[Treatment=="letrozole"&Cohort=="Discovery"]$dynamic_class3))$table)%>% group_by(Reference) %>%
  mutate(prop = Freq/sum(Freq))
tableL2 <- data.frame(confusionMatrix(as.factor(summarizedCohorts[Treatment=="letrozole"&Cohort=="Validation"]$Pred), as.factor(summarizedCohorts[Treatment=="letrozole"&Cohort=="Validation"]$dynamic_class3))$table)%>% group_by(Reference) %>%
  mutate(prop = Freq/sum(Freq))
tableR1 <- data.frame(confusionMatrix(as.factor(summarizedCohorts[Treatment!="letrozole"&Cohort=="Discovery"]$Pred), as.factor(summarizedCohorts[Treatment!="letrozole"&Cohort=="Discovery"]$dynamic_class3))$table)%>% group_by(Reference) %>%
  mutate(prop = Freq/sum(Freq))
tableR2 <- data.frame(confusionMatrix(as.factor(summarizedCohorts[Treatment!="letrozole"&Cohort=="Validation"]$Pred), as.factor(summarizedCohorts[Treatment!="letrozole"&Cohort=="Validation"]$dynamic_class3))$table)%>% group_by(Reference) %>%
  mutate(prop = Freq/sum(Freq))

tableR12 <- data.table( data.frame(confusionMatrix(as.factor(summarizedCohorts[Treatment!="letrozole"]$Pred), as.factor(summarizedCohorts[Treatment!="letrozole"]$dynamic_class3))$table)%>% group_by(Reference) %>%
  mutate(prop = Freq/sum(Freq)) )
tableL12 <- data.table( data.frame(confusionMatrix(as.factor(summarizedCohorts[Treatment=="letrozole"]$Pred), as.factor(summarizedCohorts[Treatment=="letrozole"]$dynamic_class3))$table)%>% group_by(Reference) %>%
  mutate(prop = Freq/sum(Freq)) )

tableR12[,TumorResponse:="Resistant"]
tableR12[Reference=="Response",TumorResponse:="Sensitive"]
tableL12[,TumorResponse:="Resistant"]
tableL12[Reference=="Response",TumorResponse:="Sensitive"]

tableR12[,MyeloidStatus:="M2"]
tableR12[Prediction=="Response",MyeloidStatus:="M1"]
tableL12[,MyeloidStatus:="M1"]
tableL12[Prediction=="Response",MyeloidStatus:="M2"]
tableL12$TumorResponse <- factor(tableL12$TumorResponse  , levels=c("Sensitive","Resistant"))
tableR12$TumorResponse <- factor(tableR12$TumorResponse  , levels=c("Sensitive","Resistant"))

jointtable<- rbind( data.table(Treatment="Combination ribociclib",tableR12),
   data.table(Treatment="Letrozole alone",tableL12) )
ggplot(data = jointtable, mapping = aes(y = TumorResponse, x = MyeloidStatus,  fill = prop)) +
  geom_tile(alpha=0.9) +  scale_fill_viridis_c(name="Proportion tumors \n with myeloid cells \n in each state",option="D")+ labs(y="Tumor response", x="Myeloid differentiation")+
  geom_text(aes(label = Freq), vjust = .5, fontface  = "bold", alpha = 1,col="black") +
  theme_classic() +facet_wrap(~Treatment)+theme(aspect.ratio=1)

ggsave(filename=paste0("/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cohort 2 results/Myeloid differentiation prediction accuracy to ribo not letrozole.png"),width=6,height=6)


# Myeloid differentiation prediction of response letrozole and ribociclib treatment
Rchisq<-chisq.test(summarizedCohorts[Treatment!="letrozole"]$MyeloidStatus,summarizedCohorts[Treatment!="letrozole"]$dynamic_class3)
Lchisq<-chisq.test(summarizedCohorts[Treatment=="letrozole"]$MyeloidStatus,summarizedCohorts[Treatment=="letrozole"]$dynamic_class3)
Rchisq
Lchisq
library(corrplot)
corrplot(Rchisq$residuals, is.cor=F,cl.lim =range(c(Rchisq$residuals,Lchisq$residuals)) )
corrplot(Lchisq$residuals, is.cor=F,cl.lim =range(c(Rchisq$residuals,Lchisq$residuals)) )




rocFun<-function(treat,cohort){
  rocobj <- roc(
    as.numeric(summarizedCohorts[Treatment==treat&Cohort==cohort]$Pred=="Response"),
    as.numeric(summarizedCohorts[Treatment==treat&Cohort==cohort]$dynamic_class3=="Response"),
    
  )
  
  #create ROC plot
  ggroc(rocobj)
  return(rocobj)
}

R1 <- ggroc( rocFun("letrozole + ribo","Discovery"),colour = 'steelblue' )+theme_classic()
R2 <- ggroc( rocFun("letrozole + ribo","Validation"),colour = 'steelblue' )+theme_classic()

L1 <- ggroc( rocFun("letrozole","Discovery"),colour = 'steelblue' )+theme_classic()
L2 <- ggroc( rocFun("letrozole","Validation"),colour = 'steelblue' )+theme_classic()
Allrocs<-rbind( data.table(Treatment="letrozole + ribo",cohort="Discovery",R1$data,auc=rocFun("letrozole + ribo","Discovery")$auc), 
       data.table(Treatment="letrozole + ribo",cohort="Validation",R2$data,auc=rocFun("letrozole + ribo","Validation")$auc),
       data.table(Treatment="letrozole",cohort="Discovery",L1$data,auc=rocFun("letrozole","Discovery")$auc), 
       data.table(Treatment="letrozole",cohort="Validation",L2$data,auc=rocFun("letrozole","Validation")$auc)
)

ggplot( Allrocs  , aes(x=1-specificity, y=sensitivity, col=Treatment, linetype=cohort,group=interaction(Treatment,cohort) ) ) + geom_path(size=1.5) +theme_classic()+
  labs(y="True positive rate (sensitivity)",x="False positive rate (1-specificity)")+
  scale_color_jco(name="Treatment",labels=c("Letrozole alone", "Combination ribociclib"))+theme(aspect.ratio=1.1)+
  scale_linetype(name="Cohort")





# 
# 
# plotTable <- tableL1 %>%
#   mutate(goodbad = ifelse(tableL1$Prediction == tableL1$Reference, "Non-response", "Response")) %>%
#   
# 
# # fill alpha relative to sensitivity/specificity by proportional outcomes within reference groups (see dplyr code above as well as original confusion matrix for comparison)
# ggplot(data = tableR1, mapping = aes(x = Reference, y = Prediction,  fill = prop)) +
#   geom_tile() +  scale_fill_viridis_c()+
#   geom_text(aes(label = Freq), vjust = .5, fontface  = "bold", alpha = 1,col="slategrey") +
#   theme_minimal() 
# ggplot(data = tableL1, mapping = aes(x = Reference, y = Prediction,  fill = prop)) +
#   geom_tile() +  scale_fill_viridis_c()+
#   geom_text(aes(label = Freq), vjust = .5, fontface  = "bold", alpha = 1,col="slategrey") +
#   theme_minimal() 
# ggplot(data = tableR2, mapping = aes(x = Reference, y = Prediction,  fill = prop)) +
#   geom_tile() +  scale_fill_viridis_c()+
#   geom_text(aes(label = Freq), vjust = .5, fontface  = "bold", alpha = 1,col="slategrey") +
#   theme_minimal() 
# 
# 
# dd1<-u_dat%>%dplyr::select(Cell.ID:file_string,V1,V2)
# 
# 
# 
# 
# 
# 
# 
# 
# ggplot( dd1,aes(-V1,V2,col=  dynamic_class3)) +  
#   geom_point(size=1.5)+scale_color_npg(name="Tumor response")+
#   theme_classic()+
#   theme(aspect.ratio=1)#axis.text=element_blank(),axis.title=element_blank(),legend.position
# 
# set.seed(1234)
# p1<-princomp ( dd1%>%dplyr::select(V1,V2) )
# u1<-umap ( dd1%>%dplyr::select(V1,V2) ,n_components=2)
# dd1$Dimension1 <- -p1$scores[,"Comp.1"]
# dd1$Dimension2 <- u1$layout[,1]
# #plot(Dimension2~V2,data=dd1)
# ggplot( dd1,aes(-V1,V2,col=  Dimension2)) +  
#   geom_point(size=1.5)+theme_classic()
# 
# ggplot( dd1,aes(y=Dimension2, x=dynamic_class3,col=dynamic_class3,fill=dynamic_class3)) +  
#   geom_violin(alpha=0.7)+  geom_point(size=1.5)+
#   scale_color_npg(name="Tumor response")+
#   theme_classic()+
#   theme(aspect.ratio=1)#axis.text=element_blank(),axis.title=element_blank(),legend.position
# 
# summardd1 <- data.table(dd1%>%group_by(Day,Patient.Study.ID,Sample,ARM,Celltype) %>%summarise(n= n(),V1=mean(V1), V2=mean(V2),Dimension1=mean(Dimension1),Dimension2=mean(Dimension2)  ))
# summardd2 <- merge(summardd1,unique(metadd%>%dplyr::select(Patient.Study.ID,prop_change,dynamic_class3))  , by="Patient.Study.ID")
# #ggplot( summardd2 , aes(y=V2,x=-V1,col=Dimension1)) +geom_point() 
# ggplot( data.table(summardd2[n>100][Day!=180]%>%group_by(Patient.Study.ID,prop_change,dynamic_class3)%>%summarise(Dimension2=weighted.mean(Dimension2,n),n=sum(n)) )[n>0]  ,
#         aes(x=-Dimension2, y=(prop_change)) )+
#   #geom_hline(yintercept=1,linetype="dashed")+
#   geom_smooth(method="gam",formula=y~s(x,k=3),level=0.8,col="black",fill="grey") + 
#   #stat_smooth(method="gam",formula=y~s(x,k=3),level=0.8,col="black",fill="grey")+
#   geom_point(aes(col=dynamic_class3,size=n)) +
#   theme_classic(base_size=22) + theme(aspect.ratio=1)+
#   coord_trans(y="exp")+
#   labs(y="Tumor size post treatment \n (relative to baseline)",x="Average myeloid M2 differentiation \n in first 14 days of treatment")+
#   scale_color_npg(name="Tumor response",labels=c("Resistant","Sensitive"))+
#   scale_size_continuous(name="Sample size")
# 
# ggsave(filename=paste0("/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cohort 2 results/Macrophage differentiation predicts tumor size response.png"),width=10,height=10)
# 
# Cohort1Extract <-function(){
#   load(file= "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/UMAP genes macrophage only ALLARMS/UMAP genes macrophage only ALLARMS.RData")
#   #load( file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Cohort2Metadata/Cohort2Metadata.RData") #metadd,cohort1metadd,annotation.file,compdataLU
#   
#   dd1<-u_dat%>%dplyr::select(Cell.ID:PhenoCelltype,V1,V2)
#   return(dd1)
# }
# ddCohort1<-Cohort1Extract()
# #p1<-princomp ( ddCohort1%>%dplyr::select(V1,V2) )
# #u1 <- umap ( ddCohort1%>%dplyr::select(V1,V2) ,n_components=2)
# #ddCohort1$Dimension1 <- -p1$scores[,"Comp.1"]
# #ddCohort1$Dimension2 <- u1$layout[,1]
# preds <- predict(u1,ddCohort1%>%dplyr::select(V1,V2)  ) 
# ddCohort1$Dimension2 <-preds[,1] 
# summardd1Cohort1 <- data.table(ddCohort1%>%group_by(Day,Patient.Study.ID,Sample,ARM,Celltype) %>%summarise(n= n(),V1=mean(V1), V2=mean(V2),Dimension2=mean(Dimension2)  ))
# summardd2Cohort1 <- merge(summardd1Cohort1,unique(cohort1metadd%>%dplyr::select(Patient.Study.ID,prop_change,dynamic_class3))  , by="Patient.Study.ID")
# #ggplot( summardd2Cohort1 , aes(y=V2,x=-V1,col=Dimension1)) +geom_point() 
# 
# ggplot( data.table(summardd2Cohort1[n>50][Day!=180]%>%group_by(Patient.Study.ID,prop_change,dynamic_class3)%>%summarise(Dimension1=NA,Dimension2=weighted.mean(Dimension2,n),n=sum(n)) )[n>0]  ,
#         aes(x=-Dimension2, y=(prop_change)) )+
#   #geom_hline(yintercept=1,linetype="dashed")+
#   geom_smooth(method="gam",formula=y~s(x,k=3),level=0.8,col="black",fill="grey") + 
#   #stat_smooth(method="gam",formula=y~s(x,k=3),level=0.8,col="black",fill="grey")+
#   geom_point(aes(col=dynamic_class3,size=n)) +
#   theme_classic(base_size=22) + theme(aspect.ratio=1)+
#   coord_trans(y="exp")+
#   labs(y="Tumor size post treatment \n (relative to baseline)",x="Average myeloid M2 differentiation \n in first 14 days of treatment")+
#   scale_color_npg(name="Tumor response",labels=c("Resistant","Sensitive"))+
#   scale_size_continuous(name="Sample size")
# 
# 
# combplot<-rbind(data.table(Cohort="Discovery",summardd2Cohort1[n>20][Day!=180]%>%group_by(Patient.Study.ID,prop_change,dynamic_class3)%>%summarise(Dimension1=NA,Dimension2=weighted.mean(Dimension2,n),n=sum(n)) )[n>0] ,
#                 data.table(Cohort="Validation", summardd2[n>100][Day!=180]%>%group_by(Patient.Study.ID,prop_change,dynamic_class3)%>%summarise(Dimension1=NA,Dimension2=weighted.mean(Dimension2,n),n=sum(n)) )[n>0] 
# )
# ggplot( combplot, aes(x=-Dimension2, y=(prop_change) ) )+
#   #geom_hline(yintercept=1,linetype="dashed")+
#   geom_smooth(method="gam",formula=y~s(x,k=3),level=0.8,col="black",fill="grey") + 
#   #stat_smooth(method="gam",formula=y~s(x,k=3),level=0.8,col="black",fill="grey")+
#   geom_point(aes(col=dynamic_class3,size=n, shape=Cohort)) +
#   theme_classic(base_size=10) + theme(aspect.ratio=1)+
#   # coord_trans(y="exp")+
#   labs(y="Tumor size post treatment \n (relative to baseline)",x="Average myeloid M2 differentiation \n in first 14 days of treatment")+
#   labs(y="Post treatment fraction of tumor remaining",x="Average myeloid M2 differentiation \n in first 14 days of treatment")+
#   scale_color_npg(name="Tumor response",labels=c("Resistant","Sensitive"))+
#   scale_size_continuous(name="Sample size")#+facet_wrap(~Cohort)
# 
# ggsave(filename=paste0("/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cohort 2 results/Cohort 1 and 2Macrophage differentiation predicts tumor size response.png"),width=5,height=5)
# 
# summary( lm(prop_change~ I(-Dimension2)+Cohort,combplot) )
# 
# combplot<-rbind(data.table(Cohort="Discovery",summardd2Cohort1[n>20][Day!=180]%>%group_by(Patient.Study.ID,prop_change,dynamic_class3)%>%summarise(Dimension1=NA,Dimension2=weighted.mean(Dimension2,n),n=sum(n)) )[n>0] ,
#                 data.table(Cohort="Validation", summardd2[n>20][Day!=180]%>%group_by(Patient.Study.ID,prop_change,dynamic_class3)%>%summarise(Dimension1=NA,Dimension2=weighted.mean(Dimension2,n),n=sum(n)) )[n>0] 
# )
# 
# combplot
# ggplot( combplot, aes(y=-Dimension2, x=(dynamic_class3) ) )+
#   geom_violin(aes(group=dynamic_class3,fill=dynamic_class3),alpha=0.25)+
#   geom_jitter(aes(col= dynamic_class3, shape= Cohort ), width=0.1,size=2.5) +
#   theme_classic(base_size=10) + theme(aspect.ratio=1)+
#   # coord_trans(y="exp")+
#   labs(y="Tumor size post treatment \n (relative to baseline)",x="Average myeloid M2 differentiation \n in first 14 days of treatment")+
#   labs(y="Post treatment fraction of tumor remaining",x="Average myeloid M2 differentiation \n in first 14 days of treatment")+
#   scale_color_npg(name="Tumor response",labels=c("Resistant","Sensitive"))+
#   scale_size_continuous(name="Sample size")#+facet_wrap(~Cohort)
# 
# 
# 
# joineddd <- rbind(
#   data.table(Cohort="Discovery" , ddCohort1%>%select(Treatment,Cell.ID,Patient.Study.ID,Day,ARM,dynamic_class3,Celltype,V1,V2) ),
#   data.table(Cohort="Vavlidation" , dd1%>%select(Treatment,Cell.ID,Patient.Study.ID,Day,ARM,dynamic_class3,Celltype,V1,V2))
# )
# 
# summardd1Cohort1 <- data.table(ddCohort1%>%group_by(Day,Patient.Study.ID,Sample,ARM,Celltype) %>%summarise(n= n(),V1=mean(V1), V2=mean(V2),Dimension2=mean(Dimension2)  ))
# 
# joineddd%>%group_by(Day,Patient.Study.ID,Sample,ARM,Celltype) %>%summarise(n= n(),V1=mean(V1), V2=mean(V2),Dimension2=mean(Dimension2)  ))
# joineddd$V2%>%hist()
# 
# 
# 
# dd1<-u_dat%>%dplyr::select(Cell.ID:file_string,V1,V2)
# 
# 
# 
