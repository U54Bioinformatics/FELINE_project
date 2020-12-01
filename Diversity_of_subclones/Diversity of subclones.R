rm(list = ls())
require(deSolve);require(data.table);require(dplyr);require(tidyr);require(ggplot2)
library(igraph)

# Patient data :: 108 patients had sufficient data to get an ok estimate of response
load(file="/Users/jason/Dropbox/FELINE/Data_share/Clinical results/FELINE clinical response U54 shared.RData") #dd4,dd3,
load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE/Classificatoin 2 phase Patient_clinical response classif smooth9.RData")
dim(Clin_resp_dd_class);dim(Clin_resp_dd)

Clin_resp_dd_classAdd<-merge(Clin_resp_dd, Clin_resp_dd_class%>%dplyr::select(Patient.Study.ID, Day,rgr_A:dynamic_class2),by=c("Patient.Study.ID", "Day"),all.x=TRUE)
Clin_resp_dd_classAdd[Patient.Study.ID=="001-122",dynamic_class:="Partial response" ]
Clin_resp_dd_classAdd[Patient.Study.ID=="001-122",dynamic_class2:="Partial response" ]
Clin_resp_dd_classAdd[Patient.Study.ID=="001-122",Dose_lab:="Constant low dose" ]
Clin_resp_dd_classAdd[Patient.Study.ID%in% c("001-137","001-116"),dynamic_class:="Stable disease" ]
Clin_resp_dd_classAdd[Patient.Study.ID%in% c("001-137","001-116"),dynamic_class2:="Stable disease" ]
Clin_resp_dd_classAdd[Patient.Study.ID%in% c("001-137","001-116"),Dose_lab:="Constant low dose" ]
#res_file_46[Patient.Study.ID%in%  unique(Clin_resp_dd_classAdd[is.na(dynamic_class)]$Patient.Study.ID) ]
Clin_resp_dd_classAdd[, dynamic_class3:="Non-response"];Clin_resp_dd_classAdd[dynamic_class2%in%c("Sustained response","Partial response"), dynamic_class3:="Response"]


res_file_46<-fread("/Users/jason/Dropbox/FELINE Project/Data_analysis/Clinical_data/FELINE_patient_1_46.clinical_data.csv")
response_dd <- data.table( unique(Clin_resp_dd_classAdd[!is.na(dynamic_class)]
                                  #Clin_resp_dd_class 
                                  %>% dplyr::select(Patient.Study.ID,prop_change,Arm,ARM,Ribo,Treatment,TreatLab,Class, Response,Burden_t0 ,BurdenTracked ,Day0, DayLastBurd,
                                                    rgr_A ,  rgr_B, classification_hclust, classification_hclust2, Type,TypeM,TypeM2, Dose_lab , dynamic_class , dynamic_class2,dynamic_class3)))
patientCode_LU <- data.table( fread("~/Dropbox/FELINE/Data_share/Modeling_Data/Patient_1_43.clinical_data.Bild_ID.TXT"))

patientCode_LU$Patient <-paste0("FEL0",1:nrow(patientCode_LU))
response_code_dd <- merge(patientCode_LU ,response_dd,by="Patient.Study.ID")
response_code_dd[,TreatCode:="RiboArm"]
response_code_dd[Treatment=="letrozole",TreatCode:="ControlArm"]
response_code_dd[,TreatCodeOrd:="a_RiboArm"]
response_code_dd[Treatment=="letrozole",TreatCodeOrd:="b_ControlArm"]


#dd0 <- fread("/Users/jason/Dropbox/FELINE Project/Data_analysis/WES/VariantCalls_pyclone_clonevol_merge_figure/input/FELINE_WES_cluster_frequency.txt")
dd0 <- fread("/Users/jason/Dropbox/FELINE Project/Data_analysis/FELINE_data_folder/VariantCalls_pyclone_clonevol_merge_figure/input/FELINE_WES_cluster_frequency.txt")

lu <- unique(dd0%>%dplyr::select(Patient, Timepoint))
lu2<- lu%>%dplyr::select(Patient)%>%unique()
lu2$pat_index<-1:nrow(lu2)
lu<- merge(lu,lu2,by="Patient")
#i<-11
dd1<-rbindlist(lapply(1:nrow(lu),function(i){
dd_pt_i <- merge(dd0,lu[i],c("Patient","Timepoint"))
g <- graph.data.frame(data.table("source"=dd_pt_i[Parent!=0]$Parent,"target"=dd_pt_i[Parent!=0]$Cluster,weight=log(dd_pt_i[Parent!=0]$Fraction /100)))
E(g)$weight <- -1*E(g)$weight
res <- exp(shortest.paths(g) * -1)
diag(res) <- 0
res1 <- round(res,digits =2)[,'1']  #res[,'1']
dd_pt_i2<-  merge(dd_pt_i, data.table(Cluster = as.numeric(names(shortest.paths(g)[,'1'])), unadjprop = res[,'1'],prop = res1/sum( res1 ) ) ,by="Cluster")
return(dd_pt_i2)
}))

dd2 <- merge(dd1,unique(response_code_dd%>%dplyr::select(Patient.Study.ID,Patient,prop_change,ARM, Ribo,Treatment,TreatLab,Dose_lab,dynamic_class3,TreatCode, TreatCodeOrd)),by="Patient")[Cluster!=1][prop>0]
dd2$ClusterID <- dense_rank( dd2$Cluster +  (max(dd2$Cluster)+1)*dd2$pat_index )

require(vegan)
H <- diversity(dd2[Patient=="FEL046"] $prop)
diversity(dd2[Patient=="FEL046"][Timepoint==0] $prop,index="simpson")
diversity(dd2[Patient=="FEL046"][Timepoint!=0] $prop,index="simpson")

diversity(dd2[Patient=="FEL046"][Timepoint==0] $prop)
diversity(dd2[Patient=="FEL046"][Timepoint!=0] $prop)

dd2[, H:=diversity(prop,index="shannon") ,by=c("Timepoint","Patient")] #increases as both the richness and the evenness increase
dd2[, Dominance:=diversity(prop,index="simpson") ,by=c("Timepoint","Patient")] #increasesdominance increases/evenness decreases
dd2[, Richness:=length(unique(ClusterID)) ,by=c("Timepoint","Patient")] #increasesdominance increases/evenness decreases
dd2[, max_t:=max(Timepoint) ,by=c("Patient")]
intis_mu <- data.table(dd2[Timepoint==0]%>%group_by(Patient)%>%dplyr::summarise(H0=mean(H),Dominance0=mean(Dominance), Richness0=mean(Richness)))

dd2.1<-merge(dd2,intis_mu,by="Patient")

ggplot( dd2[Patient=="FEL046"] ,aes(y=prop,x=Timepoint,fill=as.factor(ClusterID),group=Patient))+geom_bar(position="stack", stat="identity")

ggplot( dd2[max_t==180] ,aes(y=H,x=Timepoint,col=dynamic_class3,group=Patient))+geom_line()+facet_wrap(~ARM)+labs(y="Tumor heterogeneity (Shannon diversity)" , x="Time")
ggplot( dd2[max_t==180] ,aes(y=Dominance,x=Timepoint,col=dynamic_class3))+geom_line(aes(group=Patient))+facet_wrap(~Ribo)+labs(y="Tumor dominance (Simpson's index)" , x="Time")+geom_smooth()
ggplot( dd2[max_t==180] ,aes(y=Dominance,x=H,col=dynamic_class3,group=Patient))+geom_line()+facet_wrap(~ARM)+labs(y="Tumor heterogeneity (Shannon diversity)" , x="Time")
ggplot( dd2[max_t==180] ,aes(y=Dominance,x=Richness,col=dynamic_class3,group=Patient))+geom_point()+facet_wrap(~Ribo)
ggplot( dd2[max_t==180] ,aes(y=Dominance,x=H,col=dynamic_class3,group=Patient))+geom_point()+facet_wrap(~Ribo)

ggplot( dd2[max_t==180] ,aes(y=Dominance,x=Timepoint,col=dynamic_class3))+geom_line(aes(group=Patient))+facet_wrap(~Ribo)+labs(y="Tumor dominance (Simpson's index)" , x="Time")+geom_smooth(method="lm")
ggplot( dd2[max_t==180] ,aes(y=Dominance,x=Timepoint,col=as.factor(Ribo)))+geom_jitter(aes(group=Patient),height=0,width=10)+facet_wrap(~dynamic_class3)+labs(y="Tumor dominance (Simpson's index)" , x="Time")+geom_smooth(method="lm")
ggplot( dd2[max_t==180] ,aes(y=Dominance,x=Timepoint,col=as.factor(Ribo)))+geom_jitter(aes(group=Patient),height=0,width=3)+labs(y="Tumor dominance (Simpson's index)" , x="Time")+geom_smooth(method="lm")
ggplot( dd2[max_t==180] ,aes(y=H,x=Timepoint,col=as.factor(Ribo)))+geom_jitter(aes(group=Patient),height=0,width=10)+labs(y="Tumor dominance (Simpson's index)" , x="Time")+geom_smooth(method="lm")+theme_classic()
ggplot( dd2[max_t==180] ,aes(y=H,x=Timepoint,col=as.factor(Ribo)))+geom_jitter(aes(group=Patient),height=0,width=10)+facet_wrap(~dynamic_class3)+labs(y="Tumor dominance (Simpson's index)" , x="Time")+geom_smooth(method="lm")+theme_classic()
ggplot( dd2[max_t==180] ,aes(y=Richness,x=Timepoint,col=as.factor(Ribo)))+geom_jitter(aes(group=Patient),height=0,width=10)+facet_wrap(~dynamic_class3)+labs(y="Tumor dominance (Simpson's index)" , x="Time")+geom_smooth(method="lm")+theme_classic()

ggplot( dd2[max_t==180] ,aes(y=Richness,x=Timepoint,col=as.factor(Ribo)))+geom_jitter(aes(group=Patient),height=0,width=10)+labs(y="Tumor dominance (Simpson's index)" , x="Time")+geom_smooth(method="lm")+theme_classic()
ggplot( dd2[max_t==180] ,aes(y=H,x=Timepoint,col=as.factor(Ribo)))+geom_jitter(aes(group=Patient),height=0,width=10)+labs(y="Tumor dominance (Simpson's index)" , x="Time")+geom_smooth(method="lm")+theme_classic()
ggplot( dd2[max_t==180] ,aes(y=Dominance,x=Timepoint,col=as.factor(Ribo)))+geom_jitter(aes(group=Patient),height=0,width=10)+labs(y="Tumor dominance (Simpson's index)" , x="Time")+geom_smooth(method="lm")+theme_classic()

require(lme4)
require(lmerTest);
in_dd<-unique(dd2[max_t==180]%>%dplyr::select(H, Richness,Dominance,Timepoint,Ribo,dynamic_class3,Patient,TreatLab))
in_ddprint <- in_dd
in_ddprint[,simpson:=Dominance]
in_ddprint[,Dominance:=1-simpson]
m1 <- lmer(H~ Timepoint*Ribo +(1|Patient)+(0+Timepoint|Patient),data=in_dd)
m1 <- lmer(H~ Timepoint*Ribo +(1|Patient),data=dd2[max_t==180])
#m1<-lmer(H~ Timepoint*Ribo +Timepoint:dynamic_class3:Ribo +(1|Patient)+(0+Timepoint|Patient),data=unique(dd2[max_t==180]%>%dplyr::select(H, Timepoint,Ribo,dynamic_class3,Patient)))
#m1<-lmer(H~ Timepoint*Ribo+Timepoint*dynamic_class3 +(1+Timepoint|Patient),data=dd2[max_t==180])
#
#m1<-lmer(H~ Timepoint*Ribo +Timepoint:dynamic_class3:Ribo +(1+Timepoint|Patient),data=dd2[max_t==180])
#m1<-lmer(H~ Timepoint*Ribo +(1+Timepoint|Patient),data=dd2[max_t==180])
plot(m1)
summary(m1)
library(merTools)
m1 <- lmer(H~ Timepoint*Ribo +(1+Timepoint|Patient),data=dd2[max_t==180])
preds <- predictInterval(m1, newdata = dd2[max_t==180], n.sims = 999,which =  "fixed")
preds2<-cbind(dd2.1[max_t==180],preds)

intis_fit<-data.table(preds2[Timepoint==0]%>%group_by(Patient)%>%dplyr::summarise(fit0=mean(fit)))

preds2.1<-merge(preds2,intis_fit,by="Patient")

library("ggsci")
library(scatterpie)
ggplot( dd2[max_t==180] ,aes(y=H,x=Timepoint,col=as.factor(Ribo)))+geom_jitter(aes(group=Patient),height=0,width=10)+facet_wrap(~dynamic_class3)+labs(y="Tumor heterogeneity (Shannon diversity)" , x="Time")+
  geom_smooth(stat="identity",data=preds2,aes(y=fit,x=Timepoint,col=as.factor(Ribo),group=as.factor(Ribo)))+theme_classic()+
  geom_ribbon(data=preds2,aes(ymax=upr ,ymin= lwr,x=Timepoint,fill=as.factor(Ribo), group=as.factor(Ribo)),alpha=0.4,col=NA)

ggplot( dd2[max_t==180] ,aes(y=H,x=Timepoint,col=as.factor(TreatLab)))+geom_jitter(aes(group=Patient),height=0,width=10)+labs(y="Tumor heterogeneity (Shannon diversity)" , x="Time")+
  geom_smooth(stat="identity",data=preds2,aes(y=fit,x=Timepoint,col=as.factor(TreatLab),group=as.factor(TreatLab)))+theme_classic(base_size=18)+
  geom_ribbon(data=preds2,aes(ymax=upr ,ymin= lwr,x=Timepoint,fill=as.factor(TreatLab), group=as.factor(TreatLab)),alpha=0.4,col=NA)+
  scale_color_jco(name="Treatment")+scale_fill_jco(name="Treatment")#+facet_wrap(~dynamic_class3)#(palette = "Set2")

ggplot(dd2[max_t==180],aes(y=exp(H),x=Timepoint,col=as.factor(TreatLab)))+geom_jitter(aes(group=Patient),height=0,width=10)+labs(y="Tumor heterogeneity (Shannon diversity)" , x="Time")+
  geom_smooth(stat="identity",data=preds2,aes(y=exp(fit),x=Timepoint,col=as.factor(TreatLab),group=as.factor(TreatLab)))+theme_classic(base_size=18)+
  geom_ribbon(data=preds2,aes(ymax=exp(upr) ,ymin= exp(lwr),x=Timepoint,fill=as.factor(TreatLab), group=as.factor(TreatLab)),alpha=0.4,col=NA)+
  scale_color_jco(name="Treatment")+scale_fill_jco(name="Treatment")#+facet_wrap(~dynamic_class3)#(palette = "Set2")


ggplot(dd2.1[max_t==180],aes(y=(H)-(H0),x=Timepoint,col=as.factor(TreatLab)))+geom_jitter(aes(group=Patient),height=0,width=10)+labs(y="Tumor heterogeneity (Shannon diversity)" , x="Time")+
  geom_smooth(stat="identity",data=preds2.1,aes(y=(fit)-(fit0),x=Timepoint,col=as.factor(TreatLab),group=as.factor(TreatLab)))+theme_classic(base_size=18)+
  geom_ribbon(data=preds2.1,aes(ymax=(upr)-(fit0) ,ymin= (lwr)-(fit0),x=Timepoint,fill=as.factor(TreatLab), group=as.factor(TreatLab)),alpha=0.4,col=NA)+
  scale_color_jco(name="Treatment")+scale_fill_jco(name="Treatment")#+facet_wrap(~dynamic_class3)#(palette = "Set2")


# ggplot( dd2[max_t==180] ,aes(y=exp(H),x=Timepoint,col=as.factor(TreatLab)))+geom_jitter(aes(group=Patient),height=0,width=10)+labs(y="Tumor heterogeneity (Shannon diversity)" , x="Time")+
#   geom_smooth(stat="identity",data=preds2,aes(y=exp(fit),x=Timepoint,col=as.factor(TreatLab),group=as.factor(TreatLab)))+theme_classic(base_size=18)+
#   geom_ribbon(data=preds2,aes(ymax=exp(upr) ,ymin= exp(lwr),x=Timepoint,fill=as.factor(TreatLab), group=as.factor(TreatLab)),alpha=0.4,col=NA)+
#   scale_color_jco(name="Treatment")+scale_fill_jco(name="Treatment")#+facet_wrap(~dynamic_class3)#(palette = "Set2")

ww<-dd2[max_t==180]%>%spread(Cluster,prop,fill =0)
setnames(ww,old=as.character(2:7),new=LETTERS[1:6])
ggplot( dd2[max_t==180] ,aes(y=exp(H),x=log(Timepoint+1),col=as.factor(TreatLab)))+
  #geom_jitter(aes(group=Patient),height=0,width=10)+
  labs(y="Tumor heterogeneity (Shannon diversity)" , x="Time")+
  geom_smooth(stat="identity",data=preds2,aes(y=exp(fit),x=log(Timepoint+1),col=as.factor(TreatLab),group=as.factor(TreatLab)))+theme_classic(base_size=18)+
  geom_ribbon(data=preds2,aes(ymax=exp(upr) ,ymin= exp(lwr),x=log(Timepoint+1),fill=as.factor(TreatLab), group=as.factor(TreatLab)),alpha=0.4,col=NA)+
  scale_color_jco(name="Treatment")+scale_fill_jco(name="Treatment")+ coord_equal()+ 
  geom_scatterpie(data=data.frame(ww), aes(y=exp(H),x=log(Timepoint+1)+(pat_index-mean(pat_index))/12, group=interaction(Patient,Timepoint), r=0.2),pie_scale = 2 , cols=LETTERS[1:6],  alpha=.8) +
  facet_wrap(~TreatLab)
+
  geom_scatterpie_legend(d$radius, x=-160, y=-55)


m2 <- lmer(Dominance~ Timepoint*Ribo +(1|Patient),data=dd2[max_t==180])
summary(m2)
m2 <- lmer(Dominance~ Timepoint*Ribo +(1+Timepoint|Patient),data=dd2[max_t==180])
predsB <- predictInterval(m2, newdata = dd2[max_t==180], n.sims = 999,which =  "fixed")
preds3<-cbind(dd2.1[max_t==180],predsB)
intis_fit3<-data.table(preds3[Timepoint==0]%>%group_by(Patient)%>%dplyr::summarise(fit0=mean(fit)))


dd3<-dd2.1[max_t==180]

dd3[,is_dom:=prop==max(prop),by=c("Patient","Timepoint")]

dd3[,is_dom0:=is_dom*(Timepoint==0),by=c("Patient","Timepoint")]
dd3[,is_dom180:=is_dom*(Timepoint==180),by=c("Patient","Timepoint")]
dd3[,Clust_dom180:=sum(is_dom180)>0,by=c("Patient","Cluster")]
dd3[,Clust_dom0:=sum(is_dom0)>0,by=c("Patient","Cluster")]

preds3.1<-merge(preds3,intis_fit3,by="Patient")


mR <- lmer(log(Richness)~ Timepoint*Ribo +(1+Timepoint|Patient),data=dd2[max_t==180])
predsC <- exp(predictInterval(mR, newdata = dd2[max_t==180], n.sims = 999,which =  "fixed"))
preds4<-cbind(dd2.1[max_t==180],predsC)
intis_fit4<-data.table(preds4[Timepoint==0]%>%group_by(Patient)%>%dplyr::summarise(fit0=mean(fit)  ))
preds4.1<-merge(preds4,intis_fit4,by="Patient")


dd3<-dd2.1[max_t==180]

dd3[,is_dom:=prop==max(prop),by=c("Patient","Timepoint")]

dd3[,is_dom0:=is_dom*(Timepoint==0),by=c("Patient","Timepoint")]
dd3[,is_dom180:=is_dom*(Timepoint==180),by=c("Patient","Timepoint")]
dd3[,Clust_dom180:=sum(is_dom180)>0,by=c("Patient","Cluster")]
dd3[,Clust_dom0:=sum(is_dom0)>0,by=c("Patient","Cluster")]

preds3.1<-merge(preds3,intis_fit3,by="Patient")



ggplot(in_dd ,aes(y=1-Dominance,x=Timepoint,col=as.factor(TreatLab)))+geom_jitter(aes(group=Patient),height=0,width=10)+labs(y="Tumor dominance (1-Simpson's index)" , x="Time")+
  geom_smooth(stat="identity",data=preds3,aes(y=1-fit,x=Timepoint,col=as.factor(TreatLab),group=as.factor(TreatLab)))+theme_classic(base_size=18)+
  geom_ribbon(data=preds3,aes(ymax=1-upr ,ymin= 1-lwr,x=Timepoint,fill=as.factor(TreatLab), group=as.factor(TreatLab)),alpha=0.4,col=NA)+
  scale_color_jco(name="Treatment")+scale_fill_jco(name="Treatment")+theme(aspect.ratio=1)#+facet_wrap(~dynamic_class3)#(palette = "Set2")

ww3<-dd3[max_t==180]%>%spread(Cluster,prop,fill =0)
setnames(ww,old=as.character(2:7),new=LETTERS[1:6])

ww3<-dd3[max_t==180]%>%spread(Clust_dom0,prop,fill =0)
setnames(ww3,old=c("FALSE","TRUE"),new=c("NonDom0","Dom0"))


ggplot( dd3[max_t==180] ,aes(y=Dominance,x=log(Timepoint+1),col=as.factor(TreatLab)))+
  #geom_jitter(aes(group=Patient),height=0,width=10)+
  labs(y="Tumor dominance (Simpson's index)" , x="Time")+
  geom_smooth(stat="identity",data=preds3,aes(y=(fit),x=log(Timepoint+1),col=as.factor(TreatLab),group=as.factor(TreatLab)))+theme_classic(base_size=18)+
  geom_ribbon(data=preds3,aes(ymax=(upr) ,ymin= (lwr),x=log(Timepoint+1),fill=as.factor(TreatLab), group=as.factor(TreatLab)),alpha=0.4,col=NA)+
  scale_color_jco(name="Treatment")+scale_fill_jco(name="Treatment")+ coord_equal()+ 
  geom_scatterpie(data=data.frame(ww3), aes(y=Dominance+(pat_index-mean(pat_index))/10,x=log(Timepoint+1)+(pat_index-mean(pat_index))/4, group=interaction(Patient,Timepoint), r=0.2),pie_scale = 2 , cols=c("NonDom0","Dom0"),  alpha=.8) +
  facet_wrap(~TreatLab)




ggplot(dd3[max_t==180],aes(y=(Dominance),x=Timepoint,col=as.factor(TreatLab)))+geom_jitter(aes(group=Patient),height=0,width=10)+labs(y="Tumor dominance (Simpson's index)" , x="Time")+
  geom_smooth(stat="identity",data=preds3.1,aes(y=(fit),x=Timepoint,col=as.factor(TreatLab),group=as.factor(TreatLab)))+theme_classic(base_size=18)+
  geom_ribbon(data=preds3.1,aes(ymax=(upr) ,ymin= (lwr),x=Timepoint,fill=as.factor(TreatLab), group=as.factor(TreatLab)),alpha=0.4,col=NA)+
  scale_color_jco(name="Treatment")+scale_fill_jco(name="Treatment")#+facet_wrap(~dynamic_class3)#(palette = "Set2")


ggplot(dd3[max_t==180],aes(y=(1-Dominance)-(1-Dominance0),x=Timepoint,col=as.factor(TreatLab)))+geom_jitter(aes(group=Patient),height=0,width=10)+labs(y="Tumor dominance (1-Simpson's index)" , x="Time")+
  geom_smooth(stat="identity",data=preds2.1,aes(y=(1-fit)-(1-fit0),x=Timepoint,col=as.factor(TreatLab),group=as.factor(TreatLab)))+theme_classic(base_size=18)+
  geom_ribbon(data=preds2.1,aes(ymax=(1-upr)-(1-fit0) ,ymin= (1-lwr)-(1-fit0),x=Timepoint,fill=as.factor(TreatLab), group=as.factor(TreatLab)),alpha=0.4,col=NA)+
  scale_color_jco(name="Treatment")+scale_fill_jco(name="Treatment")#+facet_wrap(~dynamic_class3)#(palette = "Set2")



ggplot(dd2.1[max_t==180],aes(y=(H)-(H0),x=Timepoint,col=as.factor(TreatLab)))+geom_jitter(aes(group=Patient),height=0,width=10)+labs(y="Tumor heterogeneity (Shannon diversity)" , x="Time")+
  geom_smooth(stat="identity",data=preds2.1,aes(y=(fit)-(fit0),x=Timepoint,col=as.factor(TreatLab),group=as.factor(TreatLab)))+theme_classic(base_size=18)+
  geom_ribbon(data=preds2.1,aes(ymax=(upr)-(fit0) ,ymin= (lwr)-(fit0),x=Timepoint,fill=as.factor(TreatLab), group=as.factor(TreatLab)),alpha=0.4,col=NA)+
  scale_color_jco(name="Treatment")+scale_fill_jco(name="Treatment")#+facet_wrap(~dynamic_class3)#(palette = "Set2")


ddFULL <- merge(dd3,dd2.1,by=names(dd3)[names(dd3)%in%names(dd2.1)])

merg2 <- preds2.1%>%mutate(metric="H",y=(H)-(H0),fitted=(fit)-(fit0),ymax=(upr)-(fit0),ymin=(lwr)-(fit0))
merg3 <- preds3.1%>%mutate(metric="iDprime",y=(1-Dominance)-(1-Dominance0),fitted=(1-fit)-(1-fit0),ymax=(1-upr)-(1-fit0),ymin=(1-lwr)-(1-fit0))
merg4 <- preds4.1%>%mutate(metric="Richness",y=(Richness)-(Richness0),fitted=(fit)-(fit0),ymax=(upr)-(fit0),ymin=(lwr)-(fit0))
#write.csv(rr,file="/Users/jason/Desktop/Shannon index values.csv")


rr<- unique(data.table(merg3%>%dplyr::select(Patient,Timepoint,ARM,Dominance,Dominance0)%>%mutate(iDominance=(1-Dominance),iDominance0=(1-Dominance0), deltaDominance=(1-Dominance)-(1-Dominance0))))[order(-deltaDominance)]
write.csv(rr,file="/Users/jason/Desktop/Dominance index values.csv")
#write.csv(rr,file="/Users/jason/Desktop/Shannon index values.csv")
preds4.1

predsFULL<-data.table(rbind(rbind(merg2,merg3),merg4))
predsFULL[,flab:="Diversity \n (Shannon index)"]
predsFULL[metric=="iDprime",flab:="Dominance \n (1-Simpson's index)"]
predsFULL[metric=="Richness",flab:="Richness"]

ggplot(predsFULL[Timepoint==180],aes(y=y,x=flab,col=as.factor(TreatLab)))+
  geom_violin(aes(group=interaction(flab,TreatLab),fill=TreatLab),alpha=0.4,col=NA)+
  geom_point(aes(group=TreatLab),position=position_dodge(width=1))+
  labs(y="Change in tumor heterogeneity over trial" , x="Heterogeneity metric")+
  geom_point(size=6,stat="identity",data=predsFULL[Timepoint==180],aes(y=fitted,x=flab,col=as.factor(TreatLab),group=interaction(flab,TreatLab)),position=position_dodge(width=1))+
  theme_classic(base_size=18)+
  geom_errorbar(data=predsFULL[Timepoint==180],aes(ymax=ymax ,ymin=ymin,x=flab,fill=as.factor(TreatLab), group=as.factor(TreatLab)),alpha=0.4,position=position_dodge(width=1))+
  scale_color_jco(name="Treatment")+scale_fill_jco(name="Treatment")+#+facet_wrap(~dynamic_class3)#(palette = "Set2")
  #scale_x_discrete(labels= c("Diversity \n (Shannon index)","Dominance \n (1-Simpson's index)","Richness")) +
  theme(aspect.ratio=1)+facet_wrap(~flab,scales="free")+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())+
  theme(legend.text =element_blank(),
        legend.title =element_blank())



ggplot(ddFULL[Timepoint==180],aes(y=y,x=Timepoint,col=as.factor(TreatLab)))+geom_jitter(aes(group=Patient),height=0,width=10)+labs(y="Tumor heterogeneity (Shannon diversity)" , x="Time")+
  geom_smooth(stat="identity",data=predsFULL[Timepoint==180],aes(y=(fit)-(fit0),x=Timepoint,col=as.factor(TreatLab),group=as.factor(TreatLab)))+theme_classic(base_size=18)+
  geom_ribbon(data=predsFULL[Timepoint==180],aes(ymax=(upr)-(fit0) ,ymin= (lwr)-(fit0),x=Timepoint,fill=as.factor(TreatLab), group=as.factor(TreatLab)),alpha=0.4,col=NA)+
  scale_color_jco(name="Treatment")+scale_fill_jco(name="Treatment")#+facet_wrap(~dynamic_class3)#(palette = "Set2")



