
rm(list = ls())
require(deSolve);require(data.table);require(dplyr);require(tidyr);require(ggplot2)
library(igraph)

cna0<-fread("~/Dropbox/FELINE Project/Data_analysis/FELINE_data_folder/VariantCalls_CNV_FACETS_gene_cnv/output/FELINE_FACETS.gene_cnv_call.short_list.txt")

# Patient data :: 108 patients had sufficient data to get an ok estimate of response
load(file="~/Dropbox/FELINE/Data_share/Clinical results/FELINE clinical response U54 shared.RData") #dd4,dd3,
load(file="~/Dropbox/Cancer_pheno_evo/data/FELINE/Classificatoin 2 phase Patient_clinical response classif smooth9.RData")
dim(Clin_resp_dd_class);dim(Clin_resp_dd)

Clin_resp_dd_classAdd<-merge(Clin_resp_dd, Clin_resp_dd_class%>%dplyr::select(Patient.Study.ID, Day,rgr_A:dynamic_class2),by=c("Patient.Study.ID", "Day"),all.x=TRUE)
Clin_resp_dd_classAdd[Patient.Study.ID=="001-122",dynamic_class:="Partial response" ]
Clin_resp_dd_classAdd[Patient.Study.ID=="001-122",dynamic_class2:="Partial response" ]
Clin_resp_dd_classAdd[Patient.Study.ID=="001-122",Dose_lab:="Constant low dose" ]
Clin_resp_dd_classAdd[Patient.Study.ID%in% c("001-137","001-116"),dynamic_class:="Stable disease" ]
Clin_resp_dd_classAdd[Patient.Study.ID%in% c("001-137","001-116"),dynamic_class2:="Stable disease" ]
Clin_resp_dd_classAdd[Patient.Study.ID%in% c("001-137","001-116"),Dose_lab:="Constant low dose" ]
Clin_resp_dd_classAdd[, dynamic_class3:="Non-response"];Clin_resp_dd_classAdd[dynamic_class2%in%c("Sustained response","Partial response"), dynamic_class3:="Response"]



res_file_46<-fread("~/Dropbox/FELINE Project/Data_analysis/Clinical_data/FELINE_patient_1_46.clinical_data.csv")
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

dd0 <- fread("/Users/jason/Dropbox/FELINE Project/Data_analysis/FELINE_data_folder/VariantCalls_pyclone_clonevol_merge_figure/input/FELINE_WES_cluster_frequency.txt")
dd2 <- unique(merge(unique(dd0%>%dplyr::select(-c(Cluster,Fraction, Parent))),
                    unique(response_code_dd%>%dplyr::select(Patient.Study.ID,Patient,prop_change,ARM, Ribo,Treatment,TreatLab,Dose_lab,dynamic_class3,TreatCode, TreatCodeOrd)),by="Patient"))

cna0[,Sample:=gsub("P","FEL0",sample)]
dd3 <- merge(cna0%>%dplyr::select(Sample,total_cn,purity,ploidy,cnlr.median,gene),dd2,by=c("Sample"))
dd3[, max_t:=max(Timepoint) ,by=c("Patient.Study.ID")]

intis_mu <- data.table(dd3[Timepoint==0]%>%group_by(Patient,gene)%>%dplyr::summarise(cnlr.median0=mean(cnlr.median),total_cn0=mean(total_cn), purity0=mean(purity)))

dd3.1<-merge(dd3,intis_mu,by=c("Patient","gene"))
write.csv(dd3.1,file="/Users/jason/Desktop/Feline CNA with metadata.csv")
g<-dd3.1$gene%>%unique
ggplot( dd3.1[gene=="ESR1"][max_t==180],aes(y=cnlr.median-cnlr.median0,x=Timepoint) )+geom_point()+geom_path(aes(group=Patient.Study.ID))+facet_wrap(~ARM)+geom_smooth(method="lm") + theme_classic()
ggplot( dd3.1[gene=="CDK6"][max_t==180],aes(y=total_cn,x=Timepoint) )+geom_point()+geom_path(aes(group=Patient.Study.ID))+facet_wrap(~ARM)+geom_smooth(method="lm")+ theme_classic()
ggplot( dd3.1[gene=="TP53"][max_t==180],aes(y=total_cn-total_cn0,x=Timepoint) )+geom_point()+geom_path(aes(group=Patient.Study.ID))+facet_wrap(~ARM)+geom_smooth(method="lm")+ theme_classic()


ggfig<-ggplot( dd3.1[gene%in%c("TP53","CDK6","ESR1","RB1")][max_t==180],aes(y=total_cn,x=Timepoint) )+geom_point()+geom_path(aes(group=Patient.Study.ID))+facet_grid(gene~ARM)+ theme_classic() #geom_smooth(method="lm")+
ggfig<-ggplot( dd3.1[max_t==180],aes(y=total_cn,x=Timepoint) )+geom_point()+geom_path(aes(group=Patient.Study.ID))+facet_grid(gene~ARM)+ theme_classic() #geom_smooth(method="lm")+
ggsave(ggfig,filename="/Users/jason/Desktop/Feline CNA.pdf",  height=25, width=8)
indd <- dd3.1[gene=="TP53"][max_t==180]
pois1<-glmer(total_cn~(Timepoint) +(1|Patient.Study.ID),family =poisson,data=indd)
summary(pois1)
indd$pred<-predict(pois1,type="response")
ggplot( indd,aes(y=total_cn,x=Timepoint) )+geom_point()+geom_path(aes(group=Patient.Study.ID))+#facet_wrap(~ARM)+
  geom_path(aes(y=pred,x=Timepoint),col="red",size=2)+ theme_classic()

ggplot( dd3.1[max_t==180],aes(y=cnlr.median,x=Timepoint) )+geom_point()+geom_path(aes(group=Patient.Study.ID))+facet_grid(gene~ARM)+geom_smooth(method="lm")+ theme_classic()
ggplot( dd3.1[max_t==180],aes(y=total_cn,x=Timepoint) )+geom_point()+geom_path(aes(group=Patient.Study.ID))+facet_grid(gene~ARM)+geom_smooth(method="lm")+ theme_classic()


