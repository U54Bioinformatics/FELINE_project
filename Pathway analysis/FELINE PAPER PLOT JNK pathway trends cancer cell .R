rm(list=ls())
require(mgcv);require(data.table);require(dplyr);require(ggplot2);require(tidyr);require(parallel);require(lme4);require(lmerTest);require(parallel)
library(effects)#require(umap)

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

response_code_dd <- merge(patientCode_LU ,response_dd,by="Patient.Study.ID")
response_code_dd[,TreatCode:="RiboArm"]
response_code_dd[Treatment=="letrozole",TreatCode:="ControlArm"]
response_code_dd[,TreatCodeOrd:="a_RiboArm"]
response_code_dd[Treatment=="letrozole",TreatCodeOrd:="b_ControlArm"]


gsea_path<-"~/Dropbox/FELINE Project/Data_analysis/scRNA/05_ssGSEA_score/Signature_c2_hallmark/results/"   #"~/Dropbox/FELINE/Data_share/Modeling_Data/pathway_seperated_files/Data_gene_count_per_celltype_model_zinbwave_ssGSEA/FEL001043/"

fit_mod<-function(DATA,ARM_I,DAY){
  dd_i <- DATA[ARM==ARM_I&Day%in%c(0,DAY)]
  responseMod <- lmer(ssGSEA ~ -1+dynamic_class2*day_fact + (1+day_fact|Patient.Study.ID), REML=FALSE,data= dd_i) ### responseMod <- lmer(ssGSEA ~ ln_prop_change*day_fact*Ribo + (1|Patient.Study.ID) + (1|day_fact), data= path_clin_long)
  #if( data.frame(anova(responseMod))[3,"Pr..F."]<0.05 ){
  #  dd_i$pred<- predict(responseMod,re.form=NA)
  #  dd_i<<-dd_i
  #  conf_ints<-data.table(
  #    expand.grid(dynamic_class2=levels(dd_i$dynamic_class2),
  #                day_fact=sort(unique(dd_i$day_fact))), 
  #    with(effect("dynamic_class2*day_fact",responseMod),cbind(lower,upper)))
  #  dd_ii <- merge(dd_i,conf_ints,by=c("dynamic_class2","day_fact"))
  #  p1<- ggplot(dd_ii[],aes(y=ssGSEA,x=dynamic_class,group=interaction(day_fact,dynamic_class),col=dynamic_class, fill=dynamic_class))+geom_violin(alpha=0.75)+
  #    facet_wrap(TreatLab~Day)+theme_classic(base_size=18)+#group=Patient.Study.ID,
  #    labs(y= paste0(path_clin_long[1]$Gene_Set, " in \n",path_clin_long[1]$Celltype) , x="Day")+
  #    theme(aspect.ratio=1)+scale_color_discrete(name="Response" )+scale_fill_discrete(name="Response" ) +
  #    geom_point(aes(y=pred))+
  #    geom_errorbar(aes(ymax=V2, ymin=V1))+
  #    geom_hline(yintercept = mean(dd_ii$ssGSEA),linetype=2,col="grey")
  #  ggsave(p1,file=paste0(paste(print_pics_loc,dd_ii[1]$Celltype,dd_ii[1]$Gene_Set,DAY,ARM_I,sep="_"),".png") ,height=10, width=30)
  #}
  output <- data.frame(summary(responseMod)$ coefficients);  output$Effect <- rownames(output)
  output$Intermediate_Clusters_ML <-dd_i[1]$Celltype   ;  output$pathwayscore <- dd_i[1]$Gene_Set
  output <- gather(output, statistic, value,"Estimate":"Pr...t..")
  output$metric <- "coefs"
  #anova info
  output2 <-   as.data.table(anova(responseMod)%>%dplyr::select(-c("NumDF","DenDF","Sum Sq","Mean Sq") ));  output2$Effect <- rownames(anova(responseMod))
  output2 <- na.omit(output2)
  output2$Intermediate_Clusters_ML <- dd_i[1]$Celltype ;    output2$pathwayscore <- dd_i[1]$Gene_Set
  output2 <- gather(output2, statistic, value,"F value":"Pr(>F)")
  output2$metric <- "signifcov"
  # join info
  outputFull <- data.table(rbind(output,output2))
  outputFull[,Day:=DAY]
  outputFull[,ARM:=ARM_I]
  return(outputFull)
}

### We will repeat for each cell type
cell_types_all <-c("Cancer_cells")
iii<- 1
#iii<- 7



  cell_type_i <- cell_types_all[iii]
  
  #load cell metadata
  file.exists("~/Dropbox/FELINE Project/Data_analysis/scRNA/01_metadata/results/FEL001046_meta_for_each_celltype/FEL001046_Cancer_cells_scRNA.metadata.txt")
  annotation.file<-paste0("~/Dropbox/FELINE Project/Data_analysis/scRNA/01_metadata/results/FEL001046_meta_for_each_celltype/FEL001046_",cell_type_i,"_scRNA.metadata.txt")
  cell_type_meta_dd <- data.table(  fread(annotation.file))[Celltype_subtype!="Low-quality cells"]
  setnames( cell_type_meta_dd,old="Sample",new="Sample_p_t")
  cell_type_meta_dd[, c("Sample", "Timepoint") := tstrsplit(Sample_p_t, "_", fixed=TRUE)]
  cell_type_meta_dd[,Day:=0]  ; cell_type_meta_dd[grepl("_M",Sample_p_t),Day:=14]   ; cell_type_meta_dd[grepl("_E",Sample_p_t),Day:=180]   
  cell_type_meta_dd[,day_fact:=as.factor(Day)]
  
  # merge response scores
  FULL_cell_type_meta_dd <- merge(cell_type_meta_dd,response_code_dd,by="Sample")        #FULL_cell_type_meta_dd[Celltype_subtype=="Fibroblasts",Celltype_subtype:="CAF-S1"]
  print(unique(FULL_cell_type_meta_dd$Celltype_subtype))
  
  # Load adjusted gsea scores
  cell_type_gsea_adj <-data.table( readRDS(file=paste0(gsea_path,"FEL001046_",cell_type_i,"_scRNA.zinbwave.normalized.ssGSEA_scores.RDS" ))) 
  setnames( cell_type_gsea_adj,old="Gene Set",new="Gene_Set") #cell_type_gsea_adj[1:5,1:5]
  #CELL_Subtype<-"Cancer cells"
  #path_i <-grep("HALLMARK_",cell_type_gsea_adj$Gene_Set)[15]
  cell_type_gsea_adj$Gene_Set[grep("HALLMARK_",cell_type_gsea_adj$Gene_Set)]
  CELL_Subtype<-"Cancer cells"
FULL_cell_type_meta_dd_CT <- FULL_cell_type_meta_dd[Celltype_subtype==CELL_Subtype]
cell_type_gsea_adj_CT<-cell_type_gsea_adj[,c("Gene_Set",FULL_cell_type_meta_dd_CT$Cell.ID),with=FALSE]  
      
      
#lll<-"SMID_BREAST_CANCER_BASAL_UP"   #II<-1;llll<-lll[1]
#lll<-"HALLMARK_ESTROGEN_RESPONSE_EARLY"
lll<-"ST_JNK_MAPK_PATHWAY"

llll<-lll[sapply(1:length(lll),function(i){ lll[i]%in%
  cell_type_gsea_adj$Gene_Set%>%unique})]
II=1
require(merTools)

pthese<-rbindlist(lapply(c(1,3),function(II){
path_long <- gather(cell_type_gsea_adj[Gene_Set==llll[II],],Cell.ID,ssGSEA,-Gene_Set)
path_clin_long <- data.table(merge( path_long , FULL_cell_type_meta_dd_CT, by= "Cell.ID") )
pathway_NAME <- path_clin_long[1,]$Gene_Set
cell_NAME <- path_clin_long[1,]$Celltype_subtype
path_clin_long[,DayLab:="Day 0"]
path_clin_long[Day==14, DayLab:="Day 14"]
path_clin_long[Day==180, DayLab:="Day 180"]
mu_ssgsea_t0<-data.table(path_clin_long[Day==0] %>% group_by(Patient.Study.ID) %>% dplyr::summarise(mu_ssGSEA_t0=mean(ssGSEA)))
path_clin_long<-merge(path_clin_long,mu_ssgsea_t0,by="Patient.Study.ID",all.x=TRUE)
path_clin_long[is.na(mu_ssGSEA_t0),mu_ssGSEA_t0:=mean(path_clin_long[Day==0]$ssGSEA)]
path_clin_long[,standardised_ssGSEA:=ssGSEA-mu_ssGSEA_t0]
path_clin_long[,ln_prop_change:=log(1+prop_change)]
path_clin_long[,exp_rgr_A:=exp(91*rgr_A)]
path_clin_long[,exp_rgr_B:=exp(89*rgr_B)]
path_clin_long[,class_rgr_A:="Response"]
path_clin_long[exp_rgr_A>0.8,class_rgr_A:="Non_response"]
path_clin_long$Patient.Study.ID <-as.factor(path_clin_long$Patient.Study.ID)

NEWXDAT2 <-rbind( data.table( unique(path_clin_long%>%dplyr::select(Patient.Study.ID,ARM,dynamic_class3)),day_fact="0"),
                  data.table( unique(path_clin_long%>%dplyr::select(Patient.Study.ID,ARM,dynamic_class3)),day_fact="14"),
                  data.table( unique(path_clin_long%>%dplyr::select(Patient.Study.ID,ARM,dynamic_class3)),day_fact="180"))
NEWXDAT2$Patient.Study.ID <-as.factor(NEWXDAT2$Patient.Study.ID)

Mod1 <- lmer(ssGSEA ~ -1+dynamic_class3*day_fact*ARM + (1+day_fact|Patient.Study.ID), REML=FALSE,data= path_clin_long) 
summary(Mod1)
# ModA <- lmer(ssGSEA ~ -1+dynamic_class3*day_fact + (1+day_fact|Patient.Study.ID), REML=FALSE,data= path_clin_long[ARM=="A"]) 
# summary(ModA)
# ModB <- lmer(ssGSEA ~ -1+dynamic_class2*day_fact + (1+day_fact|Patient.Study.ID), REML=FALSE,data= path_clin_long[ARM=="B"]) 
# summary(ModB)
# ModB <- lmer(ssGSEA ~ day_fact + (1+day_fact|Patient.Study.ID), REML=FALSE,data= path_clin_long[ARM=="B"][Day!=14]) 
# summary(ModB)
# ModC <- lmer(ssGSEA ~ -1+dynamic_class3*day_fact + (1+day_fact|Patient.Study.ID), REML=FALSE,data= path_clin_long[ARM=="C"]) 
# summary(ModC)
# ModResp <- lmer(ssGSEA ~ day_fact + (1+day_fact|Patient.Study.ID), REML=FALSE,data= path_clin_long[dynamic_class3=="Response"][ARM!="A"]) 
# summary(ModResp)
# ModBResp <- lmer(ssGSEA ~ day_fact*dynamic_class3 + (1+day_fact|Patient.Study.ID), REML=FALSE,data= path_clin_long[][ARM=="B"]) 
# summary(ModBResp)

PRED.lme4 <-  data.table(y = predict(Mod1, newdata=NEWXDAT2,re.form = NA) ,NEWXDAT2)
PRED.lme4 <-  data.table(y = predict(Mod1, newdata=NEWXDAT2,re.form = NA) ,
                         yi = predict(Mod1, newdata=NEWXDAT2,re.form =~(1+day_fact|Patient.Study.ID)) ,NEWXDAT2)

cinfreg<-data.table( predictInterval(Mod1,newdata=PRED.lme4, which ="fixed", level = 0.95) )

PRED.lme4<-data.table( data.table(cbind(PRED.lme4,cinfreg)) %>% group_by(ARM,day_fact,dynamic_class3)%>%dplyr::mutate(lwr=min(fit),upr=max(fit),fit=mean(fit) ))
PRED.lme4<-data.table(PRED.lme4 %>% group_by(ARM,day_fact,dynamic_class3)%>%dplyr::mutate(lwr=min(fit),upr=max(fit),fit=mean(fit) ))




PRED.lme4[,pathwayscore:=pathway_NAME]
#PRED.lme4[,PathType:=unique_path_annot[pathwayscore==pathway_NAME]$pathtype]
PRED.lme4[,Intermediate_Clusters_ML:=FULL_cell_type_meta_dd_CT[1]$Celltype]
PRED.lme4[,is_day0:=0];PRED.lme4[day_fact=="0",is_day0:=1]
PRED.lme4[, y_init:=sum(is_day0*yi), by=c("Patient.Study.ID","ARM","dynamic_class3")]
PRED.lme4[, mu_init:=sum(is_day0*y), by=c("Patient.Study.ID","ARM","dynamic_class3")]

PRED.lme4[Patient.Study.ID=="001-101"]
PRED.lme4[, upr1:=quantile(yi-y_init,probs =0.8), by=c("day_fact","ARM","dynamic_class3")]
PRED.lme4[, lwr1:=quantile(yi-y_init,probs =0.2), by=c("day_fact","ARM","dynamic_class3")]
PRED.lme4[,dum_day:=0]
PRED.lme4[day_fact==14,dum_day:=1]
PRED.lme4[day_fact==180,dum_day:=2]
}))


PRED.lme4[,Treat:="Letrozole alone"]
PRED.lme4[ARM=="B",Treat:="Intermittent high \n dose ribociclib"]
PRED.lme4[ARM=="C",Treat:="Continuous low \n dose ribociclib"]
PRED.lme4$Treat <- factor(PRED.lme4$Treat, levels = c("Letrozole alone", "Intermittent high \n dose ribociclib", "Continuous low \n dose ribociclib"))

#PRED.lme4[, mu_i_init:=sum(is_day0*yi), by=c("Patient.Study.ID","ARM","dynamic_class3","pathwayscore")]
# ggplot(PRED.lme4,
#        aes(y=y-mu_init,x=dum_day,group=dynamic_class3,col=dynamic_class3))+
#   geom_line(aes(group=interaction(pathwayscore,dynamic_class3,ARM)),size=2,alpha=1)+#geom_line(aes(y=yi-y_init,group=interaction(pathwayscore,dynamic_class3,Patient.Study.ID)),size=1,alpha=1)+
#   geom_ribbon(aes(ymax=upr-mu_init,ymin=lwr-mu_init,fill=dynamic_class3),alpha=0.4,col=NA)+
#   facet_grid(.~ARM,scales="free") + theme_classic(base_size=20)+
#   scale_color_discrete(name="Response",labels=c("Non-responder","Responder"))+
#   scale_fill_discrete(name="Response",labels=c("Non-responder","Responder"))+
#   scale_x_continuous(name="Day",breaks=0:2,labels=c(0,14,180))+
#   ylab(paste0(llll[II])) + theme(aspect.ratio=1)+
#   geom_line(aes(y=yi-mu_i_init,x=dum_day,group=Patient.Study.ID),linetype=2)

PRED.lme4[, mu_i_init:=sum(is_day0*yi), by=c("Patient.Study.ID","ARM","dynamic_class3","pathwayscore")]
ggplot(PRED.lme4,
       aes(y=y-mu_init,x=dum_day,group=dynamic_class3,col=dynamic_class3))+
  geom_line(aes(group=interaction(pathwayscore,dynamic_class3,ARM)),size=2,alpha=1)+#geom_line(aes(y=yi-y_init,group=interaction(pathwayscore,dynamic_class3,Patient.Study.ID)),size=1,alpha=1)+
  geom_ribbon(aes(ymax=upr1,ymin=lwr1,fill=dynamic_class3),alpha=0.4,col=NA)+
  facet_grid(.~Treat,scales="free") + theme_classic(base_size=20)+
  scale_color_discrete(name="Response",labels=c("Non-responder","Responder"))+
  scale_fill_discrete(name="Response",labels=c("Non-responder","Responder"))+
  scale_x_continuous(name="Day",breaks=0:2,labels=c(0,14,180))+
  ylab(paste0(llll[II])) + theme(aspect.ratio=1)+
  geom_line(aes(y=yi-mu_i_init,x=dum_day,group=Patient.Study.ID),linetype=2)+ theme(aspect.ratio=1)

ggplot(PRED.lme4,
       aes(y=y-mu_init,x=dum_day,group=dynamic_class3,col=dynamic_class3))+
  geom_line(aes(group=interaction(pathwayscore,dynamic_class3,ARM)),size=2,alpha=1)+#geom_line(aes(y=yi-y_init,group=interaction(pathwayscore,dynamic_class3,Patient.Study.ID)),size=1,alpha=1)+
  geom_ribbon(aes(ymax=upr1,ymin=lwr1,fill=dynamic_class3),alpha=0.4,col=NA)+
  facet_grid(.~Treat,scales="free") + theme_classic(base_size=20)+
  scale_color_discrete(name="Response",labels=c("Non-responder","Responder"))+
  scale_fill_discrete(name="Response",labels=c("Non-responder","Responder"))+
  scale_x_continuous(name="Day",breaks=0:2,labels=c(0,14,180))+
  ylab(paste0(llll[II])) + theme(aspect.ratio=1)+
  geom_line(aes(y=yi-mu_i_init,x=dum_day,group=Patient.Study.ID),linetype=2)+ theme(aspect.ratio=1)+
  theme(
    strip.background = element_blank(),
    strip.text = element_blank()
  )+
  theme(
    axis.title=element_blank(),
    axis.text=element_blank(),
    #axis.ticks=element_blank(),
    legend.position = "none"
  )




PRED.lme4[, mu_i_init:=sum(is_day0*yi), by=c("Patient.Study.ID","ARM","dynamic_class3","pathwayscore")]
ggplot(PRED.lme4,
       aes(y=y-mu_init,x=dum_day,group=dynamic_class3,col=dynamic_class3))+
  geom_line(aes(group=interaction(pathwayscore,dynamic_class3,ARM)),size=2,alpha=1)+#geom_line(aes(y=yi-y_init,group=interaction(pathwayscore,dynamic_class3,Patient.Study.ID)),size=1,alpha=1)+
  geom_ribbon(aes(ymax=upr,ymin=lwr,fill=dynamic_class3),alpha=0.4,col=NA)+
  facet_grid(.~ARM,scales="free") + theme_classic(base_size=20)+
  scale_color_discrete(name="Response",labels=c("Non-responder","Responder"))+
  scale_fill_discrete(name="Response",labels=c("Non-responder","Responder"))+
  scale_x_continuous(name="Day",breaks=0:2,labels=c(0,14,180))+
  ylab(paste0(llll[II])) + theme(aspect.ratio=1)+
  geom_line(aes(y=yi-mu_i_init,x=dum_day,group=Patient.Study.ID),linetype=2)


# ggplot(path_clin_long[Patient.Study.ID=="001-129"],aes(y=ssGSEA,x=as.factor(Day)) )+geom_point()+  facet_grid(.~ARM,scales="free") 
# PRED.lme4[yi==max(yi)]
# pthese
library(ggplot2)
library(scales)
library("ggsci")

# Get the colors with 3 classes
cols <- hue_pal()(2)

pthese[, mu_i_init:=sum(is_day0*yi), by=c("Patient.Study.ID","ARM","dynamic_class3","pathwayscore")]
pthese[Patient.Study.ID=="001-101"]

pthese[,PATH:="Estrogen pathway"]

ggplot(pthese,
       aes(y=y-mu_init,x=dum_day,group=dynamic_class3,col=dynamic_class3))+
  geom_line(aes(group=interaction(pathwayscore,dynamic_class3,ARM)),size=2,alpha=1)+#geom_line(aes(y=yi-y_init,group=interaction(pathwayscore,dynamic_class3,Patient.Study.ID)),size=1,alpha=1)+
  geom_line(linetype=2,size=0.5,aes(y=yi-mu_i_init,x=dum_day, group=interaction(pathwayscore,dynamic_class3,ARM,Patient.Study.ID)),alpha=1)+#geom_line(aes(y=yi-y_init,group=interaction(pathwayscore,dynamic_class3,Patient.Study.ID)),size=1,alpha=1)+
  
  geom_ribbon(aes(ymax=upr,ymin=lwr,group=interaction(pathwayscore,dynamic_class3),fill=dynamic_class3),alpha=0.4,col=NA)+
  facet_grid(PATH~ARM,scales="free") + theme_classic(base_size=20)+  
  scale_color_manual(values = (cols),name="Response",labels=c("Non-responder","Responder"))+
  scale_fill_manual(values = (cols),name="Response",labels=c("Non-responder","Responder"))+
  scale_x_continuous(name="Day",breaks=0:2,labels=c(0,14,180))+
  ylab("Pathway activity") + theme(aspect.ratio=1)+
  theme(
    strip.background = element_blank(),
    strip.text = element_blank()
  )+
  theme(
    axis.title=element_blank(),
    axis.text=element_blank(),
    axis.ticks=element_blank(),
    legend.position = "none"
  )

pthese[,PATH:="JNK"]
pthese[pathwayscore=="BIOCARTA_ERK_PATHWAY",PATH:="ERK"]


ggplot(pthese[PATH=="JNK"],
       aes(y=y-mu_init,x=dum_day,group=dynamic_class3,col=dynamic_class3))+
  geom_line(aes(group=interaction(pathwayscore,dynamic_class3,ARM)),size=2,alpha=1)+#geom_line(aes(y=yi-y_init,group=interaction(pathwayscore,dynamic_class3,Patient.Study.ID)),size=1,alpha=1)+
  geom_line(linetype=2,size=0.5,aes(y=yi-mu_i_init,x=dum_day, group=interaction(pathwayscore,dynamic_class3,ARM,Patient.Study.ID)),alpha=1)+#geom_line(aes(y=yi-y_init,group=interaction(pathwayscore,dynamic_class3,Patient.Study.ID)),size=1,alpha=1)+
  
  geom_ribbon(aes(ymax=upr,ymin=lwr,group=interaction(pathwayscore,dynamic_class3),fill=dynamic_class3),alpha=0.4,col=NA)+
  facet_grid(PATH~ARM,scales="free") + theme_classic(base_size=20)+  
  scale_color_manual(values = (cols),name="Response",labels=c("Non-responder","Responder"))+
  scale_fill_manual(values = (cols),name="Response",labels=c("Non-responder","Responder"))+
  scale_x_continuous(name="Day",breaks=0:2,labels=c(0,14,180))+
  ylab("Pathway activity") + theme(aspect.ratio=1)+
  theme(
    strip.background = element_blank(),
    strip.text = element_blank()
  )+
  theme(
    axis.title=element_blank(),
    axis.text=element_blank(),
    axis.ticks=element_blank(),
    legend.position = "none"
  )

pthese[,Treat:="Letrozole alone"]
pthese[ARM=="B",Treat:="Intermittent high \n dose ribociclib"]
pthese[ARM=="C",Treat:="Continuous low \n dose ribociclib"]
pthese$Treat <- factor(pthese$Treat, levels = c("Letrozole alone", "Intermittent high \n dose ribociclib", "Continuous low \n dose ribociclib"))

ggplot(pthese,
       aes(y=y-mu_init,x=dum_day,group=dynamic_class3,col=pathwayscore))+geom_hline(yintercept=0,linetype=2,col="grey",size=1.25)+
  geom_line(aes(group=interaction(pathwayscore,dynamic_class3,ARM)),size=2,alpha=1)+#geom_line(aes(y=yi-y_init,group=interaction(pathwayscore,dynamic_class3,Patient.Study.ID)),size=1,alpha=1)+
  geom_line(linetype=2,size=0.5,aes(y=yi-mu_i_init,x=dum_day, group=interaction(pathwayscore,dynamic_class3,ARM,Patient.Study.ID)),size=2,alpha=1)+#geom_line(aes(y=yi-y_init,group=interaction(pathwayscore,dynamic_class3,Patient.Study.ID)),size=1,alpha=1)+
  geom_ribbon(aes(ymax=upr,ymin=lwr,group=interaction(pathwayscore,dynamic_class3),fill=pathwayscore),alpha=0.4,col=NA)+
  facet_grid(dynamic_class3~Treat,scales="free") + theme_classic(base_size=20)+  
  # scale_fill_brewer(palette = "Paired",name="Pathway",labels=c("ERK activity","JNK activity"))+
  # scale_color_brewer(palette = "Paired",name="Pathway",labels=c("ERK activity","JNK activity"))+
  scale_color_jco(name="Pathway",labels=c("ERK activity","JNK activity"))+
  scale_fill_jco(name="Pathway",labels=c("ERK activity","JNK activity"))+
  scale_x_continuous(name="Day",breaks=0:2,labels=c(0,14,180))+
  ylab("Pathway activity") + theme(aspect.ratio=1)#+
# theme(
#   strip.background = element_blank(),
#  strip.text = element_blank()
# )

ggplot(pthese[PATH=="ERK"],
       aes(y=y-mu_init,x=dum_day,group=dynamic_class3,col=dynamic_class3))+
  geom_line(aes(group=interaction(pathwayscore,dynamic_class3,ARM)),size=2,alpha=1)+#geom_line(aes(y=yi-y_init,group=interaction(pathwayscore,dynamic_class3,Patient.Study.ID)),size=1,alpha=1)+
  geom_line(linetype=2,size=0.5,aes(y=yi-mu_i_init,x=dum_day, group=interaction(pathwayscore,dynamic_class3,ARM,Patient.Study.ID)),alpha=1)+#geom_line(aes(y=yi-y_init,group=interaction(pathwayscore,dynamic_class3,Patient.Study.ID)),size=1,alpha=1)+
  
  geom_ribbon(aes(ymax=upr,ymin=lwr,group=interaction(pathwayscore,dynamic_class3),fill=dynamic_class3),alpha=0.4,col=NA)+
  facet_grid(PATH~ARM,scales="free") + theme_classic(base_size=20)+  
  scale_color_manual(values = (cols),name="Response",labels=c("Non-responder","Responder"))+
  scale_fill_manual(values = (cols),name="Response",labels=c("Non-responder","Responder"))+
  scale_x_continuous(name="Day",breaks=0:2,labels=c(0,14,180))+
  ylab("Pathway activity") + theme(aspect.ratio=1)+
  theme(
    strip.background = element_blank(),
    strip.text = element_blank()
  )+
  theme(
    axis.title=element_blank(),
    axis.text=element_blank(),
    axis.ticks=element_blank(),
    legend.position = "none"
  )
