require(mgcv);require(data.table);require(dplyr);require(ggplot2);require(tidyr);require(parallel);require(lme4);require(lmerTest);require(parallel)
library(effects)#require(umap)
rm(list=ls())

cell_types_all <-c("Cancer_cells")
iii<- 1
load( file="/Users/jason/Dropbox/FELINE Project/Data_analysis/scRNA/15_subclone_interaction/ESR1Subclonal ABC.RData")#ESR1Subclonal
ESR1Subclonal<-data.table(ESR1Subclonal)
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
FULL_cell_type_meta_dd <- ESR1Subclonal#merge(cell_type_meta_dd,response_code_dd,by="Sample")        #FULL_cell_type_meta_dd[Celltype_subtype=="Fibroblasts",Celltype_subtype:="CAF-S1"]
print(unique(FULL_cell_type_meta_dd$Celltype_subtype))

# Load adjusted gsea scores
gsea_path<-"~/Dropbox/FELINE Project/Data_analysis/scRNA/05_ssGSEA_score/Signature_c2_hallmark/results/"   #"~/Dropbox/FELINE/Data_share/Modeling_Data/pathway_seperated_files/Data_gene_count_per_celltype_model_zinbwave_ssGSEA/FEL001043/"
cell_type_gsea_adj <-data.table( readRDS(file=paste0(gsea_path,"FEL001046_",cell_type_i,"_scRNA.zinbwave.normalized.ssGSEA_scores.RDS" ))) 
setnames( cell_type_gsea_adj,old="Gene Set",new="Gene_Set") #cell_type_gsea_adj[1:5,1:5]

# lu<- unique(FULL_cell_type_meta_dd%>%dplyr::select(Subclone, Patient.Study.ID, ARM,     Response))
# 
# xx=1
# 
# #tic <- Sys.time()
# output <-
#   rbindlist(lapply(1:nrow(lu),function(xx){
#     # subset metadata for one subclonal population
#     FULL_cell_type_meta_dd_CT <- merge(FULL_cell_type_meta_dd,lu[xx],by=c("Subclone", "Patient.Study.ID", "ARM" ,     "Response"))
#     # subset gsea data for one subclonal population
#     cell_type_gsea_adj_CT<-cell_type_gsea_adj[,c("Gene_Set",FULL_cell_type_meta_dd_CT$Cell.ID),with=FALSE]  
#     
#     # Reshape for analysis :wide->long ->wide other way around (transposed)
#     path_long <- gather(cell_type_gsea_adj_CT,Cell.ID,ssGSEA,-Gene_Set)
#     path_clin_long <- data.table(merge( path_long , FULL_cell_type_meta_dd_CT, by= "Cell.ID") )
#     path_clin_wider <- data.table(spread( path_clin_long , Gene_Set, ssGSEA) )
#     # select gene sets to correlate
#     indep_vars<- path_clin_wider%>%dplyr::select(-c(Cell.ID, Subclone, Patient.Study.ID, ARM, Day,      Response ,     Sub_catagory ,Expression))
#     # correlate, order by coef, add metadata and update name of column for correlation coef
#     C1<-data.table(lu[xx],as.data.table(cor(as.matrix(indep_vars), path_clin_wider$Expression),keep.rownames = T)[order(-abs(V1))])
#     setnames(C1,old="V1",new="ESR1_correlation")
#   return(C1)
#   }))


lu<- unique(FULL_cell_type_meta_dd%>%dplyr::select(Subclone, Patient.Study.ID, ARM,   Day,  Response))

output <-
  rbindlist(lapply(1:nrow(lu),function(xx){
    # subset metadata for one subclonal population
    FULL_cell_type_meta_dd_CT <- merge(FULL_cell_type_meta_dd,lu[xx],by=c("Subclone", "Patient.Study.ID", "ARM" , "Day",    "Response"))
    # subset gsea data for one subclonal population
    cell_type_gsea_adj_CT<-cell_type_gsea_adj[,c("Gene_Set",FULL_cell_type_meta_dd_CT$Cell.ID),with=FALSE]  
    
    # Reshape for analysis :wide->long ->wide other way around (transposed)
    path_long <- gather(cell_type_gsea_adj_CT,Cell.ID,ssGSEA,-Gene_Set)
    path_clin_long <- data.table(merge( path_long , FULL_cell_type_meta_dd_CT, by= "Cell.ID") )
    path_clin_wider <- data.table(spread( path_clin_long , Gene_Set, ssGSEA) )
    # select gene sets to correlate
    indep_vars<- path_clin_wider%>%dplyr::select(-c(Cell.ID, Subclone, Patient.Study.ID, ARM, Day,      Response ,     Sub_catagory ,Expression))
    # correlate, order by coef, add metadata and update name of column for correlation coef
    C1<-data.table(lu[xx],as.data.table(cor(as.matrix(indep_vars), path_clin_wider$Expression),keep.rownames = T)[order(-abs(V1))])
    setnames(C1,old="V1",new="ESR1_correlation")
    return(C1)
  }))
save(output, file="/Users/jason/Dropbox/FELINE Project/Data_analysis/scRNA/15_subclone_interaction/Correlated subclonal ESR1 ABC 2.RData")#ESR1Subclonal
load( file="/Users/jason/Dropbox/FELINE Project/Data_analysis/scRNA/15_subclone_interaction/Correlated subclonal ESR1 ABC 2.RData")#ESR1Subclonal

output[,Ribo:="1"]
output[ARM!="A",Ribo:="0"]

s1<-output[rn%in%c("SMID_BREAST_CANCER_BASAL_UP","CREIGHTON_ENDOCRINE_THERAPY_RESISTANCE_5")][is.finite(ESR1_correlation)&abs(ESR1_correlation)<0.9][order(ESR1_correlation)]
write.csv(s1,file= "/Users/jason/Desktop/SR1 correlation with basal and endocrine resistance phenotypes.csv")


d1<-output[rn%in%c("SMID_BREAST_CANCER_BASAL_UP")][is.finite(ESR1_correlation)&abs(ESR1_correlation)<0.9][order(-ESR1_correlation)]
d1%>%group_by()
m1<-lm(ESR1_correlation~-1+Ribo,d1)
d1$pred<-predict(m1,se.fit=TRUE)$fit
d1$ucl<-predict(m1,se.fit=TRUE)$fit + predict(m1,se.fit=TRUE)$se.fit
d1$lcl<-predict(m1,se.fit=TRUE)$fit - predict(m1,se.fit=TRUE)$se.fit
summary(m1)


x1<-ggplot(d1,
           aes(y=ESR1_correlation,x=ARM,group=ARM))+geom_hline(yintercept=0, linetype="dashed")+geom_violin(aes(col=ARM,fill=ARM))+
  geom_jitter(width=0.05,alpha=0.5)+theme_classic()+  geom_errorbar(aes(ymax=ucl,ymin=lcl))+geom_point(size=3,aes(y=pred))+
  labs(y="Correlation of pathway activity \n with ESR1 expression",x="Treatment" )+
  theme(
    legend.position="none",
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+
  # geom_hline(yintercept=logit(1e-3+1/ 4.3 ))
  theme(axis.title=element_blank(),axis.text=element_blank())
ggsave(x1,filename = "/Users/jason/Desktop/BLANK ESR1 smid basal correl by arm.pdf", width = 9, height = 9)






d1<-output[rn%in%c("CREIGHTON_ENDOCRINE_THERAPY_RESISTANCE_5")][is.finite(ESR1_correlation)&abs(ESR1_correlation)<0.9][order(-ESR1_correlation)]
m1<-lm(ESR1_correlation~-1+Ribo,d1)
d1$pred<-predict(m1,se.fit=TRUE)$fit
d1$ucl<-predict(m1,se.fit=TRUE)$fit + predict(m1,se.fit=TRUE)$se.fit
d1$lcl<-predict(m1,se.fit=TRUE)$fit - predict(m1,se.fit=TRUE)$se.fit
summary(m1)


x2<-ggplot(d1,
           aes(y=ESR1_correlation,x=ARM,group=ARM))+geom_hline(yintercept=0, linetype="dashed")+geom_violin(aes(col=ARM,fill=ARM))+
  geom_jitter(width=0.05,alpha=0.5)+theme_classic()+  geom_errorbar(aes(ymax=ucl,ymin=lcl))+geom_point(size=3,aes(y=pred))+
  labs(y="Correlation of pathway activity \n with ESR1 expression",x="Treatment" )+
  theme(
    legend.position="none",
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+
  # geom_hline(yintercept=logit(1e-3+1/ 4.3 ))
  theme(axis.title=element_blank(),axis.text=element_blank())
ggsave(x2,filename = "/Users/jason/Desktop/BLANK ESR1 CREIGHTON_ENDOCRINE_THERAPY_RESISTANCE_5 correl by arm.pdf", width = 9, height = 9)






ggplot(output[rn%in%c("SMID_BREAST_CANCER_BASAL_UP")][is.finite(ESR1_correlation)&abs(ESR1_correlation)<0.9][order(-ESR1_correlation)],
       aes(y=ESR1_correlation,x=ARM,group=ARM))+geom_hline(yintercept=0, linetype="dashed")+geom_violin(aes(col=ARM,fill=ARM))+geom_boxplot()+
  #facet_wrap(~ARM) +
  theme_classic()+
  theme(
    legend.position="none",
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+
  # geom_hline(yintercept=logit(1e-3+1/ 4.3 ))
  theme(axis.title=element_blank(),axis.text=element_blank())
ggsave(Smid1,filename = "/Users/jason/Desktop/BLANK ESR1 expression vs smid basal by arm.pdf", width = 9, height = 9)


ggplot(output[rn%in%c("SMID_BREAST_CANCER_BASAL_UP","CREIGHTON_ENDOCRINE_THERAPY_RESISTANCE_5")][is.finite(ESR1_correlation)&abs(ESR1_correlation)<0.9][order(-ESR1_correlation)],
       aes(x=ESR1_correlation))+geom_histogram()+facet_wrap(~rn)






fit_mod<-function(DATA,ARM_I){
  dd_i <- DATA[ARM==ARM_I]
  responseMod <- lmer(ssGSEA ~ day_fact + (1+day_fact|Subclone), REML=FALSE,data= dd_i) ### responseMod <- lmer(ssGSEA ~ ln_prop_change*day_fact*Ribo + (1|Patient.Study.ID) + (1|day_fact), data= path_clin_long)
  
  output <- data.frame(summary(responseMod)$ coefficients);  output$Effect <- rownames(output)
  output$Intermediate_Clusters_ML <-"Cancer cells"   ;  output$pathwayscore <- dd_i[1]$Gene_Set
  output <- gather(output, statistic, value,"Estimate":"Pr...t.."); output$metric <- "coefs"
  #anova info
  output2 <-   as.data.table(anova(responseMod)%>%dplyr::select(-c("NumDF","DenDF","Sum Sq","Mean Sq") ));  output2$Effect <- rownames(anova(responseMod))
  output2 <- na.omit(output2)
  output2$Intermediate_Clusters_ML <- "Cancer cells" ;    output2$pathwayscore <- dd_i[1]$Gene_Set
  output2 <- gather(output2, statistic, value,"F value":"Pr(>F)"); output2$metric <- "signifcov"
  # join info
  outputFull <- data.table(rbind(output,output2));   outputFull[,Day:="All"] ;   outputFull[,ARM:=ARM_I]
  return(outputFull)
}

lu<- data.table(unique(FULL_cell_type_meta_dd%>%group_by(Patient.Study.ID)%>%mutate(u=length(unique(Day)))%>%dplyr::select( Patient.Study.ID, ARM,u)))[u>1];  lu[,u:=NULL]
output <- rbindlist(lapply(1:nrow(lu)
                           ,function(xx){
                             # subset metadata for one subclonal population
                             FULL_cell_type_meta_dd_CT <- merge(FULL_cell_type_meta_dd,lu[xx],by=names(lu[xx]))
                             # subset gsea data for one subclonal population
                             cell_type_gsea_adj_CT<-cell_type_gsea_adj[,c("Gene_Set",FULL_cell_type_meta_dd_CT$Cell.ID),with=FALSE]  
                             
                             # Reshape for analysis :wide->long ->wide other way around (transposed)
                             path_long <- gather(cell_type_gsea_adj_CT,Cell.ID,ssGSEA,-Gene_Set)
                             path_clin_long <- data.table(merge( path_long , FULL_cell_type_meta_dd_CT, by= "Cell.ID") )
                             path_clin_wider <- data.table(spread( path_clin_long , Gene_Set, ssGSEA) )
                             # select gene sets to correlate
                             # correlate, order by coef, add metadata and update name of column for correlation coef
                             pthwy<-path_clin_long$Gene_Set%>%unique
                             output_CT <- rbindlist(#mc
                               lapply(1:length(pthwy)    
                                      , function(path_i){
                                        path_long2<- path_clin_long[Gene_Set==pthwy[path_i] ]
                                        pathway_NAME <- path_long2[1,]$Gene_Set
                                        cell_NAME <- "Cancer cells"
                                        
                                        path_long2[,DayLab:="Day 0"]
                                        path_long2[Day==14, DayLab:="Day 14"]
                                        path_long2[Day==180, DayLab:="Day 180"]
                                        path_long2[,day_fact:="0"]
                                        path_long2[Day==14, day_fact:="14"]
                                        path_long2[Day==180, day_fact:="180"]
                                        
                                        mu_ssgsea_t0 <- data.table(path_long2[Day==0]%>%group_by(Patient.Study.ID)%>%summarise(mu_ssGSEA_t0=mean(ssGSEA)))
                                        path_long2<-merge(path_long2,mu_ssgsea_t0,by="Patient.Study.ID",all.x=TRUE)
                                        path_long2[is.na(mu_ssGSEA_t0),mu_ssGSEA_t0:=mean(path_long2[Day==min(Day)]$ssGSEA)]
                                        path_long2[,standardised_ssGSEA:=ssGSEA-mu_ssGSEA_t0]
                                        
                                        
                                        
                                        All14<- tryCatch({
                                          All14<-fit_mod(path_long2,ARM_I=path_long2[1]$ARM)   
                                        },error=function(cond) { data.table(Effect=NA, Intermediate_Clusters_ML=NA, pathwayscore=NA, statistic =NA,      value=NA,    metric=NA,Day=NA ,ARM=NA) })
                                        All14[,Patient.Study.ID:=path_long2[1,]$Patient.Study.ID]
                                        All14[,Response:=path_long2[1,]$Response]
                                        return(All14)
                                      }
                                      #,mc.cores =detectCores()-1
                               ))
                           }))



save(output, file="/Users/jason/Dropbox/FELINE Project/Data_analysis/scRNA/15_subclone_interaction/Trend in patient specific expression.RData")#ESR1Subclonal
load( file="/Users/jason/Dropbox/FELINE Project/Data_analysis/scRNA/15_subclone_interaction/Trend in patient specific expression.RData")#ESR1Subclonal
out_coefs <- spread( output[grep("day_",Effect)][metric=="coefs"] ,statistic,value)
out_coefsraw<- spread( output[grep("day_",Effect)][metric=="coefs"] ,statistic,value)
N_pathways <- length(out_coefs$pathwayscore%>%unique)
out_coefs$adj.pval <- c(unlist(
  lapply(seq(1,nrow(out_coefs) - 1000,by=1000),function(i){
    p.adjust(out_coefs$Pr...t..[i:(i+999)], method ="holm", n = N_pathways)
  })
),p.adjust(out_coefs$Pr...t..[(max(seq(1,nrow(out_coefs)-1000,by=1000)) + 1000):nrow(out_coefs)], method ="holm", n =N_pathways))



out_coefs_sig <- out_coefs[adj.pval<0.05][!grep("_DN",pathwayscore)]
pp<-sort(unique(out_coefs_sig[ARM=="C"]$Patient.Study.ID))
ntot <-length(pp)
out_coefs_sig[Patient.Study.ID==pp[8]][order(-abs(t.value))][1:30]
ARMC_paths_det<-data.table( out_coefs_sig[ARM=="C"][Effect=="day_fact180"] %>%group_by(pathwayscore) %>%summarise(n=length(Estimate),ntot=ntot,med_t=mean(abs(t.value))) )[order(-n,-med_t)][n/ntot>0.3]
ARMC_paths_det[,ARM:="C"]
na.omit(ARMC_paths_det[1:30])


SC_summar <- data.table( out_coefsraw[pathwayscore%in%c("SMID_BREAST_CANCER_BASAL_UP","CREIGHTON_ENDOCRINE_THERAPY_RESISTANCE_5")]%>%
                           dplyr::select(pathwayscore,ARM, Patient.Study.ID  ,    Response,Effect,df,Estimate,Std..Error ,  t.value, Pr...t.. )
)[order(pathwayscore,-t.value)]
write.csv(SC_summar,file="/Users/jason/Desktop/SC_summar.csv")
SC_summar[pathwayscore=="SMID_BREAST_CANCER_BASAL_UP"][ Pr...t..<0.05]
SC_summar[ Pr...t..<0.05]
out_coefs_sig[pathwayscore%in%c("SMID_BREAST_CANCER_BASAL_UP","CREIGHTON_ENDOCRINE_THERAPY_RESISTANCE_5")]

# subset gsea data for one subclonal population
path_long <- gather(cell_type_gsea_adj[Gene_Set=="SMID_BREAST_CANCER_BASAL_UP"][,c("Gene_Set",FULL_cell_type_meta_dd$Cell.ID),with=FALSE],Cell.ID,ssGSEA,-Gene_Set)

path_clin_long <- data.table(merge( path_long , FULL_cell_type_meta_dd, by= "Cell.ID") )
#path_clin_wider <- data.table(spread( path_clin_long , Gene_Set, ssGSEA) )

ggplot(path_clin_long,aes(x=ssGSEA,y=Expression))+geom_point()+facet_wrap(Day~ARM)+  geom_smooth(method="gam",formula=y~s(x,k=4),aes(col=as.factor(ARM)))+theme_classic()
ggplot(path_clin_long,aes(x=ssGSEA,y=Expression))+geom_point()+  geom_smooth(method="gam",formula=y~s(x,k=4),aes(col=as.factor(ARM)))+theme_classic()
ggplot(path_clin_long[Day==180],aes(x=ssGSEA,y=Expression,col=as.factor(Day)  ))+
  geom_point()+facet_wrap(~ARM)+
  geom_smooth(formula=y~s(x,k=4),method.args=list(gamma=5),col="black")+theme_classic()


ggplot(path_clin_long[ARM=="C"],aes(x=ssGSEA,y=Expression  ))+
  geom_point(size=0.3,aes(col=as.factor(Day),fill=as.factor(Day)),pch=21)+#facet_wrap(Subclone~Patient.Study.ID)+
  geom_smooth(formula=y~s(x,k=4),method.args=list(gamma=5),col="black")+theme_classic()

ggplot(path_clin_long[ARM=="C"]%>%mutate(ssGSEAround=round(ssGSEA*200)/200)%>%group_by(ARM,ssGSEAround,Day,Subclone, Patient.Study.ID)%>% dplyr::summarise(mu_ESR1=mean(Expression)), aes(x=ssGSEAround,y=mu_ESR1  ,group= as.factor(Day),col=as.factor(Day),fill=as.factor(Day) ))+
  geom_point(size=1,pch=21)+facet_wrap(~Day)+
  geom_smooth(method="gam",formula=y~s(x,k=4),col="black")+theme_classic()


Smid1<- ggplot(path_clin_long%>%mutate(ssGSEAround=round(ssGSEA*250)/250)%>%group_by(ARM,ssGSEAround,Day,Subclone, Patient.Study.ID)%>% dplyr::summarise(mu_ESR1=mean(Expression)), 
               aes(x=(ssGSEAround-min(ssGSEAround)),y=(mu_ESR1) ,group= as.factor(ARM),fill=as.factor(ARM) ))+
  geom_jitter(size=2.5,pch=21,alpha=1,height=0,width=0.01)+#facet_wrap(~Day)+
  geom_smooth(method="gam",formula=y~s(x,k=4),aes(col=as.factor(ARM)))+theme_classic()+
  labs(y="ESR1 expression",x="Smid Breast cancer basal signature" )+
  theme(
    legend.position="none",
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+
  # geom_hline(yintercept=logit(1e-3+1/ 4.3 ))
  theme(axis.title=element_blank(),axis.text=element_blank())

ggsave(Smid1,filename = "/Users/jason/Desktop/BLANK ESR1 expression vs smid basal by arm.pdf", width = 9, height = 9)






# subset gsea data for one subclonal population
path_long <- gather(cell_type_gsea_adj[Gene_Set=="CREIGHTON_ENDOCRINE_THERAPY_RESISTANCE_5"][,c("Gene_Set",FULL_cell_type_meta_dd$Cell.ID),with=FALSE],Cell.ID,ssGSEA,-Gene_Set)

path_clin_long <- data.table(merge( path_long , FULL_cell_type_meta_dd, by= "Cell.ID") )
#path_clin_wider <- data.table(spread( path_clin_long , Gene_Set, ssGSEA) )

ggplot(path_clin_long,aes(x=ssGSEA,y=Expression))+geom_point()+facet_wrap(Day~ARM)+  geom_smooth(method="gam",formula=y~s(x,k=4),aes(col=as.factor(ARM)))+theme_classic()
ggplot(path_clin_long,aes(x=ssGSEA,y=Expression))+geom_point()+  geom_smooth(method="gam",formula=y~s(x,k=4),aes(col=as.factor(ARM)))+theme_classic()
ggplot(path_clin_long[Day==180],aes(x=ssGSEA,y=Expression,col=as.factor(Day)  ))+
  geom_point()+facet_wrap(~ARM)+
  geom_smooth(formula=y~s(x,k=4),method.args=list(gamma=5),col="black")+theme_classic()


ggplot(path_clin_long[ARM=="C"],aes(x=ssGSEA,y=Expression  ))+
  geom_point(size=0.3,aes(col=as.factor(Day),fill=as.factor(Day)),pch=21)+#facet_wrap(Subclone~Patient.Study.ID)+
  geom_smooth(formula=y~s(x,k=4),method.args=list(gamma=5),col="black")+theme_classic()

ggplot(path_clin_long[ARM=="C"]%>%mutate(ssGSEAround=round(ssGSEA*200)/200)%>%group_by(ARM,ssGSEAround,Day,Subclone, Patient.Study.ID)%>% dplyr::summarise(mu_ESR1=mean(Expression)), aes(x=ssGSEAround,y=mu_ESR1  ,group= as.factor(Day),col=as.factor(Day),fill=as.factor(Day) ))+
  geom_point(size=1,pch=21)+facet_wrap(~Day)+
  geom_smooth(method="gam",formula=y~s(x,k=4),col="black")+theme_classic()


Creighton1<- ggplot(path_clin_long%>%mutate(ssGSEAround=round(ssGSEA*250)/250)%>%group_by(ARM,ssGSEAround,Day,Subclone, Patient.Study.ID)%>% dplyr::summarise(mu_ESR1=mean(Expression)), 
                    aes(x=(ssGSEAround-min(ssGSEAround)),y=(mu_ESR1) ,group= as.factor(ARM),fill=as.factor(ARM) ))+
  geom_jitter(size=2.5,pch=21,alpha=1,height=0,width=0.01)+#facet_wrap(~Day)+
  geom_smooth(method="gam",formula=y~s(x,k=4),aes(col=as.factor(ARM)))+theme_classic()+
  labs(y="ESR1 expression",x="Creighton endocrine therapy resistance 5 signature" )+
  theme(
    legend.position="none",
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+
  # geom_hline(yintercept=logit(1e-3+1/ 4.3 ))
  theme(axis.title=element_blank(),axis.text=element_blank())

ggsave(Creighton1,filename = "/Users/jason/Desktop/BLANK ESR1 expression vs Creighton5 by arm.pdf", width = 9, height = 9)







pp<-sort(unique(out_coefs_sig[ARM=="B"]$Patient.Study.ID))
ntot <-length(pp)
ARMB_paths_det<-data.table( out_coefs_sig[ARM=="B"][Effect=="day_fact180"] %>%group_by(pathwayscore) %>%summarise(n=length(Estimate),ntot=ntot,med_t=mean(abs(t.value))) )[order(-n,-med_t)][n/ntot>0.3]
ARMB_paths_det[,ARM:="B"]
na.omit(ARMB_paths_det[1:30])

pp<-sort(unique(out_coefs_sig[ARM=="A"]$Patient.Study.ID))
ntot <-length(pp)
ARMA_paths_det<-data.table( out_coefs_sig[ARM=="A"][Effect=="day_fact180"] %>%group_by(pathwayscore) %>%summarise(n=length(Estimate),ntot=ntot,med_t=mean(abs(t.value))) )[order(-n,-med_t)][n/ntot>0.3]
ARMA_paths_det[,ARM:="A"]
na.omit(ARMA_paths_det[1:30])

HALLMARK_TNFA_SIGNALING_VIA_NFKB
ZHOU_INFLAMMATORY_RESPONSE_FIMA_UP
TIAN_TNF_SIGNALING_NOT_VIA_NFKB
HALLMARK_INTERFERON_ALPHA_RESPONSE
REACTOME_INTERFERON_ALPHA_BETA_SIGNALING
sort(unique(out_coefs_sig[ARM=="C"]$Patient.Study.ID))
out_coefs_sig[Patient.Study.ID=="001-102"][order(-abs(t.value))][1:30] #esr1 loss , est res, erbb
out_coefs_sig[Patient.Study.ID=="001-118"][order(-abs(t.value))][1:30] #TNFA , immune ERBB, stemness
out_coefs_sig[Patient.Study.ID=="001-132"][order(-abs(t.value))][1:30] #esr loss immune up
out_coefs_sig[Patient.Study.ID=="001-123"][order(-abs(t.value))][1:30]  #immune , EMT
out_coefs_sig[Patient.Study.ID=="001-140"][order(-abs(t.value))][1:30]  #esr1 loss , est res, wnt/EMT, MAPK, ERBB, EGF
out_coefs_sig[Patient.Study.ID=="2972-007-701"][order(-abs(t.value))][1:30] #esr1 loss , immune , EMT
out_coefs_sig[Patient.Study.ID=="001-143"][order(-abs(t.value))][1:30] # EGF TNFA, immune, stemness

sort(unique(out_coefs_sig[ARM=="B"]$Patient.Study.ID))
out_coefs_sig[Patient.Study.ID=="001-101"][order(-abs(t.value))][1:30] # Basal , FGF, ESR loss, immune
out_coefs_sig[Patient.Study.ID=="001-113"][order(-abs(t.value))][1:30] # Basal ,ESR loss, immune
out_coefs_sig[Patient.Study.ID=="001-124"][order(-abs(t.value))][1:30] #esr1 loss , est res
out_coefs_sig[Patient.Study.ID=="001-125"][order(-abs(t.value))][1:30] #imm, est res
out_coefs_sig[Patient.Study.ID=="001-129"][order(-abs(t.value))][1:30] # Basal ,ESR loss, 
out_coefs_sig[Patient.Study.ID=="001-142"][order(-abs(t.value))][1:30] #ESR loss, basal, MAPK, imm

out_coefs_sig[Patient.Study.ID=="001-145"][order(-abs(t.value))][1:30] #basal
out_coefs_sig[Patient.Study.ID=="2972-006-601"][order(-abs(t.value))][1:30] #ESR loss,basal,est res,
out_coefs_sig[Patient.Study.ID=="2972-007-703"][order(-abs(t.value))][1:30] 

out_coefs_sig[pathwayscore=="ABBUD_LIF_SIGNALING_1_DN"]
load(file="/Users/jason/Dropbox/FELINE Project/Data_analysis/scRNA/15_subclone_interaction/Correlated subclonal ESR1 ABC 2.RData")#ESR1Subclonal

output[!grep("_DN", output$Gene_Set)][ESR1_correlation>-0.99][ESR1_correlation< -0.4]$ESR1_correlation%>%hist(breaks=100)
#cell_type_gsea_adj[]
#SMID_BREAST_CANCER_BASAL_UP


output[Day!=0][!grep("_DN", rn)][ARM=="C"][Patient.Study.ID=="001-102"][Subclone=="Cluster3"][ESR1_correlation>-0.99][ESR1_correlation< -0.3][order(-abs(ESR1_correlation))]
output[ARM=="C"][Day==180][xx]

data.table(output[ARM=="C"][Day==180]%>%group_by(rn)%>%dplyr::summarise(ESR1_correlation=median(ESR1_correlation)))[ESR1_correlation< 0][order(-abs(ESR1_correlation))][1:10]
xx<-15
#FULL_cell_type_meta_dd_CT <- merge(FULL_cell_type_meta_dd,lu[ARM=="C"][Patient.Study.ID=="001-102"][Subclone=="Cluster3"][Day==180],by=c("Subclone", "Patient.Study.ID", "ARM" , "Day",    "Response"))
FULL_cell_type_meta_dd_CT <- merge(FULL_cell_type_meta_dd,lu[ARM=="C"][Day==180][xx],by=c("Subclone", "Patient.Study.ID", "ARM" , "Day",    "Response"))
# subset gsea data for one subclonal population
cell_type_gsea_adj_CT <- cell_type_gsea_adj[,c("Gene_Set",FULL_cell_type_meta_dd_CT$Cell.ID),with=FALSE]  

# Reshape for analysis :wide->long ->wide other way around (transposed)
path_long <- gather(cell_type_gsea_adj_CT, Cell.ID,ssGSEA,-Gene_Set)
path_clin_long <- data.table(merge( path_long , FULL_cell_type_meta_dd_CT, by= "Cell.ID") )
path_clin_wider <- data.table(spread( path_clin_long , Gene_Set, ssGSEA) )
merge(output,path_clin_wider[1]%>%dplyr::select(Subclone,Patient.Study.ID),by=c("Subclone","Patient.Study.ID"))[ESR1_correlation>-0.99][ESR1_correlation< 0][order(-abs(ESR1_correlation))][1:15]#[ESR1_correlation< -0.3]
# merge(output[!grep("_DN", output$rn)][ESR1_correlation>-0.99][ESR1_correlation< -0.3],
#       FULL_cell_type_meta_dd[Day==180]%>%dplyr::select(Subclone,Day, Patient.Study.ID, ARM),
#       by=c("Subclone", "Patient.Study.ID", "ARM"  ))

ncells<-data.table(FULL_cell_type_meta_dd%>%group_by(Patient.Study.ID,Subclone,Day,ARM,Response)%>%dplyr::summarise(ncells=length(Cell.ID)))
ncells[Day==180][ncells>20]
out_1<-output[Day==180][ESR1_correlation>-0.99][ESR1_correlation< -0.2 ][!grep("_DN", rn)][order(-abs(ESR1_correlation))]
out_1b<-merge(out_1, ncells[Day==180][ncells>20], by=c("Patient.Study.ID", "Subclone", "Day", "ARM", "Response"))
out_1b$ESR1_correlation%>%hist
out_1b[ARM!="A"][order(ESR1_correlation)][1:30]

out_2<-data.table(out_1b%>%group_by(Subclone,Patient.Study.ID,ARM, Day, Response,ncells)%>%slice(1:10))[order(ARM,Patient.Study.ID,ESR1_correlation)]
out_2[ARM!="A"]
out_1_sum <- data.table(out_1%>%group_by(Subclone,Patient.Study.ID)%>%dplyr::summarise(med_ESR1_correlation=median(ESR1_correlation)))
median(out_1_sum$med_ESR1_correlation)
data.table(out_2[Day!=0][!grep("_DN", rn)][ARM!="A"][ESR1_correlation>-0.99][ESR1_correlation< -0.3][order(-abs(ESR1_correlation))]%>%group_by(rn)%>%summarise(n=length(rn)))[order(-n)][][1:25]

data.table(output[Day!=0][!grep("_DN", rn)][ARM=="A"][ESR1_correlation>-0.99][ESR1_correlation< -0.3][order(-abs(ESR1_correlation))]%>%group_by(rn)%>%summarise(n=length(rn)))[order(-n)][][1:25]
data.table(output[Day!=0][!grep("_DN", rn)][ARM!="A"][ESR1_correlation>-0.99][ESR1_correlation< -0.3][order(-abs(ESR1_correlation))]%>%group_by(rn)%>%summarise(n=length(rn)))[order(-n)][][1:25]
data.table(output[Day!=0][!grep("_DN", rn)][ARM=="B"][ESR1_correlation>-0.99][ESR1_correlation< -0.3][order(-abs(ESR1_correlation))]%>%group_by(rn)%>%summarise(n=length(rn)))[order(-n)][][1:25]
data.table(output[Day!=0][!grep("_DN", rn)][ARM=="C"][ESR1_correlation>-0.99][ESR1_correlation< -0.3][order(-abs(ESR1_correlation))]%>%group_by(rn)%>%summarise(n=length(rn)))[order(-n)][][1:25]

out_1b%>%group_by(ARM)%>%summarise(length(unique(paste(Patient.Study.ID,Subclone) )))
#HALLMARK_G2M_CHECKPOINT

ggplot(path_clin_wider,aes(SMID_BREAST_CANCER_BASAL_UP,Expression))+geom_point()+geom_smooth(method="lm")
ggplot(path_clin_wider,aes(HALLMARK_APOPTOSIS,Expression))+geom_point()+geom_smooth(method="lm")
ggplot(path_clin_wider,aes(KEGG_P53_SIGNALING_PATHWAY,Expression))+geom_point()+geom_smooth(method="lm")
ggplot(path_clin_wider,aes(BROWNE_INTERFERON_RESPONSIVE_GENES,Expression))+geom_point()+geom_smooth(method="lm")
#ggplot( ttt[Day==180][ARM=="C"][Patient.Study.ID=="001-102"],aes(BIOCARTA_CELLCYCLE_PATHWAY,ESR1,col=HALLMARK_INTERFERON_GAMMA_RESPONSE))+geom_point()+facet_wrap(~Subclone)

# cell_type_gsea_adj$Gene_Set[grep("HALLMARK_",cell_type_gsea_adj$Gene_Set)]


rm(list=c("path_clin_long","path_clin_wider","path_long","cell_type_gsea_adj_CT","cell_type_gsea_adj"))


rm(list=ls())
#basal receptor genes
basal_genelist <- c("ADM", "ARTN","CCL18","CCL2","CCL20","CCL5","CCL7","CCL8","CCN2","CX3CL1","CXCL1","CXCL10","CXCL11","CXCL5","CXCL8","DEFB1","DKK1",
                    "EDN1","EDN2","FGF9","GAL","HBEGF","IL32","LEFTY2","LTBP1","MIA","NRTN","PDGFRA","REG1A","S100A6","SPP1","TNFRSF11B","VEGFA",
                    "WNT1","WNT5B","WNT6","KRT5","KRT17","FZD6",
                    "TGFB1","TGFB2",
                    "INHBB","FGF13","FGF12","GDF15","BMP7","GMFB",#CREIGHTON_ENDOCRINE_THERAPY_RESISTANCE_5
                    "BDNF","MIA","TNFSF15","CXCL15","CCN2","CCN3","LIF")





require(mgcv);require(data.table);require(dplyr);require(ggplot2);require(tidyr);require(parallel);require(lme4);require(lmerTest);require(parallel)
library(effects)#require(umap)
load( file="/Users/jason/Dropbox/FELINE Project/Data_analysis/scRNA/15_subclone_interaction/Correlated subclonal ESR1 ABC 2.RData")#ESR1Subclonal
library(TSP)
require(mclust)
#ttt<-data.table(readRDS(file="~/Dropbox/FELINE Project/Data_analysis/scRNA/15_subclone_interaction/FEL011046_data_gene_pathway.RDS"))
ttt<-data.table(readRDS(file="~/Dropbox/FELINE Project/Data_analysis/scRNA/15_subclone_interaction/FEL011046_data_gene_pathway.v3.RDS"))
#ttt2<-data.table(gather(ttt,Variable,Expression,ERBB4:YOSHIMURA_MAPK8_TARGETS_UP))
ttt2<-data.table(gather(ttt,Variable,Expression, -c(Cell.ID, Patient.Study.ID, ARM, Day , Response, Subclone, Phase)))
gene_pathway_lu<-data.table(read.csv(file="~/Dropbox/FELINE Project/Data_analysis/scRNA/15_subclone_interaction/gene and pathway catag AB.csv")%>%dplyr::select(-row.ID))
ttt3 <- merge(ttt2,gene_pathway_lu, by="Variable",all.x=TRUE)

#possible transformation  
ttt3[,sign_nn := 1* (Expression>0)*1 + -1*( (Expression<=0)*1)]
#require("VGAM")#lamW") #muhat_i=100  #ttt3[,Expression2:=muhat_i * exp(1 + lambertW((-muhat_i + (Expression/sign_nn)^2/2      )/(muhat_i* exp(1) )))]
ttt3[,Expression2:=Expression]
ttt3[VariableType=="pathway",Expression2:=Expression]


rm(list="ttt2")
Basal_gene_list <- read.csv(file="~/Desktop/Basal genes.csv")
Basal_gene_list$cytokine_growthfactor_gene
Basal_gene_list$SMID_BREAST_CANCER_BASAL_UP
basal_genelist


output[Day==180][ARM=="C"][Patient.Study.ID=="001-102"]
ggplot( ttt[Day==180][ARM=="C"][Patient.Study.ID=="001-102"],aes(HALLMARK_ESTROGEN_RESPONSE_LATE,HALLMARK_ANDROGEN_RESPONSE))+geom_point()+facet_wrap(~Subclone)
ggplot( ttt[Day==180][ARM=="C"][Patient.Study.ID=="001-102"],aes(HALLMARK_ESTROGEN_RESPONSE_LATE,CCL28))+geom_point()+facet_wrap(~Subclone)
require("VGAM")#lamW")
require(umap)
signal_recept_tranUnSlice <- (data.table( ttt3[Day==180][Variable%in%c("ESR1",basal_genelist,as.character(Basal_gene_list$SMID_BREAST_CANCER_BASAL_UP,Basal_gene_list$CREIGHTON_ENDOCRINE_THERAPY_RESISTANCE_5))][!ARM%in%c("A")]%>%
                                            dplyr::select(-c(VariableType,Negative.regulator,Sub_catagory,Expression2,sign_nn,Major_catagory))%>%spread(Variable,Expression)
))
signal_recept_tranSlice <- na.omit(data.table(signal_recept_tranUnSlice %>% group_by(Patient.Study.ID,ARM, Day,Subclone,Response) %>% dplyr::slice(1:100)))

in_dd<- signal_recept_tranSlice%>%dplyr::select(-c(Cell.ID, Patient.Study.ID, ARM, Day,  Response, Subclone, Phase))
set.seed(1234)
umap_md3 <- umap(in_dd, n_neighbors=25, n_components=5)
#umap_md3 <- umap(signal_recept_tranSlice%>%dplyr::select(-c(Cell.ID, Patient.Study.ID, ARM, Day,  Response, Subclone, Phase)), n_neighbors=10, n_components=2)
ggplot(   data.table(signal_recept_tranSlice,umap_md3$layout) , aes(V1,V2,col=Subclone)) + geom_point(size=2)+facet_wrap(~Patient.Study.ID)
ggplot(   data.table(signal_recept_tranSlice,umap_md3$layout) , aes(V1,V2,col=Response,shape=Subclone))+geom_point(size=2)+theme_classic()
ggplot(   data.table(signal_recept_tranSlice,umap_md3$layout) , aes(V1,V2,col=Patient.Study.ID,shape=Subclone))+geom_point(size=2)+theme_classic()
ggplot(   data.table(signal_recept_tranSlice,umap_md3$layout) , aes(V1,V2,col=Day,shape=Subclone))+geom_point(size=2)+theme_classic()


u_datFULL <- data.table(signal_recept_tranUnSlice,
                        predict(umap_md3,
                                data=signal_recept_tranUnSlice %>%
                                  dplyr::select(-c(Cell.ID, Patient.Study.ID, ARM, Day,  Response, Subclone, Phase))   )   )



ggplot(   u_datFULL, aes(V1,V2,col=Response,shape=Subclone))+geom_point(size=2)+theme_classic()
ggplot(   u_datFULL, aes(V1,V2,col=log(CD44-min(CD44)) ))+geom_point(size=2)+theme_classic()


signal_recept_tranUnSliceAllt <- (data.table( ttt3[ARM%in%c("A","B","C")][Variable%in%c("ESR1",basal_genelist,as.character(Basal_gene_list$SMID_BREAST_CANCER_BASAL_UP))]%>%
                                                dplyr::select(-c(VariableType,Negative.regulator,Sub_catagory,Expression2,sign_nn,Major_catagory))%>%spread(Variable,Expression)
))

project_this<- signal_recept_tranUnSliceAllt%>%dplyr::select(-c(Cell.ID, Patient.Study.ID, ARM, Day,  Response, Subclone, Phase))

rm(list=c("ttt","ttt2"))

u_datFULL_ALLtime <- data.table(signal_recept_tranUnSliceAllt, predict(umap_md3, data=project_this  )   )

corV1<- as.data.table(cor(u_datFULL%>%dplyr::select(-c(Cell.ID, Patient.Study.ID, ARM, Day,  Response, Subclone, Phase,V1:V5)),
                          u_datFULL%>%dplyr::select(V1:V5) ),keep.rownames=TRUE)
corV1[order(-abs(V1))][1:10]
corV1[order(-abs(V2))][1:10]
pcaOverall<- princomp(u_datFULL_ALLtime%>%dplyr::select(V1,V5))
u_datFULL_ALLtimeDims<-data.table(u_datFULL_ALLtime, dimension1=pcaOverall$scores[,1],dimension2=pcaOverall$scores[,2])

ggplot(   u_datFULL_ALLtimeDims, aes(V1, V2, col=Response, shape=Subclone))+geom_point(size=2)+theme_classic()+facet_wrap(Response~Day)
corV2<- as.data.table(cor(u_datFULL_ALLtimeDims%>%dplyr::select(-c(Cell.ID, Patient.Study.ID, ARM, Day,  Response, Subclone, Phase,V1:V5,dimension1:dimension2)),
                          u_datFULL_ALLtimeDims%>%dplyr::select(dimension1:dimension2) ),keep.rownames=TRUE)
corV2[order(-abs(dimension1))][1:10]
corV2[order(-abs(dimension2))][1:10]
ggplot(  u_datFULL_ALLtimeDims , aes(y=ESR1,x=as.factor(Day),group=Day))+geom_violin()+facet_wrap(ARM~Response)+theme_classic()
ggplot(  u_datFULL_ALLtimeDims , aes(y=TGFB2,x=as.factor(Day),group=Day))+geom_violin()+facet_wrap(ARM~Response)+theme_classic()
ggplot(  u_datFULL_ALLtimeDims , aes(y=VEGFA,x=as.factor(Day),group=Day))+geom_violin()+facet_wrap(ARM~Response)+theme_classic()


save(u_datFULL_ALLtimeDims, file="/Users/jason/Dropbox/FELINE Project/Data_analysis/scRNA/15_subclone_interaction/Dim reduced basal.RData")

fit_mod

fit_mod<-function(DATA,ARM_I,DAY){
  dd_i <- DATA[ARM==ARM_I&Day%in%c(0,DAY)]
  responseMod <- lmer(Expression ~ -1+Response*day_fact + (1+day_fact|Patient.Study.ID), REML=FALSE,data= dd_i) ### responseMod <- lmer(ssGSEA ~ ln_prop_change*day_fact*Ribo + (1|Patient.Study.ID) + (1|day_fact), data= path_clin_long)
  output <- data.frame(summary(responseMod)$ coefficients);  output$Effect <- rownames(output)
  output$Subclone <- dd_i[1]$Subclone ;   output$Patient.Study.ID <- dd_i[1]$Patient.Study.ID ;   output$Variable <- dd_i[1]$Variable
  output <- gather(output, statistic, value,"Estimate":"Pr...t..")
  output$metric <- "coefs"
  #anova info
  output2 <-   as.data.table(anova(responseMod)%>%dplyr::select(-c("NumDF","DenDF","Sum Sq","Mean Sq") ));  output2$Effect <- rownames(anova(responseMod))
  output2 <- na.omit(output2)
  output2$Subclone <- dd_i[1]$Subclone ;   output2$Patient.Study.ID <- dd_i[1]$Patient.Study.ID ;   output2$Variable <- dd_i[1]$Variable
  output2 <- gather(output2, statistic, value,"F value":"Pr(>F)")
  output2$metric <- "signifcov"
  # join info
  outputFull <- data.table(rbind(output,output2))
  outputFull[,Day:=DAY]
  outputFull[,ARM:=ARM_I]
  return(outputFull)
}

todo<-unique(c("ESR1",basal_genelist,as.character(Basal_gene_list$SMID_BREAST_CANCER_BASAL_UP,Basal_gene_list$CREIGHTON_ENDOCRINE_THERAPY_RESISTANCE_5)))
ttt3[Variable=="ESR1"]
output_CT <- rbindlist(#mc
  lapply(todo
         #1:length(todo)
         , function(path){
           path_clin_long <-  ttt3[Variable==path]
           pathway_NAME <- path  ;   cell_NAME <- "cancer"
           path_clin_long[,DayLab:="Day 0"];   path_clin_long[Day==14, DayLab:="Day 14"];    path_clin_long[Day==180, DayLab:="Day 180"]
           path_clin_long[,day_fact:=as.factor(Day)]
 
           All14<- tryCatch({
             A14<-fit_mod(path_clin_long,ARM_I="A",DAY=14)   
             B14<-fit_mod(path_clin_long,ARM_I="B",DAY=14)   
             C14<-fit_mod(path_clin_long,ARM_I="C",DAY=14)   
             A180<-fit_mod(path_clin_long,ARM_I="A",DAY=180)   
             B180<-fit_mod(path_clin_long,ARM_I="B",DAY=180)   
             C180<-fit_mod(path_clin_long,ARM_I="C",DAY=180)  
             
             All14<-rbind(A14,B14,C14,A180,B180,C180)
           },error=function(cond) { data.table(Effect=NA,Subclone=NA, Patient.Study.ID=NA, Variable=NA, statistic =NA,      value=NA,    metric=NA,Day=NA ,ARM=NA) })
           return(All14)
         }
  ))

names(u_datFULL_ALLtimeDims)[names(u_datFULL_ALLtimeDims)%in%as.character(Basal_gene_list$cytokine_growthfactor_gene)]
names(u_datFULL_ALLtimeDims)[names(u_datFULL_ALLtimeDims)%in%as.character(Basal_gene_list$SMID_BREAST_CANCER_BASAL_UP)]
names(u_datFULL_ALLtimeDims)[names(u_datFULL_ALLtimeDims)%in%as.character(Basal_gene_list$CREIGHTON_ENDOCRINE_THERAPY_RESISTANCE_5)]

ggplot(  u_datFULL_ALLtimeDims , aes(y=(FGF12),x=as.factor(Day),group=Day))+geom_violin()+facet_wrap(ARM~Response)+theme_classic()


C_A<-cor(u_datFULL_ALLtimeDims$ESR1,u_datFULL_ALLtimeDims%>%dplyr::select(one_of(todo)) )
as.data.table(t(C_A),keep.rownames = T)[V1<0][order(V1)]
ggplot(  u_datFULL_ALLtimeDims , aes(y=((FZD6)),x=as.factor(Day),group=Day))+geom_violin()+facet_wrap(ARM~Response)+theme_classic()

ggplot(  u_datFULL_ALLtimeDims , aes(y=(FGF13),x=as.factor(Day),group=Day))+geom_violin()+facet_wrap(ARM~Response)+theme_classic()
ggplot(  u_datFULL_ALLtimeDims , aes(y=(FGF12),x=as.factor(Day),group=Day))+geom_violin()+facet_wrap(ARM~Response)+theme_classic()
ggplot(  u_datFULL_ALLtimeDims , aes(y=((INHBB)),x=as.factor(Day),group=Day))+geom_violin()+facet_wrap(ARM~Response)+theme_classic()
ggplot(  u_datFULL_ALLtimeDims , aes(y=((LTBP1)),x=as.factor(Day),group=Day))+geom_violin()+facet_wrap(ARM~Response)+theme_classic()
ggplot(  u_datFULL_ALLtimeDims , aes(y=((TGFB1)),x=as.factor(Day),group=Day))+geom_violin()+facet_wrap(ARM~Response)+theme_classic()
ggplot(  u_datFULL_ALLtimeDims , aes(y=((TGFB2)),x=as.factor(Day),group=Day))+geom_violin()+facet_wrap(ARM~.)+theme_classic()
ggplot(  u_datFULL_ALLtimeDims , aes(y=((VEGFA)),x=as.factor(Day),group=Day))+geom_violin()+facet_wrap(ARM~.)+theme_classic()
ggplot(  u_datFULL_ALLtimeDims , aes(y=((CDKN2A)),x=as.factor(Day),group=Day))+geom_violin()+facet_wrap(ARM~.)+theme_classic()
ggplot(  u_datFULL_ALLtimeDims , aes(y=((IFRD1)),x=as.factor(Day),group=Day))+geom_violin()+facet_wrap(ARM~.)+theme_classic()


ggplot(  u_datFULL_ALLtimeDims , aes(y=log(1+-min(COL2A1)+COL2A1),x=as.factor(Day),group=Day))+geom_violin()+facet_wrap(ARM~Response)+theme_classic()
ggplot(  u_datFULL_ALLtimeDims , aes(y=log(1+-min(COL2A1)+COL2A1),x=as.factor(Day),group=Day))+geom_violin()+facet_wrap(ARM~Response)+theme_classic()


ggplot(  u_datFULL_ALLtimeDims[TGFB2<6] , aes(y=(TGFB2),x=as.factor(Day),group=Day))+geom_violin()+facet_wrap(ARM~Response)+theme_classic()

ggplot(  u_datFULL_ALLtimeDims , aes(y=log(1+-min(RRAS2)+RRAS2),x=as.factor(Day),group=Day))+geom_violin()+facet_wrap(ARM~Response)+theme_classic()
ggplot(  u_datFULL_ALLtimeDims , aes(y=log(1+-min(RASA2)+RASA2),x=as.factor(Day),group=Day))+geom_violin()+facet_wrap(ARM~Response)+theme_classic()
ggplot(  u_datFULL_ALLtimeDims , aes(y=ESR1,x=as.factor(Day),fill=Response,group=interaction(Day,Response)))+geom_violin()+facet_wrap(Response~ARM)+theme_classic()
ggplot(  u_datFULL_ALLtimeDims , aes(y=INHBB,x=as.factor(Day),fill=Response,group=interaction(Day,Response)))+geom_violin()+facet_wrap(Response~ARM)+theme_classic()
ggplot(  u_datFULL_ALLtimeDims , aes(y=log(1+-min(SLC7A5)+SLC7A5),x=as.factor(Day),group=Day))+geom_violin()+facet_wrap(ARM~Response)+theme_classic()

ggplot(  u_datFULL_ALLtimeDims , aes(y=log(1+-min(TCF7L1)+TCF7L1),x=as.factor(Day),group=Day))+geom_violin()+facet_wrap(ARM~Response)+theme_classic()
ggplot(  u_datFULL_ALLtimeDims , aes(y=log(1+-min(CDKN2A)+CDKN2A),x=as.factor(Day),group=Day))+geom_violin()+facet_wrap(ARM~Response)+theme_classic()
FZD6
ggplot(  u_datFULL_ALLtimeDims , aes(y=log(1+-min(ITPKB)+ITPKB),x=as.factor(Day),group=Day))+geom_violin()+facet_wrap(ARM~Response)+theme_classic()
SLC7A5
#save(output_CT,ttt3, file="/Users/jason/Dropbox/FELINE Project/Data_analysis/scRNA/15_subclone_interaction/temporal trend basal genes.RData")
output_CT[,Subclone:=NULL]
output_CT[][metric=="coefs"&statistic=="Pr...t.."][value<0.05][order(value)]
output_CT[][metric=="coefs"&statistic=="Pr...t.."][value<0.05][order(value)][Day==180]
output_CT[][metric=="coefs"&statistic=="Pr...t.."][value<0.05][order(value)][Day==180][ARM!="A"]
output_CT[][metric=="coefs"&statistic=="Pr...t.."][value<0.05][order(value)][Day==180][ARM=="C"][grep("day_fact",Effect)]
output_CT[][metric=="coefs"&statistic=="Pr...t.."][value<0.05][order(value)][Day==180][ARM=="B"][grep("day_fact",Effect)]
output_CT[][metric=="coefs"&statistic=="Pr...t.."][value<0.05][order(value)][Day==180][ARM!="A"][ARM=="C"]

data.table(output_CT[][metric=="coefs"&statistic=="Pr...t.."][value<0.01][order(value)][Day==180][ARM!="A"][grep("day_fact",Effect)]%>%group_by(Variable)%>%dplyr::summarise(n=length(Effect)))[order(-n)]

output_CT[][metric=="coefs"&statistic=="Pr...t.."][order(value)][Day==180][grep("day_fact",Effect)][Variable=="CDKN2A"]

ggplot(  u_datFULL_ALLtimeDims , aes(ESR1,TGFB2,col=log(ESR1-min(ESR1))))+geom_point(size=0.5)+facet_wrap(ARM~Day)+theme_classic()

ggplot(  u_datFULL_ALLtimeDims , aes(V1,V2,col=Patient.Study.ID,shape=Subclone))+geom_point(size=2)+facet_wrap(Response~Day)+theme_classic()
ggplot(  u_datFULL_ALLtimeDims , aes(dimension1,dimension2,col=Response))+geom_point(size=0.5)+facet_wrap(ARM~Day)+theme_classic()
KRT15
ggplot(  u_datFULL_ALLtimeDims , aes(V1,V2,col=log(PTPN14-min(PTPN14))))+geom_point(size=0.5)+facet_wrap(ARM~Day)+theme_classic()
ggplot(  u_datFULL_ALLtimeDims , aes(V1,V2,col=log(ESR1-min(ESR1))))+geom_point(size=0.5)+facet_wrap(ARM~Day)+theme_classic()

require(mclust)
u_datFULL_ALLtimeDims$classification<- NULL





cell_types_all <-c("Cancer_cells")
iii<- 1
load( file="/Users/jason/Dropbox/FELINE Project/Data_analysis/scRNA/15_subclone_interaction/ESR1Subclonal ABC.RData")#ESR1Subclonal
ESR1Subclonal<-data.table(ESR1Subclonal)
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
FULL_cell_type_meta_dd <- ESR1Subclonal#merge(cell_type_meta_dd,response_code_dd,by="Sample")       
print(unique(FULL_cell_type_meta_dd$Celltype_subtype))

# Load adjusted gsea scores
rm(list=c("ESR1Subclonal","cell_type_meta_dd"))
gene_path<-"~/Dropbox/FELINE Project/Data_analysis/scRNA/02_count/FEL001046/"
# Load adjusted gsea scores
cell_type_gsea_adj <-data.table( readRDS(file=paste0(gene_path,"FEL001046_",cell_type_i,"_scRNA.zinbwave.normalized.RDS" ))) 
setnames( cell_type_gsea_adj,old="Gene.ID",new="Gene_Set") #cell_type_gsea_adj[1:5,1:5]
cell_type_gsea_adj$Gene_Set[grep("ESR1",cell_type_gsea_adj$Gene_Set)]


tic <- Sys.time()
output <-
  rbindlist(lapply(unique(FULL_cell_type_meta_dd$Celltype),function(CELL_Subtype){
    
    FULL_cell_type_meta_dd_CT <- FULL_cell_type_meta_dd[Celltype==CELL_Subtype]
    cell_type_gsea_adj_CT<-cell_type_gsea_adj[,c("Gene_Set",FULL_cell_type_meta_dd_CT$Cell.ID),with=FALSE]  
    output_CT <- rbindlist(#mc
      lapply(#grep("HALLMARK_",cell_type_gsea_adj$Gene_Set) 
        1:nrow(cell_type_gsea_adj)
        , function(path_i,TYPE){
          
          path_long <- gather(cell_type_gsea_adj_CT[path_i,],Cell.ID,ssGSEA,-Gene_Set)
          path_clin_long <- data.table(merge( path_long , FULL_cell_type_meta_dd_CT, by= "Cell.ID") )
          pathway_NAME <- path_clin_long[1,]$Gene_Set
          cell_NAME <- path_clin_long[1,]$Celltype_subtype
          
          path_clin_long[,DayLab:="Day 0"]
          path_clin_long[Day==14, DayLab:="Day 14"]
          path_clin_long[Day==180, DayLab:="Day 180"]
          
          mu_ssgsea_t0<-data.table(path_clin_long[Day==0]%>%group_by(Patient.Study.ID)%>%summarise(mu_ssGSEA_t0=mean(ssGSEA)))
          path_clin_long<-merge(path_clin_long,mu_ssgsea_t0,by="Patient.Study.ID",all.x=TRUE)
          path_clin_long[is.na(mu_ssGSEA_t0),mu_ssGSEA_t0:=mean(path_clin_long[Day==0]$ssGSEA)]
          path_clin_long[,standardised_ssGSEA:=ssGSEA-mu_ssGSEA_t0]
          path_clin_long[,ln_prop_change:=log(1+prop_change)]
          path_clin_long[,exp_rgr_A:=exp(91*rgr_A)]
          path_clin_long[,exp_rgr_B:=exp(89*rgr_B)]
          
          path_clin_long[,class_rgr_A:="Response"]
          path_clin_long[exp_rgr_A>0.8,class_rgr_A:="Non_response"]
          path_clin_long[,ssGSEA_tran:=log(ssGSEA+10)]
          
          All14<- tryCatch({
            A14<-fit_mod(path_clin_long,ARM_I="A",DAY=14)   
            B14<-fit_mod(path_clin_long,ARM_I="B",DAY=14)   
            C14<-fit_mod(path_clin_long,ARM_I="C",DAY=14)   
            
            A180<-fit_mod(path_clin_long,ARM_I="A",DAY=180)   
            B180<-fit_mod(path_clin_long,ARM_I="B",DAY=180)   
            C180<-fit_mod(path_clin_long,ARM_I="C",DAY=180)   
            
            All<-rbind(A14,B14,C14,A180,B180,C180)
          },error=function(cond) { data.table(Effect=NA, Intermediate_Clusters_ML=NA, pathwayscore=NA, statistic =NA,      value=NA,    metric=NA,Day=NA ,ARM=NA) })
          return(All)
        },TYPE=CELL_Subtype
        #,mc.cores =detectCores()-1
      ))
  }))
toc <- Sys.time()
execTime<- toc-tic



