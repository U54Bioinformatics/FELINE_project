rm(list=ls())
require(data.table);require(dplyr);require(tidyr);require(ggplot2)
require(boot);require("compositions");require(vegan);require(ggsci)
require(mgcv);require(lme4);require(lmerTest);require(parallel)

# macrophage annotations
load( file=paste0("~/Dropbox/Cancer_pheno_evo/data/FELINE2/PhenotypesAllArmsCohort2/CPMPhenotpyeLandscape_C2ProjectC2RevisednewMacrophagesDC.RData"))#C1umap,u_dat,DAY,cell_types_all,ARMS,Subtype,
load( file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Cohort2Metadata/Cohort2Metadata.RData") #metadd,cohort1metadd,annotation.file,compdataLU

dd1 <- u_dat %>% dplyr::select(Cell.ID:file_string,V1,V2)

Cohort1Extract <- function(){
  load(file= "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/UMAP genes macrophage only ALLARMS/UMAP genes macrophage only ALLARMS.RData")
  dd1<-u_dat%>%dplyr::select(Cell.ID:PhenoCelltype,V1,V2)
  return(dd1)
}
ddCohort1 <- Cohort1Extract()

ddCohort2 <- dd1
ddCohort2[,Treatment:="letrozole + ribo"]
ddCohort2[ARM=="A",Treatment:="letrozole"]


joineddd <- rbind(
  data.table(Cohort="Discovery" ,    ddCohort1%>%dplyr::select(Treatment,Cell.ID,Sample,Patient.Study.ID,Day,ARM,dynamic_class3,Celltype,Celltype_subtype,V1,V2) ),
  data.table(Cohort="Validation" ,  ddCohort2%>%dplyr::select(Treatment,Cell.ID,Sample,Patient.Study.ID,Day,ARM,dynamic_class3,Celltype,Celltype_subtype,V1,V2))
)

# all cell type annotations
dd0 <-rbind(
  data.table(Cohort="Validation", metadd %>%dplyr::select(Cell.ID,Celltype,Celltype_subtype,dynamic_class3,Patient.Study.ID,Day,ARM)) , 
  data.table(Cohort="Discovery",cohort1metadd%>%dplyr::select(Cell.ID,Celltype,Celltype_subtype,dynamic_class3,Patient.Study.ID,Day,ARM) ))

dd0[,Treatmentlab:= "Combination ribociclib"]
dd0[ARM=="A",Treatmentlab:= "Letrozole alone"]
dd0[Cell.ID%in% joineddd[Celltype=="Macrophages"][Celltype_subtype=="Macrophages"][V2>0]$Cell.ID, Celltype_subtype:= "M2 macrophages" ]
dd0[Cell.ID%in% joineddd[Celltype=="Macrophages"][Celltype_subtype=="Macrophages"][V2<0]$Cell.ID, Celltype_subtype:= "M1 macrophages" ]
dd0[Celltype=="Macrophages"]$Celltype_subtype%>%table()
dd0[Celltype_subtype=="Endothelial cells",Celltype_subtype:="Vas-Endo"]

dd0$CellAnnot1<-dd0$Celltype
dd0$CellAnnot2<-dd0$Celltype_subtype
dd0[Celltype%in%c("Cancer cells","Normal epithelial cells"),CellAnnot1:="Epithelial cells"]

# Total cell count per sample= sampling effort
dd0[,totalCount:= length(Cell.ID) ,by=c("Patient.Study.ID","Day","Cohort","dynamic_class3","Treatmentlab")]




dd_wbc <- data.table(  read.csv("/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE/Blood work CBC tidy.csv"))
dd_wbc <- data.table(  read.csv("/Users/jason/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Supplementary Data/FELINE blood work.csv"))

dd_wbc[Visit.Name=="EOT Tumor progression"]
dd_wbc[Visit.Name=="EOT Tumor progression",Visit.Name:="EOT Days 8-15 of cycle 6"]
unique(dd_wbc$Visit.Name)
lu_visName<- data.table(Visit.Name=c("Screening","Day 1 of cycle 2 (+/-3 days)",
                                     "Day 1 of cycle 3 (+/-3 days)","Day 1 of cycle 4 (+/-3 days)",
                                     "Day 1 of cycle 5 (+/-3 days)","Day 1 of cycle 6 (+/-3 days)","EOT Days 8-15 of cycle 6"),
                        Day=c(0,  31,  61 , 91 ,121, 151, 180)
)

ddM <- merge(dd_wbc,lu_visName,by="Visit.Name")
ddM[,Hematology.Sample.Collection.Date:=NULL]
ddM[,Data.Entry.Date:=NULL]
ddM[,Visit.Name:=NULL]

ggplot(ddM[][order(Day)], aes(y=log(1+White.Blood.Cells..WBC..Blood.K.UL), x=Day))+geom_path(aes(group=Patient.Study.ID))


ddM_resp <- merge(unique( dd0%>%dplyr::select(Patient.Study.ID,dynamic_class3,Treatmentlab)),
                  ddM,
                  by=c("Patient.Study.ID"))[!is.na( Lymphocytes..Absolute.Lymph.Count..Blood.K.UL)]#,all=TRUE)

ddM_resp[is.na( Neutrophils..Absolute.Neutrophil.Count..Blood.K.UL )]

WBCdata<-ddM_resp[! is.na(dynamic_class3)][order(Patient.Study.ID,Day)]
WBCdata$Patient.Study.ID%>%unique()
save(WBCdata, file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/Figure6/PeripheralWBC.RData") 

g1<- gam(I(White.Blood.Cells..WBC..Blood.K.UL)*1000~s(Day,k=4,by=as.factor(Treatmentlab)), data=WBCdata)
summary(g1)
plot(g2)
g2<- gam(I(log(White.Blood.Cells..WBC..Blood.K.UL*1000))~s(sqrt(Day),k=4,by=as.factor(Treatmentlab)), data=WBCdata)
summary(g2)

m1<- lm(I(log(White.Blood.Cells..WBC..Blood.K.UL*1000))~Treatmentlab, data=WBCdata[Day==0])
summary(m1)

ggplot(ddM_resp[! is.na(dynamic_class3)][order(Day)][], aes(y=(White.Blood.Cells..WBC..Blood.K.UL)*1000, x=Day, col=Treatmentlab))+
  geom_point(position=position_dodge(width=10),size=2.5) +
  geom_smooth(method = "gam",  se = T, formula=y~s(x, k=3),
              method.args = list(family = "poisson"))+
  coord_trans(y="log")+
  theme_classic(base_size=20)+theme(aspect.ratio=1)+
  labs(y="Peripheral white blood cells (per ul)")+
  scale_color_manual(name="Treatment", values=rev(pal_jco("default")(2)))


ggplot(ddM_resp[! is.na(dynamic_class3)][order(Day)][], aes(y=(White.Blood.Cells..WBC..Blood.K.UL), x=sqrt(Day), fill=Treatmentlab,col=Treatmentlab))+
  geom_smooth(method = "gam",  se = T, alpha=0.6, formula=y~s(x, k=3))+
  geom_point()+
  coord_trans(y="log",x = scales::trans_new("square", function(x) x^2, "sqrt"))+
  theme_classic(base_size=26)+theme(aspect.ratio=1)+
  labs(y="Peripheral white blood cells (K/ul)",x="Day")+
  scale_color_manual(name="Treatment", values=rev(pal_jco("default")(2)))+
  scale_fill_manual(name="Treatment", values=rev(pal_jco("default")(2))) +
  scale_y_continuous(breaks=c(2,4,8,16,32),labels=c(2,4,8,16,32))+
  scale_x_continuous(breaks=sqrt(c(0,90,180)),labels=c(0,90,180))
paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Ribo and Letrozole Peripheral white blood cell count by treatment over time.png"),height=10,width=10)

ggplot(ddM_resp[! is.na(dynamic_class3)][order(Day)][], aes(y=(White.Blood.Cells..WBC..Blood.K.UL), x=sqrt(Day), fill=Treatmentlab,col=Treatmentlab))+
  geom_smooth(method = "gam",  se = T, alpha=0.6, formula=y~s(x, k=3))+
  geom_point(size=3.5)+
  coord_trans(y="log",x = scales::trans_new("square", function(x) x^2, "sqrt"))+
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


ggplot(ddM_resp[! is.na(dynamic_class3)][order(Day)][], aes(y=(White.Blood.Cells..WBC..Blood.K.UL), x=sqrt(Day), fill=dynamic_class3,col=dynamic_class3))+
  geom_smooth(method = "gam",  se = T, alpha=0.6, formula=y~s(x, k=3))+
  geom_point()+
  coord_trans(y="log",x = scales::trans_new("square", function(x) x^2, "sqrt"))+
  theme_classic(base_size=26)+theme(aspect.ratio=1)+
  labs(y="Peripheral white blood cells (K/ul)",x="Day")+
  scale_color_npg(name="Response", labels=c("Resistant","Sensitive"))+
  scale_fill_npg(name="Response", labels=c("Resistant","Sensitive")) +
  scale_y_continuous(breaks=c(2,4,8,16,32),labels=c(2,4,8,16,32))+
  scale_x_continuous(breaks=sqrt(c(0,90,180)),labels=c(0,90,180))+
  facet_wrap(~Treatmentlab, ncol=1)
#ggsave(paste0(paperfile,"Ribo and Letrozole Peripheral white blood cell count by response and treatment over time.png"),height=10,width=10)




ggplot(ddM_resp[order(Day)][Monocytes..Absolute.Monocyte.Count..Blood..K.UL<5], aes(y=(Monocytes..Absolute.Monocyte.Count..Blood..K.UL), x=sqrt(Day), fill=dynamic_class3,col=dynamic_class3))+
  geom_point()+
  coord_trans(y="log",x = scales::trans_new("square", function(x) x^2, "sqrt"))+
  theme_classic(base_size=26)+theme(aspect.ratio=1)+
  labs(y="Peripheral white blood cells (K/ul)",x="Day")+
  scale_color_npg(name="Response", labels=c("Resistant","Sensitive"))+
  scale_fill_npg(name="Response", labels=c("Resistant","Sensitive")) +
  scale_y_continuous(breaks=c(2,4,8,16,32),labels=c(2,4,8,16,32))+
  scale_x_continuous(breaks=sqrt(c(0,90,180)),labels=c(0,90,180))+
  facet_wrap(~Treatmentlab, ncol=1)
#ggsave(paste0(paperfile,"Ribo and Letrozole Peripheral white blood cell count by response and treatment over time.png"),height=10,width=10)

boxplot( ddM_resp$Lymphocytes..Absolute.Lymph.Count..Blood.K.UL )

ggplot(ddM_resp[order(Day)][], aes(y=(Lymphocytes..Absolute.Lymph.Count..Blood.K.UL), x=sqrt(Day), fill=Treatmentlab,col=Treatmentlab))+
  geom_smooth(method = "gam",  se = T, alpha=0.6, formula=y~s(x, k=3))+
  geom_point()+
  coord_trans(y="log",x = scales::trans_new("square", function(x) x^2, "sqrt"))+
  theme_classic(base_size=26)+theme(aspect.ratio=1)+
  labs(y="Peripheral lymphocytes (K/ul)",x="Day")+
  scale_color_manual(name="Treatment", values=rev(pal_jco("default")(2)))+
  scale_fill_manual(name="Treatment", values=rev(pal_jco("default")(2))) +
  scale_y_continuous(breaks=c(2,4,8,16,32),labels=c(2,4,8,16,32))+
  scale_x_continuous(breaks=sqrt(c(0,90,180)),labels=c(0,90,180))
paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Ribo and Letrozole Peripheral lymphocyte count by treatment over time.png"),height=10,width=10)

ggplot(ddM_resp[order(Day)][Lymphocytes..Absolute.Lymph.Count..Blood.K.UL<25], aes(y=(Lymphocytes..Absolute.Lymph.Count..Blood.K.UL), x=sqrt(Day), fill=Treatmentlab,col=Treatmentlab))+
  geom_smooth(method = "gam",  se = T, alpha=0.6, formula=y~s(x, k=3))+
  geom_point()+
  coord_trans(y="log",x = scales::trans_new("square", function(x) x^2, "sqrt"))+
  theme_classic(base_size=26)+theme(aspect.ratio=1)+
  labs(y="Peripheral lymphocytes (K/ul)",x="Day")+
  scale_color_manual(name="Treatment", values=rev(pal_jco("default")(2)))+
  scale_fill_manual(name="Treatment", values=rev(pal_jco("default")(2))) +
  scale_y_continuous(breaks=c(2,4,8,16,32),labels=c(2,4,8,16,32))+
  scale_x_continuous(breaks=sqrt(c(0,90,180)),labels=c(0,90,180))
paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Ribo and Letrozole Peripheral lymphocyte count without outliers by treatment over time.png"),height=10,width=10)



ggplot(ddM_resp[order(Day)][], aes(y=(Monocytes..Absolute.Monocyte.Count..Blood..K.UL), x=sqrt(Day), fill=Treatmentlab,col=Treatmentlab))+
  geom_smooth(method = "gam",  se = T, alpha=0.6, formula=y~s(x, k=3))+
  geom_point()+
  coord_trans(y="log",x = scales::trans_new("square", function(x) x^2, "sqrt"))+
  theme_classic(base_size=26)+theme(aspect.ratio=1)+
  labs(y="Peripheral monocytes (K/ul)",x="Day")+
  scale_color_manual(name="Treatment", values=rev(pal_jco("default")(2)))+
  scale_fill_manual(name="Treatment", values=rev(pal_jco("default")(2))) +
  scale_y_continuous(breaks=c(2,4,8,16,32),labels=c(2,4,8,16,32))+
  scale_x_continuous(breaks=sqrt(c(0,90,180)),labels=c(0,90,180))
#ggsave(paste0(paperfile,"Ribo and Letrozole Peripheral monocytes count without outliers by treatment over time.png"),height=10,width=10)


ggplot(ddM_resp[order(Day)][Lymphocytes..Absolute.Lymph.Count..Blood.K.UL<10], aes(y=(Lymphocytes..Absolute.Lymph.Count..Blood.K.UL), x=sqrt(Day), fill=dynamic_class3,col=dynamic_class3))+
  geom_point()+
  coord_trans(y="log", x = scales::trans_new("square", function(x) x^2, "sqrt"))+
  theme_classic(base_size=26)+theme(aspect.ratio=1)+
  labs(y="Peripheral white blood cells (K/ul)",x="Day")+
  scale_color_npg(name="Response", labels=c("Resistant","Sensitive"))+
  scale_fill_npg(name="Response", labels=c("Resistant","Sensitive")) +
  scale_y_continuous(breaks=c(2,4,8,16,32),labels=c(2,4,8,16,32))+
  scale_x_continuous(breaks=sqrt(c(0,90,180)),labels=c(0,90,180))+
  facet_wrap(~Treatmentlab, ncol=1)




ggplot(ddM_resp[! is.na(dynamic_class3)][order(Day)][White.Blood.Cells..WBC..Blood.K.UL<20], aes(y=(White.Blood.Cells..WBC..Blood.K.UL)*1000, x=Day, col=Treatmentlab))+
  geom_point(position=position_dodge(width=10),size=2.5) +
  geom_smooth(method = "gam",  se = T, formula=y~s(x, k=3),
            method.args = list(family = "poisson"))+
  coord_trans(y="log")+
  theme_classic(base_size=20)+
  labs(y="Peripheral white blood cells (per ul)")


ggplot(ddM_resp[! is.na(dynamic_class3)][order(Day)], aes(y=(Lymphocytes..Absolute.Lymph.Count..Blood.K.UL)*1000, x=Day, col=dynamic_class3))+
  geom_point(position=position_dodge(width=10),size=2.5) +
  geom_smooth(method = "gam",  se = T, formula=y~s(x, k=3),
              method.args = list(family = "poisson"))+
  coord_trans(y="log")+
  theme_classic()+
  labs(y="Lymphocytes cells (per ul)")+
  facet_wrap(~Treatmentlab)

ggplot(ddM_resp[! is.na(dynamic_class3)][order(Day)], aes(y=(Monocytes..Absolute.Monocyte.Count..Blood..K.UL)*1000, x=Day, col=Treatmentlab))+
  geom_point(position=position_dodge(width=10),size=2.5) +
  geom_smooth(method = "gam",  se = T, formula=y~s(x, k=3),
              method.args = list(family = "poisson"))+
  coord_trans(y="log")+
  theme_classic()+
  labs(y="Lymphocytes cells (per ul)")+
  facet_wrap(~dynamic_class3)



ggplot(ddM_resp[!is.na(dynamic_class3)][order(Day)], aes(y=log(White.Blood.Cells..WBC..Blood.K.UL), x=Day, col=dynamic_class3))+geom_point()+
  geom_smooth(method="gam", formula = y ~ s(x, k=5))+
  facet_wrap(~Treatmentlab, ncol=1)

ggplot(ddM_resp[!is.na(dynamic_class3)][order(Day)], aes(y=log(Lymphocytes..Absolute.Lymph.Count..Blood.K.UL), x=Day, col=Treatmentlab))+geom_point()+
  geom_smooth(method="gam", formula = y ~ s(x, k=5))

ggplot(ddM_resp[!is.na(dynamic_class3)][order(Day)], aes(y=log(Lymphocytes..Absolute.Lymph.Count..Blood.K.UL), x=Day, col=Treatmentlab))+geom_point()+
  geom_smooth(method="gam", formula = y ~ s(x, k=5))

ggplot(ddM_resp[!is.na(dynamic_class3)][order(Day)], aes(y=log(Lymphocytes..Absolute.Lymph.Count..Blood.K.UL/Monocytes..Absolute.Monocyte.Count..Blood..K.UL), x=Day, col=dynamic_class3))+geom_point()+
  geom_smooth(method="gam", formula = y ~ s(x, k=5))+
  facet_wrap(~Treatmentlab, ncol=1)

ddM_resplong <- data.table( ddM_resp%>%gather(Var, Val, White.Blood.Cells..WBC..Blood.K.UL:Basophils..Absolute.Basophil.Count..Blood.K.UL) )
ddM_resplong[, Count:=1000*Val]
ddM_resplong[,InitialCount:= sum(Count* (Day==0), na.rm=T), by=c("Var","Patient.Study.ID")]
ddM_resplong[, RelativelnCount:=log(Count)-log(InitialCount)]


ggplot(ddM_resplong[! is.na(dynamic_class3)][order(Day)][Var=="Lymphocytes..Absolute.Lymph.Count..Blood.K.UL"][Count<10000], aes(y=(Count), x=Day, col=dynamic_class3))+
  geom_point(position=position_dodge(width=10),size=2.5) +
  geom_boxplot(alpha=0.6,aes(fill=dynamic_class3,group=interaction(Day,dynamic_class3))) +
  facet_wrap(~Treatmentlab)+
    geom_smooth(method = "gam",  se = T, formula= y~s(x, k=4),
                method.args = list(family = "poisson"))+
    coord_trans(y="log")
  
  









