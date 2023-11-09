rm(list=ls())
require(mgcv);require(data.table);require(dplyr);require(ggplot2);require(tidyr);require(lme4);require(lmerTest);require(parallel)
library(effects);require(umap);require(Rfast);require(ider);library("dendextend");library(ggdendro);require(ggsci);require(viridis);require("Rdimtools")

load( file=paste0("~/Dropbox/Cancer_pheno_evo/data/FELINE2/PhenotypesAllArmsCohort2/CPMPhenotpyeLandscape_C2ProjectC2RevisednewMacrophagesDC.RData"))#C1umap,u_dat,DAY,cell_types_all,ARMS,Subtype,
load( file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Cohort2Metadata/Cohort2Metadata.RData") #metadd,cohort1metadd,annotation.file,compdataLU

u_dat[,Treatmentlab:="Combination ribociclib"]
u_dat[ARM=="A", Treatmentlab:="Letrozole alone"]

dd1<-u_dat[ARM!="A"]%>%dplyr::select(Cell.ID:file_string,V1,V2)
ggplot( dd1,aes(-V1,V2,col=  dynamic_class3)) +  
  geom_point(size=1.5)+scale_color_npg(name="Tumor response")+
  theme_classic()+
  theme(aspect.ratio=1)#axis.text=element_blank(),axis.title=element_blank(),legend.position

ggplot( u_dat[Day!=180],aes(-V1,V2,col=  dynamic_class3)) +  
  geom_point(size=0.5)+scale_color_npg(name="Tumor response")+
  theme_classic()+
  guides(col=guide_legend(override.aes=list(size=4)))+
  theme(aspect.ratio=1)+#axis.text=element_blank(),axis.title=element_blank(),legend.position
  facet_wrap(~Treatmentlab)+labs(y="Phenotype Umap 2",x="Phenotype Umap 1")


M1 <- lm(V2~dynamic_class3*Treatmentlab , data=u_dat[Day!=180])
M1 <- lm(V2~dynamic_class3 , data=u_dat[ARM!="A"][Day!=180])
summary(M1)

M1r <- lmer(V2~dynamic_class3+(1|Patient.Study.ID) , data=u_dat[ARM!="A"][Day!=180])
summary(M1r)

set.seed(1234)
p1<-princomp ( dd1%>%dplyr::select(V1,V2) )
u1<-umap ( dd1%>%dplyr::select(V1,V2) ,n_components=2)
dd1$Dimension1 <- -p1$scores[,"Comp.1"]
dd1$Dimension2 <- u1$layout[,1]
#plot(Dimension2~V2,data=dd1)
ggplot( dd1,aes(-V1,V2,col=  Dimension2)) +  
  geom_point(size=1.5)+theme_classic()

ggplot( dd1,aes(y=Dimension2, x=dynamic_class3,col=dynamic_class3,fill=dynamic_class3)) +  
  geom_violin(alpha=0.7)+  geom_point(size=1.5)+
  scale_color_npg(name="Tumor response")+
  theme_classic()+
  theme(aspect.ratio=1)#axis.text=element_blank(),axis.title=element_blank(),legend.position

summardd1 <- data.table(dd1%>%group_by(Day,Patient.Study.ID,Sample,ARM,Celltype) %>%summarise(n= n(),V1=mean(V1), V2=mean(V2),Dimension1=mean(Dimension1),Dimension2=mean(Dimension2)  ))
summardd2 <- merge(summardd1,unique(metadd%>%dplyr::select(Patient.Study.ID,prop_change,dynamic_class3))  , by="Patient.Study.ID")
#ggplot( summardd2 , aes(y=V2,x=-V1,col=Dimension1)) +geom_point() 
ggplot( data.table(summardd2[n>100][Day!=180]%>%group_by(Patient.Study.ID,prop_change,dynamic_class3)%>%summarise(Dimension2=weighted.mean(Dimension2,n),n=sum(n)) )[n>0]  ,
        aes(x=-Dimension2, y=(prop_change)) )+
  #geom_hline(yintercept=1,linetype="dashed")+
  geom_smooth(method="gam",formula=y~s(x,k=3),level=0.8,col="black",fill="grey") + 
  #stat_smooth(method="gam",formula=y~s(x,k=3),level=0.8,col="black",fill="grey")+
  geom_point(aes(col=dynamic_class3,size=n)) +
  theme_classic(base_size=22) + theme(aspect.ratio=1)+
  coord_trans(y="exp")+
  labs(y="Tumor size post treatment \n (relative to baseline)",x="Average myeloid M2 differentiation \n in first 14 days of treatment")+
  scale_color_npg(name="Tumor response",labels=c("Resistant","Sensitive"))+
  scale_size_continuous(name="Sample size")

ggsave(filename=paste0("/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cohort 2 results/Macrophage differentiation predicts tumor size response.png"),width=10,height=10)

Cohort1Extract <-function(){
  load(file= "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/UMAP genes macrophage only ALLARMS/UMAP genes macrophage only ALLARMS.RData")
  #load( file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Cohort2Metadata/Cohort2Metadata.RData") #metadd,cohort1metadd,annotation.file,compdataLU
  
  dd1<-u_dat[ARM!="A"]%>%dplyr::select(Cell.ID:PhenoCelltype,V1,V2)
  return(dd1)
}
ddCohort1<-Cohort1Extract()
#p1<-princomp ( ddCohort1%>%dplyr::select(V1,V2) )
#u1 <- umap ( ddCohort1%>%dplyr::select(V1,V2) ,n_components=2)
#ddCohort1$Dimension1 <- -p1$scores[,"Comp.1"]
#ddCohort1$Dimension2 <- u1$layout[,1]
preds <- predict(u1,ddCohort1%>%dplyr::select(V1,V2)  ) 
ddCohort1$Dimension2 <-preds[,1] 
summardd1Cohort1 <- data.table(ddCohort1%>%group_by(Day,Patient.Study.ID,Sample,ARM,Celltype) %>%summarise(n= n(),V1=mean(V1), V2=mean(V2),Dimension2=mean(Dimension2)  ))
summardd2Cohort1 <- merge(summardd1Cohort1,unique(cohort1metadd%>%dplyr::select(Patient.Study.ID,prop_change,dynamic_class3))  , by="Patient.Study.ID")
#ggplot( summardd2Cohort1 , aes(y=V2,x=-V1,col=Dimension1)) +geom_point() 

ggplot( data.table(summardd2Cohort1[n>50][Day!=180]%>%group_by(Patient.Study.ID,prop_change,dynamic_class3)%>%summarise(Dimension1=NA,Dimension2=weighted.mean(Dimension2,n),n=sum(n)) )[n>0]  ,
        aes(x=-Dimension2, y=(prop_change)) )+
  #geom_hline(yintercept=1,linetype="dashed")+
  geom_smooth(method="gam",formula=y~s(x,k=3),level=0.8,col="black",fill="grey") + 
  #stat_smooth(method="gam",formula=y~s(x,k=3),level=0.8,col="black",fill="grey")+
  geom_point(aes(col=dynamic_class3,size=n)) +
  theme_classic(base_size=22) + theme(aspect.ratio=1)+
  coord_trans(y="exp")+
  labs(y="Tumor size post treatment \n (relative to baseline)",x="Average myeloid M2 differentiation \n in first 14 days of treatment")+
  scale_color_npg(name="Tumor response",labels=c("Resistant","Sensitive"))+
  scale_size_continuous(name="Sample size")


combplot<-rbind(data.table(Cohort="Discovery",summardd2Cohort1[n>20][Day!=180]%>%group_by(Patient.Study.ID,prop_change,dynamic_class3)%>%summarise(Dimension1=NA,Dimension2=weighted.mean(Dimension2,n),n=sum(n)) )[n>0] ,
                data.table(Cohort="Validation", summardd2[n>100][Day!=180]%>%group_by(Patient.Study.ID,prop_change,dynamic_class3)%>%summarise(Dimension1=NA,Dimension2=weighted.mean(Dimension2,n),n=sum(n)) )[n>0] 
)
ggplot( combplot, aes(x=-Dimension2, y=(prop_change) ) )+
  #geom_hline(yintercept=1,linetype="dashed")+
  geom_smooth(method="gam",formula=y~s(x,k=3),level=0.8,col="black",fill="grey") + 
  #stat_smooth(method="gam",formula=y~s(x,k=3),level=0.8,col="black",fill="grey")+
  geom_point(aes(col=dynamic_class3,size=n, shape=Cohort)) +
  theme_classic(base_size=10) + theme(aspect.ratio=1)+
  # coord_trans(y="exp")+
  labs(y="Tumor size post treatment \n (relative to baseline)",x="Average myeloid M2 differentiation \n in first 14 days of treatment")+
  labs(y="Post treatment fraction of tumor remaining",x="Average myeloid M2 differentiation \n in first 14 days of treatment")+
  scale_color_npg(name="Tumor response",labels=c("Resistant","Sensitive"))+
  scale_size_continuous(name="Sample size")#+facet_wrap(~Cohort)

ggsave(filename=paste0("/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cohort 2 results/Cohort 1 and 2Macrophage differentiation predicts tumor size response.png"),width=5,height=5)

summary( lm(prop_change~ I(-Dimension2)+Cohort,combplot) )

combplot<-rbind(data.table(Cohort="Discovery",summardd2Cohort1[n>20][Day!=180]%>%group_by(Patient.Study.ID,prop_change,dynamic_class3)%>%summarise(Dimension1=NA,Dimension2=weighted.mean(Dimension2,n),n=sum(n)) )[n>0] ,
                data.table(Cohort="Validation", summardd2[n>20][Day!=180]%>%group_by(Patient.Study.ID,prop_change,dynamic_class3)%>%summarise(Dimension1=NA,Dimension2=weighted.mean(Dimension2,n),n=sum(n)) )[n>0] 
)

combplot
ggplot( combplot, aes(y=-Dimension2, x=(dynamic_class3) ) )+
  geom_violin(aes(group=dynamic_class3,fill=dynamic_class3),alpha=0.25)+
  geom_jitter(aes(col= dynamic_class3, shape= Cohort ), width=0.1,size=2.5) +
  theme_classic(base_size=10) + theme(aspect.ratio=1)+
  # coord_trans(y="exp")+
  labs(y="Tumor size post treatment \n (relative to baseline)",x="Average myeloid M2 differentiation \n in first 14 days of treatment")+
  labs(y="Post treatment fraction of tumor remaining",x="Average myeloid M2 differentiation \n in first 14 days of treatment")+
  scale_color_npg(name="Tumor response",labels=c("Resistant","Sensitive"))+
  scale_size_continuous(name="Sample size")#+facet_wrap(~Cohort)



joineddd <- rbind(
  data.table(Cohort="Discovery" , ddCohort1%>%select(Treatment,Cell.ID,Patient.Study.ID,Day,ARM,dynamic_class3,Celltype,V1,V2) ),
 data.table(Cohort="Vavlidation" , dd1%>%select(Treatment,Cell.ID,Patient.Study.ID,Day,ARM,dynamic_class3,Celltype,V1,V2))
)

summardd1Cohort1 <- data.table(ddCohort1%>%group_by(Day,Patient.Study.ID,Sample,ARM,Celltype) %>%summarise(n= n(),V1=mean(V1), V2=mean(V2),Dimension2=mean(Dimension2)  ))

joineddd%>%group_by(Day,Patient.Study.ID,Sample,ARM,Celltype) %>%summarise(n= n(),V1=mean(V1), V2=mean(V2),Dimension2=mean(Dimension2)  ))
joineddd$V2%>%hist()



dd1<-u_dat%>%dplyr::select(Cell.ID:file_string,V1,V2)



