rm(list=ls())
require(data.table)
require(dplyr)
require(tidyr)
require(ggplot2)
require(boot)
require("compositions")
require(vegan)
require(ggsci)
#library(microbiome)
#library(phyloseq)

## Start by getting the M1-M2 phenotypes of myeloid cells
require(mgcv);require(data.table);require(dplyr);require(ggplot2);require(tidyr);require(lme4);require(lmerTest);require(parallel)
library(effects);require(umap);require(Rfast);require(ider);library("dendextend");library(ggdendro);require(ggsci);require(viridis);require("Rdimtools")
library(caret)
library(pROC)

load( file=paste0("~/Dropbox/Cancer_pheno_evo/data/FELINE2/PhenotypesAllArmsCohort2/CPMPhenotpyeLandscape_C2ProjectC2RevisednewMacrophagesDC.RData"))#C1umap,u_dat,DAY,cell_types_all,ARMS,Subtype,
load( file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Cohort2Metadata/Cohort2Metadata.RData") #metadd,cohort1metadd,annotation.file,compdataLU

dd1 <- u_dat %>% dplyr::select(Cell.ID:file_string,V1,V2)

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
  data.table(Cohort="Discovery" ,    ddCohort1%>%select(Treatment,Cell.ID,Sample,Patient.Study.ID,Day,ARM,dynamic_class3,Celltype,Celltype_subtype,V1,V2) ),
  data.table(Cohort="Validation" ,  ddCohort2%>%select(Treatment,Cell.ID,Sample,Patient.Study.ID,Day,ARM,dynamic_class3,Celltype,Celltype_subtype,V1,V2))
)


dd0 <-rbind(
  data.table(Cohort="Validation", metadd %>%select(Cell.ID,Celltype,Celltype_subtype,dynamic_class3,Patient.Study.ID,Day,ARM)) , 
  data.table(Cohort="Discovery",cohort1metadd%>%select(Cell.ID,Celltype,Celltype_subtype,dynamic_class3,Patient.Study.ID,Day,ARM) ))

dd0[,Treatmentlab:= "Combination ribociclib"]
dd0[ARM=="A",Treatmentlab:= "Letrozole alone"]
dd0[Cell.ID%in% joineddd[Celltype=="Macrophages"][Celltype_subtype=="Macrophages"][V2>0]$Cell.ID, Celltype_subtype:= "M2 macrophages" ]
dd0[Cell.ID%in% joineddd[Celltype=="Macrophages"][Celltype_subtype=="Macrophages"][V2<0]$Cell.ID, Celltype_subtype:= "M1 macrophages" ]
dd0[Celltype=="Macrophages"]$Celltype_subtype%>%table()
dd0[Celltype_subtype=="Endothelial cells",Celltype_subtype:="Vas-Endo"]

dd0$CellAnnot1<-dd0$Celltype
dd0$CellAnnot2<-dd0$Celltype_subtype
dd0[Celltype%in%c("Cancer cells","Normal epithelial cells"),CellAnnot1:="Epithelial cells"]
dd0[Celltype%in%c("Cancer cells","Normal epithelial cells"),CellAnnot2:="Epithelial cells"]


# Total cell count per sample= sampling effort
dd0[,totalCount:= length(Cell.ID) ,by=c("Patient.Study.ID","Day","Cohort","dynamic_class3","Treatmentlab")]

# count of each cell type's abundance and frequency
countTable0 <- data.table( dd0 %>% group_by(Patient.Study.ID,Day,Cohort,Treatmentlab,dynamic_class3,
                                            CellAnnot2,totalCount)%>% summarise(count=n() , frac=n()/unique(totalCount) ) ) 
# ggplot(countTable0[CellAnnot2 %in% c("M1 macrophages", "M2 macrophages") ], 
#        aes(log(frac), x=log(1+Day), fill=dynamic_class3,col=dynamic_class3)) + theme_classic(base_size=16)+
#   geom_violin(aes(group=interaction(Day,dynamic_class3)),alpha=0.6,position=position_dodge(width=0.9))+
#   geom_point(size=2,position=position_dodge(width=0.9)) + 
#   facet_grid(Treatmentlab~CellAnnot2) + theme(axis.text.x=element_text(angle=90)) +
#   geom_smooth(method="gam", formula=y~s(x,k=3),se=F)

countTable0wideM1M2<-countTable0[CellAnnot2 %in% unique(dd0[CellAnnot1=="Macrophages"]$CellAnnot2) ]%>% select(-count)%>%spread(CellAnnot2,frac,fill=0) 
#countTable0wideM1M2<-countTable0[CellAnnot2 %in% c("M1 macrophages", "M2 macrophages") ]%>% select(-count)%>%spread(CellAnnot2,frac,fill=0) 
names(countTable0wideM1M2) <- gsub(" ","_",names(countTable0wideM1M2))
countTable0wideM1M2[,M1M2ratio:= M1_macrophages/(  M1_macrophages + M2_macrophages ) ]
countTable0wideM1M2[,Myeloidtot := ( DC + M1_macrophages + M2_macrophages + Monocytes)  ]
#countTable0wideM1M2[,M1M2ratio:= M1_macrophages/(M1_macrophages+M2_macrophages) ]
#countTable0wideM1M2[,Myeloidtot := (M1_macrophages+M2_macrophages) ]

ggplot(countTable0wideM1M2, 
       aes(M1M2ratio, x=log(1+Day), fill=dynamic_class3,col=dynamic_class3)) + theme_classic(base_size=16)+
  geom_violin(aes(group=interaction(Day,dynamic_class3)),alpha=0.6,position=position_dodge(width=0.9))+
  geom_point(aes(shape=Cohort,group=interaction(Day,dynamic_class3)),size=2,position=position_dodge(width=0.9)) + 
  facet_grid(~Treatmentlab) + theme(aspect.ratio=1,axis.text.x=element_text(angle=90)) +
  geom_smooth(method="gam", formula=y~s(x,k=3),se=F,size=2)+
  scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"),guide="none")+
  scale_fill_npg(name="Tumor response", labels=c("Resistant", "Sensitive"),guide="none")+
  scale_x_continuous(name="Day", breaks= log(1+c(0,14,180)),labels=c(0,14,180))


countTable0wideM1M2[,Success:=1]
countTable0wideM1M2[dynamic_class3=="Non-response",Success:=0]
countTable0wideM1M2<-countTable0wideM1M2[is.finite(M1M2ratio)]
countTable0wideM1M2[,x:=M1M2ratio ]
countTable0wideM1M2[,y:=log(Myeloidtot) ]


glmmodAll <-glm(Success~ M1M2ratio*Treatmentlab,family="binomial",data=countTable0wideM1M2[])
summary(glmmodAll)

glmmodR <-glm(Success~ x*y,family="binomial",data=countTable0wideM1M2[Treatmentlab!="Letrozole alone"])
glmmodL <-glm(Success~ x*y,family="binomial",data=countTable0wideM1M2[Treatmentlab=="Letrozole alone"])
glmmodRs <-glm(Success~ M1M2ratio,family="binomial",data=countTable0wideM1M2[Treatmentlab!="Letrozole alone"][Day<180])
glmmodLs <-glm(Success~ M1M2ratio,family="binomial",data=countTable0wideM1M2[Treatmentlab=="Letrozole alone"][Day<180])
#countTable0wideM1M2$predglm <- as.vector( predict(glmmod,data=countTable0wideM1M2,type="response")>0.5 )
summary(glmmodRs)
summary(glmmodLs)


# a more fancy plot with the boundary added
boundary <- function(model, data, class = NULL, predict_type = "class",
                     resolution = 100, showgrid = TRUE, ...) {
  if(!is.null(class)) cl <- data[,class] else cl <- 1
  data <- data.frame( model$data)[,names(model$data)%in% names(coef(model)) ] 
  k <- length(unique(cl))
  # make grid of x values
  ys <- seq(min(data$y), max(data$y), length.out = resolution)
  
  # get model coefficients to predict boundary
  b0 = coef(model)["(Intercept)"]
  bx = coef(model)["x"]
  by = coef(model)["y"]
  bxy = coef(model)["x:y"]
  xbound <- -(b0 + by*ys )/(bx + bxy*ys) 
  model$data$ClassPrediction <-  predict(model,type="response")>0.5
  # data to plot classification boundary
  boundary<- data.table(y=ys,x=xbound)[x>=min(data$x)& x<=max(data$x)] 
  # plot
  ggplot(model$data,aes(y=y,x=x))+
    geom_point(aes(col=as.factor(Class),shape= as.factor(1*ClassPrediction)))+
    scale_color_discrete(name="Cell lineage")+
    scale_shape_discrete(name="Prediction")+
    geom_path(data= boundary ) + theme_classic() +
    labs(y="UMAP 2",x="UMAP 1")
  
  return(list(data= model$data, boundary=boundary))
  
}

outR <- boundary(model=glmmodR, 
                data=data.frame(countTable0wideM1M2[Treatmentlab!="Letrozole alone"]),
                #class="Class",
                predict_type="response",
                resolution = 1000)




ggplot(countTable0wideM1M2, 
       aes(log(Myeloidtot), x=log(1+Day), fill=dynamic_class3,col=dynamic_class3)) + theme_classic(base_size=16)+
  geom_violin(aes(group=interaction(Day,dynamic_class3)),alpha=0.6,position=position_dodge(width=0.9))+
  geom_point(aes(shape=Cohort),size=2,position=position_dodge(width=0.9)) + 
  facet_grid(~Treatmentlab) + theme(aspect.ratio=1,axis.text.x=element_text(angle=90)) +
  geom_smooth(method="gam", formula=y~s(x,k=3),se=F,size=2)+
  scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"),guide="none")+
  scale_fill_npg(name="Tumor response", labels=c("Resistant", "Sensitive"),guide="none")+
  scale_x_continuous(name="Day", breaks= log(1+c(0,14,180)),labels=c(0,14,180)) 

ggplot(countTable0wideM1M2, 
       aes(M1M2ratio, x=log(1+Day), fill=dynamic_class3,col=dynamic_class3)) + theme_classic(base_size=16)+
  geom_violin(aes(group=interaction(Day,dynamic_class3)),alpha=0.6,position=position_dodge(width=0.9))+
  geom_point(aes(group=interaction(Day,dynamic_class3),shape=Cohort),size=2,position=position_dodge(width=0.9)) + 
  facet_grid(~Treatmentlab) + theme(aspect.ratio=1,axis.text.x=element_text(angle=90)) +
  geom_smooth(method="gam",formula=y~s(x,k=3),se=F,size=2)+
  scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"),guide="none")+
  scale_fill_npg(name="Tumor response", labels=c("Resistant", "Sensitive"),guide="none")+
  scale_x_continuous(name="Day", breaks= log(1+c(0,14,180)),labels=c(0,14,180)) 


ggplot(countTable0wideM1M2,#[Day!=180], 
       aes(y=log(Myeloidtot), x=M1M2ratio)) + theme_classic(base_size=18)+
  #geom_violin(aes(group=interaction(Day,dynamic_class3)),alpha=0.6,position=position_dodge(width=0.9))+
  geom_point(aes(fill=dynamic_class3,col=dynamic_class3,shape=as.factor(Day)),size=3) + 
  #facet_wrap(Cohort~Treatmentlab) + 
  facet_wrap(~Treatmentlab,ncol=1) + 
  theme(aspect.ratio=1 ,panel.spacing=unit(20,"pt"))+
  #geom_smooth(method="gam", formula=y~s(x,k=3),se=F,size=2)+
  scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  scale_shape_discrete(name="Day")+
  scale_y_continuous(name="Myeloid cell fraction",labels=c(0.001,0.01,0.1,0.5),breaks=log(c(0.001,0.01,0.1,0.5)))+
  labs(x="Proportion myeloid cells in M1 state")+
  geom_path(data=outR$boundary,aes(x=x,y=y ) ,col="black",fill="black",size=2,show.legend=FALSE)

ggsave(filename=paste0("/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cohort 2 results/Composition Myeloid differentiation not abundance predicts response to ribo not letrozole.png"),width=8,height=6)


lm(M1M2ratio~dynamic_class3,countTable0wideM1M2[Day!=180][Treatmentlab=="Combination ribociclib"] )%>%summary()
lm(M1M2ratio~dynamic_class3,countTable0wideM1M2[Day!=180][Treatmentlab!="Combination ribociclib"] )%>%summary()


ggplot(countTable0wideM1M2[Day!=180], 
       aes(x= M1M2ratio, fill= dynamic_class3, col= dynamic_class3 )) + theme_classic(base_size= 26)+
  geom_density(aes(group= interaction(dynamic_class3)),alpha=0.6) + 
  facet_wrap(~ Treatmentlab,ncol=1) + theme(aspect.ratio= 1) +
  #geom_smooth(method="gam", formula=y~s(x,k=3),se=F,size=2)+
  scale_color_npg(name= "Tumor response", labels= c("Resistant", "Sensitive"))+
  scale_fill_npg(name= "Tumor response", labels= c("Resistant", "Sensitive"))+
  labs(x="Proportion myeloid cells in M1 state",y="Probability density") + geom_rug(size=2)+theme(legend.position="bottom")
ggsave(filename=paste0("/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cohort 2 results/Early treatment Composition abundance of M1 myeloid differentiation links to response to ribo not letrozole.png"),width=10,height=7)
paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Ribo and Letrozole Early treatment M1 abundance high in response to ribo not letrozole.png"),height=10,width=10)


ggplot(countTable0wideM1M2[Day!=180], 
       aes(y= (M1M2ratio), x= dynamic_class3, fill= dynamic_class3 )) + theme_classic(base_size= 26)+
  #geom_violin(aes(group= interaction(dynamic_class3)),alpha=0.6) + 
  geom_boxplot()+
  geom_boxplot(outlier.colour=NA, position= position_dodge() ,col="black")+
  stat_boxplot(geom="errorbar",position=position_dodge(1.75),width=0.5)+#geom_smooth(method="lm")+
  facet_wrap(~ Treatmentlab,ncol=1) + theme(aspect.ratio= 1) +
  #geom_smooth(method="gam", formula=y~s(x,k=3),se=F,size=2)+
  scale_fill_npg(name= "Tumor response", labels= c("Resistant", "Sensitive"))+
  labs(y="Proportion myeloid cells in M1 state",x="Tumor response") + 
  geom_point(size=1.5, col="black")+theme(legend.position="none")

ggsave(paste0(paperfile,"Ribo and Letrozole Early treatment M1 abundance high in response to ribo not letrozole boxplot.png"),height=10,width=10)

ggplot(countTable0wideM1M2[Day!=180], 
       aes(y= (M1M2ratio), x= dynamic_class3, fill= dynamic_class3 )) + theme_classic(base_size= 26)+
  #geom_violin(aes(group= interaction(dynamic_class3)),alpha=0.6) + 
  geom_boxplot()+
  geom_boxplot(outlier.colour=NA, position= position_dodge() ,col="black")+
  stat_boxplot(geom="errorbar",position=position_dodge(1.75),width=0.5)+#geom_smooth(method="lm")+
  facet_wrap(~ Treatmentlab,ncol=1) + theme(aspect.ratio= 1) +
  #geom_smooth(method="gam", formula=y~s(x,k=3),se=F,size=2)+
  scale_fill_npg(name= "Tumor response", labels= c("Resistant", "Sensitive"))+
  labs(y="Proportion myeloid cells in M1 state",x="Tumor response") + 
  geom_point(size=1.5, col="black")+theme(legend.position="none")+
  theme(axis.title=element_blank(),  axis.text=element_blank(),strip.text = element_blank(), legend.position="none")
#ggsave(paste0(paperfile,"BLANK Ribo and Letrozole Early treatment M1 abundance high in response to ribo not letrozole boxplot.png.png"),height=10,width=10)




ggplot(countTable0wideM1M2[Day!=180], 
       aes(x= M1M2ratio, fill= dynamic_class3, col= dynamic_class3 )) + theme_classic(base_size= 16)+
  geom_density(aes(group= interaction(dynamic_class3)),alpha=0.4) + 
  facet_grid(~ Treatmentlab) + theme(aspect.ratio= 1) +
  #geom_smooth(method="gam", formula=y~s(x,k=3),se=F,size=2)+
  scale_color_npg(name= "Tumor response", labels= c("Resistant", "Sensitive"))+
  scale_fill_npg(name= "Tumor response", labels= c("Resistant", "Sensitive"))+
  labs(x="Proportion myeloid cells in M1 state",y="Density") + geom_rug(size=2) +
  theme(axis.title=element_blank(),  axis.text=element_blank(),strip.text = element_blank(), legend.position="none")
#ggsave(paste0(paperfile,"BLANK Ribo and Letrozole Early treatment M1 abundance high in response to ribo not letrozole.png"),height=10,width=10)





#ggsave(paste0(paperfile,"BLANK Ribo and Letrozole Hallmark IFN gamma versus Tumor wide IL15.png"),height=10,width=10)



psimplif<-countTable0wideM1M2[Day!=180]
psimplif[,M2Hi:=T]
psimplif[M1M2ratio<0.5,M2Hi:=F]
psimplif1 <-data.table( psimplif%>%group_by(dynamic_class3,M2Hi,Treatmentlab)%>%summarise(n=length(M2Hi)) )
psimplif1[,prop:=n/sum(n),by=c("dynamic_class3","Treatmentlab")]
ggplot(psimplif1[M2Hi==T], 
       aes(y= prop, x= dynamic_class3, fill= M2Hi )) + theme_classic(base_size= 16)+
  geom_bar(stat="identity")+
  #geom_density(aes(group= interaction(dynamic_class3)),alpha=0.4) + 
  facet_grid(~ Treatmentlab) + theme(aspect.ratio= 1) +
  #geom_smooth(method="gam", formula=y~s(x,k=3),se=F,size=2)+
  scale_color_npg(name= "Tumor response", labels= c("Resistant", "Sensitive"))+
  scale_fill_npg(name= "Tumor response", labels= c("Resistant", "Sensitive"))+
  labs(x="Proportion myeloid cells in M1 state",y="Frequency density")+geom_rug(size=2)
ggsave(filename=paste0("/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cohort 2 results/Early treatment Composition abundance of M1 myeloid differentiation links to response to ribo not letrozole.png"),width=10,height=7)



ggplot(countTable0wideM1M2,#[Day!=180], 
       aes(x= M1M2ratio, fill= dynamic_class3, col= dynamic_class3 )) + theme_classic(base_size= 16)+
  geom_density(aes(group= interaction(dynamic_class3,Day)),alpha=0.4) + 
  facet_grid(Day~ Treatmentlab) + theme(aspect.ratio= 1) +
  #geom_smooth(method="gam", formula=y~s(x,k=3),se=F,size=2)+
  scale_color_npg(name= "Tumor response", labels= c("Resistant", "Sensitive"))+
  scale_fill_npg(name= "Tumor response", labels= c("Resistant", "Sensitive"))+
  labs(x="Proportion myeloid cells in M1 state",y="Frequency density")+geom_rug(size=2)
ggsave(filename=paste0("/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cohort 2 results/Composition abundance of M1 myeloid differentiation links to response to ribo not letrozole.png"),width=10,height=7)


ggplot(countTable0wideM1M2,#[Day!=180], 
       aes(x= M1M2ratio, fill= dynamic_class3, col= dynamic_class3 )) + theme_classic(base_size= 16)+
  geom_histogram(position = "dodge",bins=20,aes(group= interaction(dynamic_class3,Day)),alpha=1) + 
  facet_grid(Day~ Treatmentlab) + theme(aspect.ratio= 1) +
  scale_color_npg(name= "Tumor response", labels= c("Resistant", "Sensitive"))+
  scale_fill_npg(name= "Tumor response", labels= c("Resistant", "Sensitive"))+
  labs(x="Proportion myeloid cells in M1 state")

summary( lm(M1M2ratio~dynamic_class3*Treatmentlab,data=countTable0wideM1M2))
summary( lm(M1M2ratio~dynamic_class3,data=countTable0wideM1M2[Treatmentlab=="Letrozole alone"]))
summary( lm(M1M2ratio~dynamic_class3,data=countTable0wideM1M2[Treatmentlab!="Letrozole alone"]))


ggplot(countTable0wideM1M2%>%gather(CellAnnot2,frac,M1_macrophages: M2_macrophages), 
       aes(log(frac), x=log(1+Day), fill=dynamic_class3,col=dynamic_class3)) + theme_classic(base_size=16)+
  geom_violin(aes(group=interaction(Day,dynamic_class3)),alpha=0.6,position=position_dodge(width=0.9))+
  geom_point(aes(shape=Cohort),size=2,position=position_dodge(width=0.9)) + 
  facet_grid(CellAnnot2~Treatmentlab) + theme(aspect.ratio=1,axis.text.x=element_text(angle=90)) +
  geom_smooth(method="gam", formula=y~s(x,k=3),se=F,size=2)+
  scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"),guide="none")+
  scale_fill_npg(name="Tumor response", labels=c("Resistant", "Sensitive"),guide="none")+
  scale_x_continuous(name="Day", breaks= log(1+c(0,14,180)),labels=c(0,14,180))

# 
# #ggplot(countTable0[Celltype=="T cells"],aes(log(1+count),x=Celltype,col=dynamic_class3))+
# #  geom_boxplot()+geom_point(position=position_dodge(width=0.75))
# ggplot(countTable0[][Celltype=="T cells"],aes(logit(frac), x= Celltype, col= dynamic_class3))+
#   geom_boxplot()+geom_point(position=position_dodge(width=0.75))+ facet_wrap(Treatmentlab~Day)
# ggplot(countTable0[Cohort=="Discovery"][Celltype=="T cells"],aes(logit(frac), x= Celltype, col= dynamic_class3))+
#   geom_boxplot()+geom_point(position=position_dodge(width=0.75))+ facet_wrap(Treatmentlab~Day)
# 
# ggplot(countTable0[Cohort=="Discovery"][Celltype=="Macrophages"],aes(logit(frac), x= Celltype, col= dynamic_class3))+
#   geom_boxplot()+geom_point(position=position_dodge(width=0.75))+ facet_wrap(Treatmentlab~Day)
# 
# ggplot(countTable0[][Celltype=="Epithelial"],aes(logit(frac), x= Celltype, col= dynamic_class3))+
#   geom_boxplot()+geom_point(position=position_dodge(width=0.75))+ facet_wrap(Treatmentlab~Day)
# 
# 
# ggplot(countTable0[][Celltype=="Macrophages"],aes(logit(frac), x= log(1+Day), fill=dynamic_class3,col= dynamic_class3))+
#   geom_violin(aes(group=interaction(Day,dynamic_class3) ),alpha=0.6)+geom_point(position=position_dodge(width=2))+ 
#   facet_wrap(Cohort~Treatmentlab, ncol=1)+theme_classic()+
#   geom_smooth(method="gam", formula=y~s(x,k=3), se=F)
# 
# ggplot(countTable0[][Celltype=="Cancer cells"],aes(logit(frac), x= log(1+Day), fill=dynamic_class3,col= dynamic_class3))+
#   geom_violin(aes(group=interaction(Day,dynamic_class3) ),alpha=0.6)+geom_point(position=position_dodge(width=2))+ 
#   facet_wrap(Cohort~Treatmentlab, ncol=1)+theme_classic()+
#   geom_smooth(method="gam", formula=y~s(x,k=3), se=F)
# 
# 
# ggplot(countTable0[][Celltype=="Cancer cells"],aes(logit(frac), x= log(1+Day), col= dynamic_class3))+
#   geom_boxplot()+geom_point(position=position_dodge(width=0.75))+ facet_wrap(Treatmentlab~Day)
# 
#

### Full cohort all samples
# fraction table spread with cell types across columns #[Timepoint==0]
#fracTable1 <- spread( countTable0[Day==180][totalCount>150]%>%dplyr::select(-count), CellAnnot2,frac,fill=0)
fracTable1 <- spread( countTable0[][]%>%dplyr::select(-count), CellAnnot2,frac,fill=0)
fracTable1$rowID <- 1:nrow(fracTable1)

# fraction matrix without metadata. name rows by ID to retain match to metadata
fracTable2 <- as.matrix( fracTable1%>%dplyr::select(-c("Patient.Study.ID", "Day", "Cohort","Treatmentlab" ,"dynamic_class3", "totalCount","rowID")) )
rownames(fracTable2)<- fracTable1$rowID

rowannot <- data.frame(fracTable1 %>% dplyr::select(rowID,Cohort, Day,dynamic_class3, Treatmentlab,totalCount))
rownames(rowannot) <- rowannot$rowID; rowannot$rowID <- NULL
pheatmap::pheatmap(logit( 1e-11+fracTable2*(.9999)),annotation_row=rowannot)





## DAY 0 analysis
# count table (raw counts) with cell types across columns #
#countTable1 <- spread( countTable0[Timepoint== 0][ARM!= "A"][] %>% dplyr::select(-frac), Celltype, count, fill= 0)
#countTable1 <- spread( countTable0[Day== 0][][] %>% dplyr::select(-frac), CellAnnot2, count, fill= 0)
countTable1 <- spread( countTable0[Day != 180][Cohort=="Discovery"][] %>% dplyr::select(-frac), CellAnnot2, count, fill= 0)
countTable1$rowID <- 1:nrow(countTable1)

# fraction matrix without metadata. name rows by ID to retain match to metadata
countTable2 <- as.matrix( countTable1 %>% dplyr::select(-c("Patient.Study.ID", "Day", "Cohort", "Treatmentlab" ,"dynamic_class3", "totalCount","rowID")) )
rownames(countTable2) <- countTable1$rowID
clrdd <- clr(countTable2)
ilrdd <- ilr(countTable2)

rowannot <- data.frame(countTable1 %>% dplyr::select(rowID, dynamic_class3,Cohort))#,Cohort,totalCount))
rownames(rowannot) <- rowannot$rowID;rowannot$rowID<-NULL
pheatmap::pheatmap(clrdd, annotation_row= rowannot)





countTable0$CellAnnot2 <- gsub("\\+","",countTable0$CellAnnot2)
countTable0$CellAnnot2 <- gsub(" ","_",countTable0$CellAnnot2)

# fraction table spread with cell types across columns #[Timepoint==0]
#fracTable1 <- spread( countTable0[Day==180][totalCount>150]%>%dplyr::select(-count), CellAnnot2,frac,fill=0)
fracTable1 <- spread( countTable0[][]%>%dplyr::select(-count), CellAnnot2,frac,fill=0)
fracTable1$rowID <- 1:nrow(fracTable1)
hist(log(fracTable1[]$totalCount))

# fraction matrix without metadata. name rows by ID to retain match to metadata
fracTable2 <- as.matrix( fracTable1%>%dplyr::select(-c("Patient.Study.ID", "Day", "Cohort","Treatmentlab" ,"dynamic_class3", "totalCount","rowID")) )
rownames(fracTable2)<- fracTable1$rowID

rowannot <- data.frame(fracTable1 %>% dplyr::select(rowID, dynamic_class3,Day, Treatmentlab,totalCount))
rownames(rowannot) <- rowannot$rowID; rowannot$rowID <- NULL
pheatmap::pheatmap(logit( 1e-10+fracTable2*0.999),annotation_row=rowannot)

distmat <- vegdist(fracTable2, method = "bray")
pheatmap::pheatmap(as.matrix(distmat),annotation_col=rowannot)

NMS <- metaMDS(distmat,distance = "bray")
set.seed(1)
u1<-umap::umap(fracTable2 , distance = "manhattan",n_components=3, n_neighbors=5)#, n_neighbors=20)
gmm1<-Mclust(u1$layout[,2:3],G=3,modelNames="EEV")

set.seed(1)
u1<-umap::umap(fracTable2 , distance = "manhattan",n_components=2, n_neighbors=10)#, n_neighbors=20)
require(mclust)
gmm1 <- Mclust(u1$layout,G=3,modelNames="EEV")
# gathr output of NMDS
umapout<-data.table( cbind( fracTable1,  u1$layout ),archetype=gmm1$classification ) #,princomp(scores(comm.bc.mds))$scores
names(umapout)<-gsub("-","_",names(umapout))
names(umapout)<-gsub(" ","_",names(umapout))
umapout$archetype<- factor(umapout$archetype,levels=c("1","3","2" ))
names(umapout)
ggplot(umapout, aes(x=V1, y=V2, 
                    col=dynamic_class3,
                    #col=logit(1e-5+Cancer_cells),#col=logit(1e-5+T_cells),
                    #col=logit(1e-5+M1_macrophages),#col=logit(1e-5+T_cells),
                    shape= (archetype)) ) + geom_point(size=2.5) + theme_classic(base_size=18) + 
  facet_wrap(Day~Treatmentlab)#+ scale_color_npg(name="Tumor response") + scale_size_continuous(name="Sample size")
ggplot(umapout, aes(x=V1, y=V2, 
                    col=archetype,
                    shape= (archetype)) ) + geom_point(size=2.5) + theme_classic(base_size=18)  #+ scale_color_npg(name="Tumor response") + scale_size_continuous(name="Sample size")


compplotdd <- umapout[]%>%gather(celltype,frac,Adipocytes:T_cells)
compplotdd$celltype<- factor(compplotdd$celltype, levels=c("Epithelial","T_cells","B_cells", "Macrophages",
                                                           "Adipocytes","Fibroblasts","Endothelial_cells","Pericytes"))

p1<- ggplot( compplotdd,
             aes(x=celltype,y=frac,group=archetype,col=archetype))+geom_point(size=2,position=position_dodge(width=0.33))+
  #facet_wrap(~archetype,ncol=1)+
  theme_classic(base_size=22)+
  labs(y="Fraction", x= "Cell type")+  scale_color_npg(name="Archetype", labels=c("Cancer \n dominated","Immune hot \n and diverse","Fibroblast/endothlial \n enriched")) + 
  theme(aspect.ratio=1,axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(legend.position="none")

p2 <- ggplot(umapout, aes(x=V1, y=V2, col=as.factor(archetype) ) )+ geom_point(size=4)+ theme_classic(base_size=22) +  #,fill=as.factor(archetype)
  scale_color_npg(name="Archetype", labels= c("Cancer \n dominated", "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched")) + 
  labs(x="Dimension A",y="Dimension B")+ theme(aspect.ratio=1)#+
#stat_ellipse(alpha=0.03,geom = "polygon",level =0.95)

ggsave(p1, filename="/Users/jason/Downloads/Feine1and2_ArchetypeCompositions.png")
ggsave(p2, filename="/Users/jason/Downloads/Feine1and2_Archetypes.png")






rm(list=ls())
# read data
dd0 <- fread("/Users/jason/Dropbox/Jason FELINE data/Feline_harmonizedAnnotation_metadata.txt")%>%
  dplyr::select(orig.ident,Patient,Cell,Timepoint,Cohort,seurat_clusters,integrated_snn_res.0.8, integrated_snn_res.0.7,Annotation, Patient.Study.ID, ARM,dynamic_class,dynamic_class2,dynamic_class3)
names(dd0)

# temp merge all epithelial cells due to classification uncertainty
dd0$Celltype<-dd0$Annotation
dd0[ grep("pithelial",Annotation), Celltype :="Epithelial"] 
# Total cell count per sample= sampling effort
dd0[,totalCount:= length(Cell) ,by=c("Patient.Study.ID","Timepoint","Cohort","dynamic_class3")]

# count of each cell type's abundance and frequency
countTable0 <- data.table( dd0%>%group_by(Patient.Study.ID,Timepoint,Cohort,ARM,dynamic_class3,Celltype,totalCount)%>% summarise(count=n() , frac=n()/unique(totalCount) ) )

# count table (raw counts) with cell types across columns #
countTable1 <- spread( countTable0[Timepoint==0][][]%>%dplyr::select(-frac),Celltype,count,fill=0)
countTable1$rowID <- 1:nrow(countTable1)

# fraction matrix without metadata. name rows by ID to retain match to metadata
countTable2 <- as.matrix( countTable1%>%dplyr::select(-c("Patient.Study.ID", "Timepoint", "Cohort","ARM" ,"dynamic_class3", "totalCount","rowID")) )
rownames(countTable2)<- countTable1$rowID

rowannot <- data.frame(countTable1%>%dplyr::select(rowID, dynamic_class3,Cohort))#,Cohort,totalCount))
rownames(rowannot) <-rowannot$rowID;rowannot$rowID<-NULL
clrdd <- clr(countTable2)
ilrdd <- ilr(countTable2)
pheatmap::pheatmap(clrdd,annotation_row=rowannot)

# fraction table spread with cell types across columns #[Timepoint==0]
fracTable1 <- spread( countTable0[Timepoint==0][][totalCount>150]%>%dplyr::select(-count),Celltype,frac,fill=0)
fracTable1$rowID <- 1:nrow(fracTable1)
hist(log(fracTable1[]$totalCount))

# fraction matrix without metadata. name rows by ID to retain match to metadata
fracTable2 <- as.matrix( fracTable1%>%dplyr::select(-c("Patient.Study.ID", "Timepoint", "Cohort","ARM" ,"dynamic_class3", "totalCount","rowID")) )
rownames(fracTable2)<- fracTable1$rowID

rowannot <- data.frame(fracTable1%>%dplyr::select(rowID,Cohort, dynamic_class3, ARM,totalCount))
rownames(rowannot) <-rowannot$rowID;rowannot$rowID<-NULL
pheatmap::pheatmap(logit( 1e-10+fracTable2),annotation_row=rowannot)

distmat <- vegdist(fracTable2, method = "bray")
pheatmap::pheatmap(as.matrix(distmat),annotation_col=rowannot)

set.seed(1)
#u1<-umap::umap(fracTable2 , distance = "manhattan",n_components=3, n_neighbors=5)#, n_neighbors=20)
#gmm1<-Mclust(u1$layout[,2:3],G=3,modelNames="EEV")

set.seed(12321)
u1<-umap::umap(fracTable2 , distance = "manhattan",n_components=2, n_neighbors=10)#, n_neighbors=20)
gmm1 <- Mclust(u1$layout,G=3,modelNames="EEV")

# gathr output of NMDS
umapout<-data.table( cbind( fracTable1,  u1$layout ),archetype=gmm1$classification ) #,princomp(scores(comm.bc.mds))$scores
names(umapout)<-gsub("-","_",names(umapout))
names(umapout)<-gsub(" ","_",names(umapout))
umapout$archetype<- factor(umapout$archetype,levels=c("2","3","1" ))

ggplot(umapout, aes(x=V1, y=V2, col= (archetype)) ) + geom_point(size=2.5) + theme_classic(base_size=18) #+ scale_color_npg(name="Tumor response") + scale_size_continuous(name="Sample size")

ggplot(umapout, aes(x=V1, y=V2, col=logit(1e-5+T_cell),shape= (archetype)) ) + geom_point(size=2.5) + theme_classic(base_size=18) #+ scale_color_npg(name="Tumor response") + scale_size_continuous(name="Sample size")
ggplot(umapout, aes(x=V1, y=V2, col=logit(1e-5+Epithelial),shape= (archetype)) ) + geom_point(size=2.5) + theme_classic(base_size=18) #+ scale_color_npg(name="Tumor response") + scale_size_continuous(name="Sample size")
ggplot(umapout, aes(x=V1, y=V2, col=log((1e-5+Fibroblast)/(1e-5+T_cell)),shape= (archetype)) ) + geom_point(size=2.5) + theme_classic(base_size=18) #+ scale_color_npg(name="Tumor response") + scale_size_continuous(name="Sample size")

compplotdd <- umapout[]%>%gather(celltype,frac,Adipocyte:Vas_Endo)
compplotdd$celltype<- factor(compplotdd$celltype, levels=c("Epithelial","T_cell","NK_cell","B_cell","Plasma_cell", "Macrophage","Monocyte","Dendritic_cell",
                                                           "Adipocyte","Fibroblast","Vas_Endo","Pericyte","Lym_Endo"))

p1 <- ggplot( compplotdd,
              aes(x=celltype,y=logit(1e-5+frac),group=archetype,col=archetype))+geom_point(size=2,position=position_dodge(width=0.33))+
  #facet_wrap(~archetype,ncol=1)+
  theme_classic(base_size=22)+
  labs(y="Fraction", x= "Cell type")+  scale_color_npg(name="Archetype", labels=c("Cancer \n dominated","Immune hot \n and diverse","Fibroblast/endothlial \n enriched")) + 
  theme(aspect.ratio=1,axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(legend.position="none")

p1 <- ggplot( compplotdd,
              aes(x=celltype,y=logit(1e-5+frac),col=archetype,group=interaction(archetype,celltype)))+
  geom_violin(size=2,position=position_dodge(width=0.33))+
  #facet_wrap(~archetype,ncol=1)+
  theme_classic(base_size=22)+
  labs(y="Fraction", x= "Cell type")+  scale_color_npg(name="Archetype", labels=c("Cancer \n dominated","Immune hot \n and diverse","Fibroblast/endothlial \n enriched")) + 
  theme(aspect.ratio=1,axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(legend.position="none")


p2 <- ggplot(umapout, aes(x=V1, y=V2, col=as.factor(archetype) ) )+ geom_point(size=4)+ theme_classic(base_size=22) +  #,fill=as.factor(archetype)
  scale_color_npg(name="Archetype", labels= c("Cancer \n dominated", "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched")) + 
  labs(y="Immune cell infiltration",x="Cancer dominance")+ theme(aspect.ratio=1)#+


p3 <- ggplot(umapout, aes(x=V1, y=V2, col=logit(1e-5+T_cell), shape=as.factor(archetype) ) )+ geom_point(size=4)+ theme_classic(base_size=22) +  #,fill=as.factor(archetype)
  scale_shape_discrete(name="Archetype", labels= c("Cancer \n dominated", "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched")) + 
  labs(y="Immune abundance",x="Cancer dominance")+ theme(aspect.ratio=1)+
  scale_color_viridis_c(option="A")
cowplot::plot_grid(p2,p3)

p3b <- ggplot(umapout, aes(x=V1, y=V2, col=logit(1e-5+Epithelial), shape=as.factor(archetype) ) )+ geom_point(size=4)+ theme_classic(base_size=22) +  #,fill=as.factor(archetype)
  scale_shape_discrete(name="Archetype", labels= c("Cancer \n dominated", "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched")) + 
  labs(y="Immune abundance",x="Cancer dominance")+ theme(aspect.ratio=1)+
  scale_color_viridis_c(option="A")
cowplot::plot_grid(p2,p3)
#stat_ellipse(alpha=0.03,geom = "polygon",level =0.95)

ggsave(p1, filename="/Users/jason/Downloads/Feine1and2_ArchetypeCompositions2_0.png")
ggsave(p2, filename="/Users/jason/Downloads/Feine1and2_Archetypes2_0.png")



cordd <- data.table(t(cor(umapout%>%dplyr::select(V1:V2),umapout%>%dplyr::select(Adipocyte:Vas_Endo))),keep.rownames = T)
cordd[order(-abs(V1))]
cordd[order(-abs(V2))]









rm(list=ls())
# read data
dd0 <- fread("/Users/jason/Dropbox/Jason FELINE data/Feline_harmonizedAnnotation_metadata.txt")%>%
  dplyr::select(orig.ident,Patient,Cell,Timepoint,Cohort,seurat_clusters,integrated_snn_res.0.8, integrated_snn_res.0.7,Annotation, Patient.Study.ID, ARM,dynamic_class,dynamic_class2,dynamic_class3)
names(dd0)

# temp merge all epithelial cells due to classification uncertainty
dd0$Celltype<-dd0$Annotation
dd0[ grep("pithelial",Annotation), Celltype :="Epithelial"] 
# Total cell count per sample= sampling effort
dd0[,totalCount:= length(Cell) ,by=c("Patient.Study.ID","Timepoint","Cohort","dynamic_class3")]

# count of each cell type's abundance and frequency
countTable0 <- data.table( dd0%>%group_by(Patient.Study.ID,Timepoint,Cohort,ARM,dynamic_class3,Celltype,totalCount)%>% summarise(count=n() , frac=n()/unique(totalCount) ) )

# count table (raw counts) with cell types across columns #
countTable1 <- spread( countTable0[][][]%>%dplyr::select(-frac),Celltype,count,fill=0)
countTable1$rowID <- 1:nrow(countTable1)

# fraction matrix without metadata. name rows by ID to retain match to metadata
countTable2 <- as.matrix( countTable1%>%dplyr::select(-c("Patient.Study.ID", "Timepoint", "Cohort","ARM" ,"dynamic_class3", "totalCount","rowID")) )
rownames(countTable2)<- countTable1$rowID

rowannot <- data.frame(countTable1%>%dplyr::select(rowID, dynamic_class3,Cohort))#,Cohort,totalCount))
rownames(rowannot) <-rowannot$rowID;rowannot$rowID<-NULL
clrdd <- clr(countTable2)
ilrdd <- ilr(countTable2)
pheatmap::pheatmap(clrdd,annotation_row=rowannot)

# fraction table spread with cell types across columns #[Timepoint==0]
fracTable1 <- spread( countTable0[][][]%>%dplyr::select(-count),Celltype,frac,fill=0)
fracTable1$rowID <- 1:nrow(fracTable1)
hist(log(fracTable1[]$totalCount))

# fraction matrix without metadata. name rows by ID to retain match to metadata
fracTable2 <- as.matrix( fracTable1%>%dplyr::select(-c("Patient.Study.ID", "Timepoint", "Cohort","ARM" ,"dynamic_class3", "totalCount","rowID")) )
rownames(fracTable2)<- fracTable1$rowID

rowannot <- data.frame(fracTable1%>%dplyr::select(rowID,Cohort, dynamic_class3, ARM,totalCount))
rownames(rowannot) <-rowannot$rowID;rowannot$rowID<-NULL
pheatmap::pheatmap(logit( 1e-10+fracTable2),annotation_row=rowannot)

distmat <- vegdist(fracTable2, method = "bray")
pheatmap::pheatmap(as.matrix(distmat),annotation_col=rowannot)


set.seed(123)
u1<-umap::umap(fracTable2 , distance = "manhattan",n_components=2, n_neighbors=10)#, n_neighbors=20)
gmm1 <- Mclust(u1$layout,G=3,modelNames="EEV")
#gmm1 <- Mclust(u1$layout,G=3,modelNames="EEI")


# gathr output of NMDS
umapout<-data.table( cbind( fracTable1,  u1$layout ),archetype=gmm1$classification ) #,princomp(scores(comm.bc.mds))$scores
names(umapout)<-gsub("-","_",names(umapout))
names(umapout)<-gsub(" ","_",names(umapout))
umapout$archetype<- factor(umapout$archetype,levels=c("1","3","2" ))


cordd <- data.table(t(cor(umapout%>%dplyr::select(V1:V2),umapout%>%dplyr::select(Adipocyte:Vas_Endo))),keep.rownames = T)
cordd[order(-abs(V1))]
cordd[order(-abs(V2))]


ggplot(umapout, aes(x=V1, y=V2, col= (archetype)) ) + geom_point(size=2.5) + theme_classic(base_size=18) #+ scale_color_npg(name="Tumor response") + scale_size_continuous(name="Sample size")
ggplot(umapout, aes(x=V1, y=V2, col= log(1+Timepoint)) ) + geom_point(size=2.5) + theme_classic(base_size=18) #+ scale_color_npg(name="Tumor response") + scale_size_continuous(name="Sample size")
ggplot(umapout, aes(x=V1, y=V2, col= (Cohort)) ) + geom_point(size=2.5) + theme_classic(base_size=18) #+ scale_color_npg(name="Tumor response") + scale_size_continuous(name="Sample size")

ggplot(umapout, aes(x=V1, y=V2, col=logit(1e-5+Fibroblast),shape= (archetype)) ) + geom_point(size=2.5) + theme_classic(base_size=18) #+ scale_color_npg(name="Tumor response") + scale_size_continuous(name="Sample size")

ggplot(umapout, aes(x=V1, y=V2, col=logit(1e-5+T_cell),shape= (archetype)) ) + geom_point(size=2.5) + theme_classic(base_size=18) #+ scale_color_npg(name="Tumor response") + scale_size_continuous(name="Sample size")
ggplot(umapout, aes(x=V1, y=V2, col=logit((1e-5+Epithelial)*0.99),shape= (archetype)) ) + geom_point(size=2.5) + theme_classic(base_size=18) #+ scale_color_npg(name="Tumor response") + scale_size_continuous(name="Sample size")
ggplot(umapout, aes(x=V1, y=V2, col=log((1e-5+Fibroblast)/(1e-5+T_cell)),shape= (archetype)) ) + geom_point(size=2.5) + theme_classic(base_size=18) #+ scale_color_npg(name="Tumor response") + scale_size_continuous(name="Sample size")

compplotdd <- umapout[]%>%gather(celltype,frac,Adipocyte:Vas_Endo)
compplotdd$celltype<- factor(compplotdd$celltype, levels=c("Epithelial","T_cell","NK_cell","B_cell","Plasma_cell", "Macrophage","Monocyte","Dendritic_cell",
                                                           "Adipocyte","Fibroblast","Vas_Endo","Pericyte","Lym_Endo"))

p1 <- ggplot( compplotdd,
              aes(x=celltype,y=logit(1e-5+frac),group=archetype,col=archetype))+geom_point(size=2,position=position_dodge(width=0.33))+
  #facet_wrap(~archetype,ncol=1)+
  theme_classic(base_size=22)+
  labs(y="Fraction", x= "Cell type")+  scale_color_npg(name="Archetype", labels=c("Cancer \n dominated","Immune hot \n and diverse","Fibroblast/endothlial \n enriched")) + 
  theme(aspect.ratio=1,axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(legend.position="none")

p1 <- ggplot( compplotdd,
              aes(x=celltype,y=logit(1e-5+frac),col=archetype,group=interaction(archetype,celltype)))+
  geom_violin(size=2,position=position_dodge(width=0.33))+
  #facet_wrap(~archetype,ncol=1)+
  theme_classic(base_size=22)+
  labs(y="Fraction", x= "Cell type")+  scale_color_npg(name="Archetype", labels=c("Cancer \n dominated","Immune hot \n and diverse","Fibroblast/endothlial \n enriched")) + 
  theme(aspect.ratio=1,axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(legend.position="none")


p2 <- ggplot(umapout, aes(y=-V1, x=V2, col=as.factor(archetype) ) )+ geom_point(size=4)+ theme_classic(base_size=22) +  #,fill=as.factor(archetype)
  scale_color_npg(name="Archetype", labels= c("Cancer \n dominated", "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched")) + 
  labs(y="Immune cell infiltration",x="Cancer dominance")+ theme(aspect.ratio=1)#+

p2 <- ggplot(umapout, aes(y=-V1, x=V2, col=as.factor(archetype) ) )+ geom_point(size=4)+ theme_classic(base_size=22) +  #,fill=as.factor(archetype)
  scale_color_npg(name="Archetype", labels= c("Cancer \n dominated", "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched")) + 
  labs(y="TME composition axis 2",x="TME composition axis 1")+ theme(aspect.ratio=1)#+
ggsave(p2, filename="/Users/jason/Downloads/Feine1and2_Archetypes3_0_newlab.png")

p2Blank <- ggplot(umapout, aes(y=-V1, x=V2, col=as.factor(archetype) ) )+ geom_point(size=8)+ theme_classic(base_size=22) +  #,fill=as.factor(archetype)
  scale_color_npg(name="Archetype", labels= c("Cancer \n dominated", "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched")) + 
  labs(y="TME composition axis 2",x="TME composition axis 1")+ theme(aspect.ratio=1)+
  theme(
    axis.text=element_blank(),axis.title.x = element_blank(),
    axis.title.y = element_blank(),legend.position="none")
ggsave(p2Blank, filename="/Users/jason/Downloads/BLANK_Feine1and2_Archetypes3_0.png")


p3 <- ggplot(umapout, aes(y=-V1, x=V2, col=logit(1e-5+T_cell), shape=as.factor(archetype) ) )+ geom_point(size=4)+ theme_classic(base_size=22) +  #,fill=as.factor(archetype)
  scale_shape_discrete(name="Archetype", labels= c("Cancer \n dominated", "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched")) + 
  labs(y="Immune abundance",x="Cancer dominance")+ theme(aspect.ratio=1)+
  scale_color_viridis_c(option="A")
cowplot::plot_grid(p2,p3)

p3b <- ggplot(umapout, aes(x=V1, y=V2, col=logit(1e-5+Epithelial), shape=as.factor(archetype) ) )+ geom_point(size=4)+ theme_classic(base_size=22) +  #,fill=as.factor(archetype)
  scale_shape_discrete(name="Archetype", labels= c("Cancer \n dominated", "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched")) + 
  labs(y="Immune abundance",x="Cancer dominance")+ theme(aspect.ratio=1)+
  scale_color_viridis_c(option="A")
cowplot::plot_grid(p2,p3)
#stat_ellipse(alpha=0.03,geom = "polygon",level =0.95)

ggsave(p1, filename="/Users/jason/Downloads/Feine1and2_ArchetypeCompositions3_0.png")
ggsave(p2, filename="/Users/jason/Downloads/Feine1and2_Archetypes3_0.png")




shannonDiv <- diversity(fracTable2,index = "shannon")

CompositionDivTable <- merge(umapout , data.table(rowID=as.numeric(names(shannonDiv)),shannonDiv=shannonDiv), by="rowID")

p4 <- ggplot(CompositionDivTable, aes(y=-V1, x=V2, shape=as.factor(archetype), col=shannonDiv ) )+ geom_point(size=4)+ theme_classic(base_size=22) +  #,fill=as.factor(archetype)
  scale_color_viridis_c(name="Shannon Diversity")+
  scale_shape_discrete(name="Archetype", labels= c("Cancer \n dominated", "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched")) + 
  labs(y="TME composition axis 2",x="TME composition axis 1")+ theme(aspect.ratio=1)#+

ggsave(p4, filename="/Users/jason/Downloads/Feline1and2_ArchetypesWithDiversity3_0.png")

p4Blank <- ggplot(CompositionDivTable, aes(y=-V1, x=V2, shape=as.factor(archetype) , col=shannonDiv) )+ geom_point(size=8)+ theme_classic(base_size=22) +  #,fill=as.factor(archetype)
  scale_color_viridis_c(name="Shannon Diversity")+
  scale_shape_discrete(name="Archetype", labels= c("Cancer \n dominated", "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched")) +  labs(y="TME composition axis 2",x="TME composition axis 1")+ theme(aspect.ratio=1)+
  theme(
    axis.text=element_blank(),axis.title.x = element_blank(),
    axis.title.y = element_blank(),legend.position="none")
ggsave(p4Blank, filename="/Users/jason/Downloads/BLANK_Feline1and2_ArchetypesWithDiversity3_0.png")


tmp<-CompositionDivTable
tmp[, res:=0]
tmp[archetype==1, res:=1]
lm(shannonDiv~res,data=tmp)%>%summary()
cohort1CancerHeterogeneity <- readxl::read_excel("/Users/jason/Dropbox/Jason FELINE data/Table S_HeterogeneityDiveristy.xlsx")%>%
  dplyr::select(Patient, Timepoint, H)
#cohort1CancerHeterogeneity$H <- as.numeric(cohort1CancerHeterogeneity$H)
mrgeCancerHet<- merge(cohort1CancerHeterogeneity,unique(dd0%>%dplyr::select(Patient.Study.ID,Patient)),by="Patient")
#mrgeCancerHet$Timepoint<- as.numeric(mrgeCancerHet$Timepoint)



TMEvsCancerDiv<- merge(CompositionDivTable,mrgeCancerHet,
                       by=c("Patient.Study.ID","Timepoint"))
TMEvsCancerDiv$H

ggplot(TMEvsCancerDiv,
       aes(x=H,y=shannonDiv, col=archetype))+geom_point()+
  geom_smooth(method="lm")+
  theme_classic()+
  labs(y="Tumor composition diversity", "Cancer subclone diversity")


p5 <- ggplot(TMEvsCancerDiv[Timepoint==0],
             aes(x=(H),y=(shannonDiv)))+geom_point(size=3)+
  geom_smooth(method="gam",size=2)+
  theme_classic(base_size = 22)+
  labs(y="Tumor composition diversity", x="Cancer subclone diversity")
ggsave(p5, filename="/Users/jason/Downloads/Feline1_PRETREAT_TMEvsCancerDiversity3_0.png")

p5Blank <- ggplot(TMEvsCancerDiv,
                  aes(x=(H),y=(shannonDiv)))+geom_point(size=5)+
  geom_smooth(method="gam",size=2)+
  theme_classic(base_size = 22)+
  theme(
    axis.text=element_blank(),axis.title.x = element_blank(),
    axis.title.y = element_blank(),legend.position="none")
ggsave(p5Blank, filename="/Users/jason/Downloads/BLANK_Feline1_TMEvsCancerDiversity3_0.png")



p5Blankpretreat <- ggplot(TMEvsCancerDiv[Timepoint==0],
                          aes(x=(H),y=(shannonDiv)))+geom_point(size=5)+
  geom_smooth(method="gam",size=2)+
  theme_classic(base_size = 22)+
  theme(
    axis.text=element_blank(),axis.title.x = element_blank(),
    axis.title.y = element_blank(),legend.position="none")
ggsave(p5Blankpretreat, filename="/Users/jason/Downloads/BLANK_Feline1_PRETREATMENT_TMEvsCancerDiversity3_0.png")

cor.test(TMEvsCancerDiv[Timepoint==180]$shannonDiv,TMEvsCancerDiv[Timepoint==180]$H)

cor.test(TMEvsCancerDiv[Timepoint==0]$shannonDiv,TMEvsCancerDiv[Timepoint==0]$H)
plot(TMEvsCancerDiv[ARM!="A"]$shannonDiv,TMEvsCancerDiv[ARM!="A"]$H)
summary( mgcv::gam(shannonDiv~H,data=TMEvsCancerDiv[]))
anova(lm(shannonDiv~H,data=TMEvsCancerDiv[]))
summary(lm(shannonDiv~H,data=TMEvsCancerDiv))
summary(lm(shannonDiv~H,data=TMEvsCancerDiv[Timepoint==0]))
summary(lm(shannonDiv~H,data=TMEvsCancerDiv[Timepoint==180]))


ggplot(TMEvsCancerDiv, aes(y=-V1, x=V2, shape=as.factor(archetype), col=H ) )+ geom_point(size=4)+ theme_classic(base_size=22) +  #,fill=as.factor(archetype)
  scale_color_viridis_c(name="Shannon Diversity")+
  scale_shape_discrete(name="Archetype", labels= c("Cancer \n dominated", "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched")) + 
  labs(y="Immune cell infiltration",x="Cancer dominance")+ theme(aspect.ratio=1)#+



ggplot(TMEvsCancerDiv,
       aes(x=(shannonDiv),y=dynamic_class3 ))+geom_point(size=3)+  # (shannonDiv)
  geom_smooth(method="gam",size=2)+
  theme_classic(base_size = 22)+
  facet_wrap(~Timepoint)+
  labs(x="Tumor composition diversity", y="Tumor response \n during treatment")

predictive_value <- CompositionDivTable[Cohort=="C1"]
predictive_value[,Outcome:=1]
predictive_value[dynamic_class3=="Non-response",Outcome:=0]

ggplot(  predictive_value,
         aes(x=(shannonDiv),y=Outcome ))+geom_point(size=3)+  # (shannonDiv)
  geom_smooth(method="gam",method.args = list(family = "binomial"),size=2)+
  theme_classic(base_size = 22)+
  facet_wrap(~Timepoint)+
  labs(x="Tumor composition diversity", y="Tumor response during treatment")

predictive_value <- TMEvsCancerDiv[Cohort=="C1"]
predictive_value[,Outcome:=1]
predictive_value[dynamic_class3=="Non-response",Outcome:=0]
predictive_value[,Ribo:=T]
predictive_value[ARM=="A",Ribo:=F]
CompositionDivTable[,Ribo:=T]
CompositionDivTable[ARM=="A",Ribo:=F]

poutArch<-ggplot(  CompositionDivTable, #CompositionDivTable[Cohort=="C1"]
                   aes(x=archetype,fill=dynamic_class3,group=interaction(archetype,dynamic_class3) ))+
  geom_bar()+  # (shannonDiv)
  #geom_bar(position = "fill")+ 
  theme_classic(base_size = 22)+
  facet_wrap(~Timepoint)+coord_flip()+
  labs(x="Tumor archetype", y="Number of tumor with each outcome")+theme(aspect.ratio=1)+scale_fill_npg(name="Tumor response",labels=c("resistant","sensitive"))+
  scale_x_discrete(name="Archetype", labels= c("Cancer \n dominated", "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched")) 
ggsave(poutArch, filename="/Users/jason/Downloads/Feine1and2_Archetypes3_0_outcomes_time.png",width=14,height=7)


ggplot(  predictive_value,
         aes(x=(H), y=Outcome ))+geom_point(size=3)+  # (shannonDiv)
  geom_smooth(method="gam",method.args = list(family = "binomial"),size=2)+
  theme_classic(base_size = 22)+
  facet_wrap(~Timepoint)+
  labs(x="Tumor composition diversity", y="Tumor response during treatment")




ggplot(umapout[Timepoint<180], aes(y=-V1, x=V2, group=Patient.Study.ID,col=as.factor(archetype) ) )+ geom_point(size=4)+
  geom_path(data=umapout[Timepoint<180],col="black",arrow = arrow(length=unit(0.20,"cm"), ends="first", type = "closed"))+ theme_classic(base_size=22) +  #,fill=as.factor(archetype)
  scale_color_npg(name="Archetype", labels= c("Cancer \n dominated", "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched")) + 
  labs(y="Immune cell infiltration",x="Cancer dominance")+ theme(aspect.ratio=1)#+


p2arr <- ggplot(umapout, aes(y=-V1, x=V2,group=Patient.Study.ID, col=as.factor(archetype) ) )+ geom_point(size=4)+ theme_classic(base_size=22) +  #,fill=as.factor(archetype)
  scale_color_npg(name="Archetype", labels= c("Cancer \n dominated", "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched")) + 
  labs(y="Immune cell infiltration",x="Cancer dominance")+ theme(aspect.ratio=1)+
  #geom_path(data=umapout[Timepoint!=14],col="black",arrow = arrow(length=unit(0.20,"cm"), ends="first", type = "closed"))+
  facet_wrap(~dynamic_class3)+
  geom_path(data=umapout[Timepoint<180],col="black",arrow = arrow(length=unit(0.20,"cm"), ends="first", type = "closed"))+
  geom_path(data=umapout[Timepoint>0],col="black",arrow = arrow(length=unit(0.20,"cm"), ends="first", type = "closed"))#+
# facet_wrap(~ARM)
ggsave(p2arr, filename="/Users/jason/Downloads/Feine1and2_Archetypes3_0_arrows.png")

ggplot(umapout, aes(y=-V1, x=V2,group=Patient.Study.ID, col=as.factor(Timepoint) ) )+ geom_point(size=4)+ theme_classic(base_size=22) +  #,fill=as.factor(archetype)
  scale_color_npg(name="Timepoint")+#, labels= c("Cancer \n dominated", "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched")) + 
  labs(y="Immune cell infiltration",x="Cancer dominance")+ theme(aspect.ratio=1)+
  facet_wrap(dynamic_class3~ARM)


ggplot(umapout[Timepoint%in%c(180)], aes(y=-V1, x=V2,group=Patient.Study.ID, col=dynamic_class3 ) )+ geom_point(size=4)+ theme_classic(base_size=22) +  #,fill=as.factor(archetype)
  scale_color_npg(name="Timepoint")+#, labels= c("Cancer \n dominated", "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched")) + 
  labs(y="Immune cell infiltration",x="Cancer dominance")+ theme(aspect.ratio=1)


Timepoint_names <- c(`0` = "Pre-treatment",`180` = "Post-treatment")
p6<- ggplot(umapout[Timepoint%in%c(0,180)], aes(y=-V1, x=V2 ) )+
  #stat_ellipse(aes(fill=as.factor(archetype)),alpha=0.2,geom = "polygon",level =0.95)+
  geom_point(aes(col=dynamic_class3),size=4)+ theme_classic(base_size=22) +  #,fill=as.factor(archetype)
  scale_color_npg(name="Tumor \n response",labels=c("Resistant","Sensitive"))+#, labels= c("Cancer \n dominated", "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched")) + 
  labs(y="Immune cell infiltration",x="Cancer dominance")+ theme(aspect.ratio=1)+
  facet_wrap(~Timepoint, labeller = as_labeller(Timepoint_names))

ggsave(p6, filename="/Users/jason/Downloads/Feine1and2_ArchetypesResponseTime.png")

p6Blankday0<- ggplot(umapout[Timepoint%in%c(0)], aes(y=-V1, x=V2 ) )+
  #stat_ellipse(aes(fill=as.factor(archetype)),alpha=0.2,geom = "polygon",level =0.95)+
  geom_point(aes(col=dynamic_class3),size=8)+ theme_classic(base_size=22) +  #,fill=as.factor(archetype)
  scale_color_npg(name="Tumor \n response",labels=c("Resistant","Sensitive"))+#, labels= c("Cancer \n dominated", "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched")) + 
  labs(y="Immune cell infiltration",x="Cancer dominance")+ theme(aspect.ratio=1)+
  theme(
    axis.text=element_blank(),axis.title.x = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.title.y = element_blank(),legend.position="none")

ggsave(p6Blankday0, filename="/Users/jason/Downloads/BLANK_Feine1and2_ArchetypesResponseTimeDay0.png")


p6Blankday180<- ggplot(umapout[Timepoint%in%c(180)], aes(y=-V1, x=V2 ) )+
  #stat_ellipse(aes(fill=as.factor(archetype)),alpha=0.2,geom = "polygon",level =0.95)+
  geom_point(aes(col=dynamic_class3),size=8)+ theme_classic(base_size=22) +  #,fill=as.factor(archetype)
  scale_color_npg(name="Tumor \n response",labels=c("Resistant","Sensitive"))+#, labels= c("Cancer \n dominated", "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched")) + 
  labs(y="Immune cell infiltration",x="Cancer dominance")+ theme(aspect.ratio=1)+
  theme(
    axis.text=element_blank(),axis.title.x = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.title.y = element_blank(),legend.position="none")

ggsave(p6Blankday180, filename="/Users/jason/Downloads/BLANK_Feine1and2_ArchetypesResponseTimeDay180.png")



average_frac1 <-data.table(compplotdd)# data.table(CompositionDivTable[]%>%
#        gather(celltype,frac,Adipocyte:Vas_Endo))
average_frac1[,celltypegrp:=celltype]
average_frac1[celltype=="Epithelial",celltypegrp:="Epithelial"]
average_frac1[celltype%in%c("Plasma_cell","B_cell","Dendritic_cell","Macrophage","Monocyte","T_cell","NK_cell"),celltypegrp:="Immune"]
average_frac1[celltype%in%c("Fibroblast","Pericyte","Vas_Endo","Lym_Endo"),celltypegrp:="Fibroblast/Endothelial"]
average_frac2 <- data.table( average_frac1%>%group_by(archetype,celltypegrp,Patient.Study.ID,Timepoint,ARM,dynamic_class3)%>%summarise(frac=sum(frac)) )
average_frac2


average_frac2 <- data.table(average_frac2%>%group_by(archetype,celltypegrp)%>%summarise(med=median(frac)))
average_frac2[,avFrac:=med/sum(med),by="archetype"]

average_frac2$celltypegrp<- factor(average_frac2$celltypegrp,levels=c("Epithelial","Immune","Fibroblast/Endothelial","Adipocyte"))
p7_Cancer<-ggplot(average_frac2[archetype=="1"], aes(x="", y=avFrac, fill=celltypegrp))+
  geom_bar(width = 1, stat = "identity")+scale_fill_discrete()+ coord_polar("y", start=0)+  
  theme_minimal()+scale_fill_npg(name="Cell type")+
  theme(
    axis.text.x=element_blank(),axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold"))#+  facet_wrap(~archetype)
p7_ImmuneDiverse<-ggplot(average_frac2[archetype=="3"], aes(x="", y=avFrac, fill=celltypegrp))+
  geom_bar(width = 1, stat = "identity")+scale_fill_discrete()+ coord_polar("y", start=0)+  
  theme_minimal()+scale_fill_npg(name="Cell type")+
  theme(
    axis.text.x=element_blank(),axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold"))#+  facet_wrap(~archetype)
p7_FibroEnd<-ggplot(average_frac2[archetype=="2"], aes(x="", y=avFrac, fill=celltypegrp))+
  geom_bar(width = 1, stat = "identity")+scale_fill_discrete()+ coord_polar("y", start=0)+  
  theme_minimal()+scale_fill_npg(name="Cell type")+
  theme(
    axis.text.x=element_blank(),axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold"))#+  facet_wrap(~archetype)

ggsave(p7_Cancer, filename="/Users/jason/Downloads/Feine1and2_ArchetypesPiechartCancer.png")
ggsave(p7_FibroEnd, filename="/Users/jason/Downloads/Feine1and2_ArchetypesPiechartFibroEnd.png")
ggsave(p7_ImmuneDiverse, filename="/Users/jason/Downloads/Feine1and2_ArchetypesPiechartImmuneDiverse.png")

p7_Cancer<-ggplot(average_frac2[archetype=="1"], aes(x="", y=avFrac, fill=celltypegrp))+
  geom_bar(width = 1, stat = "identity")+scale_fill_discrete()+ coord_polar("y", start=0)+  
  theme_minimal()+scale_fill_npg(name="Cell type")+
  theme(legend.position = "none",
        axis.text.x=element_blank(),axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        axis.ticks = element_blank(),
        plot.title=element_text(size=14, face="bold"))#+  facet_wrap(~archetype)
p7_ImmuneDiverse<-ggplot(average_frac2[archetype=="3"], aes(x="", y=avFrac, fill=celltypegrp))+
  geom_bar(width = 1, stat = "identity")+scale_fill_discrete()+ coord_polar("y", start=0)+  
  theme_minimal()+scale_fill_npg(name="Cell type")+
  theme(legend.position = "none",
        axis.text.x=element_blank(),axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        axis.ticks = element_blank(),
        plot.title=element_text(size=14, face="bold"))#+  facet_wrap(~archetype)
p7_FibroEnd<-ggplot(average_frac2[archetype=="2"], aes(x="", y=avFrac, fill=celltypegrp))+
  geom_bar(width = 1, stat = "identity")+scale_fill_discrete()+ coord_polar("y", start=0)+  
  theme_minimal()+scale_fill_npg(name="Cell type")+
  theme(legend.position = "none",
        axis.text.x=element_blank(),axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        axis.ticks = element_blank(),
        plot.title=element_text(size=14, face="bold"))#+  facet_wrap(~archetype)
ggsave(p7_Cancer, filename="/Users/jason/Downloads/BLANK_Feine1and2_ArchetypesPiechartCancer.png")
ggsave(p7_FibroEnd, filename="/Users/jason/Downloads/BLANK_Feine1and2_ArchetypesPiechartFibroEnd.png")
ggsave(p7_ImmuneDiverse, filename="/Users/jason/Downloads/BLANK_Feine1and2_ArchetypesPiechartImmuneDiverse.png")

# average_frac2 <- data.table(umapout[archetype==2]%>%gather(celltype,frac,Adipocyte:Vas_Endo)%>%group_by(celltype)%>%summarise(med=inv.logit(mean(logit(1e-5+frac)))))
# average_frac2[,avFrac:=med/sum(med)]
# average_frac3 <- data.table(umapout[archetype==3]%>%gather(celltype,frac,Adipocyte:Vas_Endo)%>%group_by(celltype)%>%summarise(med=inv.logit(mean(logit(1e-5+frac)))))
# average_frac3[,avFrac:=med/sum(med)]
# p1<- ggplot(average_frac1, aes(x="", y=avFrac, fill=celltype))+
#   geom_bar(width = 1, stat = "identity")+scale_fill_discrete()+ coord_polar("y", start=0)+  blank_theme +
#   theme(axis.text.x=element_blank())
# p2<- ggplot(average_frac2, aes(x="", y=avFrac, fill=celltype))+
#   geom_bar(width = 1, stat = "identity")+scale_fill_discrete()+ coord_polar("y", start=0)+  blank_theme +
#   theme(axis.text.x=element_blank())
# p3<- ggplot(average_frac3, aes(x="", y=avFrac, fill=celltype))+
#   geom_bar(width = 1, stat = "identity")+scale_fill_discrete()+ coord_polar("y", start=0)+  blank_theme +
#   theme(axis.text.x=element_blank())
# 









# fraction table spread with cell types across columns #[Timepoint==0]
fracTableFull1 <- spread( countTable0[][ARM!="A"][]%>%dplyr::select(-count),Celltype,frac,fill=0)
fracTableFull1$rowID <- 1:nrow(fracTableFull1)
hist(log(fracTableFull1[]$totalCount))

# fraction matrix without metadata. name rows by ID to retain match to metadata
fracTableFull2 <- as.matrix( fracTableFull1%>%dplyr::select(-c("Patient.Study.ID", "Timepoint", "Cohort","ARM" ,"dynamic_class3", "totalCount","rowID")) )
rownames(fracTableFull2)<- fracTableFull1$rowID

rowannotFull <- data.frame(fracTableFull1%>%dplyr::select(rowID,Cohort, Timepoint,dynamic_class3, ARM,totalCount))
rownames(rowannotFull) <-rowannotFull$rowID;rowannotFull$rowID<-NULL

ttt <- logit( 1e-5+fracTableFull2)
ttt[is.nan(ttt)] <- logit(1-1e-5)
for(i in 1:ncol(ttt)){
  ttt[,i]<- scale(ttt[,i])
}
pheatmap::pheatmap( ttt , annotation_row= rowannotFull )


distmatFull <- vegdist(fracTableFull2, method = "bray")
pheatmap::pheatmap(as.matrix(distmatFull),annotation_col=rowannotFull)


# upreds <- predict(data=fracTableFull2,u1)
# plot(upreds)
# plot(u1$layout)
# 
# 
# 
# NMS <- metaMDS(distmat,
#                distance = "bray",
# )
# set.seed(1)
# u1FULL<-umap::umap(fracTableFull2 , distance = "manhattan",n_components=3, n_neighbors=5)#, n_neighbors=20)
# gmm1<-Mclust(u1$layout[,2:3],G=3,modelNames="EEV")
# 
# u1<-umap::umap(fracTable2 , distance = "manhattan",n_components=2, n_neighbors=10)#, n_neighbors=20)
# require(mclust)
# gmm1 <- Mclust(u1$layout,G=3,modelNames="EEV")
# # gathr output of NMDS
# umapout<-data.table( cbind( fracTable1,  u1$layout ),archetype=gmm1$classification ) #,princomp(scores(comm.bc.mds))$scores
# names(umapout)<-gsub("-","_",names(umapout))
# names(umapout)<-gsub(" ","_",names(umapout))
# umapout$archetype<- factor(umapout$archetype,levels=c("1","3","2" ))
# 
# ggplot(umapout, aes(x=V1, y=V2, col=logit(1e-5+T_cell),shape= (archetype)) ) + geom_point(size=2.5) + theme_classic(base_size=18) #+ scale_color_npg(name="Tumor response") + scale_size_continuous(name="Sample size")
# 
# compplotdd <- umapout[]%>%gather(celltype,frac,Adipocyte:Vas_Endo)
# compplotdd$celltype<- factor(compplotdd$celltype, levels=c("Epithelial","T_cell","NK_cell","B_cell","Plasma_cell", "Macrophage","Monocyte","Dendritic_cell",
#                                                            "Adipocyte","Fibroblast","Vas_Endo","Pericyte","Lym_Endo"))
# 
# p1<- ggplot( compplotdd,
#              aes(x=celltype,y=frac,group=archetype,col=archetype))+geom_point(size=2,position=position_dodge(width=0.33))+
#   #facet_wrap(~archetype,ncol=1)+
#   theme_classic(base_size=22)+
#   labs(y="Fraction", x= "Cell type")+  scale_color_npg(name="Archetype", labels=c("Cancer \n dominated","Immune hot \n and diverse","Fibroblast/endothlial \n enriched")) + 
#   theme(aspect.ratio=1,axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   theme(legend.position="none")
# 
# p2 <- ggplot(umapout, aes(x=V1, y=V2, col=as.factor(archetype) ) )+ geom_point(size=4)+ theme_classic(base_size=22) +  #,fill=as.factor(archetype)
#   scale_color_npg(name="Archetype", labels= c("Cancer \n dominated", "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched")) + 
#   labs(x="Dimension A",y="Dimension B")+ theme(aspect.ratio=1)#+
# #stat_ellipse(alpha=0.03,geom = "polygon",level =0.95)
# 
# ggsave(p1, filename="/Users/jason/Downloads/Feine1and2_ArchetypeCompositions.png")
# ggsave(p2, filename="/Users/jason/Downloads/Feine1and2_Archetypes.png")
# 
#    
# 
# 
# 
# 
# 
# 
# 
# ggplot( umapout[]%>%gather(celltype,frac,Adipocyte:Vas_Endo)%>%
#           group_by(celltype,archetype)%>%summarise(frac=inv.logit(mean(logit(1e-5+frac)))),
#         aes(x=celltype,y=(1e-5+frac),col=archetype))+geom_point()+
#   facet_wrap(~archetype,ncol=1)+theme_classic()
# 
#   
# 
# ggplot(umapout, aes(x=V2, y=V3, col=as.factor(archetype) ) )+ geom_point(size=4)+ theme_classic(base_size=22) + 
#   scale_color_npg(name="Archetype", labels=c("Cancer /n dominated","Immune /n hot","Fibroblast/endothlial \n enriched")) + labs(x="Dimension A",y="Dimension B")
# 
# ggplot(umapout, aes(x=V2, y=V3, col=as.factor(archetype) ) )+ geom_point(size=4)+ theme_classic(base_size=22) + 
#   scale_color_npg(name="Archetype", labels=c("Cancer /n dominated","Immune /n hot","Fibroblast/endothlial \n enriched")) + labs(x="Dimension A",y="Dimension B")
# 
# ggplot(umapout, aes(x=V2, y=V3, col=logit(1e-5+T_cell),shape= rev(dynamic_class3)) ) + geom_point(size=2.5) + theme_classic(base_size=18) #+ scale_color_npg(name="Tumor response") + scale_size_continuous(name="Sample size")
# ggplot(umapout, aes(x=V2, y=V3, col=logit(1e-5+Adipocyte),shape= rev(dynamic_class3)) ) + geom_point(size=2.5) + theme_classic(base_size=18) #+ scale_color_npg(name="Tumor response") + scale_size_continuous(name="Sample size")
# ggplot(umapout, aes(x=V2, y=V3, col=logit(1e-5+Macrophage),shape= rev(dynamic_class3)) ) + geom_point(size=2.5) + theme_classic(base_size=18) #+ scale_color_npg(name="Tumor response") + scale_size_continuous(name="Sample size")
# ggplot(umapout, aes(x=V2, y=V3, col=logit(1e-5+Fibroblast),shape= rev(dynamic_class3)) ) + geom_point(size=2.5) + theme_classic(base_size=18) #+ scale_color_npg(name="Tumor response") + scale_size_continuous(name="Sample size")
# 
# average_frac1 <- data.table(umapout[archetype==1]%>%gather(celltype,frac,Adipocyte:Vas_Endo)%>%group_by(celltype)%>%summarise(med=inv.logit(mean(logit(1e-5+frac)))))
# average_frac1[,avFrac:=med/sum(med)]
# average_frac2 <- data.table(umapout[archetype==2]%>%gather(celltype,frac,Adipocyte:Vas_Endo)%>%group_by(celltype)%>%summarise(med=inv.logit(mean(logit(1e-5+frac)))))
# average_frac2[,avFrac:=med/sum(med)]
# average_frac3 <- data.table(umapout[archetype==3]%>%gather(celltype,frac,Adipocyte:Vas_Endo)%>%group_by(celltype)%>%summarise(med=inv.logit(mean(logit(1e-5+frac)))))
# average_frac3[,avFrac:=med/sum(med)]
# blank_theme <- theme_minimal()+
#   theme(
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank(),
#     panel.border = element_blank(),
#     panel.grid=element_blank(),
#     axis.ticks = element_blank(),
#     plot.title=element_text(size=14, face="bold"))
# p1<- ggplot(average_frac1, aes(x="", y=avFrac, fill=celltype))+
#   geom_bar(width = 1, stat = "identity")+scale_fill_discrete()+ coord_polar("y", start=0)+  blank_theme +
#   theme(axis.text.x=element_blank())
# p2<- ggplot(average_frac2, aes(x="", y=avFrac, fill=celltype))+
#   geom_bar(width = 1, stat = "identity")+scale_fill_discrete()+ coord_polar("y", start=0)+  blank_theme +
#   theme(axis.text.x=element_blank())
# p3<- ggplot(average_frac3, aes(x="", y=avFrac, fill=celltype))+
#   geom_bar(width = 1, stat = "identity")+scale_fill_discrete()+ coord_polar("y", start=0)+  blank_theme +
#   theme(axis.text.x=element_blank())
# 
# 
# cordd <- data.table(t(cor(umapout%>%dplyr::select(V1:V3),umapout%>%dplyr::select(Adipocyte:Vas_Endo))),keep.rownames = T)
# cordd[order(-abs(V1))]
# cordd[order(-abs(V2))]
# cordd[order(-abs(V3))]
# 
# 
# ggplot(umapout, aes(x=Epithelial, y=V3) ) + geom_point() 
# 
# # Plot first 2 dimensions and colour by response group
# ggplot(umapout, aes(x=V2, y=V3, col=as.factor(archetype) ) )+ geom_point(size=3)+ theme_classic(base_size=22) + scale_color_npg(name="Archetype") + labs(x="Dimension A",y="Dimension B")
# 
# 
# ggplot(umapout, aes(x=V3, y=V1, col=dynamic_class3,size=totalCount) ) + geom_point() + theme_classic(base_size=18) + scale_color_npg(name="Tumor response") + scale_size_continuous(name="Sample size")
# ggplot(umapout, aes(x=V1, y=V2, col=logit(1e-5+T_cell+NK_cell),shape= ARM,size=totalCount) ) + geom_point() + theme_classic(base_size=18) #+ scale_color_npg(name="Tumor response") + scale_size_continuous(name="Sample size")
# ggplot(umapout, aes(x=V1, y=V2, col=logit(1e-5+Epithelial),shape= ARM,size=totalCount) ) + geom_point() + theme_classic(base_size=18) #+ scale_color_npg(name="Tumor response") + scale_size_continuous(name="Sample size")
# 
# ggplot(umapout, aes(x=V1, y=V3, col=Cohort) ) + geom_point(size=2.5) + theme_classic(base_size=18) + scale_color_npg(name="Tumor response") + scale_size_continuous(name="Sample size")
# 
# ggplot(umapout, aes(x=V1, y=V3, col=logit(1e-5+Epithelial),shape= rev(dynamic_class3)) ) + geom_point(size=2.5) + theme_classic(base_size=18) #+ scale_color_npg(name="Tumor response") + scale_size_continuous(name="Sample size")
# ggplot(umapout, aes(x=V2, y=V3, col=logit(1e-5+T_cell),shape= rev(dynamic_class3)) ) + geom_point(size=2.5) + theme_classic(base_size=18) #+ scale_color_npg(name="Tumor response") + scale_size_continuous(name="Sample size")
# ggplot(umapout, aes(x=V1, y=V3, col=logit(1e-5+Fibroblast),shape= rev(dynamic_class3)) ) + geom_point(size=2.5) + theme_classic(base_size=18) + scale_color_viridis_c()#+ scale_color_npg(name="Tumor response") + scale_size_continuous(name="Sample size")
# ggplot(umapout, aes(x=V1, y=V3, col=logit(1e-5+Epithelial),shape= rev(dynamic_class3)) ) + geom_point(size=2.5) + theme_classic(base_size=18) + scale_color_viridis_c()#+ scale_color_npg(name="Tumor response") + scale_size_continuous(name="Sample size")
# 
# ggplot(umapout, aes(x=V1, y=V3, col=logit(1e-5+Vas_Endo+Pericyte),shape= rev(dynamic_class3)) ) + geom_point(size=2.5) + theme_classic(base_size=18) + scale_color_viridis_c()#+ scale_color_npg(name="Tumor response") + scale_size_continuous(name="Sample size")
# ggplot(umapout, aes(x=V1, y=V3, col=logit(1e-5+Pericyte),shape= rev(dynamic_class3)) ) + geom_point(size=2.5) + theme_classic(base_size=18) + scale_color_viridis_c()#+ scale_color_npg(name="Tumor response") + scale_size_continuous(name="Sample size")
# ggplot(umapout, aes(x=V1, y=V3, col=logit(1e-5+Adipocyte),shape= rev(dynamic_class3)) ) + geom_point(size=2.5) + theme_classic(base_size=18) + scale_color_viridis_c()#+ scale_color_npg(name="Tumor response") + scale_size_continuous(name="Sample size")
# 
# 
# ggplot(umapout, aes(x=V1, y=V3, col=logit(1e-5+Dendritic_cell),shape= ARM,size=totalCount) ) + geom_point() + theme_classic(base_size=18) #+ scale_color_npg(name="Tumor response") + scale_size_continuous(name="Sample size")
# ggplot(umapout, aes(x=V1, y=V3, col=logit(1e-5+Fibroblast),shape= ARM,size=totalCount) ) + geom_point() + theme_classic(base_size=18) #+ scale_color_npg(name="Tumor response") + scale_size_continuous(name="Sample size")
# ggplot(umapout, aes(x=V1, y=V3, col=logit(1e-5+Vas_Endo),shape= ARM,size=totalCount) ) + geom_point() + theme_classic(base_size=18) #+ scale_color_npg(name="Tumor response") + scale_size_continuous(name="Sample size")
# 
# # Perform NMDS using a bray curtis dissimilarity metric as a measure of the compostititonal differences between biopsies 
# # As we converted to fractions above, this is equivalent to the Manhattan distance between the sample composition (this has a nice theoretical basis).
# set.seed(1)
# comm.bc.mds <- metaMDS(fracTable2,  dist = "bray",k = 2,
#                        maxit = 999, 
#                        trymax = 500,
#                        wascores = TRUE)
# 
# # Check that the stress is below 0.1 and re-run with higher k if not
# comm.bc.mds$stress
# 
# # Produces a Shepards diagram
# #correlation statistics indicating the fit between ordiantion distances and observed dissimilarities
# stressplot(comm.bc.mds) 
# 
# # gathr output of NMDS
# nmdsout<-data.table( cbind( fracTable1%>%dplyr::select(rowID,Patient.Study.ID:totalCount),  scores(comm.bc.mds) ) ) #,princomp(scores(comm.bc.mds))$scores
# 
# # Plot first 2 dimensions and colour by response group
# ggplot(nmdsout, aes(x=NMDS1, y=NMDS2, col=dynamic_class3,size=totalCount) ) + geom_point()+theme_classic(base_size=18) + scale_color_npg(name="Tumor response") + scale_size_continuous(name="Sample size")
# 
# # Look at what kinds of composition the NMDS axes are associated with
# corvals <- cor(scores(comm.bc.mds), fracTable2 )
# corvals[1,][ order(-abs(corvals[1,]))][1:10]
# corvals[2,][ order(-abs(corvals[2,]))][1:10]
# corvals[3,][ order(-abs(corvals[3,]))][1:10]
# 
# # Which of this small number of compositional variables are best describing response.. maybe we can find out what they represent biologically later
# corResponse <- cor(scores(comm.bc.mds), 1*(nmds1$dynamic_class3=="Response") )
# 
# # Plot some dimensions that look importnat and how the composition changes over time
# ggplot(nmdsout, aes(x= NMDS6, y= NMDS5, col=dynamic_class3, size=samplesize_it)) + 
#   geom_point()+
#   theme_classic(base_size=18)+
#   scale_color_npg(name="Tumor response")+ 
#   scale_size_continuous(name="Sample size") +
#   facet_wrap(~Day) +theme(aspect.ratio=1)
# 
# 
# 
# 
# NMS$stress # Check that the stress is below 0.1 and re-run with higher k if not
# goodness(NMS)
# stressplot(NMS) # Produces a Shepards diagram
# plot(NMS, "sites")   # Produces distance 
# 
# grep("Epithelial",dd0$Annotation[1:13])
