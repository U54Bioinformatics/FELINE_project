rm(list=ls())
require(mgcv);require(data.table);require(dplyr);require(ggplot2);require(tidyr);require(parallel);require(lme4);require(lmerTest);require(parallel)
library(effects)#require(umap)
library(TSP)
require(mclust)
#ttt<-data.table(readRDS(file="~/Dropbox/FELINE Project/Data_analysis/scRNA/15_subclone_interaction/FEL011046_data_gene_pathway.RDS"))
ttt<-data.table(readRDS(file="~/Dropbox/FELINE Project/Data_analysis/scRNA/15_subclone_interaction/FEL011046_data_gene_pathway.v1.RDS"))
ttt2<-data.table(gather(ttt,Variable,Expression,ERBB4:YOSHIMURA_MAPK8_TARGETS_UP))
gene_pathway_lu<-data.table(read.csv(file="~/Dropbox/FELINE Project/Data_analysis/scRNA/15_subclone_interaction/gene and pathway catag AB.csv")%>%dplyr::select(-row.ID))
ttt3 <- merge(ttt2,gene_pathway_lu, by="Variable")
require("VGAM")#lamW")
#possible transformation  
ttt3[,sign_nn := 1* (Expression>0)*1 + -1*( (Expression<=0)*1)]
ttt3[,Expression2:=Expression]
ttt3[VariableType=="pathway",Expression2:=Expression]


signal_recept_tranUnSlice <- na.omit(data.table( ttt3[Day==180][VariableType=="gene"][Major_catagory%in%c("Cell cycle regulator")]%>%
                                                   dplyr::select(-c(VariableType,Negative.regulator,Sub_catagory,Expression2,sign_nn,Major_catagory))%>%spread(Variable,Expression)
))
signal_recept_tranSlice <- na.omit(data.table(signal_recept_tranUnSlice %>% group_by(Patient.Study.ID,ARM, Day,Subclone,Response) %>% dplyr::slice(1:100)))

in_dd<- signal_recept_tranSlice%>%dplyr::select(-c(Cell.ID, Patient.Study.ID, ARM, Day,  Response, Subclone, Phase))
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


signal_recept_tranUnSliceAllt <- na.omit(data.table( ttt3[ARM%in%c("A","B","C")][VariableType=="gene"][Major_catagory%in%c("Cell cycle regulator")]%>%
                                                       dplyr::select(-c(VariableType,Negative.regulator,Sub_catagory,Expression2,sign_nn,Major_catagory))%>%spread(Variable,Expression)
))

project_this<- signal_recept_tranUnSliceAllt%>%dplyr::select(-c(Cell.ID, Patient.Study.ID, ARM, Day,  Response, Subclone, Phase))

rm(list=c("ttt","ttt2"))

u_datFULL_ALLtime <- data.table(signal_recept_tranUnSliceAllt, predict(umap_md3, data=project_this  )   )

corV1<- as.data.table(cor(u_datFULL%>%dplyr::select(-c(Cell.ID, Patient.Study.ID, ARM, Day,  Response, Subclone, Phase,V1:V5)),
                          u_datFULL%>%dplyr::select(V1:V5) ),keep.rownames=TRUE)
corV1[order(-abs(V1))][1:10]
corV1[order(-abs(V5))][1:10]
pcaOverall<- princomp(u_datFULL_ALLtime%>%dplyr::select(V1,V5))
u_datFULL_ALLtimeDims<-data.table(u_datFULL_ALLtime, dimension1=pcaOverall$scores[,1],dimension2=pcaOverall$scores[,2])

ggplot(   u_datFULL_ALLtime, aes(V1, V2, col=Response, shape=Subclone))+geom_point(size=2)+theme_classic()+facet_wrap(Response~Day)
corV2<- as.data.table(cor(u_datFULL_ALLtimeDims%>%dplyr::select(-c(Cell.ID, Patient.Study.ID, ARM, Day,  Response, Subclone, Phase,V1:V5,dimension1:dimension2)),
                          u_datFULL_ALLtimeDims%>%dplyr::select(dimension1:dimension2) ),keep.rownames=TRUE)
corV2[order(-abs(dimension1))][1:10]
corV2[order(-abs(dimension2))][1:10]

ggplot(  u_datFULL_ALLtimeDims , aes(dimension1,dimension2,col=Patient.Study.ID,shape=Subclone))+geom_point(size=2)+facet_wrap(Response~Day)+theme_classic()
ggplot(  u_datFULL_ALLtimeDims , aes(dimension1,dimension2,col=Response))+geom_point(size=0.5)+facet_wrap(Response~Day)+theme_classic()


require(mclust)
u_datFULL_ALLtimeDims$classification<- NULL


# cluster
set.seed(1234)
mclust.mod<-Mclust(u_datFULL_ALLtimeDims[ARM=="B"]%>%  dplyr::select(V1:V2),G=30);
set.seed(1234)
mclust.mod2<-Mclust(u_datFULL_ALLtimeDims[ARM=="B"]%>%  dplyr::select(V1:V2),G=4);
u_datFULL_ALLtimeDims$classification <- predict(mclust.mod,newdata=u_datFULL_ALLtimeDims%>%  dplyr::select(V1:V2))$classification # mclust.mod$classification
u_datFULL_ALLtimeDims$classification2 <- predict(mclust.mod2,newdata=u_datFULL_ALLtimeDims%>%  dplyr::select(V1:V2))$classification # mclust.mod$classification
#mclust.mod<-Mclust(u_datFULL%>%  dplyr::select(V1:V2),G=20)
# Solve traveling salesman problem
library(TSP)
tsp <- TSP(dist(t(mclust.mod$parameters$mean), method = "euclidean")); tour <- solve_TSP(tsp)

tdf0 <- data.table(tour=1:max(mclust.mod$classification),classification=as.integer(tour))[order(classification)]
tdf  <- data.table(tdf0, t(mclust.mod$parameters$mean))[order(tour)]
#tdf <- data.table(tour,classification=1:max(mclust.mod$classification), t(mclust.mod$parameters$mean))[order(tour)]
tdf <- data.table(tdf%>%mutate( V1e=lead(V1),         V2e=lead(V2)))
tdf[is.na(V1e),V1e:=tdf[tour==1]$V1]
tdf[is.na(V2e),V2e:=tdf[tour==1]$V2]



ggplot(  merge(u_datFULL_ALLtimeDims ,tdf%>%dplyr::select( tour ,classification),by=c("classification")) , aes(V1,V2))+
  geom_point(aes(col=as.factor(tour)),size=2)+facet_wrap(Response~Day)+theme_classic()+facet_wrap(~Patient.Study.ID)+
  geom_point(data=tdf, size=3, aes(x=V1  ,      y=  V2,col=as.factor(tour)))+geom_segment(data=tdf,aes(x=V1  ,xend= V1e,     y=  V2,yend= V2e ),col="black")+
  scale_color_discrete(name="cell cycle stage")#,col=as.factor(tour)


ggplot(  merge(u_datFULL_ALLtimeDims ,tdf%>%dplyr::select( tour ,classification),by=c("classification")) , aes(V1,V2))+
  geom_point(aes(col=as.factor(classification2)),size=2)+theme_classic()

ggplot(  merge(u_datFULL_ALLtimeDims ,tdf%>%dplyr::select( tour ,classification),by=c("classification")) , aes(V1,V2))+
  geom_point(aes(col=as.factor(Response)),size=2)+theme_classic()+facet_wrap(Day~ARM)

ggplot(  merge(u_datFULL_ALLtimeDims ,tdf%>%dplyr::select( tour ,classification),by=c("classification")) , aes(V1,V2))+
  geom_point(aes(col=as.factor(tour)),size=2)+#facet_wrap(Response~Day)+theme_classic()+facet_wrap(Day~ARM)+
  geom_point(data=tdf, size=3, aes(x=V1  ,      y=  V2,col=as.factor(tour)))+geom_segment(data=tdf,aes(x=V1  ,xend= V1e,     y=  V2,yend= V2e ),col="black")+
  scale_color_discrete(name="cell cycle stage")#,col=as.factor(tour)





ggplot(  merge(u_datFULL_ALLtimeDims ,tdf%>%dplyr::select( tour ,classification),by=c("classification")) , aes(y=CDK6,x=tour))+
  geom_point(aes(col=as.factor(tour)),size=2)+
  coord_polar()+
  #facet_wrap(Response~Day)+
  theme_classic()+
  facet_wrap(~Patient.Study.ID)

ggplot(  merge(u_datFULL_ALLtimeDims ,tdf%>%dplyr::select( tour ,classification),by=c("classification")) , aes(y=CDK6,x=tour))+
  geom_point(aes(col=as.factor(Response)),size=2)+
  coord_polar()+
  facet_wrap(ARM~Day)+
  theme_classic()
ggplot(  merge(u_datFULL_ALLtimeDims ,tdf%>%dplyr::select( tour ,classification),by=c("classification")) , aes(y=CCND1,x=tour))+
  geom_point(aes(col=as.factor(Response)),size=2)+
  coord_polar()+
  facet_wrap(ARM~Day)+
  theme_classic()
ggplot(  merge(u_datFULL_ALLtimeDims ,tdf%>%dplyr::select( tour ,classification),by=c("classification")) , aes(y=CDKN1B,x=tour))+
  geom_point(aes(col=as.factor(Response)),size=2)+
  coord_polar()+
  facet_wrap(ARM~Day)+
  theme_classic()


ggplot(  merge(u_datFULL_ALLtimeDims ,tdf%>%dplyr::select( tour ,classification),by=c("classification")) , aes(y=CDK6,x=tour))+
  geom_point(aes(col=as.factor(tour)),size=2)+
  coord_polar()+
  facet_wrap(Response~Day)+
  theme_classic()+
  geom_smooth(aes(group=Response,col=Response),method="gam",formula= y~s(x,bs="cc",k=4))#+#


u_datFULL_ALLtimeDims2<- merge(u_datFULL_ALLtimeDims ,tdf%>%dplyr::select( tour ,classification),by=c("classification"))
u_datFULL_ALLtimeDims2[, classification3:= "G1"]
u_datFULL_ALLtimeDims2[tour<(0.24*30),classification3:= "G0"]
u_datFULL_ALLtimeDims2[tour>(0.7*30),classification3:= "S/G2"]

# save(u_datFULL_ALLtimeDims,u_datFULL_ALLtimeDims2,tsp,mclust.mod,mclust.mod2,tdf, file="/Users/jason/Dropbox/FELINE Project/Data_analysis/scRNA/15_subclone_interaction/Dim reduced signalling and transduction.RData")
load( file="/Users/jason/Dropbox/FELINE Project/Data_analysis/scRNA/15_subclone_interaction/Dim reduced signalling and transduction.RData")

xyz<- spread(ttt3[Variable%in%c("ESR1","MAPK10","MAPK8","FGFR2","FOS","JUNB")]%>%dplyr::select(Cell.ID,Variable,  Expression),Variable,  Expression)

ggplot( merge(u_datFULL_ALLtimeDims2[ARM=="C"][Day==180],xyz,by="Cell.ID")%>%group_by(ARM,Day,tour,classification3,Response,Patient.Study.ID,Subclone)%>%dplyr::summarise(CDK6=mean(CDK6))
        
        
        , aes(y=CDK6,x=tour/max(tour)))+
  geom_point(aes(fill=as.factor(classification3)),pch=21,size=2)+
  coord_polar()+
  facet_wrap(Patient.Study.ID~Day)+
  theme_minimal()+
  geom_smooth(aes(group=Response,col=Response),size=1.5,method="gam",formula= y~s(x,bs="cc",k=5),se=T)+scale_x_continuous(breaks=seq(0,1,by=0.1))

ggplot( merge(u_datFULL_ALLtimeDims2,xyz,by="Cell.ID")%>%group_by(ARM,Day,tour,classification3,Response,Patient.Study.ID,Subclone)%>%dplyr::summarise(CDK6=mean(CDK6))
        , aes(y=CDK6,x=tour/max(tour)))+
  #geom_point(aes(fill=as.factor(classification3)),pch=21,size=2)+
  coord_polar()+
  facet_wrap(ARM~Day)+
  theme_minimal()+
  geom_smooth(aes(group=Response,col=Response),size=1.5,method="gam",formula= y~s(x,bs="cc",k=5),se=T)+scale_x_continuous(breaks=seq(0,1,by=0.1))

ggplot( merge(u_datFULL_ALLtimeDims2,xyz,by="Cell.ID")%>%group_by(ARM,Day,tour,classification3,Response,Patient.Study.ID,Subclone)%>%dplyr::summarise(MAPK8=mean(MAPK8),MAPK10=mean(MAPK10),CDK6=mean(CDK6),ESR1=mean(ESR1) )
        , aes(y=ESR1,x=tour/max(tour)))+
  # geom_point(aes(fill=as.factor(classification3)),pch=21,size=2)+
  coord_polar()+
  facet_wrap(ARM~Day)+
  theme_minimal()+
  geom_smooth(aes(group=Response,col=Response),size=1.5,method="gam",formula= y~s(x,bs="cc",k=5),se=T)+scale_x_continuous(breaks=seq(0,1,by=0.1))

ggplot( merge(u_datFULL_ALLtimeDims2,xyz,by="Cell.ID")%>%group_by(ARM,Day,tour,classification3,Response,Patient.Study.ID,Subclone)%>%dplyr::summarise(MAPK8=mean(MAPK8),MAPK10=mean(MAPK10),CDK6=mean(CDK6),ESR1=mean(ESR1) )
        , aes(y=MAPK8,x=tour/max(tour)))+  coord_polar()+
  facet_wrap(ARM~Day)+
  theme_minimal()+
  geom_smooth(aes(group=Response,col=Response),size=1.5,method="gam",formula= y~s(x,bs="cc",k=7),se=T)+scale_x_continuous(breaks=seq(0,1,by=0.1))

# The palette with black:
cbbPalette <- c(  "#56B4E9", "#E69F00","#D55E00")#   "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
pa<-ggplot( merge(u_datFULL_ALLtimeDims2,xyz,by="Cell.ID")%>%group_by(ARM,Day,tour,classification3,Response,Patient.Study.ID,Subclone)%>%dplyr::summarise(MAPK8=mean(MAPK8),MAPK10=mean(MAPK10),CDK6=mean(CDK6),ESR1=mean(ESR1) ), 
            aes(y=MAPK8,x=(tour-1)/max(tour-1)))+  coord_polar()+
  facet_wrap(ARM~.)+
  theme_minimal()+
  geom_smooth(aes(group=as.factor(Day),col=as.factor(Day),linetype=as.factor(Day) ),size=2.5,method="gam",formula= y~s(x,bs="cc",k=5),se=F)+
  scale_color_manual(name="Day",values=cbbPalette[1:3])+labs(y="",x="")+
  #scale_color_brewer(name="Day",palette="Set2")+labs(y="",x="")+
  scale_linetype_discrete(name="Day")+theme(aspect.ratio = 1)+
  theme(
    legend.position="none",
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+
  # geom_hline(yintercept=logit(1e-3+1/ 4.3 ))
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())+scale_x_continuous(breaks=seq(0,1,by=0.1))+xlim(0,1)
ggsave(pa,filename = "/Users/jason/Dropbox/FELINE Project/Manuscript/Final Main Figures/Blank figures/JNK1 cell cycle expression by arm and time.png", width = 12, height = 4)


pb<-ggplot( merge(u_datFULL_ALLtimeDims2,xyz,by="Cell.ID")%>%group_by(ARM,Day,tour,classification3,Response,Patient.Study.ID,Subclone)%>%dplyr::summarise(MAPK8=mean(MAPK8),MAPK10=mean(MAPK10),CDK6=mean(CDK6),ESR1=mean(ESR1) ,FOS=mean(FOS),JUNB=mean(JUNB)), 
            aes(y=FOS,x=(tour-1)/max(tour-1)))+  coord_polar()+
  facet_wrap(ARM~.)+
  theme_minimal()+
  geom_smooth(aes(group=as.factor(Day),col=as.factor(Day),linetype=as.factor(Day) ),size=2.5,method="gam",formula= y~s(x,bs="cc",k=5),se=F)+
  scale_color_manual(name="Day",values=cbbPalette[1:3])+labs(y="",x="")+
  #scale_color_brewer(name="Day",palette="Set2")+labs(y="",x="")+
  scale_linetype_discrete(name="Day")+theme(aspect.ratio = 1)+
  theme(
    legend.position="none",
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+
  # geom_hline(yintercept=logit(1e-3+1/ 4.3 ))
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())+scale_x_continuous(breaks=seq(0,1,by=0.1))+xlim(0,1)
ggsave(pb,filename = "/Users/jason/Dropbox/FELINE Project/Manuscript/Final Main Figures/Blank figures/FOS cell cycle expression by arm and time.png", width = 12, height = 4)



pb<- ggplot( merge(u_datFULL_ALLtimeDims2,xyz,by="Cell.ID")%>%group_by(ARM,Day,tour,classification3,Response,Patient.Study.ID,Subclone)%>%dplyr::summarise(MAPK8=mean(MAPK8),MAPK10=mean(MAPK10),CDK6=mean(CDK6),ESR1=mean(ESR1) ,FOS=median(FOS),JUNB=mean(JUNB)), 
             aes(y=JUNB,x=(tour-1)/max(tour-1)))+  coord_polar()+
  facet_wrap(ARM~.)+
  theme_minimal()+
  geom_smooth(aes(group=as.factor(Day),col=as.factor(Day),linetype=as.factor(Day) ),size=2.5,method="gam",formula= y~s(x,bs="cc",k=4),se=F)+
  scale_color_manual(name="Day",values=cbbPalette[1:3])+labs(y="",x="")+
  #scale_color_brewer(name="Day",palette="Set2")+labs(y="",x="")+
  scale_linetype_discrete(name="Day")+theme(aspect.ratio = 1)+
  theme(
    legend.position="none",
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+
  # geom_hline(yintercept=logit(1e-3+1/ 4.3 ))
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())+scale_x_continuous(breaks=seq(0,1,by=0.1))+xlim(0,1)+ylim(c(-0.22,0.2))
ggsave(pb,filename = "/Users/jason/Dropbox/FELINE Project/Manuscript/Final Main Figures/Blank figures/JUNB cell cycle expression by arm and time.png", width = 12, height = 4)

ggplot( merge(u_datFULL_ALLtimeDims2,xyz,by="Cell.ID")%>%group_by(ARM,Day,tour,classification3)%>%dplyr::summarise(MAPK8=mean(MAPK8),MAPK10=mean(MAPK10),CDK6=mean(CDK6),ESR1=mean(ESR1) ), 
        aes(y=ESR1,x=tour/max(tour)))+  coord_polar()+
  facet_wrap(ARM~.)+ylim(c(-5,3))+
  geom_smooth(aes(group=as.factor(Day),col=as.factor(Day),linetype=as.factor(Day) ),size=2.5,method="gam",formula= y~s(x,bs="cc",k=5),se=F)


pc<- ggplot( merge(u_datFULL_ALLtimeDims2,xyz,by="Cell.ID")%>%group_by(ARM,Day,tour,classification3)%>%dplyr::summarise(MAPK8=mean(MAPK8),MAPK10=mean(MAPK10),CDK6=mean(CDK6),ESR1=mean(ESR1) ), 
             aes(y=ESR1,x=(tour-1)/max(tour-1)))+  coord_polar()+
  facet_wrap(ARM~.)+
  theme_minimal()+
  geom_smooth(aes(group=as.factor(Day),col=as.factor(Day),linetype=as.factor(Day) ),size=2.5,method="gam",formula= y~s(x,bs="cc",k=5),se=F)+
  scale_color_manual(name="Day",values=cbbPalette[1:3])+labs(y="",x="")+
  scale_linetype_discrete(name="Day")+theme(aspect.ratio = 1)+
  theme(
    legend.position="none",
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+
  # geom_hline(yintercept=logit(1e-3+1/ 4.3 ))
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())+scale_x_continuous(breaks=seq(0.1,1,by=0.1))+xlim(0,1)+ylim(c(-6,3.95))
ggsave(pc,filename = "/Users/jason/Dropbox/FELINE Project/Manuscript/Final Main Figures/Blank figures/ESR1 cell cycle expression by arm and time.png", width = 12, height = 4)

pd<-ggplot( merge(u_datFULL_ALLtimeDims2,xyz,by="Cell.ID")%>%group_by(ARM,Day,tour,classification3,Response,Patient.Study.ID,Subclone)%>%dplyr::summarise(MAPK8=mean(MAPK8),MAPK10=mean(MAPK10),CDK6=mean(CDK6),ESR1=mean(ESR1) ), 
            aes(y=CDK6,x=(tour-1)/max(tour-1)))+  coord_polar()+
  facet_wrap(ARM~.)+
  theme_minimal()+
  geom_smooth(aes(group=as.factor(Day),col=as.factor(Day),linetype=as.factor(Day) ),size=2.5,method="gam",formula= y~s(x,bs="cc",k=5),se=F)+
  scale_color_manual(name="Day",values=cbbPalette[1:3])+labs(y="",x="")+
  scale_linetype_discrete(name="Day")+theme(aspect.ratio = 1)+
  theme(
    legend.position="none",
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+
  # geom_hline(yintercept=logit(1e-3+1/ 4.3 ))
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())+scale_x_continuous(breaks=seq(0.1,1,by=0.1))+xlim(0,1)
ggsave(pd,filename = "/Users/jason/Dropbox/FELINE Project/Manuscript/Final Main Figures/Blank figures/CDK6 cell cycle expression by arm and time.png", width = 12, height = 4)



ggplot( merge(u_datFULL_ALLtimeDims2,xyz,by="Cell.ID")%>%group_by(ARM,Day,tour,classification3,Response,Patient.Study.ID,Subclone)%>%dplyr::summarise(MAPK8=mean(MAPK8),MAPK10=mean(MAPK10),CDK6=mean(CDK6),ESR1=mean(ESR1) ,CDKN1B=median(CDKN1B),CDKN2A=mean(CDKN2A)), 
        aes(y=CDKN1B,x=(tour-1)/max(tour-1)))+  coord_polar()+
  facet_wrap(ARM~.)+
  theme_minimal()+
  geom_smooth(aes(group=as.factor(Day),col=as.factor(Day),linetype=as.factor(Day) ),size=2.5,method="gam",formula= y~s(x,bs="cc",k=5),se=F)+
  scale_color_manual(name="Day",values=cbbPalette[1:3])+labs(y="",x="")+
  scale_linetype_discrete(name="Day")+theme(aspect.ratio = 1)+
  theme(
    legend.position="none",
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+
  # geom_hline(yintercept=logit(1e-3+1/ 4.3 ))
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())+scale_x_continuous(breaks=seq(0.1,1,by=0.1))+xlim(0,1)



pe<-ggplot( merge(u_datFULL_ALLtimeDims2,xyz,by="Cell.ID")%>%group_by(ARM,Day,tour,classification3,Response,Patient.Study.ID,Subclone)%>%dplyr::summarise(MAPK8=mean(MAPK8),MAPK10=mean(MAPK10),CDK6=mean(CDK6),ESR1=mean(ESR1) ,CDKN1B=median(CDKN1B),CDKN2A=mean(CDKN2A)), 
            aes(y=CDKN2A,x=(tour-1)/max(tour-1)))+  coord_polar()+
  facet_wrap(ARM~.)+
  theme_minimal()+
  geom_smooth(aes(group=as.factor(Day),col=as.factor(Day),linetype=as.factor(Day) ),size=2.5,method="gam",formula= y~s(x,bs="cc",k=5),se=F)+
  scale_color_manual(name="Day",values=cbbPalette[1:3])+labs(y="",x="")+
  scale_linetype_discrete(name="Day")+theme(aspect.ratio = 1)+
  theme(
    legend.position="none",
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+
  # geom_hline(yintercept=logit(1e-3+1/ 4.3 ))
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())+scale_x_continuous(breaks=seq(0.1,1,by=0.1))+xlim(0,1)

ggsave(pe,filename = "/Users/jason/Dropbox/FELINE Project/Manuscript/Final Main Figures/Blank figures/CDKN2A cell cycle expression by arm and time.png", width = 12, height = 4)

ggplot( merge(u_datFULL_ALLtimeDims2,xyz,by="Cell.ID")%>%group_by(ARM,Day,tour,classification3,Response,Patient.Study.ID,Subclone)%>%dplyr::summarise(MAPK8=mean(MAPK8),MAPK10=mean(MAPK10),CDK6=mean(CDK6),ESR1=mean(ESR1) ,CDKN1B=median(CDKN1B),CDKN2A=mean(CDKN2A)), 
        aes(y=CDKN2A,x=(tour-1)/max(tour-1)))+  coord_polar()+
  facet_wrap(ARM~.)+
  theme_minimal()+
  geom_smooth(aes(group=as.factor(Day),col=as.factor(Day),linetype=as.factor(Day) ),size=2.5,method="gam",formula= y~s(x,bs="cc",k=5),method.args=list(gamma=1.25),se=T)+
  scale_color_manual(name="Day",values=cbbPalette[1:3])+labs(y="",x="")+
  scale_linetype_discrete(name="Day")+theme(aspect.ratio = 1)+
  theme(
    legend.position="none",
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+
  # geom_hline(yintercept=logit(1e-3+1/ 4.3 ))
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())+scale_x_continuous(breaks=seq(0.1,1,by=0.1))+xlim(0,1)

gamddA<-data.table(merge(u_datFULL_ALLtimeDims2,xyz,by="Cell.ID")%>%group_by(ARM,Day,tour,classification3,Response,Patient.Study.ID,Subclone)%>%dplyr::summarise(MAPK8=mean(MAPK8),MAPK10=mean(MAPK10),CDK6=mean(CDK6),ESR1=mean(ESR1) ,CDKN1B=median(CDKN1B),CDKN2A=mean(CDKN2A))
)[ARM=="A"]
gamddA[, x:=(tour-1)/max(tour-1)]
gamA <- gam(CDKN2A~as.factor(Day) ,data=gamddA)
gamA <- gam(CDKN2A~as.factor(Day)+s(x,by=as.factor(Day) ,k=5, bs="cs"),data=gamddA)
#plot(gam1)
summary(gamA)#[4]

gamddB<-data.table(merge(u_datFULL_ALLtimeDims2,xyz,by="Cell.ID")%>%group_by(ARM,Day,tour,classification3,Response,Patient.Study.ID,Subclone)%>%dplyr::summarise(MAPK8=mean(MAPK8),MAPK10=mean(MAPK10),CDK6=mean(CDK6),ESR1=mean(ESR1) ,CDKN1B=median(CDKN1B),CDKN2A=mean(CDKN2A))
)[ARM=="B"]
gamddB[, x:=(tour-1)/max(tour-1)]
gamB <- gam(CDKN2A~as.factor(Day) ,data=gamddB)
gamB <- gam(CDKN2A~as.factor(Day)+s(x,by=as.factor(Day) ,k=5, bs="cs"),data=gamddB)
#plot(gam1)
summary(gamB)


gamddA<-data.table(merge(u_datFULL_ALLtimeDims2,xyz,by="Cell.ID")%>%group_by(ARM,Day,tour,classification3,Response,Patient.Study.ID,Subclone)%>%dplyr::summarise(MAPK8=mean(MAPK8),MAPK10=mean(MAPK10),CDK6=mean(CDK6),ESR1=mean(ESR1) ,CDKN1B=median(CDKN1B),CDKN2A=mean(CDKN2A))
)[ARM=="A"]
gamddA[, x:=(tour-1)/max(tour-1)]
gamA <- gam(CDKN2A~as.factor(Day) ,data=gamddA)
gamA <- gam(CDKN2A~as.factor(Day)+s(x,by=as.factor(Day) ,k=5, bs="cs"),data=gamddA)
#plot(gam1)
summary(gamA)
gam1 <- gam(CDKN2A~as.factor(Day)+s(x,by=as.factor(Day) ,k=5, bs="cs"),data=gamdd)
#plot(gam1)
summary(gam1)#[4]


str(summary(gam1))
gamdd$pred<-predict(gam1)

gam1 <- gam(CDKN2A~as.factor(Day)+s(x,by=as.factor(Day) ,k=5, bs="cs"),data=gamdd)


pf <-ggplot( merge(u_datFULL_ALLtimeDims2,xyz,by="Cell.ID")%>%group_by(ARM,Day,tour,classification3,Response,Patient.Study.ID,Subclone)%>%dplyr::summarise(MAPK8=mean(MAPK8),MAPK10=mean(MAPK10),CDK6=mean(CDK6),ESR1=mean(ESR1) ,CDKN1B=median(CDKN1B),CDKN2A=mean(CDKN2A),RB1=mean(RB1)), 
             aes(y=RB1,x=(tour-1)/max(tour-1)))+  coord_polar()+
  facet_wrap(ARM~.)+
  theme_minimal()+
  geom_smooth(aes(group=as.factor(Day),col=as.factor(Day),linetype=as.factor(Day) ),size=2.5,method="gam",formula= y~s(x,bs="cc",k=5),se=F)+
  scale_color_manual(name="Day",values=cbbPalette[1:3])+labs(y="",x="")+
  scale_linetype_discrete(name="Day")+theme(aspect.ratio = 1)+
  theme(
    legend.position="none",
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+
  # geom_hline(yintercept=logit(1e-3+1/ 4.3 ))
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())+scale_x_continuous(breaks=seq(0.1,1,by=0.1))+xlim(0,1)
ggsave(pf,filename = "/Users/jason/Dropbox/FELINE Project/Manuscript/Final Main Figures/Blank figures/RB1 cell cycle expression by arm and time.png", width = 12, height = 4)








tempn <-data.table( merge(u_datFULL_ALLtimeDims2,xyz,by="Cell.ID")[Day==0]%>%group_by(Patient.Study.ID)%>%dplyr::summarise(mu_ESR1=mean(ESR1) ) )

esr1_dd<-data.table(merge( merge(u_datFULL_ALLtimeDims2,xyz,by="Cell.ID"),tempn,by= "Patient.Study.ID")%>%mutate(R=1*(Response=="Responder"),N=1*(Response!="Responder")))
esr1_dd[,y:=ESR1-mu_ESR1]
esr1_dd[,lnDay:=log(1+Day)]
esr1_dd[,ARMResponse:=as.factor(paste(ARM,Response,sep="_"))]
esr1_dd[,FacPatient.Study.ID:=as.factor(Patient.Study.ID)]
esr1_dd[,FacARM:=as.factor(ARM)]

gammod<-gam(y~s(lnDay, by=FacARM, k=3)+s(FacPatient.Study.ID,bs="re"),gamma=0.1,
            data=esr1_dd)
summary(gammod)
lmmod<-lmer(y~lnDay*FacARM+(1|FacPatient.Study.ID),
            data=esr1_dd)
summary(lmmod)
esr1_dd%>%dplyr(ARM,Response,)
newy <- data.table(expand.grid(ARM=unique(esr1_dd$ARM),
                               #Response=unique(esr1_dd$Response),
                               lnDay=seq(min(esr1_dd$lnDay),max(esr1_dd$lnDay),length=100)  ,Patient.Study.ID=NA))
newy[,ARMResponse:=as.factor(paste(ARM,Response,sep="_"))]
newy[,FacPatient.Study.ID:=as.factor(Patient.Study.ID)]
newy[,FacARM:=as.factor(ARM)]
newy$Pred <- predict(gammod,newdata=newy,exclude="s(FacPatient.Study.ID)",newdata.guaranteed=TRUE)
tempd <- newy[lnDay==min(lnDay)]%>%mutate(Pred0=Pred)%>%dplyr::select(-c(Pred,lnDay))
newy2 <- merge(newy,tempd,by=names(tempd)[names(tempd)%in%names(newy)])
newy2$ARM<- as.character(newy2$ARM)
esr1_dd$ARM<- as.character(esr1_dd$ARM)

esr1_dd[(ESR1-esr1_dd) <quantile(ESR1-mu_ESR1,probs = 0.9) & (ESR1-esr1_dd) >quantile(ESR1-mu_ESR1,probs = 0.1),by=c("ARM","lnDay")]

esr1_dd[,uq:=quantile(ESR1-mu_ESR1,probs = 0.9) ,by=c("ARM","lnDay")]
esr1_dd[,lq:=quantile(ESR1-mu_ESR1,probs = 0.1) ,by=c("ARM","lnDay")]


esr1_dd[ (ESR1-mu_ESR1) > lq & (ESR1-mu_ESR1) <uq]
#erplot <- ggplot( newy2,aes(y=Pred-Pred0,x=lnDay,group=Response,col=Response ))+  geom_line(size=2)+facet_wrap(~ARM)
erplotN <- ggplot( newy2,aes(y=Pred-Pred0,x=exp(lnDay)-1,group=ARM,col=ARM ))+  geom_line(size=4) + theme_classic(base_size=20)+
  #geom_point(data=esr1_dd,aes(y=ESR1-mu_ESR1,x=Day))+
  geom_violin(data=esr1_dd[ (ESR1-mu_ESR1) > lq & (ESR1-mu_ESR1) <uq],aes(y=ESR1-mu_ESR1,x=exp(lnDay),fill=ARM,col=ARM,group=interaction(ARM,Day)),alpha=0.5,trim=TRUE) +
  labs(y="ESR1 expression", x="Day")
ggsave(erplotN,filename = "/Users/jason/Desktop/ESR1 expression over time by arm.pdf", width = 9, height = 9)

erplotN2 <- ggplot( newy2,aes(y=Pred-Pred0,x=(lnDay)-1,group=ARM,col=ARM ))+  geom_line(size=4) + theme_classic(base_size=20)+
  #geom_point(data=esr1_dd,aes(y=ESR1-mu_ESR1,x=Day))+
  geom_violin(data=esr1_dd[ (ESR1-mu_ESR1) > lq & (ESR1-mu_ESR1) <uq],aes(y=ESR1-mu_ESR1,x=(lnDay),fill=ARM,col=ARM,group=interaction(ARM,Day)),alpha=0.5,trim=TRUE) +
  labs(y="ESR1 expression", x="Day")+ scale_x_continuous(breaks=log(1+c(0,14,180)) ,labels=c(0,14,180))
ggsave(erplotN2,filename = "/Users/jason/Desktop/ESR1 expression over time by arm logday.pdf", width = 9, height = 9)

erplotNB <- ggplot( newy2,aes(y=Pred-Pred0,x=(1+lnDay),group=ARM,col=ARM ))+  geom_line(size=4) + theme_classic(base_size=20)+
  #geom_point(data=esr1_dd,aes(y=ESR1-mu_ESR1,x=Day))+
  geom_violin(data=esr1_dd[ (ESR1-mu_ESR1) > lq & (ESR1-mu_ESR1) <uq],aes(y=ESR1-mu_ESR1,x=(1+lnDay),fill=ARM,col=ARM,group=interaction(ARM,Day)),alpha=0.5,trim=TRUE) +
  labs(y="ESR1 expression", x="Day")+ scale_x_continuous(breaks=log(1+c(0,14,180)) ,labels=c(0,14,180))+
  theme(
    legend.position="none",
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+
  # geom_hline(yintercept=logit(1e-3+1/ 4.3 ))
  theme(axis.title=element_blank(),axis.text=element_blank())+facet_wrap(~ARM)
#ggsave(erplotNB,filename = "/Users/jason/Desktop/BLANK ESR1 expression over time by arm logday.pdf", width = 9, height = 9)
ggsave(erplotNB,filename = "/Users/jason/Desktop/BLANK ESR1 expression over time by arm logday by arm.pdf", width = 15, height = 5)






erplot <- ggplot( data.table(merge( merge(u_datFULL_ALLtimeDims2,xyz,by="Cell.ID"),tempn,by= "Patient.Study.ID")%>%mutate(R=1*(Response=="Responder"),N=1*(Response!="Responder"))), 
                  aes(y=ESR1-mu_ESR1,x=0.5*R+log(1+Day),group=Response,col=Response ))+  
  facet_wrap(~ARM)+
  theme_classic(base_size=18)+
  geom_violin(aes(fill=Response,group=interaction(Day,Response)),alpha=0.5)+
  geom_jitter(aes(col=Response,group=Response),size=0.35,position = position_jitter(width = .05))+#scale_color_viridis_c(option="B")+
  scale_x_continuous(name="Day",breaks=log(1+c(0,14,180)),labels=c(0,14,180))+
  ylab("ESR1 expression")+
  geom_smooth(method="gam",formula= y~s(x,k=3),se=T,alpha=0.5)+theme(aspect.ratio=1)

ggsave(erplot,filename = "/Users/jason/Desktop/ESR1 expression over time by arm and response.png", width = 21, height = 7)



ggplot( merge( merge(u_datFULL_ALLtimeDims2,xyz,by="Cell.ID"),tempn,by= "Patient.Study.ID")%>%mutate(R=1*(Response=="Responder"),N=1*(Response!="Responder")), 
        aes(y=ESR1-mu_ESR1,x=R+log(1+Day)))+  
  facet_wrap(Response~ARM)+
  theme_classic()+
  geom_violin(aes(group=interaction(Day)))+
  geom_jitter(aes(col=ESR1-mu_ESR1),width=0.25,size=0.5)+#scale_color_viridis_c(option="B")+
  scale_x_continuous(name="ESR1 expression",breaks=log(1+c(0,14,180)),labels=c(0,14,180))+
  
  geom_smooth(method="lm",aes(fill=Response ))

