rm(list=ls())
require(mgcv);require(data.table);require(dplyr);require(ggplot2);require(tidyr);require(parallel);require(lme4);require(lmerTest);require(parallel)
library(effects)#require(umap)
library(TSP)
require(mclust)
ttt<-data.table(readRDS(file="~/Dropbox/FELINE Project/Data_analysis/scRNA/15_subclone_interaction/FEL011046_data_gene_pathway.v1.RDS"))
ttt2<-data.table(gather(ttt,Variable,Expression,ERBB4:YOSHIMURA_MAPK8_TARGETS_UP))
gene_pathway_lu<-data.table(read.csv(file="~/Dropbox/FELINE Project/Data_analysis/scRNA/15_subclone_interaction/gene and pathway catag AB.csv")%>%dplyr::select(-row.ID))
ttt3 <- merge(ttt2,gene_pathway_lu, by="Variable")
require("VGAM")#lamW")
ttt3[,sign_nn := 1* (Expression>0)*1 + -1*( (Expression<=0)*1)]
ttt3[,Expression2:=Expression]
ttt3[VariableType=="pathway",Expression2:=Expression]


signal_recept_tranUnSlice <- na.omit(data.table( ttt3[Day==180][VariableType=="gene"][Major_catagory%in%c("Cell cycle regulator")]%>%
                                                   dplyr::select(-c(VariableType,Negative.regulator,Sub_catagory,Expression2,sign_nn,Major_catagory))%>%spread(Variable,Expression)
))
signal_recept_tranSlice <- na.omit(data.table(signal_recept_tranUnSlice %>% group_by(Patient.Study.ID,ARM, Day,Subclone,Response) %>% dplyr::slice(1:100)))

in_dd<- signal_recept_tranSlice%>%dplyr::select(-c(Cell.ID, Patient.Study.ID, ARM, Day,  Response, Subclone, Phase))
umap_md3 <- umap(in_dd, n_neighbors=25, n_components=5)
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
#set.seed(1234)
#mclust.mod2<-Mclust(u_datFULL_ALLtimeDims[ARM=="B"]%>%  dplyr::select(V1:V2),G=4);
#mclust.mod<-Mclust(u_datFULL_ALLtimeDims%>%  dplyr::select(V1:V2),G=20);
#mclust.mod<-Mclust(u_datFULL_ALLtimeDims[ARM=="B"]%>%  dplyr::select(V1:V2),G=100);
#mclust.mod<-Mclust(u_datFULL_ALLtimeDims[Day==180]%>%  dplyr::select(V1:V2),G=20);
u_datFULL_ALLtimeDims$classification <- predict(mclust.mod,newdata=u_datFULL_ALLtimeDims%>%  dplyr::select(V1:V2))$classification # mclust.mod$classification
u_datFULL_ALLtimeDims$classification2 <- predict(mclust.mod2,newdata=u_datFULL_ALLtimeDims%>%  dplyr::select(V1:V2))$classification # mclust.mod$classification
#mclust.mod<-Mclust(u_datFULL%>%  dplyr::select(V1:V2),G=20)
# Solve traveling salesman problem
library(TSP)
tsp <- TSP(dist(t(mclust.mod$parameters$mean), method = "euclidean")); tour <- solve_TSP(tsp)

tdf0 <- data.table(tour=1:max(mclust.mod$classification),classification=as.integer(tour))[order(classification)]
tdf  <- data.table(tdf0, t(mclust.mod$parameters$mean))[order(tour)]
tdf <- data.table(tdf%>%mutate( V1e=lead(V1),         V2e=lead(V2)))
tdf[is.na(V1e),V1e:=tdf[tour==1]$V1]
tdf[is.na(V2e),V2e:=tdf[tour==1]$V2]



ggplot(  merge(u_datFULL_ALLtimeDims ,tdf%>%dplyr::select( tour ,classification),by=c("classification")) , aes(V1,V2))+
  geom_point(aes(col=as.factor(tour)),size=2)+facet_wrap(Response~Day)+theme_classic()+facet_wrap(~Patient.Study.ID)+
  geom_point(data=tdf, size=3, aes(x=V1  ,      y=  V2,col=as.factor(tour)))+geom_segment(data=tdf,aes(x=V1  ,xend= V1e,     y=  V2,yend= V2e ),col="black")+
  scale_color_discrete(name="cell cycle stage")+
  labs(y="Cell cycle Umap dimension 2", x="Cell cycle Umap dimension 1")#,col=as.factor(tour)

ggplot(  merge(u_datFULL_ALLtimeDims ,tdf%>%dplyr::select( tour ,classification),by=c("classification")) , aes(V1,V2))+
  geom_point(aes(col=as.factor(tour)),size=0.1)+theme_classic()+
  geom_point(data=tdf, size=3, aes(x=V1  ,      y=  V2,col=as.factor(tour)))+geom_segment(data=tdf,aes(x=V1  ,xend= V1e,     y=  V2,yend= V2e ),col="black")+
  scale_color_discrete(name="cell cycle stage")+
  labs(y="Cell cycle Umap dimension 2", x="Cell cycle Umap dimension 1")+theme(aspect.ratio=1)#,col=as.factor(tour)

ggplot(  merge(u_datFULL_ALLtimeDims ,tdf%>%dplyr::select( tour ,classification),by=c("classification")) , aes(x=as.factor(tour),y=1))+
  geom_jitter(aes(col=as.factor(tour)),size=0.1)+theme_classic()+coord_polar()+
  theme(axis.title = element_blank(),
        panel.border = element_blank(),
        legend.key = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid  = element_blank())+ylim(c(0,2))+theme(legend.position="none",aspect.ratio=1)

geom_point(data=tdf, size=3, aes(x=V1  ,      y=  V2,col=as.factor(tour)))+geom_segment(data=tdf,aes(x=V1  ,xend= V1e,     y=  V2,yend= V2e ),col="black")+
  scale_color_discrete(name="cell cycle stage")+
  labs(y="Cell cycle Umap dimension 2", x="Cell cycle Umap dimension 1")#,col=as.factor(tour)





ggplot(  merge(u_datFULL_ALLtimeDims ,tdf%>%dplyr::select( tour ,classification),by=c("classification")) , aes(V1,V2))+
  geom_point(aes(col=as.factor(classification2)),size=2)+theme_classic()

ggplot(  merge(u_datFULL_ALLtimeDims ,tdf%>%dplyr::select( tour ,classification),by=c("classification")) , aes(V1,V2))+
  geom_point(aes(col=as.factor(Response)),size=2)+theme_classic()+facet_wrap(Day~ARM)

ggplot(  merge(u_datFULL_ALLtimeDims ,tdf%>%dplyr::select( tour ,classification),by=c("classification")) , aes(V1,V2))+theme_classic()+
  geom_point(aes(col=as.factor(tour)),size=2)+#facet_wrap(Response~Day)+theme_classic()+facet_wrap(Day~ARM)+
  geom_point(data=tdf, size=3, aes(x=V1  ,      y=  V2,col=as.factor(tour)))+geom_segment(data=tdf,aes(x=V1  ,xend= V1e,     y=  V2,yend= V2e ),col="black")+
  scale_color_discrete(name="cell cycle stage")+#,col=as.factor(tour)
  labs(y="Cell cycle Umap dimension 2", x="Cell cycle Umap dimension 1")#,col=as.factor(tour)





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
  geom_smooth(aes(group=as.factor(Day),col=as.factor(Day),linetype=as.factor(Day) ),size=2.5,method="gam",formula= y~s(x,bs="cc",k=7),se=F)+
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

ggplot( merge(u_datFULL_ALLtimeDims2,xyz,by="Cell.ID")%>%group_by(ARM,Day,tour,classification3,Response,Patient.Study.ID,Subclone)%>%dplyr::summarise(JUNB=mean(JUNB),MAPK8=mean(MAPK8),MAPK10=mean(MAPK10),CDK6=mean(CDK6),ESR1=mean(ESR1) ,CDKN1B=median(CDKN1B),CDKN2A=mean(CDKN2A)), 
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

#MAPK8 JUNB CDKN2A CDK6
gamddA<-data.table(merge(u_datFULL_ALLtimeDims2,xyz,by="Cell.ID")%>%group_by(ARM,Day,tour,classification3,Response,Patient.Study.ID,Subclone)%>%dplyr::summarise(MAPK8=mean(MAPK8),MAPK10=mean(MAPK10),CDK6=mean(CDK6),ESR1=mean(ESR1) ,CDKN1B=median(CDKN1B),CDKN2A=mean(CDKN2A))
)[ARM=="A"]
gamddA[, x:=(tour-1)/max(tour-1)]
gamA <- gam(MAPK8~as.factor(Day) ,data=gamddA)
gamA <- gam(MAPK8~as.factor(Day)+s(x,by=as.factor(Day) ,k=5, bs="cs"),data=gamddA)
#plot(gam1)
summary(gamA)#[4]

gamddB<-data.table(merge(u_datFULL_ALLtimeDims2,xyz,by="Cell.ID")%>%group_by(ARM,Day,tour,classification3,Response,Patient.Study.ID,Subclone)%>%dplyr::summarise(MAPK8=mean(MAPK8),MAPK10=mean(MAPK10),CDK6=mean(CDK6),ESR1=mean(ESR1) ,CDKN1B=median(CDKN1B),CDKN2A=mean(CDKN2A))
)[ARM=="B"]
gamddB[, x:=(tour-1)/max(tour-1)]
gamB <- gam(MAPK8~as.factor(Day) ,data=gamddB)
gamB <- gam(MAPK8~as.factor(Day)+s(x,by=as.factor(Day) ,k=5, bs="cs"),data=gamddB)
#plot(gam1)
summary(gamB)


gamddA<-data.table(merge(u_datFULL_ALLtimeDims2,xyz,by="Cell.ID")%>%group_by(ARM,Day,tour,classification3,Response,Patient.Study.ID,Subclone)%>%dplyr::summarise(JUNB=mean(JUNB),MAPK8=mean(MAPK8),MAPK10=mean(MAPK10),CDK6=mean(CDK6),ESR1=mean(ESR1) ,CDKN1B=median(CDKN1B),CDKN2A=mean(CDKN2A))
)[ARM!="A"]
gamddA[, x:=(tour-1)/max(tour-1)]
gamA <- gam(CDK6~as.factor(Day) ,data=gamddA)
gamA <- gam(MAPK8~as.factor(Day)+s(x,by=as.factor(Day) ,k=5, bs="cs"),data=gamddA)
#plot(gam1)
summary(gamA)
gam1 <- gam(MAPK8~as.factor(Day)+s(x,by=as.factor(Day) ,k=5, bs="cs"),data=gamdd)
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
#ab_dat <-data.table(esr1_dd%>%group_by(Patient.Study.ID,ARM,Response,Day)%>%dplyr::summarise(ESR1=mean(ESR1)))[order(Day,ARM,Response)]
#write.csv( ab_dat[Day==180] , file="/Users/jason/Desktop/expected estrogen expression at end of trial.csv")
ggplot(ab_dat[Day==180],aes(x=ESR1,group=Patient.Study.ID,fill=Response))+geom_histogram()+facet_wrap(~ARM)
ggplot(ab_dat,aes(y=ESR1,x=Day,group=Patient.Study.ID,col=Response))+geom_point()+facet_wrap(~ARM)+geom_line()
newy <- data.table(expand.grid(ARM=unique(esr1_dd$ARM),
                               #Response=unique(esr1_dd$Response),
                               lnDay=seq(min(esr1_dd$lnDay),max(esr1_dd$lnDay),length=100)  ,Patient.Study.ID=NA))
#newy[,ARMResponse:=as.factor(paste(ARM,Response,sep="_"))]
newy[,FacPatient.Study.ID:=as.factor(Patient.Study.ID)]
newy[,FacARM:=as.factor(ARM)]
newy$Pred <- predict(gammod,newdata=newy,exclude="s(FacPatient.Study.ID)",newdata.guaranteed=TRUE)
tempd <- newy[lnDay==min(lnDay)]%>%mutate(Pred0=Pred)%>%dplyr::select(-c(Pred,lnDay))
newy2 <- merge(newy,tempd,by=names(tempd)[names(tempd)%in%names(newy)])
newy2$ARM<- as.character(newy2$ARM)
esr1_dd$ARM<- as.character(esr1_dd$ARM)

esr1_dd[(ESR1-mu_ESR1) <quantile(ESR1-mu_ESR1,probs = 0.9) & (ESR1-mu_ESR1) >quantile(ESR1-mu_ESR1,probs = 0.1),by=c("ARM","lnDay")]

esr1_dd[,uq:=quantile(ESR1-mu_ESR1,probs = 0.9) ,by=c("ARM","lnDay")]
esr1_dd[,lq:=quantile(ESR1-mu_ESR1,probs = 0.1) ,by=c("ARM","lnDay")]


esr1_dd[ (ESR1-mu_ESR1) > lq & (ESR1-mu_ESR1) <uq]
#erplot <- ggplot( newy2,aes(y=Pred-Pred0,x=lnDay,group=Response,col=Response ))+  geom_line(size=2)+facet_wrap(~ARM)
erplotN <- ggplot( newy2,aes(y=Pred-Pred0,x=exp(lnDay)-1,group=ARM,col=ARM ))+  geom_line(size=4) + theme_classic(base_size=20)+
  #geom_point(data=esr1_dd,aes(y=ESR1-mu_ESR1,x=Day))+
  geom_violin(data=esr1_dd[ (ESR1-mu_ESR1) > lq & (ESR1-mu_ESR1) <uq],aes(y=ESR1-mu_ESR1,x=exp(lnDay),fill=ARM,col=ARM,group=interaction(ARM,Day)),alpha=0.5,trim=TRUE) +
  labs(y="ESR1 expression", x="Day")
#ggsave(erplotN,filename = "/Users/jason/Desktop/ESR1 expression over time by arm.pdf", width = 9, height = 9)

erplotN2 <- ggplot( newy2,aes(y=Pred-Pred0,x=(lnDay)-1,group=ARM,col=ARM ))+  geom_line(size=4) + theme_classic(base_size=20)+
  #geom_point(data=esr1_dd,aes(y=ESR1-mu_ESR1,x=Day))+
  geom_violin(data=esr1_dd[ (ESR1-mu_ESR1) > lq & (ESR1-mu_ESR1) <uq],aes(y=ESR1-mu_ESR1,x=(lnDay),fill=ARM,col=ARM,group=interaction(ARM,Day)),alpha=0.5,trim=TRUE) +
  labs(y="ESR1 expression", x="Day")+ scale_x_continuous(breaks=log(1+c(0,14,180)) ,labels=c(0,14,180))
#ggsave(erplotN2,filename = "/Users/jason/Desktop/ESR1 expression over time by arm logday.pdf", width = 9, height = 9)

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
#ggsave(erplotNB,filename = "/Users/jason/Desktop/BLANK ESR1 expression over time by arm logday by arm.pdf", width = 15, height = 5)

newy2[,Treat:="Letrozole alone"]
newy2[ARM=="B",Treat:="Intermittent high \n dose ribociclib"]
newy2[ARM=="C",Treat:="Continuous low \n dose ribociclib"]
newy2$Treat <- factor(newy2$Treat, levels = c("Letrozole alone", "Intermittent high \n dose ribociclib", "Continuous low \n dose ribociclib"))
esr1_dd[,Treat:="Letrozole alone"]
esr1_dd[ARM=="B",Treat:="Intermittent high \n dose ribociclib"]
esr1_dd[ARM=="C",Treat:="Continuous low \n dose ribociclib"]
esr1_dd$Treat <- factor(esr1_dd$Treat, levels = c("Letrozole alone", "Intermittent high \n dose ribociclib", "Continuous low \n dose ribociclib"))

ggplot( newy2,aes(y=Pred-Pred0,x=(1+lnDay),group=ARM,col=ARM ))+  geom_line(size=4) + theme_classic(base_size=20)+
  #geom_point(data=esr1_dd,aes(y=ESR1-mu_ESR1,x=Day))+
  geom_violin(data=esr1_dd[ (ESR1-mu_ESR1) > lq & (ESR1-mu_ESR1) <uq],aes(y=ESR1-mu_ESR1,x=(1+lnDay),fill=ARM,col=ARM,group=interaction(ARM,Day)),alpha=0.5,trim=TRUE) +
  labs(y="ESR1 expression", x="Day")+ scale_x_continuous(breaks=log(1+c(0,14,180)) ,labels=c(0,14,180))+
  theme(aspect.ratio=1, legend.position="none")+  # geom_hline(yintercept=logit(1e-3+1/ 4.3 ))
  facet_wrap(~Treat)




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

#             ,size=2.5,method="gam",formula= y~s(x,k=3),se=T,alpha=0.5)+
# 
# theme(
#   legend.position="none",
#   strip.background = element_blank(),
#   strip.text.x = element_blank()
# )+
# # geom_hline(yintercept=logit(1e-3+1/ 4.3 ))
# theme(axis.title=element_blank(),
#       axis.text=element_blank(),
#       axis.ticks=element_blank())




ggplot( merge(u_datFULL_ALLtimeDims2,xyz,by="Cell.ID")%>%group_by(ARM,Day,Patient.Study.ID,Response,tour,classification3,Subclone)%>%dplyr::summarise(MAPK8=mean(MAPK8),MAPK10=mean(MAPK10),CDK6=mean(CDK6),ESR1=mean(ESR1) ), 
        aes(y=ESR1,x=log(1+Day)))+  
  facet_wrap(ARM~.)+
  theme_classic()+
  geom_point()+
  geom_smooth(aes(group=as.factor(Response),col=as.factor(Response) ),size=2.5,method="gam",formula= y~s(x,k=3),se=T,alpha=0.5)+
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






zzSUB[, classification3:= "G1"]
zzSUB[tour<(0.24*30),classification3:= "G0"]
zzSUB[tour>(0.7*30),classification3:= "S/G2"]
u_datFULL_ALLtimeDims
require(boot)
zz<- data.table(merge(merge(u_datFULL_ALLtimeDims,xyz,by="Cell.ID") ,tdf%>%dplyr::select( tour ,classification),by=c("classification"))%>%
                  group_by(Response,Subclone,Patient.Study.ID,Day,ARM)%>%
                  mutate(n_sampled=length(Cell.ID))%>%
                  group_by(classification,classification2,tour,Response,Patient.Study.ID,Subclone,Day,ARM)%>%
                  dplyr::summarise(
                    prop_tour=length(tour)/max(n_sampled),
                    MAPK10=median(MAPK10),ESR1=median(ESR1),RB1=median(RB1),CDK6=median(CDK6),CDK7=median(CDK7),E2F3=median(E2F3),CDKN2A=median(CDKN2A),CDKN1B=median(CDKN1B),CCND1=median(CCND1),CCND3=median(CCND3)
                  ))

# zz<- data.table(merge(u_datFULL_ALLtimeDims ,tdf%>%dplyr::select( tour ,classification),by=c("classification"))%>%
#                   group_by(Response,Subclone,Patient.Study.ID,Day,ARM)%>%
#                   mutate(n_sampled=length(Cell.ID))%>%
#                   group_by(classification,classification2,tour,Response,Patient.Study.ID,Subclone,Day,ARM)%>%
#                   dplyr::summarise(
#                     prop_tour=length(tour)/max(n_sampled),
#                     RB1=median(RB1),CDK6=median(CDK6),CDK7=median(CDK7),E2F3=median(E2F3),CDKN2A=median(CDKN2A),CDKN1B=median(CDKN1B),CCND1=median(CCND1),CCND3=median(CCND3)
#                     ))
zzSUB<-zz[prop_tour<0.5]# [!(classification2==2&tour>0.4*30)&!(classification2==1&tour>0.57*30)&!(classification2==1&tour<0.24*30)] 
zzSUB<-zz[prop_tour<0.5][Day==180]# [!(classification2==2&tour>0.4*30)&!(classification2==1&tour>0.57*30)&!(classification2==1&tour<0.24*30)] 
zzSUB[, classification3:= "G1"]
zzSUB[tour<(0.24*30),classification3:= "G0"]
zzSUB[tour>(0.7*30),classification3:= "S/G2"]
#zzSUB[tour>(0.57*30),classification3:= "S/G2"]
zzSUB$classification3<- as.factor(zzSUB$classification3)
zzSUB$classification2<- as.factor(zzSUB$classification2)
zzSUB[,ARMfact:=as.factor(ARM)]
zzSUB[,ARMRespfact:=as.factor(paste0(ARM,Response))]

#zzSUB<-zzSUB[ARMRespfact=="BResponder"]

ggplot( zzSUB , aes(x=tour /max(tour) ,      y=  logit(1e-3+prop_tour), col=as.factor(classification2)))+
  coord_polar()+
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  geom_point(x=0,y=logit(1e-3),pch=3,size=3,col="black")+
  geom_jitter(size=2,width = 0.005, height =0)  +
  facet_wrap(ARM~Response,ncol=2)

ggplot( zz[prop_tour<0.5] , aes(x=tour /max(tour) ,      y=  log(1e-3+ESR1-min(ESR1)), col=as.factor(classification2)))+
  coord_polar()+
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  geom_point(x=0,y=logit(1e-3),pch=3,size=3,col="black")+
  geom_jitter(size=2,width = 0.005, height =0) 


ggplot( zz[tour>(0.7*30)] , aes(x=CDK6 ,      y=  ESR1, col=as.factor(classification2)))+geom_point()+facet_wrap(~ARM)
coord_polar()+
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  geom_point(x=0,y=logit(1e-3),pch=3,size=3,col="black")+
  geom_jitter(size=2,width = 0.005, height =0)+facet_wrap(~ARM)
#f  



prop_g2 <- data.table(zzSUB%>%group_by(Patient.Study.ID, Subclone, Day,ARM,Response,ARMRespfact)%>%summarise(ESR1mu=mean(ESR1), p=sum(classification3=="S/G2")/length(classification3)   ))
ggplot( prop_g2,
        aes(x=ESR1mu ,      y=  (p), col=as.factor(ARM) ))+
  geom_point()+facet_wrap(Day~ARM)



zzingle_pred<- data.table(zzSUB%>%group_by(tour,classification3,ARMfact,Response,ARMRespfact)%>%summarise(prop_tour=median(prop_tour)))
zzingle_pred$pred<- predict( gam(formula= (prop_tour)~s(tour,by=ARMRespfact,bs="cc",k=18),data=zzSUB),newdata=zzingle_pred)
zzingle_pred[,pred:= pred/sum(pred),by="ARMRespfact"]
zzingle_pred[pred<=0,pred:= 0]
ggplot( zzingle_pred,
        aes(x=tour /max(tour) ,      y=  (pred), fill=as.factor(classification3) ))+
  geom_bar( stat="identity")+
  coord_polar()+
  theme_minimal(base_size=18)+
  labs(x="Cell cycle phase",y="Proportion of cells")+theme(aspect.ratio=1)+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())+
  scale_fill_brewer(name="Stage",palette="Dark2")+#,palette="Dark2",labels=c("G1","G0","G2/S"))+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+facet_wrap(ARMfact~Response,ncol=2)




zzSUB<-zz[prop_tour<0.5][Day==0]# [!(classification2==2&tour>0.4*30)&!(classification2==1&tour>0.57*30)&!(classification2==1&tour<0.24*30)] 
zzSUB[, classification3:= "G1"]
zzSUB[tour<(0.24*30),classification3:= "G0"]
zzSUB[tour>(0.7*30),classification3:= "S/G2"]
#zzSUB[tour>(0.57*30),classification3:= "S/G2"]
zzSUB$classification3<- as.factor(zzSUB$classification3)
zzSUB$classification2<- as.factor(zzSUB$classification2)
zzSUB[,ARMfact:=as.factor(ARM)]
zzingle_pred<- data.table(zzSUB%>%group_by(tour,classification3,ARMfact)%>%summarise(prop_tour=median(prop_tour)))
zzingle_pred$pred<- predict( gam(formula= prop_tour~s(tour,bs="cc",k=12),data=zzSUB),newdata=zzingle_pred)
zzingle_pred$pred <- zzingle_pred$pred/sum(zzingle_pred$pred)
zzingle_pred[pred<=0,pred:= 0]
ggplot( zzingle_pred,
        aes(x=tour /max(tour) ,      y=  (pred), fill=as.factor(classification3) ))+
  geom_bar( stat="identity")+
  coord_polar()+
  theme_minimal(base_size=18)+
  labs(x="Cell cycle phase",y="Proportion of cells")+theme(aspect.ratio=1)+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())+
  scale_fill_brewer(name="Stage",palette="Dark2")+#,palette="Dark2",labels=c("G1","G0","G2/S"))+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )#+facet_wrap(~ARMfact,ncol=1)





ggplot( zzSUB , aes(x=tour /max(tour) ,      y=  logit(1e-3+prop_tour)))+
  coord_polar()+
  geom_jitter(aes( col=as.factor(classification3)),size=2,width = 0.004, height =0)+  
  #facet_wrap(ARM~Response,ncol=2)+
  geom_smooth(size=1.5,method="gam",formula= y~s(x,bs="cc",k=12),se=T,col="black")+
  theme_minimal(base_size=18)+
  labs(x="Cell cycle phase",y="Proportion of cells")+theme(aspect.ratio=1)+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())+
  scale_color_brewer(name="Stage",palette="Dark2")+#,palette="Dark2",labels=c("G1","G0","G2/S"))+
  #scale_color_discrete(name="Stage")+scale_fill_discrete(name="Stage")+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )#+ # geom_hline(yintercept=logit(1e-3+1/ 2.9 ))+# geom_point(x=0,y=logit(1e-3),pch=3,size=3,col="black")

osc1<-ggplot( zzSUB[quantile(RB1, 0.95) > RB1 &quantile(RB1, 0.05) <RB1] , aes(x=(tour-1) /max(tour-1) ,      y=   (1e-3+RB1) ))+
  coord_polar()+
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  geom_point(x=0,y=logit(1e-3),pch=3,size=3,col="black")+
  geom_jitter(aes( col=as.factor(classification3)),size=2.5,width = 0.005, height =0)+  
  #facet_wrap(ARM~Response,ncol=2)+
  geom_smooth(size=1.5,method="gam",formula= y~s(x,bs="cc",k=5),se=T,col="black")+
  theme_minimal(base_size=18)+
  labs(x="Cell cycle phase",y="Proportion of cells")+theme(aspect.ratio=1)+theme(axis.title=element_blank(),
                                                                                 axis.text=element_blank(),
                                                                                 axis.ticks=element_blank())+
  scale_color_brewer(name="Stage",palette="Dark2")+
  #scale_color_discrete(name="Stage")+scale_fill_discrete(name="Stage")+
  theme(
    legend.position="none",
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+scale_x_continuous(breaks=seq(0.1,1,by=0.1))+xlim(0,1)#+
# geom_hline(yintercept=logit(1e-3+1/ 4.3 ))



osc2<-ggplot( zzSUB[quantile(CDK6, 0.95) > CDK6 &quantile(CDK6, 0.05) <CDK6] , aes(x=(tour-1) /max(tour-1) ,      y=   (1e-3+CDK6) ))+
  coord_polar()+
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  geom_point(x=0,y=logit(1e-3),pch=3,size=3,col="black")+
  geom_jitter(aes(col=as.factor(classification3)),size=2.5,width = 0.005, height =0)+  
  #facet_wrap(ARM~Response,ncol=2)+
  geom_smooth(size=1.5,method="gam",formula= y~s(x,bs="cc",k=5),se=T,col="black")+
  theme_minimal(base_size=18)+
  labs(x="Cell cycle phase",y="Proportion of cells")+theme(aspect.ratio=1)+theme(axis.title=element_blank(),
                                                                                 axis.text=element_blank(),
                                                                                 axis.ticks=element_blank())+
  scale_color_brewer(name="Stage",palette="Dark2",labels=c("c","a","b"))+
  #scale_color_discrete(name="Stage")+scale_fill_discrete(name="Stage")+
  theme(
    legend.position="none",
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+
  geom_hline(yintercept=logit(1e-3+1/ 4 ))+scale_x_continuous(breaks=seq(0.1,1,by=0.1))+xlim(0,1)


osc3<-ggplot( zzSUB[quantile(CCND3, 0.95) > CCND3 &quantile(CCND3, 0.05) <CCND3] , aes(x=(tour-1) /max(tour-1) ,      y=   (1e-3+CCND3)))+
  coord_polar()+
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  #geom_point(x=0,y=logit(1e-3),pch=3,size=3,col="black")+
  geom_jitter(aes(col=as.factor(classification3) ),size=2.5,width = 0.005, height =0)+  
  #facet_wrap(ARM~Response,ncol=2)+
  geom_smooth(size=1.5,method="gam",formula= y~s(x,bs="cc",k=5),se=T,col="black")+
  theme_minimal(base_size=18)+
  labs(x="Cell cycle phase",y="Proportion of cells")+theme(aspect.ratio=1)+theme(axis.title=element_blank(),
                                                                                 axis.text=element_blank(),
                                                                                 axis.ticks=element_blank())+
  scale_color_brewer(name="Stage",palette="Dark2")+
  #scale_color_discrete(name="Stage")+scale_fill_discrete(name="Stage")+
  theme(
    legend.position="none",
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+scale_x_continuous(breaks=seq(0.1,1,by=0.1))+xlim(0,1)#+
#geom_hline(yintercept=logit(1e-3+1/ 5 ))



osc4<-ggplot( zzSUB[quantile(E2F3, 0.95) > E2F3 &quantile(E2F3, 0.05) <E2F3] , aes(x=(tour-1) /max(tour-1) ,      y=   (1e-3+E2F3) ))+
  coord_polar()+
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  geom_point(x=0,y=logit(1e-3),pch=3,size=3,col="black")+
  geom_jitter(aes( col=as.factor(classification3)), size=2.5,width = 0.005, height =0)+  
  #facet_wrap(ARM~Response,ncol=2)+
  geom_smooth(size=1.5,method="gam",formula= y~s(x,bs="cc",k=8),se=T,col="black")+
  theme_minimal(base_size=18)+
  labs(x="Cell cycle phase",y="Proportion of cells")+theme(aspect.ratio=1)+theme(axis.title=element_blank(),
                                                                                 axis.text=element_blank(),
                                                                                 axis.ticks=element_blank())+
  scale_color_brewer(name="Stage",palette="Dark2",labels=c("c","a","b"))+
  #scale_color_discrete(name="Stage")+scale_fill_discrete(name="Stage")+
  theme(
    legend.position="none",
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+
  geom_hline(yintercept=logit(1e-3+1/ 5 ))+scale_x_continuous(breaks=seq(0.1,1,by=0.1))+xlim(0,1)

ggplot( zzSUB[quantile(RB1, 0.95) > RB1 &quantile(RB1, 0.05) <RB1] , aes(x=(tour-1) /max(tour-1) ,      y=   (1e-3+RB1) ))+
  coord_polar()+
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  geom_point(x=0,y=logit(1e-3),pch=3,size=3,col="black")+
  geom_jitter(aes( col=as.factor(tour-1)),size=2.5,width = 0.005, height =0)+  
  #facet_wrap(ARM~Response,ncol=2)+
  geom_smooth(size=1.5,method="gam",formula= y~s(x,bs="cc",k=5),se=T,col="black")+
  theme_minimal(base_size=18)+
  labs(x="Cell cycle phase",y="Proportion of cells")+theme(aspect.ratio=1)+theme(axis.title=element_blank(),
                                                                                 axis.text=element_blank(),
                                                                                 axis.ticks=element_blank())+
  # scale_color_brewer(name="Stage",palette="Dark2")+
  #scale_color_discrete(name="Stage")+scale_fill_discrete(name="Stage")+
  theme(
    legend.position="none",
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+scale_x_continuous(breaks=seq(0.1,1,by=0.1))+xlim(0,1)#+
# geom_hline(yintercept=logit(1e-3+1/ 4.3 ))




ggplot( )+
  #coord_polar()+
  geom_smooth(data=zzSUB[quantile(RB1, 0.95) > RB1 &quantile(RB1, 0.05) <RB1] , aes(x=(tour-1) /max(tour-1) ,      y=   (1e-3+RB1) ),size=1.5,method="gam",formula= y~s(x,bs="cc",k=5),se=F,col="blue")+
  geom_smooth(data=zzSUB[quantile(CDK6, 0.95) > CDK6 &quantile(CDK6, 0.05) <CDK6] , aes(x=(tour-1) /max(tour-1) ,      y=   (1e-3+CDK6) ),size=1.5,method="gam",formula= y~s(x,bs="cc",k=5),se=F,col="orange")+
  geom_smooth(data=zzSUB[quantile(E2F3, 0.95) > E2F3 &quantile(E2F3, 0.05) <E2F3] , aes(x=(tour-1) /max(tour-1) ,      y=   (1e-3+E2F3) ),size=1.5,method="gam",formula= y~s(x,bs="cc",k=8),se=F,col="green")+
  geom_smooth(data=zzSUB[quantile(CDKN2A, 0.95) > CDKN2A &quantile(CDKN2A, 0.05) <CDKN2A] , aes(x=(tour-1) /max(tour-1) ,      y=   (1e-3+CDKN2A-mean(1e-3+CDKN2A))/sd(1e-3+CDKN2A) ),size=1.5,method="gam",formula= y~s(x,bs="cc",k=12),se=F,col="red")+
  theme_minimal(base_size=18)


ggplot( )+
  coord_polar()+
  geom_smooth(data=zzSUB[quantile(RB1, 0.95) > RB1 &quantile(RB1, 0.05) <RB1] , aes(x=(tour-1) /max(tour-1) ,      y=   (1e-3+RB1-mean(1e-3+RB1))/sd(1e-3+RB1) )   ,size=2.5,method="gam",formula= y~s(x,bs="cc",k=5),se=F,col="blue")+
  geom_smooth(data=zzSUB[quantile(CDK6, 0.95) > CDK6 &quantile(CDK6, 0.05) <CDK6] , aes(x=(tour-1) /max(tour-1) ,      y=  (1e-3+CDK6-mean(1e-3+CDK6))/sd(1e-3+CDK6)  ),size=2.5,method="gam",formula= y~s(x,bs="cc",k=5),se=F,col="orange")+
  geom_smooth(data=zzSUB[quantile(E2F3, 0.95) > E2F3 &quantile(E2F3, 0.05) <E2F3] , aes(x=(tour-1) /max(tour-1) ,      y=   (1e-3+E2F3-mean(1e-3+E2F3))/sd(1e-3+E2F3)   ),size=2.5,method="gam",formula= y~s(x,bs="cc",k=8),se=F,col="green")+
  geom_smooth(data=zzSUB[quantile(CDKN2A, 0.95) > CDKN2A &quantile(CDKN2A, 0.05) <CDKN2A] , aes(x=(tour-1) /max(tour-1) ,      y=   3*(1e-3+CDKN2A-mean(1e-3+CDKN2A))/sd(1e-3+CDKN2A) ),size=2.5,method="gam",formula= y~s(x,bs="cc",k=10),se=F,col="red")+
  theme_minimal(base_size=18)+
  theme(aspect.ratio=1)+theme(axis.title=element_blank(),
                              axis.text=element_blank(),
                              axis.ticks=element_blank())+
  theme(
    legend.position="none",
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+scale_x_continuous(breaks=seq(0.1,1,by=0.1))+xlim(0,1)


ggplot( zzSUB[quantile(CDK6, 0.95) > CDK6 &quantile(CDK6, 0.05) <CDK6] , aes(x=(tour-1) /max(tour-1) ,      y=   (1e-3+CDK6) ))+
  coord_polar()+
  geom_smooth(size=1.5,method="gam",formula= y~s(x,bs="cc",k=5),se=F,col="orange")+
  theme_minimal(base_size=18)

ggplot( zzSUB[quantile(E2F3, 0.95) > CDK6 &quantile(E2F3, 0.05) <E2F3] , aes(x=(tour-1) /max(tour-1) ,      y=   (1e-3+E2F3) ))+
  coord_polar()+
  geom_smooth(size=1.5,method="gam",formula= y~s(x,bs="cc",k=5),se=F,col="green")+
  theme_minimal(base_size=18)

ggplot( zzSUB[quantile(CDKN2A, 0.95) > CDK6 &quantile(CDKN2A, 0.05) <CDKN2A] , aes(x=(tour-1) /max(tour-1) ,      y=   (1e-3+CDKN2A) ))+
  coord_polar()+
  geom_smooth(size=1.5,method="gam",formula= y~s(x,bs="cc",k=15),se=F,col="red")+
  theme_minimal(base_size=18)


ss<-6
ggsave(osc1,filename = "/Users/jason/Desktop/RB1 cell cycle expression.png", width = ss, height = ss)
ggsave(osc2,filename = "/Users/jason/Desktop/CDK6 cell cycle expression.png", width = ss, height = ss)
ggsave(osc3,filename = "/Users/jason/Desktop/CCND3 cell cycle expression.png", width = ss, height = ss)
ggsave(osc4,filename = "/Users/jason/Desktop/E2F3 cell cycle expression.png", width = ss, height = ss)



temp_dd<-merge(u_datFULL_ALLtimeDims ,tdf%>%dplyr::select( tour ,classification),by=c("classification"))
temp_dd[, classification3:= "G1"]
temp_dd[tour<(0.24*30),classification3:= "G0"]
temp_dd[tour>(0.7*30),classification3:= "S/G2"]

zznew<- data.table(temp_dd%>%
                     group_by(Response,Subclone,Patient.Study.ID,Day,ARM)%>%
                     mutate(n_sampled=length(Cell.ID))%>%
                     group_by(classification3,Response,Patient.Study.ID,Subclone,Day,ARM)%>%
                     dplyr::summarise(
                       prop_stage=length(tour)/max(n_sampled),
                       RB1=median(RB1),CDK6=median(CDK6),CDK7=median(CDK7),E2F3=median(E2F3),CDKN2A=median(CDKN2A),CDKN1B=median(CDKN1B),CCND1=median(CCND1),CCND3=median(CCND3)
                     ))
ggplot( zznew      , aes(x=log(1+Day) ,      y=  prop_stage,group=as.factor(Response),col=as.factor(Response)))+
  #coord_polar()+
  #geom_hline(yintercept=logit(1e-3+1/ 5 ))+
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  #geom_point(x=0,y=logit(1e-3),pch=3,size=3,col="black")+
  facet_wrap(ARM~classification3)+
  geom_jitter(width = 0.005, height =0)+geom_smooth()+theme_classic()

qq1<-ggplot( zznew[ARM!="A"][prop_stage!=1]      , aes(x=log(1+Day) ,      y=  prop_stage,group=(Response),fill=as.factor(Response),col=as.factor(Response)))+
  facet_wrap(~classification3)+
  geom_point(size=1.7)+geom_smooth(method="gam",method.args=list(family="quasibinomial"),
                                   
                                   formula=y~s(x,k=3),alpha=0.5,size=1.2)+theme_classic()+labs(y="Proportion of cells in stage",x="Day")+
  scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))+
  scale_fill_discrete(name="Response")+
  scale_color_discrete(name="Response")+
  theme_classic(base_size=18)+theme(aspect.ratio=1)
ggsave(qq1,filename = "/Users/jason/Desktop/prop cells at diff cell cycle stages.png", width = 15, height = 5)

qq2<-ggplot( zznew[classification3=="S/G2"][prop_stage!=1]      , aes(x=log(1+Day) ,      y=  prop_stage,group=(Response),fill=as.factor(Response),col=as.factor(Response)))+
  facet_wrap(~ARM)+
  geom_point(size=2)+geom_smooth(method="gam",method.args=list(family="quasibinomial",gamma=0.4),
                                 
                                 formula=y~s(x,k=3),alpha=0.5,size=1.2)+theme_classic()+labs(y="Proportion of cells in stage",x="Day")+
  scale_x_continuous(breaks=log(1+c(0,14,180)),labels=c(0,14,180))+
  scale_fill_discrete(name="Response")+
  scale_color_discrete(name="Response")+
  theme_classic(base_size=18)+theme(aspect.ratio=1)+
  theme(aspect.ratio=1,
        strip.background = element_blank(),
        strip.text.x = element_blank()
  ) +
  theme(axis.title=element_blank(), axis.text=element_blank(),legend.title=element_blank(),legend.text=element_blank()) +theme(legend.position = "none")

ggsave(qq2,filename = "/Users/jason/Dropbox/FELINE Project/Manuscript/Final Main Figures/Blank figures/G2 phase over time across arms.png", width = 15, height = 5)

glm(prop_stage~ log(1+Day)  , family="quasibinomial" ,data= zznew[classification3=="S/G2"][prop_stage!=1][ARM!="A"])%>%summary()
glm(prop_stage~ log(1+Day)*Response  , family="quasibinomial" ,data= zznew[classification3=="S/G2"][prop_stage!=1][ARM=="B"])%>%summary()
glm(prop_stage~ log(1+Day)*Response  , family="quasibinomial" ,data= zznew[classification3=="S/G2"][prop_stage!=1][ARM=="A"])%>%summary()


ggplot( zznew      , aes(x=classification3 ,      y=  logit(1e-3+prop_stage),group=as.factor(Day),col=as.factor(Day)))+
  #coord_polar()+
  #geom_hline(yintercept=logit(1e-3+1/ 5 ))+
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  #geom_point(x=0,y=logit(1e-3),pch=3,size=3,col="black")+
  facet_wrap(ARM~Response,ncol=2)+
  geom_jitter(data=zznew,size=2,aes(x=classification3 ,      y=  logit(1e-3+prop_stage),group=as.factor(Day),col=as.factor(Day)),
              width = 0.005, height =0)+#facet_wrap(Response~Day)+
  theme_void(base_size=18)+
  geom_smooth(aes(group=as.factor(Day),col=as.factor(Day),fill=as.factor(Day)),size=1.5,method="gam",formula= y~s(x,bs="cc",k=9),se=T)+
  labs(x="Cell cycle phase",y="Proportion of cells")+theme(aspect.ratio=1)+theme(axis.title=element_blank(),
                                                                                 axis.text=element_blank(),
                                                                                 axis.ticks=element_blank())+
  scale_color_discrete(name="Day")+
  scale_fill_discrete(name="Day")+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )




ggplot( zzSUB[quantile(CDKN2A, 0.95) > CDKN2A &quantile(CDKN2A, 0.05) <CDKN2A] , aes(x=tour /max(tour) ,      y=   (1e-3+CDKN2A) ))+
  coord_polar()+
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  geom_point(x=0,y=logit(1e-3),pch=3,size=3,col="black")+
  geom_jitter(aes( col=as.factor(classification3)),size=2,width = 0.005, height =0)+  
  #facet_wrap(ARM~Response,ncol=2)+
  geom_smooth(size=1.5,method="gam",formula= y~s(x,bs="cc",k=5),se=T,col="black")+
  theme_minimal(base_size=18)+
  labs(x="Cell cycle phase",y="Proportion of cells")+theme(aspect.ratio=1)+theme(axis.title=element_blank(),
                                                                                 axis.text=element_blank(),
                                                                                 axis.ticks=element_blank())+
  scale_color_brewer(name="Stage",palette="Dark2",labels=c("c","a","b"))+
  #scale_color_discrete(name="Stage")+scale_fill_discrete(name="Stage")+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+
  geom_hline(yintercept=logit(1e-3+1/ 5 ))

ggplot( zzSUB[quantile(CDKN1B, 0.95) > CDKN1B &quantile(CDKN1B, 0.05) <CDKN1B] , aes(x=tour /max(tour) ,      y=   (1e-3+CDKN1B)))+
  coord_polar()+
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  geom_point(x=0,y=logit(1e-3),pch=3,size=3,col="black")+
  geom_jitter(aes( col=as.factor(classification3)),size=2,width = 0.005, height =0)+  
  #facet_wrap(ARM~Response,ncol=2)+
  geom_smooth(size=1.5,method="gam",formula= y~s(x,bs="cc",k=5),se=T,col="black")+
  theme_minimal(base_size=18)+
  labs(x="Cell cycle phase",y="Proportion of cells")+theme(aspect.ratio=1)+theme(axis.title=element_blank(),
                                                                                 axis.text=element_blank(),
                                                                                 axis.ticks=element_blank())+
  scale_color_brewer(name="Stage",palette="Dark2",labels=c("c","a","b"))+
  #scale_color_discrete(name="Stage")+scale_fill_discrete(name="Stage")+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+
  geom_hline(yintercept=logit(1e-3+1/ 5 ))

ggplot( zzSUB[quantile(CDK7, 0.95) > CDK7 &quantile(CDK7, 0.05) <CDK7] , aes(x=tour /max(tour) ,      y=   (1e-3+CDK7)))+
  coord_polar()+
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  geom_point(x=0,y=logit(1e-3),pch=3,size=3,col="black")+
  geom_jitter(aes( col=as.factor(classification3)),size=2,width = 0.005, height =0)+  
  #facet_wrap(ARM~Response,ncol=2)+
  geom_smooth(size=1.5,method="gam",formula= y~s(x,bs="cc",k=5),se=T,col="black")+
  theme_minimal(base_size=18)+
  labs(x="Cell cycle phase",y="Proportion of cells")+theme(aspect.ratio=1)+theme(axis.title=element_blank(),
                                                                                 axis.text=element_blank(),
                                                                                 axis.ticks=element_blank())+
  scale_color_brewer(name="Stage",palette="Dark2",labels=c("c","a","b"))+
  #scale_color_discrete(name="Stage")+scale_fill_discrete(name="Stage")+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )








zza<-data.table(zz%>%
                  group_by(classification,tour,Response,Day,ARM)%>%
                  dplyr::summarise(prop_tour=median(prop_tour)))

ggplot( zz[prop_tour<0.5]       , aes(x=tour /max(tour) ,      y=  logit(1e-3+prop_tour),group=as.factor(Day),col=as.factor(Day)))+
  coord_polar()+
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  geom_point(x=0,y=logit(1e-3),pch=3,size=3,col="black")+
  facet_wrap(ARM~Response,ncol=2)+
  geom_jitter(data=zza,size=2,aes(x=tour /max(tour) ,      y=  logit(1e-3+prop_tour),group=as.factor(Day),col=as.factor(Day)),
              width = 0.005, height =0)+#facet_wrap(Response~Day)+
  theme_void(base_size=18)+
  geom_smooth(aes(group=as.factor(Day),col=as.factor(Day),fill=as.factor(Day)),size=1.5,method="gam",formula= y~s(x,bs="cc",k=7),se=T)+
  labs(x="Cell cycle phase",y="Proportion of cells")+theme(aspect.ratio=1)+theme(axis.title=element_blank(),
                                                                                 axis.text=element_blank(),
                                                                                 axis.ticks=element_blank())+
  scale_color_discrete(name="Day")+
  scale_fill_discrete(name="Day")+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )


ggplot( zz[prop_tour<0.5]      , aes(x=tour /max(tour) ,      y=  logit(1e-3+prop_tour),group=as.factor(Day),col=as.factor(Day)))+
  coord_polar()+
  geom_hline(yintercept=logit(1e-3+1/ 5 ))+
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  geom_point(x=0,y=logit(1e-3),pch=3,size=3,col="black")+
  facet_wrap(ARM~Response,ncol=2)+
  geom_jitter(data=zza,size=2,aes(x=tour /max(tour) ,      y=  logit(1e-3+prop_tour),group=as.factor(Day),col=as.factor(Day)),
              width = 0.005, height =0)+#facet_wrap(Response~Day)+
  theme_void(base_size=18)+
  geom_smooth(aes(group=as.factor(Day),col=as.factor(Day),fill=as.factor(Day)),size=1.5,method="gam",formula= y~s(x,bs="cc",k=9),se=T)+
  labs(x="Cell cycle phase",y="Proportion of cells")+theme(aspect.ratio=1)+theme(axis.title=element_blank(),
                                                                                 axis.text=element_blank(),
                                                                                 axis.ticks=element_blank())+
  scale_color_discrete(name="Day")+
  scale_fill_discrete(name="Day")+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )


zzb<-data.table(zz%>%
                  group_by(classification2,Response,Day,ARM)%>%
                  dplyr::summarise(prop_tour=median(prop_tour)))


ggplot( zzb[prop_tour<0.5][classification2==3]       , aes(x=Day ,      y=  logit(1e-3+prop_tour),col=as.factor(Day)))+
  #coord_polar()+
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  geom_point(x=0,y=logit(1e-3),pch=3,size=3,col="black")+
  facet_wrap(ARM~Response,ncol=2)+
  geom_jitter(size=2,width = 0.00, height =0)+#facet_wrap(Response~Day)+
  theme_void(base_size=18)+
  geom_smooth(aes(group=as.factor(Day),col=as.factor(Day),fill=as.factor(Day)),size=1.5,method="gam",formula= y~s(x,bs="cc",k=8),se=T)+
  labs(x="Cell cycle phase",y="Proportion of cells")+theme(aspect.ratio=1)+theme(axis.title=element_blank(),
                                                                                 axis.text=element_blank(),
                                                                                 axis.ticks=element_blank())+
  scale_color_discrete(name="Day")+
  scale_fill_discrete(name="Day")+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )


load(file="/Users/jason/Dropbox/FELINE Project/Data_analysis/scRNA/15_subclone_interaction/MAPK major axis proj.RData")
setnames(MAPK_Major_axis, old=c("V1","V2"), new=c("mapk_axisV1","mapk_axisV2"))
#MAPK_Major_axis,which_plot,corV3,subset_u_datFULL_ALLtimeDims,
mapk_ttt<-data.table(ttt3[Variable%in%c("MAP3K5","MAPK8","MAPK10","MAP3K1","MAP3K14","MAP4K3")]%>%dplyr::select(Cell.ID,Variable,Expression)   %>%spread(Variable,Expression)%>%dplyr::select(Cell.ID,MAP3K5,MAPK8,MAPK10,MAP3K1,MAP3K14,MAP4K3))
mapk_ttt<-data.table(ttt3[Major_catagory%in%c("Signal transduction","Cell surface receptor")][VariableType=="gene"]%>%dplyr::select(Cell.ID,Variable,Expression)   %>%spread(Variable,Expression)  )#%>%dplyr::select(Cell.ID,MAP3K5,MAPK8,MAPK10,MAP3K1,MAP3K14,MAP4K3))

mrg_umap_mapk_tsp  <-merge(u_datFULL_ALLtimeDims2,mapk_ttt,by="Cell.ID") 
mrg_umap_mapk_tsp2 <- merge(  mrg_umap_mapk_tsp, MAPK_Major_axis, by=c("Cell.ID", "Patient.Study.ID", "ARM", "Day",      "Response", "Subclone"))

save(mrg_umap_mapk_tsp2,zz,MAPK_Major_axis,file="/Users/jason/Dropbox/FELINE Project/Data_analysis/scRNA/15_subclone_interaction/merging rtk mapk and cell cycle.RData")

ggplot(  mrg_umap_mapk_tsp2[abs(V1_disc)<2], aes(y=CDK6 ,      x=  V1_disc,col=Response))+
  #coord_polar()+
  #geom_jitter(size=2)+
  theme_classic(base_size=18)+
  geom_smooth(aes(group=interaction(Response),col=Response,fill=Response),size=1.5,method="gam",formula=y~s(x,k=6),se=T)+
  theme(aspect.ratio=1)+
  facet_wrap(ARM~Day) + theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) +labs(x="JNK activation score" , x= "CDK6")+
  geom_point(data=mrg_umap_mapk_tsp2[abs(V1_disc)<2]%>%group_by(Response, ARM,Day,V1_disc)%>%summarise(CDK6=mean(CDK6)), aes(y=CDK6 ,      x=  V1_disc,col=Response) )

ggplot( mrg_umap_mapk_tsp2[tour%in%c(29)][Day==180],aes(y=log(CDK6), x=V1_disc,col=tour))+geom_point()+facet_grid(ARM~Response)+
  geom_smooth(aes(group=interaction(tour),col=tour,fill=tour),size=1.5,method="gam",formula=y~s(x,k=3),se=F)



ggplot(  mrg_umap_mapk_tsp2, aes(x=tour ,   y=  V2_disc,col=as.factor(Day),fill=as.factor(Day)))+
  coord_polar()+
  # geom_jitter(size=2)+
  theme_minimal(base_size=18)+
  geom_smooth(aes(group=interaction(Day)),size=1.5,method="gam",formula=y~s(x,k=3,bs="cc"),se=T)+
  theme(aspect.ratio=1)+
  facet_wrap(ARM~Response,ncol=2) 




zzMAPK<- data.table(mrg_umap_mapk_tsp2%>%
                      group_by(Response,Subclone,Patient.Study.ID,Day,ARM)%>%
                      mutate(n_sampled=length(Cell.ID))%>%
                      group_by(classification,classification2,classification3,tour,Response,Patient.Study.ID,Day,ARM)%>%# removed subclone grouping
                      dplyr::summarise(
                        prop_tour=length(tour)/max(n_sampled),
                        RB1=median(RB1),CDK6=median(CDK6),CDK7=median(CDK7),E2F3=median(E2F3),CDKN2A=median(CDKN2A),CDKN1B=median(CDKN1B),CCND1=median(CCND1),CCND3=median(CCND3),
                        AKT3=median(AKT3),FOS=median(FOS),JUNB=median(JUNB),MAPK10=median(MAPK10),MAPK3=median(MAPK3),MAP3K1=median(MAP3K1),V1_disc=median(V1_disc),
                        ESR1=median(ESR1),ESR2=median(ESR2),ERBB2=median(ERBB2),ERBB3=median(ERBB3),ERBB4=median(ERBB4) ,FGFR1=median(FGFR1),FGFR2=mean(FGFR2),RORA=median(RORA),RORC=median(RORC),ROR2=median(ROR2) )
)

ggplot( zzMAPK[prop_tour<0.5] , aes(x=tour /max(tour) ,      y=  logit(prop_tour), col=as.factor(classification3)))+
  coord_polar()+
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  geom_point(x=0,y=logit(1e-3),pch=3,size=3,col="black")+
  geom_jitter(size=2,width = 0.005, height =0)  
#facet_wrap(ARM~Response,ncol=2)

zzMAPKSUB<-zzMAPK[prop_tour<0.5]
zzMAPKSUB$classification3<- as.factor(zzMAPKSUB$classification3)
zzMAPKSUB$classification2<- as.factor(zzMAPKSUB$classification2)
zzMAPKSUB[,ARMfact:=as.factor(ARM)]
zzMAPKSUB[,ARMRespfact:=as.factor(paste0(ARM,Response))]
#[quantile(V1_disc, 0.95) > V1_disc &quantile(V1_disc, 0.05) <V1_disc]
ggplot( zzMAPKSUB[Day==180] , aes(x=tour /max(tour) ,      y=   (1e-3+V1_disc)))+
  coord_polar()+
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  geom_point(x=0,y=logit(1e-3),pch=3,size=3,col="black")+
  geom_jitter(aes( col=as.factor(classification3)),size=2,width = 0.005, height =0)+  
  facet_wrap(ARM~Response,ncol=2)+
  geom_smooth(size=1.5,method="gam",formula= y~s(x,bs="cc",k=5),se=T,col="black")+
  theme_minimal(base_size=18)+
  labs(x="Cell cycle phase",y="Proportion of cells")+theme(aspect.ratio=1)+theme(axis.title=element_blank(),
                                                                                 axis.text=element_blank(),
                                                                                 axis.ticks=element_blank())+
  scale_color_brewer(name="Stage",palette="Dark2",labels=c("c","a","b"))+
  #scale_color_discrete(name="Stage")+scale_fill_discrete(name="Stage")+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+
  geom_hline(yintercept=logit(1e-3+1/ 5 ))

ggplot( zzMAPKSUB[Day==0] , aes(x=tour /max(tour) ,      y=   (1e-3+FGFR2)))+
  coord_polar()+
  geom_jitter(aes( col=as.factor(classification3)),size=2,width = 0.005, height =0)+  
  facet_wrap(ARM~Response,ncol=2)+
  geom_smooth(size=1.5,method="gam",formula= y~s(x,bs="cc",k=5),se=T,col="black")+
  theme_minimal(base_size=18)+
  labs(x="Cell cycle phase",y="Proportion of cells")+theme(aspect.ratio=1)+theme(axis.title=element_blank(),
                                                                                 axis.text=element_blank(),
                                                                                 axis.ticks=element_blank())+
  scale_color_brewer(name="Stage",palette="Dark2",labels=c("c","a","b"))+
  #scale_color_discrete(name="Stage")+scale_fill_discrete(name="Stage")+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

ggplot( zzMAPKSUB[Day==180] , aes(x=tour /max(tour) ,      y=   (1e-3+FGFR2)))+
  coord_polar()+
  geom_jitter(aes( col=as.factor(classification3)),size=2,width = 0.005, height =0)+  
  facet_wrap(ARM~Response,ncol=2)+
  geom_smooth(size=1.5,method="gam",formula= y~s(x,bs="cc",k=5),se=T,col="black")+
  theme_minimal(base_size=18)+
  labs(x="Cell cycle phase",y="Proportion of cells")+theme(aspect.ratio=1)+theme(axis.title=element_blank(),
                                                                                 axis.text=element_blank(),
                                                                                 axis.ticks=element_blank())+
  scale_color_brewer(name="Stage",palette="Dark2",labels=c("c","a","b"))+
  #scale_color_discrete(name="Stage")+scale_fill_discrete(name="Stage")+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

ggplot( zzMAPKSUB[Day==180] , aes(x=tour /max(tour) ,      y=   (1e-3+ESR1)))+
  coord_polar()+
  geom_jitter(aes( col=as.factor(classification3)),size=2,width = 0.005, height =0)+  
  facet_wrap(ARM~Response,ncol=2)+
  geom_smooth(size=1.5,method="gam",formula= y~s(x,bs="cc",k=5),se=T,col="black")+
  theme_minimal(base_size=18)+
  labs(x="Cell cycle phase",y="Proportion of cells")+theme(aspect.ratio=1)+theme(axis.title=element_blank(),
                                                                                 axis.text=element_blank(),
                                                                                 axis.ticks=element_blank())+
  scale_color_brewer(name="Stage",palette="Dark2",labels=c("c","a","b"))+
  #scale_color_discrete(name="Stage")+scale_fill_discrete(name="Stage")+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )


ggplot( zzMAPKSUB[Day==180] , aes(x=tour /max(tour) ,      y=   (1e-3+ERBB4)))+
  coord_polar()+
  geom_jitter(aes( col=as.factor(classification3)),size=2,width = 0.005, height =0)+  
  facet_wrap(ARM~Response,ncol=2)+
  geom_smooth(size=1.5,method="gam",formula= y~s(x,bs="cc",k=5),se=T,col="black")+
  theme_minimal(base_size=18)+
  labs(x="Cell cycle phase",y="Proportion of cells")+theme(aspect.ratio=1)+theme(axis.title=element_blank(),
                                                                                 axis.text=element_blank(),
                                                                                 axis.ticks=element_blank())+
  scale_color_brewer(name="Stage",palette="Dark2",labels=c("c","a","b"))+
  #scale_color_discrete(name="Stage")+scale_fill_discrete(name="Stage")+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

ggplot( zzMAPKSUB[Day==180] , aes(x=tour /max(tour) ,      y=   (1e-3+AKT3)))+
  coord_polar()+
  geom_jitter(aes( col=as.factor(classification3)),size=2,width = 0.005, height =0)+  
  facet_wrap(ARM~Response,ncol=2)+
  geom_smooth(size=1.5,method="gam",formula= y~s(x,bs="cc",k=5),se=T,col="black")+
  theme_minimal(base_size=18)+
  labs(x="Cell cycle phase",y="Proportion of cells")+theme(aspect.ratio=1)+theme(axis.title=element_blank(),
                                                                                 axis.text=element_blank(),
                                                                                 axis.ticks=element_blank())+
  scale_color_brewer(name="Stage",palette="Dark2",labels=c("c","a","b"))+
  #scale_color_discrete(name="Stage")+scale_fill_discrete(name="Stage")+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )


ggplot( zzMAPKSUB[Day==180] , aes(x=tour /max(tour) ,      y=   (1e-3+ESR2)))+
  coord_polar()+
  geom_jitter(aes( col=as.factor(classification3)),size=2,width = 0.005, height =0)+  
  facet_wrap(ARM~Response,ncol=2)+
  geom_smooth(size=1.5,method="gam",formula= y~s(x,bs="cc",k=5),se=T,col="black")+
  theme_minimal(base_size=18)+
  labs(x="Cell cycle phase",y="Proportion of cells")+theme(aspect.ratio=1)+theme(axis.title=element_blank(),
                                                                                 axis.text=element_blank(),
                                                                                 axis.ticks=element_blank())+
  scale_color_brewer(name="Stage",palette="Dark2",labels=c("c","a","b"))+
  #scale_color_discrete(name="Stage")+scale_fill_discrete(name="Stage")+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

zzMAPKSUB[,FGFR2round:=round(FGFR2)]

ggplot( zzMAPKSUB , aes(x=(V1_disc)    ))+
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  geom_density(aes( fill=as.factor(Day)),alpha = 0.5)+facet_wrap(Response~ARM)

ggplot( zzMAPKSUB[Day==180] , aes(x=(V1_disc)    ))+
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  geom_density(aes( fill=as.factor(Response)),alpha = 0.5)+facet_wrap(Day~ARM)


mrg_umap_mapk_tsp2[,V1_disc2:=round(mapk_axisV1*35)/35]


ggplot( mrg_umap_mapk_tsp2[Day==180]%>%group_by(V1_disc)%>%summarise(FGFR2=median(FGFR2),ESR1=median(ESR1),CDK6=median(CDK6)), aes(y=log(FGFR2-min(FGFR2)) ,      x=   (V1_disc),col=log(ESR1-min(ESR1)) )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  #facet_wrap(ARM~.,ncol=3)+
  geom_smooth(size=1.5, method="gam",formula= y~s(x,k=4 ),se=T,col="firebrick",fill="firebrick")+
  geom_point(size=2)


ggplot( mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9][ARM!="A"]%>%group_by(V1_disc2)%>%summarise(FGFR2=mean(FGFR2),ESR1=median(ESR1),CDK6=mean(CDK6)), aes(y=(CDK6) ,      x=   (V1_disc2),col=log(FGFR2-min(FGFR2)) )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  #facet_wrap(ARM~.,ncol=3)+
  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=4), se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+labs(y="CDK6 ",x="JNK activation score")+
  scale_color_viridis(option="B",begin = 0.1,name="FGFR2")+
  theme(aspect.ratio=1,
        strip.background = element_blank(),
        strip.text.x = element_blank()
  )


ggplot( mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9]%>%group_by(ARM,V1_disc2)%>%summarise(FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)), aes(y=(CDK6) ,      x=   (V1_disc2),col=log(ESR1-min(ESR1)) )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(~ARM,ncol=3)+
  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=4), se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+labs(y="CDK6 ",x="JNK activation score")+
  scale_color_viridis(option="B",begin = 0.1,name="ESR1")+
  theme(aspect.ratio=1,
        strip.background = element_blank(),
        strip.text.x = element_blank()
  )

ggplot( mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9]%>%group_by(ARM,V1_disc2)%>%summarise(FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)), aes(y=(CDK6) ,      x=   (V1_disc2),col=log(FGFR2-min(FGFR2)) )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(~ARM,ncol=3)+
  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=4), se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+labs(y="CDK6 ",x="JNK activation score")+
  scale_color_viridis(option="B",begin = 0.1,name="FGFR2")+
  theme(aspect.ratio=1,
        strip.background = element_blank(),
        strip.text.x = element_blank()
  )


ggplot( mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9]%>%group_by(Response,ARM,V1_disc2)%>%summarise(FGFR2=mean(FGFR2),ESR1=mean(ESR1),ERBB4=mean(ERBB4),CDK6=mean(CDK6),RB1=mean(RB1),E2F3=mean(E2F3)), aes(y=(CDK6) ,      x=   (ERBB4),col=log(ESR1-min(ESR1)) )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(Response~ARM,ncol=3)+
  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=4), se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+labs(y="CDK6 ",x="ERBB4 score")+
  scale_color_viridis(option="B",begin = 0.1,name="ESR1")+
  theme(aspect.ratio=1,
        strip.background = element_blank(),
        strip.text.x = element_blank()
  )

ggplot( mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9]%>%group_by(Response,ARM,V1_disc2)%>%summarise(FGFR2=mean(FGFR2),ESR1=mean(ESR1),ERBB4=mean(ERBB4),CDK6=mean(CDK6),RB1=mean(RB1),E2F3=mean(E2F3)), aes(y=(CDK6) ,      x=   (ESR1),col=log(ERBB4-min(ERBB4)) )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(Response~ARM,ncol=3)+
  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=4), se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+labs(y="CDK6 ",x="ERBB4 score")+
  scale_color_viridis(option="B",begin = 0.1,name="ESR1")+
  theme(aspect.ratio=1,
        strip.background = element_blank(),
        strip.text.x = element_blank()
  )


ggplot( mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9]%>%group_by(Response,ARM,V1_disc2)%>%summarise(FGFR2=mean(FGFR2),ESR1=mean(ESR1),ERBB4=mean(ERBB4),CDK6=mean(CDK6),RB1=mean(RB1),E2F3=mean(E2F3)), aes(y=(CDK6) ,      x=   (ESR1),col=Response,group=Response )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(~ARM,ncol=3)+
  #  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=4), se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+labs(y="CDK6 ",x="ESR1 score")+
  # scale_color_viridis(option="B",begin = 0.1,name="ESR1")+
  theme(aspect.ratio=1,
        strip.background = element_blank(),
        strip.text.x = element_blank()
  )

ggplot( mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9]%>%group_by(Response,ARM,V1_disc2)%>%summarise(FGFR2=mean(FGFR2),ESR1=mean(ESR1),ERBB4=mean(ERBB4),CDK6=mean(CDK6),RB1=mean(RB1),E2F3=mean(E2F3)), aes(y=(CDK6) ,      x=   (ERBB4),col=Response,group=Response )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(~ARM,ncol=3)+
  #  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=4), se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+labs(y="CDK6 ",x="ERBB4 score")+
  # scale_color_viridis(option="B",begin = 0.1,name="ESR1")+
  theme(aspect.ratio=1,
        strip.background = element_blank(),
        strip.text.x = element_blank()
  )

ggplot( mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9]%>%group_by(Response,ARM,V1_disc2)%>%summarise(FGFR2=mean(FGFR2),ESR1=mean(ESR1),ERBB4=mean(ERBB4),CDK6=mean(CDK6),RB1=mean(RB1),E2F3=mean(E2F3)), aes(y=(CDK6) ,      x=   (FGFR2),col=Response,group=Response )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(~ARM,ncol=3)+
  #  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=4), se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+labs(y="CDK6 ",x="FGFR2 score")+
  # scale_color_viridis(option="B",begin = 0.1,name="ESR1")+
  theme(aspect.ratio=1,
        strip.background = element_blank(),
        strip.text.x = element_blank()
  )


ggplot( mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9][ARM!="A"]%>%group_by(ARM,Response,V1_disc2)%>%summarise(FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)), aes(y=(FGFR2) ,      x=   (ESR1) ,col=Response,group=Response,fill=Response)) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(ARM~.,ncol=3)+
  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=3), se=T)+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)#+#labs(y="CDK6 ",x="JNK activation score")+
#scale_color_viridis(option="B",begin = 0.1,name="ESR1")


ggplot( mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9][ARM=="A"]%>%group_by(V1_disc2)%>%summarise(FGFR2=mean(FGFR2),ESR1=median(ESR1),CDK6=mean(CDK6)), aes(y=(CDK6) ,      x=   (V1_disc2),col=log(FGFR2-min(FGFR2)) )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  #facet_wrap(ARM~.,ncol=3)+
  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=4), se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+labs(y="CDK6 ",x="JNK activation score")+
  scale_color_viridis(option="B",begin = 0.1)



########
test1<- zzMAPKSUB[ARM=="B"]
test1[,fgfrtoERS1:= (FGFR2-min(FGFR2))/(ESR1-min(ESR1))]
ggplot( test1[fgfrtoERS1<1.5] , aes(x=fgfrtoERS1 ,      y=   (V1_disc)))+
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  # facet_wrap(ARM~.,ncol=3)+
  geom_smooth(size=1.5, method="gam",formula= y~s(x,k=3 ),se=T,col="firebrick",fill="firebrick")+
  geom_point( col="firebrick2",size=2)+  
  theme_classic(base_size=18)+
  labs(x="FGFR2 expression",y="JNK activation score")+theme(aspect.ratio=1)




ggplot( zzMAPKSUB[Day==180][FGFR2>-1.4] , aes(x=log(FGFR2+1.5) ,      y=   (V1_disc)))+
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  # facet_wrap(ARM~.,ncol=3)+
  geom_smooth(size=1.5, method="gam",formula= y~s(x,k=4 ),se=T,col="firebrick",fill="firebrick")+
  geom_point( col="firebrick2",size=2)+  
  theme_classic(base_size=18)+
  labs(x="FGFR2 expression",y="JNK activation score")+theme(aspect.ratio=1)


ggplot( zzMAPKSUB[Day==180][FGFR2>-1.4] , aes(x=log(FGFR2+1.5) ,      y=   (V1_disc)))+
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  # facet_wrap(ARM~.,ncol=3)+
  geom_smooth(size=1.5, method="gam",formula= y~s(x,k=4 ),se=T,col="firebrick",fill="firebrick")+
  geom_point( col="firebrick2",size=2)+  
  theme_classic(base_size=18)+
  labs(x="FGFR2 expression",y="JNK activation score")+theme(aspect.ratio=1)


ggplot( zzMAPKSUB%>%group_by(Response, Patient.Study.ID, Subclone, Day ,ARM) %>%summarise(V1_disc=mean(V1_disc),FGFR2=mean(FGFR2)), aes(x=log(FGFR2-min(FGFR2)) ,     
                                                                                                                                        y=   (V1_disc)))+
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  geom_jitter(aes( col=ARM),size=2,width = 0.05, height =0)+  
  # facet_wrap(ARM~.,ncol=3)+
  geom_smooth(size=1.5, method="gam",formula= y~s(x,k=7 ),se=T,col="black")+
  theme_classic(base_size=18)+
  labs(x="FGFR2 expression",y="JNK activation score")+theme(aspect.ratio=1)+#theme(axis.title=element_blank(),
  #     axis.text=element_blank(),
  #     axis.ticks=element_blank())+
  scale_color_brewer(name="Stage",palette="Dark2")

ggplot( mrg_umap_mapk_tsp2[Day==180] , aes(x=(FGFR2) ,      y=   (V1_disc)))+
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  geom_jitter(aes( col=ARM),size=2,width = 0.05, height =0)+  
  # facet_wrap(ARM~.,ncol=3)+
  geom_smooth(size=1.5, method="gam",formula= y~s(x,k=4 ),se=T,col="black")+
  theme_classic(base_size=18)+
  labs(x="FGFR2 expression",y="JNK activation score")+theme(aspect.ratio=1)+#theme(axis.title=element_blank(),
  #     axis.text=element_blank(),
  #     axis.ticks=element_blank())+
  scale_color_brewer(name="Stage",palette="Dark2")




ggplot( zzMAPKSUB[Day==180] , aes(x=(V1_disc) ,      y=   (CDK6),group=tour ,col=classification3,fill=classification3))+
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  geom_jitter(size=2,width = 0.005, height =0)+  
  #facet_wrap(~tour)+
  geom_smooth(size=1.5,method="gam",formula= y~s(x,k=4),se=T)+
  theme_minimal(base_size=18)+
  labs(x="Cell cycle phase",y="Proportion of cells")+theme(aspect.ratio=1)+theme(axis.title=element_blank(),
                                                                                 axis.text=element_blank(),
                                                                                 axis.ticks=element_blank())+
  scale_color_brewer(name="Stage",palette="Dark2")+scale_fill_brewer(name="Stage",palette="Dark2")





ggplot( zzMAPKSUB[Day==180] , aes(y=log(1+CDKN2A-min(CDKN2A)) ,      x=   (V1_disc)))+
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  geom_jitter(aes( col=as.factor(classification3)),size=2,width = 0.005, height =0)+  
  facet_wrap(ARM~Day,ncol=3)+
  geom_smooth(size=1.5,method="gam",formula= y~s(x,k=3),se=T,col="black")+
  theme_minimal(base_size=18)+
  labs(x="Cell cycle phase",y="Proportion of cells")+theme(aspect.ratio=1)+theme(axis.title=element_blank(),
                                                                                 axis.text=element_blank(),
                                                                                 axis.ticks=element_blank())+
  scale_color_brewer(name="Stage",palette="Dark2")


ggplot( zzMAPKSUB[Day==180] , aes(y=ERBB4 ,      x=   (V1_disc)))+
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  geom_jitter(aes( col=as.factor(classification3)),size=2,width = 0.005, height =0)+  
  facet_wrap(ARM~Response,ncol=2)+
  geom_smooth(size=1.5,method="gam",formula= y~s(x,k=3),se=T,col="black")+
  theme_minimal(base_size=18)+
  labs(x="Cell cycle phase",y="Proportion of cells")+theme(aspect.ratio=1)+theme(axis.title=element_blank(),
                                                                                 axis.text=element_blank(),
                                                                                 axis.ticks=element_blank())+
  scale_color_brewer(name="Stage",palette="Dark2")

ggplot( zzMAPKSUB[Day==180] , aes(y=FGFR2 ,      x=   (V1_disc)))+
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  geom_jitter(aes( col=as.factor(classification3)),size=2,width = 0.005, height =0)+  
  facet_wrap(ARM~Response,ncol=2)+
  geom_smooth(size=1.5,method="gam",formula= y~s(x,k=3),se=T,col="black")+
  theme_minimal(base_size=18)+
  labs(x="Cell cycle phase",y="Proportion of cells")+theme(aspect.ratio=1)+theme(axis.title=element_blank(),
                                                                                 axis.text=element_blank(),
                                                                                 axis.ticks=element_blank())+
  scale_color_brewer(name="Stage",palette="Dark2")


ggplot( zzMAPKSUB[Day==180] , aes(y=ESR1 ,      x=   (V1_disc)))+
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  geom_jitter(aes( col=as.factor(classification3)),size=2,width = 0.005, height =0)+  
  facet_wrap(ARM~Response,ncol=2)+
  geom_smooth(size=1.5,method="gam",formula= y~s(x,k=5),se=T,col="black")+
  theme_minimal(base_size=18)+
  labs(x="Cell cycle phase",y="Proportion of cells")+theme(aspect.ratio=1)+theme(axis.title=element_blank(),
                                                                                 axis.text=element_blank(),
                                                                                 axis.ticks=element_blank())+
  scale_color_brewer(name="Stage",palette="Dark2")



ggplot( zzMAPKSUB[Day==180] , aes(y=CDKN2A ,      x=   (V1_disc)))+
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  geom_jitter(aes( col=as.factor(classification3)),size=2,width = 0.005, height =0)+  
  facet_wrap(ARM~Response,ncol=2)+
  geom_smooth(size=1.5,method="gam",formula= y~s(x,k=5),se=T,col="black")+
  theme_minimal(base_size=18)+
  labs(x="Cell cycle phase",y="Proportion of cells")+theme(aspect.ratio=1)+theme(axis.title=element_blank(),
                                                                                 axis.text=element_blank(),
                                                                                 axis.ticks=element_blank())+
  scale_color_brewer(name="Stage",palette="Dark2")



ggplot( zzMAPKSUB[Day==180] , aes(y=CDK6 ,      x=   (1e-3+V1_disc)))+
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  geom_jitter(aes( col=as.factor(classification3)),size=2,width = 0.005, height =0)+  
  facet_wrap(ARM~Response,ncol=2)+
  geom_smooth(size=1.5,method="gam",formula= y~s(x,k=5),se=T,col="black")+
  theme_minimal(base_size=18)+
  labs(x="Cell cycle phase",y="Proportion of cells")+theme(aspect.ratio=1)+theme(axis.title=element_blank(),
                                                                                 axis.text=element_blank(),
                                                                                 axis.ticks=element_blank())+
  scale_color_brewer(name="Stage",palette="Dark2",labels=c("c","a","b"))+
  #scale_color_discrete(name="Stage")+scale_fill_discrete(name="Stage")+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+
  geom_hline(yintercept=logit(1e-3+1/ 5 ))
#ggplot(  mrg_umap_mapk_tsp2, aes(y=MAPK3 ,   x=  V1_disc,col=as.factor(Day),fill=as.factor(Day)))+
# geom_point(size=2)+geom_smooth(aes(group=interaction(Day)),size=1.5,method="gam",formula=y~s(x,k=3),se=T)
#ggplot(  mrg_umap_mapk_tsp2, aes(y=MAPK10 ,   x=  V1_disc,col=as.factor(Day),fill=as.factor(Day)))+
# geom_point(size=2)+geom_smooth(aes(group=interaction(Day)),size=1.5,method="gam",formula=y~s(x,k=3),se=T)


ggplot(  mrg_umap_mapk_tsp2, aes(y=CDK6 ,   x=  V1_disc,col=as.factor(Day),fill=as.factor(Day)))+
  #coord_polar()+
  #geom_jitter(size=2)+facet_wrap(Day~ARM)+
  theme_classic(base_size=18)+
  geom_smooth(aes(group=interaction(Day)),size=1.5,method="gam",formula=y~s(x,k=3),se=T)+
  theme(aspect.ratio=1)+
  facet_wrap(ARM~Response,ncol=2) 


ggplot(  mrg_umap_mapk_tsp2, aes(y=CDK6 ,      x=  V1,col=Response))+
  #coord_polar()+
  geom_jitter(size=2)+facet_wrap(Day~ARM)+
  theme_classic(base_size=18)+
  geom_smooth(aes(group=interaction(Response),col=Response,fill=Response),size=1.5,method="loess",se=F)+
  theme(aspect.ratio=1)+
  facet_wrap(ARM~Day)   

ggplot(  mrg_umap_mapk_tsp2, aes(y=V2 ,      x=  V1,col=Response))+
  #coord_polar()+
  geom_jitter(size=2)+facet_wrap(Day~ARM)+
  theme_classic(base_size=18)+
  geom_smooth(aes(group=interaction(Response),col=Response,fill=Response),size=1.5,method="loess",se=F)+
  theme(aspect.ratio=1)+
  facet_wrap(ARM~Day)   


ggplot(  mrg_umap_mapk_tsp2[Day==180], aes(y=V2 ,      x=  V1,col=Response))+
  #coord_polar()+
  geom_jitter(size=2)+facet_wrap(ARM~Response)+
  theme_classic(base_size=18)+
  geom_smooth(aes(group=interaction(Response),col=Response,fill=Response),size=1.5,method="loess",se=F)+
  theme(aspect.ratio=1)+
  facet_wrap(ARM~Day)   
#mmm <- data.table(mrg_umap_mapk_tsp2%>%gather(var,val,-c(classification:Phase,V1:dimension2,tour))%>%
#                    group_by(ARM,Patient.Study.ID,classification,tour,Day,Response,var)%>%
#                    dplyr::summarise(val=median(val)))

mmm <- data.table(mrg_umap_mapk_tsp2%>%gather(var,val,-c(Cell.ID:Phase,V1:dimension2,tour))%>%
                    group_by(ARM,Patient.Study.ID,classification,tour,Day,Response,var)%>%
                    dplyr::summarise(val=median(val)))

m0<-lm(CDK6~    CDKN2A*V1_disc, data=mrg_umap_mapk_tsp2[ARM=="B"&Day==180])
ggplot(  mrg_umap_mapk_tsp2, aes(y=Response ,      x=  tour,col=Phase))+
  # coord_polar()+
  geom_jitter(size=2)+facet_wrap(Day~ARM)+
  theme_classic(base_size=18)+
  geom_smooth(aes(group=interaction(Response),col=Response,fill=Response),size=1.5,method="gam",formula=y~s(x,k=3),se=F)+
  theme(aspect.ratio=1)+
  facet_wrap(ARM~Day)   

mrg_umap_mapk_tsp3<-data.table(mrg_umap_mapk_tsp2%>% group_by(ARM,Day,Response,V1_disc)%>%
                                 dplyr::summarise(CDK6=mean(CDK6),CDKN2A=mean(CDKN2A)))

ggplot(  mrg_umap_mapk_tsp3, aes(y=CDK6 ,      x=  V1_disc,col=Response))+
  #coord_polar()+
  geom_jitter(size=2)+facet_wrap(Day~ARM)+
  theme_classic(base_size=18)+
  geom_smooth(aes(group=interaction(Response),col=Response,fill=Response),size=1.5,method="gam",formula=y~s(x,k=3),se=F)+
  theme(aspect.ratio=1)+
  facet_wrap(ARM~Day)   
ggplot(  mrg_umap_mapk_tsp3, aes(y=CDKN2A ,      x=  V1_disc,col=Response))+
  #coord_polar()+
  geom_jitter(size=2)+facet_wrap(Day~ARM)+
  theme_classic(base_size=18)+
  geom_smooth(aes(group=interaction(Response),col=Response,fill=Response),size=1.5,method="gam",formula=y~s(x,k=3),se=F)+
  theme(aspect.ratio=1)+
  facet_wrap(ARM~Day)   
ggplot(  mrg_umap_mapk_tsp3, aes(y=CDK6 ,      x=  V1_disc,col=CDKN2A))+
  #coord_polar()+
  geom_jitter(size=2)+facet_wrap(Day~ARM)+
  theme_classic(base_size=18)+
  geom_smooth(aes(group=interaction(Response),col=Response,fill=Response),size=1.5,method="gam",formula=y~s(x,k=3),se=F)+
  theme(aspect.ratio=1)+
  facet_wrap(ARM~Day) 

cor(mrg_umap_mapk_tsp2$mapk_axisV1,mrg_umap_mapk_tsp2%>%dplyr::select(CCND1:V2))
ggplot(  mrg_umap_mapk_tsp2, aes(y=CDK6 ,      x=  V1_disc,col=Response))+
  #coord_polar()+
  geom_jitter(size=2)+facet_wrap(Day~ARM)+
  theme_classic(base_size=18)+
  geom_smooth(aes(group=interaction(Response),col=Response,fill=Response),size=1.5,method="loess",se=F)+
  theme(aspect.ratio=1)+
  facet_wrap(ARM~Day)   



ggplot(  mmm[var=="CDK6"], aes(x=tour /max(tour) ,      y=  val,col=Response))+
  coord_polar()+
  geom_jitter(size=2)+facet_wrap(Day~ARM)+
  theme_classic(base_size=18)+
  geom_smooth(aes(group=interaction(Response),col=Response,fill=Response),size=1.5,method="loess",se=F)+
  theme(aspect.ratio=1)+
  facet_wrap(ARM~Day)   


ggplot(  mmm[var=="CDK6"], aes(x=tour /max(tour) ,      y=  val,col=Response))+
  coord_polar()+
  geom_jitter(size=2)+facet_wrap(Day~ARM)+
  theme_classic(base_size=18)+
  geom_smooth(aes(group=interaction(Response),col=Response,fill=Response),size=1.5,method="loess",se=F)+
  theme(aspect.ratio=1)+
  facet_wrap(ARM~Day)   



ggplot(  mmm[var=="CDK6"], aes(x=tour /max(tour) ,      y=  val,col=Response))+
  coord_polar()+
  geom_jitter(size=2)+facet_wrap(Day~ARM)+
  theme_classic(base_size=18)+
  geom_smooth(aes(group=interaction(Response),col=Response,fill=Response),size=1.5,method="loess",se=F)+
  theme(aspect.ratio=1)+
  facet_wrap(ARM~Day)              

ggplot(  mmm[var=="RB1"], aes(x=tour /max(tour) ,      y=  val,col=Response))+
  geom_point(size=2)+facet_wrap(Day~ARM)+
  theme_classic(base_size=18)+
  geom_smooth(aes(group=interaction(Patient.Study.ID,Response),col=Response,fill=Response),size=1.5,method="loess",se=F)+
  theme(aspect.ratio=1)+
  facet_wrap(ARM~Day)                

ggplot(  mmm[var=="MAPK10"], aes(x=tour /max(tour) ,      y=  val,col=Response))+
  geom_point(size=2)+facet_wrap(Day~ARM)+
  theme_classic(base_size=18)+
  geom_smooth(aes(group=interaction(Response),col=Response,fill=Response),size=1.5,method="loess",se=F)+
  theme(aspect.ratio=1)+
  facet_wrap(ARM~Day)                


ggplot(  mrg_umap_mapk_tsp, aes(x=MAPK10,      y=  CCND1))+
  geom_point(size=2)+
  theme_classic(base_size=18)+
  geom_smooth(aes(group=interaction(Response),col=Response,fill=Response),size=1.5,method="lm",se=F)+
  theme(aspect.ratio=1)+
  facet_wrap(ARM~Day)                




ggplot(  mmm[var=="CCND1"], aes(x=tour /max(tour) ,      y=  val,col=Response))+
  geom_point(size=2)+facet_wrap(Day~ARM)+
  theme_classic(base_size=18)+
  geom_smooth(aes(group=interaction(Patient.Study.ID,Response),col=Response,fill=Response),size=1.5,method="loess",se=F)+
  theme(aspect.ratio=1)+
  facet_wrap(ARM~Day)  


ggplot(  mmm[var=="CCNH"], aes(x=tour /max(tour) ,      y=  val,col=Response))+
  geom_point(size=2)+facet_wrap(Day~ARM)+
  theme_classic(base_size=18)+
  geom_smooth(aes(group=interaction(Patient.Study.ID,Response),col=Response,fill=Response),size=1.5,method="loess",se=F)+
  theme(aspect.ratio=1)+
  facet_wrap(ARM~Day)  





