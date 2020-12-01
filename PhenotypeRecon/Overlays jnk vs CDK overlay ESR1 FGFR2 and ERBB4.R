load(file="/Users/jason/Dropbox/FELINE Project/Data_analysis/scRNA/15_subclone_interaction/merging rtk mapk and cell cycle.RData")

ggplot( mrg_umap_mapk_tsp2[Day==180]%>%group_by(V1_disc)%>%summarise(FGFR2=median(FGFR2),ESR1=median(ESR1),CDK6=median(CDK6)), aes(y=log(FGFR2-min(FGFR2)) ,      x=   (V1_disc),col=log(ESR1-min(ESR1)) )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  #facet_wrap(ARM~.,ncol=3)+
  geom_smooth(size=1.5, method="gam",formula= y~s(x,k=4 ),se=T,col="firebrick",fill="firebrick")+
  geom_point(size=2)

require(viridis)
# ggplot( data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9]%>%group_by(ARM,V1_disc2)%>%summarise(FGFR2=mean(FGFR2),ESR1=median(ESR1),CDK6=mean(CDK6)))[CDK6<2.4], 
#         aes(y=(CDK6) ,      x=   (V1_disc2),col=log(FGFR2-min(FGFR2)) )) +
#   #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
#   facet_wrap(~ARM,ncol=3)+
#   geom_smooth(size=2.5, method="gam",formula= y~s(x,k=4), se=T, col="firebrick", fill="firebrick")+
#   geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+labs(y="CDK6 ",x="JNK activation score")+
#   scale_color_viridis(option="B",begin = 0.1,name="FGFR2")+
#   theme(aspect.ratio=1,
#         strip.background = element_blank(),
#         strip.text.x = element_blank()
#   )

ggplot(ttt[Day==180],aes(y=CDK6,x=Patient.Study.ID,group=Patient.Study.ID))+facet_wrap(~ARM,ncol=1)+geom_violin()+ theme(axis.text.x = element_text(angle = 90))
ff<-mrg_umap_mapk_tsp2#[ARM!="A"]
ff[,CDK6_mu:=0]
ff[Patient.Study.ID%in%c("001-118","001-119","001-143","2972-005-202","2972-005-202"),CDK6_mu:=1]

mean_ESR1lev <- data.table(ff %>%group_by(CDK6_mu,Patient.Study.ID,ARM,Day,Response)%>%summarise(ESR1=median(ESR1) ))

pc <- mrg_umap_mapk_tsp2[ARM=="C"]$Patient.Study.ID%>%unique()

####GOOD
bb<-"steelblue"
newfin1<-ggplot( data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9]
                   [!Patient.Study.ID%in%c(pc[c(2,7,11)],"001-137","2972-005-202")]#[!Patient.Study.ID%in%c("001-118")]   #6
                   #[!Patient.Study.ID%in%c("001-118","001-119","001-143")]
                   %>%group_by(ARM,V1_disc2)%>%summarise(ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4], 
        aes(y=(CDK6) ,      x=   (V1_disc2),col=log(1+ESR1-min(ESR1)) )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(~ARM,ncol=3)+
  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=6), method.args=list(gamma=4),se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+labs(y="CDK6 ",x="JNK activation score")+
  scale_color_viridis(option="B",begin = 0.1,name="ESR1")+
  theme(aspect.ratio=1,
        strip.background = element_blank(),
        strip.text.x = element_blank()
  ) +
  geom_smooth(data=
                data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9][Patient.Study.ID%in%c(pc[c(2,7,11)],"001-137","2972-005-202")]#6
                           %>%group_by(ARM,V1_disc2)%>%summarise(ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4],size=2.5, method="gam",formula= y~s(x,k=9), method.args=list(gamma=10), se=T, col=bb, fill=bb)+
  
  geom_point(data=
               data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9][Patient.Study.ID%in%c(pc[c(2,7,11)],"001-137","2972-005-202") ]#6
                          %>%group_by(ARM,V1_disc2)%>%summarise(ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4],shape=2#pch=21,aes(fill=ERBB4)
  )
newfin1

#ggsave(newfin1, filename = "/Users/jason/Desktop/ESR1 overlay of CDK vs JNK.png",height=3,width=9)

#ggsave(newfin1, filename = "/Users/jason/Dropbox/FELINE Project/Manuscript/Final Main Figures/Blank figures/ESR1 overlay of CDK vs JNK.png",height=3,width=9)

newfin2<-ggplot( data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9]
                            [!Patient.Study.ID%in%c(pc[c(2,7,11)],"001-137","2972-005-202")]#[!Patient.Study.ID%in%c("001-118")]   #6
                            #[!Patient.Study.ID%in%c("001-118","001-119","001-143")]
                            %>%group_by(ARM,V1_disc2)%>%summarise(ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4], 
                 aes(y=(CDK6) ,      x=   (V1_disc2),col=log(1+ERBB4-min(ERBB4)) )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(~ARM,ncol=3)+
  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=6), method.args=list(gamma=4),se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+labs(y="CDK6 ",x="JNK activation score")+
  scale_color_viridis(option="B",begin = 0.1,name="ERBB4")+
  theme(aspect.ratio=1,
        strip.background = element_blank(),
        strip.text.x = element_blank()
  ) +
  geom_smooth(data=
                data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9][Patient.Study.ID%in%c(pc[c(2,7,11)],"001-137","2972-005-202")]#6
                           %>%group_by(ARM,V1_disc2)%>%summarise(ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4],size=2.5, method="gam",formula= y~s(x,k=9), method.args=list(gamma=10), se=T, col="firebrick2", fill="firebrick2")+
  
  geom_point(data=
               data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9][Patient.Study.ID%in%c(pc[c(2,7,11)],"001-137","2972-005-202")]#6
                          %>%group_by(ARM,V1_disc2)%>%summarise(ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4],shape=2,size=2#pch=21,aes(fill=ERBB4)
  )
ggsave(newfin2, filename = "/Users/jason/Desktop/ERBB4 overlay of CDK vs JNK.png",height=3,width=9)

newfin3<-ggplot( data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9]
                            [!Patient.Study.ID%in%c(pc[c(2,7,11)],"001-137","2972-005-202")]#[!Patient.Study.ID%in%c("001-118")]   #6
                            #[!Patient.Study.ID%in%c("001-118","001-119","001-143")]
                            %>%group_by(ARM,V1_disc2)%>%summarise(ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4], 
                 aes(y=(CDK6) ,      x=   (V1_disc2),col=log(1+FGFR2-min(FGFR2)) )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(~ARM,ncol=3)+
  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=6), method.args=list(gamma=4),se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+labs(y="CDK6 ",x="JNK activation score")+
  scale_color_viridis(option="B",begin = 0.1,name="FGFR2")+
  theme(aspect.ratio=1,
        strip.background = element_blank(),
        strip.text.x = element_blank()
  ) +
  geom_smooth(data=
                data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9][Patient.Study.ID%in%c(pc[c(2,7,11)],"001-137","2972-005-202")]#6
                           %>%group_by(ARM,V1_disc2)%>%summarise(ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4],size=2.5, method="gam",formula= y~s(x,k=9), method.args=list(gamma=10), se=T, col="firebrick2", fill="firebrick2")+
  
  geom_point(data=
               data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9][Patient.Study.ID%in%c(pc[c(2,7,11)],"001-137","2972-005-202")]#6
                          %>%group_by(ARM,V1_disc2)%>%summarise(ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4],shape=2,size=2#pch=21,aes(fill=ERBB4)
  )
ggsave(newfin3, filename = "/Users/jason/Desktop/FGFR2 overlay of CDK vs JNK.png",height=3,width=9)

# BLANKS

Bnewfin1<-ggplot( data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9]
                            [!Patient.Study.ID%in%c(pc[c(2,7,11)],"001-137","2972-005-202")]#[!Patient.Study.ID%in%c("001-118")]   #6
                            #[!Patient.Study.ID%in%c("001-118","001-119","001-143")]
                            %>%group_by(ARM,V1_disc2)%>%summarise(ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4], 
                 aes(y=(CDK6) ,      x=   (V1_disc2),col=log(1+ESR1-min(ESR1)) )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(~ARM,ncol=3)+
  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=6), method.args=list(gamma=4),se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+labs(y="CDK6 ",x="JNK activation score")+
  scale_color_viridis(option="B",begin = 0.1,name="ESR1")+
  theme(aspect.ratio=1,
        strip.background = element_blank(),
        strip.text.x = element_blank()
  ) +
  geom_smooth(data=
                data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9][Patient.Study.ID%in%c(pc[c(2,7,11)],"001-137","2972-005-202")]#6
                           %>%group_by(ARM,V1_disc2)%>%summarise(ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4],size=2.5, method="gam",formula= y~s(x,k=9), method.args=list(gamma=10), se=T, col=bb, fill=bb)+
  
  geom_point(data=
               data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9][Patient.Study.ID%in%c(pc[c(2,7,11)],"001-137","2972-005-202") ]#6
                          %>%group_by(ARM,V1_disc2)%>%summarise(ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4],shape=2,stroke=1#pch=21,aes(fill=ERBB4)
  )+theme(aspect.ratio=1,
          strip.background = element_blank(),
          strip.text.x = element_blank()
  )+ 
  theme(axis.title=element_blank(), axis.text=element_blank(),legend.title=element_blank(),legend.text=element_blank()) +theme(legend.position = "none")

ggsave(Bnewfin1, filename = "/Users/jason/Dropbox/FELINE Project/Manuscript/Final Main Figures/Blank figures/ESR1 overlay of CDK vs JNK.png",height=3,width=9)
#ggsave(newfin1, filename = "/Users/jason/Dropbox/FELINE Project/Manuscript/Final Main Figures/Blank figures/ESR1 overlay of CDK vs JNK.png",height=3,width=9)

Bnewfin2<-ggplot( data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9]
                            [!Patient.Study.ID%in%c(pc[c(2,7,11)],"001-137","2972-005-202")]#[!Patient.Study.ID%in%c("001-118")]   #6
                            #[!Patient.Study.ID%in%c("001-118","001-119","001-143")]
                            %>%group_by(ARM,V1_disc2)%>%summarise(ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4], 
                 aes(y=(CDK6) ,      x=   (V1_disc2),col=log(1+ERBB4-min(ERBB4)) )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(~ARM,ncol=3)+
  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=6), method.args=list(gamma=4),se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+labs(y="CDK6 ",x="JNK activation score")+
  scale_color_viridis(option="B",begin = 0.1,end = 0.885,name="ERBB4")+
  theme(aspect.ratio=1,
        strip.background = element_blank(),
        strip.text.x = element_blank()
  ) +
  geom_smooth(data=
                data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9][Patient.Study.ID%in%c(pc[c(2,7,11)],"001-137","2972-005-202")]#6
                           %>%group_by(ARM,V1_disc2)%>%summarise(ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4],size=2.5, method="gam",formula= y~s(x,k=9), method.args=list(gamma=10), se=T, col=bb, fill=bb)+
  
  geom_point(data=
               data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9][Patient.Study.ID%in%c(pc[c(2,7,11)],"001-137","2972-005-202")]#6
                          %>%group_by(ARM,V1_disc2)%>%summarise(ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4],shape=2,size=2,stroke=1#pch=21,aes(fill=ERBB4)
  )+theme(aspect.ratio=1,
          strip.background = element_blank(),
          strip.text.x = element_blank()
  )+ 
  theme(axis.title=element_blank(), axis.text=element_blank(),legend.title=element_blank(),legend.text=element_blank()) +theme(legend.position = "none")

ggsave(Bnewfin2, filename = "/Users/jason/Dropbox/FELINE Project/Manuscript/Final Main Figures/Blank figures/ERBB4 overlay of CDK vs JNK.png",height=3,width=9)

Bnewfin3<-ggplot( data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9]
                            [!Patient.Study.ID%in%c(pc[c(2,7,11)],"001-137","2972-005-202")]#[!Patient.Study.ID%in%c("001-118")]   #6
                            #[!Patient.Study.ID%in%c("001-118","001-119","001-143")]
                            %>%group_by(ARM,V1_disc2)%>%summarise(ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4], 
                 aes(y=(CDK6) ,      x=   (V1_disc2),col=log(1+FGFR2-min(FGFR2)) )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(~ARM,ncol=3)+
  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=6), method.args=list(gamma=4),se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+labs(y="CDK6 ",x="JNK activation score")+
  scale_color_viridis(option="B",begin = 0.1,name="FGFR2")+
  theme(aspect.ratio=1,
        strip.background = element_blank(),
        strip.text.x = element_blank()
  ) +
  geom_smooth(data=
                data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9][Patient.Study.ID%in%c(pc[c(2,7,11)],"001-137","2972-005-202")]#6
                           %>%group_by(ARM,V1_disc2)%>%summarise(ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4],size=2.5, method="gam",formula= y~s(x,k=9), method.args=list(gamma=10), se=T, col=bb, fill=bb)+
  
  geom_point(data=
               data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9][Patient.Study.ID%in%c(pc[c(2,7,11)],"001-137","2972-005-202")]#6
                          %>%group_by(ARM,V1_disc2)%>%summarise(ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4],shape=2,size=2,stroke=1#pch=21,aes(fill=ERBB4)
  )+theme(aspect.ratio=1,
          strip.background = element_blank(),
          strip.text.x = element_blank()
  )+ 
  theme(axis.title=element_blank(), axis.text=element_blank(),legend.title=element_blank(),legend.text=element_blank()) +theme(legend.position = "none")

ggsave(Bnewfin3, filename = "/Users/jason/Dropbox/FELINE Project/Manuscript/Final Main Figures/Blank figures/FGFR2 overlay of CDK vs JNK.png",height=3,width=9)


legend <- cowplot::get_legend(newfin1+
                                theme(axis.title=element_blank(), axis.text=element_blank(),legend.title=element_blank(),legend.text=element_blank()) )
library(grid)
library(gridExtra) 
grid.newpage()
grid.draw(legend)


legend <- cowplot::get_legend(newfin2+
                                theme(axis.title=element_blank(), axis.text=element_blank(),legend.title=element_blank(),legend.text=element_blank()) )
library(grid)
library(gridExtra) 
grid.newpage()
grid.draw(legend)



legend <- cowplot::get_legend(newfin3+
                                theme(axis.title=element_blank(), axis.text=element_blank(),legend.title=element_blank(),legend.text=element_blank()) )
library(grid)
library(gridExtra) 
grid.newpage()
grid.draw(legend)





mrg_umap_mapk_tsp2 [, CDK6_CNA:="Normal"]
mrg_umap_mapk_tsp2 [ Patient.Study.ID%in%c(pc[c(2,7,11)],"001-137","2972-005-202") , CDK6_CNA:="Gain" ]
yy<- data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9]
               %>%
                 group_by(CDK6_CNA,ARM,V1_disc2)%>%summarise(ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4]
yy[,CDK6_CNA_ARM:=paste0(CDK6_CNA, ARM)]
gm<-gam(CDK6~s(V1_disc2,by=as.factor(CDK6_CNA_ARM),k=6),data=yy[ARM=="B"])
cor(predict(gm),yy[ARM=="B"]$CDK6)

plot(gm)
summary(gm)
lm1<-lm(CDK6~V1_disc2*CDK6_CNA,data=yy[ARM!="A"])

lm1<-lm(CDK6~V1_disc2,data=yy[CDK6_CNA!="Normal"][ARM!="A"])
lm1<-lm(CDK6~ESR1*CDK6_CNA,data=yy[ARM!="A"])
lm1<-lm(CDK6~ESR1,data=yy[CDK6_CNA!="Normal"][ARM!="A"])
summary(lm1)
lm1<-lm(ESR1~V1_disc2,data=yy[CDK6_CNA=="Normal"][ARM!="A"])
summary(lm1)
summary(gm)

lm1<-lm(ERBB4~ARM*CDK6_CNA,data=yy[ARM!="A"])



c(pc[c(2,7,11)],"001-137","2972-005-202")

erbb4data<-ttt3[Variable=="ERBB4"]
ggplot( erbb4data[Patient.Study.ID%in%c(pc[c(2,7,11)],"001-137","2972-005-202")],aes(y=Expression- median(erbb4data$Expression),
                                                                                     x=log(1+Day),group=Day))+geom_violin()+facet_wrap(~Patient.Study.ID)

mean(erbb4data$Expression)

FGFR2data<-ttt3[Variable=="FGFR2"]
FGFR2data[,CNA:=FALSE]
FGFR2data[Patient.Study.ID%in%c("2972-007-703","001-145"   ,"001-135"   ,"001-101"   ,"001-102","001-117")]%>%dplyr::select(Patient.Study.ID, ARM, CNA)%>%unique()
FGFR2data[ARM=="B"]%>%dplyr::select(Patient.Study.ID, ARM)%>%unique()
(FGFR2data[Patient.Study.ID%in%
            unique(FGFR2data[ARM=="B"]$Patient.Study.ID)
            ]%>%dplyr::select(Patient.Study.ID, ARM, CNA)%>%unique())[order(Patient.Study.ID)]

ggplot( FGFR2data[ARM=="B"],aes(y=Expression- median(FGFR2data$Expression),
             x=log(1+Day),group=Day))+geom_violin()+facet_wrap(~Patient.Study.ID)

FGFR2data[Patient.Study.ID%in%c(pc[c(2,7,11)],"001-137","2972-005-202"),CNA:=TRUE]
ggplot( FGFR2data[Patient.Study.ID%in%c(pc[c(2,7,11)],"001-137","2972-005-202")],aes(y=Expression- median(FGFR2data$Expression),
                                                                                     x=log(1+Day),group=Day))+geom_violin()+facet_wrap(~Patient.Study.ID)

ggplot( FGFR2data[][ARM!="A"],aes(y=Expression- median(FGFR2data$Expression),
                                                                                     x=log(1+Day),group=Day,fill=CNA))+geom_violin()+facet_wrap(~Patient.Study.ID)

ggplot( FGFR2data[][],aes(y=Expression- median(FGFR2data$Expression),
                                  x=log(1+Day),group=Day,fill=CNA))+geom_violin()+facet_wrap(~Patient.Study.ID)



mean(erbb4data$Ex
     
     
     
ggplot(mean_ESR1lev, aes(y=ESR1,x=Day,group=Patient.Study.ID))+geom_point()+geom_line()+facet_wrap(Response~ARM)

ggplot( data.table(ff[Day==180][abs(V1_disc2)<1.9]%>%group_by(CDK6_mu,Patient.Study.ID,ARM,V1_disc2)%>%summarise(FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)) 
)[CDK6<2.4], 
aes(col=(CDK6) ,  group= Patient.Study.ID,   x=   (V1_disc2),y=log(1+(ESR1)-min((ESR1))) )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(CDK6_mu~ARM,ncol=7)+
  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=4), se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+ylim(c(1.5,3.5))


ggplot( data.table(ff[Day==180][abs(V1_disc2)<1.9]%>%group_by(CDK6_mu,Patient.Study.ID,ARM,V1_disc2)%>%summarise(FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)) 
)[CDK6<2.4], 
aes(col=(CDK6) ,  group= Patient.Study.ID,   x=   (V1_disc2),y=log(1+(ESR1)-min((ESR1))) )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(CDK6_mu~ARM,ncol=7)+
  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=4), se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+ylim(c(1.5,3.5))

#,meanlnESR1=mean(log(1+ESR1-min(ESR1)))
ggplot( data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9][]%>%group_by(Patient.Study.ID,ARM,V1_disc2)%>%summarise(FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)) %>%group_by(Patient.Study.ID,ARM)%>%
                     mutate(meanESR1=mean(ESR1) )
)[CDK6<2.4], 
aes(col=(Patient.Study.ID) , group=Patient.Study.ID,  x=   (V1_disc2),y=ESR1-meanESR1 )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(~ARM,ncol=3,scales="free_y")+
  geom_smooth(size=2.5, method="lm",formula= y~s(x,k=4),method.args=list(gamma=7), se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)


ggplot( data.table(ff[Day==180][abs(V1_disc2)<1.9]%>%group_by(CDK6_mu,Patient.Study.ID,Subclone,ARM,V1_disc2)%>%summarise(FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4], 
        aes(y=(CDK6) ,  group= interaction(Subclone,Patient.Study.ID),   x=   (V1_disc2),col=log(ESR1-min(ESR1)) )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(CDK6_mu~ARM,ncol=3)+
  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=4), se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)

ggplot( data.table(ff[Day==180][abs(V1_disc)<1.9]%>%group_by(CDK6_mu,Patient.Study.ID,Subclone,ARM,V1_disc)%>%summarise(FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4], 
        aes(y=(CDK6) ,  group= interaction(ARM,Subclone,Patient.Study.ID),  
            x=   (V1_disc),col=log(1+ESR1-min(ESR1)) )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(~ARM,ncol=3)+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+
  
  geom_smooth(aes( group= interaction(ARM,Patient.Study.ID) ) ,size=2.5, method="gam",formula= y~s(x,k=3), se=T, col="firebrick", fill="firebrick")

ggplot( data.table(ff[ARM=="C"][Day==180][abs(V1_disc2)<1.9]%>%group_by(CDK6_mu,Patient.Study.ID,Subclone,ARM,V1_disc2)%>%summarise(FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4], 
        aes(y=(CDK6) ,  group= interaction(ARM,Subclone,Patient.Study.ID),  
            x=   (V1_disc2),col=log(1+ESR1-min(ESR1)) )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(~Patient.Study.ID,ncol=3)+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+
  
  geom_smooth(aes( group= interaction(ARM,Patient.Study.ID) ) ,size=2.5, method="gam",formula= y~s(x,k=3), se=T, col="firebrick", fill="firebrick")+
  scale_color_viridis(option="B",begin = 0.1,name="ESR1")


ggplot( data.table(ff[ARM=="C"][Day==180][abs(V1_disc2)<1.9]%>%group_by(CDK6_mu,Patient.Study.ID,Subclone,ARM,V1_disc2)%>%summarise(ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4], 
        aes(y=(CDK6) ,  group= interaction(ARM,Subclone,Patient.Study.ID),  
            x=   (V1_disc2),col=log(1+ERBB4-min(ERBB4)) )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(~Patient.Study.ID,ncol=3)+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+
  
  geom_smooth(aes( group= interaction(ARM,Patient.Study.ID) ) ,size=2.5, method="gam",formula= y~s(x,k=3), se=T, col="firebrick", fill="firebrick")+
  scale_color_viridis(option="B",begin = 0.1,name="ERBB4")



ggplot( data.table(ff[Day==180][abs(V1_disc2)<1.9][ARM=="C"]%>%group_by(CDK6_mu,Patient.Study.ID,Response,ARM,V1_disc2)%>%summarise(FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4], 
        aes(y=(CDK6) ,  group= Patient.Study.ID,   x=   (V1_disc2),col=CDK6_mu )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(Response~Patient.Study.ID,ncol=3)+
  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=4), se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)




ggplot( data.table(ff[Day==180][abs(V1_disc2)<1.9]%>%group_by(CDK6_mu,ARM,V1_disc2)%>%summarise(FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDKN2A=mean(CDKN2A),CDK6=mean(CDK6)))[CDK6<2.4], 
        aes(col=(CDK6_mu) ,  group=CDK6_mu,    x=   (V1_disc2),y=log(ESR1-min(ESR1)) )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(~ARM,ncol=3)+
  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=3), se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)




+labs(y="CDK6 ",x="JNK activation score")+
  scale_color_viridis(option="B",begin = 0.1,name="ESR1")+
  theme(aspect.ratio=1,
        strip.background = element_blank(),
        strip.text.x = element_blank()
  ) +
  theme(axis.title=element_blank(), axis.text=element_blank(),legend.title=element_blank(),legend.text=element_blank()) +theme(legend.position = "none")


ggplot( data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9][!(ARM=="C"&!Patient.Study.ID%in%c("001-118","001-119","001-143"))]%>%group_by(ARM,V1_disc2)%>%summarise(FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4], 
        aes(y=(CDK6) ,      x=   (V1_disc2),col=log(ESR1-min(ESR1)) )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(~ARM,ncol=3)+
  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=4), se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+labs(y="CDK6 ",x="JNK activation score")+
  scale_color_viridis(option="B",begin = 0.1,name="ESR1")



ggplot( data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9][!Patient.Study.ID%in%c("001-118","001-119","001-143")]%>%group_by(Patient.Study.ID,ARM,V1_disc2)%>%summarise(FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4], 
        aes(y=(CDK6) ,   group=Patient.Study.ID,   x=   (V1_disc2),col=log(ESR1-min(ESR1)) )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(~ARM,ncol=3)+
  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=4), se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+labs(y="CDK6 ",x="JNK activation score")+
  scale_color_viridis(option="B",begin = 0.1,name="ESR1")


ggplot( data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9][]%>%group_by(Patient.Study.ID,ARM,V1_disc2)%>%summarise(FGFR2=mean(FGFR2),ESR1=mean(ESR1),ERBB4=mean(ERBB4),CDK6=mean(CDK6)))[CDK6<2.4], 
        aes(y=(CDK6) ,   group=Patient.Study.ID,   x=   (V1_disc2),col=log(ERBB4-min(ERBB4)) )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(~ARM,ncol=3)+
  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=4), se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+labs(y="CDK6 ",x="JNK activation score")+
  scale_color_viridis(option="B",begin = 0.1,name="ERBB4")

ggplot( data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9][]%>%group_by(Patient.Study.ID,ARM,V1_disc2)%>%summarise(FGFR2=mean(FGFR2),ESR1=mean(ESR1),ERBB4=mean(ERBB4),CDK6=mean(CDK6)))[CDK6<2.4], 
        aes(y=(CDK6) ,   group=Patient.Study.ID,   x=   (V1_disc2),col=log(FGFR2-min(FGFR2)) )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(~ARM,ncol=3)+
  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=4), se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+labs(y="CDK6 ",x="JNK activation score")+
  scale_color_viridis(option="B",begin = 0.1,name="FGFR2")

ggplot( data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9][]%>%group_by(Patient.Study.ID,ARM,Response,V1_disc2)%>%summarise(MAPK8=mean(MAPK8),MAPK10=mean(MAPK10),FGFR2=mean(FGFR2),ESR1=mean(ESR1),ERBB4=mean(ERBB4),CDK6=mean(CDK6),CDKN2A=mean(CDKN2A) ))[CDK6<2.4], 
        aes(y=(CDK6) ,   group=Patient.Study.ID,   x=   (V1_disc2),col=Response )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(~ARM,ncol=3)+
  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=4), se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+labs(y="CDK6 ",x="JNK activation score")#+
# scale_color_viridis(option="B"
ggplot( data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9][]%>%group_by(Patient.Study.ID,ARM,Response,V1_disc2)%>%summarise(RB1=mean(RB1),MAPK8=mean(MAPK8),MAPK10=mean(MAPK10),FGFR2=mean(FGFR2),ESR1=mean(ESR1),ERBB4=mean(ERBB4),CDK6=mean(CDK6),CDKN2A=mean(CDKN2A) ))[CDK6<2.4], 
        aes(y=log(ERBB4-min(ERBB4)) ,   group=Patient.Study.ID,   x=   log(ESR1-min(ESR1)),col=Response )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(~ARM,ncol=3)+
  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=4), se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)#+#labs(y="CDK6 ",x="JNK activation score")#+
ggplot( data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9][]%>%group_by(Patient.Study.ID,ARM,Response,V1_disc2)%>%summarise(RB1=mean(RB1),MAPK8=mean(MAPK8),MAPK10=mean(MAPK10),FGFR2=mean(FGFR2),ESR1=mean(ESR1),ERBB4=mean(ERBB4),CDK6=mean(CDK6),CDKN2A=mean(CDKN2A) ))[CDK6<2.4], 
        aes(y=log(FGFR2-min(FGFR2)) ,   group=Patient.Study.ID,   x=   log(ESR1-min(ESR1)),col=Response )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(~ARM,ncol=3)+
  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=4), se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)#+#labs(y="CDK6 ",x="JNK activation score")#+

ggplot( data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9][]%>%group_by(Patient.Study.ID,ARM,Response,V1_disc2)%>%summarise(RORC=mean(RORC),RORA=mean(RORA),RB1=mean(RB1),MAPK8=mean(MAPK8),MAPK10=mean(MAPK10),FGFR2=mean(FGFR2),ESR1=mean(ESR1),ERBB4=mean(ERBB4),CDK6=mean(CDK6),CDKN2A=mean(CDKN2A) ))[CDK6<2.4], 
        aes(y=log(RORC-min(RORC)) ,   group=Patient.Study.ID,   x=   log(ESR1-min(ESR1)),col=Response )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(~ARM,ncol=3)+
  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=4), se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)#+#labs(y="CDK6 ",x="JNK activation score")#+

ggplot( data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9][!Patient.Study.ID%in%c("001-118","001-119","001-143")]%>%group_by(Patient.Study.ID,ARM,V1_disc2)%>%summarise(RORC=mean(RORC),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4], 
        aes(y=(CDK6) ,    group=Patient.Study.ID,  x=   (V1_disc2),col=log(ESR1-min(ESR1)) )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(~ARM,ncol=3)+
  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=4), se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+labs(y="CDK6 ",x="JNK activation score")+
  scale_color_viridis(option="B",begin = 0.1,name="RORC")


ggplot( data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9][]%>%group_by(Patient.Study.ID,Response,ARM,V1_disc2)%>%summarise(RORC=mean(RORC),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4], 
        aes(y=log(ESR1-min(ESR1)) ,    group=Patient.Study.ID,  x=   (V1_disc2),fill=Response,col=log(ESR1-min(ESR1)) )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(~ARM,ncol=3)+
  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=4), se=T, col="firebrick",alpha=0.8)+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+
  scale_color_viridis(option="B",begin = 0.1)



ggplot( data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9][ARM=="C"]
                   [!Patient.Study.ID%in%c("001-118")]
                   #[!Patient.Study.ID%in%c("001-118","001-119","001-143")]
                   %>%group_by(ARM,Patient.Study.ID,V1_disc2)%>%summarise(ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4], 
        aes(col=(CDK6) ,    group=Patient.Study.ID,  x=   (V1_disc2),y=log(ESR1-min(ESR1)) )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(~Patient.Study.ID,ncol=3)+
  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=4), se=T, col="firebrick", fill="firebrick")





# scale_color_viridis(option="B",begin = 0.1,name="MAPK8")
ggplot( data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9]
                   [!Patient.Study.ID%in%c("001-118","001-143","2972-007-701")]
                   #[!Patient.Study.ID%in%c("001-118","001-119","001-143")]
                   %>%group_by(ARM,V1_disc2)%>%summarise(ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4], 
        aes(y=(CDK6) ,      x=   (V1_disc2),col=log(ESR1-min(ESR1)) )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(~ARM,ncol=3)+
  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=4), se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+labs(y="CDK6 ",x="JNK activation score")+
  scale_color_viridis(option="B",begin = 0.1,name="ESR1")+
  theme(aspect.ratio=1,
        strip.background = element_blank(),
        strip.text.x = element_blank()
  ) +
  theme(axis.title=element_blank(), axis.text=element_blank(),legend.title=element_blank(),legend.text=element_blank()) +theme(legend.position = "none")




ggplot( data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9]
                   [!Patient.Study.ID%in%pc[c(2,7)]]#[!Patient.Study.ID%in%c("001-118")]   #6
                   #[!Patient.Study.ID%in%c("001-118","001-119","001-143")]
                   %>%group_by(ARM,V1_disc2)%>%summarise(FOS=mean(FOS),ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4], 
        aes(y=(CDK6) ,      x=   (V1_disc2),col=log(1+FOS-min(FOS)) )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(~ARM,ncol=3)+
  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=4), se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+labs(y="CDK6 ",x="JNK activation score")+
  scale_color_viridis(option="B",begin = 0.1,name="FOS")+
  theme(aspect.ratio=1,
        strip.background = element_blank(),
        strip.text.x = element_blank()
  ) +
  geom_smooth(data=
                data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9][Patient.Study.ID%in%pc[c(2,7,6)]]#6
                           %>%group_by(ARM,V1_disc2)%>%summarise(FOS=mean(FOS),ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4],size=2.5, method="gam",formula= y~s(x,k=3), se=T, col="firebrick2", fill="firebrick2")+
  
  geom_point(data=
               data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9][Patient.Study.ID%in%pc[c(2,7,6)]]#6
                          %>%group_by(ARM,V1_disc2)%>%summarise(FOS=mean(FOS),ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4],shape=2#pch=21,aes(fill=ERBB4)
  )


ggplot( data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9]
                   [!Patient.Study.ID%in%pc[c(2,7)]]#[!Patient.Study.ID%in%c("001-118")]   #6
                   #[!Patient.Study.ID%in%c("001-118","001-119","001-143")]
                   %>%group_by(ARM,V1_disc2)%>%summarise(JUNB=mean(JUNB),ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4], 
        aes(y=(CDK6) ,      x=   (V1_disc2),col=log(1+FGFR2-min(FGFR2)) )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(~ARM,ncol=3)+
  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=4), se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+labs(y="CDK6 ",x="JNK activation score")+
  scale_color_viridis(option="B",begin = 0.1,name="FGFR2")+
  theme(aspect.ratio=1,
        strip.background = element_blank(),
        strip.text.x = element_blank()
  ) +
  geom_smooth(data=
                data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9][Patient.Study.ID%in%pc[c(2,7,6)]]#6
                           %>%group_by(ARM,V1_disc2)%>%summarise(FGFR2=mean(FGFR2),ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4],size=2.5, method="gam",formula= y~s(x,k=3), se=T, col="firebrick2", fill="firebrick2")+
  
  geom_point(data=
               data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9][Patient.Study.ID%in%pc[c(2,7,6)]]#6
                          %>%group_by(ARM,V1_disc2)%>%summarise(FGFR2=mean(FGFR2),ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4],shape=2#pch=21,aes(fill=ERBB4)
  )


ggplot( data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9]
                   [!Patient.Study.ID%in%pc[c(2,7)]]#[!Patient.Study.ID%in%c("001-118")]   #6
                   #[!Patient.Study.ID%in%c("001-118","001-119","001-143")]
                   %>%group_by(ARM,V1_disc2)%>%summarise(CDKN2A=mean(CDKN2A),ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4], 
        aes(y=(CDK6) ,      x=   (V1_disc2),col=log(1+CDKN2A-min(CDKN2A)) )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(~ARM,ncol=3)+
  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=4), se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+labs(y="CDK6 ",x="JNK activation score")+
  scale_color_viridis(option="B",begin = 0.1,name="CDKN2A")+
  theme(aspect.ratio=1,
        strip.background = element_blank(),
        strip.text.x = element_blank()
  ) +
  geom_smooth(data=
                data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9][Patient.Study.ID%in%pc[c(2,7,6)]]#6
                           %>%group_by(ARM,V1_disc2)%>%summarise(CDKN2A=mean(CDKN2A),ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4],size=2.5, method="gam",formula= y~s(x,k=3), se=T, col="firebrick2", fill="firebrick2")+
  
  geom_point(data=
               data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9][Patient.Study.ID%in%pc[c(2,7,6)]]#6
                          %>%group_by(ARM,V1_disc2)%>%summarise(CDKN2A=mean(CDKN2A),ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4],shape=2#pch=21,aes(fill=ERBB4)
  )


ggplot( data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9]
                   [!Patient.Study.ID%in%pc[c(2,7)]]#[!Patient.Study.ID%in%c("001-118")]   #6
                   #[!Patient.Study.ID%in%c("001-118","001-119","001-143")]
                   %>%group_by(ARM,V1_disc2)%>%summarise(RORC=mean(RORC),ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4], 
        aes(y=(CDK6) ,      x=   (V1_disc2),col=log(1+RORC-min(RORC)) )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(~ARM,ncol=3)+
  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=4), se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+labs(y="CDK6 ",x="JNK activation score")+
  scale_color_viridis(option="B",begin = 0.1,name="RORC")+
  theme(aspect.ratio=1,
        strip.background = element_blank(),
        strip.text.x = element_blank()
  ) +
  geom_smooth(data=
                data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9][Patient.Study.ID%in%pc[c(2,7,6)]]#6
                           %>%group_by(ARM,V1_disc2)%>%summarise(RORC=mean(RORC),ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4],size=2.5, method="gam",formula= y~s(x,k=3), se=T, col="firebrick2", fill="firebrick2")+
  
  geom_point(data=
               data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9][Patient.Study.ID%in%pc[c(2,7,6)]]#6
                          %>%group_by(ARM,V1_disc2)%>%summarise(RORC=mean(RORC),ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4],shape=2#pch=21,aes(fill=ERBB4)
  )




ggplot( data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9]
                   [!Patient.Study.ID%in%pc[c(2,6,7)]]#[!Patient.Study.ID%in%c("001-118")]
                   #[!Patient.Study.ID%in%c("001-118","001-119","001-143")]
                   %>%group_by(ARM,V1_disc2)%>%summarise(ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4], 
        aes(y=(CDK6) ,      x=   (V1_disc2),col=log(1+ERBB4-min(ERBB4)) )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(~ARM,ncol=3)+
  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=4), se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+labs(y="CDK6 ",x="JNK activation score")+
  scale_color_viridis(option="B",begin = 0.1,name="ERBB4")+
  theme(aspect.ratio=1,
        strip.background = element_blank(),
        strip.text.x = element_blank()
  ) +
  
  geom_point(data=
               data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9][Patient.Study.ID%in%pc[c(2,6,7)]]
                          %>%group_by(ARM,V1_disc2)%>%summarise(ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4],shape=2#pch=21,aes(fill=ERBB4)
  )#+  scale_fill_viridis(option="D",begin = 0.1,name="ERBB4")

ggplot( data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9]
                   [!Patient.Study.ID%in%pc[c(2,6,7)]]#[!Patient.Study.ID%in%c("001-118")]
                   #[!Patient.Study.ID%in%c("001-118","001-119","001-143")]
                   %>%group_by(ARM,V1_disc2)%>%summarise(RORA=mean(RORA),ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4], 
        aes(y=(CDK6) ,      x=   (V1_disc2),col=log(1+RORA-min(RORA)) )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(~ARM,ncol=3)+
  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=4), se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+labs(y="CDK6 ",x="JNK activation score")+
  scale_color_viridis(option="B",begin = 0.1,name="RORA")+
  theme(aspect.ratio=1,
        strip.background = element_blank(),
        strip.text.x = element_blank()
  ) +
  
  geom_point(data=
               data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9][Patient.Study.ID%in%pc[c(2,6,7)]]
                          %>%group_by(ARM,V1_disc2)%>%summarise(RORA=mean(RORA),ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4],shape=2#pch=21,aes(fill=ERBB4)
  )#+  scale_fill_viridis(option="D",b


ggplot( data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9]
                   [Patient.Study.ID%in%pc[c(2,6,11,7)]]#[!Patient.Study.ID%in%c("001-118")]
                   #[!Patient.Study.ID%in%c("001-118","001-119","001-143")]
                   %>%group_by(ARM,V1_disc2)%>%summarise(ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4], 
        aes(y=(CDK6) ,      x=   (V1_disc2),col=log(ESR1-min(ESR1)) )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(~ARM,ncol=3)+
  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=4), se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+labs(y="CDK6 ",x="JNK activation score")+
  scale_color_viridis(option="B",begin = 0.1,name="ESR1")+
  theme(aspect.ratio=1,
        strip.background = element_blank(),
        strip.text.x = element_blank()
  ) 
pc[c(2,6,11,7)]
ggplot( data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9]
                   [!Patient.Study.ID%in%c("001-118")]
                   #[!Patient.Study.ID%in%c("001-118","001-119","001-143")]
                   %>%group_by(Day,ARM,Patient.Study.ID,V1_disc2)%>%summarise(ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4], 
        aes(x=(CDK6) ,group=Patient.Study.ID,      col=   log(ESR1-min(ESR1)),y=log(ERBB4-min(ERBB4)) )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(Day~ARM,ncol=3)+
 # geom_smooth(size=2.5, method="gam",formula= y~s(x,k=4), se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+
  scale_color_viridis(option="B",begin = 0.1)
 
gmdd<- data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9][!Patient.Study.ID%in%c("001-118")]%>%group_by(Day,ARM,Patient.Study.ID,V1_disc2)%>%summarise(ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4]
gmdd[,lnERBB4:=log(1+ERBB4-min(ERBB4))]
gmdd[,lnESR1:=log(1+ESR1-min(ESR1))]
gmB<- gam(CDK6~ te(lnERBB4 ,V1_disc2,k=c(3,3)),       data=gmdd[ARM=="B"])
vis.gam( gmB,theta=150 )
gmC<- gam(CDK6~ te(lnERBB4 ,V1_disc2,k=c(3,3)),       data=gmdd[ARM=="C"])
vis.gam( gmC ,theta=150)
gmA<- gam(CDK6~ te(lnERBB4 ,V1_disc2,k=c(4,4)),       data=gmdd[ARM=="A"])
vis.gam( gmA ,theta=50)

gmB<- gam(CDK6~ s(lnERBB4,k=4) + s(V1_disc2,k=4),       data=gmdd[ARM=="B"])
plot( gmB )
gmC<- gam(CDK6~ s(lnERBB4,k=3) + s(V1_disc2,k=3),       data=gmdd[ARM=="C"][!Patient.Study.ID%in%c("001-118")])
plot( gmC )

gmC<- gam(CDK6~ s(lnERBB4,k=3) + s(V1_disc2,k=3),       data=gmdd[ARM=="A"])
plot( gmC )
11
6
7
pc[c(2,6,11,7)]
gmC<- gam(lnESR1~  s(V1_disc2,k=6,by=as.factor(Patient.Study.ID)),       data=gmdd[ARM=="C"][Patient.Study.ID%in%pc[7]])
plot( gmC )
fin1<-ggplot( data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9]%>%group_by(ARM,V1_disc2)%>%summarise(FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4], 
        aes(y=(CDK6) ,      x=   (V1_disc2),col=log(ESR1-min(ESR1)) )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(~ARM,ncol=3)+
  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=4), se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+labs(y="CDK6 ",x="JNK activation score")+
  scale_color_viridis(option="B",begin = 0.1,name="ESR1")+
  theme(aspect.ratio=1,
        strip.background = element_blank(),
        strip.text.x = element_blank()
  ) +
  theme(axis.title=element_blank(), axis.text=element_blank(),legend.title=element_blank(),legend.text=element_blank()) +theme(legend.position = "none")

ggsave(fin1, filename = "/Users/jason/Dropbox/FELINE Project/Manuscript/Final Main Figures/Blank figures/ESR1 overlay of CDK vs JNK.png",height=3,width=9)


fin2<- ggplot( data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9]%>%group_by(ARM,V1_disc2)%>%summarise(ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=median(ESR1),CDK6=mean(CDK6)))[CDK6<2.4], 
        aes(y=(CDK6) ,      x=   (V1_disc2),col=log(FGFR2-min(FGFR2)) )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(~ARM,ncol=3)+
  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=4), se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+labs(y="CDK6 ",x="JNK activation score")+
  scale_color_viridis(option="B",begin = 0.1,name="FGFR2")+
  theme(aspect.ratio=1,
        strip.background = element_blank(),
        strip.text.x = element_blank()
  ) +
  theme(axis.title=element_blank(), axis.text=element_blank(),legend.title=element_blank(),legend.text=element_blank()) +theme(legend.position = "none")

ggsave(fin2, filename = "/Users/jason/Dropbox/FELINE Project/Manuscript/Final Main Figures/Blank figures/FGFR2 overlay of CDK vs JNK.png",height=3,width=9)

fin3<- ggplot( data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9]%>%group_by(ARM,V1_disc2)%>%summarise(ERBB4=median(ERBB4),FGFR2=mean(FGFR2),ESR1=median(ESR1),CDK6=mean(CDK6)))[CDK6<2.4], 
        aes(y=(CDK6) ,      x=   (V1_disc2),col=log(ERBB4-min(ERBB4)) )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(~ARM,ncol=3)+
  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=4), se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+labs(y="CDK6 ",x="JNK activation score")+
  scale_color_viridis(option="B",begin = 0.1,name="ERBB4")+
  theme(aspect.ratio=1,
        strip.background = element_blank(),
        strip.text.x = element_blank()
  )+ 
  theme(axis.title=element_blank(), axis.text=element_blank(),legend.title=element_blank(),legend.text=element_blank()) +theme(legend.position = "none")

ggsave(fin3, filename = "/Users/jason/Dropbox/FELINE Project/Manuscript/Final Main Figures/Blank figures/ERBB4 overlay of CDK vs JNK.png",height=3,width=9)








ggplot( mrg_umap_mapk_tsp2[Day==180][(V1_disc2)> -1.8&(V1_disc2)< 2]%>%group_by(ARM,V1_disc2)%>%summarise(ERBB4=median(ERBB4),FGFR2=mean(FGFR2),ESR1=median(ESR1),CDK6=mean(CDK6)), 
        aes(y=(CDK6) ,      x=   (V1_disc2),col=log(ERBB4-min(ERBB4)) )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(~ARM,ncol=3)+
  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=4), se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+labs(y="CDK6 ",x="JNK activation score")+
  scale_color_viridis(option="B",begin = 0.1,name="ERBB4")+
  theme(aspect.ratio=1,
        strip.background = element_blank(),
        strip.text.x = element_blank()
  )


ggplot( mrg_umap_mapk_tsp2[Day==180][(V1_disc2)> -1.8&(V1_disc2)< 2]%>%group_by(ARM,V1_disc2)%>%summarise(ERBB4=median(ERBB4),FGFR2=mean(FGFR2),ESR1=median(ESR1),CDK6=mean(CDK6)), 
        aes(y=(CDK6) ,      x=   (V1_disc2),col=log(ESR1-min(ESR1)) )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(~ARM,ncol=3)+
  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=4), se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+labs(y="CDK6 ",x="JNK activation score")+
  scale_color_viridis(option="B",begin = 0.1,name="ESR1")+
  theme(aspect.ratio=1,
        strip.background = element_blank(),
        strip.text.x = element_blank()
  )



ggplot( mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9]%>%group_by(ARM,V1_disc2)%>%summarise(ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=median(ESR1),CDK6=mean(CDK6)), 
        aes(y=log(ERBB4-min(ERBB4)) ,      x= log(FGFR2-min(FGFR2)),col=  (CDK6) )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(~ARM,ncol=3)+
  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=4), se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+#labs(y="CDK6 ",x="JNK activation score")+
  scale_color_viridis(option="B",begin = 0.1)+
  theme(aspect.ratio=1,
        strip.background = element_blank(),
        strip.text.x = element_blank()
  )

# 
# ggplot( data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9]
#                    [!Patient.Study.ID%in%pc[c(2,6,7)]]#[!Patient.Study.ID%in%c("001-118")]
#                    #[!Patient.Study.ID%in%c("001-118","001-119","001-143")]
#                    %>%group_by(ARM,V1_disc2)%>%summarise(ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4], 
#         aes(y=(CDK6) ,      x=   (V1_disc2),col=log(1+ESR1-min(ESR1)) )) +
#   #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
#   facet_wrap(~ARM,ncol=3)+
#   geom_smooth(size=2.5, method="gam",formula= y~s(x,k=4), se=T, col="firebrick", fill="firebrick")+
#   geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+labs(y="CDK6 ",x="JNK activation score")+
#   scale_color_viridis(option="B",begin = 0.1,name="ESR1")+
#   theme(aspect.ratio=1,
#         strip.background = element_blank(),
#         strip.text.x = element_blank()
#   ) +
#   
#   geom_point(data=
#   data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9][Patient.Study.ID%in%pc[c(2,6,7)]]
#              %>%group_by(ARM,V1_disc2)%>%summarise(ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4],shape=2#pch=21,aes(fill=ERBB4)
#    )#+  scale_fill_viridis(option="D",begin = 0.1,name="ERBB4")
# 
# ggplot( data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9]
#                    [!Patient.Study.ID%in%pc[c(2,6,7)]]#[!Patient.Study.ID%in%c("001-118")]
#                    #[!Patient.Study.ID%in%c("001-118","001-119","001-143")]
#                    %>%group_by(ARM,V1_disc2)%>%summarise(ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4], 
#         aes(y=(CDK6) ,      x=   (V1_disc2),col=log(1+FGFR2-min(FGFR2)) )) +
#   #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
#   facet_wrap(~ARM,ncol=3)+
#   geom_smooth(size=2.5, method="gam",formula= y~s(x,k=4), se=T, col="firebrick", fill="firebrick")+
#   geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+labs(y="CDK6 ",x="JNK activation score")+
#   scale_color_viridis(option="B",begin = 0.1,name="ESR1")+
#   theme(aspect.ratio=1,
#         strip.background = element_blank(),
#         strip.text.x = element_blank()
#   ) +
#   
#   geom_point(data=
#                data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9][Patient.Study.ID%in%pc[c(2,6,7)]]
#                           %>%group_by(ARM,V1_disc2)%>%summarise(ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4],shape=2#pch=21,aes(fill=ERBB4)
#   )#+  scal

ggplot( data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9]
                   [!Patient.Study.ID%in%pc[c(2,7)]]#[!Patient.Study.ID%in%c("001-118")]   #6
                   #[!Patient.Study.ID%in%c("001-118","001-119","001-143")]
                   %>%group_by(ARM,V1_disc2)%>%summarise(ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4], 
        aes(y=(CDK6) ,      x=   (V1_disc2),col=log(1+ERBB4-min(ERBB4)) )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(~ARM,ncol=3)+
  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=4), se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+labs(y="CDK6 ",x="JNK activation score")+
  scale_color_viridis(option="B",begin = 0.1,name="ERBB4")+
  theme(aspect.ratio=1,
        strip.background = element_blank(),
        strip.text.x = element_blank()
  ) +
  geom_smooth(data=
                data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9][Patient.Study.ID%in%pc[c(2,7,6)]]#6
                           %>%group_by(ARM,V1_disc2)%>%summarise(ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4],size=2.5, method="gam",formula= y~s(x,k=3), se=T, col="firebrick2", fill="firebrick2")+
  
  geom_point(data=
               data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9][Patient.Study.ID%in%pc[c(2,7,6)]]#6
                          %>%group_by(ARM,V1_disc2)%>%summarise(ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4],shape=2#pch=21,aes(fill=ERBB4)
  )

ggplot( data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9][ARM=="C"]
                   [!Patient.Study.ID%in%pc[c(2,7)]]#[!Patient.Study.ID%in%c("001-118")]   #6
                   #[!Patient.Study.ID%in%c("001-118","001-119","001-143")]
                   %>%group_by(ARM,V1_disc2,Patient.Study.ID)%>%summarise(ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4], 
        aes( group=Patient.Study.ID,    x=   (V1_disc2),y=log(1+ESR1-min(ESR1)) )) +
  #geom_point(data=zz[prop_tour<0.5]%>%group_by(Response,Day,ARM,tour)%>%dplyr::summarise(prop_tour=median(prop_tour)) )
  facet_wrap(~Patient.Study.ID,ncol=3)+
  geom_smooth(size=2.5, method="gam",formula= y~s(x,k=4), se=T, col="firebrick", fill="firebrick")+
  geom_point(size=2)+  theme_classic(base_size=18)+theme(aspect.ratio=1)+labs(y="CDK6 ",x="JNK activation score")+
  scale_color_viridis(option="B",begin = 0.1,name="ESR1")+
  theme(aspect.ratio=1,
        strip.background = element_blank(),
        strip.text.x = element_blank()
  ) +
  geom_smooth(data=
                data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9][Patient.Study.ID%in%pc[c(2,7,6)]]#6
                           %>%group_by(ARM,V1_disc2)%>%summarise(ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4],size=2.5, method="gam",formula= y~s(x,k=3), se=T, col="firebrick2", fill="firebrick2")+
  
  geom_point(data=
               data.table(mrg_umap_mapk_tsp2[Day==180][abs(V1_disc2)<1.9][Patient.Study.ID%in%pc[c(2,7,6)]]#6
                          %>%group_by(ARM,V1_disc2)%>%summarise(ERBB4=mean(ERBB4),FGFR2=mean(FGFR2),ESR1=mean(ESR1),CDK6=mean(CDK6)))[CDK6<2.4],shape=2#pch=21,aes(fill=ERBB4)
  )

