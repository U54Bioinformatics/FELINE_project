rm(list=ls())
require(mgcv);require(data.table);require(dplyr);require(ggplot2);require(tidyr);require(parallel);require(lme4);require(lmerTest);require(parallel)
library(effects)#require(umap)
library(TSP)
require(mclust);require(umap)
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

signal_recept_tranUnSlice <- na.omit(data.table( ttt3[Day==180][VariableType=="gene"][Major_catagory%in%c("Signal transduction")]%>%
                                                   dplyr::select(-c(VariableType,Negative.regulator,Sub_catagory,Expression2,sign_nn,Major_catagory))%>%spread(Variable,Expression)
))
signal_recept_tranSlice <- na.omit(data.table(signal_recept_tranUnSlice %>% group_by(Patient.Study.ID,ARM, Day,Subclone,Response) %>% dplyr::slice(1:100)))

in_dd<- signal_recept_tranSlice%>%dplyr::select(-c(Cell.ID, Patient.Study.ID, ARM, Day,  Response, Subclone, Phase))
umap_md3 <- umap(in_dd, n_neighbors=25, n_components=5)
#umap_md3 <- umap(signal_recept_tranSlice%>%dplyr::select(-c(Cell.ID, Patient.Study.ID, ARM, Day,  Response, Subclone, Phase)), n_neighbors=10, n_components=2)

#ggplot(   data.table(signal_recept_tranSlice,umap_md3$layout) , aes(V1,V2,col=Subclone)) + geom_point(size=2)+facet_wrap(~Patient.Study.ID)
ggplot(   data.table(signal_recept_tranSlice,umap_md3$layout) , aes(V1,V2,col=Response,shape=Subclone))+geom_point(size=2)+theme_classic()+facet_wrap(~ARM)
ggplot(   data.table(signal_recept_tranSlice,umap_md3$layout) , aes(V1,V2,col=Patient.Study.ID,shape=Subclone))+geom_point(size=2)+theme_classic()
ggplot(   data.table(signal_recept_tranSlice,umap_md3$layout) , aes(V1,V2,col=Day,shape=Subclone))+geom_point(size=2)+theme_classic()

corV1<- as.data.table(cor(signal_recept_tranSlice%>%dplyr::select(-c(Cell.ID, Patient.Study.ID, ARM, Day,  Response, Subclone, Phase)),
                          umap_md3$layout[,]),keep.rownames=TRUE)
corV1[order(-abs(V1))]
corV1[order(-abs(V3))][1:20]

u_datFULL <- data.table(signal_recept_tranUnSlice,
                        predict(umap_md3,
                                data=signal_recept_tranUnSlice %>%
                                  dplyr::select(-c(Cell.ID, Patient.Study.ID, ARM, Day,  Response, Subclone, Phase))   )   )

corV1 <- as.data.table(cor(u_datFULL%>%dplyr::select(-c(Cell.ID, Patient.Study.ID, ARM, Day,  Response, Subclone, Phase,V1:V5)),
                           u_datFULL%>%dplyr::select(V1:V5) ),keep.rownames=TRUE)
corV1[order(-abs(V1))][1:5]
corV1[order(-abs(V5))][1:5]
ggplot(   u_datFULL, aes(V1,V5,col=Response,shape=Subclone))+geom_point(size=2)+theme_classic()+facet_wrap(~ARM)
ggplot(  u_datFULL , aes(V1,V2,col=MAPK10,shape=Subclone))+geom_point(size=2)+facet_wrap(Response~Day)+theme_classic()+scale_color_viridis()
ggplot(  u_datFULL , aes(V1,V2,col=MAP3K14,shape=Subclone))+geom_point(size=2)+facet_wrap(Response~Day)+theme_classic()+scale_color_viridis()
ggplot(  u_datFULL , aes(V1,V2,col=FOS,shape=Subclone))+geom_point(size=2)+facet_wrap(Response~Day)+theme_classic()+scale_color_viridis()
ggplot(  u_datFULL , aes(MAPK10,FOS,col=Response,shape=Subclone))+geom_point(size=2)+facet_wrap(Response~Day)+theme_classic()+scale_color_viridis()

summary(lm(V1~V2,data=data.table(signal_recept_tranSlice,umap_md3$layout)))



signal_recept_tran_all <- ttt3[ARM%in%c("A","B","C")][Day==180 ][VariableType=="gene"][Major_catagory%in%c("Signal transduction")]%>%dplyr::select(-c(VariableType,Negative.regulator,Sub_catagory,Expression2,sign_nn,Major_catagory))%>%spread(Variable,Expression)


lu <- expand.grid(wchArm=c("A","B","C"),wchRes=c("Non-responder","Responder"))




signal_recept_tranUnSliceAllt <- na.omit(data.table( ttt3[ARM%in%c("A","B","C")][VariableType=="gene"][Major_catagory%in%c("Signal transduction")]%>%
                                                       dplyr::select(-c(VariableType,Negative.regulator,Sub_catagory,Expression2,sign_nn,Major_catagory))%>%spread(Variable,Expression)
))

project_this <- signal_recept_tranUnSliceAllt%>%dplyr::select(-c(Cell.ID, Patient.Study.ID, ARM, Day,  Response, Subclone, Phase))

#rm(list=c("ttt","ttt2"))

u_datFULL_ALLtime <- data.table(signal_recept_tranUnSliceAllt, predict(umap_md3, data=project_this  )   )
pcaOverall<- princomp(u_datFULL_ALLtime%>%dplyr::select(V1,V5))
u_datFULL_ALLtimeDims<-data.table(u_datFULL_ALLtime, dimension1=pcaOverall$scores[,1],dimension2=pcaOverall$scores[,2])
#u_datFULL_ALLtime[1:1000]%>%dplyr::select(V1:V5)%>%pairs

corV1<- as.data.table(cor(u_datFULL%>%dplyr::select(-c(Cell.ID, Patient.Study.ID, ARM, Day,  Response, Subclone, Phase,V1:V5)),
                          u_datFULL%>%dplyr::select(V1:V5) ),keep.rownames=TRUE)
corV1[order(-abs(V1))][1:10]
corV1[order(-abs(V2))][1:10]
corV1[order(-abs(V3))][1:10]
corV1[order(-abs(V5))][1:10]

ggplot(   u_datFULL_ALLtime, aes(V1, V2, col=Response))+geom_point(size=0.85,alpha=0.75)+theme_classic()+facet_wrap(ARM~Day)
ggplot(   u_datFULL_ALLtime, aes(V1, V2, col=Response))+geom_point(size=0.85,alpha=0.75)+theme_classic(base_size=22)+facet_wrap(~Day)

ggplot(   u_datFULL_ALLtimeDims[ARM!="A"], aes(dimension1, dimension2, col=Response))+geom_point(size=0.85,alpha=0.75)+theme_classic(base_size=22)+facet_wrap(ARM~Day)
ggplot(   u_datFULL_ALLtimeDims, aes(V1, -V2))+geom_hex(bins = 100)+  scale_fill_viridis(name="Cell count")+theme_classic(base_size=22)+labs(y="Umap dimension 2", x="Umap dimension 1") +theme(aspect.ratio = 1)+
  facet_wrap(~ARM)

ggplot(   u_datFULL_ALLtimeDims, aes(V1, V5, col=Response))+geom_point(size=0.85,alpha=0.75)+theme_classic(base_size=22)+facet_wrap(ARM~Day)
ggplot(   u_datFULL_ALLtimeDims, aes(V2, V5, col=Response))+geom_point(size=0.85,alpha=0.75)+theme_classic(base_size=22)+facet_wrap(ARM~Day)
ggplot(   u_datFULL_ALLtimeDims, aes(V1, V2, col=Response))+geom_point(size=0.5,alpha=0.75)+theme_classic(base_size=22)+facet_wrap(ARM~Day)
summary(lm(V1~V2,data=u_datFULL_ALLtimeDims))
pcaOverall<- princomp(u_datFULL_ALLtimeDims%>%dplyr::select(AKT1:MAPKBP1))
plot( cumsum(pcaOverall$sdev^2/ sum(pcaOverall$sdev^2)),type="b")
#save(u_datFULL_ALLtimeDims, umap_md3,file="/Users/jason/Dropbox/FELINE Project/Data_analysis/scRNA/15_subclone_interaction/Dim reduced transduction only 3 bproj.RData")
load(file="/Users/jason/Dropbox/FELINE Project/Data_analysis/scRNA/15_subclone_interaction/Dim reduced transduction only 3 bproj.RData")


corV3<- as.data.table(cor(u_datFULL_ALLtimeDims%>%dplyr::select(-c(Cell.ID, Patient.Study.ID, ARM, Day,  Response, Subclone, Phase,V1:dimension2)),
                          u_datFULL_ALLtimeDims%>%dplyr::select(V1:dimension2) ),keep.rownames=TRUE)
corV3[order(-abs(dimension1))][1:10]
corV3[order(-abs(dimension2))][1:10]
corV3[order(-abs(V1))][1:10]
corV3[order(-abs(V2))][1:10]
corV3[order(-abs(V5))][1:10]

corVA<- as.data.table(cor(u_datFULL_ALLtimeDims[ARM=="A"]%>%dplyr::select(V1:dimension2) ,
                          1*(u_datFULL_ALLtimeDims[ARM=="A"]$Response=="Responder") ),
                      keep.rownames=TRUE)
corVB<- as.data.table(cor(u_datFULL_ALLtimeDims[ARM=="B"]%>%dplyr::select(V1:dimension2) ,
                          1*(u_datFULL_ALLtimeDims[ARM=="B"]$Response=="Responder") ),
                      keep.rownames=TRUE)
corVC<- as.data.table(cor(u_datFULL_ALLtimeDims[ARM=="C"]%>%dplyr::select(V1:dimension2) ,
                          1*(u_datFULL_ALLtimeDims[ARM=="C"]$Response=="Responder") ),
                      keep.rownames=TRUE)

pcaOverall<- princomp(u_datFULL_ALLtime%>%dplyr::select(V1,V2))
u_datFULL_ALLtimeDims<-data.table(u_datFULL_ALLtime, dimension1=pcaOverall$scores[,1],dimension2=pcaOverall$scores[,2])
u_datFULL_ALLtimeDims[,V1_disc:=(round(V1*10)/10) ]
u_datFULL_ALLtimeDims[,V2_disc:=(round(V2*3)/3) ]
MAPK_Major_axis<-data.table(u_datFULL_ALLtimeDims%>%dplyr::select(Cell.ID,Patient.Study.ID,ARM,Day,Response,Subclone,V1,V2,V1_disc,    V2_disc))


ggplot(   u_datFULL_ALLtimeDims[abs(V1_disc)<=2], 
          aes(y=V2_disc, x=V1_disc))+#
  geom_point()
ggplot(   u_datFULL_ALLtimeDims[abs(V1_disc)<=2]%>%group_by(V1_disc,V2_disc)%>%dplyr::summarise(freq=length(Cell.ID)), 
          aes(y=V2_disc, x=V1_disc,fill=freq))+#
  geom_tile(size=1.1,alpha=1)



extra_genes<-spread( ttt3[Variable%in%c("ESR1","ERBB4","FGFR2")]%>%dplyr::select(Cell.ID,Variable,Expression),Variable,Expression)
which_plot<-unique( c(corV3[order(-abs(V1))][1:10]$rn,  corV3[order(-abs(V2))][1:10]$rn) )
which_plot<-which_plot[which_plot!="MAPK1" ]


subset_u_datFULL_ALLtimeDims <- data.table(u_datFULL_ALLtimeDims%>% #[ARM!="A"]
                                             dplyr::select(c("Patient.Study.ID", "ARM", "Day",  "Response", "V1_disc"),which_plot)%>%
                                             gather(var,val,-c(Patient.Study.ID, ARM, Day,  Response,V1_disc)) )


summar_u_datFULL_ALLtimeDims <-data.table(subset_u_datFULL_ALLtimeDims   %>% group_by(Patient.Study.ID, ARM, Day,  Response,V1_disc,var)%>%
                                            dplyr::summarise(val_mu=mean(val),ncell=length(val)))

subset_u_datFULL_ALLtimeDims[var=="FOS"][Day==180][Response=="Responder"][V1_disc==median(V1_disc)]
subset_u_datFULL_ALLtimeDims[var=="FOS"][Day==180][Response=="Non-responder"][V1_disc==median(V1_disc)]

#save(MAPK_Major_axis,which_plot,corV3,subset_u_datFULL_ALLtimeDims,file="/Users/jason/Dropbox/FELINE Project/Data_analysis/scRNA/15_subclone_interaction/MAPK major axis proj.RData")

summar_u_datFULL_ALLtimeDims2 <- data.table( subset_u_datFULL_ALLtimeDims%>% group_by( Response,ARM, Day,V1_disc)%>%
                                               summarise(ncell=length(val))%>%spread(Response,ncell)%>%mutate(prop_resp=Responder/(Responder+`Non-responder` )))
summar_u_datFULL_ALLtimeDims2[is.na(Responder),Responder:=0]
summar_u_datFULL_ALLtimeDims2[is.na(`Non-responder`),`Non-responder`:=0]
summar_u_datFULL_ALLtimeDims2[,prop_resp:=Responder/(Responder+`Non-responder` )]
#   group_by( Response,ARM, Day,V1_disc)
#summar_u_datFULL_ALLtimeDims3 <-  data.table(summar_u_datFULL_ALLtimeDims           %>%
#               group_by( ARM, Day,V1_disc)%>%  
#               mutate(nsamp=sum(ncell) )%>%
#               group_by( Response,ARM, Day,V1_disc)%>%
#               mutate(prop_resp=ncell/nsamp)  
#             )

summar_u_datFULL_ALLtimeDims4<- merge(summar_u_datFULL_ALLtimeDims2,
                                      u_datFULL_ALLtimeDims%>%group_by(ARM,Day)%>%summarise(s=sum(Response=="Responder")/sum(Response%in%c("Responder","Non-responder"))),
                                      by=c("ARM",     "Day"))
summar_u_datFULL_ALLtimeDims4[,prop_resp2:=prop_resp/(s*2),by=c("ARM","Day")]
summar_u_datFULL_ALLtimeDims4[,prop_resp2:=prop_resp]
summar_u_datFULL_ALLtimeDims4[ARM=="A",prop_resp2:=0.25*prop_resp/s]

draw<- data.table(summar_u_datFULL_ALLtimeDims4[Day==180]%>%group_by(ARM)%>%mutate(sum_prop=sum(prop_resp),mean_prop=mean(prop_resp)))#[(`Non-responder` +Responder) >1400]
draw%>%group_by(ARM)%>%mutate(sum_prop=sum(prop_resp),mean_prop=mean(prop_resp))
ggplot(   draw[(`Non-responder` +Responder) >200], aes(x=V1_disc,y=(prop_resp2),col=ARM,fill=ARM))+#
  geom_point(size=2,alpha=1)+
  facet_wrap(~ARM,scales="free_y")+
  theme_classic(base_size=26)+theme(aspect.ratio=1)+ scale_color_brewer(palette = "Dark2")+ scale_fill_brewer(palette = "Dark2") +
  #+ scale_color_grey(start = 0.8, end = 0.2) + scale_fill_grey(start = 0.8, end = 0.2) +
  geom_smooth(method="gam",formula=y~s(x,k=5), method.args=list(family="quasibinomial", gamma=6)) +labs(y="Proportion of cells from responders", x="MAPK activation phenotpye" )+
  theme(legend.position="none")+
  theme(legend.position="none")+theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

ggplot(   draw[(`Non-responder` +Responder) >200], aes(x=V1_disc,y=(prop_resp2),col=ARM,fill=ARM))+#
  geom_point(size=2,alpha=1)+
  facet_wrap(~ARM)+
  theme_classic(base_size=26)+theme(aspect.ratio=1)+ scale_color_brewer(palette = "Dark2")+ scale_fill_brewer(palette = "Dark2") +
  #+ scale_color_grey(start = 0.8, end = 0.2) + scale_fill_grey(start = 0.8, end = 0.2) +
  geom_smooth(method="gam",formula=y~s(x,k=5), method.args=list(family="quasibinomial", gamma=6)) +labs(y="Proportion of cells from responders", x="MAPK activation phenotpye" )+
  theme(legend.position="none")+
  theme(legend.position="none")+theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )
ggplot(   summar_u_datFULL_ALLtimeDims4[Day==180][(`Non-responder` +Responder) >1400], aes(x=V1_disc,y=(prop_resp),col=ARM,fill=ARM))+#
  geom_point(size=2,alpha=1)+
  facet_wrap(~ARM,scales="free_y")+
  theme_classic(base_size=26)+theme(aspect.ratio=1)+ scale_color_brewer(palette = "Dark2")+ scale_fill_brewer(palette = "Dark2") +
  #+ scale_color_grey(start = 0.8, end = 0.2) + scale_fill_grey(start = 0.8, end = 0.2) +
  geom_smooth(method="gam",formula=y~s(x,k=4), method.args=list(gamma=1.55,family="quasibinomial")) +labs(y="Proportion of cells from responders", x="JNK activation phenotpye" )+
  theme(legend.position="none")+theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )


ggplot(   summar_u_datFULL_ALLtimeDims4[Day==180], aes(x=V1_disc,y=(prop_resp),col=ARM,fill=ARM))+#
  geom_point(size=2,alpha=1)+
  facet_wrap(~ARM)+
  theme_classic(base_size=26)+theme(aspect.ratio=1)+ scale_color_brewer(palette = "Dark2")+ scale_fill_brewer(palette = "Dark2") +
  #+ scale_color_grey(start = 0.8, end = 0.2) + scale_fill_grey(start = 0.8, end = 0.2) +
  geom_smooth(method="gam",formula=y~s(x,k=4), method.args=list(gamma=1.55,family="quasibinomial")) +labs(y="Proportion of cells from responders", x="JNK activation phenotpye" )+
  theme(legend.position="none")+theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+xlim(c(-2,2))


summar_u_datFULL_ALLtimeDims3 <- data.table(merge(summar_u_datFULL_ALLtimeDims,summar_u_datFULL_ALLtimeDims2%>%dplyr::select(ARM ,Day  ,  V1_disc,prop_resp),by=c("ARM","Day","V1_disc")))
#save(summar_u_datFULL_ALLtimeDims3,summar_u_datFULL_ALLtimeDims2 ,file="/Users/jason/Dropbox/FELINE Project/Data_analysis/scRNA/15_subclone_interaction/Visualising MAPK transduction major axis.RData")
#load(file="/Users/jason/Dropbox/FELINE Project/Data_analysis/scRNA/15_subclone_interaction/Visualising MAPK transduction major axis.RData")
#summar_u_datFULL_ALLtimeDims3,summar_u_datFULL_ALLtimeDims2
ggplot(   summar_u_datFULL_ALLtimeDims2[ARM!="A"][Day==180], aes(x=V1_disc,y=(prop_resp),col=ARM,fill=ARM))+#
  geom_point(size=2,alpha=1)+
  facet_grid(~ARM)+
  theme_classic(base_size=22)+theme(aspect.ratio=1)+ scale_color_brewer(palette = "Dark2")+ scale_fill_brewer(palette = "Dark2") +
  #+ scale_color_grey(start = 0.8, end = 0.2) + scale_fill_grey(start = 0.8, end = 0.2) +
  geom_smooth(method="gam", method.args=list(family="quasibinomial")) +labs(y="Proportion of cells from responders", x="MAPK activation phenotpye" )+
  theme(legend.position="none")

ggplot(   summar_u_datFULL_ALLtimeDims3[var!="V2_disc"][ARM!="A"][Day==180], aes(y=val_mu, x=V1_disc, group=Response,col=Response))+#
  geom_point(size=1.1,alpha=1)+
  facet_grid(ARM~var)+
  theme_classic(base_size=22)+
  geom_smooth(method="gam")+theme(aspect.ratio=1)+
  labs(y="MAPK gene expression", x="MAPK activation phenotpye")+theme(aspect.ratio=1)

ggplot(   summar_u_datFULL_ALLtimeDims3[var!="V2_disc"][Day==180], 
          aes(y=val_mu, x=V1_disc, group=Response,col=Response))+#
  geom_point(size=1.1,alpha=1)+
  facet_wrap(~var,scales="free")+
  theme_classic(base_size=22)+
  geom_smooth(method="gam")+theme(aspect.ratio=1)+
  labs(y="MAPK gene expression", x="JNK activation phenotpye")+theme(aspect.ratio=1)


pl_sub<-summar_u_datFULL_ALLtimeDims3[ARM!="A"][Day==180][var%in%c("MAP3K1","MAP3K14","MAP3K5","MAP4K3","MAPK10","MAPK8")]
pl_sub[,slope:="neg"]
pl_sub[var%in%c("MAP3K5","MAPK10","MAPK8"),slope:="pos"]
pl_sub$var <- factor(pl_sub$var , levels = c("MAP3K5","MAPK10","MAPK8",
                                             "MAP3K1","MAP3K14","MAP4K3"
))

ggplot(   pl_sub[!(var=="MAPK3"&val_mu>0.5)&!(var=="MAPK3"&Patient.Study.ID=="001-002")], aes(y=val_mu, x=V1_disc, group=var,col=var))+#
  geom_point(size=1.1,alpha=1)+
  facet_wrap(~var,nrow=2,scales="free_y")+
  theme_classic(base_size=22)+
  geom_smooth(method="gam",formula=y~s(x,k=4))+theme(aspect.ratio=1)+
  labs(y="MAPK gene expression", x="MAPK activation phenotpye")+ theme(legend.position="none")+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )


ggplot(   pl_sub[var=="MAPK3"][!(var=="MAPK3"&val_mu>0.5)], aes(y=val_mu, x=V1_disc, group=var,col=Patient.Study.ID))+#
  geom_point(size=1.1,alpha=1)+
  facet_wrap(~Patient.Study.ID,nrow=2,scales="free_y")+
  theme_classic(base_size=22)+
  geom_smooth(method="gam",formula=y~s(x,k=4))+theme(aspect.ratio=1)+
  labs(y="MAPK gene expression", x="MAPK activation phenotpye")+ theme(legend.position="none")



#library(factoextra)
subset_u_datFULL_ALLtimeDims <- data.table(merge(u_datFULL_ALLtimeDims,extra_genes,by="Cell.ID")%>% #[ARM!="A"]
                                             dplyr::select(c("Patient.Study.ID", "ARM", "Day",  "Response", "V1_disc"),which_plot,"ESR1","ERBB4","FGFR2","JUNB")%>%
                                             gather(var,val,-c(Patient.Study.ID, ARM, Day,  Response,V1_disc)) )


summar_u_datFULL_ALLtimeDims <-data.table(subset_u_datFULL_ALLtimeDims   %>% group_by(Patient.Study.ID, ARM, Day,  Response,V1_disc,var)%>%
                                            dplyr::summarise(val_mu=mean(val),ncell=length(val)))

subset_u_datFULL_ALLtimeDims[var=="FOS"][Day==180][Response=="Responder"][V1_disc==median(V1_disc)]
subset_u_datFULL_ALLtimeDims[var=="FOS"][Day==180][Response=="Non-responder"][V1_disc==median(V1_disc)]

#save(MAPK_Major_axis,which_plot,corV3,subset_u_datFULL_ALLtimeDims,file="/Users/jason/Dropbox/FELINE Project/Data_analysis/scRNA/15_subclone_interaction/MAPK major axis proj.RData")

summar_u_datFULL_ALLtimeDims2 <- data.table( subset_u_datFULL_ALLtimeDims%>% group_by( Response,ARM, Day,V1_disc)%>%
                                               summarise(ncell=length(val))%>%spread(Response,ncell)%>%mutate(prop_resp=Responder/(Responder+`Non-responder` )))
summar_u_datFULL_ALLtimeDims2[is.na(Responder),Responder:=0]
summar_u_datFULL_ALLtimeDims2[is.na(`Non-responder`),`Non-responder`:=0]
summar_u_datFULL_ALLtimeDims2[,prop_resp:=Responder/(Responder+`Non-responder` )]
#   group_by( Response,ARM, Day,V1_disc)
#summar_u_datFULL_ALLtimeDims3 <-  data.table(summar_u_datFULL_ALLtimeDims           %>%
#               group_by( ARM, Day,V1_disc)%>%  
#               mutate(nsamp=sum(ncell) )%>%
#               group_by( Response,ARM, Day,V1_disc)%>%
#               mutate(prop_resp=ncell/nsamp)  
#             )

summar_u_datFULL_ALLtimeDims4<- merge(summar_u_datFULL_ALLtimeDims2,
                                      u_datFULL_ALLtimeDims%>%group_by(ARM,Day)%>%summarise(s=sum(Response=="Responder")/sum(Response%in%c("Responder","Non-responder"))),
                                      by=c("ARM",     "Day"))
summar_u_datFULL_ALLtimeDims4[,prop_resp2:=prop_resp/(s*2),by=c("ARM","Day")]
summar_u_datFULL_ALLtimeDims4[,prop_resp2:=prop_resp]
summar_u_datFULL_ALLtimeDims4[ARM=="A",prop_resp2:=0.25*prop_resp/s]


summar_u_datFULL_ALLtimeDims3 <- data.table(merge(summar_u_datFULL_ALLtimeDims,summar_u_datFULL_ALLtimeDims2%>%dplyr::select(ARM ,Day  ,  V1_disc,prop_resp),by=c("ARM","Day","V1_disc")))
#save(summar_u_datFULL_ALLtimeDims3,summar_u_datFULL_ALLtimeDims2 ,file="/Users/jason/Dropbox/FELINE Project/Data_analysis/scRNA/15_subclone_interaction/Visualising MAPK transduction major axis.RData")

#library("pheatmap")
plt<- summar_u_datFULL_ALLtimeDims3[var!="V2_disc"][ARM!="A"][Day==180]
#pheatmap(spread(plt[ARM=="B"] ,var,     val_mu)%>%select(-c(ARM:prop_resp)), scale = "row")
#pheatmap(cor(spread(plt[ARM=="B"] ,var,     val_mu)%>%select(-c(ARM:prop_resp)) ,method="kendall") , scale = "row")
CM <- cor(spread(plt[ARM=="B"] ,var,     val_mu)%>%dplyr::select(-c(ARM:prop_resp)) ,method="kendall")

par(mfrow=c(1,1))
hmod<-hclust(dist(CM))
CM2<-CM[c("MAPK8","MAPKAP1","FGFR2","MAP3K5","MAPK10","ERBB4","ESR1","MAP2K4","FOS","JUNB","MAP3K8","MAP3K14","MAP3K2","MAP3K20","MAP2K5","MAP3K1","MAP4K3"),
        c("MAPK8","MAPKAP1","FGFR2","MAP3K5","MAPK10","ERBB4","ESR1","MAP2K4","FOS","JUNB","MAP3K8","MAP3K14","MAP3K2","MAP3K20","MAP2K5","MAP3K1","MAP4K3")]

fviz_dend(hmod, cex = 1.5, main = "",xlab = "", ylab = "", sub = "",k = 2,rect=TRUE,rect_fill=TRUE,k_colors = c("#FC4E07","#2E9FDF"), #, "#E7B800", "#00AFBB"),
          color_labels_by_k = TRUE,  # color labels by groups
          ggtheme = theme_void()    )

CM <- cor(spread(plt[ARM=="B"] ,var,     val_mu)%>%dplyr::select(-c(ARM:prop_resp,FOS,       JUNB)) ,method="kendall")

par(mfrow=c(1,1))
hmod<-hclust(dist(CM))
CM2<-CM[c( "MAPK8","FGFR2","MAPK10","MAP3K5","MAPKAP1","ERBB4","ESR1","MAP3K8","MAP3K14","MAP3K2","MAP3K1","MAP4K3","MAP3K20", "MAP2K4","MAP2K5"),
        c( "MAPK8","FGFR2","MAPK10","MAP3K5","MAPKAP1","ERBB4","ESR1","MAP3K8","MAP3K14","MAP3K2","MAP3K1","MAP4K3","MAP3K20", "MAP2K4","MAP2K5")]

fviz_dend(hmod, cex = 1.5, main = "",xlab = "", ylab = "", sub = "",k = 2,rect=TRUE,rect_fill=TRUE,k_colors = c("#FC4E07","#2E9FDF"), #, "#E7B800", "#00AFBB"),
          color_labels_by_k = TRUE,  # color labels by groups
          ggtheme = theme_void()    )
#library("corrplot")
# 9 x 7
corrplot::corrplot(CM2, order = "original"         ,diag=FALSE, addrect = 3,type = "lower",method= "color",tl.cex = 2,
                   col=colorRampPalette(c("blue","white","red"))(200))


corrplot::corrplot(CM2, order = "original"         ,diag=FALSE, addrect = 3,type = "lower",method= "color",tl.cex = 2,
                   col=rev(colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                              "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                              "#4393C3", "#2166AC", "#053061"))(100)) )

ggplot(   summar_u_datFULL_ALLtimeDims3[var!="V2_disc"][ARM!="A"][Day==180], aes(y=val_mu,x=prop_resp,group=Response,col=Response))+#
  geom_point(size=1.1,alpha=1)+
  facet_grid(ARM~var)+
  theme_classic(base_size=22)+
  geom_smooth(method="gam")



ggplot(   summar_u_datFULL_ALLtimeDims2[ARM!="A"][Day==180], aes(y=val_mu,x=V1_disc,group=var,col=var))+#
  geom_point(size=1.1,alpha=1)+facet_wrap(ARM~Response)+
  theme_classic(base_size=22)+
  geom_smooth(method="loess")




