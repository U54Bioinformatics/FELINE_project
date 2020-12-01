

## Default plot method
plot(out)
MAPK_Major_axis<-data.table(u_datFULL_ALLtimeDims%>%dplyr::select(Cell.ID,Patient.Study.ID,ARM,Day,Response,Subclone,V1,V2,V1_disc,    V2_disc))
#save(MAPK_Major_axis,which_plot,corV3,subset_u_datFULL_ALLtimeDims,file="/Users/jason/Dropbox/FELINE Project/Data_analysis/scRNA/15_subclone_interaction/MAPK major axis proj.RData")

corV3<- as.data.table(cor(u_datFULL_ALLtimeDims%>%dplyr::select(-c(Cell.ID, Patient.Study.ID, ARM, Day,  Response, Subclone, Phase,V1:dimension2)),
                          u_datFULL_ALLtimeDims%>%dplyr::select(V1:dimension2) ),keep.rownames=TRUE)
corV3[order(-abs(dimension1))][1:10]
corV3[order(-abs(dimension2))][1:10]
corV3[order(-abs(V1))][1:10]
corV3[order(-abs(V2))][1:10]
corV3[order(-abs(V5))][1:10]

u_datFULL_ALLtimeDims[,V1_disc:=(round(V1*10)/10) ]
u_datFULL_ALLtimeDims[,V2_disc:=(round(V2*3)/3) ]

which_plot<-unique( c(corV3[order(-abs(V1))][1:10]$rn,  corV3[order(-abs(V2))][1:10]$rn) )
subset_u_datFULL_ALLtimeDims <- data.table(u_datFULL_ALLtimeDims%>% #[ARM!="A"]
                                             dplyr::select(c("Patient.Study.ID", "ARM", "Day",  "Response", "V1_disc"),which_plot)%>%
                                             gather(var,val,-c(Patient.Study.ID, ARM, Day,  Response,V1_disc)) )
summar_u_datFULL_ALLtimeDims <-data.table(subset_u_datFULL_ALLtimeDims   %>% group_by(Patient.Study.ID, ARM, Day,  Response,V1_disc,var)%>%
                                            dplyr::summarise(val_mu=mean(val),ncell=length(val)))

subset_u_datFULL_ALLtimeDims[var=="FOS"][Day==180][Response=="Responder"][V1_disc==median(V1_disc)]
subset_u_datFULL_ALLtimeDims[var=="FOS"][Day==180][Response=="Non-responder"][V1_disc==median(V1_disc)]


summar_u_datFULL_ALLtimeDims2 <- data.table( subset_u_datFULL_ALLtimeDims%>% group_by( Response,ARM, Day,V1_disc)%>%
                                               summarise(ncell=length(val))%>%spread(Response,ncell)%>%mutate(prop_resp=Responder/(Responder+`Non-responder` )))
summar_u_datFULL_ALLtimeDims2[is.na(Responder),Responder:=0]
summar_u_datFULL_ALLtimeDims2[is.na(`Non-responder`),`Non-responder`:=0]
summar_u_datFULL_ALLtimeDims2[,prop_resp:=Responder/(Responder+`Non-responder` )]
summar_u_datFULL_ALLtimeDims2[ARM=="A"][Day==180]$prop_resp*  1/(sum(summar_u_datFULL_ALLtimeDims2[ARM=="A"][Day==180]$prop_resp)/length( summar_u_datFULL_ALLtimeDims2[ARM=="A"][Day==180]))
summar_u_datFULL_ALLtimeDims2[ARM=="C"][Day==180]$prop_resp*  1/(sum(summar_u_datFULL_ALLtimeDims2[ARM=="C"][Day==180]$prop_resp)/length( summar_u_datFULL_ALLtimeDims2[ARM=="C"][Day==180]))
summar_u_datFULL_ALLtimeDims2[ARM=="A"][Day==180]$prop_resp2
1/(sum(summar_u_datFULL_ALLtimeDims2[ARM=="C"][Day==180]$prop_resp)/length( summar_u_datFULL_ALLtimeDims2[ARM=="C"][Day==180]))
summar_u_datFULL_ALLtimeDims2[ARM=="A"][Day==180]$prop_resp*  1/(mean(summar_u_datFULL_ALLtimeDims2[ARM=="A"][Day==180]$prop_resp))

#summar_u_datFULL_ALLtimeDims2[,prop_resp2:=0.5*prop_resp /(mean(prop_resp)),by=c("ARM","Day")]
summar_u_datFULL_ALLtimeDims2[,prop_resp2:=prop_resp/(0.2+mean(prop_resp))/2,by=c("ARM","Day")]
summar_u_datFULL_ALLtimeDims2[,mu_prop_resparmday:=(mean(prop_resp)),by=c("ARM","Day")]

summar_u_datFULL_ALLtimeDims4<- merge(summar_u_datFULL_ALLtimeDims2,
                                      u_datFULL_ALLtimeDims%>%group_by(ARM,Day)%>%summarise(s=sum(Response=="Responder")/sum(Response%in%c("Responder","Non-responder"))),
                                      by=c("ARM",     "Day"))
summar_u_datFULL_ALLtimeDims4[,prop_resp2:=prop_resp/(s*2),by=c("ARM","Day")]
summar_u_datFULL_ALLtimeDims4[,prop_resp2:=prop_resp]
summar_u_datFULL_ALLtimeDims4[ARM=="A",prop_resp2:=0.25*prop_resp/s]

ggplot(   summar_u_datFULL_ALLtimeDims4[Day==180], aes(x=V1_disc,y=(prop_resp2),col=ARM,fill=ARM))+#
  geom_point(size=2,alpha=1)+
  facet_grid(~ARM)+
  theme_classic(base_size=26)+theme(aspect.ratio=1)+ scale_color_brewer(palette = "Dark2")+ scale_fill_brewer(palette = "Dark2") +
  #+ scale_color_grey(start = 0.8, end = 0.2) + scale_fill_grey(start = 0.8, end = 0.2) +
  geom_smooth(method="gam",formula=y~s(x,k=4), method.args=list(family="quasibinomial")) +labs(y="Proportion of cells from responders", x="MAPK activation phenotpye" )+
  theme(legend.position="none")
u_datFULL_ALLtimeDims%>%group_by(ARM,Day)%>%summarise(s=sum(Response=="Responder")/sum(Response%in%c("Responder","Non-responder")))


ggplot(   summar_u_datFULL_ALLtimeDims4[Day==180], aes(x=V1_disc,y=0.9*(prop_resp2),col=ARM,fill=ARM))+#
  geom_point(size=2,alpha=1)+
  facet_grid(~ARM)+
  theme_classic(base_size=26)+theme(aspect.ratio=1)+ scale_color_brewer(palette = "Dark2")+ scale_fill_brewer(palette = "Dark2") +
  #+ scale_color_grey(start = 0.8, end = 0.2) + scale_fill_grey(start = 0.8, end = 0.2) +
  geom_smooth(method="gam",formula=y~s(x,k=4), method.args=list(family="quasibinomial")) +labs(y="Proportion of cells from responders", x="MAPK activation phenotpye" )+
  theme(legend.position="none")
u_datFULL_ALLtimeDims%>%group_by(ARM,Day)%>%summarise(s=sum(Response=="Responder")/sum(Response%in%c("Responder","Non-responder")))


#   group_by( Response,ARM, Day,V1_disc)
#summar_u_datFULL_ALLtimeDims3 <-  data.table(summar_u_datFULL_ALLtimeDims           %>%
#               group_by( ARM, Day,V1_disc)%>%  
#               mutate(nsamp=sum(ncell) )%>%
#               group_by( Response,ARM, Day,V1_disc)%>%
#               mutate(prop_resp=ncell/nsamp)  
#             )


summar_u_datFULL_ALLtimeDims3 <- data.table(merge(summar_u_datFULL_ALLtimeDims,summar_u_datFULL_ALLtimeDims2%>%dplyr::select(ARM ,Day  ,  V1_disc,prop_resp),by=c("ARM","Day","V1_disc")))
#save(summar_u_datFULL_ALLtimeDims3,summar_u_datFULL_ALLtimeDims2 ,file="/Users/jason/Dropbox/FELINE Project/Data_analysis/scRNA/15_subclone_interaction/Visualising MAPK transduction major axis.RData")
load(file="/Users/jason/Dropbox/FELINE Project/Data_analysis/scRNA/15_subclone_interaction/Visualising MAPK transduction major axis.RData")

ggplot(   summar_u_datFULL_ALLtimeDims2[Day==180], aes(x=V1_disc,y=(prop_resp),col=ARM,fill=ARM))+#
  geom_point(size=2,alpha=1)+
  facet_grid(~ARM)+
  theme_classic(base_size=26)+theme(aspect.ratio=1)+ scale_color_brewer(palette = "Dark2")+ scale_fill_brewer(palette = "Dark2") +
  #+ scale_color_grey(start = 0.8, end = 0.2) + scale_fill_grey(start = 0.8, end = 0.2) +
  geom_smooth(method="gam", method.args=list(family="quasibinomial")) +labs(y="Proportion of cells from responders", x="MAPK activation phenotpye" )+
  theme(legend.position="none")

ggplot(   summar_u_datFULL_ALLtimeDims3[var!="V2_disc"][ARM!="A"][Day==180], aes(y=val_mu, x=V1_disc, group=Response,col=Response))+#
  geom_point(size=1.1,alpha=1)+
  facet_grid(ARM~var)+
  theme_classic(base_size=22)+
  geom_smooth(method="gam")+theme(aspect.ratio=1)+
  labs(y="MAPK gene expression", x="MAPK activation phenotpye")+theme(aspect.ratio=1)

pl_sub<-summar_u_datFULL_ALLtimeDims3[ARM!="A"][Day==180][var%in%c("MAP3K1","MAP3K14","MAP3K5","MAP4K3","MAPK10","MAPK8")]
pl_sub[,slope:="neg"]
pl_sub[var%in%c("MAP3K5","MAPK10","MAPK8"),slope:="pos"]
pl_sub$var <- factor(pl_sub$var , levels = c("MAP3K5","MAPK10","MAPK8",
                                             "MAP3K1","MAP3K14","MAP4K3"
))

ggplot(   pl_sub, aes(y=val_mu, x=V1_disc, group=var,col=var))+#
  geom_point(size=1.1,alpha=1)+
  facet_wrap(~var,nrow=2)+
  theme_classic(base_size=22)+
  geom_smooth(method="gam",formula=y~s(x,k=5))+theme(aspect.ratio=1)+
  labs(y="MAPK gene expression", x="MAPK activation phenotpye")+ theme(legend.position="none")


pl_sub<-summar_u_datFULL_ALLtimeDims3[][Day==180][var%in%c("MAP3K1","MAP3K14","MAP3K5","MAP4K3","MAPK10","MAPK8")]
pl_sub[,slope:="neg"]
pl_sub[var%in%c("MAP3K5","MAPK10","MAPK8"),slope:="pos"]
pl_sub$var <- factor(pl_sub$var , levels = c("MAP3K5","MAPK10","MAPK8",
                                             "MAP3K1","MAP3K14","MAP4K3"
))
ggplot(   pl_sub[Day==180], aes(y=log(1+ncell), x=V1_disc, group=interaction(ARM,Response),col=ARM))+#
  geom_point(size=1.1,alpha=1)+
  facet_wrap(~var,nrow=2)+
  theme_classic(base_size=22)+
  geom_smooth(method="gam",formula=y~s(x,k=9))+theme(aspect.ratio=1)+
  labs(y="MAPK gene expression", x="MAPK activation phenotpye")+ theme(legend.position="none")


ggplot(   pl_sub, aes(y=val_mu, x=V1_disc, group=var,col=ARM))+#
  geom_point(size=1.1,alpha=1)+
  facet_wrap(~var,nrow=2)+
  theme_classic(base_size=22)+
  geom_smooth(method="gam",formula=y~s(x,k=5))+theme(aspect.ratio=1)+
  labs(y="MAPK gene expression", x="MAPK activation phenotpye")+ theme(legend.position="none")



library(factoextra)
library(corrplot)

plt<- summar_u_datFULL_ALLtimeDims3[var!="V2_disc"][Day==180]
#pheatmap(spread(plt[ARM=="B"] ,var,     val_mu)%>%select(-c(ARM:prop_resp)), scale = "row")
#pheatmap(cor(spread(plt[ARM=="B"] ,var,     val_mu)%>%select(-c(ARM:prop_resp)) ,method="kendall") , scale = "row")
CM<- cor(spread(plt ,var,     val_mu)%>%dplyr::select(-c(ARM:prop_resp)) ,method="kendall")
hmod<-hclust(dist(CM))
#7 x 14
fviz_dend(hmod, cex = 1.5, main = "",xlab = "", ylab = "", sub = "",k = 2,rect=TRUE,rect_fill=TRUE,k_colors = c("#2E9FDF", "#FC4E07"), #, "#E7B800", "#00AFBB"),
          color_labels_by_k = TRUE,  # color labels by groups
          ggtheme = theme_void()    )
#8x8
corrplot:corrplot(CM, order = "hclust",diag=FALSE, addrect = 3,type = "lower",method= "color",tl.cex=1.25)




#library("pheatmap")
plt<- summar_u_datFULL_ALLtimeDims3[var!="V2_disc"][ARM=="C"][Day==180]
#pheatmap(spread(plt[ARM=="B"] ,var,     val_mu)%>%select(-c(ARM:prop_resp)), scale = "row")
#pheatmap(cor(spread(plt[ARM=="B"] ,var,     val_mu)%>%select(-c(ARM:prop_resp)) ,method="kendall") , scale = "row")
CM<- cor(spread(plt[ARM=="C"] ,var,     val_mu)%>%select(-c(ARM:prop_resp)) ,method="kendall")
hmod<-hclust(dist(CM))
#7 x 14
fviz_dend(hmod, cex = 1.5, main = "",xlab = "", ylab = "", sub = "",k = 2,rect=TRUE,rect_fill=TRUE,k_colors = c( "#FC4E07","#2E9FDF"), #, "#E7B800", "#00AFBB"),
          color_labels_by_k = TRUE,  # color labels by groups
          ggtheme = theme_void()    )
#8x8
corrplot(CM, order = "hclust",diag=FALSE, addrect = 3,type = "lower",method= "color")


plt<- summar_u_datFULL_ALLtimeDims3[var!="V2_disc"][ARM=="A"][Day==180]

CM2<- cor(spread(plt[ARM=="A"] ,var,     val_mu)%>%select(-c(ARM:prop_resp)) ,method="kendall")
hmod2<-hclust(dist(CM2))

fviz_dend(hmod2, cex = 1.5, main = "",xlab = "", ylab = "", sub = "",k = 2,rect=TRUE,rect_fill=TRUE,k_colors = c("#2E9FDF", "#FC4E07"), #, "#E7B800", "#00AFBB"),
          color_labels_by_k = TRUE,  # color labels by groups
          ggtheme = theme_void()    )

corrplot(CM2, order = "hclust",diag=FALSE, addrect = 3,type = "lower",method= "color")

