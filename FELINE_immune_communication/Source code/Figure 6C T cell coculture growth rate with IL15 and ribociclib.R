rm(list=ls())
require(data.table);require(dplyr);require(ggplot2);require(tidyr);
require(mgcv);
require(ggsci)

savloc<- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/Figure6/"
alldd <- data.table(read.csv(file=paste0(savloc,"SourceData_Figure6_TcellCocultureDrugImpact.csv")))

# order label information
alldd[,Cellcatlab:=  paste(CellLine,Lineagelab,sep=" ")]

alldd$Cellcatlab <- factor(alldd$Cellcatlab , levels=c(
  "CAMA1 Resistant" , "MDAMB134 Sensitive",
  "CAMA1 Sensitive","MDAMB134 Resistant"
))

alldd$Cellcatlab <- factor(alldd$Cellcatlab , levels=c(
  "CAMA1 Resistant" , "CAMA1 Sensitive",
  "MDAMB134 Resistant","MDAMB134 Sensitive"
))

alldd$IL15Level <- factor(alldd$IL15Level , levels=c(
  "IL-15: 0 ngmL","IL-15: 0.1 ngmL" ,"IL-15: 0.25 ngmL" , "IL-15: 0.5 ngmL","IL-15: 0.75 ngmL",
  "IL-15: 1 ngmL","IL-15: 2.5 ngmL","IL-15: 5 ngmL","IL-15: 10 ngmL"  
))

# order composition
alldd$Composition <- factor(alldd$Composition , levels=c(  "Cancer and T cells","Cancer alone" ))
alldd$InitialTcellFraction<- round(alldd$InitialTcellFraction,2)

alldd[,MaxHour:=max(Hour), by=c("Expt","CellLine","RepID","Rep")]
alldd[,MinHour:=min(Hour), by=c("Expt","CellLine","RepID","Rep")]

ggplot(alldd[][InitialTcellFraction%in%c(0,0.2)][Ribociclib_uM%in%c(0,1)][IL15_ngmL%in%c(0,1,5)], aes(x = Hour, y= log(1+Count), col= Treatlab)) +
  geom_point(size=3,aes(shape=Composition,
                        col = Treatlab)) +
  geom_path(aes(linetype=Composition,
                group=interaction(RepID,Cellcatlab)))+
  facet_grid(Cellcatlab~paste0("IL-15:",IL15_ngmL,"(ngmL)"))+ 
  #stat_summary(geom="line", fun.data= mean_se, aes(col= Ribociclib_uM)) +
  #stat_summary(geom="errorbar", fun.data= mean_se, aes(col= Ribociclib_uM)) +
  labs(x= "Time (Hours)", y= "Cancer cell number (log(1+x))")+
  scale_color_jco(name="Treatment")+  scale_fill_jco(name="Treatment")+
  theme_classic(base_size=26)+
  theme(aspect.ratio=1)
#ggsave(file=paste0(svloc0,"TcellcocultureTimeCourse_IL15_ribo.png"), dpi=320,width=16,height=16)

ggplot(alldd[][InitialTcellFraction%in%c(0,0.2)][Ribociclib_uM%in%c(0,1)][IL15_ngmL%in%c(0,1,5)], aes(x = Hour, y= log(1+Count), col= Treatlab)) +
  geom_point(size=3,aes(shape=Composition,
                        col = Treatlab)) +
  facet_grid(Cellcatlab~paste0("IL-15:",IL15_ngmL,"(ngmL)"))+ 
  geom_smooth(aes(fill=(Treatlab),linetype=Composition,group=interaction(Treatlab,Composition)
  ),
  method="gam",alpha=0.4, formula=y~s(x,k=3))+
  labs(x= "Time (Hours)", y= "Cancer cell number (log(1+x))")+
  scale_color_jco(name="Treatment")+  scale_fill_jco(name="Treatment")+
  theme_classic(base_size=26)+
  theme(aspect.ratio=1)



alldd[, IL15lab:=paste0("IL-15:",IL15_ngmL,"(ng/mL)") ]
alldd$IL15lab <- factor(alldd$IL15lab ,levels=c(
  "IL-15:0(ng/mL)","IL-15:0.5(ng/mL)", "IL-15:0.75(ng/mL)","IL-15:1(ng/mL)","IL-15:2.5(ng/mL)","IL-15:5(ng/mL)"   
))
ggplot(alldd[][InitialTcellFraction%in%c(0,0.2)][Ribociclib_uM%in%c(0,1)][], aes(x = Hour, y= log(1+Count), col= Treatlab)) +
  geom_point(size=3,aes(shape=Composition,
                        col = Treatlab)) +
  facet_grid(IL15lab~Cellcatlab)+ 
  geom_smooth(aes(fill=(Treatlab),linetype=Composition,group=interaction(Treatlab,Composition)
  ),
  method="gam",alpha=0.4, formula=y~s(x,k=3))+
  labs(x= "Time (Hours)", y= "Cancer cell number (log(1+x))")+
  scale_color_jco(name="Treatment")+  scale_fill_jco(name="Treatment")+
  theme_classic(base_size=26)+
  theme(aspect.ratio=1)
#ggsave(file=paste0(svloc0,"TcellcocultureTimeCourseSmoothed_IL15_riboFull.png"), dpi=320,width=20.5,height=20.5)

FS18dd<-alldd[][InitialTcellFraction%in%c(0,0.2)][Ribociclib_uM%in%c(0,1)]
savloc<-"/Users/jason/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/SI data/FigureS18/"
write.csv(FS18dd,file=paste0(savloc,"SourceData_FigureS18_CocultureCancerTcellGrowthIL15Ribo.csv"))
# Calculate average growth rate of cancer population throughout treatment
alldd[,Count0:=sum((Hour==MinHour)*Count,na.rm=T),by=c("Expt","CellLine","RepID","Rep") ]
Finaldd <- alldd[Ribociclib_uM%in%c(0,1)][Hour==MaxHour]
Finaldd[,RGR:=(log(Count)-log(Count0))/(MaxHour-MinHour) ]
Finaldd<-Finaldd[!is.na(RGR)]

# Create modified variables for fitting
Finaldd[,sqrt_IL15:=sqrt(as.numeric(as.character(IL15_ngmL)))]
Finaldd[,FactorCellLineID:=as.factor(Cellcatlab)]

# Image of all the data
ggplot(Finaldd[InitialTcellFraction%in%c(0,0.2)],
       aes( y= RGR, x = sqrt_IL15, col = Treatlab,fill=Treatlab, 
            group=interaction(Treatlab,Composition,
                              InitialTcellFraction,CellLine,Lineagelab)
       )) +
  geom_smooth(aes(linetype=Composition),method="gam",alpha=0.4, formula=y~s(x,k=3))+
  geom_point(aes(col = Treatlab, shape=Composition,group=RepID),size=3) +
  facet_grid(CellLine~Lineagelab, scale="free")+theme_classic(base_size=26)+
  scale_color_jco(name="Treatment")+  scale_fill_jco(name="Treatment")+
  labs(y="Cancer cell growth rate (mean(cells/cell/day))", x="IL-15 (ngmL)")+ 
  theme(aspect.ratio=1)+scale_linewidth(range=c(1,2))+scale_size( range = c(1.5, 4.5))

#ggsave(file=paste0(svloc0,"SI TcellcocultureRGR_IL15_ribo_FitRAWgrowth.png"), dpi=320,width=16,height=16)
FS19dd<-Finaldd[InitialTcellFraction%in%c(0,0.2)]
savloc<-"/Users/jason/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/SI data/FigureS19/"
write.csv(FS19dd,file=paste0(savloc,"SourceData_FigureS19_MonoandCoculture_CancerGrowthRateIL15Ribo.csv"))

####################################
# Analysis of the data
#Composition = catag(InitialTcellFraction)   # FactorCellLineID=interaction(CellLine,Lineagelab)
statsdata <- Finaldd %>% dplyr::select(RepID,Rep,RGR,sqrt_IL15,Treatlab,Composition,InitialTcellFraction,Cellcatlab,CellLine,Lineagelab,FactorCellLineID)
statsdata[,RiboTreated:=0]
statsdata[Treatlab=="Ribociclib", RiboTreated:=1]
statsdata[,FacTreatlab:=as.factor(Treatlab)]
statsdata[, Ribo__sqrt_IL15:=RiboTreated*sqrt_IL15]#*(InitialTcellFraction==0)]
statsdata[, Ribo__InitialTcellFraction:=RiboTreated*InitialTcellFraction]

statsdata[, isCoculture:=0]
statsdata[InitialTcellFraction>0, isCoculture:=1]
statsdata[, isMonoculture:=0]
statsdata[InitialTcellFraction==0, isMonoculture:=1]

statsdata[, isControl:=0]
statsdata[Treatlab!="Ribociclib", isControl:=1]


####### GAM model and prediction
ModOut<-lapply(unique(statsdata$FactorCellLineID),function(x){
  cellXanalysis <- statsdata[FactorCellLineID==x][InitialTcellFraction%in%c(0,0.2)]
  
  cellXanalysis[, Ribo__sqrt_IL15:=RiboTreated*sqrt_IL15]#*(InitialTcellFraction==0)]
  cellXanalysis[, isCoculture_sqrt_IL15:=isCoculture*sqrt_IL15]
  cellXanalysis[, Ribo__isCoculture_sqrt_IL15:=RiboTreated*isCoculture*sqrt_IL15]
  
  gam0 <- gam(RGR~ RiboTreated+  # Monoculture ribo effect
                s(sqrt_IL15, k=3,bs="ts")+  # Monoculture IL-15 effect
                s(Ribo__sqrt_IL15, k=3,bs="ts")+ # Monoculture ribo + IL-15 synergy
                RiboTreated:isCoculture + # Coculture impact on ribo effect (change in ribo effect in coculutre without IL15)
                s(isCoculture_sqrt_IL15, k=3,bs="ts")+ # Coculture impact on IL15 effect
                s(Ribo__isCoculture_sqrt_IL15, k=3,bs="ts"),# Coculture impact on ribo + IL-15 synergy
              data=cellXanalysis)
  # Predict using model
  cellXanalysis[ , preds:=predict(gam0, type="response")]
  cellXanalysis[ , ucl:=predict(gam0, type="response")+predict(gam0,se=T, type="response")$se.fit]
  cellXanalysis[ , lcl:=predict(gam0, type="response")-predict(gam0,se=T, type="response")$se.fit]
  
  # Extract terms
  gam_terms <- cbind("(Intercept)"=coef(gam0)["(Intercept)"] ,predict(gam0,type="terms" )) #plot(rowSums(gam_terms),predict(gam0))
  colnames(gam_terms)<-paste0("gam",c("Intercept", "MonoRiboEffect", "MonoIL15Effect","MonoRiboIL15synergy",
                                      "CoRiboEffect", "CoIL15Effect","CoRiboIL15synergy"))
  cellXanalysisTerms <- cbind(cellXanalysis,gam_terms)
  
  # Make fine scale predictions
  lu<- data.table(sqrt_IL15= unique(c(seq(min(cellXanalysisTerms$sqrt_IL15),max(cellXanalysisTerms$sqrt_IL15),length=50),
                                      unique(cellXanalysisTerms$sqrt_IL15)))
  )
  meta2<-unique(cellXanalysis%>%dplyr::select(Lineagelab,CellLine,RiboTreated,Treatlab,Composition,isCoculture))
  
  predsdd<- rbindlist(lapply(1:nrow(lu),function(i){
    data.table(lu[i], meta2)  
  }))
  
  predsdd[, Ribo__sqrt_IL15:=RiboTreated*sqrt_IL15]
  predsdd[, isCoculture_sqrt_IL15:=isCoculture*sqrt_IL15]
  predsdd[, Ribo__isCoculture_sqrt_IL15:=RiboTreated*isCoculture*sqrt_IL15]
  
  
  predsdd[ , preds:=predict(gam0, newdata=predsdd ,type="response")]
  predsdd[ , ucl:=predict(gam0, newdata=predsdd ,type="response")+predict(gam0,se=T, newdata=predsdd ,type="response")$se.fit]
  predsdd[ , lcl:=predict(gam0, newdata=predsdd ,type="response")-predict(gam0,se=T, newdata=predsdd ,type="response")$se.fit]
  
  list(gam0,cellXanalysisTerms,meta2,predsdd)
  
})


allgammods<- lapply(ModOut,"[[",1)

allXanalysisTerms<- rbindlist(lapply(ModOut,"[[",2))
allpredsdd<- rbindlist(lapply(ModOut,"[[",4))

# Goodness of fit
ggplot(allXanalysisTerms,aes( y = preds, x= RGR,col = Treatlab,fill=Treatlab, 
                              group=interaction(Treatlab,Composition,InitialTcellFraction,CellLine,Lineagelab))) +
  geom_abline(linetype="dashed")+
  geom_smooth(aes(linewidth=(InitialTcellFraction/10)),method="gam",se=F,alpha=0.6, formula=y~s(x,k=3))+
  geom_point(aes(col = Treatlab, shape=Composition,group=RepID)) +
  geom_point(aes(col = Treatlab, shape=Composition,group=RepID,size=InitialTcellFraction)) +
  facet_grid(CellLine~Lineagelab)+
  theme_classic(base_size=26)+
  scale_color_jco(name="Treatment")+  scale_fill_jco(name="Treatment")+
  labs(y="Predicted rowth rate \n (mean(cells/cell/day))", x="Observed rowth rate \n (mean(cells/cell/day))")+ 
  theme(aspect.ratio=1)+scale_linewidth(range=c(1,2))+scale_size( range = c(1.5, 4.5))


# Predictions overlayed on data

allXanalysisTerms[,Cellcatlab:=  paste(CellLine,Lineagelab,sep=" ")]
allpredsdd[,Cellcatlab:=  paste(CellLine,Lineagelab,sep=" ")]

normvals<-data.table(allpredsdd[Treatlab=="Control"][Composition=="Cancer alone"]%>%group_by(Lineagelab,CellLine,sqrt_IL15)%>%summarise(predsNorm=mean(preds)))
allXanalysisTermsB<-merge(allXanalysisTerms,normvals,by=c("Lineagelab", "CellLine" , "sqrt_IL15"))
allpredsddB<-merge(allpredsdd,normvals,by=c("Lineagelab", "CellLine" , "sqrt_IL15"))

normvals<-data.table(allpredsdd[Treatlab=="Control"][Composition=="Cancer alone"][sqrt_IL15==0]%>%group_by(Lineagelab,CellLine)%>%summarise(predsNorm=mean(preds)))
allXanalysisTermsB<-merge(allXanalysisTerms,normvals,by=c("Lineagelab", "CellLine" ))
allpredsddB<-merge(allpredsdd,normvals,by=c("Lineagelab", "CellLine" ))

normvals<-data.table(allpredsdd[][Composition=="Cancer alone"][sqrt_IL15==0]%>%group_by(Treatlab,Lineagelab,CellLine)%>%summarise(predsNorm=mean(preds)))
allXanalysisTermsB<-merge(allXanalysisTerms,normvals,by=c("Lineagelab", "CellLine","Treatlab" ))
allpredsddB<-merge(allpredsdd,normvals,by=c("Lineagelab", "CellLine","Treatlab"  ))

svloc00<-"~/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"

ggplot(allXanalysisTermsB,
       aes( y= RGR-predsNorm, x = sqrt_IL15, col = Treatlab,fill=Treatlab, 
            group=interaction(Treatlab,Composition,
                              CellLine,Lineagelab)
       )) +
  geom_ribbon(data=allpredsddB ,alpha=0.2,aes(col=NA,y=preds-predsNorm,ymax=ucl-predsNorm,ymin=lcl-predsNorm))+
  geom_line(data=allpredsddB ,aes(linetype=Composition,y=preds-predsNorm))+#scale_linewidth(range=c(1,2))+scale_size( range = c(1.5, 4.5))+
  geom_point(size=3,alpha=1,aes(col = Treatlab, shape=Composition,group=RepID)) +
  facet_grid(CellLine~Lineagelab, scale="free")+theme_classic(base_size=26)+
  scale_color_jco(name="Treatment")+  scale_fill_jco(name="Treatment")+
  scale_x_continuous(labels = c(0,0.5,1,2,4), breaks=sqrt(c(0,0.5,1,2,4)) )+
  labs(y="Cancer growth rate \n (relative to ribociclib treatment \n monoculture mean of IL-15 control)", x="IL-15 (ngmL)")+ 
  theme(aspect.ratio=1)
#ggsave(file= paste0(svloc00,"TcellcocultureRGR_IL15_ribo_FitRelativegrowth.png"), dpi=320,width=11,height=11)

ggplot(allXanalysisTermsB,
       aes( y= RGR-predsNorm, x = sqrt_IL15, col = Treatlab,fill=Treatlab, 
            group=interaction(Treatlab,Composition,
                              CellLine,Lineagelab)
       )) +
  geom_ribbon(data=allpredsddB ,alpha=0.2,aes(col=NA,y=preds-predsNorm,ymax=ucl-predsNorm,ymin=lcl-predsNorm))+
  geom_line(data=allpredsddB ,aes(linetype=Composition,y=preds-predsNorm))+#scale_linewidth(range=c(1,2))+scale_size( range = c(1.5, 4.5))+
  geom_point(size=3,alpha=1,aes(col = Treatlab, shape=Composition,group=RepID)) +
  facet_grid(CellLine~Lineagelab, scale="free")+theme_classic(base_size=26)+
  scale_color_jco(name="Treatment")+  scale_fill_jco(name="Treatment")+
  scale_x_continuous(labels = c(0,0.5,1,2,4), breaks=sqrt(c(0,0.5,1,2,4)) )+
  labs(y="Cancer growth rate \n (relative to ribociclib treatment \n monoculture mean of IL-15 control)", x="IL-15 (ngmL)")+ 
  theme(aspect.ratio=1)+
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        strip.background = element_blank(),
        strip.text= element_blank(),
        legend.position = "none")
#ggsave(file=paste0(svloc00,"BLANK_TcellcocultureRGR_IL15_ribo_FitRelativegrowth.png"), dpi=320,width=8,height=8)




resout<-allXanalysisTermsB%>%dplyr::select(CellLine,Lineagelab,Treatlab,Composition,RepID ,Rep,
                                           RiboTreated, sqrt_IL15, InitialTcellFraction,
                                           RGR,predsNorm,preds,ucl,lcl)
write.csv(resout,file=paste0(savloc,"Outputs/","SourceData_Figure6_TcellCocultureDrugImpact_Out.csv"))




ggplot(allXanalysisTermsB,
       aes( y= RGR-predsNorm, x = sqrt_IL15, col = Treatlab,fill=Treatlab, 
            group=interaction(Treatlab,Composition,
                              CellLine,Lineagelab)
       )) +
  geom_ribbon(data=allpredsddB ,alpha=0.2,aes(col=NA,y=preds-predsNorm,ymax=ucl-predsNorm,ymin=lcl-predsNorm))+
  geom_line(data=allpredsddB ,aes(linetype=Cellcatlab,y=preds-predsNorm))+#scale_linewidth(range=c(1,2))+scale_size( range = c(1.5, 4.5))+
  geom_point(size=3,alpha=1,aes(col = Treatlab, shape=Cellcatlab,group=RepID)) +
  facet_grid(
    CellLine~Composition, scale="free")+
  theme_classic(base_size=26)+
  scale_color_jco(name="Treatment")+  scale_fill_jco(name="Treatment")+
  scale_x_continuous(labels = c(0,0.5,1,2,4), breaks=sqrt(c(0,0.5,1,2,4)) )+
  labs(y="Cancer growth rate", x="IL-15 (ngmL)")+ 
  theme(aspect.ratio=1)



ggplot(allXanalysisTerms,
       aes( y= RGR, x = sqrt_IL15, col = Treatlab,fill=Treatlab, 
            group=interaction(Treatlab,Composition,
                              CellLine,Lineagelab)
       )) +
  geom_ribbon(data=allpredsdd ,alpha=0.2,aes(col=NA,y=preds,ymax=ucl,ymin=lcl))+
  geom_line(data=allpredsdd ,aes(linetype=Composition,y=preds))+#scale_linewidth(range=c(1,2))+scale_size( range = c(1.5, 4.5))+
  geom_point(alpha=1,aes(col = Treatlab,size=Composition, shape=Composition,group=RepID)) +
  facet_grid(CellLine~Lineagelab, scale="free")+theme_classic(base_size=26)+
  scale_color_jco(name="Treatment")+  scale_fill_jco(name="Treatment")+
  scale_x_continuous(labels = c(0,0.5,1,2,4), breaks=sqrt(c(0,0.5,1,2,4)) )+
  labs(y="Cancer growth rate", x="IL-15 (ng/mL)")+ 
  theme(aspect.ratio=1)+
  scale_size_manual(values=c(4,2))+
  theme(axis.text.x=element_text(size=rel(0.9)))

#ggsave(file=paste0(svloc00,"SI TcellcocultureRGR_IL15_ribo_Fit.png"), dpi=320,width=11,height=11)


ggplot(allXanalysisTerms,
       aes( y= RGR, x = sqrt_IL15, col = Treatlab,fill=Treatlab, 
            group=interaction(Treatlab,Composition,
                              CellLine,Lineagelab)
       )) +
  geom_ribbon(data=allpredsdd ,alpha=0.2,aes(col=NA,y=preds,ymax=ucl,ymin=lcl))+
  geom_line(data=allpredsdd ,aes(linetype=Composition,y=preds))+
  geom_point(alpha=1,aes(col = Treatlab,size=Composition, shape=Composition,group=RepID)) +
  facet_grid(CellLine~Lineagelab, scale="free")+theme_classic(base_size=26)+
  scale_color_jco(name="Treatment")+  scale_fill_jco(name="Treatment")+
  scale_x_continuous(labels = c(0,0.5,1,2,4), breaks=sqrt(c(0,0.5,1,2,4)) )+
  labs(y="Cancer growth rate", x="IL-15 (ngmL)")+ 
  theme(aspect.ratio=1)+
  scale_size_manual(values=c(4,2))+
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        strip.background = element_blank(),
        strip.text= element_blank(),
        legend.position = "none")

ggsave(file=paste0(svloc,"BLANK_TcellcocultureRGR_IL15_ribo_Fit.png"), dpi=320,width=8,height=8)



# stats
ss1<-rbindlist( lapply(1:length(allgammods), function(g){
  data.table( ModOut[[g]][[2]][1] %>%dplyr::select(Cellcatlab:FactorCellLineID),
              data.table(summary(allgammods[[g]])$s.table,keep.rownames = T)[rn=="s(isCoculture_sqrt_IL15)"]
  )
}) )

ss2<-rbindlist( lapply(1:length(allgammods), function(g){
  data.table( ModOut[[g]][[2]][1] %>%dplyr::select(Cellcatlab:FactorCellLineID),
              data.table(summary(allgammods[[g]])$s.table,keep.rownames = T)[rn=="s(Ribo__isCoculture_sqrt_IL15)"]
  )
}) )

#(IL-15 effect in coculutre:CAMA-1 resistant:edf=2.00, F=1323.00,p<1e-16,CAMA-1 sensitive:edf=2.00, F=543.25,p<1e-16,MDAMB134 resistant:edf=1.25, F=173.77,p<1e-16,MDAMB134 sensitive:edf=1.70, F=684.73,p<1e-16)
#(ribociclib impact on IL-15 effect in coculutre:CAMA-1 resistant:edf=1.88, F=56.15,p<1e-16,CAMA-1 sensitive:edf=1.00, F=9.61,p=2.7e-5,MDAMB134 resistant:edf=1.13, F=6.22,p=5.5e-16,MDAMB134 sensitive:edf=1.94, F=51.74,p<e-16)

ss<- rbind(ss1,ss2)
#write.csv(ss,    "/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/temporary/TcellcocultureRGR_IL15_ribo_Stats.csv")




ggplot(allXanalysisTerms[isCoculture==T][Treatlab=="Ribociclib"],
       aes( gamCoRiboEffect + #gamCoIL15Effect + 
              gamCoRiboIL15synergy,
            x = sqrt_IL15, col=Lineagelab,
            group=interaction(CellLine,Lineagelab)
       )) +
  geom_hline(linetype="dashed",yintercept=0)+
  geom_point(size=3,aes(shape=CellLine)) +
  geom_path(size=1.5 ) +
  theme_classic(base_size=26)+
  scale_x_continuous(labels = c(0,0.5,1,2,4), breaks=sqrt(c(0,0.5,1,2,4)) )+
  labs(y="Ribociclib effect on \n cocultured cancer growth rate", x="IL-15 (ngmL)")+ 
  theme(aspect.ratio=1)+
  scale_color_npg(name="Ribociclib \n resistance")

ggsave(file=paste0(svloc,"TcellcocultureRGR_IL15_ribo_RiboEffect.png"), dpi=320,width=9,height=9)


ggplot(allXanalysisTerms[isCoculture==T][Treatlab=="Ribociclib"],
       aes( gamCoRiboIL15synergy,
            x = sqrt_IL15, col=Lineagelab,
            group=interaction(CellLine,Lineagelab)
       )) +
  geom_hline(linetype="dashed",yintercept=0)+
  geom_point(size=3,aes(shape=CellLine)) +
  geom_path(size=1.5 ) +
  theme_classic(base_size=26)+
  scale_x_continuous(labels = c(0,0.5,1,2,4), breaks=sqrt(c(0,0.5,1,2,4)) )+
  labs(y="Ribociclib effect on \n cocultured cancer growth rate", x="IL-15 (ngmL)")+ 
  theme(aspect.ratio=1)+
  scale_color_npg(name="Ribociclib \n resistance")


