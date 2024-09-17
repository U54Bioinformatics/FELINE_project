rm(list=ls())

### Load Ligand Receptor database list of Ramilowski et al 2015
LRpairsFiltered <- data.table(read.csv( file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/SI data/FigureS4/SourceData_Figure S4 LigandReceptorPairsRamilowski2015.csv"))
LRgenelist <- unique( c(LRpairsFiltered$HPMR.Receptor, LRpairsFiltered$HPMR.Ligand) )

### Load GeneOntology Cytokine Receptor List
cytokineReceptors <- read.csv(  file= "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/SI data/FigureS4/SourceData_Figure S4 GeneOntologyCytokineReceptors.csv")$x

### Discovery/validation cohort communications between cell types
savloc <- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/Figure5/"
TcellSubtypeCommCCI_all <- data.table(read.csv( file=paste0(savloc,"SourceData_Figure5_MyeloidSubtypetoTcellCommunication.csv")) )

### Discovery cohort communications between cell types
TcellSubtypeCommCCI <- data.table(TcellSubtypeCommCCI_all [Cohort== "Discovery"]%>% dplyr::select(-Cohort))

subsetInflam <- TcellSubtypeCommCCI[Pair.Name %in% LRpairsFiltered[
  HPMR.Receptor%in% cytokineReceptors
  ]$Pair.Name][LigandPhenoCelltype=="Macrophages"]
subsetInflam[,startval:=sum((Day==0)*log(1+scaleTransduction))/sum(Day==0), by=c("Treatmentlab","Pair.Name","Patient.Study.ID","ReceptorPhenoCelltype")]
subsetInflam[ ,startvalAv:=median(startval,na.rm=T), by=c("Treatmentlab","Pair.Name","dynamic_class3","ReceptorPhenoCelltype")]
subsetInflam[,"Tumorresponse":= "Resistant"]
subsetInflam[dynamic_class3=="Response","Tumorresponse":= "Sensitive"]

pairstotest <- unique(subsetInflam$Pair.Name)
resultout <- rbindlist(lapply( 1:length(pairstotest),function(x){
  print (x)
  res<-tryCatch({
    m1 <- lmer( log(1+scaleTransduction) ~ log(1+Day) * dynamic_class3 +(1|Patient.Study.ID), data=   subsetInflam[ARM!="A"][Pair.Name== pairstotest[x] ])
    data.table( Pair.Name=pairstotest[x],  data.table(summary(m1)$coef,keep.rownames=T) )
  },
  error=function(x){
    return(NULL)
  })
  res
}))
setnames(resultout,old=c("t value","Pr(>|t|)"),new=c("tval","pval"))
resultout$adj.pval <- p.adjust(p=resultout$pval,method="fdr")#,n=length(unique(resultout$Pair.Name)))

pathsR <- c("F13A1_ITGA4" ,  "VCAM1_ITGA4","VCAM1_ITGB1","F13A1_ITGB1","ICAM3_ITGB2",  
           "ADAM12_ITGB1",  "COL18A1_ITGB1", "THBS2_ITGA4" ,  "ICAM2_ITGB2")
pathsA <-c("IL15_IL2RB" ,   "IL18_IL18R1" ,"IL15_IL2RA"  ,  "IL18_IL18RAP","IL15_IL15RA")
paths<- c(pathsA,pathsR)
plotthis<- subsetInflam[ Pair.Name%in%paths]
plotthis[, Pair.Name2 :=gsub("_","-",Pair.Name) ]
plotthis[, Function := "Activation" ]
plotthis[Pair.Name%in%pathsR, Function := "Recruitment" ]
plotthis$Pair.Name2 <- factor(plotthis$Pair.Name2 , levels= gsub("_","-", c(sort(pathsA) ,
                                                                            sort(pathsR) )))                                              

ggplot( plotthis[Day==180]%>%group_by(Pair.Name2)%>%mutate(y=scale( log(1 + scaleTransduction))) , aes(y= y, Tumorresponse,group=interaction(Tumorresponse,Pair.Name2),col= Pair.Name2,fill= Pair.Name2)) +
  theme_classic(base_size=18) +
  geom_boxplot(alpha = 0.6,position=position_dodge(width=1))+
  #geom_smooth(alpha=1,col= "black",method="gam",formula=y~s(x,k=3),se=T) +
  theme(aspect.ratio=1)+
  geom_point(position=position_dodge(width=1))+ facet_wrap(Treatmentlab~Function)+
  geom_vline(xintercept = 1.5,linetype="dashed") +
  scale_color_discrete(name="Communication \n pathway")+
  scale_fill_discrete(name="Communication \n pathway")+
  labs(y="Myeloid to T cell \n communication post treatment",x="Tumor response")

plotthis[,Daylab:=paste0("Day ",Day)]
discdat <-plotthis[Function=="Activation"][Day==180]%>%group_by(Pair.Name2)%>%mutate(y=scale( log(1 + scaleTransduction))) 
load(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/temporary/DiscoveryMyeloidTcellActivation.RData")
Activationdiscdat <- data.table(Cohort="Discovery", discdat )

discdat <-plotthis[Function=="Recruitment"][Day==180]%>%group_by(Pair.Name2)%>%mutate(y=scale( log(1 + scaleTransduction)))
#save(discdat,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/temporary/DiscoveryMyeloidTcellRecruitment.RData")
load(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/temporary/DiscoveryMyeloidTcellRecruitment.RData")
Recruitmentdiscdat <- data.table(Cohort="Discovery", discdat )



# Stats
pairstotest <- c("IL15_IL2RA","IL15_IL2RB","IL18_IL18R1","IL18_IL18RAP", "IL15_IL15RA" )

activeoutRiboDisc <- rbindlist(lapply( 1:length(pairstotest),function(x){
  print (x)
  res<-tryCatch({
    m1 <- lmer( log(1+scaleTransduction) ~ log(1+Day) * dynamic_class3 +(1|Patient.Study.ID), data=   subsetInflam[ARM!="A"][Pair.Name== pairstotest[x] ])
    data.table( Pair.Name=pairstotest[x],  data.table(summary(m1)$coef,keep.rownames=T) )
  },
  error=function(x){
    return(NULL)
  })
  setnames(res,old=c("t value","Pr(>|t|)"),new=c("tval","pval"))
  res
}))[rn=="log(1 + Day):dynamic_class3Response"]

activeoutLetrozoleDisc <- rbindlist(lapply( 1:length(pairstotest),function(x){
  print (x)
  res<-tryCatch({
    m1 <- lmer( log(1+scaleTransduction) ~ log(1+Day) * dynamic_class3 +(1|Patient.Study.ID), data=   subsetInflam[ARM=="A"][Pair.Name== pairstotest[x] ])
    data.table( Pair.Name=pairstotest[x],  data.table(summary(m1)$coef,keep.rownames=T) )
  },
  error=function(x){
    return(NULL)
  })
  setnames(res,old=c("t value","Pr(>|t|)"),new=c("tval","pval"))
  res
}))[rn=="log(1 + Day):dynamic_class3Response"]



pairstotest <- c("ADAM12_ITGB1" , "COL18A1_ITGB1", "F13A1_ITGA4","F13A1_ITGB1","ICAM2_ITGB2"  , "ICAM3_ITGB2","VCAM1_ITGA4","VCAM1_ITGB1","THBS2_ITGA4")
recruitoutRiboDisc <- rbindlist(lapply( 1:length(pairstotest),function(x){
  print (x)
  res<-tryCatch({
    m1 <- lmer( log(1+scaleTransduction) ~ log(1+Day) * dynamic_class3 +(1|Patient.Study.ID), data=   subsetInflam[ARM!="A"][Pair.Name== pairstotest[x] ])
    data.table( Pair.Name=pairstotest[x],  data.table(summary(m1)$coef,keep.rownames=T) )
  },
  error=function(x){
    return(NULL)
  })
  setnames(res,old=c("t value","Pr(>|t|)"),new=c("tval","pval"))
  res
}))[rn=="log(1 + Day):dynamic_class3Response"]

recruitoutLetrozoleDisc <- rbindlist(lapply( 1:length(pairstotest),function(x){
  print (x)
  res<-tryCatch({
    m1 <- lmer( log(1+scaleTransduction) ~ log(1+Day) * dynamic_class3 +(1|Patient.Study.ID), data=   subsetInflam[ARM=="A"][Pair.Name== pairstotest[x] ])
    data.table( Pair.Name=pairstotest[x],  data.table(summary(m1)$coef,keep.rownames=T) )
  },
  error=function(x){
    return(NULL)
  })
  setnames(res,old=c("t value","Pr(>|t|)"),new=c("tval","pval"))
  res
}))[rn=="log(1 + Day):dynamic_class3Response"]

recruitoutRiboDisc[pval<0.05]
recruitoutLetrozoleDisc[pval<0.05]
activeoutRiboDisc[pval<0.05]
activeoutLetrozoleDisc[pval<0.05]

totalstatsDiscovery <- rbind( data.table(Cohort="Discovery",Treatment="Combination ribociclib" ,Function="Recruitment" , recruitoutRiboDisc),
                              data.table(Cohort="Discovery",Treatment="Letrozole alone" ,Function="Recruitment" , recruitoutLetrozoleDisc),
                              data.table(Cohort="Discovery",Treatment="Combination ribociclib" ,Function="Activation" , activeoutRiboDisc),
                              data.table(Cohort="Discovery",Treatment="Letrozole alone" ,Function="Activation" , activeoutLetrozoleDisc) )
#write.csv(totalstatsDiscovery,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/temporary/DiscoveryStatsMyeloidTcellComm.csv")





#######
ActivComms <- c("IL15_IL15RA", "IL15_IL2RA","IL15_IL2RB","IL18_IL18R1","IL18_IL18RAP")
resultoutActivDiscRibo <- rbindlist(lapply( 1:length(ActivComms),function(x){
  print (x)
  res<-tryCatch({
    m1 <- lmer( log(1+scaleTransduction) ~ log(1+Day) * dynamic_class3 +(1|Patient.Study.ID), data=   subsetInflam[ARM!="A"][Pair.Name== ActivComms[x] ])
    data.table( Pair.Name=ActivComms[x],  data.table(summary(m1)$coef,keep.rownames=T) )
  },
  error=function(x){
    return(NULL)
  })
  setnames(res,old=c("t value","Pr(>|t|)"),new=c("tval","pval"))
  res
}))
resultoutActivDiscRibo$adj.pval <- p.adjust(p=resultoutActivDiscRibo$pval,method="fdr")#,n=length(unique(resultout$Pair.Name)))
resultoutActivDiscRibo[rn=="log(1 + Day):dynamic_class3Response"]
#######

## Validation cohort communication between cell types 
TcellSubtypeCommCCI <- data.table(TcellSubtypeCommCCI_all [Cohort== "Validation"]%>%select(-Cohort))


# incdicator of notable macrophage enrichement and ERRB pathway activation indicative of apoptosis induction
TcellSubtypeCommCCI[,MacroEnriched :=F]
TcellSubtypeCommCCI[Patient.Study.ID=="001-167",MacroEnriched :=T]

subsetInflam <- TcellSubtypeCommCCI[ARM!="A"][Pair.Name %in% LRpairsFiltered[
  HPMR.Receptor%in% cytokineReceptors
  ]$Pair.Name][LigandPhenoCelltype=="Macrophages"]

subsetInflam[,startval:=sum((Day==0)*log(1+scaleTransduction))/sum(Day==0), by=c("Pair.Name","Patient.Study.ID","ReceptorPhenoCelltype")]
subsetInflam[ ,startvalAv:=median(startval,na.rm=T), by=c("Pair.Name","dynamic_class3","ReceptorPhenoCelltype")]



pairstotest <- unique(subsetInflam$Pair.Name)

resultout <- rbindlist(lapply( 1:length(pairstotest),function(x){
  print (x)
  res<-tryCatch({
    m1 <- lmer( log(1+scaleTransduction) ~ log(1+Day) * dynamic_class3 +(1|Patient.Study.ID), data=   subsetInflam[ARM!="A"][Pair.Name== pairstotest[x] ])
    data.table( Pair.Name=pairstotest[x],  data.table(summary(m1)$coef,keep.rownames=T) )
  },
  error=function(x){
    return(NULL)
  })
  res
}))
setnames(resultout,old=c("t value","Pr(>|t|)"),new=c("tval","pval"))
resultout$adj.pval <- p.adjust(p=resultout$pval,method="fdr")#,n=length(unique(resultout$Pair.Name)))



plotthis <- TcellSubtypeCommCCI[][Pair.Name %in%c("IL15_IL15RA","IL15_IL2RA","IL15_IL2RB","IL18_IL18R1" ,"IL18_IL18RAP")    ][LigandPhenoCelltype=="Macrophages"][MacroEnriched==F]
plotthis[,Tumorresponse:= "Resistant"]
plotthis[dynamic_class3=="Response",Tumorresponse:= "Sensitive"]

plotthis[, Pair.Name2 :=gsub("_","-",Pair.Name) ]
plotthis[, Function := "Activation" ]
plotthis$Pair.Name2 <- factor(plotthis$Pair.Name2 , levels= c("IL15-IL15RA","IL15-IL2RA","IL15-IL2RB","IL18-IL18R1" ,"IL18-IL18RAP") )

plotthis[,Treatmentlab:= "Combination ribociclib"]
plotthis[ARM=="A",Treatmentlab:= "Letrozole alone"]

plotthis[, y:=  scale( log(1 + scaleTransduction)),by=Pair.Name2]
plotthis[,y0:=sum(y*(Day==0),na.rm=T)/sum(Day==0) , by=c("Patient.Study.ID","Pair.Name2")]
plotthis[,deltay:=y-y0 , by=c("Patient.Study.ID","Pair.Name2")]
# 

validdat <-plotthis[Function=="Activation"][Day==180]
#save(validdat,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/temporary/ValidationMyeloidTcellActivation.RData")
load(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/temporary/ValidationMyeloidTcellActivation.RData")
Activationvaliddat <- data.table(Cohort="Validation", validdat )

ggplot( plotthis[Function=="Activation"][Day==180][],#%>%group_by(Pair.Name2)%>%mutate(y=scale( log(1 + scaleTransduction))) , 
        aes(y= deltay#y
            , Tumorresponse,group=interaction(Tumorresponse,Pair.Name2),col= Pair.Name2,fill= Pair.Name2)) +
  theme_classic(base_size=28) +
  stat_boxplot(geom ='errorbar',position=position_dodge(width=1),col="black") + 
  geom_boxplot(alpha = 0.6,col="black",position=position_dodge(width=1))+
  #geom_smooth(alpha=1,col= "black",method="gam",formula=y~s(x,k=3),se=T) +
  theme(aspect.ratio=1)+
  geom_point(position=position_dodge(width=1))+ facet_wrap(~Treatmentlab, ncol=1)+
  geom_vline(xintercept = 1.5,linetype="dashed") +
  #scale_color_discrete(name="Communication \n pathway")+
  scale_colour_brewer(name="Communication \n pathway", palette = "YlGn") +
  scale_fill_brewer(name="Communication \n pathway", palette = "YlGn") +
  labs(y="Myeloid to CD8+ T cell \n activation communications post treatment",x="Tumor response")

paperfile<- "/Users/jason/Dropbox/FELINE Project (1)/Manuscript  Feline immune communication/Figures Communication Project/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Validation IL communicaiton Ribo and letrozole Myeloid to T cell post treatment activation communications Validation C2.png"),height=10,width=10)



plotthis <- TcellSubtypeCommCCI[][Pair.Name %in%c("ADAM12_ITGB1","COL18A1_ITGB1","F13A1_ITGA4","F13A1_ITGB1","ICAM2_ITGB2","ICAM3_ITGB2","THBS2_ITGA4","VCAM1_ITGA4","VCAM1_ITGB1")    ][LigandPhenoCelltype=="Macrophages"][MacroEnriched==F]
plotthis[,"Tumorresponse":= "Resistant"]
plotthis[dynamic_class3=="Response","Tumorresponse":= "Sensitive"]

plotthis[, Pair.Name2 :=gsub("_","-",Pair.Name) ]
plotthis[, Function := "Recruitment" ]
plotthis$Pair.Name2 <- factor(plotthis$Pair.Name2 , levels= c("ADAM12-ITGB1","COL18A1-ITGB1","F13A1-ITGA4","F13A1-ITGB1","ICAM2-ITGB2","ICAM3-ITGB2","THBS2-ITGA4","VCAM1-ITGA4","VCAM1-ITGB1") )

plotthis[,Treatmentlab:= "Combination ribociclib"]
plotthis[ARM=="A",Treatmentlab:= "Letrozole alone"]

plotthis[, y:=  scale( log(1 + scaleTransduction)),by=Pair.Name2]
plotthis[,y0:=sum(y*(Day==0))/sum(Day==0) , by=c("Patient.Study.ID","Pair.Name2")]
plotthis[,deltay:=y-y0 , by=c("Patient.Study.ID","Pair.Name2")]
# 
validdat <-plotthis[][Day==180]#[!(Patient.Study.ID=="001-167")]
#save(validdat,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/temporary/ValidationMyeloidTcellRecruitment.RData")
load(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/temporary/ValidationMyeloidTcellRecruitment.RData")
Recruitmentvaliddat <- data.table(Cohort="Validation", validdat )

ggplot( plotthis[][Day==180],#[!(Patient.Study.ID=="001-167")],#%>%group_by(Pair.Name2)%>%mutate(y=scale( log(1 + scaleTransduction))) , 
        aes(y= deltay#y
            , Tumorresponse,group=interaction(Tumorresponse,Pair.Name2),col= Pair.Name2,fill= Pair.Name2)) +
  theme_classic(base_size=28) +
  stat_boxplot(geom ='errorbar',position=position_dodge(width=1),col="black") + 
  geom_boxplot(alpha = 0.6,col="black",position=position_dodge(width=1))+
  #geom_smooth(alpha=1,col= "black",method="gam",formula=y~s(x,k=3),se=T) +
  theme(aspect.ratio=1)+
  geom_point(position=position_dodge(width=1))+ facet_wrap(~Treatmentlab, ncol=1)+
  geom_vline(xintercept = 1.5,linetype="dashed") +
  #scale_color_discrete(name="Communication \n pathway")+
  scale_colour_brewer(name="Communication \n pathway", palette = "YlGn") +
  scale_fill_brewer(name="Communication \n pathway", palette = "YlGn") +
  labs(y="Myeloid to CD8+ T cell \n activation communications post treatment",x="Tumor response")


# Stats
pairstotest <- c("IL15_IL2RA","IL15_IL2RB","IL18_IL18R1","IL18_IL18RAP", "IL15_IL15RA" )

activeoutRiboValid <- rbindlist(lapply( 1:length(pairstotest),function(x){
  print (x)
  res<-tryCatch({
    m1 <- lmer( log(1+scaleTransduction) ~ log(1+Day) * dynamic_class3 +(1|Patient.Study.ID), data=   subsetInflam[Patient.Study.ID!="001-167"][ARM!="A"][Pair.Name== pairstotest[x] ])
    data.table( Pair.Name=pairstotest[x],  data.table(summary(m1)$coef,keep.rownames=T) )
  },
  error=function(x){
    return(NULL)
  })
  setnames(res,old=c("t value","Pr(>|t|)"),new=c("tval","pval"))
  res
}))[rn=="log(1 + Day):dynamic_class3Response"]

activeoutLetrozoleValid <- rbindlist(lapply( 1:length(pairstotest),function(x){
  print (x)
  res<-tryCatch({
    m1 <- lmer( log(1+scaleTransduction) ~ log(1+Day) * dynamic_class3 +(1|Patient.Study.ID), data=   subsetInflam[Patient.Study.ID!="001-167"][ARM=="A"][Pair.Name== pairstotest[x] ])
    data.table( Pair.Name=pairstotest[x],  data.table(summary(m1)$coef,keep.rownames=T) )
  },
  error=function(x){
    return(NULL)
  })
  setnames(res,old=c("t value","Pr(>|t|)"),new=c("tval","pval"))
  res
}))[rn=="log(1 + Day):dynamic_class3Response"]



pairstotest <- c("ADAM12_ITGB1" , "COL18A1_ITGB1", "F13A1_ITGA4","F13A1_ITGB1","ICAM2_ITGB2"  , "ICAM3_ITGB2","VCAM1_ITGA4","VCAM1_ITGB1","THBS2_ITGA4")
recruitoutRiboValid <- rbindlist(lapply( 1:length(pairstotest),function(x){
  print (x)
  res<-tryCatch({
    m1 <- lmer( log(1+scaleTransduction) ~ log(1+Day) * dynamic_class3 +(1|Patient.Study.ID), data=   subsetInflam[Patient.Study.ID!="001-167"][ARM!="A"][Pair.Name== pairstotest[x] ])
    data.table( Pair.Name=pairstotest[x],  data.table(summary(m1)$coef,keep.rownames=T) )
  },
  error=function(x){
    return(NULL)
  })
  setnames(res,old=c("t value","Pr(>|t|)"),new=c("tval","pval"))
  res
}))[rn=="log(1 + Day):dynamic_class3Response"]

recruitoutLetrozoleValid<- rbindlist(lapply( 1:length(pairstotest),function(x){
  print (x)
  res<-tryCatch({
    m1 <- lmer( log(1+scaleTransduction) ~ log(1+Day) * dynamic_class3 +(1|Patient.Study.ID), data=   subsetInflam[Patient.Study.ID!="001-167"][ARM=="A"][Pair.Name== pairstotest[x] ])
    data.table( Pair.Name=pairstotest[x],  data.table(summary(m1)$coef,keep.rownames=T) )
  },
  error=function(x){
    return(NULL)
  })
  setnames(res,old=c("t value","Pr(>|t|)"),new=c("tval","pval"))
  res
}))[rn=="log(1 + Day):dynamic_class3Response"]

recruitoutRiboValid[pval<0.05]
recruitoutLetrozoleValid[pval<0.05]
recruitoutLetrozoleValid[pval>0.05]
activeoutRiboValid[pval<0.05]
activeoutLetrozoleValid[pval<0.05]

totalstatsValidation <- rbind( data.table(Cohort="Validation",Treatment="Combination ribociclib" ,Function="Recruitment" , recruitoutRiboValid),
                               data.table(Cohort="Validation",Treatment="Letrozole alone" ,Function="Recruitment" , recruitoutLetrozoleValid),
                               data.table(Cohort="Validation",Treatment="Combination ribociclib" ,Function="Activation" , activeoutRiboValid),
                               data.table(Cohort="Validation",Treatment="Letrozole alone" ,Function="Activation" , activeoutLetrozoleValid) )
write.csv(totalstatsValidation,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/temporary/ValidationStatsMyeloidTcellComm.csv")

rm(list=ls())
totalstatsValidation<- data.table(read.csv(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/temporary/ValidationStatsMyeloidTcellComm.csv"))
totalstatsDiscovery<- data.table( read.csv(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/temporary/DiscoveryStatsMyeloidTcellComm.csv"))
allstats<-rbind(totalstatsDiscovery,totalstatsValidation)
allstats[,X:=NULL]
write.csv(allstats,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/Figure5/DiscoveryandValidationStatsMyeloidTcellComm.csv")





######## Combine all results

Activationalldat <- rbind(
  Activationdiscdat%>%dplyr::mutate(Communication=y)%>%dplyr::select(Function,Cohort,Patient.Study.ID, Day, dynamic_class3, ARM,    Pair.Name,Pair.Name2,key_,LigandPhenoCelltype,ReceptorPhenoCelltype,Communication,Treatmentlab,Tumorresponse),
  Activationvaliddat%>%dplyr::mutate(Communication=deltay)%>%dplyr::select(Function,Cohort,Patient.Study.ID, Day, dynamic_class3, ARM,    Pair.Name,key_,Pair.Name2,LigandPhenoCelltype,ReceptorPhenoCelltype,Communication,Treatmentlab,Tumorresponse)
)[!is.na(Communication)]
Activationalldat[, Daylab:=paste0("Day ",Day)]

#save(Activationalldat, file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/Figure5/DiscoveryandValidationMyeloidTcellActivation.RData")
savloc<- "/Users/jason/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/Figure5/"
write.csv(Activationalldat, 
          file=paste0(savloc,"MyeloidtoTcellActivation_Output.csv"))

ggplot( Activationalldat , aes(y= Communication, Tumorresponse,group=interaction(Tumorresponse,Pair.Name2),col= Pair.Name2,fill= Pair.Name2)) +
  theme_classic(base_size=26) +
  stat_boxplot(geom ='errorbar',position=position_dodge(width=1),col="black") + 
  geom_boxplot(alpha = 0.6,col="black",position=position_dodge(width=1))+
  #geom_smooth(alpha=1,col= "black",method="gam",formula=y~s(x,k=3),se=T) +
  theme(aspect.ratio=1)+
  geom_point(position=position_dodge(width=1))+ facet_grid(Treatmentlab~Cohort)+
  geom_vline(xintercept = 1.5,linetype="dashed") +
  #scale_color_discrete(name="Communication \n pathway")+
  scale_colour_brewer(name="Communication \n pathway", palette = "YlGn") +
  scale_fill_brewer(name="Communication \n pathway", palette = "YlGn") +
  labs(y="Myeloid to T cell \n activation communications post treatment",x="Tumor response")

ggplot( Activationalldat , aes(y= Communication, x=Tumorresponse, group=interaction(Tumorresponse,Pair.Name2), col= Pair.Name2, fill= Pair.Name2)) +
  theme_classic(base_size=28) +
  stat_boxplot(geom ='errorbar',position=position_dodge(width=1),col="black") + 
  geom_boxplot(alpha = 0.6,col="black",position=position_dodge(width=1))+
  #geom_smooth(alpha=1,col= "black",method="gam",formula=y~s(x,k=3),se=T) +
  theme(aspect.ratio=1)+
  geom_point(position=position_dodge(width=1))+ facet_wrap(Treatmentlab~Cohort, ncol=2)+
  geom_vline(xintercept = 1.5,linetype="dashed") +
  #scale_color_discrete(name="Communication \n pathway")+
  scale_colour_brewer(name="Communication \n pathway", palette = "YlGn") +
  scale_fill_brewer(name="Communication \n pathway", palette = "YlGn") +
  labs(y="Myeloid to T cell \n activation communication post treatment",x="Tumor response") +
  theme(axis.title=element_blank(),  axis.text=element_blank(),strip.text = element_blank(),legend.text=element_blank(),legend.title=element_blank())

Activationalldat[,Cohort2:=" "]
Activationalldat[Cohort=="Validation",Cohort2:="  "]
ggplot( Activationalldat , aes(y= Communication, Tumorresponse,group=interaction(Tumorresponse,Pair.Name2),col= Pair.Name2,fill= Pair.Name2)) +
  theme_classic(base_size=26) +
  stat_boxplot(geom ='errorbar',position=position_dodge(width=1),col="black") + 
  geom_boxplot(alpha = 0.6,col="black",position=position_dodge(width=1))+
  theme(aspect.ratio=1)+
  geom_point(position=position_dodge(width=1))+ facet_wrap(~paste0(Treatmentlab,Cohort2))+
  geom_vline(xintercept = 1.5,linetype="dashed") +
  scale_colour_brewer(name="Communication \n pathway", palette = "YlGn") +
  scale_fill_brewer(name="Communication \n pathway", palette = "YlGn") +
  labs(y="Myeloid to T cell \n activation communications post treatment",x="Tumor response")




Recruitmentalldat <- rbind(
  Recruitmentdiscdat%>%dplyr::mutate(Communication=y)%>%dplyr::select(Function,Cohort,Patient.Study.ID, Day, dynamic_class3, ARM,    Pair.Name,Pair.Name2,key_,LigandPhenoCelltype,ReceptorPhenoCelltype,Communication,Treatmentlab,Tumorresponse),
  Recruitmentvaliddat%>%dplyr::mutate(Communication=deltay)%>%dplyr::select(Function,Cohort,Patient.Study.ID, Day, dynamic_class3, ARM,    Pair.Name,key_,Pair.Name2,LigandPhenoCelltype,ReceptorPhenoCelltype,Communication,Treatmentlab,Tumorresponse)
)[!is.na(Communication)]
Recruitmentalldat[, Daylab:=paste0("Day ",Day)]
#save(Recruitmentalldat, file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/Figure5/DiscoveryandValidationMyeloidTcellRecruitment.RData")
savloc<- "/Users/jason/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/Figure5/Outputs/"
write.csv(Recruitmentalldat, file=paste0(savloc,"MyeloidtoTcellRecruitment_Output.csv"))

ggplot( Recruitmentalldat , aes(y= Communication, Tumorresponse,group=interaction(Tumorresponse,Pair.Name2),col= Pair.Name2,fill= Pair.Name2)) +
  theme_classic(base_size=26) +
  stat_boxplot(geom ='errorbar',position=position_dodge(width=1),col="black") + 
  geom_boxplot(alpha = 0.6,col="black",position=position_dodge(width=1))+
  theme(aspect.ratio=1)+
  geom_point(position=position_dodge(width=1))+ facet_grid(Treatmentlab~Cohort)+
  geom_vline(xintercept = 1.5,linetype="dashed") +
  scale_colour_brewer(name="Communication \n pathway", palette = "YlGn") +
  scale_fill_brewer(name="Communication \n pathway", palette = "YlGn") +
  labs(y="Myeloid to T cell \n recruitment communications post treatment",x="Tumor response")

ggplot( Recruitmentalldat , aes(y= Communication, x=Tumorresponse, group=interaction(Tumorresponse,Pair.Name2), col= Pair.Name2, fill= Pair.Name2)) +
  theme_classic(base_size=28) +
  stat_boxplot(geom ='errorbar',position=position_dodge(width=1),col="black") + 
  geom_boxplot(alpha = 0.6,col="black",position=position_dodge(width=1))+
  theme(aspect.ratio=1)+
  geom_point(position=position_dodge(width=1))+ facet_wrap(Treatmentlab~Cohort, ncol=2)+
  geom_vline(xintercept = 1.5,linetype="dashed") +
  scale_colour_brewer(name="Communication \n pathway", palette = "YlGn") +
  scale_fill_brewer(name="Communication \n pathway", palette = "YlGn") +
  labs(y="Myeloid to T cell \n recruitment communication post treatment",x="Tumor response") +
  theme(axis.title=element_blank(),  axis.text=element_blank(),strip.text = element_blank(),legend.text=element_blank(),legend.title=element_blank())

Recruitmentalldat[,Cohort2:=" "]
Recruitmentalldat[Cohort=="Validation",Cohort2:="  "]
ggplot( Recruitmentalldat , aes(y= Communication, Tumorresponse,group=interaction(Tumorresponse,Pair.Name2),col= Pair.Name2,fill= Pair.Name2)) +
  theme_classic(base_size=26) +
  stat_boxplot(geom ='errorbar',position=position_dodge(width=1),col="black") + 
  geom_boxplot(alpha = 0.6,col="black",position=position_dodge(width=1))+
  theme(aspect.ratio=1)+
  geom_point(position=position_dodge(width=1))+ facet_wrap(~paste0(Treatmentlab,Cohort2))+
  geom_vline(xintercept = 1.5,linetype="dashed") +
  scale_colour_brewer(name="Communication \n pathway", palette = "YlGn") +
  scale_fill_brewer(name="Communication \n pathway", palette = "YlGn") +
  labs(y="Myeloid to T cell \n recruitment communications post treatment",x="Tumor response")
