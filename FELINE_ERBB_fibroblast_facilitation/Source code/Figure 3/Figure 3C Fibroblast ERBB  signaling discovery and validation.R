rm(list=ls())
require(mgcv);require(data.table);require(dplyr);require(ggplot2);require(tidyr);require(lme4);require(lmerTest);require(parallel)
library(effects);require(umap)
require(Rfast);require(ider)
library("dendextend");library(ggdendro);require(ggsci);require(viridis)
require("Rdimtools")

# Define data location
SourceDataLoc <- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript Feline ERBB facilitation/Molecular Systems Biology/Revision and submission folder/Source Data/"
Intermediateloc <- paste0(SourceDataLoc,"Figure 3/")

Discovery_FibroExpression <- data.table(read.csv( file= paste0(Intermediateloc,"SourceData_Figure_3C_DiscoveryFibroblastERBBLigands.csv" ))%>%
                                          gather(var,val,NRG1:ADAM17))
#Discretize to 10 fibroblast mesenchymal differentiation levels
Discovery_FibroExpression[,emtlevel:= round(-V3*1.15)/1.15]#length(unique(Discovery_FibroExpression$emtlevel))
Discovery_FibroExpression[,egfproliflevel:= round(V1*1.925)/1.925]#length(unique(Discovery_FibroExpression$egfproliflevel))

#Aggregate ligand expression across genes for each cell (assume equal weights)
Discovery_FibroExpression[,scaleval:=scale( log(1+val) ) ,by=var]
Discovery_FibroSignalInput <- data.table( Discovery_FibroExpression%>%group_by(Cell.ID,Timepoint,Day,
                                                           Patient.Study.ID)%>%dplyr::summarise(emtlevel=emtlevel[1], egfproliflevel= egfproliflevel[1], lnvalmu= mean( scaleval )  )  )
#ggplot(Discovery_FibroExpression[][var=="NRG1"],aes(y= scaleval ,x=emtlevel))+geom_boxplot(aes(group=emtlevel)) +geom_smooth(method="gam", formula= y~s(x,k=3))+facet_wrap(~Timepoint,scales="free_y")
#ggplot(Discovery_FibroExpression[val>0][var=="HBEGF"],aes(y= scaleval ,x=emtlevel))+geom_point() +geom_smooth(method="gam", formula= y~s(x,k=4))+facet_wrap(~var,scales="free_y")

# Annotate differentiation states and cohort 
Discovery_FibroSignalInput[,egfproliflevelB:="Quiescent"]
Discovery_FibroSignalInput[egfproliflevel==T,egfproliflevelB:="EGF proliferative"]
Discovery_FibroSignalInput$egfproliflevelB <- factor(Discovery_FibroSignalInput$egfproliflevelB , levels= c("Quiescent","EGF proliferative"  ) )
Discovery_catLU<-data.table(emtlevel=sort( unique(Discovery_FibroSignalInput$emtlevel)) , EMTlevel= formatC(  1:length(unique(Discovery_FibroSignalInput$emtlevel))   , width = 2, format = "d", flag = "0")  )
Discovery_FibroSignalInput <- merge(Discovery_FibroSignalInput,Discovery_catLU, by="emtlevel")
Discovery_FibroSignalInput[,Cohort:="Discovery cohort"]



## Repeat for validation cohort
Validation_FibroExpression <- data.table(read.csv( file= paste0(Intermediateloc,"SourceData_Figure_3C_ValidationFibroblastERBBLigands.csv" ))%>%
                                          gather(var,val,ANXA1:TGFA))

#Discretize to 10 fibroblast mesenchymal differentiation levels
Validation_FibroExpression[,emtlevel:= round(V3*1.84)/1.84]
unique(Validation_FibroExpression$emtlevel)%>%length()
Validation_FibroExpression[,egfproliflevel:= round(V1*1.925)/1.925]

#Aggregate ligand expression across genes for each cell (assume equal weights)
Validation_FibroExpression[,scaleval:=scale( log(1+val) ) ,by=var]
Validation_FibroSignalInput <- data.table( Validation_FibroExpression%>%group_by(Cell.ID,Timepoint,Day,
                                                                               Patient.Study.ID)%>%dplyr::summarise(emtlevel=emtlevel[1], egfproliflevel= egfproliflevel[1], lnvalmu= mean( scaleval )  )  )
# Annotate differentiation states and cohort 
Validation_FibroSignalInput[,egfproliflevelB:="Quiescent"]
Validation_FibroSignalInput[egfproliflevel==T,egfproliflevelB:="EGF proliferative"]
Validation_FibroSignalInput$egfproliflevelB <- factor(Validation_FibroSignalInput$egfproliflevelB , levels= c("Quiescent","EGF proliferative"  ) )
Validation_catLU<-data.table(emtlevel=sort( unique(Validation_FibroSignalInput$emtlevel)) , EMTlevel= formatC(  1:length(unique(Validation_FibroSignalInput$emtlevel))   , width = 2, format = "d", flag = "0")  )
Validation_FibroSignalInput <- merge(Validation_FibroSignalInput,Validation_catLU, by="emtlevel")
Validation_FibroSignalInput[,Cohort:="Validation cohort"]

nrow(Discovery_FibroExpression)+nrow(Validation_FibroExpression)

# Aggregate by fibroblast subtype within each sample
plotallSummar <- data.table(rbind(Discovery_FibroSignalInput,Validation_FibroSignalInput)%>%
                              group_by(emtlevel,EMTlevel,#PlotID,
                                       Cohort,Timepoint,Day,
                                       Patient.Study.ID)%>%
                              dplyr::summarise(lnvalmu1=mean(lnvalmu),
                                               lnvalmu2=median(lnvalmu)))

ggplot(plotallSummar,aes(y= lnvalmu1 ,x=(as.numeric(EMTlevel)) ))+
  theme_classic(base_size=26)+theme(aspect.ratio=1)+ # or aspect.ratio=0.71
  labs(x="Fibroblast myCAF differentiation level",y="ERBB ligand signaling (scaled)")+
  geom_boxplot(aes(group=interaction(Cohort,EMTlevel)    ,fill= as.numeric(EMTlevel) ),position = position_dodge(.7),outlier.color=NA) +
  stat_boxplot(aes(group=interaction(Cohort,EMTlevel)    ,fill= as.numeric(EMTlevel),col=Cohort),geom = "errorbar", position = position_dodge(.7),width =0.75,size=2) +
  geom_point(data=plotallSummar[],aes(group=interaction(Cohort,EMTlevel)  ,col=Cohort  ,fill= as.numeric(EMTlevel) ),pch=21,position=position_dodge(.7),size=2.5) +
  scale_fill_gradient(guide=F,low =pal_aaas()(6)[3] ,high =pal_aaas()(6)[1]) + 
  scale_color_manual(values=c("slategray","grey")) + 
  theme(legend.position = "none")+
  #lims(y=c(-0.5,2))+
  scale_x_continuous(breaks = 1:10)

# Aggegate within highly and lowly differentiated cell populations within each tumor
plotallSummar[,EMTlevel2:=(as.numeric(EMTlevel)>7)]
plotthisSummary <- data.table(plotallSummar %>% group_by(EMTlevel2, Cohort, Day, Patient.Study.ID )%>%summarise(lnvalmu1=median(lnvalmu1)))

# Statistical analysis
summary(lm(lnvalmu1~EMTlevel2,plotthisSummary[Cohort=="Discovery cohort"]))
summary(lm(lnvalmu1~EMTlevel2,plotthisSummary[Cohort=="Validation cohort"]))

# Frequency of detection of fibroblasts high/low mesenchymal subtypes within tumor samples in the discover and validation cohort
plotthisSummary[Cohort=="Discovery cohort"]%>%group_by(EMTlevel2)%>%summarise(n())
plotthisSummary[Cohort=="Validation cohort"]%>%group_by(EMTlevel2)%>%summarise(n())

plotthisSummary[Cohort=="Discovery cohort"]%>%group_by(EMTlevel2)%>%summarise(length(unique(Patient.Study.ID)))
plotthisSummary[Cohort=="Validation cohort"]%>%group_by(EMTlevel2)%>%summarise(length(unique(Patient.Study.ID)))

plotthisSummary[Cohort=="Discovery cohort"][EMTlevel2==T]%>%select(Patient.Study.ID,Day)
plotthisSummary[Cohort=="Discovery cohort"][EMTlevel2==F]%>%select(Patient.Study.ID,Day)
plotthisSummary[Cohort=="Validation cohort"][EMTlevel2==T]%>%select(Patient.Study.ID,Day)
plotthisSummary[Cohort=="Validation cohort"][EMTlevel2==F]%>%select(Patient.Study.ID,Day)
plotthisSummary[Cohort=="Discovery cohort"]%>%select(Patient.Study.ID,Day)
plotthisSummary[Cohort=="Validation cohort"] %>%select(Patient.Study.ID,Day)

# Visualize
ggplot(plotthisSummary,aes(y= lnvalmu1 , x=Cohort ))+
  theme_classic(base_size=26)+theme(aspect.ratio=1)+ # or aspect.ratio=0.71
  labs(x="Cohort",y="ERBB ligand signaling (scaled)")+
  geom_boxplot(aes(group=interaction(Cohort,EMTlevel2)    ,fill= as.numeric(EMTlevel2) ),position = position_dodge(.7),width=0.5,outlier.color=NA) +
  stat_boxplot(aes(group=interaction(Cohort,EMTlevel2)    ,fill= as.numeric(EMTlevel2)),geom = "errorbar", position = position_dodge(.7),width =0.75) +
  geom_point(data=plotthisSummary[],aes(group=interaction(Cohort,EMTlevel2)  ,fill= as.numeric(EMTlevel2) ), pch=21,position=position_dodge(.7),size=2.5) +
  scale_fill_gradient(guide=F,low =pal_aaas()(6)[3] ,high =pal_aaas()(6)[1]) + 
  #scale_color_manual(values=c("slategray","grey")) + 
  #theme(legend.position = "none")+
  #lims(y=c(-0.5,2))+
  scale_x_discrete(labels=c("Discovery","Validation"))
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Discovery and Validation Fibroblast differentiation increases ERBBligandexpression of fibrobastsRevised.pdf",height=6,width=6,dpi = 320)

