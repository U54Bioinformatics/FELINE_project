rm(list=ls()); require(data.table);require(dplyr);require(ggplot2);require(tidyr);
library(plater)
library(scales)
require(lmerTest)


fileLoc<-"/Users/jason/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Fibroblast_Experiments/coculture/20230223_Feng Afatinib Fulvestrant CAMA1 MCF7 T47D/RawDataQC/"
fileNms <-paste0(fileLoc,list.files(fileLoc)[grepl(".csv",list.files(fileLoc))])

# Read Multiple Plates
dd <- data.table(read_plates(
  files =  fileNms[],
  plate_names = c("Plate_1", "Plate_2", "Plate_3", "Plate_4", "Plate_5", "Plate_6"),
  well_ids_column = "Well",   sep = "," 
))

names(dd)<-gsub("YFP i" ,"Sens_I",names(dd))
names(dd)<-gsub("CFP i" ,"Res_I",names(dd))
names(dd)<-gsub("red i" ,"Fibro_I",names(dd))
dd$PercentFBS <-1

dd[, Resistance := "Sensitive"]
dd[!grep("Sens", CancerCell), Resistance := "RiboR"]

#Rep identifiers
dd[,RepID:=1:nrow(dd)]
dd[,Rep:=1:length(RepID), by=c("CancerCell","Resistance","Fibroblasts","PercentFBS","Fulvestrant_nM","Afatinib_uM")]

# First we'll continue to work with florescence measured total intensity; filtering ATP conten
dat2 <- dd %>% dplyr::select(c("RepID","Rep","Resistance","CancerCell":"Fulvestrant_nM","Sens_Intensity_24hr":"Sens_Intensity_336hr","Fibro_Intensity_24hr":"Fibro_Intensity_336hr", "Res_Intensity_24hr":"Res_Intensity_336hr"))


# Convert to Long format 
dat_long <- data.table(gather(dat2, Time, Measurement, 'Sens_Intensity_24hr':'Res_Intensity_336hr'))
#dat_long$CancerCell%>%unique()

# Split one column with 3 pieces of information into 3 different columns
dat_long[, c("CellType", "  ","Timehr") := tstrsplit(Time, "_", fixed=TRUE)]
dat_long[,"  ":=NULL]

# get numeric time
dat_long[, Hour := as.numeric( gsub("[^0-9]","",Time) )]
dat_long[,Timehr:=NULL]
dat_long[,Time:=NULL]

# adjust labeling to match Fengs standard curve data
dat_long[, Variable:=paste(gsub(" ","_",CellType),"Intensity",sep="_") ]
setnames(dat_long,old=c("Measurement"),new=c("Value"))
inputdd <- data.table(dat_long%>%dplyr::select(-c("CellType"))%>%spread(Variable,Value))


# load regression model fits obtained from Feng's standard curve mapping data (translate intensity to counts )
getCountModels <- function(){
  # Load the data
  dd <- data.table( read.csv(file="~/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Fibroblast_Experiments/coculture/221110_FengCellNumberIntensity/RawDataQC_PC/FengCellNumberIntensities_Parsed_Clean.csv") )
  # Additional columns to subset either by cell state (fibrobalst sensitive or resistant cancer) or cell type (lineage name)
  dd[, Celltype := "Fibroblast"]
  dd[grep("Sens_", CancerCell), Celltype := "Sens"]
  dd[grep("RiboR_",CancerCell), Celltype := "RiboR"]
  
  # Add column labelling different cancer cell types
  dd[, CelltypeName :=  "MRC5"]
  dd[grep("CAMA1",CancerCell), CelltypeName :=  "CAMA1"]
  dd[grep("T47D",CancerCell), CelltypeName :=  "T47D"]
  dd[grep("MCF7",CancerCell), CelltypeName :=  "MCF7"]
  
  # Alter intensity varaible names to match that of Feng's full dataset
  names(dd) <- gsub("_24hr","", names(dd))
  names(dd)[names(dd)=="Res_Intensity"] <-  "RiboR_Intensity"
  
  # Use linear regression to predict cell number from fluorescent intensity
  SensLm <- lm(Cell_Number~ -1 + (Sens_Intensity):CancerCell,   data = dd[Celltype =="Sens"] ) # plot(SensLm) # coef(SensLm)
  ResLm <- lm(Cell_Number~ -1 + (RiboR_Intensity):CancerCell,   data = dd[Celltype =="RiboR"] )#  plot(ResLm) # coef(ResLm)
  FibroLm <- lm(Cell_Number~-1 + (Fibro_Intensity),  data = dd[Celltype =="Fibroblast"] )# plot(FibroLm)# coef(FibroLm)
  return(list("SensLm"=SensLm,"ResLm"=ResLm,"FibroLm"=FibroLm))
}
CountModels <- getCountModels()

# Predict counts for each cell type
names(inputdd)[names(inputdd)=="Res_Intensity"] <-  "RiboR_Intensity"
inputdd[,Sens_Count:= 0]
inputdd[,Res_Count:= 0]
inputdd[,Sensitivity:=Resistance]
inputdd[Resistance=="Sensitive",Sensitivity:="Sens"]

inputdd[CancerCell=="CAMA1_riboR-CER2",CancerCell:="CAMA1_RiboR_Cer2"]
inputdd[CancerCell=="MCF7_riboR-CER2",CancerCell:="MCF7_RiboR_Cer2"]
inputdd[CancerCell=="T47D_riboR-CER2",CancerCell:="T47D_RiboR_Cer2"]
inputdd[Sensitivity=="Sens", Sens_Count:= predict(CountModels$SensLm,newdata=inputdd[Sensitivity=="Sens"])]
inputdd[Sensitivity=="RiboR", Res_Count:= predict(CountModels$ResLm,newdata=inputdd[Sensitivity=="RiboR"])]

# remove intensity now that we have counts
inputdd[,RiboR_Intensity :=NULL]
inputdd[,Sens_Intensity:=NULL]

# create a single total count column instead of one count column per cell type and remove other redundant data from object 
inputdd[,Count:=Sens_Count + Res_Count ]
inputdd[Ratio_CancerCelltoFibro==0, Count:=0]

# Estimate: Intensity per fibroblast cell assuming 1000 fibbroblasts at 1st timepoint
FibroIntenPerCell <- (inputdd[Ratio_CancerCelltoFibro==0][Hour==24][Afatinib_uM==0][Fulvestrant_nM==0]$Fibro_Intensity)%>%median()/1000
inputdd[,FibroCount:=Fibro_Intensity/FibroIntenPerCell ]
inputdd[,Sens_Count:=NULL]
inputdd[,Res_Count:=NULL]
inputdd[, "CancerCellName" := tstrsplit(CancerCell, "_", fixed=TRUE)[1]]

inputdd[,CancerCell:=NULL]
inputdd[Fibroblasts=="None"]$FibroCount%>%hist()
inputdd[Fibroblasts=="None",FibroCount:=0 ]
inputdd[,Fibroblasts:=gsub("-",".",Fibroblasts)]
inputdd[,FibroblastsPresent:=T]
inputdd[Fibroblasts=="None",FibroblastsPresent:=F]

setnames(inputdd, old="Ratio_CancerCelltoFibro", new="InitialCancerFraction")
inputdd[,CancerPresent:=T]
inputdd[InitialCancerFraction==0,CancerPresent:=F]
inputdd[,Afatinib_uMLab:=paste0("Afatinib: ",Afatinib_uM , " uM")]
inputdd[Afatinib_uM==0,Afatinib_uMLab:="Control"]
inputdd$Afatinib_uMLab  <- factor(inputdd$Afatinib_uMLab  ,  levels=c("Control","Afatinib: 2.5 uM" ) )
inputdd[,Fulvestrant_nMLab:=paste0("Fulvestrant: ",Fulvestrant_nM , " nM")]
inputdd[Fulvestrant_nM==0,Fulvestrant_nMLab:="Control"]
inputdd$Fulvestrant_nMLab  <- factor(inputdd$Fulvestrant_nMLab  ,  levels=c("Control","Fulvestrant: 5 nM" ) )


inputdd[,ResistanceLab:="Resistant"]
inputdd[Resistance=="Sensitive",ResistanceLab:="Sensitive"]



inputdd[, CompositionLab:="Fibroblasts alone"]
inputdd[InitialCancerFraction>0&InitialCancerFraction<1&Fibroblasts=="EMT.MRC5_C2", CompositionLab:="Cancer and primed fibroblast"]
inputdd[InitialCancerFraction>0&InitialCancerFraction<1&Fibroblasts=="MRC5_C2", CompositionLab:="Cancer and unprimed fibroblast"]
inputdd[InitialCancerFraction==1, CompositionLab:="Cancer alone"]
inputdd$CompositionLab <- factor(inputdd$CompositionLab , levels=c("Fibroblasts alone","Cancer alone","Cancer and unprimed fibroblast","Cancer and primed fibroblast"))
inputdd[,maxHour:=max(Hour), by="RepID"]
inputdd[, initialCount:= sum((Hour==24)*Count), by=RepID]
inputdd[, initialFibroCount:= sum((Hour==24)*FibroCount), by=RepID]


plotthis00 <-rbindlist(lapply(unique(inputdd[Fibroblasts!="EMT.MRC5_C2"][Fulvestrant_nM%in%c(0,5)][Afatinib_uM%in%c(0,2.5)]$RepID),function(i){
  cat(i) 
  tmp<- inputdd[RepID==i][Fibroblasts!="EMT.MRC5_C2"][Fulvestrant_nM%in%c(0,5)][Afatinib_uM%in%c(0,2.5)]
  tmp[Count==0,Count:=1]
  tmp[FibroCount==0,FibroCount:=1]
  if(tmp$CompositionLab[1]=="Fibroblasts alone"){
    RGR=NA
  }else{
    RGR = coef(summary(lm(log(Count)~ Hour, data=tmp)))["Hour","Estimate"]
  }  
  if(tmp$CompositionLab[1]=="Cancer alone"){
    RGRfibro=NA
  }else{
    RGRfibro= coef(summary(lm(log(FibroCount)~ Hour, data=tmp[is.finite(FibroCount)])))["Hour","Estimate"]
  }
  out<- data.table(unique(tmp%>%
                            select(RepID,Rep,Resistance,Fibroblasts,Cytation,Total_Cell_Number,InitialCancerFraction,Afatinib_uM,Fulvestrant_nM,Sensitivity,
                                CancerCellName,FibroblastsPresent,CancerPresent,Afatinib_uMLab,Fulvestrant_nMLab,ResistanceLab,CompositionLab,maxHour))
                   ,RGR=RGR,RGRfibro=RGRfibro)
}))


# Part 3 = Facilitation of cancer increases with fibroblast number

FfreqDat <- plotthis00[Afatinib_uM==0][Fibroblasts!="EMT.MRC5_C2"][][InitialCancerFraction>0]

FfreqDat[,RGRalone:=sum(RGR*(CompositionLab=="Cancer alone"))/sum(CompositionLab=="Cancer alone"),by=c("Fulvestrant_nM","CancerCellName","ResistanceLab")]
FfreqDat[, wrapLab:= "DMSO"]
FfreqDat[Afatinib_uM==2.5&Fulvestrant_nM==0, wrapLab:= "Afatinib (2.5uM)"]
FfreqDat[Afatinib_uM==0&Fulvestrant_nM==5, wrapLab:= "Fulvestrant (5nM)"]
FfreqDat[Afatinib_uM==2.5&Fulvestrant_nM==5, wrapLab:= "Combination"]
FfreqDat$wrapLab <- factor(FfreqDat$wrapLab , levels=c("DMSO","Fulvestrant (5nM)","Afatinib (2.5uM)","Combination"))

collabs<-(ggsci::pal_aaas()(6)[c(2,3,1)])

ggplot(FfreqDat, aes(y= RGR-RGRalone, x= 100*(1-InitialCancerFraction)))+
  geom_smooth(method="gam",col="black",fill="black", formula=y~s(x,k=3,bs="cr"),alpha=0.2,se=T)+
  stat_boxplot(aes(fill=CompositionLab,col=CompositionLab, group= interaction(InitialCancerFraction,CancerCellName,ResistanceLab,Fulvestrant_nM,CompositionLab)),geom="errorbar",position=position_dodge(width=8))+
  geom_boxplot(col="black",aes(fill=CompositionLab,col=CompositionLab, group= interaction(InitialCancerFraction,CancerCellName,ResistanceLab,Fulvestrant_nM,CompositionLab)),
               width=6,position=position_dodge(width=8))+
  geom_point(aes(fill=CompositionLab,col=CompositionLab),pch=21,col="black",
             position=position_dodge(width=8))+
  theme_classic(base_size=26)+
  facet_grid(wrapLab~.)+theme(aspect.ratio=1)+
  labs(y="Cancer growth rate \n (compared to monoculture)",x="Fibroblast coculture frequency \n(% of fixed cancer population)")+
  scale_color_manual(name="",values=collabs[-2])+
  scale_fill_manual(name="",values=collabs[-2])+
  scale_x_continuous(breaks=(100*unique(1-FfreqDat$InitialCancerFraction))  )+
  theme(aspect.ratio=1,legend.position = "bottom")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
  theme(legend.position="none")

#ggsave(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Fibroblast_Experiments/coculture/20230223_Feng Afatinib Fulvestrant CAMA1 MCF7 T47D/Modeling_JG/MUTUALISM FINALIZE part3 cancer growth rate improves under more frequent fibroblast coculture DMSO and fulv treatment.pdf", width=9, height=12, dpi=320)
ggsave("~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/MUTUALISM FINALIZE part3 cancer growth rate improves under more frequent fibroblast coculture DMSO and fulv treatment.pdf",height=9,width=12,dpi = 320)

# Stats
require(mgcv) 
require(lmerTest)
dmod<-gam( I(RGR-RGRalone)~ s(I(100*(1-InitialCancerFraction)),k=3),data=FfreqDat[wrapLab=="DMSO"])
fmod<-gam( I(RGR-RGRalone)~ s(I(100*(1-InitialCancerFraction)),k=3),data=FfreqDat[wrapLab!="DMSO"])

FfreqDat[,CancerCellNameSensitivity:= as.factor(paste0(CancerCellName,Sensitivity))]
dmodlmer<-lmer( I(RGR-RGRalone)~ I(100*(1-InitialCancerFraction))+
                  (-1+InitialCancerFraction|CancerCellNameSensitivity ),data=FfreqDat[wrapLab=="DMSO"])
summary(dmod)
summary(fmod)
###### Final stats
fmodlmer<-lmer( I(RGR-RGRalone)~ I(100*(1-InitialCancerFraction))+
                  (-1+InitialCancerFraction|CancerCellNameSensitivity ),data=FfreqDat[wrapLab!="DMSO"])

fmodlm<-lm( I(RGR-RGRalone)~ I(100*(1-InitialCancerFraction)) ,data=FfreqDat[wrapLab!="DMSO"])
summary(dmodlmer)
summary(fmodlmer)
###### 



# save( FfreqDat , file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/Figure5/Cancer growth increase with fibroblast freq DMSO and fulv.RData")
# write.csv( FfreqDat , file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/Figure5/Fig5c/Cancer growth increase with fibroblast freq DMSO and fulv.csv")











# Part 1/2 = Facilitation of cancer and fibroblasts
plotthis0 <- plotthis00[Fibroblasts!="EMT.MRC5_C2"][Fulvestrant_nM%in%c(0,5)][Afatinib_uM%in%c(0,2.5)][InitialCancerFraction%in%c(0,0.5,1)]
plotthis0[, wrapLab:= "DMSO"]
plotthis0[Afatinib_uM==2.5&Fulvestrant_nM==0, wrapLab:= "Afatinib (2.5uM)"]
plotthis0[Afatinib_uM==0&Fulvestrant_nM==5, wrapLab:= "Fulvestrant (5nM)"]
plotthis0[Afatinib_uM==2.5&Fulvestrant_nM==5, wrapLab:= "Combination"]
plotthis0$wrapLab <- factor(plotthis0$wrapLab , levels=c("DMSO","Fulvestrant (5nM)","Afatinib (2.5uM)","Combination"))

plotthis0[CompositionLab=="Fibroblasts alone"][wrapLab=="DMSO"]
plotthis0[is.infinite(RGRfibro),RGRfibro:=NA]
plotthis0[is.infinite(RGR),RGR:=NA]

plotthis0[Fulvestrant_nM==0][Afatinib_uM==0][CompositionLab=="Cancer alone"]$RGR%>%mean()
#Cancer alone 0.002480009
plotthis0[Fulvestrant_nM==0][Afatinib_uM==0][CompositionLab=="Cancer alone"][Sensitivity=="RiboR"]$RGR%>%mean()
#MONOCULTURE Resitant Cancer  0.002328489
plotthis0[Fulvestrant_nM==0][Afatinib_uM==0][CompositionLab=="Cancer and unprimed fibroblast"][Sensitivity=="RiboR"]$RGR%>%mean()
#COCULUTRE Resitant Cancer 0.002495303
##beta_C=(0.003303814/0.002328489 -1)/1000


plotthis0[Fulvestrant_nM==0][Afatinib_uM==0][CompositionLab=="Fibroblasts alone"][]$RGRfibro%>%max()
#MONOCULTURE Fibroblast  1e-5 -0.002825221 (set slightly positive)
plotthis0[Fulvestrant_nM==0][Afatinib_uM==0][CompositionLab=="Cancer and unprimed fibroblast"][]$RGRfibro%>%max()
#COCULUTRE Fibroblast  0.002077994 
##beta_M=(0.002077994/1e-5 -1)/1000


plotthis0[Fulvestrant_nM==5][Afatinib_uM==0][CompositionLab=="Cancer alone"][Sensitivity=="RiboR"]$RGR%>%mean()
# fulv treated cancer monoculture: 0.001138915
#X = (0.002328489/0.001138915) -1


plotthis <- data.table(plotthis0%>% group_by(wrapLab,ResistanceLab,Fulvestrant_nM,Afatinib_uM,CompositionLab,CancerCellName)%>%
                         dplyr::summarise(RGR=mean(RGR,na.rm=T),
                                          RGRfibro=mean(RGRfibro,na.rm=T)))

refRGR <- data.table( plotthis[CompositionLab=="Cancer alone"][] %>% select(wrapLab,ResistanceLab,CancerCellName,RGR))
refRGRfibro <- data.table( plotthis[CompositionLab=="Fibroblasts alone"][] %>% select(wrapLab,RGRfibro))
setnames(refRGR,old=c("RGR"), new=c("DMSOaloneRGR"))
setnames(refRGRfibro,old=c("RGRfibro"), new=c("DMSOaloneRGRfibro"))

plotthis1 <- merge(merge(plotthis0,refRGR,by=c("ResistanceLab", "CancerCellName","wrapLab"  )),
                   refRGRfibro%>%group_by(wrapLab)%>%summarise(DMSOaloneRGRfibro=mean(DMSOaloneRGRfibro)),
                   by="wrapLab")

collabs<-(ggsci::pal_aaas()(6)[c(2,3,1)])
plotthis1[, FibroComLab:=paste0( CancerCellName," \n ",ResistanceLab)]
plotthis1[, CancerComLab:=paste0( CancerCellName," \n ",ResistanceLab)]
plotthis1[CompositionLab=="Cancer alone", FibroComLab:="No fibroblast"]
plotthis1[CompositionLab=="Fibroblasts alone", FibroComLab:="Alone"]
plotthis1$FibroComLab<- factor(plotthis1$FibroComLab, 
                               levels=c("Alone","No fibroblast" ,
                                        "CAMA1 \n Resistant","CAMA1 \n Sensitive",
                                        "MCF7 \n Resistant" ,"MCF7 \n Sensitive",
                                        "T47D \n Resistant","T47D \n Sensitive"))



fibroblastplot <-ggplot(plotthis1[Afatinib_uM==0][FibroComLab!="No fibroblast"],aes(y= (RGRfibro-DMSOaloneRGRfibro),#/abs(DMSOaloneRGRfibro),
                                                                                    x= FibroComLab, 
                                                                                    group= interaction(FibroComLab),col=CompositionLab,fill=CompositionLab))+
  stat_boxplot(geom="errorbar")+
  geom_boxplot(col="black")+theme_classic()+
  geom_point(pch=21,col="black")+
  theme_classic(base_size=26)+theme(aspect.ratio = 0.3)+
  #coord_trans(y="log")+
  facet_grid( wrapLab~., scales="free_x") +
  scale_x_discrete(labels=c("Alone","CAMA1 (CRR)","CAMA1 (ETR)",
                                "MCF7 (CRR)" ,"MCF7 (ETR)",
                                "T47D (CRR)","T47D (ETR)") )+
  scale_color_manual(name="",labels=c("Fibroblast \n alone","Cancer and \n fibroblasts"),values=collabs[-1])+
  scale_fill_manual(name="",labels=c("Fibroblast \n alone","Cancer and \n fibroblasts"),values=collabs[-1])+
  labs(x="Fibroblast coculture", y="Fibroblast growth rate \n (compared to monoculture)")+
  #scale_x_discrete(labels=c("Cancer \n alone", "Cancer + \n unprimed fibroblasts", "Cancer + \n primed fibroblasts"))+
  theme(aspect.ratio=1,legend.position = "bottom")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1) ,
        axis.text.x = element_text(vjust=0.5,angle=90))+
  theme(legend.position = "none")

plotthis1$CancerComLab<- factor(plotthis1$CancerComLab, 
                               levels=c("CAMA1 \n Resistant","CAMA1 \n Sensitive",
                                        "MCF7 \n Resistant" ,"MCF7 \n Sensitive",
                                        "T47D \n Resistant","T47D \n Sensitive"))


plotthis1[CompositionLab=="Cancer and unprimed fibroblast",CompositionLab2:="Fibroblast \n coculutre"]
plotthis1[,CompositionLab2:="Alone"]
cancerplot <- ggplot(plotthis1[Afatinib_uM==0][FibroComLab!="Alone"],aes(y= (RGR-DMSOaloneRGR),
                                                                         x= interaction(ResistanceLab,CancerCellName,CompositionLab) ,#CompositionLab,#CancerComLab, 
                                                                         group= interaction(ResistanceLab,CancerCellName,CompositionLab),col=CompositionLab,fill=CompositionLab))+
  stat_boxplot(geom="errorbar")+
  geom_boxplot(col="black")+theme_classic()+
  geom_point(position=position_dodge(width=1),pch=21,col="black")+
  theme_classic(base_size=26)+theme(aspect.ratio = 0.3)+
  #coord_trans(y="log")+
  facet_grid( wrapLab~., scales="free_y") +
  scale_color_manual(name="",labels=c("Cancer \n alone","Cancer and \n fibroblasts"),values=collabs[-2])+
  scale_fill_manual(name="",labels=c("Cancer \n alone","Cancer and \n fibroblasts"),values=collabs[-2])+
  labs(x="Cancer cell line", y="Cancer growth rate \n (compared to monoculture)")+
  scale_x_discrete(labels=rep(c("CAMA1 (CRR)","CAMA1 (ETR)",
                            "MCF7 (CRR)" ,"MCF7 (ETR)",
                            "T47D (CRR)","T47D (ETR)") ,2))+
  theme(aspect.ratio=1,legend.position = "bottom")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text.x = element_text(vjust=0.5,angle=90))+
  theme(legend.position = "none")


fibroblastplot
#ggsave(fibroblastplot,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Fibroblast_Experiments/coculture/20230223_Feng Afatinib Fulvestrant CAMA1 MCF7 T47D/Modeling_JG/MUTUALISM FINALIZED part1 fibroblast growth rate change under coculture DMSO and fulv treatmentlong.pdf", width=12, height=12, dpi=320)
ggsave(fibroblastplot,file="~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/MUTUALISM FINALIZED part1 fibroblast growth rate change under coculture DMSO and fulv treatmentlong.pdf",height=10,width=12,dpi = 320)

cancerplot 
#ggsave(cancerplot,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Fibroblast_Experiments/coculture/20230223_Feng Afatinib Fulvestrant CAMA1 MCF7 T47D/Modeling_JG/MUTUALISM FINALIZED part1 cancer growth rate change under coculture DMSO and fulv treatmentlong.pdf", width=12, height=12, dpi=320)
ggsave(cancerplot,file="~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/MUTUALISM FINALIZED part1 cancer growth rate change under coculture DMSO and fulv treatmentlong.pdf",height=10,width=12,dpi = 320)


### Cancer stats
plotthis1[,CancerCellNameResistance:=paste0(CancerCellName,ResistanceLab)]
#Cancer facilitation: DMSO stats
rbindlist(lapply(unique(plotthis1$CancerCellNameResistance)[],function(iii){
  data.table(Fulvestrant_nMLab="Control",
             CancerCellNameResistance=iii,coef(summary(lm( I(RGR-DMSOaloneRGR)~-1+CompositionLab,
                                                           plotthis1[Afatinib_uM==0][FibroComLab!="Alone"][CancerCellNameResistance==iii][Fulvestrant_nMLab=="Control"]))),keep.rownames = T)[rn!="CompositionLabCancer alone"]
  
}))[order(CancerCellNameResistance)]

#Cancer facilitation: Fulv stats
FresCancer<-rbindlist(lapply(unique(plotthis1$CancerCellNameResistance)[],function(iii){
  data.table(Fulvestrant_nMLab="Fulvestrant: 5 nM",
             CancerCellNameResistance=iii,coef(summary(lm( I(RGR-DMSOaloneRGR)~-1+CompositionLab,
                                                           plotthis1[Afatinib_uM==0][FibroComLab!="Alone"][CancerCellNameResistance==iii][Fulvestrant_nMLab=="Fulvestrant: 5 nM"]
             ))),keep.rownames = T)[rn!="CompositionLabCancer alone"]
  
}))

rbindlist(lapply(unique(plotthis1$CancerCellNameResistance)[],function(iii){
  data.table(Fulvestrant_nMLab="Control",CancerCellNameResistance=iii,
             coef(summary(lm( I(RGRfibro-DMSOaloneRGRfibro)~ -1+FibroComLab,
                              plotthis1[Afatinib_uM==0][FibroComLab!="No fibroblast"][CancerCellNameResistance==iii|FibroComLab=="Alone"][Fulvestrant_nMLab=="Control"]))),keep.rownames = T)[rn!="FibroComLabAlone"]
  
}))

##### Fibro facilitated by cancer
summary(lmer( I(RGRfibro-DMSOaloneRGRfibro)~ -1+CompositionLab+(1|CancerCellNameResistance),
              plotthis1[Afatinib_uM==0][FibroComLab!="No fibroblast"][][Fulvestrant_nMLab=="Control"]))

summary(lmer( I(RGRfibro-DMSOaloneRGRfibro)~ -1+CompositionLab+(1|CancerCellNameResistance),
              plotthis1[Afatinib_uM==0][FibroComLab!="No fibroblast"][][Fulvestrant_nMLab=="Fulvestrant: 5 nM"]))

##### Cancer facilitated by Fibro
summary(lmer( I(RGR-DMSOaloneRGR)~ -1+CompositionLab+(1|CancerCellNameResistance),
              plotthis1[Afatinib_uM==0][FibroComLab!="Fibroblasts alone"][][Fulvestrant_nMLab=="Control"]))

summary(lmer( I(RGR-DMSOaloneRGR)~ -1+CompositionLab+(1|CancerCellNameResistance),
              plotthis1[Afatinib_uM==0][FibroComLab!="Fibroblasts alone"][][Fulvestrant_nMLab=="Fulvestrant: 5 nM"]))
####

out<- plotthis1 [Afatinib_uM==0]
# save( out, file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/Figure5/Cancer and fibroblast facilitation DMSO and fulv.RData")

#
cancerout <- plotthis1[Afatinib_uM==0][FibroComLab!="Alone"]
fibroout <-  plotthis1[Afatinib_uM==0][FibroComLab!="No fibroblast"]

write.csv(fibroout,"/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/Figure5/Fig5a/Fibroblast growth facilitated by cancer cells.csv")
write.csv(cancerout,"/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/Figure5/Fig5b/Cancer growth facilitated by cancer cells.csv")






