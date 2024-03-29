rm(list=ls()); require(data.table);require(dplyr);require(ggplot2);require(tidyr);
library(plater)
library(scales)

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
dat_long$CancerCell%>%unique()


# Split one column with 3 pieces of information into 3 different columns
dat_long[, c("CellType", "  ","Timehr") := tstrsplit(Time, "_", fixed=TRUE)]
dat_long[,"  ":=NULL]

# get numeric time
dat_long[, Hour := as.numeric( gsub("[^0-9]","",Time) )]
dat_long[,Timehr:=NULL]
dat_long[,Time:=NULL]
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
                            dplyr::select(RepID,Rep,Resistance,Fibroblasts,Cytation,Total_Cell_Number,InitialCancerFraction,Afatinib_uM,Fulvestrant_nM,Sensitivity,
                                   CancerCellName,FibroblastsPresent,CancerPresent,Afatinib_uMLab,Fulvestrant_nMLab,ResistanceLab,CompositionLab,maxHour))
                   ,RGR=RGR,RGRfibro=RGRfibro)
}))


plotthis0 <- plotthis00[Fibroblasts!="EMT.MRC5_C2"][Fulvestrant_nM%in%c(0,5)][Afatinib_uM%in%c(0,2.5)][InitialCancerFraction%in%c(0,0.5,1)]

plotthis0[, wrapLab:= "DMSO"]
plotthis0[Afatinib_uM==2.5&Fulvestrant_nM==0, wrapLab:= "Afatinib (2.5uM)"]
plotthis0[Afatinib_uM==0&Fulvestrant_nM==5, wrapLab:= "Fulvestrant (5nM)"]
plotthis0[Afatinib_uM==2.5&Fulvestrant_nM==5, wrapLab:= "Combination"]
plotthis0$wrapLab <- factor(plotthis0$wrapLab , levels=c("DMSO","Fulvestrant (5nM)","Afatinib (2.5uM)","Combination"))
plotthis0[CompositionLab=="Fibroblasts alone"][wrapLab=="DMSO"]
plotthis0[is.infinite(RGRfibro),RGRfibro:=NA]
plotthis0[is.infinite(RGR),RGR:=NA]

plotthis <- data.table(plotthis0%>% group_by(wrapLab,ResistanceLab,Fulvestrant_nM,Afatinib_uM,CompositionLab,CancerCellName)%>%
                         dplyr::summarise(RGR=mean(RGR,na.rm=T),
                                          RGRfibro=mean(RGRfibro,na.rm=T)))

refRGR <- data.table( plotthis[CompositionLab=="Cancer alone"][wrapLab=="DMSO"] %>% dplyr::select(ResistanceLab,CancerCellName,RGR))
refRGRful <- data.table( plotthis[CompositionLab=="Cancer alone"][wrapLab=="Fulvestrant (5nM)"] %>% dplyr::select(ResistanceLab,CancerCellName,RGR))
refRGRfibro <- data.table( plotthis[CompositionLab=="Fibroblasts alone"][wrapLab=="DMSO"] %>% dplyr::select(RGRfibro))
refRGRfulfibro <- data.table( plotthis[CompositionLab=="Fibroblasts alone"][wrapLab=="Fulvestrant (5nM)"] %>% dplyr::select(RGRfibro))
setnames(refRGR,old=c("RGR"), new=c("DMSOaloneRGR"))
setnames(refRGRful,old=c("RGR"), new=c("fulvaloneRGR"))
setnames(refRGRfibro,old=c("RGRfibro"), new=c("DMSOaloneRGRfibro"))
setnames(refRGRfulfibro,old=c("RGRfibro"), new=c("fulvaloneRGRfibro"))

plotthis1 <- merge(
  data.table(merge(plotthis0,refRGR,by=c("ResistanceLab", "CancerCellName"  )),"DMSOaloneRGRfibro"=mean(refRGRfibro$DMSOaloneRGRfibro), "fulvaloneRGRfibro"=mean(refRGRfulfibro$fulvaloneRGRfibro))
  , refRGRful,by=c("ResistanceLab", "CancerCellName"  ))

collabs<-(ggsci::pal_aaas()(6)[c(2,3,1)])


plotthis1[, FibroComLab:=paste0("Fibroblast + \n ", ResistanceLab," ",CancerCellName)]
plotthis1[, FibroComLab:=paste0( ResistanceLab," \n ",CancerCellName)]

plotthis1[CompositionLab=="Cancer alone", FibroComLab:="No fibroblast"]
plotthis1[CompositionLab=="Fibroblasts alone", FibroComLab:="Alone"]
plotthis1$FibroComLab<- factor(plotthis1$FibroComLab, 
                               levels=c("Alone","No fibroblast" ,
                                        "Resistant \n CAMA1","Sensitive \n CAMA1",
                                        "Resistant \n MCF7" ,"Sensitive \n MCF7",
                                        "Resistant \n T47D","Sensitive \n T47D"))



# ##"Fibroblast + \n Resistant T47D" ,"Fibroblast + \n Sensitive T47D" ,
# fibroblastplot <-ggplot(plotthis1[][FibroComLab!="No fibroblast"],aes(y= (RGRfibro-DMSOaloneRGRfibro),#/abs(DMSOaloneRGRfibro),
#                                                                       x= FibroComLab, 
#                                                                       group= interaction(FibroComLab),col=CompositionLab,fill=CompositionLab))+
#   geom_hline(yintercept=0, linetype="dashed")+
#   stat_boxplot(geom="errorbar")+
#   geom_boxplot(col="black")+theme_classic()+
#   geom_point(pch=21,col="black")+
#   theme_classic(base_size=26)+theme(aspect.ratio = 0.3)+
#   #coord_trans(y="log")+
#   facet_grid( Fulvestrant_nMLab~Afatinib_uMLab , scales="free_x") +
#   scale_color_manual(name="",values=collabs[-1])+
#   scale_fill_manual(name="",values=collabs[-1])+
#   labs(x="Fibroblast coculture", y="Fibroblast growth rate \n (relative to DMSO monoculutre)")+
#   #scale_x_discrete(labels=c("Cancer \n alone", "Cancer + \n unprimed fibroblasts", "Cancer + \n primed fibroblasts"))+
#   theme(aspect.ratio=1,legend.position = "bottom")+
#   theme(panel.border = element_rect(color = "black", fill = NA, size = 1) ,
#         axis.text.x = element_text(vjust=0.5,angle=90))#+
# 
# cancerplot <- ggplot(plotthis1[][FibroComLab!="Alone"],aes(y= (RGR-DMSOaloneRGR),
#                                                            x= interaction(ResistanceLab,CancerCellName), 
#                                                            group= interaction(ResistanceLab,CancerCellName,CompositionLab),col=CompositionLab,fill=CompositionLab))+
#   geom_hline(yintercept=0, linetype="dashed")+
#   stat_boxplot(geom="errorbar")+
#   geom_boxplot(col="black")+theme_classic()+
#   geom_point(position=position_dodge(width=1),pch=21,col="black")+
#   theme_classic(base_size=26)+theme(aspect.ratio = 0.3)+
#   #coord_trans(y="log")+
#   facet_grid(Fulvestrant_nMLab~Afatinib_uMLab , scales="free_y") +
#   scale_color_manual(name="",values=collabs[-2])+
#   scale_fill_manual(name="",values=collabs[-2])+
#   labs(x="Cancer coculture", y="Cancer growth rate \n (relative to DMSO monoculutre)")+
#   scale_x_discrete(labels=c("Resistant \n CAMA1","Sensitive \n CAMA1",
#                             "Resistant \n MCF7" ,"Sensitive \n MCF7",
#                             "Resistant \n T47D","Sensitive \n T47D"))+
#   theme(aspect.ratio=1,legend.position = "bottom")+
#   theme(panel.border = element_rect(color = "black", fill = NA, size = 1),axis.text.x = element_text(vjust=0.5,angle=90))#+
# 
# fibroblastplot
# #ggsave(fibroblastplot,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Fibroblast_Experiments/20230223_Feng Afatinib Fulvestrant CAMA1 MCF7 T47D/Modeling_JG/BLOCK MUTUALISM part1 fibroblast growth rate under coculture afatinib and fulv treatment.pdf", width=12, height=12, dpi=320)
# 
# cancerplot 
# #ggsave(cancerplot,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Fibroblast_Experiments/20230223_Feng Afatinib Fulvestrant CAMA1 MCF7 T47D/Modeling_JG/BLOCK MUTUALISM part1 cancer growth rate under coculture afatinib and fulv treatment.pdf", width=12, height=12, dpi=320)


plotthis1[,labs2:="DMSO"]
plotthis1[wrapLab=="Afatinib (2.5uM)",labs2:="Afatinib (2.5uM)"]

plotthis1[wrapLab=="Fulvestrant (5nM)",labs2:="Fulvestrant (5nM)"]
plotthis1[wrapLab=="Combination",labs2:="Afatinib (2.5uM) + \n Fulvestrant (5 nM)"]
plotthis1$labs2<-factor(plotthis1$labs2,levels=c("DMSO","Afatinib (2.5uM)","Fulvestrant (5nM)","Afatinib (2.5uM) + \n Fulvestrant (5 nM)"))

fibroblastplot <-ggplot(plotthis1[Fulvestrant_nM==0][FibroComLab!="No fibroblast"],aes(y= (RGRfibro-DMSOaloneRGRfibro),#/abs(DMSOaloneRGRfibro),
                                                                                       x= FibroComLab, 
                                                                                       group= interaction(FibroComLab),col=CompositionLab,fill=CompositionLab))+
  geom_hline(yintercept=0, linetype="dashed")+
  stat_boxplot(geom="errorbar")+
  geom_boxplot(col="black")+theme_classic()+
  geom_point(pch=21,col="black")+
  theme_classic(base_size=26)+theme(aspect.ratio = 0.3)+
  #coord_trans(y="log")+
  facet_grid(~labs2 , scales="free_x") +
  scale_color_manual(name="",values=collabs[-1])+
  scale_fill_manual(name="",values=collabs[-1])+
  scale_x_discrete(labels=c("Alone","CAMA1 (CRR)","CAMA1 (ETR)",
                            "MCF7 (CRR)" ,"MCF7 (ETR)",
                            "T47D (CRR)","T47D (ETR)"))+
  labs(x="Fibroblast coculture", y="Fibroblast growth rate \n (compared to DMSO monoculture)")+
  #scale_x_discrete(labels=c("Cancer \n alone", "Cancer + \n unprimed fibroblasts", "Cancer + \n primed fibroblasts"))+
  theme(aspect.ratio=1,legend.position = "none")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1) ,
        axis.text.x = element_text(vjust=0.5,angle=90))#+

# cancerplot <- ggplot(plotthis1[Fulvestrant_nM==0][FibroComLab!="Alone"],aes(y= (RGR-DMSOaloneRGR),
#                                                                             x= interaction(ResistanceLab,CancerCellName), 
#                                                                             group= interaction(ResistanceLab,CancerCellName,CompositionLab),col=CompositionLab,fill=CompositionLab))+
#   geom_hline(yintercept=0, linetype="dashed")+
#   stat_boxplot(geom="errorbar")+
#   geom_boxplot(col="black")+theme_classic()+
#   geom_point(position=position_dodge(width=1),pch=21,col="black")+
#   theme_classic(base_size=26)+theme(aspect.ratio = 0.3)+
#   #coord_trans(y="log")+
#   facet_grid(~labs2 , scales="free_y") +
#   scale_color_manual(name="",values=collabs[-2])+
#   scale_fill_manual(name="",values=collabs[-2])+
#   labs(x="Cancer cell line", y="Cancer growth rate \n (relative to DMSO monoculutre)")+
#   scale_x_discrete(labels=c("Resistant \n CAMA1","Sensitive \n CAMA1",
#                             "Resistant \n MCF7" ,"Sensitive \n MCF7",
#                             "Resistant \n T47D","Sensitive \n T47D"))+
#   theme(aspect.ratio=1,legend.position = "bottom")+
#   theme(panel.border = element_rect(color = "black", fill = NA, size = 1),axis.text.x = element_text(vjust=0.5,angle=90))#+

plotthis1$CancerComLab<- factor(plotthis1$CancerComLab, 
                                levels=c("CAMA1 \n Resistant","CAMA1 \n Sensitive",
                                         "MCF7 \n Resistant" ,"MCF7 \n Sensitive",
                                         "T47D \n Resistant","T47D \n Sensitive"))


plotthis1[CompositionLab=="Cancer and unprimed fibroblast",CompositionLab2:="Fibroblast \n coculutre"]
plotthis1[,CompositionLab2:="Alone"]

cancerplot <- ggplot(plotthis1[Fulvestrant_nM==0][FibroComLab!="Alone"],aes(y= (RGR-DMSOaloneRGR),
                                                                         x= interaction(ResistanceLab,CancerCellName,CompositionLab) ,#CompositionLab,#CancerComLab, 
                                                                         group= interaction(ResistanceLab,CancerCellName,CompositionLab),col=CompositionLab,fill=CompositionLab))+
  geom_hline(yintercept=0, linetype="dashed")+
  stat_boxplot(geom="errorbar")+
  geom_boxplot(col="black")+theme_classic()+
  geom_point(position=position_dodge(width=1),pch=21,col="black")+
  theme_classic(base_size=26)+theme(aspect.ratio = 0.3)+
  #coord_trans(y="log")+
  #facet_grid( wrapLab~., scales="free_y") +
  facet_grid(~labs2 , scales="free_x") +
  scale_color_manual(name="",labels=c("Cancer \n alone","Cancer and \n unprimed fibroblasts"),values=collabs[-2])+
  scale_fill_manual(name="",labels=c("Cancer \n alone","Cancer and \n unprimed fibroblasts"),values=collabs[-2])+
  labs(x="Cancer cell line", y="Cancer growth rate \n (compared to monoculture)")+
  scale_x_discrete(labels=rep(c("CAMA1 (CRR)","CAMA1 (ETR)",
                                "MCF7 (CRR)" ,"MCF7 (ETR)",
                                "T47D (CRR)","T47D (ETR)") ,2))+
  theme(aspect.ratio=1,legend.position = "none")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text.x = element_text(vjust=0.5,angle=90))#+


fibroblastplot
#ggsave(fibroblastplot,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Fibroblast_Experiments/coculture/20230223_Feng Afatinib Fulvestrant CAMA1 MCF7 T47D/Modeling_JG/BLOCK MUTUALISM FINALIZED part2 fibroblast growth rate under coculture afatinib treatment.pdf", width=12, height=12, dpi=320)
#ggsave(fibroblastplot,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Fibroblast_Experiments/coculture/20230223_Feng Afatinib Fulvestrant CAMA1 MCF7 T47D/Modeling_JG/BLOCK MUTUALISM FINALIZED part2 fibroblast growth rate under coculture afatinib treatmentB.pdf", width=10, height=10, dpi=320)
ggsave(fibroblastplot,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Fibroblast_Experiments/coculture/20230223_Feng Afatinib Fulvestrant CAMA1 MCF7 T47D/Modeling_JG/BLOCK MUTUALISM FINALIZED part2 fibroblast growth rate under coculture afatinib treatmentB.pdf", width=10, height=10, dpi=320)

cancerplot 
#ggsave(cancerplot,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Fibroblast_Experiments/coculture/20230223_Feng Afatinib Fulvestrant CAMA1 MCF7 T47D/Modeling_JG/BLOCK MUTUALISM FINALIZED part2 cancer growth rate under coculture afatinib treatment.pdf", width=12, height=12, dpi=320)
#ggsave(cancerplot,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Fibroblast_Experiments/coculture/20230223_Feng Afatinib Fulvestrant CAMA1 MCF7 T47D/Modeling_JG/BLOCK MUTUALISM FINALIZED part2 cancer growth rate under coculture afatinib treatmentB.pdf", width=10, height=10, dpi=320)


#### Stats
##### Fibro facilitated by cancer
require(lmerTest)
plotthis1[,CancerCellNameResistance:=paste0(CancerCellName,ResistanceLab)]

summary(lmer( I(RGRfibro-DMSOaloneRGRfibro)~ -1+CompositionLab+(1|CancerCellNameResistance),
              plotthis1[Fulvestrant_nM==0][FibroComLab!="No fibroblast"][Afatinib_uMLab=="Control"]))
#Est=5.7e-3,se=6.4e-4,df=25,t=9.0,p=2.5e-9
summary(lmer( I(RGRfibro-DMSOaloneRGRfibro)~ CompositionLab+(1|CancerCellNameResistance),
              plotthis1[Fulvestrant_nM==0][FibroComLab!="No fibroblast"][Afatinib_uMLab=="Afatinib: 2.5 uM"]))


summary(lmer( I(RGRfibro-DMSOaloneRGRfibro)~ Afatinib_uMLab+(1|CancerCellNameResistance),
              plotthis1[Fulvestrant_nM==0][(CompositionLab=="Cancer and unprimed fibroblast") ]))
#Est=-0.0049,se=0.00075,df=29,t=-6.68,p=2.5e-7


summary(lmer( I(RGRfibro-DMSOaloneRGRfibro)~ CompositionLab+(1|CancerCellNameResistance),
              plotthis1[Fulvestrant_nM==0][
                (CompositionLab=="Fibroblasts alone"&Afatinib_uMLab=="Control" ) |
                (CompositionLab=="Cancer and unprimed fibroblast"&Afatinib_uMLab=="Afatinib: 2.5 uM" ) ]))
summary(lmer( I(RGRfibro-DMSOaloneRGRfibro)~ Afatinib_uMLab+(1|CancerCellNameResistance),
                             plotthis1[Fulvestrant_nM==0][CompositionLab=="Fibroblasts alone"]))
#Est=-5.0e-3,se=1.5e-3,df=14,t=-3.35,p=0.0047 

#summary(lmer( I(RGRfibro-DMSOaloneRGRfibro)~ -1+ CompositionLab*Afatinib_uMLab+(1|CancerCellNameResistance),
#              plotthis1[Fulvestrant_nM==0][FibroComLab!="No fibroblast"][]))

##### Cancer facilitated by Fibro
summary(lmer( I(RGR-DMSOaloneRGR)~ CompositionLab+(1|CancerCellNameResistance),
              plotthis1[Fulvestrant_nM==0][FibroComLab!="Alone"][Afatinib_uMLab=="Control"]))
# Est=5.45e-4,se=1.63e-4,df=29,t=3.33,p=0.00235

summary(lmer( I(RGR-DMSOaloneRGR)~ CompositionLab+(1|CancerCellNameResistance),
              plotthis1[Fulvestrant_nM==0][FibroComLab!="Alone"][][Afatinib_uMLab=="Afatinib: 2.5 uM"]))
# Est=4.4e-4,se=9.59e-4,df=29,t=-0.46,p=0.65
####

out<- plotthis1 [Fulvestrant_nM==0]
# save( out, file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/Figure6/Cancer and fibroblast blocking facilitation DMSO and afatinib.RData")
#write.csv( out, file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/Figure6/Cancer and fibroblast blocking facilitation DMSO and afatinib.csv")









fibroblastplot <-ggplot(plotthis1[Fulvestrant_nM==5][FibroComLab!="No fibroblast"],aes(y= (RGRfibro-fulvaloneRGRfibro),#/abs(DMSOaloneRGRfibro),
                                                                                       x= FibroComLab, 
                                                                                       group= interaction(FibroComLab),col=CompositionLab,fill=CompositionLab))+
  geom_hline(yintercept=0, linetype="dashed")+
  stat_boxplot(geom="errorbar")+
  geom_boxplot(col="black")+theme_classic()+
  geom_point(pch=21,col="black")+
  theme_classic(base_size=26)+theme(aspect.ratio = 0.3)+
  #coord_trans(y="log")+
  facet_grid(~labs2 , scales="free_x") +
  scale_color_manual(name="",values=collabs[-1])+
  scale_fill_manual(name="",values=collabs[-1])+
  labs(x="Fibroblast coculture", y="Fibroblast growth rate \n (relative to DMSO monoculutre)")+
  #scale_x_discrete(labels=c("Cancer \n alone", "Cancer + \n unprimed fibroblasts", "Cancer + \n primed fibroblasts"))+
  theme(aspect.ratio=1,legend.position = "bottom")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1) ,
        axis.text.x = element_text(vjust=0.5,angle=90))#+

cancerplot <- ggplot(plotthis1[Fulvestrant_nM==5][FibroComLab!="Alone"],aes(y= (RGR-fulvaloneRGR),
                                                                            x= interaction(ResistanceLab,CancerCellName), 
                                                                            group= interaction(ResistanceLab,CancerCellName,CompositionLab),col=CompositionLab,fill=CompositionLab))+
  geom_hline(yintercept=0, linetype="dashed")+
  stat_boxplot(geom="errorbar")+
  geom_boxplot(col="black")+theme_classic()+
  geom_point(position=position_dodge(width=1),pch=21,col="black")+
  theme_classic(base_size=26)+theme(aspect.ratio = 0.3)+
  #coord_trans(y="log")+
  facet_grid(~labs2 , scales="free_y") +
  scale_color_manual(name="",values=collabs[-2])+
  scale_fill_manual(name="",values=collabs[-2])+
  labs(x="Cancer cell line", y="Cancer growth rate \n (relative to DMSO monoculutre)")+
  scale_x_discrete(labels=c("Resistant \n CAMA1","Sensitive \n CAMA1",
                            "Resistant \n MCF7" ,"Sensitive \n MCF7",
                            "Resistant \n T47D","Sensitive \n T47D"))+
  theme(aspect.ratio=1,legend.position = "bottom")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),axis.text.x = element_text(vjust=0.5,angle=90))#+

fibroblastplot
#ggsave(fibroblastplot,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Fibroblast_Experiments/20230223_Feng Afatinib Fulvestrant CAMA1 MCF7 T47D/Modeling_JG/BLOCK MUTUALISM part2 fibroblast growth rate under coculture afatinib + fulv treatment.pdf", width=12, height=12, dpi=320)

cancerplot 
#ggsave(cancerplot,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Fibroblast_Experiments/20230223_Feng Afatinib Fulvestrant CAMA1 MCF7 T47D/Modeling_JG/BLOCK MUTUALISM part2 cancer growth rate under coculture afatinib + fulv  treatment.pdf", width=12, height=12, dpi=320)




plotthis1[, cancerIDcode:=paste(ResistanceLab,CancerCellName,sep="_")]


# statsFibro<-rbindlist(lapply( as.character(unique(plotthis1[!FibroComLab%in%c("Alone","No fibroblast")]$FibroComLab)) , function(i){
#   data.table(Celltype="Fibroblast",Coculture=i, coef(summary(lm( RGRfibro-fulvaloneRGRfibro~ Afatinib_uMLab*CompositionLab,
#                                                                  data= plotthis1[][Fulvestrant_nM==5][FibroComLab%in%c("Alone",i)]
#   ))) , keep.rownames = T)
# }))
statsFibro<-rbindlist(lapply( as.character(unique(plotthis1[!FibroComLab%in%c("Alone","No fibroblast")]$FibroComLab)) , function(i){
  data.table(Celltype="Fibroblast",Coculture=i, coef(summary(lm( RGRfibro-DMSOaloneRGRfibro~ Afatinib_uMLab*CompositionLab,
                                                                 data= plotthis1[][Fulvestrant_nM==0][FibroComLab%in%c("Alone",i)]
  ))) , keep.rownames = T)
}))
setnames(statsFibro,old=c("Std. Error","t value","Pr(>|t|)"), new=c("Std.Error","tvalue","pvalue"))
statsFibro[pvalue<0.05][rn=="Afatinib_uMLabAfatinib: 2.5 uM:CompositionLabCancer and unprimed fibroblast"]
statsFibro[pvalue<0.05][rn=="CompositionLabCancer and unprimed fibroblast"]
statsFibro[rn=="Afatinib_uMLabAfatinib: 2.5 uM"] [pvalue<0.05]
statsFibro[rn=="CompositionLabCancer and unprimed fibroblast"][pvalue<0.05]

#Resistant CAMA1:est=0.0049,se=0.0012,t=3.90,p=0.001;Sensitive CAMA1:est=0.0052,se=0.0012,t=4.34,p=0.0004;Resistant MCF7:est=0.0029,se=0.0012,t=2.47,p=0.024;Sensitive MCF7:est=0.0034,se=0.0014,t=2.54,p=0.020;Resistant T47D:est=0.0037,se=0.0012,t=3.08,p=0.006;Sensitive T47D:est=0.0037,se=0.0013,t=2.92,p=0.009


plotthis1[,CancerCellNameFibcomp:=cancerIDcode]
plotthis1[CancerPresent==F,CancerCellNameFibcomp:="None"]
summary(lm( RGRfibro-DMSOaloneRGRfibro~ Afatinib_uMLab,
            data= plotthis1[][Fulvestrant_nM==0][CancerCellNameFibcomp=="None"]))
#Afatinib_uMLabAfatinib: 2.5 uM -0.0054069  0.0009858  -5.485 8.04e-05 ***


statsFibro<-rbindlist(lapply( as.character(unique(plotthis1[!FibroComLab%in%c("Alone","No fibroblast")]$CancerCellNameFibcomp)) , function(i){
  data.table(Celltype="Fibroblast",Coculture=i, coef(summary(lm( RGRfibro-DMSOaloneRGRfibro~ Afatinib_uMLab,
                                                                 data= plotthis1[CancerCellNameFibcomp==i][Fulvestrant_nM==0][!FibroComLab%in%c("No fibroblast")]
  ))) , keep.rownames = T)
}))
setnames(statsFibro,old=c("Std. Error","t value","Pr(>|t|)"), new=c("Std.Error","tvalue","pvalue"))
statsFibro[pvalue<0.05][rn=="Afatinib_uMLabAfatinib: 2.5 uM:CompositionLabCancer and unprimed fibroblast"]
statsFibro[pvalue<0.05][rn=="CompositionLabCancer and unprimed fibroblast"]
statsFibro[pvalue<0.05][rn=="Afatinib_uMLabAfatinib: 2.5 uM"]

Resistant CAMA1:est=-0.0031,se=0.001,t=-3.01,p=0.040;Sensitive CAMA1:est=-0.0039,se=0.001,t=-5.65,p=0.048;Resistant MCF7:est=-0.0037,se=0.0003,t=-11.89,p=0.0003;Sensitive MCF7:est=-0.0084se=0.002,t=-4.94,p=0.0078;Resistant T47D:est=-0.0014,se=0.0004,t=-3.30,p=0.03;Sensitive T47D:est=-0.0045,se=0.001,t=-3.76,p=0.02

# statsFibroB<-rbindlist(lapply( as.character(unique(plotthis1[!FibroComLab%in%c("Alone","No fibroblast")]$FibroComLab)) , function(i){
#   data.table(Celltype="Fibroblast",Coculture=i, coef(summary(lm( RGRfibro-DMSOaloneRGRfibro~ Afatinib_uMLab,
#                                                                  data= plotthis1[][Fulvestrant_nM==0][FibroComLab%in%c(i)]
#   ))) , keep.rownames = T)
# }))
# setnames(statsFibroB,old=c("Std. Error","t value","Pr(>|t|)"), new=c("Std.Error","tvalue","pvalue"))
# statsFibroB[pvalue<0.05][rn=="Afatinib_uMLabAfatinib: 2.5 uM"]
# 


# statsCancer <- rbindlist(lapply( as.character(unique(plotthis1[!FibroComLab%in%c("Alone","No fibroblast")]$cancerIDcode)) , function(i){
#   data.table(Celltype="Cancer",cancerIDcode=i, coef(summary(lm( I(RGR-fulvaloneRGR) ~ Afatinib_uMLab*CompositionLab,
#                                                                  data= plotthis1[][Fulvestrant_nM==5][!FibroComLab%in%c("Alone")][cancerIDcode==i]
# 
#   ))) , keep.rownames = T)
# }))
# setnames(statsCancer,old=c("Std. Error","t value","Pr(>|t|)"), new=c("Std.Error","tvalue","pvalue"))
# statsCancer[pvalue<0.05][rn!="(Intercept)"]
# statsCancer[cancerIDcode=="Sensitive_T47D"]
# statsCancer[cancerIDcode=="Resistant_MCF7"]
#statsCancer[cancerIDcode=="Resistant_CAMA1"]

statsCancer <- rbindlist(lapply( as.character(unique(plotthis1[!FibroComLab%in%c("Alone","No fibroblast")]$cancerIDcode)) , function(i){
  data.table(Celltype="Cancer",cancerIDcode=i, coef(summary(lm( I(RGR-DMSOaloneRGR) ~ Afatinib_uMLab*CompositionLab,
                                                                data= plotthis1[][Fulvestrant_nM==0][!FibroComLab%in%c("Alone")][cancerIDcode==i]
                                                                
  ))) , keep.rownames = T)
}))
setnames(statsCancer,old=c("Std. Error","t value","Pr(>|t|)"), new=c("Std.Error","tvalue","pvalue"))
statsCancer[pvalue<0.05][rn!="(Intercept)"]
statsCancer[cancerIDcode=="Sensitive_T47D"]
statsCancer[cancerIDcode=="Resistant_MCF7"]
statsCancer[cancerIDcode=="Resistant_CAMA1"]


statsCancerA <- rbindlist(lapply( as.character(unique(plotthis1[!FibroComLab%in%c("Alone","No fibroblast")]$cancerIDcode)) , function(i){
  KW <-     kruskal.test( I(RGR-DMSOaloneRGR)~ Afatinib_uMLab,    data= plotthis1[][Fulvestrant_nM==0][FibroComLab%in%c("No fibroblast")][cancerIDcode==i])
  est=data.table(coef(summary(lm( I(RGR-DMSOaloneRGR)~ Afatinib_uMLab, data= plotthis1[][Fulvestrant_nM==0][FibroComLab%in%c("No fibroblast")][cancerIDcode==i]))))$Estimate[-1]
  data.table(Celltype="Cancer",cancerIDcode=i, method=KW$method, "Kruskal-Wallis chi-squared"=KW$statistic ,est=est,df=KW$parameter,pvalue=KW$p.value)
}))

statsCancerB <- rbindlist(lapply( as.character(unique(plotthis1[!FibroComLab%in%c("Alone","No fibroblast")]$cancerIDcode)) , function(i){
  KW <-     kruskal.test( I(RGR-DMSOaloneRGR)~ Afatinib_uMLab,    data= plotthis1[][Fulvestrant_nM==0][!FibroComLab%in%c("Alone","No fibroblast")][cancerIDcode==i])
  est=data.table(coef(summary(lm( I(RGR-DMSOaloneRGR)~ Afatinib_uMLab, data= plotthis1[][Fulvestrant_nM==0][FibroComLab%in%c("Alone","No fibroblast")][cancerIDcode==i]))))$Estimate[-1]
  data.table(Celltype="Cancer",cancerIDcode=i, method=KW$method, "Kruskal-Wallis chi-squared"=KW$statistic ,est=est,df=KW$parameter,pvalue=KW$p.value)
}))

statsCancer<- rbind(data.table(Coculture="Cancer and Fibroblasts",statsCancerB),data.table(Coculture="Cancer Alone",statsCancerA))
statsCancer[Coculture=="Cancer and Fibroblasts"][pvalue<0.05]
statsCancer[Coculture=="Cancer and Fibroblasts"][pvalue>0.05]
statsCancer[Coculture=="Cancer Alone"][pvalue<0.05]

# coculture p={Resistant CAMA1=0.04630159,Sensitive CAMA1=0.04953461,Resistant MCF7=0.04630159,Sensitive MCF7=0.04953461,Resistant T47D=0.04630159,Sensitive T47D=0.04953461}
#  p={Resistant CAMA1=0.04953461,Sensitive CAMA1=0.04953461,Resistant MCF7=0.82725935,Sensitive MCF7=0.04953461,Resistant T47D=0.51269076,Sensitive T47D=0.04953461}

statsCancerA <- rbindlist(lapply( as.character(unique(plotthis1[!FibroComLab%in%c("Alone","No fibroblast")]$cancerIDcode)) , function(i){
  KW <-     kruskal.test( I(RGR-fulvaloneRGR)~ Afatinib_uMLab,    data= plotthis1[][Fulvestrant_nM==5][FibroComLab%in%c("No fibroblast")][cancerIDcode==i])
  est=data.table(coef(summary(lm( I(RGR-fulvaloneRGR)~ Afatinib_uMLab, data= plotthis1[][Fulvestrant_nM==5][FibroComLab%in%c("No fibroblast")][cancerIDcode==i]))))$Estimate[-1]
  data.table(Celltype="Cancer",cancerIDcode=i, method=KW$method, "Kruskal-Wallis chi-squared"=KW$statistic ,est=est,df=KW$parameter,pvalue=KW$p.value)
}))

statsCancerB <- rbindlist(lapply( as.character(unique(plotthis1[!FibroComLab%in%c("Alone","No fibroblast")]$cancerIDcode)) , function(i){
  KW <-     kruskal.test( I(RGR-fulvaloneRGR)~ Afatinib_uMLab,    data= plotthis1[][Fulvestrant_nM==5][!FibroComLab%in%c("Alone","No fibroblast")][cancerIDcode==i])
  lmmod=data.table(coef(summary(lm( I(RGR-fulvaloneRGR)~ Afatinib_uMLab, data= plotthis1[][Fulvestrant_nM==5][FibroComLab%in%c("Alone","No fibroblast")][cancerIDcode==i]))))
  est<-lmmod$Estimate[-1]
  unlist( lmmod[2,"Pr(>|t|)"] )
  data.table(Celltype="Cancer",cancerIDcode=i, method=KW$method, "Kruskal-Wallis chi-squared"=KW$statistic ,est=est,df=KW$parameter,pvalue=KW$p.value)
}))
statsCancer<- rbind(data.table(Coculture="Cancer and Fibroblasts",statsCancerB),data.table(Coculture="Cancer Alone",statsCancerA))
#setnames(statsCancer,old=c("Std. Error","t value","Pr(>|t|)"), new=c("Std.Error","tvalue","pvalue"))
statsCancer[pvalue<0.05]
statsCancer[!pvalue<0.05]




summary(lmer( I(RGR-DMSOaloneRGR)~ Afatinib_uMLab+(1|cancerIDcode),
              data= plotthis1[][Fulvestrant_nM==0][!FibroComLab%in%c("Alone","No fibroblast")]
))


summary(lmer( I(RGR)~ Afatinib_uMLab*CompositionLab+(1|cancerIDcode),
              data= plotthis1[][Fulvestrant_nM==0][!FibroComLab%in%c("Alone")]
))


data.table( coef(summary(lm( RGRfibro-DMSOaloneRGRfibro~ Afatinib_uMLab*CompositionLab,
                             data= plotthis1[][Fulvestrant_nM==0][FibroComLab%in%c("Alone","Resistant \n CAMA1")]
))) , keep.rownames = T)

summary(lm( RGRfibro-DMSOaloneRGRfibro~ Afatinib_uMLab*CompositionLab,
            data= plotthis1[][Fulvestrant_nM==0][FibroComLab%in%c("Alone","Sensitive \n CAMA1")]
))

summary(lm( RGRfibro-DMSOaloneRGRfibro~ Afatinib_uMLab*CompositionLab,
            data= plotthis1[][Fulvestrant_nM==0][FibroComLab%in%c("Alone","Resistant \n MCF7")]
))

summary(lm( RGRfibro-DMSOaloneRGRfibro~ Afatinib_uMLab*CompositionLab,
            data= plotthis1[][Fulvestrant_nM==0][FibroComLab%in%c("Alone","Sensitive \n MCF7")]
))

summary(lm( RGRfibro-DMSOaloneRGRfibro~ Afatinib_uMLab*CompositionLab,
            data= plotthis1[][Fulvestrant_nM==0][FibroComLab%in%c("Alone","Resistant \n T47D")]
))

summary(lm( RGRfibro-DMSOaloneRGRfibro~ Afatinib_uMLab*CompositionLab,
            data= plotthis1[CancerCellName%in%c("T47D")][Fulvestrant_nM==0][FibroComLab%in%c("Alone","Sensitive \n T47D")]
))
