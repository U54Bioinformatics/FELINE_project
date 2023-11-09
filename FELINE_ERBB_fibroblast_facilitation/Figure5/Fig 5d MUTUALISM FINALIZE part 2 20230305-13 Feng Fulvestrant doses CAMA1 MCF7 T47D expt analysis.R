#  20230305-13_Feng Fulvestrant dose CAMA1 MCF7 T47D

rm(list=ls()); require(data.table);require(dplyr);require(ggplot2);require(tidyr);
library(plater)

fileLoc<-"/Users/jason/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Fibroblast_Experiments/coculture/20230305-13_Feng  Fulvestrant dose CAMA1 MCF7 T47D/RawDataQC/"
paste0(fileLoc,list.files(fileLoc)[grepl(".csv",list.files(fileLoc))])
nplates<- length(paste0(fileLoc,list.files(fileLoc)[grepl(".csv",list.files(fileLoc))]))
# Read Multiple Plates
dd <- data.table(read_plates(
  files =  paste0(fileLoc,list.files(fileLoc)[grepl(".csv",list.files(fileLoc))]),
  plate_names = paste0("Plate_", 1:nplates),
  well_ids_column = "Well",   sep = "," 
))[!is.na( CancerCell)]
dd<-dd[CancerCell!="CAMA1_FulR-C2"]

names(dd)<-gsub("YFP i" ,"YFP_i",names(dd))
names(dd)<-gsub("CFP i" ,"CFP_i",names(dd))
names(dd)<-gsub("red i" ,"red_i",names(dd))
names(dd)<-gsub("YFP_i" ,"Sens_I",names(dd))
names(dd)<-gsub("CFP_i" ,"Res_I",names(dd))
names(dd)<-gsub("red_i" ,"Fibro_I",names(dd))

dd[,CancerCell:=gsub("-","_",CancerCell)]

ggplot( dd[CancerCell=="MCF7_riboR_CER2"][is.finite(Res_Intensity_336hr)],#[Ratio_CancerCelltoFibro>0.7],
        aes(y=log(Res_Intensity_336hr), x=Fibroblasts))+
  geom_point()+ facet_grid(fulvestrant_nM~Ratio_CancerCelltoFibro,scale="free")

dd[, Resistance := "Sensitive"]
dd[grep("riboR_CER2", CancerCell), Resistance := "RiboR"]
dd[grep("FulR_C2", CancerCell), Resistance := "FulR"]
dd[,PercentFBS:=1]
dd[,Cytation:=NULL]
initialTotal_Cell_Number <- dd$Total_Cell_Number%>%unique()
dd[,Total_Cell_Number:=NULL]
dd[,RepID:=1:nrow(dd)]
dd[,Rep:=1:length(RepID), by=c("Plate","PercentFBS","CancerCell","Fibroblasts","Ratio_CancerCelltoFibro","fulvestrant_nM")]
dd[CancerCell=="MCF7_riboR_CER2"][is.finite(Res_Intensity_336hr)][fulvestrant_nM==0][Ratio_CancerCelltoFibro%in%c(0,1)]
ggplot( dd[CancerCell=="MCF7_riboR_CER2"][is.finite(Res_Intensity_336hr)][fulvestrant_nM==0],
        #[Ratio_CancerCelltoFibro==0.5] ,
        aes(y=log(Res_Intensity_336hr), x=Fibroblasts))+
  geom_point()+ facet_grid(CancerCell~Ratio_CancerCelltoFibro,scale="free")


dat2 <- dd %>% dplyr::select(c("RepID","Rep","Resistance","PercentFBS","CancerCell":"fulvestrant_nM","Sens_Intensity_24hr":"Sens_Intensity_336hr","Fibro_Intensity_24hr":"Fibro_Intensity_336hr", "Res_Intensity_24hr":"Res_Intensity_336hr"))
# Convert to Long format 
dat_long <- data.table(gather(dat2, Time, Measurement, 'Sens_Intensity_24hr':'Res_Intensity_336hr'))
dat_long$CancerCell%>%unique()
# remove minus signs from strings
dat_long$CancerCell<-gsub("CER2_","CER2" ,  gsub("-","_",dat_long$CancerCell))
dat_long$CancerCell<-gsub(" ","" ,  dat_long$CancerCell)
dat_long$Fibroblasts<- gsub("-","_",dat_long$Fibroblasts)
# Split one column with 3 pieces of information into 3 different columns
dat_long[, c("CellType", "  ","Timehr") := tstrsplit(Time, "_", fixed=TRUE)]
dat_long[,"  ":=NULL]
dat_long[, c("CellLine", "Sensitivity","FluorLabel") := tstrsplit(CancerCell, "_", fixed=TRUE)]
dat_long[Sensitivity=="riboR",Sensitivity:="RiboR"]

# get numeric time
dat_long[, Hour := as.numeric( gsub("[^0-9]","",Time) )]
dat_long[,Timehr:=NULL]

# adjust labeling to match Fengs standard curve data
dat_long[, Variable:=gsub(" ","_",CellType) ]
dat_long[, Variable:=paste(Variable,"Intensity",sep="_") ]
dat_long[,CellType:=NULL]
dat_long[,Time:=NULL]
dat_long[,FluorLabel:=NULL]
setnames(dat_long,old="Measurement",new="Value")

# As we are dealing with monocultures there is no fluoresence if R cells in S culture and vica versa. 
dat_long[(Sensitivity=="Sens"&Variable=="Res_Intensity" ),Value:=0 ]
dat_long[(Sensitivity=="RiboR"&Variable=="Sens_Intensity" ),Value:=0 ]
dat_long[Ratio_CancerCelltoFibro==0&Variable=="Res_Intensity" ,Value:=0 ]
dat_long[Ratio_CancerCelltoFibro==0&Variable=="Sens_Intensity" ,Value:=0 ]
dat_long[CancerCell== "MCF7_riboR_CER2",CancerCell:="MCF7_RiboR_Cer2"]
dat_long[CancerCell== "CAMA1_riboR_CER2",CancerCell:="CAMA1_RiboR_Cer2"]
dat_long[, Fibroblasts:=gsub("_C2","",Fibroblasts) ]
dat_long[, Fibroblasts:=gsub("_CER2","",Fibroblasts) ]

ggplot(dat_long[fulvestrant_nM==0][][Variable!="Fibro_Intensity"][Ratio_CancerCelltoFibro%in%c(1,0.8)][], 
       aes(y= log10(Value), x= Hour,col=Fibroblasts,fill=Fibroblasts))+
  geom_point()+   theme_classic(base_size =26)+theme(aspect.ratio = 1)+
  facet_grid(Sensitivity~CellLine, scales="free_y") +
  geom_smooth(method="gam", formula=y~s(x,k=4), se= T) +
  labs(y="Cancer cell intensity")


# Translate intensity to counts
inputdd <- data.table(dat_long%>%spread(Variable,Value))
inputdd2 <- inputdd
inputdd2[,Cancer_Intensity:=Res_Intensity+ Sens_Intensity]



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

inputdd[,CancerCell:=gsub("T47D_riboR_CER2","T47D_RiboR_Cer2",CancerCell)]
inputdd[,CancerCell:=gsub("CAMA1_FulR_C2","CAMA1_RiboR_Cer2",CancerCell)]

inputdd[,Sens_Count:= 0]
inputdd[,Res_Count:= 0]
inputdd[Sensitivity=="Sens", Sens_Count:= predict(CountModels$SensLm,newdata=inputdd[Sensitivity=="Sens"])]
inputdd[Sensitivity!="Sens", Res_Count:= predict(CountModels$ResLm,newdata=inputdd[Sensitivity!="Sens"])]

# remove intensity now that we have counts
inputdd[,RiboR_Intensity :=NULL]
inputdd[,Sens_Intensity:=NULL]

# create a single total count column instead of one count column per cell type and remove other redundant data from object 
inputdd[,Count:=Sens_Count + Res_Count ]
inputdd[Ratio_CancerCelltoFibro==0, Count:=0]

FibroIntenPerCell <- (inputdd[Ratio_CancerCelltoFibro==min(Ratio_CancerCelltoFibro)][Hour==24][fulvestrant_nM==0]$Fibro_Intensity)%>%median()/1000
inputdd[,FibroCount:=Fibro_Intensity/FibroIntenPerCell ]
inputdd[,Sens_Count:=NULL]
inputdd[,Res_Count:=NULL]
inputdd[,CancerCell:=NULL]
inputdd[Fibroblasts=="None"]$FibroCount%>%hist()
inputdd[Fibroblasts=="None",FibroCount:=0 ]
inputdd[,Fibroblasts:=gsub("-",".",Fibroblasts)]
inputdd[,FibroblastsPresent:=T]
inputdd[Fibroblasts=="None",FibroblastsPresent:=F]

setnames(inputdd, old="Ratio_CancerCelltoFibro", new="InitialCancerFraction")
inputdd[,CancerPresent:=T]
inputdd[InitialCancerFraction==0,CancerPresent:=F]
inputdd[,fulvestrant_nMLab:=paste0("Fulvestrant: ",fulvestrant_nM , " uM")]
inputdd[fulvestrant_nM==0,fulvestrant_nMLab:="Control"]
inputdd$fulvestrant_nMLab  <- factor(inputdd$fulvestrant_nMLab  ,  levels=c("Control","Fulvestrant: 1 uM", "Afatinib: 5 uM",   "Afatinib: 25 uM" ) )
inputdd[,ResistanceLab:="Resistant"]
inputdd[Resistance=="Sensitive",ResistanceLab:="Sensitive"]
inputdd[, initialCount:= sum((Hour==24)*Count), by= RepID]
inputdd[, initialFibroCount:= sum((Hour==24)*FibroCount), by= RepID]

ggplot(inputdd[Resistance!="Sensitive"][][][is.finite(Count)][InitialCancerFraction%in%c(0.5,1)][], 
       aes(y= log(Count/initialCount), x= Hour, 
           group= interaction(Fibroblasts),col=Fibroblasts,fill=Fibroblasts))+
  geom_point(pch=21,col="black")+theme_classic()+
  geom_path(aes(group=RepID))+
  facet_grid(paste0(CellLine,Sensitivity)~paste0("Fulvestrant (",fulvestrant_nM,"nM)"),scales="free_y")

require(scales)
hexcol <- hue_pal()(3)
inputdd[, CompositionLab:="Cancer alone"]
inputdd[InitialCancerFraction<1&Fibroblasts=="EMT_MRC5", CompositionLab:="Cancer and \n primed fibroblast"]
inputdd[InitialCancerFraction<1&Fibroblasts=="MRC5", CompositionLab:="Cancer and \n unprimed fibroblast"]
inputdd$CompositionLab <- factor(inputdd$CompositionLab , levels=c("Cancer alone","Cancer and \n unprimed fibroblast","Cancer and \n primed fibroblast"))

ggplot(inputdd[][Hour%in%c(264,336)][][is.finite(Count)][InitialCancerFraction%in%c(0.5,1)][], 
       aes(y= log2(Count/initialCount), x= sqrt(fulvestrant_nM), 
           group= interaction(CompositionLab),col=CompositionLab,fill=CompositionLab))+
  geom_jitter(pch=21,col="black", height=0,width=0.05)+theme_classic(base_size=26)+
  theme(aspect.ratio=1)+geom_smooth(method="gam", formula=y~s(x,k=4),se=T)+
  #geom_path(aes(group=RepID))+
  facet_grid(CellLine~ResistanceLab,scales="free_y")+
  scale_color_manual(name="Composition",values=rev(hexcol))+
  scale_fill_manual(name="Composition",values=rev(hexcol))+
  labs(y="Cancer growth \n (log2 fold change during treatment)",x="Fulvestrant (nM)")+
  scale_x_continuous(breaks=sqrt(unique(inputdd$fulvestrant_nM)) , labels=(unique(inputdd$fulvestrant_nM)) )

#ggsave(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Fibroblast_Experiments/20230305-13_Feng  Fulvestrant dose CAMA1 MCF7 T47D/Modeling_JG/Cancer abundance in mon vs coculture across fulv doses and cell types.pdf", width=11, height=11,dpi=320)

show_col(ggsci::pal_aaas()(6))
collabs<-(ggsci::pal_aaas()(6)[c(2,3,1)])

plotthis00 <-rbindlist(lapply(unique(inputdd[Fibroblasts%in%c("None","MRC5")][is.finite(Count)][InitialCancerFraction%in%c(0.5,1)]$RepID),function(i){
  cat(i)
  tmp<- inputdd[RepID==i][Fibroblasts%in%c("None","MRC5")][is.finite(Count)][InitialCancerFraction%in%c(0.5,1)]
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
                            select(RepID,Rep,Resistance,PercentFBS,Fibroblasts,,InitialCancerFraction,fulvestrant_nM,CellLine,Sensitivity,
                                   FibroblastsPresent,CancerPresent,fulvestrant_nMLab,ResistanceLab,CompositionLab))
                   ,RGR=RGR,RGRfibro=RGRfibro)
}))



#plotthis<-inputdd[Fibroblasts%in%c("None","MRC5")][Hour%in%c(336)][][is.finite(Count)][InitialCancerFraction%in%c(0.5,1)][]
#plotthis[,RGR:=log(Count/initialCount)/Hour]

plotthis<-plotthis00
plotthis[,FibroPresent:="Added"]
plotthis[Fibroblasts=="None",FibroPresent:="None"]

plotthis
plotthis[fulvestrant_nM==0][CellLine=="CAMA1"][Sensitivity=="Sens"]
plotthis[,RGRalone:=sum(RGR*(CompositionLab=="Cancer alone"))/sum(CompositionLab=="Cancer alone"),by=c("fulvestrant_nM","CellLine","Sensitivity")]

ggplot(plotthis, aes(y= (RGR-RGRalone), x= sqrt(fulvestrant_nM), 
                     group= interaction(CompositionLab), col= CompositionLab, fill= CompositionLab))+
  stat_boxplot(geom="errorbar",aes(  group= interaction(CellLine, ResistanceLab,fulvestrant_nM,CompositionLab)))+
  geom_boxplot(col="black",aes(  group= interaction(CellLine, ResistanceLab,fulvestrant_nM,CompositionLab)))+
  geom_point(aes(  group= interaction(CellLine, ResistanceLab,fulvestrant_nM,CompositionLab)),pch=21, col="black", position=position_dodge(width=1)) +
  theme_classic(base_size= 26)+
  theme(aspect.ratio=1) + 
  geom_smooth(method="gam", formula=y~s(x,k=3,bs="cr"),se=T)+
  scale_color_manual(name="",values=collabs[-2],labels=c("Cancer \n alone","Cancer and \n fibroblasts"))+
  scale_fill_manual(name="",values=collabs[-2],labels=c("Cancer \n alone","Cancer and \n fibroblasts"))+
  labs(y="Cancer growth rate \n (compared to treatment monoculture)",x="Fulvestrant (nM)")+
  scale_x_continuous(breaks=sqrt(unique(inputdd$fulvestrant_nM)) , labels=(unique(inputdd$fulvestrant_nM)) )+
  theme(aspect.ratio=1,legend.position = "none")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))#,axis.text.x = element_text(vjust=0.5,angle=90))#+
ggsave( file="~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/FACILITATION Cancer growth rate in coculture across fulv doses relative to monoculutre across cancer cell lines.pdf",height=8.3,width=7,dpi = 320)
ggsave( file="~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/FACILITATION Cancer growth rate in coculture across fulv doses relative to monoculutre across cancer cell lines.pdf",height=7,width=7,dpi = 320)


##### Cancer facilitated by Fibro
plotthis[,CancerCellNameResistance:=paste0(CellLine,ResistanceLab)]

summary(lmer( I(RGR-RGRalone)~ sqrt(fulvestrant_nM)*CompositionLab+(1+sqrt(fulvestrant_nM)|CancerCellNameResistance),
              plotthis))
####


outFulv<- plotthis
# save( outFulv, file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/Figure5/Cancer facilitation under increasing fulv.RData")

#write.csv( outFulv, file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/Figure5/Fig5d/Cancer facilitation under increasing fulv.csv")
