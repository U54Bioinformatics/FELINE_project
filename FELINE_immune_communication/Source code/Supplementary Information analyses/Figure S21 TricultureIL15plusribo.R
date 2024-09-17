rm(list=ls()); require(data.table);require(dplyr);require(ggplot2);require(tidyr);
library(colorspace)
require(segmented)
fileLocation <- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/SI data/FigureS21/"
saveName <- "SourceData_FigureS21_TricultureIL15plusribo.csv"
# Load the data
dd <- data.table(read.csv(paste0(fileLocation, saveName)))[!is.na(CancerCell)]
dd[,RepID:=1:nrow(dd)]
dd[,Rep:=1:length(RepID), by=c("Plate","CancerCell","Resistance","Cancer_Number","THP1_M0_Number","Tcell_Number","Ribociclib_uM","IL15_ngmL")]
dd$CancerCell%>%unique()
dd[CancerCell=="CAMA1_SensV2",CancerCell:="CAMA1_Sens_V2"]
dd[CancerCell=="CAMA1_RiboRCer2",CancerCell:="CAMA1_RiboR_Cer2"]
dd[CancerCell=="MCF7_SensV2",CancerCell:="MCF7_Sens_V2"]
dd[CancerCell=="MCF7_RiboRCer2",CancerCell:="MCF7_RiboR_Cer2"]

ggplot(dd[CancerCell=="None"],aes(y=log(ATP_288hr),x=log(Area_288hr),
                                  col=sqrt(IL15_ngmL),
                                  shape=as.factor(Ribociclib_uM)))+geom_point()
tcellATPdat<-dd[CancerCell=="None"]
ggplot(tcellATPdat, aes(y= log(ATP_288hr/Area_288hr), x= sqrt(IL15_ngmL),
                                  col=as.factor(Ribociclib_uM))) + geom_point()
TcellATPmod <- lm(log(ATP_288hr/Area_288hr) ~ sqrt(IL15_ngmL) * Ribociclib_uM , data= tcellATPdat)
summary(TcellATPmod)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                   -3.19486    0.03561 -89.726  < 2e-16 ***
#  sqrt(IL15_ngmL)                2.05869    0.05630  36.567  < 2e-16 ***
#  Ribociclib_uM                  0.10060    0.05036   1.998   0.0519 .  
#sqrt(IL15_ngmL):Ribociclib_uM -0.53654    0.07962  -6.739 2.76e-08 ***
  
Tcellnewdat<- expand.grid(IL15_ngmL=seq(0,1,by=0.01),Ribociclib_uM=c(0,1))
TcellY<-predict(TcellATPmod,newdata=Tcellnewdat,se.fit=T)
Tcellnewdat$preds<-TcellY$fit
Tcellnewdat$ucl<-TcellY$fit+1.96*TcellY$se.fit
Tcellnewdat$lcl<-TcellY$fit-1.96*TcellY$se.fit

ggplot(tcellATPdat, aes(y= (ATP_288hr/Area_288hr), x= (IL15_ngmL),
                       col= as.factor(Ribociclib_uM), fill= as.factor(Ribociclib_uM)))+
  theme_classic(base_size = 26)+
  theme(aspect.ratio=1)+
  geom_ribbon(data=Tcellnewdat,aes(y=exp(preds),ymax=exp(ucl), ymin=exp(lcl)),alpha=0.6,col=NA)+
  geom_line(data=Tcellnewdat,aes(y=exp(preds)))+
  geom_point()+
  labs(y="T cell ATP \n (per unit spheroid area)",x="IL-15 (ng/mL)")+
  scale_color_jco(name="Treatment",labels=c("Control","Ribociclib (1uM)"))+  scale_fill_jco(name="Treatment",labels=c("Control","Ribociclib (1uM)"))
svloc0<-"/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Tcell_Experiments/Exp60/Modelling_JG/"
ggsave(file= paste0(svloc0, "TcellATPperunitarea.png"), dpi=320, width=11, height=11)


# First we'll continue to work with florescence measured total intensity; filtering ATP conten
dat2 <- dd[CancerCell!="None"] %>% dplyr::select(c("RepID","Rep","CancerCell":"IL15_ngmL",
                               names(dd) [c( grep("YFP_",names(dd)), grep("TexasRed_",names(dd)))]
                               #"Sens_Intensity_0hr":"Sens_Intensity_146hr", "Res_Intensity_0hr":"Res_Intensity_146hr")
) )


# Convert to Long format 
dat_long <- data.table(gather(dat2, Time, Measurement, 
                              names(dd) [c( grep("YFP_",names(dd)), grep("TexasRed_",names(dd)))]
                              #'Sens_Intensity_0hr':'Res_Intensity_200hr'
))
dat_long$Ribociclib_uM <- as.factor(dat_long$Ribociclib_uM)
dat_long$IL15_ngmL <- as.factor(dat_long$IL15_ngmL)
#dat_long$Ratio_Cancer_Tcell <- as.factor(dat_long$Ratio_Cancer_Tcell)

# Split one column with 3 pieces of information into 3 different columns
dat_long[, c("CellLine", "Sensitivity","Fluoresence") := tstrsplit(CancerCell, "_", fixed=TRUE)]
dat_long[,Fluoresence:=NULL]

# get numeric time
dat_long[, Hour := as.numeric( gsub("[^0-9]","",Time) )]
# adjust labeling to match Fengs standard curve data
dat_long[,CellType:="Cancer"]
dat_long[grep("TexasRed", Time),CellType:="Macrophage"]

dat_long[, Variable:="Sens_Intensity" ]
dat_long[grep("TexasRed_",Time), Variable:="Res_Intensity" ]
dat_long[,Time:=NULL]
setnames(dat_long,old="Measurement",new="Value")

# As we are dealing with monocultures there is no fluoresence if R cells in S culture and vica versa. 
dat_long[(Sensitivity=="Sens"&Variable=="Res_Intensity" ),Value:=0 ]
dat_long[(Sensitivity=="RiboR"&Variable=="Sens_Intensity" ),Value:=0 ]

# Translate intensity to counts
inputdd <- data.table(dat_long[CellType=="Cancer"][THP1_M0_Number%in%c(0,2500)]%>%spread(Variable,Value))


# load regression model fits obtained from Feng's standard curve mapping data (translate intensity to counts )
getCountModels <- function(){
  # Load the data
  dd <- data.table( read.csv(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/Lab Fibroblast/FengCellNumberIntensities_Parsed_Clean.csv") )
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
inputdd[Sensitivity=="Sens", Sens_Count:= predict(CountModels$SensLm,newdata=inputdd[Sensitivity=="Sens"])]
inputdd[Sensitivity=="RiboR", Res_Count:= predict(CountModels$ResLm,newdata=inputdd[Sensitivity=="RiboR"])]
#inputdd[Sensitivity=="Sens", Sens_Count:= 5.937333e-06*Sens_Intensity]
#inputdd[Sensitivity=="RiboR", Res_Count:= 5.937333e-06*RiboR_Intensity]
#ggplot(inputdd,aes(Sens_Intensity,Sens_Count))+geom_point()

# remove intensity now that we have counts
inputdd[,RiboR_Intensity :=NULL]
inputdd[,Sens_Intensity:=NULL]

# create a single total count column instead of one count column per cell type and remove other redundant data from object 
inputdd[,Count:=Sens_Count+Res_Count ]
inputdd[,Sens_Count:=NULL]
inputdd[,Res_Count:=NULL]
inputdd[,CancerCell:=NULL]
#inputdd[,FluorLabel:=NULL]

# create a log count variable
inputdd[,lnCount:= log(Count)]

# Calculate Initial Tcell Fraction (note  Ratio_Cancer_Tcell==0 indicates 100% cancer not 0 )
#inputdd[,InitialTcellFraction:=Fraction_Tcell]
#inputdd[Ratio_Cancer_Tcell==0,InitialTcellFraction:=0]

# labels for treatment and cell lineage on plots
inputdd[,Treatlab:="Control"]
inputdd[Ribociclib_uM!=0,Treatlab:="Ribociclib"]

inputdd[,Lineagelab:="Sensitive"]
inputdd[Sensitivity=="RiboR",Lineagelab:="Resistant"]
inputdd$Sensitivity <- factor(inputdd$Sensitivity , levels=c("Sens","RiboR"))


inputdd[,Composition:="Cancer + Macrophages + T cells"]
inputdd[THP1_M0_Number==0 &Tcell_Number==0 ,Composition:="Cancer alone"]
inputdd[THP1_M0_Number!=0 &Tcell_Number==0 ,Composition:="Cancer + Macrophages"]
inputdd[THP1_M0_Number==0 &Tcell_Number!=0 ,Composition:="Cancer + T cells"]

inputdd[,Cytokinetreatment:="No IL-15"]
inputdd[as.numeric(as.character(IL15_ngmL))>0,Cytokinetreatment:="IL-15"]

expt60dd <- inputdd
expt60dd[ ,Expt:="expt60"]
expt60dd[ ,IL15Level:= paste0("IL-15: ",IL15_ngmL," ngmL")]

alldd<-expt60dd%>%dplyr::select(Expt,CellLine,RepID,Rep,Count,lnCount,Hour,Ribociclib_uM,IL15_ngmL,Treatlab,Lineagelab,Composition,Tcell_Number,THP1_M0_Number,IL15Level)
alldd[,Cellcatlab:=  paste(CellLine,Lineagelab,sep=" ")]
alldd$Cellcatlab <- factor(alldd$Cellcatlab , levels=c( "CAMA1 Sensitive"))
alldd$IL15Level <- factor(alldd$IL15Level , levels=c( "IL-15: 0 ngmL","IL-15: 0.1 ngmL" , "IL-15: 0.5 ngmL","IL-15: 1 ngmL" ))

alldd$Composition <- factor(alldd$Composition , levels=c(  "Cancer alone"  ,"Cancer + Macrophages","Cancer + T cells", "Cancer + Macrophages + T cells"  ))
alldd[,MaxHour:=max(Hour), by=c("Expt","CellLine","RepID","Rep")]
alldd[,MinHour:=min(Hour), by=c("Expt","CellLine","RepID","Rep")]

svloc0<-"/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Tcell_Experiments/Exp60/Modelling_JG/"
  #svloc0<-"/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
ggplot(alldd[][][Ribociclib_uM%in%c(0,1)][], aes(x = Hour, y= log(1+Count), col= Treatlab)) +
  geom_point(size=3,aes(shape=Composition,
                        col = Treatlab)) +
  geom_path(aes(linetype=Composition,
                group=interaction(RepID,Cellcatlab)))+
  facet_grid(Composition~paste0("IL-15:",IL15_ngmL,"(ngmL)"))+ 
  #stat_summary(geom="line", fun.data= mean_se, aes(col= Ribociclib_uM)) +
  #stat_summary(geom="errorbar", fun.data= mean_se, aes(col= Ribociclib_uM)) +
  labs(x= "Time (Hours)", y= "Cancer cell number (log(1+x))")+
  scale_color_jco(name="Treatment")+  scale_fill_jco(name="Treatment")+
  theme_classic(base_size=26)+
  theme(aspect.ratio=1)
ggsave(file=paste0(svloc0,"TcellMacrophagetricultureTimeCourse_IL15_ribo.png"), dpi=320,width=18,height=18)


# Calculate average growth rate of cancer population throughout treatment
alldd[,Count0:=sum((Hour==MinHour)*Count,na.rm=T),by=c("Expt","CellLine","RepID","Rep") ]
Finaldd <- alldd[Ribociclib_uM%in%c(0,1)][Hour==MaxHour]
Finaldd[,RGR:=(log(Count)-log(Count0))/(MaxHour-MinHour) ]
Finaldd<-Finaldd[!is.na(RGR)]

# Create modified variables for fitting
Finaldd[,sqrt_IL15:=sqrt(as.numeric(as.character(IL15_ngmL)))]
Finaldd[,FactorCellLineID:=as.factor(Cellcatlab)]

# Image of all the data
ggplot(Finaldd[],
       aes( y= RGR, x = sqrt_IL15, col = Treatlab,fill=Treatlab, 
            group=interaction(Treatlab,Composition,
                              CellLine,Lineagelab)
       )) +
  #geom_smooth(aes(linetype=(Composition)),method="gam",se=F,alpha=0.6, formula=y~s(x,k=3))+
  geom_smooth(method="gam",se=T,alpha=0.6, formula=y~s(x,k=3))+
  geom_point(aes(col = Treatlab, shape=Composition,group=RepID)) +
  facet_grid(CellLine~Composition, scale="free")+theme_classic(base_size=6)+
  scale_color_jco(name="Treatment")+  scale_fill_jco(name="Treatment")+
  #scale_y_continuous(labels = (1:10)^2*2000, breaks=log((1:10)^2*2000))+
  #scale_y_continuous(labels = c(1000,2000,4000,8000,16000,32000,64000), breaks=log10(c(1000,2000,4000,8000,16000,32000,64000)))+
  labs(y="Cancer cell growth rate (mean(cells/cell/day))", x="IL-15 (ngmL)")+ 
  theme(aspect.ratio=1)+scale_linewidth(range=c(1,2))+scale_size( range = c(1.5, 4.5))


####################################
# Analysis of the data


#Composition = catag(InitialTcellFraction)   # FactorCellLineID=interaction(CellLine,Lineagelab)
statsdata <- Finaldd %>% dplyr::select(RepID,Rep,RGR,sqrt_IL15,Treatlab,Composition,Cellcatlab,CellLine,Lineagelab,FactorCellLineID)
statsdata[,RiboTreated:=0]
statsdata[Treatlab=="Ribociclib", RiboTreated:=1]
statsdata[,FacTreatlab:=as.factor(Treatlab)]
statsdata[, Ribo__sqrt_IL15:=RiboTreated*sqrt_IL15]#*(InitialTcellFraction==0)]
#statsdata[, Ribo__InitialTcellFraction:=RiboTreated*InitialTcellFraction]

statsdata[, isTCoculture:=0]
statsdata[Composition%in%c("Cancer + T cells","Cancer + Macrophages + T cells"), isTCoculture:=1]
statsdata[, isMCoculture:=0]
statsdata[Composition%in%c("Cancer + Macrophages","Cancer + Macrophages + T cells"), isMCoculture:=1]
statsdata[, isMonoculture:=0]
statsdata[Composition=="Cancer alone", isMonoculture:=1]

statsdata[, isControl:=0]
statsdata[Treatlab!="Ribociclib", isControl:=1]


#######
  cellXanalysis <- statsdata#[FactorCellLineID==x][]
  
  cellXanalysis[, Ribo__sqrt_IL15:=RiboTreated*sqrt_IL15]#*(InitialTcellFraction==0)]
  cellXanalysis[, isTCoculture_sqrt_IL15:=isTCoculture*sqrt_IL15]
  cellXanalysis[, isMCoculture_sqrt_IL15:=isMCoculture*sqrt_IL15]
  cellXanalysis[, Ribo__isTCoculture_sqrt_IL15:=RiboTreated*isTCoculture*sqrt_IL15]
  cellXanalysis[, Ribo__isMCoculture_sqrt_IL15:=RiboTreated*isMCoculture*sqrt_IL15]
  cellXanalysis[, isTCoculture__isMCoculture_sqrt_IL15:=isTCoculture*isMCoculture*sqrt_IL15]
  
  gam0 <- gam(RGR~ RiboTreated+ # Monoculture ribo effect
                isTCoculture+
                isMCoculture + #main coculutre effects
                s(sqrt_IL15, k=3,bs="ts")+  # Monoculture IL-15 effect
                s(Ribo__sqrt_IL15, k=3,bs="ts")+ # Monoculture ribo + IL-15 synergy
                RiboTreated:isTCoculture + # TCoculture impact on ribo effect (change in ribo effect in coculutre without IL15)
                RiboTreated:isMCoculture + # MCoculture impact on ribo effect (change in ribo effect in coculutre without IL15)
                isMCoculture:isTCoculture +
                s(isTCoculture_sqrt_IL15, k=3,bs="ts")+ # TCoculture impact on IL15 effect
                s(isMCoculture_sqrt_IL15, k=3,bs="ts")+ # MCoculture impact on IL15 effect
                
                s(Ribo__isTCoculture_sqrt_IL15, k=3,bs="ts")+# TCoculture impact on ribo + IL-15 synergy
                s(Ribo__isMCoculture_sqrt_IL15, k=3,bs="ts")+# MCoculture impact on ribo + IL-15 synergy
                s(isTCoculture__isMCoculture_sqrt_IL15, k=3,bs="ts")+
                RiboTreated:isTCoculture:isMCoculture,# Coculture impact on ribo + IL-15 synergy
              data=cellXanalysis)
  # Predict using model
  cellXanalysis[ , preds:=predict(gam0, type="response")]
  cellXanalysis[ , ucl:=predict(gam0, type="response")+predict(gam0,se=T, type="response")$se.fit]
  cellXanalysis[ , lcl:=predict(gam0, type="response")-predict(gam0,se=T, type="response")$se.fit]
  
  # Extract terms
  gam_terms <- cbind("(Intercept)"=coef(gam0)["(Intercept)"] ,predict(gam0,type="terms" )) #plot(rowSums(gam_terms),predict(gam0))
  #colnames(gam_terms)<-paste0("gam",c("Intercept", "MonoRiboEffect", "TCoRiboEffect","MCoRiboEffect",
   #                                   "MonoIL15Effect","MonoRiboIL15synergy",
   #                                    "TCoIL15Effect","TCoRiboIL15synergy",
  #                                    "MCoIL15Effect","MCoRiboIL15synergy"))
  colnames(gam_terms)<-paste0("gam",c("Intercept", "MonoRiboEffect","MainTCoculture","MainMCoculture",
                                      "TCoRiboEffect","MCoRiboEffect", "TMTriculture","TMTriRiboEffect",
                                      "MonoIL15Effect","MonoRiboIL15synergy",
                                      "TCoIL15Effect","MCoIL15Effect",
                                      "TCoRiboIL15synergy", "MCoRiboIL15synergy",
                                      "TMTriIL15Effect"))
  cellXanalysisTerms <- cbind(cellXanalysis,gam_terms)
  
  # Make fine scale predictions
  lu<- data.table(sqrt_IL15= unique(c(seq(min(cellXanalysisTerms$sqrt_IL15),max(cellXanalysisTerms$sqrt_IL15),length=50),
                                      unique(cellXanalysisTerms$sqrt_IL15)))
  )
  meta2<-unique(cellXanalysis%>%dplyr::select(Lineagelab,CellLine,RiboTreated,Treatlab,Composition,isTCoculture,isMCoculture))
  
  predsdd<- rbindlist(lapply(1:nrow(lu),function(i){
    data.table(lu[i], meta2)  
  }))
  
  predsdd[, Ribo__sqrt_IL15:=RiboTreated*sqrt_IL15]
  predsdd[, isTCoculture_sqrt_IL15:=isTCoculture*sqrt_IL15]
  predsdd[, isMCoculture_sqrt_IL15:=isMCoculture*sqrt_IL15]
  predsdd[, Ribo__isTCoculture_sqrt_IL15:=RiboTreated*isTCoculture*sqrt_IL15]
  predsdd[, Ribo__isMCoculture_sqrt_IL15:=RiboTreated*isMCoculture*sqrt_IL15]
  predsdd[, isTCoculture__isMCoculture_sqrt_IL15:=isTCoculture*isMCoculture*sqrt_IL15]
  
  
  predsdd[ , preds:=predict(gam0, newdata=predsdd ,type="response")]
  predsdd[ , ucl:=predict(gam0, newdata=predsdd ,type="response")+predict(gam0,se=T, newdata=predsdd ,type="response")$se.fit]
  predsdd[ , lcl:=predict(gam0, newdata=predsdd ,type="response")-predict(gam0,se=T, newdata=predsdd ,type="response")$se.fit]
  
  ModOut<-  list(gam0,cellXanalysisTerms,meta2,predsdd)
allgammods<- gam0
allXanalysisTerms<- cellXanalysisTerms
allpredsdd<- predsdd

summary(gam0)
# Goodness of fit
ggplot(allXanalysisTerms,aes( y = preds, x= RGR,col = Treatlab,fill=Treatlab, 
                              group=interaction(Treatlab,Composition,CellLine,Lineagelab))) +
  geom_abline(linetype="dashed")+
  geom_smooth(method="gam",se=F,alpha=0.6, formula=y~s(x,k=3))+
  geom_point(aes(col = Treatlab, shape=Composition,group=RepID)) +
  geom_point(aes(col = Treatlab, shape=Composition,group=RepID)) +
  facet_grid(CellLine~Lineagelab)+
  theme_classic(base_size=6)+
  scale_color_jco(name="Treatment")+  scale_fill_jco(name="Treatment")+
  labs(y="Predicted rowth rate \n (mean(cells/cell/day))", x="Observed rowth rate \n (mean(cells/cell/day))")+ 
  theme(aspect.ratio=1)+scale_linewidth(range=c(1,2))+scale_size( range = c(1.5, 4.5))

#lm stats
summary(lm(RGR~ (RiboTreated+ sqrt_IL15+isTCoculture+isMCoculture)^4,
           data=cellXanalysis))
summary(lm(RGR~ (RiboTreated+ sqrt_IL15+isTCoculture+isMCoculture)^2+
             RiboTreated:sqrt_IL15:isTCoculture +RiboTreated:sqrt_IL15:isMCoculture ,
           data=cellXanalysis))

# Predictions overlayed on data

normvals<-data.table(allpredsdd[][Composition=="Cancer alone"][sqrt_IL15==0]%>%group_by(Treatlab,Lineagelab,CellLine)%>%summarise(predsNorm=mean(preds)))
#normvals<-data.table(allpredsdd[][][sqrt_IL15==0]%>%group_by(Composition,Treatlab,Lineagelab,CellLine)%>%summarise(predsNorm=mean(preds)))
allXanalysisTermsB<-merge(allXanalysisTerms,normvals,by=c("Lineagelab", "CellLine","Treatlab" ))
allpredsddB<-merge(allpredsdd,normvals,by=c("Lineagelab", "CellLine","Treatlab"   ))
allXanalysisTermsB$Composition<- factor(allXanalysisTermsB$Composition, levels=c("Cancer alone","Cancer + T cells","Cancer + Macrophages","Cancer + Macrophages + T cells"))
allpredsddB$Composition<- factor(allpredsddB$Composition, levels=c("Cancer alone","Cancer + T cells","Cancer + Macrophages","Cancer + Macrophages + T cells"))

ggplot(allXanalysisTermsB,
       aes( y= RGR-predsNorm, x = sqrt_IL15, col = Treatlab,fill=Treatlab, 
            group=interaction(Treatlab,Composition,
                              CellLine,Lineagelab)
       )) +
  #geom_smooth(data=allXanalysisTermsB ,method="gam", se=T,formula=y~s(x,k=3),alpha=0.2,col=NA,aes(y=RGR-predsNorm))+
  geom_ribbon(data=allpredsddB ,alpha=0.2,col=NA,aes(y=preds-predsNorm,ymax=ucl-predsNorm,ymin=lcl-predsNorm))+
  geom_line(data=allpredsddB ,aes(y=preds-predsNorm))+#scale_linewidth(range=c(1,2))+scale_size( range = c(1.5, 4.5))+
  geom_point(size=3,alpha=1,aes(col = Treatlab, shape=Composition,group=RepID)) +
  facet_wrap(~Composition)+theme_classic(base_size=26)+
  scale_color_jco(name="Treatment")+  scale_fill_jco(name="Treatment")+
  #scale_x_continuous(labels = unique(allXanalysisTerms$sqrt_IL15)^2, breaks=unique(allXanalysisTerms$sqrt_IL15) )+
  scale_x_continuous(labels = c(0,0.1,0.5,1), breaks=sqrt(c(0,0.1,0.5,1)) )+
  #scale_x_continuous(labels = c(1000,2000,4000,8000,16000,32000,64000), breaks=log10(c(1000,2000,4000,8000,16000,32000,64000)))+
  labs(y="Cancer growth rate \n (relative to ribociclib treatment \n monoculture mean of IL-15 control)", x="IL-15 (ngmL)")+ 
  theme(aspect.ratio=1)+
  scale_shape_manual(values=c("triangle","circle","square","diamond"))
  
ggsave(file= paste0(svloc0,"Normalized TcellMacrophageTricultureRGR_IL15_ribo_FitRelativegrowth.png"), dpi=320,width=17,height=17)


ggplot(allXanalysisTermsB,
       aes( y= RGR, x = sqrt_IL15, col = Treatlab,fill=Treatlab, 
            group=interaction(Treatlab,Composition,
                              CellLine,Lineagelab)
       )) +
  #geom_smooth(data=allXanalysisTermsB ,method="gam", se=T,formula=y~s(x,k=3),alpha=0.2,col=NA,aes(y=RGR-predsNorm))+
  geom_ribbon(data=allpredsddB ,alpha=0.2,col=NA,aes(y=preds,ymax=ucl,ymin=lcl))+
  geom_line(data=allpredsddB ,aes(y=preds))+#scale_linewidth(range=c(1,2))+scale_size( range = c(1.5, 4.5))+
  geom_point(size=3,alpha=1,aes(col = Treatlab, shape=Composition,group=RepID)) +
  facet_wrap(~Composition)+theme_classic(base_size=26)+
  scale_color_jco(name="Treatment",labels=c("Control","Ribociclib (1uM)"))+  scale_fill_jco(name="Treatment",labels=c("Control","Ribociclib (1uM)"))+
  #scale_x_continuous(labels = unique(allXanalysisTerms$sqrt_IL15)^2, breaks=unique(allXanalysisTerms$sqrt_IL15) )+
  scale_x_continuous(labels = c(0,0.1,0.5,1), breaks=sqrt(c(0,0.1,0.5,1)) )+
  #scale_x_continuous(labels = c(1000,2000,4000,8000,16000,32000,64000), breaks=log10(c(1000,2000,4000,8000,16000,32000,64000)))+
  labs(y="Cancer growth rate", x="IL-15 (ngmL)")+ 
  theme(aspect.ratio=1)+
  scale_shape_manual(values=c("triangle","circle","square","diamond"))

ggsave(file= paste0(svloc0,"TcellMacrophageTricultureGAMRGR_IL15_ribo_FitRelativegrowth.png"), dpi=320,width=16,height=16)





allXanalysisTerms[,Cellcatlab:=  paste(CellLine,Lineagelab,sep=" ")]
allpredsdd[,Cellcatlab:=  paste(CellLine,Lineagelab,sep=" ")]

normvals<-data.table(allpredsdd[Treatlab=="Control"][Composition=="Cancer alone"]%>%group_by(Lineagelab,CellLine,sqrt_IL15)%>%summarise(predsNorm=mean(preds)))
allXanalysisTermsB<-merge(allXanalysisTerms,normvals,by=c("Lineagelab", "CellLine" , "sqrt_IL15"))
allpredsddB<-merge(allpredsdd,normvals,by=c("Lineagelab", "CellLine" , "sqrt_IL15"))

normvals<-data.table(allpredsdd[Treatlab=="Control"][Composition=="Cancer alone"][sqrt_IL15==0]%>%group_by(Lineagelab,CellLine)%>%summarise(predsNorm=mean(preds)))
allXanalysisTermsB<-merge(allXanalysisTerms,normvals,by=c("Lineagelab", "CellLine" ))
allpredsddB<-merge(allpredsdd,normvals,by=c("Lineagelab", "CellLine" ))

#normvals<-data.table(allpredsdd[][Composition=="Cancer alone"][sqrt_IL15==0]%>%group_by(Treatlab,Lineagelab,CellLine)%>%summarise(predsNorm=mean(preds)))
#allXanalysisTermsB<-merge(allXanalysisTerms,normvals,by=c("Lineagelab", "CellLine","Treatlab" ))
#allpredsddB<-merge(allpredsdd,normvals,by=c("Lineagelab", "CellLine","Treatlab"  ))

#svloc<-"~/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"

ggplot(allXanalysisTermsB,
       aes( y= RGR-predsNorm, x = sqrt_IL15, col = Treatlab,fill=Treatlab, 
            group=interaction(Treatlab,Composition,
                              CellLine,Lineagelab)
       )) +
  geom_ribbon(data=allpredsddB ,alpha=0.2,col=NA,aes(y=preds-predsNorm,ymax=ucl-predsNorm,ymin=lcl-predsNorm))+
  geom_line(data=allpredsddB ,aes(y=preds-predsNorm))+#scale_linewidth(range=c(1,2))+scale_size( range = c(1.5, 4.5))+
  geom_point(size=3,alpha=1,aes(col = Treatlab,group=RepID)) +
  facet_wrap(~Composition)+theme_classic(base_size=26)+
  scale_color_jco(name="Treatment")+  scale_fill_jco(name="Treatment")+
  #scale_x_continuous(labels = unique(allXanalysisTerms$sqrt_IL15)^2, breaks=unique(allXanalysisTerms$sqrt_IL15) )+
  scale_x_continuous(labels = c(0,0.1,0.5,1), breaks=sqrt(c(0,0.1,0.5,1)) )+
  #scale_x_continuous(labels = c(1000,2000,4000,8000,16000,32000,64000), breaks=log10(c(1000,2000,4000,8000,16000,32000,64000)))+
  labs(y="Cancer growth rate \n (relative to monoculture mean \n of ribociclib and IL-15 control)", x="IL-15 (ngmL)")+ 
  theme(aspect.ratio=1)
#ggsave(file= paste0(svloc0,"TcellMacrophageTricultureRGR_IL15_ribo_FitRelativegrowth.png"), dpi=320,width=15,height=15)


allXanalysisTermsB$Composition<- factor(allXanalysisTermsB$Composition, levels=c("Cancer alone","Cancer + T cells","Cancer + Macrophages","Cancer + Macrophages + T cells"))
allpredsddB$Composition<- factor(allpredsddB$Composition, levels=c("Cancer alone","Cancer + T cells","Cancer + Macrophages","Cancer + Macrophages + T cells"))

ggplot(allXanalysisTermsB,
       aes( y= RGR#-predsNorm
            , x = sqrt_IL15, col = Treatlab,fill=Treatlab, 
            group=interaction(Treatlab,Composition,
                              CellLine,Lineagelab)
       )) +
  geom_ribbon(data=allpredsddB ,alpha=0.2,aes(col=NA,y=preds#-predsNorm
                                              ,ymax=ucl#-predsNorm
                                              ,ymin=lcl#-predsNorm
                                              ))+
  geom_line(data=allpredsddB ,aes(#linetype=(Composition=="Cancer alone"),
                                  y=preds#-predsNorm
                                  ))+#scale_linewidth(range=c(1,2))+scale_size( range = c(1.5, 4.5))+
  geom_point(size=3,alpha=1,aes(col = Treatlab, shape=Composition,group=RepID)) +
  facet_wrap(~Composition)+theme_classic(base_size=26)+
  scale_shape_manual(values=c("triangle","circle","square","diamond"))+
  scale_color_jco(name="Treatment")+  scale_fill_jco(name="Treatment")+
  scale_x_continuous(labels = c(0,0.1,0.5,1), breaks=sqrt(c(0,0.1,0.5,1)) )+
  labs(y="Cancer growth rate \n (relative to ribociclib treatment \n monoculture mean of IL-15 control)", x="IL-15 (ngmL)")+ 
  theme(aspect.ratio=1)+
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        strip.background = element_blank(),
        strip.text= element_blank(),
        legend.position = "none")
#ggsave(file=paste0(svloc0,"BLANK_TcellMacrophageTricultureRGR_IL15_ribo_FitRelativegrowth.png"), dpi=320,width=8,height=8)





ggplot(allXanalysisTermsB,
       aes( y= RGR-predsNorm, x = sqrt_IL15, col = Treatlab,fill=Treatlab, 
            group=interaction(Treatlab,Composition,
                              CellLine,Lineagelab)
       )) +
  geom_ribbon(data=allpredsddB ,col=NA,alpha=0.2,aes(y=preds-predsNorm,ymax=ucl-predsNorm,ymin=lcl-predsNorm))+
  geom_line(data=allpredsddB ,aes(linetype=Cellcatlab,y=preds-predsNorm))+#scale_linewidth(range=c(1,2))+scale_size( range = c(1.5, 4.5))+
  geom_point(size=3,alpha=1,aes(col = Treatlab, shape=Cellcatlab,group=RepID)) +
  facet_grid(
    #CellLine~Lineagelab
    CellLine~Composition, scale="free")+
  theme_classic(base_size=26)+
  scale_color_jco(name="Treatment")+  scale_fill_jco(name="Treatment")+
  #scale_x_continuous(labels = unique(allXanalysisTerms$sqrt_IL15)^2, breaks=unique(allXanalysisTerms$sqrt_IL15) )+
  scale_x_continuous(labels = c(0,0.5,1,2,4), breaks=sqrt(c(0,0.5,1,2,4)) )+
  #scale_x_continuous(labels = c(1000,2000,4000,8000,16000,32000,64000), breaks=log10(c(1000,2000,4000,8000,16000,32000,64000)))+
  labs(y="Cancer growth rate", x="IL-15 (ngmL)")+ 
  theme(aspect.ratio=1)



ggplot(allXanalysisTerms,
       aes( y= RGR, x = sqrt_IL15, col = Treatlab,fill=Treatlab, 
            group=interaction(Treatlab,Composition,
                              CellLine,Lineagelab)
       )) +
  geom_ribbon(data=allpredsdd ,alpha=0.2,col=NA,aes(y=preds,ymax=ucl,ymin=lcl))+
  geom_line(data=allpredsdd ,aes(y=preds))+#scale_linewidth(range=c(1,2))+scale_size( range = c(1.5, 4.5))+
  geom_point(alpha=1,aes(col = Treatlab,group=RepID)) +
  facet_wrap(~Composition)+theme_classic(base_size=26)+
  scale_color_jco(name="Treatment")+  scale_fill_jco(name="Treatment")+
  #scale_x_continuous(labels = unique(allXanalysisTerms$sqrt_IL15)^2, breaks=unique(allXanalysisTerms$sqrt_IL15) )+
  scale_x_continuous(labels = c(0,0.1,0.5,1), breaks=sqrt(c(0,0.1,0.5,1)) )+
  #scale_x_continuous(labels = c(1000,2000,4000,8000,16000,32000,64000), breaks=log10(c(1000,2000,4000,8000,16000,32000,64000)))+
  labs(y="Cancer growth rate", x="IL-15 (ngmL)")+ 
  theme(aspect.ratio=1)+
  scale_size_manual(values=c(4,2))+
  theme(axis.text.x=element_text(size=rel(0.9)))

ggsave(file=paste0(svloc0,"TcellMacrophageTricultureRGR_IL15_ribo_Fit.png"), dpi=320,width=14,height=14)

#4 cell lines
#6 IL-15 concentrations
#2 culture compositions
#2 ribocilib doses
#3 replicates

unique(allXanalysisTerms%>%dplyr::select(CellLine,Lineagelab,FacTreatlab,sqrt_IL15,Treatlab,Composition))
unique(allXanalysisTerms%>%dplyr::select(CellLine,Lineagelab,FacTreatlab,sqrt_IL15,Treatlab,Composition))$CellLine%>%table()
unique(allXanalysisTerms%>%dplyr::select(CellLine,Lineagelab,FacTreatlab,sqrt_IL15,Treatlab,Composition))$sqrt_IL15%>%table()

ggplot(allXanalysisTerms,
       aes( y= RGR, x = sqrt_IL15, col = Treatlab,fill=Treatlab, 
            group=interaction(Treatlab,Composition,
                              CellLine,Lineagelab)
       )) +
  geom_ribbon(data=allpredsdd ,alpha=0.2,col=NA,aes(y=preds,ymax=ucl,ymin=lcl))+
  geom_line(data=allpredsdd ,aes(y=preds))+#scale_linewidth(range=c(1,2))+scale_size( range = c(1.5, 4.5))+
  geom_point(alpha=1,aes(col = Treatlab,group=RepID)) +
  facet_grid(CellLine~Lineagelab, scale="free")+theme_classic(base_size=26)+
  scale_color_jco(name="Treatment")+  scale_fill_jco(name="Treatment")+
  #scale_x_continuous(labels = unique(allXanalysisTerms$sqrt_IL15)^2, breaks=unique(allXanalysisTerms$sqrt_IL15) )+
  scale_x_continuous(labels = c(0,0.5,1,2,4), breaks=sqrt(c(0,0.5,1,2,4)) )+
  #scale_x_continuous(labels = c(1000,2000,4000,8000,16000,32000,64000), breaks=log10(c(1000,2000,4000,8000,16000,32000,64000)))+
  labs(y="Cancer growth rate", x="IL-15 (ngmL)")+ 
  theme(aspect.ratio=1)+
  scale_size_manual(values=c(4,2))+
  #theme(axis.text.x=element_text(size=rel(0.9)))+
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        strip.background = element_blank(),
        strip.text= element_blank(),
        legend.position = "none")

# stats
data.table(summary(allgammods)$s.table,keep.rownames = T)
data.table(summary(allgammods)$p.table,keep.rownames = T)

summary(lm(RGR~ (RiboTreated+ sqrt_IL15+isTCoculture+isMCoculture)^2+
             RiboTreated:sqrt_IL15:isTCoculture +RiboTreated:sqrt_IL15:isMCoculture ,
           data=cellXanalysis))

ss1<-data.table(summary(allgammods)$s.table,keep.rownames = T)
data.table(summary(allgammods)$p.table,keep.rownames = T)
 
#(IL-15 drives decreased cancer growth (increased cancer killing effect) in T coculture:CAMA-1 resistant:edf=1.92, F=18.14,p=2.62e-7
#(IL-15 drives decreased cancer growth in M coculture:CAMA-1 resistant:edf=1.27, F=9.50,p=2.46e-5
# NS effect of IL-15 on cancer monoculture :edf=6.686e-05, F=0.000,p=0.75
# Ribo Treatment decreases growth of monocultures: est=-6.51,t=-6.92,p= 9.62e-10
# Greater effect of Ribo Treatment to decrease cancer growth when cocultured with macropahges: est=-5.34,t=-3.69,p= 0.00041

#other Higher order interactions
#s(isTCoculture_sqrt_IL15)               1.924e+00      2 18.139 2.62e-07 ***
 # s(isMCoculture_sqrt_IL15)               1.269e+00      2  9.496 2.46e-05 ***
 # s(Ribo__isTCoculture_sqrt_IL15)         1.885e+00      2  7.804 0.000346 ***
 # s(Ribo__isMCoculture_sqrt_IL15)         9.480e-01      2 12.216 7.24e-07 ***

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
  #scale_x_continuous(labels = unique(allXanalysisTerms$sqrt_IL15)^2, breaks=unique(allXanalysisTerms$sqrt_IL15) )+
  scale_x_continuous(labels = c(0,0.5,1,2,4), breaks=sqrt(c(0,0.5,1,2,4)) )+
  
  labs(y="Ribociclib effect on \n cocultured cancer growth rate", x="IL-15 (ngmL)")+ 
  theme(aspect.ratio=1)+
  scale_color_npg(name="Ribociclib \n resistance")

ggsave(file=paste0(svloc,"TcellcocultureRGR_IL15_ribo_RiboEffect.png"), dpi=320,width=9,height=9)


