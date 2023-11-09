rm(list=ls()); require(data.table);require(dplyr);require(ggplot2);require(tidyr);
library(colorspace)
require(segmented)
require(ggsci)
# Load the data
dd <- data.table( read.csv(file="~/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Tcell_Experiments/Exp52/RawDataQC_PC/Exp52_Plate152-155_MDAMB134_Parsed_13Kcutoff.csv"))

dd[,RepID:=1:nrow(dd)]
dd[,Rep:=1:length(RepID), by=c("Plate","CancerCell","Resistance","Fraction_Tcell","TCellPrep","Ribociclib_uM","IL15_ngmL")]
dd$CancerCell%>%unique()
dd[CancerCell=="MDAMB134_SensV2",CancerCell:="MDAMB134_Sens_V2"]
dd[CancerCell=="MDAMB134_RiboRCer2",CancerCell:="MDAMB134_RiboR_Cer2"]

# First we'll continue to work with florescence measured total intensity; filtering ATP conten
dat2 <- dd %>% dplyr::select(c("RepID","Rep","CancerCell":"IL15_ngmL",
                               names(dd) [c( grep("Sens_Intensity_",names(dd)), grep("Res_Intensity_",names(dd)))]
                               #"Sens_Intensity_0hr":"Sens_Intensity_146hr", "Res_Intensity_0hr":"Res_Intensity_146hr")
) )


# Convert to Long format 
dat_long <- data.table(gather(dat2, Time, Measurement, 
                              names(dd) [c( grep("Sens_Intensity_",names(dd)), grep("Res_Intensity_",names(dd)))]
))
dat_long$Ribociclib_uM <- as.factor(dat_long$Ribociclib_uM)
dat_long$IL15_ngmL <- as.factor(dat_long$IL15_ngmL)

# Split one column with 3 pieces of information into 3 different columns
dat_long[, c("CellLine", "Sensitivity","Fluoresence") := tstrsplit(CancerCell, "_", fixed=TRUE)]
dat_long[,Fluoresence:=NULL]

# get numeric time
dat_long[, Hour := as.numeric( gsub("[^0-9]","",Time) )]
# adjust labeling to match Fengs standard curve data
dat_long[, Variable:="Sens_Intensity" ]
dat_long[grep("Res_Intensity",Time), Variable:="Res_Intensity" ]
dat_long[,Time:=NULL]
setnames(dat_long,old="Measurement",new="Value")

# As we are dealing with monocultures there is no fluoresence if R cells in S culture and vica versa. 
dat_long[(Sensitivity=="Sens"&Variable=="Res_Intensity" ),Value:=0 ]
dat_long[(Sensitivity=="RiboR"&Variable=="Sens_Intensity" ),Value:=0 ]

# Translate intensity to counts
inputdd <- data.table(dat_long%>%spread(Variable,Value))

# load regression model fits obtained from Feng's standard curve mapping data (translate intensity to counts )
getCountModels <- function(){
  # Load the data
  dd <- data.table( read.csv(file="~/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Tcell_Experiments/Exp53/RawDataQC_PC/Exp53_MDAMB134_LY2_CellNumber_Parsed.csv") )
  
  # Additional columns to subset either by cell state (fibrobalst sensitive or resistant cancer) or cell type (lineage name)
  dd[grep("Sensitive", Resistance), Celltype := "Sens"]
  dd[grep("Ribo_Resistant",Resistance), Celltype := "RiboR"]
  # Add column labelling different cancer cell types
  dd[, CelltypeName :=  "MDAMB134"]
  dd[grep("LY2",CancerCell), CelltypeName := "LY2"]
  dd[CancerCell=="MDAMB134_SensV2",CancerCell:="MDAMB134_Sens_V2"]
  # Alter intensity varaible names to match that of Feng's full dataset
  names(dd) <- gsub("_24hr","", names(dd))
  names(dd)[names(dd)=="Res_Int"] <-  "RiboR_Intensity"
  names(dd)[names(dd)=="Sens_Int"] <-  "Sens_Intensity"
  names(dd)[names(dd)=="CellNumber"] <-  "Cell_Number"
  
  # Use linear regression to predict cell number from fluorescent intensity
  SensLm <- lm(Cell_Number~ -1 + (Sens_Intensity):CancerCell,   data = dd[Celltype =="Sens"] ) # plot(SensLm) # coef(SensLm)
  ResLm <- lm(Cell_Number~ -1 + (RiboR_Intensity):CancerCell,   data = dd[Celltype =="RiboR"] )#  plot(ResLm) # coef(ResLm)
  return(list("SensLm"=SensLm,"ResLm"=ResLm))
}
CountModels <- getCountModels()

# Predict counts for each cell type
names(inputdd)[names(inputdd)=="Res_Intensity"] <-  "RiboR_Intensity"
inputdd[,Sens_Count:= 0]
inputdd[,Res_Count:= 0]
inputdd[Sensitivity=="Sens", Sens_Count:= predict(CountModels$SensLm,newdata=inputdd[Sensitivity=="Sens"])]
inputdd[Sensitivity=="RiboR", Res_Count:= predict(CountModels$ResLm,newdata=inputdd[Sensitivity=="RiboR"])]

# remove intensity now that we have counts
inputdd[,RiboR_Intensity :=NULL]
inputdd[,Sens_Intensity:=NULL]

# create a single total count column instead of one count column per cell type and remove other redundant data from object 
inputdd[,Count:=Sens_Count+Res_Count ]
inputdd[,Sens_Count:=NULL]
inputdd[,Res_Count:=NULL]
inputdd[,CancerCell:=NULL]

# create a log count variable
inputdd[,lnCount:= log(Count)]

# Calculate Initial Tcell Fraction (note  Ratio_Cancer_Tcell==0 indicates 100% cancer not 0 )
inputdd[,InitialTcellFraction:=Fraction_Tcell]

# labels for treatment and cell lineage on plots
inputdd[,Treatlab:="Control"]
inputdd[Ribociclib_uM!=0,Treatlab:="Ribociclib"]

inputdd[,Lineagelab:="Sensitive"]
inputdd[Sensitivity=="RiboR",Lineagelab:="Resistant"]
inputdd$Sensitivity <- factor(inputdd$Sensitivity , levels=c("Sens","RiboR"))

inputdd[,Composition:="Cancer and T cells"]
inputdd[InitialTcellFraction==0,Composition:="Cancer alone"]

inputdd[,Cytokinetreatment:="No IL-15"]
inputdd[as.numeric(as.character(IL15_ngmL))>0,Cytokinetreatment:="IL-15"]

expt52dd <- inputdd#[TCellPrep != "Stained"][][][IL15_ngmL%in%c(0,1,5)]
#save(expt52dd,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Tcell_Experiments/Exp52/Modeling_JG/data/INTENSITYSubsetParsedprocessed_data.RData")
save(expt52dd,file="~/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Tcell_Experiments/Exp52/Modeling_JG/data/INTENSITYSubsetParsedprocessed_dataFinal.RData")



