rm(list=ls())
require(data.table);require(ggplot2);require(mclust);require(png);require(EBImage)
require(dplyr);require(tidyr);require(reader );require(scales)

## Step 1 : extract the distribution of intensity values per pixel across images
# Folder path (where are the files)
filepath0 <- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Tcell_Experiments/Exp52/RawDataQC_PC/RawImages/"

# Load metadata containing the image specific peak intensity values. Reformat to have one peak intensity values per row 
meta <- data.table(read.csv("/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Tcell_Experiments/Exp52/RawDataQC_PC/Exp52_Plate152-155_MDAMB134PeakIntParsedFINAL.csv"))
meta[,RepID:=1:nrow(meta)]
meta[,Rep:=1:length(RepID), by=c("Plate","CancerCell","Resistance","Fraction_Tcell","TCellPrep","Ribociclib_uM","IL15_ngmL")]
orgA <- data.table(meta%>% 
                     dplyr::select(c("Plate":"IL15_ngmL","RepID", "Rep", names(meta) [grep("FileName",names(meta))])) %>% 
                     gather( FileHour,FileName, names(meta) [grep("FileName",names(meta))]) )

orgB <- data.table(meta %>% 
                     dplyr::select(c("Plate":"IL15_ngmL","RepID", "Rep", names(meta) [grep("PeakInt",names(meta))])) %>% 
                     gather( PeakIntHour,PeakInt, names(meta) [grep("PeakInt",names(meta))]) )

orgA[, Hour := as.numeric( gsub("[^0-9]", "", FileHour) )]
orgA[, FileHour := NULL]
orgB[, Hour := as.numeric( gsub("[^0-9]", "", PeakIntHour) )]
orgB[, PeakIntHour := NULL]

meta <- data.table( merge(orgA, orgB, by=c("Hour", names( orgA %>% dplyr::select("Plate":"IL15_ngmL","RepID", "Rep")) ) ) )[!is.na(FileName)][!is.na(PeakInt)]


# Define number of pixels to sample per image 
npix <- 1000

# Run through images collecting a sample of grweyscale intensities and re-scaling to the raw intensity values
Intens <- rbindlist(lapply(1:nrow(meta), function(i){
  cat(i);cat("  ")
  #Find image file corresponding to row i of metadata
  Fi_loc <- find.file(meta[i]$FileName , dir = filepath0, dirs =  paste0(filepath0, list.files(filepath0),"/") )
  # Read image i and get its peak intensity
  img_i <- readImage(Fi_loc)   # display(img_i)
  PeakInt_i <- meta[i]$PeakInt
  # Extract foreground, rescale intensity and collect metadata
  dd <-  data.table(meta[i], Intensity= as.vector(PeakInt_i*(img_i)), PixelID= 1:length(img_i)  )
  set.seed(1)
  if( nrow(dd)>npix ) { dd <-   dd[PixelID%in%sample(1:length(img_i),npix)]  }
  # logtransform foreground intensity (backtransformed from greyscale)
  dd[ ,lnIntensity:= log(Intensity)]
  return(dd)
}))

# Can save this Intens output if needed
save(Intens,file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Tcell_Experiments/Exp52/Modeling_JG/data/PixIntensDist.RData")


## Step 2 : explore the distribution of intensity values per pixel across images

# View distribution of intensity values. Do we see evidence of multi-modality (multiple distinct peaks)? 
ggplot(Intens, aes(x=(Intensity)  ))+ 
  geom_histogram(bins=100)+
  theme_classic()+facet_wrap(~CancerCell, scales="free_x")+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))  


ggplot(Intens, aes(x=Intensity, group=interaction(Fraction_Tcell,IL15_ngmL ,Ribociclib_uM,Hour  )))+ 
  geom_density()+
  theme_classic()+facet_grid(Hour~Resistance, scales="free_x")+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))  


# Focus on RiboR cells  and look across T cell coculture fractions
ggplot(Intens[CancerCell=="MDAMB134_RiboRCer2"][Resistance=="RiboR"], aes(x=(Intensity), group=interaction(Fraction_Tcell,IL15_ngmL ,Ribociclib_uM,Hour  )))+ 
  geom_density()+
  theme_classic()+facet_grid(Hour~paste0("T cell \n",100*Fraction_Tcell,"%"), scales="free_x")+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))  

# Examine high T cell fraction cocultures with IL-15 T cell activation and no ribo treatment
ggplot(Intens[Fraction_Tcell==0.2][IL15_ngmL==5][Ribociclib_uM==0][CancerCell=="MDAMB134_RiboRCer2"], 
       aes(x=log(Intensity) , col=as.factor(Hour),group=Hour  ))+ 
  geom_density()+
  theme_classic()+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))  




## Step 2 : identify threshold values for labelled resistant and sensitive cancer cells using the distribution of intensity values per pixel across images

## RiboR cells
# Select a threshold intensity by looking across images
Intens_RiboR <-Intens[CancerCell=="MDAMB134_RiboRCer2"][is.finite(lnIntensity)]
gmm_RiboR <-Mclust(Intens_RiboR$lnIntensity, G=4)
Intens_RiboR$Class <- as.factor(predict(gmm_RiboR)$classification)


# set cutoff based on max value of non cancer clusters  (x 1.05 to give some buffer (perhaps needs to be higher? ))
ln_cutoff_RiboR<-max(Intens_RiboR[Class!=which.max(gmm_RiboR $parameters$mean)]$lnIntensity )
Intens_RiboR[, Threshold:=exp(ln_cutoff_RiboR)]

ggplot(Intens_RiboR, aes(x=lnIntensity, fill=Class, group=Class))+ geom_histogram(bins=100)+theme_classic()+
  geom_vline(xintercept = ln_cutoff_RiboR)
ggplot(Intens_RiboR, aes(x=lnIntensity))+ geom_histogram(bins=100)+theme_classic()+
  geom_vline(xintercept = ln_cutoff_RiboR)

ggplot(Intens_RiboR, 
       aes(x=log(Intensity)  ))+ 
  geom_density(fill="grey")+
  theme_classic()+facet_wrap(~CancerCell) +
  geom_density(data=Intens_RiboR[Fraction_Tcell==0.5][IL15_ngmL==5][Ribociclib_uM==0], alpha=0.4,fill="red")+
  geom_vline(xintercept = ln_cutoff_RiboR)


ggplot(Intens_RiboR[Fraction_Tcell==0.5][IL15_ngmL==5][Ribociclib_uM==0], 
       aes(x=log(Intensity) , col=as.factor(Hour),group=Hour  ))+ 
  geom_density()+
  theme_classic()+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))  +
  geom_vline(xintercept = ln_cutoff_RiboR)

cutoff_RiboR <- exp(ln_cutoff_RiboR)




## Sens cells:
#Select a threshold intensity by looking across images
Intens_Sens <-Intens[CancerCell=="MDAMB134_SensV2"][is.finite(lnIntensity)]
gmm_Sens <-Mclust(Intens_Sens$lnIntensity, G=4)
Intens_Sens$Class <- as.factor(predict(gmm_Sens)$classification)


# set cutoff based on max value of non cancer clusters  (x 1.05 to give some buffer (perhaps needs to be higher? ))
ln_cutoff_Sens<-max(Intens_Sens[Class!=which.max(gmm_Sens $parameters$mean)]$lnIntensity )
Intens_Sens[, Threshold:=exp(ln_cutoff_Sens)]

ggplot(Intens_Sens, aes(x=lnIntensity, fill=Class, group=Class))+ geom_histogram(bins=100)+theme_classic()+
  geom_vline(xintercept = ln_cutoff_Sens)

ggplot(Intens_Sens, 
       aes(x=log(Intensity)  ))+ 
  geom_density(fill="grey")+
  theme_classic()+facet_wrap(~CancerCell) +
  geom_density(data=Intens_Sens[Fraction_Tcell==0.5][IL15_ngmL==5][Ribociclib_uM==0], alpha=0.4,fill="red")+
  geom_vline(xintercept = ln_cutoff_Sens)


ggplot(Intens_Sens[Fraction_Tcell==0.5][IL15_ngmL==5][Ribociclib_uM==0], 
       aes(x=(Intensity) , col=as.factor(Hour),group=Hour  ))+ 
  geom_density()+
  theme_classic()+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))  +
  geom_vline(xintercept = cutoff_Sens)

cutoff_Sens <- exp(ln_cutoff_Sens)



## Final threshold values
data.table(CancerCell=c("MDAMB134_RiboRCer2","MDAMB134_SensV2"), Resistance=c("RiboR","Sensitive"), Threshold=c(cutoff_RiboR ,cutoff_Sens))

rbind(Intens_RiboR,Intens_Sens)

ggplot(rbind(Intens_RiboR,Intens_Sens), aes(x=Intensity))+ geom_density(aes(fill=CancerCell),alpha=0.4)+theme_classic()+
  facet_wrap(~CancerCell,scales="free_x") +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))  +
  geom_vline(aes(xintercept = Threshold))+theme(aspect.ratio=1,legend.position = "none")+
  labs(y="Intensity", x="Probability density")


ggplot(rbind(Intens_RiboR,Intens_Sens), aes(y=Intensity, x=""))+ geom_violin(aes(fill=CancerCell),alpha=0.4)+theme_classic()+
  facet_wrap(~CancerCell,scales="free_y") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))  +
  geom_hline(aes(yintercept = Threshold))+theme(aspect.ratio=1,legend.position = "none")+
  labs(y="Intensity", x="Cell type")

#ggsave("MDAMB134 intensity thresholds.png", dpi=320)
chk <- meta[Ribociclib_uM==0][CancerCell=="MDAMB134_RiboRCer2"][Hour==max(Hour)][Fraction_Tcell==0.2][IL15_ngmL==1][Rep==1]
chk <- meta[Ribociclib_uM==0][CancerCell=="MDAMB134_RiboRCer2"][Hour==49][Fraction_Tcell==0.2][IL15_ngmL==1][Rep==1]

Fi_loc <- find.file(chk$FileName , dir = filepath0, dirs =  paste0(filepath0, list.files(filepath0), "/") )

PeakInt_i <- chk$PeakInt

# Read image i and get its peak intensity
img_i <- readImage(Fi_loc)   # 
img_strong <- img_i*(img_i>(cutoff_RiboR/PeakInt_i))
img_manual <- img_i*(img_i>(exp(9)/PeakInt_i)) 

img_otsu <- img_i*(img_i> otsu(img_i))
display( combine( list(img_i, img_strong,(img_i-img_strong) )))
display( combine( list( img_strong, img_otsu )))
display( combine( list( img_strong, img_manual )))

# display( (img_otsu-img_strong))


rm(list=ls())

load(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Tcell_Experiments/Exp52/Modeling_JG/data/PixIntensDist.RData")

IntensMDA<- Intens
load(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Tcell_Experiments/Exp51/Modeling_JG/data/PixIntensDist.RData")
IntensCAM<- Intens

Intens <- rbind(data.table(CellLine="MDA-MB-134", IntensMDA), data.table(CellLine="CAMA-1", IntensCAM))
Intens[Resistance=="RiboR",Resistance:="Ribociclib \n resistant" ]


Intens[CancerCell=="CAMA1_SensV2",cutoff:= 2000]
Intens[CancerCell=="CAMA1_RiboRCer2",cutoff:= 5550]
Intens[CancerCell=="MDAMB134_SensV2",cutoff:= 1207]
Intens[CancerCell=="MDAMB134_RiboRCer2",cutoff:= 13868]


ggplot(Intens[Intensity>0], aes(y=(Intensity) ,x="", group=Resistance,fill=Resistance  ))+ 
  geom_violin(scale ="width")+
  theme_classic(base_size=26)+facet_wrap(Resistance~CellLine, scales="free")+labs(x="Cell line",y="Pixel intensity")+
  geom_hline(aes(yintercept=cutoff))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))  +theme(aspect.ratio=1,legend.position="none")

ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/SI threshold fluorescent intensity.png", dpi=320, width=12, height=12)

ggplot(Intens[Intensity>0], aes(y=(Intensity) ,x=Resistance, group=Resistance,fill=Resistance  ))+ 
  geom_violin(scale ="width")+
  theme_classic()+facet_wrap(CellLine~Resistance, scales="free")+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))  +theme(aspect.ratio=1,legend.position="none")

thredd<-Intens[Intensity>0]
savloc<- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/SI data/FigureS25/"
write.csv(thredd, file=paste0(savloc,"SourceData_FigureS25_DistributionSpheroidImageFluorescenceIntensities.csv"))
