rm(list=ls())
require(data.table); require(dplyr); require(ggplot2); require(tidyr)
require(RColorBrewer)

load(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/temporary/DiscoveryBootstrap.RData")
DiscoveryBootstrap <- data.table(Cohort="Discovery",plotdd)
load(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/temporary/ValidationBootstrap.RData")
ValidationBootstrap <- data.table(Cohort="Validation",plotdd)

plotdd <- rbind(DiscoveryBootstrap,ValidationBootstrap)

myPalette <- colorRampPalette(rev(brewer.pal(11,"Spectral")))
plotdd[,Cohort2:="  "]
plotdd[Cohort=="Validation",Cohort2:=" "]

ggplot(plotdd,
       aes(y= ReceptorPhenoCelltype, x= LigandPhenoCelltype, fill= obs.z2)) + 
  theme_classic(base_size= 25) +
  geom_tile() + 
  facet_grid(Treatmentlab~paste0(Cohort2,Daylab) ) + 
  theme(axis.text.x= element_text(angle= 90,vjust=0.5)) + 
  scale_fill_gradientn(colours=myPalette(100), 
                       name="Communication \n in resistant \n relative to \n sensitive \n tumors") + 
  theme(aspect.ratio = 1) +
  labs(y= "Signal receiver cell type", x= "Signal sender cell type")

paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Discovery and Validation Ribo and Letrozole Bootstrap communication Both D0and 180  Significant differences in cell type communication.png"),height=10,width=20, dpi=320)

#plotdd
#save(plotdd,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/Figure2/Discovery and Validation Bootstrap Communication.RData")


rm(list=ls())

load(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Bootstrap communication results/Ribo and Letrozole Cohort 2 Bootstrap communication results.RData")
pltdd <- alldata[adjusted_pval < 0.05][LigandPhenoCelltype!= "B cells"][ReceptorPhenoCelltype!= "B cells"]
pltdd <- alldata[Day!=14][adjusted_pval <=0.05][LigandPhenoCelltype!= "B cells"][ReceptorPhenoCelltype!= "B cells"]
pltdd[,obs.z2:=scale(obs.z,center=F), by="Treat"]
pltdd[,obs.z2:=(obs.z-min(obs.z))/max((obs.z+min(obs.z))), by="Treat"]

Czmax<- max(pltdd[Day!=14][Treat=="CombinationRibo"][adjusted_pval<0.05]$obs.z)
Lzmax<- max(pltdd[Day!=14][Treat!="CombinationRibo"][adjusted_pval<0.05]$obs.z)
Cd <- pltdd[Day!=14][Treat=="CombinationRibo"][adjusted_pval<0.05]
Ld <- pltdd[Day!=14][Treat!="CombinationRibo"][adjusted_pval<0.05]
Ld[,obs.z2:= scale(obs.z,center=F), by="Treat"]
Cd[,obs.z2:= scale(obs.z,center=F), by="Treat"]

plotdd <- rbind(Cd, Ld)
plotdd[, Treatmentlab:= "Combination ribociclib"]
plotdd[Treat=="LetrozoleAlone", Treatmentlab:= "Letrozole alone"]
plotdd[, Daylab:= paste0("Day ", Day) ]
plotdd[ReceptorPhenoCelltype=="Normal epithelial cells", ReceptorPhenoCelltype:="Diploid epithelial cells"]
plotdd[LigandPhenoCelltype=="Normal epithelial cells", LigandPhenoCelltype:="Diploid epithelial cells"]
plotdd[LigandPhenoCelltype=="Macrophages",LigandPhenoCelltype:= "Myeloid cells"]
plotdd[ReceptorPhenoCelltype=="Macrophages",ReceptorPhenoCelltype:= "Myeloid cells"]

plotdd$LigandPhenoCelltype <-factor(plotdd$LigandPhenoCelltype, 
                                    levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Diploid epithelial cells","Cancer cells", "Myeloid cells","CD8+ T cells","CD4+ T cells"))
plotdd$ReceptorPhenoCelltype <-factor(plotdd$ReceptorPhenoCelltype, 
                                      levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Diploid epithelial cells","Cancer cells", "Myeloid cells","CD8+ T cells","CD4+ T cells") )

plotdd[ReceptorPhenoCelltype== "Myeloid cells"][Treat=="CombinationRibo"][Day==0]

ggplot(plotdd[
  ][
    ][Day%in%c(0,180)],#[adjusted_pval<0.001],
  aes(y= ReceptorPhenoCelltype, x= LigandPhenoCelltype, fill= obs.z2)) + 
  theme_classic(base_size= 19) +
  geom_tile() + 
  facet_wrap(Treatmentlab ~ Daylab) + 
  theme(axis.text.x= element_text(angle= 90,vjust=0.5)) + 
  scale_fill_gradientn(colours=myPalette(100), 
                       name="Communication \n in resistant \n relative to \n sensitive \n tumors") + 
  theme(aspect.ratio = 1) +
  labs(y= "Signal receiver cell type", x= "Signal sender cell type")

paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#paperfile<- "/Users/jason/Dropbox/FELINE Project (1)/Manuscript  Feline immune communication/Figures Communication Project/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Ribo and Letrozole Cohort 2 validation Both D0and 180  Significant differences in cell type communication.png"),height=10,width=10)


ggplot(plotdd,
       aes(y= ReceptorPhenoCelltype, x= LigandPhenoCelltype, fill= obs.z2)) + 
  theme_classic(base_size= 19) +
  geom_tile() + 
  facet_wrap(Treatmentlab ~ Daylab) + 
  theme(axis.text.x= element_text(angle= 90,vjust=0.5)) + 
  scale_fill_gradientn(colours=myPalette(100), 
                       name="Communication \n in resistant \n relative to \n sensitive \n tumors") + 
  theme(aspect.ratio = 1) +
  labs(y= "Signal receiver cell type", x= "Signal sender cell type")+
  theme(axis.title=element_blank(),  axis.text=element_blank(),axis.text.x=element_blank(),strip.text = element_blank(),legend.title=element_blank(),  legend.text=element_blank())


