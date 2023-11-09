rm(list=ls())
# install.packages("lmerTest")
require(lmerTest)
require(dplyr)
require(data.table)
require(ggplot2)
require(tidyr)

ddfile<-"/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/doi_10.5061_dryad.03n60__v1/DS2/DS2_datafile.tsv"
#Gray_data_raw_dose_response.csv
dd0 <- fread(ddfile)
names(dd0) <- gsub(" ","_",names(dd0))
names(dd0)[6] <-"Treatment"

dd0[,is.Focal:=F]
dd0[Cell_Name%in%c("CAMA-1","T47D","MCF7"),is.Focal:=T]

# Recreate a version of fig2
# ggplot( data.table( dd0[is.finite(GR50)]  %>%group_by(  Treatment)%>%
#                       summarise(
#                         GR50=median(GR50,na.rm=T),
#                         GRmax=median(GRmax,na.rm=T)) )[GR50<100] ,aes(y=GRmax,x=log10(GR50))
#         
# )+ geom_point()


# Select experimental results for cells treated with "(Z)-4-Hydroxytamoxifen"
dd <- dd0[Perturbagen_Target== "ESR1" ]#[Treatment=="(Z)-4-Hydroxytamoxifen"]

# Average across replicates of the same cell line
dd1<-data.table( dd%>%group_by(Cell_HMS_LINCS_ID , Cell_Name,is.Focal,Clinical_Subtype,Transcriptional_Subtype,Treatment,Perturbagen_Target,Conc_Unit)%>%
  summarise(
    GR50=median(GR50,na.rm=T),
    GEC50=median(GEC50,na.rm=T),
    GEC50=median(GEC50,na.rm=T),
    GRmax=median(GRmax,na.rm=T),
    GRinf=median(GRinf,na.rm=T),
  
    GR_Hill_Coefficient=median(GR_Hill_Coefficient,na.rm=T),
    GR_AOC=median(GR_AOC,na.rm=T),
    GR_r2=median(GR_r2,na.rm=T),
    Nominal_Division_Rate=median(Nominal_Division_Rate,na.rm=T)) )
 
#non-malignant (NM) and hormone-receptor-positive cell lines (HR+)

dd1[is.Focal==T]$GR_AOC
mean(dd1[!Clinical_Subtype%in%c("TNBC","NM")][Transcriptional_Subtype=="Luminal"]$GR_AOC)
dd1$is.Focal <- factor(dd1$is.Focal, levels= c("TRUE","FALSE"))


ggplot(dd1[], aes(x=GR_AOC,fill=is.Focal))+ geom_histogram(binwidth=0.01)+ theme_classic()+
  scale_fill_aaas(name="Cell line", labels=c("Focal","Non-focal"))+facet_wrap(~Clinical_Subtype)

ggplot(dd1[!Clinical_Subtype%in%c("TNBC","NM")][Transcriptional_Subtype=="Luminal"], aes(x=GR_AOC,fill=is.Focal))+ geom_histogram(binwidth=0.07)+ theme_classic(base_size=28)+
  scale_fill_aaas(name="Cell line", labels=c("Focal","Non-focal"))+
  geom_vline(xintercept =mean(dd1[!Clinical_Subtype%in%c("TNBC","NM")][Transcriptional_Subtype=="Luminal"]$GR_AOC) )+
  labs(x="Endocrine sensitivity", y="Frequency")+
  theme(aspect.ratio = 1)
ggsave(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/CAMAMCFT47DEndocrineresistance.png", dpi=320, height=8,width=8)

  

