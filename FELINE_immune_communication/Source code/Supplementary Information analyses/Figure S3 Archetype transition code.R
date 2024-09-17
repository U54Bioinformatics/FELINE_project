rm(list=ls())
require(data.table)
require(dplyr)
require(tidyr)
require(ggplot2)
require(boot)
require("compositions")
require(vegan)
require(ggsci)
require(mclust)

## Start by getting the M1-M2 phenotypes of myeloid cells
require(mgcv);require(lme4);require(lmerTest);require(parallel)
library(effects);require(umap);require(Rfast);require(ider);library("dendextend");library(ggdendro);require(ggsci);require(viridis);require("Rdimtools")
library(caret)
library(pROC)

# Load cell annotation data
data.loc<- "/Users/jason/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/Figure2/"
data.name <-"SourceData_Figure2_ArchetypeAnalysisOutput.csv"
allumapout <- data.table(read.csv(file=paste0(savlocFig2, data.name)))

## Calculate state transition probabilities given the pre/post treatment states of each tumor
# Extract first and last timepoint states and keep track of treatment and response of each patient's tumor
transitionsInput<-na.omit(allumapout[Day!=14]%>% 
                            mutate(DayClass = paste0("Day",Day),
                                   archetypeClass = paste0("Archetype",archetype)) %>%
                            dplyr::select(Treatmentlab,dynamic_class3,Patient.Study.ID,DayClass,archetypeClass)%>% spread(DayClass,archetypeClass))

# Explore end point frequencies of different archetypes
transitionsInput[dynamic_class3=="Response"]$Day180%>%table()
transitionsInput[dynamic_class3!="Response"]$Day180%>%table()

# Calculate transition rates between states
transitionsprobs <- data.table( data.table( table( transitionsInput%>%dplyr::select(Day0,Day180,Treatmentlab,dynamic_class3 )) )  %>%group_by(Treatmentlab,dynamic_class3,Day0) %>%mutate(prop=N/sum(N)) )
umapout[Treatmentlab=="Letrozole alone"]%>%dplyr::select(dynamic_class3,archetype)%>%table()
transitionsprobs[,From:="Immune hot \n and diverse"]
transitionsprobs[Day0=="Archetype2",From:="Fibroblast/endothlial \n enriched"]
transitionsprobs[Day0=="Archetype3",From:="Cancer \n dominated"]
transitionsprobs[,To:="Immune hot \n and diverse"]
transitionsprobs[Day180=="Archetype2",To:="Fibroblast/endothlial \n enriched"]
transitionsprobs[Day180=="Archetype3",To:="Cancer \n dominated"]
transitionsprobs[prop==0,prop:=NA]

ggplot(transitionsprobs,aes(x=From,y=To, alpha=(prop) ))+geom_tile() + theme_classic(base_size=12)+ facet_wrap(Treatmentlab~dynamic_class3) + 
  scale_alpha_continuous(name="Transition probability",na.value=0)+theme(aspect.ratio=1,axis.text.x=element_text(angle=90,vjust = 0.5))

#paperfile<- "/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication"
#ggsave(paste0(paperfile,"Ribo and Letrozole State transition probability.png"),height=8,width=8)

ggplot(transitionsprobs,aes(x=From,y=To, alpha=(prop) ))+geom_tile() + theme_classic(base_size=12)+ facet_wrap(Treatmentlab~dynamic_class3) + 
  scale_alpha_continuous(name="Transition probability",na.value=0)+theme(aspect.ratio=1,axis.text.x=element_text(angle=90,vjust = 0.5))+
  theme(axis.title=element_blank(),  axis.text=element_blank(),strip.text = element_blank(),legend.title=element_blank(),legend.text=element_blank() )
#ggsave(paste0(paperfile,"BLANK Ribo and Letrozole Composition State transition probability.png"),height=10,width=10)

#savlocS3<-"/Users/jason/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/SI data/"
#write.csv(transitionsInput, file=paste0(savlocS3,"Archetype transition data.csv"))



