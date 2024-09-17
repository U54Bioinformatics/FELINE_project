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
dd0 <- data.table(read.csv( file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/Figure2/SourceData_Figure2_CellAnnotations.csv"))

# count of each cell type's abundance and frequency
countTable0 <- data.table( dd0 %>% group_by(Patient.Study.ID,Day,Cohort,Treatmentlab,dynamic_class3,
                                            CellAnnot2,totalCount)%>% summarise(count=n() , frac=n()/unique(totalCount) ) ) 
countTableFeline1 <- countTable0[Cohort=="Discovery"]

# reformat data rows filling missing cell types with zero values and adjusting column names 
abundDD <- data.table(countTable0 %>% dplyr::select(-count)%>%spread(CellAnnot2,frac,fill=0))
names(abundDD) <- gsub("\\+", "_", names(abundDD))
names(abundDD) <- gsub("\\-", "_", names(abundDD))
names(abundDD) <- gsub(" ", "_", names(abundDD))


#ggplot( abundDD[Cohort=="Discovery"] , aes(y=logit( 1e-5+ +Tregs+CD4__T_cells+CD8__T_cells+NK_cells ),x=log(1+Day),col=dynamic_class3))+geom_point() +geom_smooth(method="lm",se=F)

eps <-1e-5

ggplot(abundDD[Cohort=="Discovery"], aes(y=logit(eps+Tregs+CD4__T_cells+CD8__T_cells+NK_cells), x=log(1+Day),group=interaction(Day,dynamic_class3), fill=dynamic_class3  ) )+  theme_classic(base_size=26) +  #,fill=as.factor(archetype)
  geom_boxplot(col="black") + geom_point(size=3,position = position_dodge(width=2)) + facet_wrap(~Treatmentlab, ncol=1) +
  scale_fill_npg(name="Tumor response", labels=c("Resistant","Sensitive"))  +
  theme(aspect.ratio=1)+#, legend.position = "none")+
  labs(y="T cell frequency", x= "Day")+
  scale_y_continuous(breaks=logit(eps+c(0.001,0.01,0.1) ), labels=c(0.001,0.01,0.1) ) +
  scale_x_continuous(breaks=log(1+c(0,14,180) ), labels=c(0,14,180) )
paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Ribo and Letrozole T cell Composition by response over time version2.png"),height=10,width=10)

ggplot(abundDD[], aes(y=logit(eps+Tregs+CD4__T_cells+CD8__T_cells+NK_cells), x=log(1+Day),group=interaction(Day,dynamic_class3), fill=dynamic_class3 ,shape= Cohort) )+  theme_classic(base_size=26) +  #,fill=as.factor(archetype)
  geom_boxplot(col="black") + geom_point(size=3,position = position_dodge(width=2)) + facet_wrap(~Treatmentlab, ncol=1) +
  scale_fill_npg(name="Tumor response", labels=c("Resistant","Sensitive"))  +
  theme(aspect.ratio=1)+#, legend.position = "none")+
  labs(y="T cell frequency", x= "Day")+
  scale_y_continuous(breaks=logit(eps+c(0.001,0.01,0.1) ), labels=c(0.001,0.01,0.1) ) +
  scale_x_continuous(breaks=log(1+c(0,14,180) ), labels=c(0,14,180) )
#ggsave(paste0(paperfile,"Ribo and Letrozole T cell Composition by response over time version2 Discovery and Validation.png"),height=10,width=10)

savloc<-"/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/Figure5/Outputs/"
write.csv(abundDD, file=paste0(savloc,"SourceData_Figure5_TcellFrequency_Output.csv"))

abundDD[,success:= round((Tregs+CD4__T_cells+CD8__T_cells+NK_cells)*totalCount)]
abundDD[,fail:= totalCount-success]


summary( glm(cbind(success,fail)   ~0+dynamic_class3*as.factor(Day) + Cohort , family="binomial", abundDD[Treatmentlab=="Combination ribociclib"]) )
summary( glm(cbind(success,fail)   ~0+dynamic_class3*as.factor(Day) + Cohort, family="binomial", abundDD[Treatmentlab=="Letrozole alone"]) )

summary( glm(cbind(success,fail)   ~dynamic_class3*as.factor(Day) + Cohort , family="binomial", abundDD[Treatmentlab=="Combination ribociclib"]) )

summary( glm(cbind(success,fail)   ~dynamic_class3*as.factor(Day) + Cohort, family="binomial", abundDD[Treatmentlab=="Letrozole alone"]) )

# 
# summary( glm(cbind(success,fail)   ~Treatmentlab*dynamic_class3 , family="binomial", umapout[Day==180]) )
# summary( glm(cbind(success,fail)   ~Treatmentlab*dynamic_class3 , family="binomial", umapout[Day==180]) )
# 
# summary( glm(cbind(success,fail)   ~dynamic_class3 , family="binomial", abundDD[Day==0][Treatmentlab=="Combination ribociclib"]) )
# summary( glm(cbind(success,fail)   ~dynamic_class3 , family="binomial", abundDD[Day==14][Treatmentlab=="Combination ribociclib"]) )
# summary( glm(cbind(success,fail)   ~dynamic_class3 , family="binomial", abundDD[Day==180][Treatmentlab=="Combination ribociclib"]) )
# 
# summary( glm(cbind(success,fail)   ~dynamic_class3+Cohort , family="binomial", abundDD[Day==0][Treatmentlab=="Letrozole alone"]) )
# summary( glm(cbind(success,fail)   ~dynamic_class3+Cohort , family="binomial", abundDD[Day==14][Treatmentlab=="Letrozole alone"]) )
# summary( glm(cbind(success,fail)   ~dynamic_class3+Cohort , family="binomial", abundDD[Day==180][Treatmentlab=="Letrozole alone"]) )

