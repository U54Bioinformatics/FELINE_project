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

# Summarise cell type frequencies in each sample
countTable0 <- data.table( dd0 %>% group_by(Patient.Study.ID,Day,Cohort,Treatmentlab,dynamic_class3,
                                            CellAnnot2,totalCount)%>% summarise(count=n() , frac=n()/unique(totalCount) ) ) 

# Extract discovery cohort data for training
countTableFeline1 <- countTable0[Cohort=="Discovery"]

# reformat data rows filling missing cell types with zero values and adjusting column names 
abundDD <- countTable0 %>% dplyr::select(-count)%>%spread(CellAnnot2,frac,fill=0)
names(abundDD) <- gsub("\\+", "_", names(abundDD))
names(abundDD) <- gsub("\\-", "_", names(abundDD))
names(abundDD) <- gsub(" ", "_", names(abundDD))

### Full cohort all samples
# fraction table spread with cell types across columns #[Timepoint==0]
#fracTable1 <- spread( countTable0[][]%>%dplyr::select(-count), CellAnnot2,frac,fill=0)
#fracTable1$rowID <- 1:nrow(fracTable1)

# fraction matrix without metadata. name rows by ID to retain match to metadata
#fracTable2 <- as.matrix( fracTable1%>%dplyr::select(-c("Patient.Study.ID", "Day", "Cohort","Treatmentlab" ,"dynamic_class3", "totalCount","rowID")) )
#rownames(fracTable2)<- fracTable1$rowID

# Construct annotation data for legend
#rowannot <- data.frame(fracTable1 %>% dplyr::select(rowID,Cohort, Day,dynamic_class3, Treatmentlab,totalCount))
#rownames(rowannot) <- rowannot$rowID; rowannot$rowID <- NULL
#pheatmap::pheatmap( logit( 1e-10+ fracTable2*0.9999),annotation_row=rowannot)
#distmat <- vegdist(fracTable2, method = "bray")
#pheatmap::pheatmap(as.matrix(distmat),annotation_col=rowannot)

###### Data to train model for discovery cohort
# fraction table spread with cell types across columns 
fracTable1 <- spread( countTable0[
  Cohort=="Discovery"
]%>%dplyr::select(-count), CellAnnot2,frac,fill=0)
fracTable1$rowID <- 1:nrow(fracTable1)

# fraction matrix without metadata. name rows by ID to retain match to metadata
fracTable2 <- as.matrix( fracTable1 %>% dplyr::select(-c("Patient.Study.ID", "Day", "Cohort","Treatmentlab" ,"dynamic_class3", "totalCount","rowID")) )
rownames(fracTable2)<- fracTable1$rowID

# add pseudo-observation to avoid logit adjust logit(0) returning -Inf and correct range to be  within 0-1
eps <- 1e-5
adjfracTable2 <- eps+ fracTable2 *(1 - 2*eps)

# Train umap model on logit transformed data (for figure S2)
set.seed(123) 
Duinput <- logit(adjfracTable2[,!colnames(adjfracTable2)%in%c("Cancer cells","Normal_epithelial_cells")]) 
u1 <- umap::umap(logit(adjfracTable2[,!colnames(adjfracTable2)%in%c("Cancer cells","Normal_epithelial_cells")]) ,n_components=2, n_neighbors=8)#, n_neighbors=20)

shannonDiv <- vegan::diversity(adjfracTable2,index = "shannon")
#Pielou’s evenness
Evenness <- shannonDiv/log(specnumber(fracTable2))
gmm1 <- Mclust(u1$layout, G=3, modelNames="EEV")

# Assess uncertainty of archetype classification when using Training data to fit umap model on logit transformed data (for figure S2)
dfuncert<- data.table(fracTable1, uncertainty= gmm1$ uncertainty)
dfuncert$pred<-predict(gam(sqrt(uncertainty)~s(sqrt(totalCount),k=4,by=as.factor(Treatmentlab)), data=dfuncert))^2
ggplot(dfuncert,aes(y=uncertainty, x=(totalCount), col=Treatmentlab ))+
  theme_classic()+
  geom_point()+
  theme(aspect.ratio=1)+
  geom_smooth(stat="identity", aes(y=pred))+
  scale_color_jco(name= "Treatment")+
  labs(y="Tumor composition archetype \n classification uncertainty", x="Biopsy sample cell count")

savloc<- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/SI data/FigureS2/"

# Source data for fig S2
savdata<- dfuncert%>%select(Patient.Study.ID, Day, Treatmentlab,    Cohort, uncertainty,totalCount)
#write.csv(savdata, file= paste0(savloc,"TumorArchetypeUncertQuant.csv") )
