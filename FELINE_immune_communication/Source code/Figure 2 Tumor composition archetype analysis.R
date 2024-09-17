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

require(mgcv);require(lme4);require(lmerTest);require(parallel)
require(effects);require(umap);require(Rfast);require(ider);require("dendextend");require(ggdendro);
require(viridis);require("Rdimtools")
require(caret)
require(pROC)


# Load cell annotation data
data.loc<- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/Figure2/"
data.name <-"SourceData_Figure2_CellAnnotations.csv"
dd0 <- data.table(read.csv( file=paste0(data.loc, data.name)))

# count of each cell type's abundance and frequency
countTable0 <- data.table( dd0 %>% group_by(Patient.Study.ID,Day,Cohort,Treatmentlab,dynamic_class3,
                                            CellAnnot2,totalCount)%>% summarise(count=n() , frac=n()/unique(totalCount) ) ) 
countTableFeline1 <- countTable0[Cohort=="Discovery"]
abundDD <- countTable0 %>% dplyr::select(-count)%>%spread(CellAnnot2,frac,fill=0)
names(abundDD) <- gsub("\\+", "_", names(abundDD))
names(abundDD) <- gsub("\\-", "_", names(abundDD))
names(abundDD) <- gsub(" ", "_", names(abundDD))

### Model training data selection and fitting
fracTable1 <- spread( countTable0[
  Cohort=="Discovery"
  ]%>%dplyr::select(-count), CellAnnot2,frac,fill=0)
fracTable1$rowID <- 1:nrow(fracTable1)

# fraction matrix without metadata. name rows by ID to retain match to metadata
fracTable2 <- as.matrix( fracTable1 %>% dplyr::select(-c("Patient.Study.ID", "Day", "Cohort","Treatmentlab" ,"dynamic_class3", "totalCount","rowID")) )
rownames(fracTable2)<- fracTable1$rowID

# adjust range
eps <- 1e-5
adjfracTable2 <- eps+ fracTable2 *(1 - 2*eps)

# Train model of compositional heterogeneity
set.seed(123) #umap version ‘0.2.3.1’ needed to get exact same dimension reduction (despite seed setting) due to algorithm being stochastic
Duinput <- logit(adjfracTable2[,!colnames(adjfracTable2)%in%c("Cancer cells","Normal_epithelial_cells")]) 
u1 <- umap::umap(logit(adjfracTable2[,!colnames(adjfracTable2)%in%c("Cancer cells","Normal_epithelial_cells")]) ,n_components=2, n_neighbors=8)#, n_neighbors=20)

# measure shannon diversity
shannonDiv <- vegan::diversity(adjfracTable2,index = "shannon")

#Pielou’s evenness
Evenness <- shannonDiv/log(specnumber(fracTable2))

# Tumors into archetypes
gmm1 <- Mclust(u1$layout, G=3, modelNames="EEV")

# Gather output of each analysis and adjust cell type names for plotting and later analysis
umapout <- data.table( cbind( fracTable1,  u1$layout,shannon=shannonDiv ,Evenness), archetype= gmm1$classification ) #,princomp(scores(comm.bc.mds))$scores
names(umapout) <- gsub("\\+", "_", names(umapout))
names(umapout) <- gsub("\\-", "_", names(umapout))
names(umapout) <- gsub(" ", "_", names(umapout))
names(umapout)
setnames(umapout, old=c("1","2"),new=c("V1","V2"))

# Explore correlation of cell type abundances with each umap dimension
cors<-data.table( t(cor(umapout%>%dplyr::select(V1,V2), umapout%>%dplyr::select(Adipocytes:Vas_Endo))) , keep.rownames =T)
cors[order(-abs(V2))]
cors[order(-abs(V1))]

# Perform Chiqu test of association of archetype state of the tumor and it's outcome
table( umapout %>% dplyr::select(archetype,dynamic_class3))%>%chisq.test()

### Construct validation dataset of all HQ cells across the discovery and validation cohorts
Vdd0 <- dd0[Quality=="HQ"]

# add treatment label column
Vdd0[,Treatmentlab:= "Combination ribociclib"]
Vdd0[ARM=="A",Treatmentlab:= "Letrozole alone"]

# rename cell annotations for subsequent analysis and make consistent with training model inputs/outputs
Vdd0[Celltype_subtype=="Endothelial cells",Celltype_subtype:="Vas-Endo"]
Vdd0$CellAnnot1<-Vdd0$Celltype
Vdd0$CellAnnot2<-Vdd0$Celltype_subtype
Vdd0[Celltype%in%c("Cancer cells","Normal epithelial cells"),CellAnnot1:="Epithelial cells"]

# Total cell count per sample= sampling effort
Vdd0[,totalCount:= length(Cell.ID) ,by=c("Patient.Study.ID","Day","Cohort","dynamic_class3","Treatmentlab")]

# Summarise cell type frequencies in each sample
VcountTable0 <- data.table( Vdd0 %>% group_by(Patient.Study.ID,Day,Cohort,Treatmentlab,dynamic_class3,
                                              CellAnnot2,totalCount)%>% summarise(count=n() , frac=n()/unique(totalCount) ) ) 

# Extract validation cohort data only
VfracTable1 <- spread( VcountTable0[      Cohort=="Validation" ]%>%dplyr::select(-count), CellAnnot2,frac,fill=0)
VfracTable1$rowID <- 1:nrow(VfracTable1)
VfracTable1$Tregs <- 0 # add back a cell type missing from validation cohort

# fraction matrix without metadata. name rows by ID to retain match to metadata
VfracTable2 <- as.matrix( VfracTable1 %>% dplyr::select(-c("Patient.Study.ID", "Day", "Cohort","Treatmentlab" ,"dynamic_class3", "totalCount","rowID")) )
rownames(VfracTable2)<- VfracTable1$rowID

# scale data as in the traning data of the discovery cohort
VadjfracTable2 <- eps + VfracTable2 *(1 - 2*eps)
Vinput <-  logit( VadjfracTable2[,!colnames(VadjfracTable2)%in%c("Cancer cells","Normal_epithelial_cells")] )[,u1$data%>%colnames()]

# Project validation cohort samples into the umpa space constructed from the training data
Vu1 <- predict(u1,data=Vinput)#colMedians(Duinput)

# Validation dataset shannon diversity
VshannonDiv <- vegan::diversity(VadjfracTable2,index = "shannon")

# Validation dataset Pielou’s evenness
VEvenness <- VshannonDiv/log(specnumber(VfracTable2))

# gather output of analysis of validation data
Vumapout <- data.table( cbind( VfracTable1,  Vu1,shannon=VshannonDiv ,Evenness=VEvenness), archetype= predict(gmm1,Vu1)$classification ) #,princomp(scores(comm.bc.mds))$scores
names(Vumapout) <- gsub("\\+", "_", names(Vumapout))
names(Vumapout) <- gsub("\\-", "_", names(Vumapout))
names(Vumapout) <- gsub(" ", "_", names(Vumapout))
setnames(Vumapout, old=c("1","2"),new=c("V1","V2"))

#### join results from the discovery (training) and validation (testing) cohorts
allumapout <- rbind( umapout,Vumapout)

# visualization
sz<-5
# Composition archetypes
ggplot(allumapout, aes(x=V1, y=V2,col=as.factor(archetype),shape=Cohort )) + geom_point(size=sz) + theme_classic(base_size=18)+theme(aspect.ratio=1)+
  scale_color_uchicago(name="Composition \n archetype", labels=c("A","B","C","D")) + #jama# jco
  labs(y="TME composition axis 2",x="TME composition axis 1" )


# Composition corelates to outcome
ggplot(allumapout[], aes(x=V1, y=V2,col=(dynamic_class3) , shape=as.factor(Day) )) + geom_point(size=sz) + theme_classic(base_size=26) +
  scale_shape_manual(name="Day",values=c(8,18,19))+
  scale_color_npg(name="Tumor response", labels=c("Resistant","Sensitive")) +
  labs(y="TME composition axis 2",x="TME composition axis 1" ) +theme(aspect.ratio=1) 

ggplot(allumapout, aes(x=V1, y=V2,col=(dynamic_class3) , shape=as.factor(Day) )) + geom_point(size=sz) + theme_classic(base_size=28) +
  scale_shape_manual(name="Day",values=c(8,18,19))+
  scale_color_npg(name="Tumor response", labels=c("Resistant","Sensitive")) +
  labs(y="TME composition axis 2",x="TME composition axis 1" ) +theme(aspect.ratio=1)+facet_wrap(~Treatmentlab, ncol=1)
#ggsave(paste0(paperfile,"Ribo and Letrozole updated LETROZOLE Composition archetypes with Tumor response by treatment.png"),height=5,width=5,dpi=320)

# Composition archetypes drivers
ggplot(allumapout, aes(x=V1, y=V2,col=logit(eps+ Cancer_cells*(1 - (100)*eps)) , shape=as.factor(archetype)) )+ geom_point(size=sz) + theme_classic(base_size=18)+
  theme(aspect.ratio=1)+labs(y="TME composition axis 2",x="TME composition axis 1" )+
  scale_color_viridis_c(name="Cancer cell \n frequency",option='C',breaks=logit(eps+c(0.005,0.05,0.5)*(1 - 100*eps)), labels=c(0.005,0.05,0.5)) +
  scale_shape_discrete(name="Composition \n archetype", labels= c( "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched", "Cancer \n dominated")) 
#ggsave(paste0(paperfile,"Ribo and Letrozole Composition archetypes with Cancer_cells frequency.png"),height=10,width=10)

ggplot(allumapout, aes(x=V1, y=V2,col=logit(eps+ Cancer_cells*(1 - (100)*eps)) , shape=as.factor(archetype)) )+ geom_point(size=1.5*sz) + theme_classic(base_size=18)+
  theme(aspect.ratio=1)+labs(y="TME composition axis 2",x="TME composition axis 1" )+
  scale_color_viridis_c(name="Cancer cell \n frequency",option='C',breaks=logit(eps+c(0.005,0.05,0.5)*(1 - 100*eps)), labels=c(0.005,0.05,0.5)) +
  scale_shape_discrete(name="Composition \n archetype", labels= c( "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched", "Cancer \n dominated"),guide=F) +
  theme(axis.title=element_blank(),  axis.text=element_blank(),strip.text = element_blank(),legend.title=element_blank(),legend.text=element_blank() )
#ggsave(paste0(paperfile,"BLANK Ribo and Letrozole Composition archetypes with Cancer_cells frequency.png"),height=10,width=10)



ggplot(allumapout, aes(x=V1, y=V2, col=logit(eps+(B_cells+Plasma_cells+CD4__T_cells+CD8__T_cells+NK_cells+Tregs+Monocytes+DC+M1_macrophages+M2_macrophages)*(1 - eps)) , shape=as.factor(archetype)) )+ geom_point(size=sz) + theme_classic(base_size=18)+
  theme(aspect.ratio=1)+labs(y="TME composition axis 2",x="TME composition axis 1" ) +
  scale_color_viridis_c(name="Immune cell \n frequency",option='B',breaks=logit(eps+c(0.001,0.01,0.1)*(1 - eps)), labels=c(0.001,0.01,0.1) ) +
  scale_shape_discrete(name="Composition \n archetype", labels= c( "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched", "Cancer \n dominated")) 
#ggsave(paste0(paperfile,"Ribo and Letrozole Composition archetypes with Immune_cells frequency.png"),height=10,width=10)

ggplot(allumapout, aes(x=V1, y=V2, col=logit(eps+(B_cells+Plasma_cells+CD4__T_cells+CD8__T_cells+NK_cells+Tregs+Monocytes+DC+M1_macrophages+M2_macrophages)*(1 - eps)) , shape=as.factor(archetype)) )+
  geom_point(size=1.5*sz) + theme_classic(base_size=18)+
  theme(aspect.ratio=1)+labs(y="TME composition axis 2",x="TME composition axis 1" ) +
  scale_color_viridis_c(name="Immune cell \n frequency",option='B',breaks=logit(eps+c(0.001,0.01,0.1)*(1 - eps)), labels=c(0.001,0.01,0.1) ) +
  scale_shape_discrete(name="Composition \n archetype", labels= c( "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched", "Cancer \n dominated"),guide=F)  +
  theme(axis.title=element_blank(),  axis.text=element_blank(),strip.text = element_blank(),legend.title=element_blank(),legend.text=element_blank() )+theme(legend.position="none")
#ggsave(paste0(paperfile,"BLANK Ribo and Letrozole Composition archetypes with Immune_cells frequency.png"),height=5,width=5)

ggplot(allumapout, aes(x=V1, y=V2, col=logit(eps+(CD4__T_cells + CD8__T_cells + NK_cells +Tregs)*(1 - eps)) , shape=as.factor(archetype)) )+ geom_point(size=sz) + theme_classic(base_size=18)+
  theme(aspect.ratio=1)+labs(y="TME composition axis 2",x="TME composition axis 1" ) +
  scale_color_viridis_c(name="T cell \n frequency",option='B',breaks=logit(eps+c(0.001,0.01,0.1)*(1 - eps)), labels=c(0.001,0.01,0.1)) +
  scale_shape_discrete(name="Composition \n archetype", labels= c( "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched", "Cancer \n dominated")) 
#ggsave(paste0(paperfile,"Ribo and Letrozole Composition archetypes with T_cells frequency.png"),height=10,width=10)

ggplot(allumapout, aes(x=V1, y=V2, col=logit(eps+(CD4__T_cells + CD8__T_cells + NK_cells +Tregs)*(1 - eps)) , shape=as.factor(archetype)) )+ geom_point(size=1.5*sz) + theme_classic(base_size=18)+
  theme(aspect.ratio=1)+labs(y="TME composition axis 2",x="TME composition axis 1" ) +
  scale_color_viridis_c(name="T cell \n frequency",option='B',breaks=logit(eps+c(0.001,0.01,0.1)*(1 - eps)), labels=c(0.001,0.01,0.1)) +
  scale_shape_discrete(name="Composition \n archetype", labels= c( "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched", "Cancer \n dominated"),guide=F) +
  theme(axis.title=element_blank(),  axis.text=element_blank(),strip.text = element_blank(),legend.title=element_blank(),legend.text=element_blank() )+theme(legend.position="none")
#ggsave(paste0(paperfile,"BLANK Ribo and Letrozole Composition archetypes with T_cells frequency.png"),height=5,width=5)

ggplot(allumapout, aes(x=V1, y=V2, col=logit(eps+(Adipocytes+Fibroblasts+Lym_Endo+Pericytes+Vas_Endo)*(1 - eps)), shape=as.factor(archetype) ) )+ geom_point(size=sz) + theme_classic(base_size=18)+
  theme(aspect.ratio=1)+labs(y="TME composition axis 2",x="TME composition axis 1" ) +
  scale_color_viridis_c(name="Stromal cell \n frequency",option='D',breaks=logit(eps+c(0.001,0.01,0.1)*(1 - eps)), labels=c(0.001,0.01,0.1)) +
  scale_shape_discrete(name="Composition \n archetype", labels= c( "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched", "Cancer \n dominated")) 
#ggsave(paste0(paperfile,"Ribo and Letrozole Composition archetypes with Stromal_cells frequency.png"),height=10,width=10)

ggplot(allumapout, aes(x=V1, y=V2, col=logit(eps+(Adipocytes+Fibroblasts+Lym_Endo+Pericytes+Vas_Endo)*(1 - eps)), shape=as.factor(archetype) ) )+ geom_point(size=1.5*sz) + theme_classic(base_size=18)+
  theme(aspect.ratio=1)+labs(y="TME composition axis 2",x="TME composition axis 1" ) +
  scale_color_viridis_c(name="Stromal cell \n frequency",option='D',breaks=logit(eps+c(0.001,0.01,0.1)*(1 - eps)), labels=c(0.001,0.01,0.1)) +
  scale_shape_discrete(name="Composition \n archetype", labels= c( "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched", "Cancer \n dominated"),guide=F) +
  theme(axis.title=element_blank(),  axis.text=element_blank(),strip.text = element_blank(),legend.title=element_blank(),legend.text=element_blank() )+theme(legend.position="none")
#ggsave(paste0(paperfile,"BLANK Ribo and Letrozole Composition archetypes with Stromal_cells frequency.png"),height=10,width=10)

ggplot(allumapout, aes(x=V1, y=V2, shape=as.factor(archetype), col=shannon ) )+ geom_point(size=sz)+ theme_classic(base_size=18) +  #,fill=as.factor(archetype)
  theme(aspect.ratio=1)+labs(y="TME composition axis 2",x="TME composition axis 1" ) +
  scale_color_viridis_c(name="Shannon \n Diversity")+
  scale_shape_discrete(name="Composition \n archetype", labels= c( "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched", "Cancer \n dominated")) +
  labs(y="TME composition axis 2",x="TME composition axis 1")+ theme(aspect.ratio=1)#+
#ggsave(paste0(paperfile,"Ribo and Letrozole Composition archetypes with Shannon Diversity.png"),height=10,width=10, dpi=320)

ggplot(allumapout, aes(x=V1, y=V2, shape=as.factor(archetype), col=shannon ) )+ geom_point(size=1.5*sz)+ theme_classic(base_size=18) +  #,fill=as.factor(archetype)
  theme(aspect.ratio=1)+labs(y="TME composition axis 2",x="TME composition axis 1" ) +
  scale_color_viridis_c(name="Shannon \n Diversity")+
  scale_shape_discrete(name="Composition \n archetype", labels= c( "Immune hot \n and diverse", "Fibroblast/endothlial \n enriched", "Cancer \n dominated"),guide=F) +
  labs(y="TME composition axis 2",x="TME composition axis 1")+ theme(aspect.ratio=1)+
  theme(axis.title=element_blank(),  axis.text=element_blank(),strip.text = element_blank(),legend.title=element_blank(),legend.text=element_blank() )+theme(legend.position="none")
#ggsave(paste0(paperfile,"BLANK Ribo and Letrozole Composition archetypes with Shannon Diversity.png"),height=5,width=5, dpi=320)

# construct variables for statistical analysis
allumapout[,isarchetype1:=F]
allumapout[archetype==1,isarchetype1:=T]
allumapout[,isarchetype3:=F]
allumapout[archetype==3,isarchetype3:=T]

#savlocFig2<-"/Users/jason/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/Figure2/"
#write.csv(allumapout, file=paste0(savlocFig2,"SourceData_Figure2_ArchetypeAnalysis_Output.csv"))


# Perform logistic regression analyses for data at each timepoint to test link between archetype and utcome of a tumor
summary( glm(I(dynamic_class3=="Response")~ as.factor(isarchetype1), data=allumapout[Day==0][Treatmentlab=="Combination ribociclib"], family="binomial") )
summary( glm(I(dynamic_class3=="Response")~ as.factor(isarchetype1), data=allumapout[Day==14][Treatmentlab=="Combination ribociclib"], family="binomial") )
summary( glm(I(dynamic_class3=="Response")~ as.factor(isarchetype1), data=allumapout[Day==180][Treatmentlab=="Combination ribociclib"], family="binomial") )

summary( glm(I(dynamic_class3=="Response")~ as.factor(isarchetype1), data=allumapout[Day==0][Treatmentlab=="Letrozole alone"], family="binomial") )
summary( glm(I(dynamic_class3=="Response")~ as.factor(isarchetype1), data=allumapout[Day==14][Treatmentlab=="Letrozole alone"], family="binomial") )
summary( glm(I(dynamic_class3=="Response")~ as.factor(isarchetype1), data=allumapout[Day==180][Treatmentlab=="Letrozole alone"], family="binomial") )

# assess link between immune cell abundance and archetpye using logistic reg model
summary( glm(cbind( I(totalCount*( B_cells+Plasma_cells+CD4__T_cells+CD8__T_cells+NK_cells+Tregs+Monocytes+DC+M1_macrophages+M2_macrophages ) ),totalCount )~ as.factor(isarchetype1), family="binomial", data=allumapout[][])  )
summary( glm(cbind( I(totalCount*( CD4__T_cells+CD8__T_cells+NK_cells+Tregs ) ),totalCount )~ as.factor(isarchetype1), family="binomial", data=allumapout[][])  )
# assess link between shannon diversity and archetpye using logistic reg model
summary( lm(log(shannon)~ as.factor(isarchetype1), data=allumapout[][]) )
summary( lm(log(shannon)~ as.factor(isarchetype3), data=allumapout[][]) )

# assess link between outcome and archetpye using logistic reg model (across time points)
m1<- glm(I(dynamic_class3=="Response") ~ 0 + as.factor(archetype), data=allumapout[], family="binomial")
summary(m1)
unique(predict(m1, type="response"))

# assess link between outcome and archetpye using logistic reg model (end point analysis)
m180<- glm(I(dynamic_class3=="Response")~0+as.factor(archetype), data=allumapout[Day==180], family="binomial")
summary(m180)
predict(m180,newdata=data.table(archetype=as.factor(1:3)), type="response")%>%unique()

# assess link between outcome and archetpye using logistic reg model (pre treatment analysis)
m0<- glm(I(dynamic_class3=="Response")~0+as.factor(archetype), data=allumapout[Day==0], family="binomial")
summary(m0)
predict(m0,newdata=data.table(archetype=as.factor(1:3)), type="response")%>%unique()


# Calculate state transition probabilities given the pre/post treatment states of each tumor
transitionsInput<-na.omit(allumapout[Day!=14]%>% 
                            mutate(DayClass = paste0("Day",Day),
                                   archetypeClass = paste0("Archetype",archetype)) %>%
                            dplyr::select(Treatmentlab,dynamic_class3,Patient.Study.ID,DayClass,archetypeClass)%>% spread(DayClass,archetypeClass))

transitionsInput[dynamic_class3=="Response"]$Day180%>%table()
transitionsInput[dynamic_class3!="Response"]$Day180%>%table()




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




# Generate heatmap of each cell types frequency in each tumor sample for the discovery and validation cohort 
# include tumor archetype annotations
fracTable1 <- spread( rbind(countTable0[Cohort=="Discovery"],VcountTable0[Cohort!="Discovery"]) %>%dplyr::select(-count), CellAnnot2,frac,fill=0)
fracTable1$rowID <- 1:nrow(fracTable1)
fracTable1<- data.table(merge(fracTable1, allumapout%>%dplyr::select(Patient.Study.ID,Day,archetype,M1_macrophages,shannon) , by=c("Patient.Study.ID","Day")))[order(archetype,-M1_macrophages,-shannon)]
fracTable1$archetype <- as.factor(fracTable1$archetype)

# fraction matrix without metadata. name rows by ID to retain match to metadata
fracTable2 <- as.matrix( fracTable1%>%dplyr::select(-c("Patient.Study.ID", "Day", "Cohort","Treatmentlab" ,"dynamic_class3", "totalCount","rowID","archetype","shannon","M1_macrophages")) )
rownames(fracTable2)<- fracTable1$rowID

# adjust range as in archetype analysis (add pseuodocount and rescale 0-1)
eps <- 1e-5
adjfracTable2<- eps+ fracTable2 *(1 - 2*eps)

# construct annotations for heatmap rows
rowannot <- data.frame(fracTable1 %>% dplyr::select(rowID, dynamic_class3,Day, Treatmentlab,archetype,Cohort))
colnames(rowannot)[colnames(rowannot)=="archetype"]<-"Archetype"
rownames(rowannot) <- rowannot$rowID; rowannot$rowID <- NULL
rowannot$Day <- as.factor(rowannot$Day)
rowannot$dynamic_class3[rowannot$dynamic_class3=="Response"]<- "Sensitive"
rowannot$dynamic_class3[rowannot$dynamic_class3=="Non-response"]<- "Resistant"
rowannot$Archetype<- as.character(rowannot$Archetype)
rowannot$Archetype[rowannot$Archetype=="1"]<- "Immune hot and diverse"
rowannot$Archetype[rowannot$Archetype=="2"]<- "Fibroblast/endothlial enriched"
rowannot$Archetype[rowannot$Archetype=="3"]<- "Cancer dominated"
rowannot$Archetype<- factor(rowannot$Archetype, levels=c("Immune hot and diverse","Fibroblast/endothlial enriched","Cancer dominated"))
setnames(rowannot,old=c("dynamic_class3","Treatmentlab"), new=c("Tumor response","Treatment"))

# construct annotations for heatmap rows
annotation_col = list( "Tumor response" = rev(pal_npg("nrc")(2) )   ,
                       Treatment=rev(pal_jco("default")(2))   ,
                       Day=c("khaki4","purple", "grey" ),
                       Cohort=c("black","grey"),
                       Archetype=c(pal_uchicago("default")(3)))
names(annotation_col[[1]] ) <- c("Sensitive","Resistant")
names(annotation_col[[2]] ) <- c("Combination ribociclib","Letrozole alone")
names(annotation_col[[3]] ) <- c("0","14","180")
names(annotation_col[[4]] ) <- c("Discovery", "Validation")
names(annotation_col[[5]] ) <- c("Immune hot and diverse","Fibroblast/endothlial enriched","Cancer dominated")

# Indicate cell orders to plot
cellord <-c("Cancer cells","DC","M1 macrophages","Monocytes","B cells","Plasma cells","M2 macrophages","Tregs","NK cells","CD4+ T cells","CD8+ T cells",
            "Normal epithelial cells","Fibroblasts"  ,"Pericytes","Vas-Endo","Adipocytes" ,"Lym-Endo" )

# reorder data to plot
y<- adjfracTable2[,cellord]

# update cell type names
colnames(y)[colnames(y)=="Normal epithelial cells"]<-"Diploid epithelial cells"
colnames(y)[colnames(y)=="DC"]<-"Dendritic cells"
colnames(y)[colnames(y)=="Vas-Endo"]<-"Vascular Endodothelial cells"
colnames(y)[colnames(y)=="Lym-Endo"]<-"Lymphatic Endodothelial cells"

cellord <-c("Cancer cells","Dendritic cells","M1 macrophages","Monocytes","B cells","Plasma cells","M2 macrophages","Tregs","NK cells","CD4+ T cells","CD8+ T cells",
            "Diploid epithelial cells","Fibroblasts"  ,"Pericytes","Vascular Endodothelial cells","Adipocytes" ,"Lymphatic Endodothelial cells" )

# Generate output data table
heatmapPlotData <- scale(logit(y  ))
#savlocFig2<-"/Users/jason/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/Figure2/"
#write.csv(heatmapPlotData, file=paste0(savlocFig2,"SourceData_Figure2_HeatmapTumorCompositionAnalysis_Output.csv"))

LogitCompData <- (logit(y  ))
#savlocFigS2<-"/Users/jason/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/SI data/FigureS2/"
#write.csv(LogitCompData, file=paste0(savlocFigS2,"SourceData_FigureS2_LogitTumorCompositionAnalysis_Output.csv"))

# Plot the heatmap of composition in each sample
#scale(logit( adjfracTable2[,cellord]))
p1<-pheatmap::pheatmap( scale(logit(y  )) ,annotation_row=rowannot,cutree_rows = 4,cutree_cols = 3,
                        annotation_colors = annotation_col,
                        border_color=NA,
                        treeheight_row=10,treeheight_col=6,
                        cluster_rows = F,cluster_cols = F,
                        labels_row=rep("",200),
                        cellheight=3,cellwidth=12)


