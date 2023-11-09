rm(list=ls())   
#require(rTensor); require(nnTensor); library(abind)
require(data.table); require(dplyr); require(ggplot2); require(tidyr)
require(igraph)

### Phenotype data
load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/PhenotypesOfAllCellTypesAllArmsCohort2/PhenotypesOfAllCellTypesAllArmsCohort2.RData")
# Encode the subtypes that have been analysed using UMAP
tmp <- data.table(allphenotypes %>% group_by(key_) %>% slice(1)) %>% dplyr::select(Celltype , Celltype_subtype) %>% unique
tmp[ , PhenoCelltype:= Celltype]
tmp[Celltype_subtype %in% c("CD4+ T cells", "Tregs"), PhenoCelltype:= "CD4+ T cells"]
tmp[Celltype_subtype %in% c("CD8+ T cells", "NK cells"), PhenoCelltype:= "CD8+ T cells"]
tmp[Celltype_subtype %in% c("Macrophages", "DC", "Monocytes"), PhenoCelltype:= "Macrophages"]
tmp[Celltype_subtype %in% c("Vas-Endo", "Lym-Endo","Endothelial cells"), PhenoCelltype:= "Endothelial cells"]
tmp[Celltype_subtype %in% c("Cancer cells"), PhenoCelltype:= "Cancer cells"]
tmp[Celltype_subtype %in% c("Normal epithelial cells"), PhenoCelltype:= "Normal epithelial cells"]
allphenotypes <- merge(allphenotypes %>% dplyr::select(-c( paste0("V", 1:5), "nCount_RNA", "nFeature_RNA",  "Platform","Sample_p_t", "file_string", "day_fact")), tmp, by= c("Celltype", "Celltype_subtype") )
allphenotypes[Patient.Study.ID=="001-171_853686",Patient.Study.ID:="001-171"]
load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE/Feline clinical input.RData")
allphenotypes <- merge(allphenotypes, unique(Clin_resp_dd_classAdd%>%dplyr::select(Patient.Study.ID,prop_change,dynamic_class)), by="Patient.Study.ID")

### Unique clusters of cells (ALL) and their umap discretization level
uu <- unique( allphenotypes %>% dplyr::select( c("Celltype","PhenoCelltype", "key_", paste0("Disc_V", 1:5) ) ) )  #"Celltype_subtype",

### Load Ligand Receptor database list of Ramilowski et al 2015
load( "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/LigandReceptor/Filtered_Human-2015-Ramilowski-LR-pairs.RData")
LRgenelist <- unique( c(LRpairsFiltered$HPMR.Receptor, LRpairsFiltered$HPMR.Ligand) )

### Load some signalling data for downstream analysis
savelocCCI <- "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/CommunicationOutputAllArmsCohort2/"

# Are we studying the per indiv or per cell type level of signalling?
perIndiv=FALSE

# Names of LR data to analyze
filenamesCCI <- list.files(savelocCCI)#%>%length

# First let's do this for the population level
#summarytable <- read.csv("/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/TensorAnalysisOutput/Ligand Receptor list fro NTD model of population cellcell communication_AB.csv")
load(file ="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/cohort2 Communication output merged/PopulationCommunicationMerged ALL ARMS cohort2.RData" )#CCI,allphenotypes, uu,perIndiv,
CCI[,Treat:="CombinationRibo"]
CCI[ARM=="A",Treat:="LetrozoleAlone"]

CCI[,scaleTransduction:=scale(TransductionMu,center=F),by=c("Pair.Name")] #"ReceptorPhenoCelltype",
CCI[!is.finite(TransductionMu),scaleTransduction:=0]
CCI[Patient.Study.ID=="001-171_853686",Patient.Study.ID:="001-171"]
CCI <-merge(CCI, unique(Clin_resp_dd_classAdd%>%dplyr::select(Patient.Study.ID,prop_change,dynamic_class)), by="Patient.Study.ID")
CCI[ prop_change < (2/3) ,dynamic_class3:="Response"]


# First let's do this for the population level
#summarytable <- read.csv("/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/TensorAnalysisOutput/Ligand Receptor list fro NTD model of population cellcell communication_AB.csv")

CCI[,scalelnTransductionMu:=scale(log(1+TransductionMu)),by=Pair.Name]

####
# generate communcation data copy to randomize 
randomizeCCI <- CCI
randomizeCCI[ , Rand_dynamic_class3:= dynamic_class3 , by= c("Treat","Pair.Name", "LigandPhenoCelltype", "Day", "ReceptorPhenoCelltype" ) ]

# calculate the average strength of communication in the actual data
summaryCCI <- data.table( CCI[][order(Treat,LigandPhenoCelltype, ReceptorPhenoCelltype ) ] %>% 
                            group_by(Treat,LigandPhenoCelltype, ReceptorPhenoCelltype, dynamic_class3, Day ) %>% 
                            #dplyr::summarise(meanSig= mean(scalelnTransductionMu) ) )
                            dplyr::summarise(meanSig=median(log(1+scaleTransduction)))  ) 
# determine the difference in the strength of interaction between the communication of each cell type in responder and non-responder tumors
summaryCCIdiff <- spread(summaryCCI, dynamic_class3, meanSig, fill= 0 ) %>% 
  setnames(old= "Non-response", new= "Nonresponse") %>% 
  mutate(diffOfResponses= Nonresponse - Response)

nRandomizations <- 1000
MultisummaryRandCCIdiff <- rbindlist( lapply(1:nRandomizations, function(x){
  randomizeCCI[ ,Rand_dynamic_class3:= sample(dynamic_class3) , by= c("Treat","Pair.Name", "LigandPhenoCelltype", "Day", "ReceptorPhenoCelltype") ]
  summaryRandCCI <- data.table( randomizeCCI[ order(Treat,LigandPhenoCelltype, ReceptorPhenoCelltype) ]%>% 
                                  group_by(Treat,LigandPhenoCelltype, ReceptorPhenoCelltype, Rand_dynamic_class3, Day) %>% 
                                  dplyr::summarise(meanSig=median(log(1+scaleTransduction)))  ) 
  summaryRandCCIdiff <- data.table( spread(summaryRandCCI, Rand_dynamic_class3, meanSig, fill= 0 ) %>%
                                      setnames(old= "Non-response", new= "Nonresponse") %>% 
                                      mutate(diffOfResponses= Nonresponse - Response),
                                    sim_id= x)
  cat(x)
  return(summaryRandCCIdiff)
}))

Randstatistics <- data.table( MultisummaryRandCCIdiff %>% group_by(Treat,LigandPhenoCelltype, ReceptorPhenoCelltype, Day) %>% 
                                dplyr::summarise(rand.meanDif= mean(diffOfResponses) , rand.sdDif= sd(diffOfResponses)  ) )

combineObs_Sim <- merge(summaryCCIdiff, Randstatistics, by= c("Treat","LigandPhenoCelltype", "ReceptorPhenoCelltype", "Day"))
combineObs_Sim[ , obs.z :=  (diffOfResponses - rand.meanDif)/rand.sdDif ]

rankPvals <- rbindlist(lapply(1:nrow(combineObs_Sim), function(x){
  whch_set <- combineObs_Sim[x]
  rank_obs <- rank(c(combineObs_Sim[Treat==whch_set$Treat][LigandPhenoCelltype==whch_set$LigandPhenoCelltype & ReceptorPhenoCelltype==whch_set$ReceptorPhenoCelltype & Day==whch_set$Day]$diffOfResponses,
                     MultisummaryRandCCIdiff[Treat==whch_set$Treat][LigandPhenoCelltype==whch_set$LigandPhenoCelltype & ReceptorPhenoCelltype==whch_set$ReceptorPhenoCelltype & Day==whch_set$Day]$diffOfResponses))[1]
  whch_set[ , rank:= rank_obs]
  #whch_set[ , pval:= rank/ (nRandomizations +  1) ]
  whch_set[ , pval:= (rank-1)/ (nRandomizations ) ]
  whch_set[ , runs:= nRandomizations ]
  whch_set[obs.z>0,pval:= 1 - (rank-1)/ (nRandomizations ) ]
  return(whch_set)
} ) )

res1 <- merge( combineObs_Sim, rankPvals %>% dplyr::select(-c(Nonresponse, Response, diffOfResponses, rand.meanDif, rand.sdDif,obs.z) ) , by=c("Treat","LigandPhenoCelltype", "ReceptorPhenoCelltype", "Day"))
res1[, adjusted_pval := p.adjust(pval, method="holm")]
#rankPvals[Day==180][LigandPhenoCelltype=="Cancer cells"][ReceptorPhenoCelltype=="CD8+ T cells"]
#res1[Day==180][LigandPhenoCelltype=="Cancer cells"][ReceptorPhenoCelltype=="CD8+ T cells"]

alldata <- merge(summaryCCIdiff,res1 %>%dplyr::select(-c(Nonresponse,Response,diffOfResponses)),by=c("Treat","LigandPhenoCelltype", "ReceptorPhenoCelltype", "Day"))
alldata$LigandPhenoCelltype <-factor(alldata$LigandPhenoCelltype, 
                                     levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Normal epithelial cells","Cancer cells","B cells","Macrophages","CD8+ T cells","CD4+ T cells"))
alldata$ReceptorPhenoCelltype <-factor(alldata$ReceptorPhenoCelltype, 
                                       levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Normal epithelial cells","Cancer cells","B cells","Macrophages","CD8+ T cells","CD4+ T cells") )

#save(res1,alldata,MultisummaryRandCCIdiff,nRandomizations,summaryCCI,file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Bootstrap communication results/Ribo and Letrozole Cohort 2 Bootstrap communication results.RData")
load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Bootstrap communication results/Ribo and Letrozole Cohort 2  Bootstrap communication results.RData")

#alldataFull <- merge(CCI[Treat== TTT] ,res1 %>%dplyr::select(-c(Nonresponse,Response)),by=c("LigandPhenoCelltype", "ReceptorPhenoCelltype", "Day"))

#ggplot(alldataFull[Day==0][][LigandPhenoCelltype=="Cancer cells"],aes(y= (scalelnTransductionMu), x= dynamic_class3,fill=adjusted_pval<0.05)) + geom_violin(alpha=0.5) + 
#  geom_point() + facet_grid(ReceptorPhenoCelltype~dynamic_class3) + theme(axis.text.x= element_text(angle= 90))
#ggplot(alldataFull[Day==180][][LigandPhenoCelltype=="Cancer cells"],aes(y= (scalelnTransductionMu), x= dynamic_class3,fill=adjusted_pval<0.05)) + geom_violin(alpha=0.5) + 
#  geom_point() + facet_grid(ReceptorPhenoCelltype~dynamic_class3) + theme(axis.text.x= element_text(angle= 90))

#ggplot(alldataFull[Day==0],aes(y=diffOfResponses,x=LigandPhenoCelltype,fill=adjusted_pval<0.05,col=adjusted_pval<0.05))+geom_violin(scale="width",alpha=0.5)+ geom_point()+facet_wrap(~Day)+theme(axis.text.x=element_text(angle=90))
pltdd <- alldata[adjusted_pval < 0.05][LigandPhenoCelltype!= "B cells"][ReceptorPhenoCelltype!= "B cells"]

require(RColorBrewer)
myPalette <- colorRampPalette(rev(brewer.pal(11,"Spectral")))

pltdd[,obs.z2:=scale(obs.z,center=F), by="Treat"]
pltdd[,obs.z2:=(obs.z-min(obs.z))/max((obs.z+min(obs.z))), by="Treat"]

Czmax<- max(pltdd[Day!=14][Treat=="CombinationRibo"][adjusted_pval<0.05]$obs.z)
Lzmax<- max(pltdd[Day!=14][Treat!="CombinationRibo"][adjusted_pval<0.05]$obs.z)
Cd <- pltdd[Day!=14][Treat=="CombinationRibo"][adjusted_pval<0.05]#[abs(obs.z)> (0.2*Czmax)]
Ld <- pltdd[Day!=14][Treat!="CombinationRibo"][adjusted_pval<0.05][abs(obs.z)> (0.2*Lzmax)]
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
ggplot(CCI[LigandPhenoCelltype=="Cancer cells"][ReceptorPhenoCelltype== "Macrophages"][Treat=="CombinationRibo"][Day==0],#[adjusted_pval<0.001],
       aes(y= scalelnTransductionMu, x= LigandPhenoCelltype, col= dynamic_class3,group=interaction(dynamic_class3,LigandPhenoCelltype))) + 
  theme_classic(base_size= 19) +
  geom_point(position=position_dodge(width=0.5)) +   geom_violin(position=position_dodge(width=0.5)) + 

  facet_wrap(~Treat ) + 
 # theme(axis.text.x= element_text(angle= 90,vjust=0.5)) + 
  #scale_fill_viridis_c(name= "Communication strength in \n resistant versus sensitive \n tumors",option= "A") + 
  scale_color_npg( 
                       #limits= c(-1,1)*max(abs(plotdd[]$obs.z2)), 
                       name="Tumor response") + 
  #name="Communication strength in \n resistant relative to sensitive \n tumors") + 
  #scale_fill_distiller(palette="RdBu", name="Communication strength in \n resistant relative to sensitive \n tumors") + #limits= c(-1,1)*max(abs(pltdd$diffOfResponses)), 
  theme(aspect.ratio = 1) +
  labs(y= "Signal receiver cell type", x= "Signal sender cell type")+
  scale_y_continuous(breaks=log(c(0.001,0.01,0.1,1,10)),labels=c(0.001,0.01,0.1,1,10))


ggplot(plotdd[ReceptorPhenoCelltype== "Myeloid cells"][Treatmentlab=="Combination ribociclib"][Day==0],#[adjusted_pval<0.001],
       aes(y= ReceptorPhenoCelltype, x= LigandPhenoCelltype, fill= obs.z)) + 
  theme_classic(base_size= 19) +
  geom_tile() + 
  facet_wrap(Treatmentlab ~ Daylab) + 
  theme(axis.text.x= element_text(angle= 90,vjust=0.5)) + 
  #scale_fill_viridis_c(name= "Communication strength in \n resistant versus sensitive \n tumors",option= "A") + 
  scale_fill_gradientn(colours=myPalette(100), 
                       #limits= c(-1,1)*max(abs(plotdd[]$obs.z2)), 
                       name="Communication \n in resistant \n relative to \n sensitive \n tumors") + 
  #name="Communication strength in \n resistant relative to sensitive \n tumors") + 
  #scale_fill_distiller(palette="RdBu", name="Communication strength in \n resistant relative to sensitive \n tumors") + #limits= c(-1,1)*max(abs(pltdd$diffOfResponses)), 
  theme(aspect.ratio = 1) +
  labs(y= "Signal receiver cell type", x= "Signal sender cell type")



paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
paperfile<- "/Users/jason/Dropbox/FELINE Project (1)/Manuscript  Feline immune communication/Figures Communication Project/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Ribo and Letrozole Both D0and 180  Significant differences in cell type communication.png"),height=10,width=10)



ggplot(plotdd,#[adjusted_pval<0.001],
       aes(y= ReceptorPhenoCelltype, x= LigandPhenoCelltype, fill= obs.z2)) + 
  theme_classic(base_size= 19) +
  geom_tile() + 
  facet_wrap(Treatmentlab ~ Daylab) + 
  theme(axis.text.x= element_text(angle= 90,vjust=0.5)) + 
  #scale_fill_viridis_c(name= "Communication strength in \n resistant versus sensitive \n tumors",option= "A") + 
  scale_fill_gradientn(colours=myPalette(100), 
                       #limits= c(-1,1)*max(abs(plotdd[]$obs.z2)), 
                       name="Communication \n in resistant \n relative to \n sensitive \n tumors") + 
  #name="Communication strength in \n resistant relative to sensitive \n tumors") + 
  #scale_fill_distiller(palette="RdBu", name="Communication strength in \n resistant relative to sensitive \n tumors") + #limits= c(-1,1)*max(abs(pltdd$diffOfResponses)), 
  theme(aspect.ratio = 1) +
  labs(y= "Signal receiver cell type", x= "Signal sender cell type")+
  theme(axis.title=element_blank(),  axis.text=element_blank(),axis.text.x=element_blank(),strip.text = element_blank(),legend.title=element_blank(),  legend.text=element_blank())




### globe plots
require(ggraph)
require(tidygraph)
ngroup <- 1
set.seed(123);  lumrg <- data.table(Pair.Name= unique(CCI$Pair.Name) , 
                                    grpvar= sample(1:ngroup, length(unique(CCI$Pair.Name) ) ,replace= T ) )



#CCI[,scalelnTransductionMu:=scale(log(1+TransductionMu)),by=Pair.Name]
CCI[,scaleTransduction2:=exp(scale(log(TransductionMu),center=F)),by=c("Pair.Name")] 
ggplot( CCI[Pair.Name=="CSF1_CSF1R"][Treat=="CombinationRibo"][LigandPhenoCelltype!="B cells"][ReceptorPhenoCelltype!="B cells"] ,
        aes(y=scalelnTransductionMu, x=dynamic_class3  )) + geom_point() + geom_violin()+
  facet_wrap(~Day) + theme_classic() +
  theme(aspect.ratio=1)

ggplot( CCI[Pair.Name=="CSF1_CSF1R"][Treat=="CombinationRibo"][LigandPhenoCelltype!="B cells"][ReceptorPhenoCelltype!="B cells"] ,
        aes(y=scalelnTransductionMu, x=scaleTransduction  )) + geom_point() 


average_muln_scaleTransduction <- data.table(CCI%>%
                                               group_by(Patient.Study.ID,LigandPhenoCelltype,ReceptorPhenoCelltype,Day,dynamic_class3,ARM,Treat) %>%
                                               #dplyr::summarise(muln_scaleTransduction=mean(log(1+scaleTransduction))) %>% 
                                               dplyr::summarise(muln_scaleTransduction=mean(scalelnTransductionMu)) %>% 
                                               group_by(LigandPhenoCelltype,ReceptorPhenoCelltype,Day,dynamic_class3,Treat) %>%
                                               dplyr::summarise(average_muln_scaleTransduction=mean(muln_scaleTransduction)))
average_muln_scaleTransduction$LigandPhenoCelltype <-factor(average_muln_scaleTransduction$LigandPhenoCelltype, 
                                                            levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Normal epithelial cells","Cancer cells","B cells","Macrophages","CD8+ T cells","CD4+ T cells"))
average_muln_scaleTransduction$ReceptorPhenoCelltype <-factor(average_muln_scaleTransduction$ReceptorPhenoCelltype, 
                                                              levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Normal epithelial cells","Cancer cells","B cells","Macrophages","CD8+ T cells","CD4+ T cells") )
ggplot( average_muln_scaleTransduction[][Treat=="CombinationRibo"][LigandPhenoCelltype!="B cells"][ReceptorPhenoCelltype!="B cells"] ,
        aes(y=average_muln_scaleTransduction, x=dynamic_class3  )) + geom_point() + geom_violin()+
  facet_wrap(~Day) + theme_classic() +
  theme(aspect.ratio=1)


initState <- data.table(average_muln_scaleTransduction[Day==0]%>%
                          group_by(LigandPhenoCelltype,ReceptorPhenoCelltype)%>% #Treat, dynamic_class3
                          dplyr::summarise(intiState=mean(average_muln_scaleTransduction)) )

average_muln_scaleTransduction2 <- merge(average_muln_scaleTransduction,initState,by=c("LigandPhenoCelltype", "ReceptorPhenoCelltype")) #,"Treat","dynamic_class3"
average_muln_scaleTransduction2[ , lnfoldchange:= (average_muln_scaleTransduction-intiState) ] 

ggplot( average_muln_scaleTransduction2[][Treat=="CombinationRibo"][LigandPhenoCelltype!="B cells"][ReceptorPhenoCelltype!="B cells"] ,
        aes(LigandPhenoCelltype,ReceptorPhenoCelltype,fill=lnfoldchange  )) + geom_tile()+
  facet_wrap(dynamic_class3~Day)+scale_fill_viridis_c(name="Log fold change \n in communication \n relative to baseline",option="B") + theme_classic()+
  labs(y="Signal receiver",x="Signal sender")+theme(aspect.ratio=1,axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot( average_muln_scaleTransduction2[][Treat!="CombinationRibo"][LigandPhenoCelltype!="B cells"][ReceptorPhenoCelltype!="B cells"] ,
        aes(LigandPhenoCelltype,ReceptorPhenoCelltype,fill=lnfoldchange  )) + geom_tile()+
  facet_wrap(dynamic_class3~Day)+scale_fill_viridis_c(name="Log fold change \n in communication \n relative to baseline",option="B") + theme_classic()+
  labs(y="Signal receiver",x="Signal sender")+theme(aspect.ratio=1,axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


ggplot( average_muln_scaleTransduction2[][Day==180][LigandPhenoCelltype!="B cells"][ReceptorPhenoCelltype!="B cells"] ,
        aes(LigandPhenoCelltype,ReceptorPhenoCelltype,fill=lnfoldchange  )) + geom_tile()+
  facet_wrap(dynamic_class3~Treat)+scale_fill_viridis_c(name="Log fold change \n in communication \n relative to baseline",option="B") + theme_classic()+
  labs(y="Signal receiver",x="Signal sender")+theme(aspect.ratio=1,axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



ggplot( average_muln_scaleTransduction2[][Day==180][LigandPhenoCelltype!="B cells"][ReceptorPhenoCelltype!="B cells"] ,aes(LigandPhenoCelltype,ReceptorPhenoCelltype,fill=lnfoldchange  )) + geom_tile()+
  facet_wrap(dynamic_class3~Treat)+scale_fill_viridis_c(name="Log fold change \n in communication \n relative to baseline",option="B") + theme_classic()+
  labs(y="Signal receiver",x="Signal sender")+theme(aspect.ratio=1,axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))




CCIcondensed <- average_muln_scaleTransduction2[LigandPhenoCelltype!="B cells"][ReceptorPhenoCelltype!="B cells"]  
# data.table(merge( lumrg, CCI, by ="Pair.Name" ) %>% 
#                           group_by(#Patient.Study.ID, 
#                          Treat,LigandPhenoCelltype, ReceptorPhenoCelltype, dynamic_class3, Day, RiboTreated, grpvar) %>%
#                        dplyr::summarise( scaleTransduction= median( log(1+scaleTransduction) ) ) )


#exp( mean( log( scaleTransduction) ) )   


CCIcondensed$LigandPhenoCelltype <-factor(CCIcondensed$LigandPhenoCelltype, 
                                          levels= sort(c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Normal epithelial cells","Cancer cells","Macrophages","CD8+ T cells","CD4+ T cells"))) #"B cells",
CCIcondensed$ReceptorPhenoCelltype <-factor(CCIcondensed$ReceptorPhenoCelltype, 
                                            levels= sort(c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Normal epithelial cells","Cancer cells","Macrophages","CD8+ T cells","CD4+ T cells") )) #"B cells",
#setnames(CCIcondensed, old=c("LigandPhenoCelltype", "ReceptorPhenoCelltype", "scaleTransduction"), new= c("from", "to", "weight"))
setnames(CCIcondensed, old=c("LigandPhenoCelltype", "ReceptorPhenoCelltype", "lnfoldchange"), new= c("from", "to", "weight"))

# specifcy color of nodes
node_col <- ggsci::pal_npg("nrc")(2)  
# circle layout
n <- length(unique(CCIcondensed$to)) -1
pts.circle <- t(sapply(1:n, function(r)c(cos(2*r*pi/n),sin(2*r*pi/n))))
NodeList <- data.table(c("Cancer cells", "CD4+ T cells" , "Adipocytes", "Fibroblasts", "Normal epithelial cells", "Pericytes"  , "Endothelial cells" ,"Macrophages", "CD8+ T cells" )  , c(0, pts.circle[,1] )  ,  c(0, pts.circle[,2] ) ) #"B cells",
#NodeList <- data.table(c("Cancer cells", "CD4+ T cells" , "Adipocytes", "Fibroblasts", "Normal epithelial cells", "Pericytes"  , "Endothelial cells" ,"Macrophages", "B cells","CD8+ T cells" )  , c(0, pts.circle[,1] )  ,  c(0, pts.circle[,2] ) )
presloc <- data.frame(NodeList)[order(NodeList$V1),]


#tmp <- CCIcondensed[RiboTreated== T]#[Day== 0]
tmp <- CCIcondensed[Treat== "CombinationRibo"][Day== 180][order( Treat,from,to)]#[Day== 0]
#tmp <- CCIcondensed[Day== 180][order( Treat,from,to)]#[Day== 0]
tmp$DayLab <- paste0("Day ", tmp$Day )
graphdd <- as_tbl_graph( tmp, directed= T)
tmp$weight <- scale(tmp$weight)
createdLayout <- create_layout(graphdd, layout= "star")
createdLayout$x <- presloc$V2
createdLayout$y <- presloc$V3


ggraph(createdLayout , aes(x= x, y= y , col= dynamic_class3 ) ) + 
  geom_edge_loop(  aes(edge_colour=weight,edge_width= weight,edge_alpha= weight) ,  arrow= arrow(length= unit(5, "mm" ) , type= "open"),
                   start_cap= circle(12, "mm" ),  end_cap= circle(12, "mm" ) ) +
  geom_edge_arc( aes(edge_colour=weight,edge_width= weight,edge_alpha= weight  ) , arrow= arrow(length= unit(5, "mm" ) , type= "open" ),
                 start_cap= circle(12, "mm" ),  end_cap= circle(12, "mm" )  ) +
  geom_node_point( col= rep( rep(node_col, each= nrow(createdLayout) ) , each= 1) , size= 25 ) + 
  geom_node_text(aes(label= name) ) + #,col= "grey80" ) + #, size= 9) +
  facet_edges(~dynamic_class3  ) + 
  theme_void(base_size= 12 ) +
  theme(aspect.ratio = 1, legend.position = "none" ) + scale_edge_color_viridis()+
  scale_edge_width(range=c(0.0,3))+
  scale_edge_alpha(range=c(0.91,0.97)) +
  scale_edge_color_viridis(option="B")



tmp <- CCIcondensed[][Day== 0][order( Treat,from,to)]
tmp$DayLab <- paste0("Day ", tmp$Day )
tmp <- data.table( tmp %>% group_by(Treat) %>% mutate(weight = scale(weight) ) )
graphdd <- as_tbl_graph( tmp, directed= T)
createdLayout <- create_layout(graphdd, layout= "star")
createdLayout$x <- presloc$V2
createdLayout$y <- presloc$V3


ggraph(createdLayout , aes(x= x, y= y , col= dynamic_class3 ) ) + 
  geom_edge_loop(  aes(edge_colour=weight,edge_width= weight,edge_alpha= weight) ,  arrow= arrow(length= unit(3, "mm" ) , type= "closed"),
                   start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" ) ) +
  geom_edge_arc( aes(edge_colour=weight,edge_width= weight,edge_alpha= weight  ) , arrow= arrow(length= unit(3, "mm" ) , type= "closed" ),
                 start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" )  ) +
  geom_node_point( col= rep( rep(node_col, each= nrow(createdLayout) ) , each= 2) , size= 20 ) + 
  geom_node_text(aes(label= name) ) + #,col= "grey80" ) + #, size= 9) +
  facet_edges(Treat~dynamic_class3  ) + 
  theme_void(base_size= 12 ) +
  theme(aspect.ratio = 1, legend.position = "none" ) + 
  scale_edge_width(range=c(0.0,2))+
  scale_edge_alpha(range=c(0.4,0.97)) +
  scale_edge_color_viridis(option="B")






tmp <- CCIcondensed[][Day== 0][order( Treat,from,to)]
tmp$DayLab <- paste0("Day ", tmp$Day )
tmp <- data.table( tmp %>% group_by(Treat) %>% mutate(weight = scale(weight) ) )
tmp[,Treatmentlab:= "Letrozole alone"]
tmp[Treat=="CombinationRibo",Treatmentlab:= "Combination ribociclib"]
tmp[,TumorResponse:="Resistant"]
tmp[dynamic_class3=="Response",TumorResponse:="Sensitive"]
tmp[from=="Normal epithelial cells",from:="Diploid \n epithelial \n cells"]
tmp[to=="Normal epithelial cells",to:="Diploid \n epithelial \n cells"]
tmp[from=="Endothelial cells",from:="Endothelial \n cells"]
tmp[to=="Endothelial cells",to:="Endothelial \n cells"]
tmp[from=="CD8+ T cells",from:="CD8+ \n T cells"]
tmp[to=="CD8+ T cells",to:="CD8+ \n T cells"]
tmp[from=="CD4+ T cells",from:="CD4+ \n T cells"]
tmp[to=="CD4+ T cells",to:="CD4+ \n T cells"]
tmp[from=="Cancer cells",from:="Cancer \n cells"]
tmp[to=="Cancer cells",to:="Cancer \n cells"]
tmp[from=="Macrophages",from:="Myeloid \n cells"]
tmp[to=="Macrophages",to:="Myeloid \n cells"]
graphdd <- as_tbl_graph( tmp, directed= T)
createdLayout <- create_layout(graphdd, layout= "star")
createdLayout$x <- presloc$V2
createdLayout$y <- presloc$V3

p0<-ggraph(createdLayout , aes(x= x, y= y , col= dynamic_class3 ) ) + 
  geom_edge_loop(  aes(edge_colour=weight,edge_width= 5*(weight),edge_alpha= (weight)) ,  arrow= arrow(length= unit(3, "mm" ) , type= "open"),
                   start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" ) ) +
  geom_edge_arc( aes(edge_colour=weight,edge_width= (weight),edge_alpha= (weight)  ) , arrow= arrow(length= unit(3, "mm" ) , type= "open" ),
                 start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" )  ) +
  #geom_node_point( col= rep( rep(node_col, each= nrow(createdLayout) ) , each= 2) , size= 20 ) + 
  geom_node_point( col= rep( rep(node_col, each= nrow(createdLayout) ) ,  2) , size= 20 ) + 
  geom_node_text(aes(label= name),size=8 ) + #,col= "grey80" ) + #, size= 9) +
  facet_edges(Treatmentlab~TumorResponse  ) + 
  theme_void(base_size= 1.5* 14 ) +
  theme(aspect.ratio = 1, legend.position = "none" ) + 
  scale_edge_width(range=c(0.0,2))+
  scale_edge_alpha(range=c(0.0,0.97)) +
  scale_edge_color_gradient2()
#scale_edge_color_viridis(option="B")
p0 + scale_edge_color_gradient2(low="white", high="black")

#ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Globe plot ggrraph/Globe plot ggrraph1.png",height=10, width = 10)


p0<-ggraph(createdLayout , aes(x= x, y= y , col= dynamic_class3 ) ) + 
  geom_edge_loop(  aes(edge_colour=weight,edge_width= (weight),edge_alpha= (weight)) ,  arrow= arrow(length= unit(3, "mm" ) , type= "open"),
                   start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" ) ) +
  geom_edge_arc( aes(edge_colour=weight,edge_width= (weight),edge_alpha= (weight)  ) , arrow= arrow(length= unit(3, "mm" ) , type= "open" ),
                 start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" )  ,strength=0.2) +
  #geom_node_point( col= rep( rep(node_col, each= nrow(createdLayout) ) , each= 2) , size= 20 ) + 
  geom_node_point( col= rep( rep(node_col, each= nrow(createdLayout) ) ,  2) , size= 20 ) + 
  geom_node_text(aes(label= name),size=6 ) + #,col= "grey80" ) + #, size= 9) +
  facet_edges(Treatmentlab~TumorResponse  ) + 
  theme_void(base_size= 1.5* 14 ) +
  theme(aspect.ratio = 1, legend.position = "none" ) + 
  scale_edge_width(range=c(0.0,2))+
  scale_edge_alpha(range=c(0.0,0.97)) +
  scale_edge_color_gradient2()
#scale_edge_color_viridis(option="B")
p0 + scale_edge_color_gradient2(low="white", high="black")

p0<-ggraph(createdLayout , aes(x= x, y= y , col= dynamic_class3 ) ) + 
  geom_edge_loop(  aes(edge_colour=weight,edge_width= (weight),edge_alpha= (weight)) ,  arrow= arrow(length= unit(3, "mm" ) , type= "open"),
                   start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" ) ) +
  geom_edge_arc( aes(edge_colour=weight,edge_width= (weight),edge_alpha= (weight)  ) , arrow= arrow(length= unit(3, "mm" ) , type= "open" ),
                 start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" ) ,strength=0.25 ) +
  #geom_node_point( col= rep( rep(node_col, each= nrow(createdLayout) ) , each= 2) , size= 20 ) + 
  geom_node_point( col= rep( rep(node_col, each= nrow(createdLayout) ) ,  2) , size= 20 ) + 
  geom_node_text(aes(label= name),size=6, col="black" ) + #,col= "grey80" ) + #, size= 9) +
  facet_edges(Treatmentlab~TumorResponse  ) + 
  theme_void(base_size= 1.5* 18) +
  theme(aspect.ratio = 1, legend.position = "none" ) + 
  scale_edge_width(range=c(0.0,2))+
  scale_edge_alpha(range=c(0.0,0.97)) +
  scale_edge_color_gradient2()
#scale_edge_color_viridis(option="B")
p0 + scale_edge_color_gradient2(low="white",high="black")+coord_cartesian(clip="off")+theme(aspect.ratio=1,
                                                                                            plot.margin=unit(c(10,30,30,30), "points"),
                                                                                            panel.spacing=unit(2, "lines") ,
                                                                                            strip.text=element_text(size=22))

#ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Globe plot ggrraph/Globe plot ggrraph 2.png",height=12, width = 10)

paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
paperfile<- "/Users/jason/Dropbox/FELINE Project (1)/Manuscript  Feline immune communication/Figures Communication Project/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Ribo and Letrozole Globe plot ggrraph.png"),height=10,width=10)
