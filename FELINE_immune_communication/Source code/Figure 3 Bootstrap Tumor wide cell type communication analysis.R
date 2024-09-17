rm(list=ls())   
require(rTensor); require(nnTensor); library(abind)
require(data.table); require(dplyr); require(ggplot2); require(tidyr)
require(igraph)
require(RColorBrewer)
myPalette <- colorRampPalette(rev(brewer.pal(11,"Spectral")))

### Load Ligand Receptor database list of Ramilowski et al 2015
LRpairsFiltered <- data.table(read.csv( file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/SI data/FigureS4/SourceData_Figure S4 LigandReceptorPairsRamilowski2015.csv"))
LRgenelist <- unique( c(LRpairsFiltered$HPMR.Receptor, LRpairsFiltered$HPMR.Ligand) )


### Load GO database list of GF receptors
#growthFactorReceptors2 <- read.csv(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/SI data/FigureS4/Figure S4 GeneOntologyGrowthFactorReceptors.csv")$x

### Load some signalling data for downstream analysis
savloc<-"/Users/jason/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/Figure3/"
CCI_all <-  read.csv(file=paste0(savloc, "SourceData_Figure3_CellCummunicationIndexTumorWideSingnaling.csv"))
CCI_discovery <- data.table(CCI_all)[Cohort=="Discovery"]
CCI_validation <-  data.table(CCI_all)[Cohort=="Validation"]




#### Discovery cohort analysis
# generate communcation data copy to randomize 
randomizeCCI_discovery <- CCI_discovery
randomizeCCI_discovery[ , Rand_dynamic_class3:= dynamic_class3 , by= c("Treat","Pair.Name", "LigandPhenoCelltype", "Day", "ReceptorPhenoCelltype" ) ]

# calculate the average strength of communication in the actual data
summaryCCI_discovery <- data.table( CCI_discovery[][order(Treat,LigandPhenoCelltype, ReceptorPhenoCelltype ) ] %>% 
                            group_by(Treat,LigandPhenoCelltype, ReceptorPhenoCelltype, dynamic_class3, Day ) %>% 
                            dplyr::summarise(meanSig=median(log(1+scaleTransduction)))  ) 
# determine the difference in the strength of interaction between the communication of each cell type in responder and non-responder tumors
summaryCCIdiff_discovery <- spread(summaryCCI_discovery, dynamic_class3, meanSig, fill= 0 ) %>% 
  setnames(old= "Non-response", new= "Nonresponse") %>% 
  mutate(diffOfResponses= Nonresponse - Response)


# For each communication pathway, day and cell type signal senders and receiver, randomize the tumor response annotations to remove the response related network structure
# Repeat many times to determine the distribtution of communication differences that may be expected by chance under the null model that communication strength do not differ between resistant and sensitive tumors.
nRandomizations <- 1000
MultisummaryRandCCIdiff_discovery <- rbindlist( lapply(1:nRandomizations, function(x){
  randomizeCCI_discovery[ ,Rand_dynamic_class3:= sample(dynamic_class3) , by= c("Treat","Pair.Name", "LigandPhenoCelltype", "Day", "ReceptorPhenoCelltype") ]
  summaryRandCCI_discovery <- data.table( randomizeCCI_discovery[ order(Treat,LigandPhenoCelltype, ReceptorPhenoCelltype) ]%>% 
                                  group_by(Treat,LigandPhenoCelltype, ReceptorPhenoCelltype, Rand_dynamic_class3, Day) %>% 
                                  dplyr::summarise(meanSig=median(log(1+scaleTransduction)))  ) 
  summaryRandCCIdiff_discovery <- data.table( spread(summaryRandCCI_discovery, Rand_dynamic_class3, meanSig, fill= 0 ) %>%
                                      setnames(old= "Non-response", new= "Nonresponse") %>% 
                                      mutate(diffOfResponses= Nonresponse - Response),
                                    sim_id= x)
  cat(x)
  return(summaryRandCCIdiff_discovery)
}))
# clculate mean and sd of randomized network connections
Randstatistics_discovery <- data.table( MultisummaryRandCCIdiff_discovery %>% group_by(Treat,LigandPhenoCelltype, ReceptorPhenoCelltype, Day) %>% 
                                dplyr::summarise(rand.meanDif= mean(diffOfResponses) , rand.sdDif= sd(diffOfResponses)  ) )

# Combine observed and randomized statistical comparisons
combineObs_Sim_discovery <- merge(summaryCCIdiff_discovery, Randstatistics_discovery, by= c("Treat","LigandPhenoCelltype", "ReceptorPhenoCelltype", "Day"))
combineObs_Sim_discovery[ , obs.z :=  (diffOfResponses - rand.meanDif)/rand.sdDif ]

# rank observed connection strengths relative to randomized and use rank to calculate a bootstrap p value 
rankPvals_discovery <- rbindlist(lapply(1:nrow(combineObs_Sim_discovery), function(x){
  whch_set <- combineObs_Sim_discovery[x]
  rank_obs <- rank(c(combineObs_Sim_discovery[Treat==whch_set$Treat][LigandPhenoCelltype==whch_set$LigandPhenoCelltype & ReceptorPhenoCelltype==whch_set$ReceptorPhenoCelltype & Day==whch_set$Day]$diffOfResponses,
                     MultisummaryRandCCIdiff_discovery[Treat==whch_set$Treat][LigandPhenoCelltype==whch_set$LigandPhenoCelltype & ReceptorPhenoCelltype==whch_set$ReceptorPhenoCelltype & Day==whch_set$Day]$diffOfResponses))[1]
  whch_set[ , rank:= rank_obs]
  whch_set[ , pval:= (rank-1)/ (nRandomizations ) ]
  whch_set[ , runs:= nRandomizations ]
  whch_set[obs.z>0,pval:= 1 - (rank-1)/ (nRandomizations ) ]
  return(whch_set)
} ) )

# merge results with data
res1_discovery <- merge( combineObs_Sim_discovery, rankPvals_discovery %>% dplyr::select(-c(Nonresponse, Response, diffOfResponses, rand.meanDif, rand.sdDif,obs.z) ) , by=c("Treat","LigandPhenoCelltype", "ReceptorPhenoCelltype", "Day"))

# perform conservative p value adjustment for multiple comparisons
res1_discovery[, adjusted_pval := p.adjust(pval, method="holm")]

# merge final results with communication summary data
alldata_discovery <- merge(summaryCCIdiff_discovery,res1_discovery %>%dplyr::select(-c(Nonresponse,Response,diffOfResponses)),by=c("Treat","LigandPhenoCelltype", "ReceptorPhenoCelltype", "Day"))
alldata_discovery$LigandPhenoCelltype <-factor(alldata_discovery$LigandPhenoCelltype, 
                                     levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Normal epithelial cells","Cancer cells","B cells","Macrophages","CD8+ T cells","CD4+ T cells"))
alldata_discovery$ReceptorPhenoCelltype <-factor(alldata_discovery$ReceptorPhenoCelltype, 
                                       levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Normal epithelial cells","Cancer cells","B cells","Macrophages","CD8+ T cells","CD4+ T cells") )

# res1<- res1_discovery;alldata <-alldata_discovery ; MultisummaryRandCCIdiff <-MultisummaryRandCCIdiff_discovery ;summaryCCI <- summaryCCI_discovery
#save(res1,alldata,MultisummaryRandCCIdiff,nRandomizations,summaryCCI,file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Bootstrap communication results/Ribo and Letrozole Bootstrap communication results.RData")
#load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Bootstrap communication results/Ribo and Letrozole Bootstrap communication results.RData")

# Detect significant statistical effect (again exclude rare B cells as too infrequently detected in samples)
pltdd_discovery <- alldata_discovery[adjusted_pval < 0.05][LigandPhenoCelltype!= "B cells"][ReceptorPhenoCelltype!= "B cells"]

# scale z statistics for comparison
pltdd_discovery[,obs.z2:=scale(obs.z,center=F), by="Treat"]
pltdd_discovery[,obs.z2:=(obs.z-min(obs.z))/max((obs.z+min(obs.z))), by="Treat"]

# identify treatment specific effect and scale z stats
Czmax<- max(pltdd_discovery[Day!=14][Treat=="CombinationRibo"][adjusted_pval<0.05]$obs.z)
Lzmax<- max(pltdd_discovery[Day!=14][Treat!="CombinationRibo"][adjusted_pval<0.05]$obs.z)
Cd_discovery <- pltdd_discovery[Day!=14][Treat=="CombinationRibo"][adjusted_pval<0.05] #[abs(obs.z)> (0.2*Czmax)]
Ld_discovery <- pltdd_discovery[Day!=14][Treat!="CombinationRibo"][adjusted_pval<0.05][abs(obs.z)> (0.2*Lzmax)]
Ld_discovery[,obs.z2:= scale(obs.z,center=F), by="Treat"]
Cd_discovery[,obs.z2:= scale(obs.z,center=F), by="Treat"]

#rejoin data and adjust labels for plotting
plotdd_discovery <- rbind(Cd_discovery, Ld_discovery)
plotdd_discovery[, Treatmentlab:= "Combination ribociclib"]
plotdd_discovery[Treat=="LetrozoleAlone", Treatmentlab:= "Letrozole alone"]
plotdd_discovery[, Daylab:= paste0("Day ", Day) ]
plotdd_discovery[ReceptorPhenoCelltype=="Normal epithelial cells", ReceptorPhenoCelltype:="Diploid epithelial cells"]
plotdd_discovery[LigandPhenoCelltype=="Normal epithelial cells", LigandPhenoCelltype:="Diploid epithelial cells"]
plotdd_discovery[LigandPhenoCelltype=="Macrophages",LigandPhenoCelltype:= "Myeloid cells"]
plotdd_discovery[ReceptorPhenoCelltype=="Macrophages",ReceptorPhenoCelltype:= "Myeloid cells"]

# reorder for plotting
plotdd_discovery$LigandPhenoCelltype <-factor(plotdd_discovery$LigandPhenoCelltype, 
                                    levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Diploid epithelial cells","Cancer cells", "Myeloid cells","CD8+ T cells","CD4+ T cells"))
plotdd_discovery$ReceptorPhenoCelltype <-factor(plotdd_discovery$ReceptorPhenoCelltype, 
                                      levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Diploid epithelial cells","Cancer cells", "Myeloid cells","CD8+ T cells","CD4+ T cells") )


## communication heatmap of bootstrap analysis results
ggplot(plotdd_discovery,
       aes(y= ReceptorPhenoCelltype, x= LigandPhenoCelltype, fill= obs.z2)) + 
  theme_classic(base_size= 19) +
  geom_tile() + 
  facet_wrap(Treatmentlab ~ Daylab) + 
  theme(axis.text.x= element_text(angle= 90,vjust=0.5)) + 
  scale_fill_gradientn(colours=myPalette(100), 
                       name="Communication \n in resistant \n relative to \n sensitive \n tumors") + 
  theme(aspect.ratio = 1) +
  labs(y= "Signal receiver cell type", x= "Signal sender cell type")

# plotdd <- plotdd_discovery
#save(plotdd,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/temporary/DiscoveryBootstrap.RData")





#### Validation cohort analysis
# generate communcation data copy to randomize 
randomizeCCI_validation <- CCI_validation
randomizeCCI_validation[ , Rand_dynamic_class3:= dynamic_class3 , by= c("Treat","Pair.Name", "LigandPhenoCelltype", "Day", "ReceptorPhenoCelltype" ) ]

# calculate the average strength of communication in the actual data
summaryCCI_validation <- data.table( CCI_validation[order(Treat,LigandPhenoCelltype, ReceptorPhenoCelltype ) ] %>% 
                            group_by(Treat,LigandPhenoCelltype, ReceptorPhenoCelltype, dynamic_class3, Day ) %>% 
                            dplyr::summarise(meanSig=median(log(1+scaleTransduction)))  ) 
# determine the difference in the strength of interaction between the communication of each cell type in responder and non-responder tumors
summaryCCIdiff_validation <- spread(summaryCCI_validation, dynamic_class3, meanSig, fill= 0 ) %>% 
  setnames(old= "Non-response", new= "Nonresponse") %>% 
  mutate(diffOfResponses= Nonresponse - Response)

nRandomizations <- 1000
MultisummaryRandCCIdiff_validation <- rbindlist( lapply(1:nRandomizations, function(x){
  randomizeCCI_validation[ ,Rand_dynamic_class3:= sample(dynamic_class3) , by= c("Treat","Pair.Name", "LigandPhenoCelltype", "Day", "ReceptorPhenoCelltype") ]
  summaryRandCCI_validation <- data.table( randomizeCCI_validation[ order(Treat,LigandPhenoCelltype, ReceptorPhenoCelltype) ]%>% 
                                  group_by(Treat,LigandPhenoCelltype, ReceptorPhenoCelltype, Rand_dynamic_class3, Day) %>% 
                                  dplyr::summarise(meanSig=median(log(1+scaleTransduction)))  ) 
  summaryRandCCIdiff_validation <- data.table( spread(summaryRandCCI_validation, Rand_dynamic_class3, meanSig, fill= 0 ) %>%
                                      setnames(old= "Non-response", new= "Nonresponse") %>% 
                                      mutate(diffOfResponses= Nonresponse - Response),
                                    sim_id= x)
  cat(x)
  return(summaryRandCCIdiff)
}))

Randstatistics_validation <- data.table( MultisummaryRandCCIdiff_validation %>% group_by(Treat,LigandPhenoCelltype, ReceptorPhenoCelltype, Day) %>% 
                                dplyr::summarise(rand.meanDif= mean(diffOfResponses) , rand.sdDif= sd(diffOfResponses)  ) )

combineObs_Sim_validation <- merge(summaryCCIdiff_validation, Randstatistics_validation, by= c("Treat","LigandPhenoCelltype", "ReceptorPhenoCelltype", "Day"))
combineObs_Sim_validation[ , obs.z :=  (diffOfResponses - rand.meanDif)/rand.sdDif ]

rankPvals_validation <- rbindlist(lapply(1:nrow(combineObs_Sim_validation), function(x){
  whch_set <- combineObs_Sim_validation[x]
  rank_obs <- rank(c(combineObs_Sim_validation[Treat==whch_set$Treat][LigandPhenoCelltype==whch_set$LigandPhenoCelltype & ReceptorPhenoCelltype==whch_set$ReceptorPhenoCelltype & Day==whch_set$Day]$diffOfResponses,
                     MultisummaryRandCCIdiff_validation[Treat==whch_set$Treat][LigandPhenoCelltype==whch_set$LigandPhenoCelltype & ReceptorPhenoCelltype==whch_set$ReceptorPhenoCelltype & Day==whch_set$Day]$diffOfResponses))[1]
  whch_set[ , rank:= rank_obs]
  whch_set[ , pval:= (rank-1)/ (nRandomizations ) ]
  whch_set[ , runs:= nRandomizations ]
  whch_set[obs.z>0,pval:= 1 - (rank-1)/ (nRandomizations ) ]
  return(whch_set)
} ) )

res1_validation <- merge( combineObs_Sim_validation, rankPvals_validation %>% dplyr::select(-c(Nonresponse, Response, diffOfResponses, rand.meanDif, rand.sdDif,obs.z) ) , by=c("Treat","LigandPhenoCelltype", "ReceptorPhenoCelltype", "Day"))
res1_validation[, adjusted_pval := p.adjust(pval, method="holm")]

alldata_validation <- merge(summaryCCIdiff_validation,res1_validation %>%dplyr::select(-c(Nonresponse,Response,diffOfResponses)),by=c("Treat","LigandPhenoCelltype", "ReceptorPhenoCelltype", "Day"))
alldata_validation$LigandPhenoCelltype <-factor(alldata_validation$LigandPhenoCelltype, 
                                     levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Normal epithelial cells","Cancer cells","B cells","Macrophages","CD8+ T cells","CD4+ T cells"))
alldata_validation$ReceptorPhenoCelltype <-factor(alldata_validation$ReceptorPhenoCelltype, 
                                       levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Normal epithelial cells","Cancer cells","B cells","Macrophages","CD8+ T cells","CD4+ T cells") )

#alldata <- alldata_validation;res1<-res1_validation;MultisummaryRandCCIdiff<-MultisummaryRandCCIdiff_validation;summaryCCI<-summaryCCI_validation
#save(res1 ,alldata ,MultisummaryRandCCIdiff,nRandomizations,summaryCCI_validation,file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Bootstrap communication results/Ribo and Letrozole Cohort 2 Bootstrap communication results.RData")
#load(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Bootstrap communication results/Ribo and Letrozole Cohort 2 Bootstrap communication results.RData")

# collect treatment effects
pltdd_validation <- alldata_validation[adjusted_pval < 0.05][LigandPhenoCelltype!= "B cells"][ReceptorPhenoCelltype!= "B cells"]
pltdd_validation <- alldata_validation[Day!=14][adjusted_pval <=0.05][LigandPhenoCelltype!= "B cells"][ReceptorPhenoCelltype!= "B cells"]
pltdd_validation[,obs.z2:=scale(obs.z,center=F), by="Treat"]
pltdd_validation[,obs.z2:=(obs.z-min(obs.z))/max((obs.z+min(obs.z))), by="Treat"]
Czmax<- max(pltdd_validation[Day!=14][Treat=="CombinationRibo"][adjusted_pval<0.05]$obs.z)
Lzmax<- max(pltdd_validation[Day!=14][Treat!="CombinationRibo"][adjusted_pval<0.05]$obs.z)
Cd_validation <- pltdd_validation[Day!=14][Treat=="CombinationRibo"][adjusted_pval<0.05]#[abs(obs.z)> (0.2*Czmax)]
Ld_validation <- pltdd_validation[Day!=14][Treat!="CombinationRibo"][adjusted_pval<0.05]#[abs(obs.z)> (0.2*Lzmax)]
Ld_validation[,obs.z2:= scale(obs.z,center=F), by="Treat"]
Cd_validation[,obs.z2:= scale(obs.z,center=F), by="Treat"]

# re gather the treatment effects together
plotdd_validation <- rbind(Cd_validation, Ld_validation)
# add labels for plotting (as in the discovery cohort)
plotdd_validation[, Treatmentlab:= "Combination ribociclib"]
plotdd_validation[Treat=="LetrozoleAlone", Treatmentlab:= "Letrozole alone"]
plotdd_validation[, Daylab:= paste0("Day ", Day) ]
plotdd_validation[ReceptorPhenoCelltype=="Normal epithelial cells", ReceptorPhenoCelltype:="Diploid epithelial cells"]
plotdd_validation[LigandPhenoCelltype=="Normal epithelial cells", LigandPhenoCelltype:="Diploid epithelial cells"]
plotdd_validation[LigandPhenoCelltype=="Macrophages",LigandPhenoCelltype:= "Myeloid cells"]
plotdd_validation[ReceptorPhenoCelltype=="Macrophages",ReceptorPhenoCelltype:= "Myeloid cells"]

# reorder factors for plotting
plotdd_validation$LigandPhenoCelltype <-factor(plotdd_validation$LigandPhenoCelltype, 
                                    levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Diploid epithelial cells","Cancer cells", "Myeloid cells","CD8+ T cells","CD4+ T cells"))
plotdd_validation$ReceptorPhenoCelltype <-factor(plotdd_validation$ReceptorPhenoCelltype, 
                                      levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Diploid epithelial cells","Cancer cells", "Myeloid cells","CD8+ T cells","CD4+ T cells") )

#plotdd <- plotdd_validation
#save(plotdd,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/temporary/ValidationBootstrap.RData")




#### Gather data together
DiscoveryBootstrap <- data.table(Cohort="Discovery",plotdd_discovery)
ValidationBootstrap <- data.table(Cohort="Validation",plotdd_validation)
plotdd <- rbind(DiscoveryBootstrap,ValidationBootstrap)
plotdd[,Cohort2:="  "]
plotdd[Cohort=="Validation",Cohort2:=" "]

ggplot(plotdd,
       aes(y= ReceptorPhenoCelltype, x= LigandPhenoCelltype, fill= obs.z2)) + 
  theme_classic(base_size= 25) +
  geom_tile() + 
  facet_grid(Treatmentlab~paste0(Cohort2,Daylab) ) + 
  theme(axis.text.x= element_text(angle= 90,vjust=0.5)) + 
  #scale_fill_viridis_c(name= "Communication strength in \n resistant versus sensitive \n tumors",option= "A") + 
  scale_fill_gradientn(colours=myPalette(100), 
                       #limits= c(-1,1)*max(abs(plotdd[]$obs.z2)), 
                       name="Communication \n in resistant \n relative to \n sensitive \n tumors") + 
  #name="Communication strength in \n resistant relative to sensitive \n tumors") + 
  #scale_fill_distiller(palette="RdBu", name="Communication \n (Resistant relative to \n sensitive tumors)") + #limits= c(-1,1)*max(abs(pltdd$diffOfResponses)), 
  theme(aspect.ratio = 1) +
  labs(y= "Signal receiver cell type", x= "Signal sender cell type")

#save(plotdd,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/Figure3/Discovery and Validation Bootstrap Communication.RData")
#load(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/Figure3/Discovery and Validation Bootstrap Communication.RData")
#plotdd
#write.csv(plotdd,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/Figure3/SourceData_Figure3_Bootstrap Communication_Out.csv")



