rm(list=ls())   
require(data.table); require(dplyr); require(ggplot2); require(tidyr);require(igraph);require(RColorBrewer);require(ggraph);require(tidygraph)


### Phenotype data
load(file="~/Dropbox/Cancer_pheno_evo/data/FELINE2/PhenotypesOfAllCellTypesAllArms/PhenotypesOfAllCellTypesAllArms.RData")
#allphenotypes,UMAPlocs ,UMAPfiles,umapDImRedloc,umapDImRedfiles,nCellTypes,
# Encode the subtypes that have been analysed using UMAP
tmp <- data.table( allphenotypes %>% group_by(key_) %>% slice(1) ) %>% dplyr::select(Celltype , Celltype_subtype) %>% unique
tmp[ , PhenoCelltype:= Celltype]
tmp[Celltype_subtype %in% c("CD4+ T cells", "Tregs"), PhenoCelltype:= "CD4+ T cells"]
tmp[Celltype_subtype %in% c("CD8+ T cells", "NK cells"), PhenoCelltype:= "CD8+ T cells"]
tmp[Celltype_subtype %in% c("Macrophages", "DC", "Monocytes"), PhenoCelltype:= "Macrophages"]
tmp[Celltype_subtype %in% c("Vas-Endo", "Lym-Endo","Endothelial cells"), PhenoCelltype:= "Endothelial cells"]
tmp[Celltype_subtype %in% c("Cancer cells"), PhenoCelltype:= "Cancer cells"]
tmp[Celltype_subtype %in% c("Normal epithelial cells"), PhenoCelltype:= "Normal epithelial cells"]
allphenotypes <- merge(allphenotypes %>% 
  dplyr::select(-c( paste0("V", 1:5), "nCount_RNA", "nFeature_RNA", "Infercnv_CNA", "Platform","Timepoint","Sample_p_t","prop_change", "
    file_string", "day_fact" , "Ribo","TreatLab", "Burden_t0" ,"BurdenTracked" ,"Day0", "DayLastBurd", "Dose_lab","TreatCode","TreatCodeOrd","dynamic_class2","rgr_A","rgr_B")), 
  tmp, by= c("Celltype", "Celltype_subtype") )

### Unique clusters of cells (ALL) and their umap discretization level
uu <- unique( allphenotypes %>% dplyr::select( c("Celltype","PhenoCelltype", "key_", paste0("Disc_V", 1:5) ) ) )  

### Load Ligand Receptor database list of Ramilowski et al 2015
load( "~/Dropbox/Cancer_pheno_evo/data/FELINE2/LigandReceptor/Filtered_Human-2015-Ramilowski-LR-pairs.RData")
LRgenelist <- unique( c(LRpairsFiltered$HPMR.Receptor, LRpairsFiltered$HPMR.Ligand) )

### Load some signalling data for downstream analysis
savelocCCI <- "~/Dropbox/Cancer_pheno_evo/data/FELINE2/CommunicationOutputAllArms/"

# Are we studying the per indiv or per cell type level of signalling?
perIndiv=FALSE

# Names of LR data to analyze
filenamesCCI <- list.files(savelocCCI)#%>%length

# First let's do this for the population level
load(file ="~/Dropbox/Cancer_pheno_evo/data/FELINE2/Communication output merged ALL ARMS/PopulationCommunicationMerged ALL ARMS.RData" )#CCI,allphenotypes, uu,perIndiv,
CCI[,Treat:="CombinationRibo"]
CCI[ARM=="A",Treat:="LetrozoleAlone"]
CCI[,scaleTransduction:=scale(TransductionMu,center=F),by=c("Pair.Name")] 
CCI[!is.finite(TransductionMu),scaleTransduction:=0]
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


# For each communication pathway, day and cell type signal senders and receiver, randomize the tumor response annotations to remove the response related network structure
# Repeatmany times to determine the distribtution of communication differences that may be expected by chance under the null model that communication strength do not differ between resistant and sensitive tumors.
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

alldata <- merge(summaryCCIdiff,res1 %>%dplyr::select(-c(Nonresponse,Response,diffOfResponses)),by=c("Treat","LigandPhenoCelltype", "ReceptorPhenoCelltype", "Day"))
alldata$LigandPhenoCelltype <-factor(alldata$LigandPhenoCelltype, 
                                     levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Normal epithelial cells","Cancer cells","B cells","Macrophages","CD8+ T cells","CD4+ T cells"))
alldata$ReceptorPhenoCelltype <-factor(alldata$ReceptorPhenoCelltype, 
                                       levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Normal epithelial cells","Cancer cells","B cells","Macrophages","CD8+ T cells","CD4+ T cells") )


#save(res1,alldata,MultisummaryRandCCIdiff,nRandomizations,summaryCCI,file="~/Dropbox/Cancer_pheno_evo/data/FELINE2/Bootstrap communication results/Ribo and Letrozole Bootstrap communication results.RData")
load(file="~/Dropbox/Cancer_pheno_evo/data/FELINE2/Bootstrap communication results/Ribo and Letrozole Bootstrap communication results.RData")



pltdd <- alldata[adjusted_pval < 0.05][LigandPhenoCelltype!= "B cells"][ReceptorPhenoCelltype!= "B cells"]
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



ggplot(plotdd, aes(y= ReceptorPhenoCelltype, x= LigandPhenoCelltype, fill= obs.z2)) + 
  theme_classic(base_size= 19) +
  geom_tile() + 
  facet_wrap(Treatmentlab ~ Daylab) + 
  theme(axis.text.x= element_text(angle= 90,vjust=0.5)) + 
  scale_fill_gradientn(colours=myPalette(100), 
                       name="Communication \n in resistant \n relative to \n sensitive \n tumors") + 
  theme(aspect.ratio = 1) +
  labs(y= "Signal receiver cell type", x= "Signal sender cell type")

paperfile<- "~/Dropbox/FELINE Project (1)/Manuscript  Feline immune communication/Figures Communication Project/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Ribo and Letrozole Both D0and 180  Significant differences in cell type communication.png"),height=10,width=10)


ggplot(plotdd, aes(y= ReceptorPhenoCelltype, x= LigandPhenoCelltype, fill= obs.z2)) + 
  theme_classic(base_size= 19) +
  geom_tile() + 
  facet_wrap(Treatmentlab ~ Daylab) + 
  theme(axis.text.x= element_text(angle= 90,vjust=0.5)) + 
  scale_fill_gradientn(colours=myPalette(100), 
                       #limits= c(-1,1)*max(abs(plotdd[]$obs.z2)), 
                       name="Communication \n in resistant \n relative to \n sensitive \n tumors") + 
  theme(aspect.ratio = 1) +
  labs(y= "Signal receiver cell type", x= "Signal sender cell type")+
  theme(axis.title=element_blank(),  axis.text=element_blank(),axis.text.x=element_blank(),strip.text = element_blank(),legend.title=element_blank(),  legend.text=element_blank())
paperfile<- "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"BLANK Ribo and Letrozole Both D0and 180  Significant differences in cell type communication.png"),height=10,width=10)


ggplot(pltdd[Treat=="CombinationRibo"],       aes(y= ReceptorPhenoCelltype, x= LigandPhenoCelltype, fill= obs.z)) + 
  theme_classic(base_size= 16) +
  geom_tile() + 
  facet_wrap(Treat ~ Day) + 
  theme(axis.text.x= element_text(angle= 90,vjust=0.5)) + 
  scale_fill_gradientn(colours=myPalette(100), limits= c(-1,1)*max(abs(pltdd[Treat=="CombinationRibo"]$obs.z)), name="Communication strength in \n resistant relative to sensitive \n tumors") + 
  theme(aspect.ratio = 1) +
  labs(y= "Signal receiver population", x= "Signal sender population")

paperfile<- "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Ribo and Letrozole RIBO Significant differences in cell type communication.png"),height=10,width=10)

ggplot(pltdd[Treat=="CombinationRibo"],       aes(y= ReceptorPhenoCelltype, x= LigandPhenoCelltype, fill= obs.z)) + 
  theme_classic(base_size= 16) +
  geom_tile() + 
  facet_wrap(Treat ~ Day) + 
  theme(axis.text.x= element_text(angle= 90,vjust=0.5)) + 
  scale_fill_gradientn(colours=myPalette(100), limits= c(-1,1)*max(abs(pltdd[Treat=="CombinationRibo"]$obs.z)), name="Communication strength in \n resistant relative to sensitive \n tumors") + 
  theme(aspect.ratio = 1) +
  labs(y= "Signal receiver population", x= "Signal sender population")+
  theme(axis.title=element_blank(), axis.text.x=element_blank(), axis.text=element_blank(),strip.text = element_blank(),legend.title=element_blank(),legend.text=element_blank())
#ggsave(paste0(paperfile,"BLANK Ribo and Letrozole RIBO Significant differences in cell type communication.png"),height=10,width=10)



ggplot(pltdd[Treat!="CombinationRibo"][adjusted_pval<0.01][abs(obs.z)>8],
       aes(y= ReceptorPhenoCelltype, x= LigandPhenoCelltype, fill= obs.z)) + 
  theme_classic(base_size= 16) +
  geom_tile() + 
  facet_wrap(Treat ~ Day) + 
  theme(axis.text.x= element_text(angle= 90)) + 
  scale_fill_gradientn(colours=myPalette(100), limits= c(-1,1)*max(abs(pltdd[Treat!="CombinationRibo"]$obs.z)), name="Communication strength in \n resistant relative to sensitive \n tumors") + 
  theme(aspect.ratio = 1) +
  labs(y= "Signal receiver population", x= "Signal sender population")
ggsave(paste0(paperfile,"Ribo and Letrozole LETROZOLE Significant differences in cell type communication.png"),height=10,width=10)

ggplot(pltdd[Treat!="CombinationRibo"][adjusted_pval<0.01][abs(obs.z)>8],
       aes(y= ReceptorPhenoCelltype, x= LigandPhenoCelltype, fill= obs.z)) + 
  theme_classic(base_size= 16) +
  geom_tile() + 
  facet_wrap(Treat ~ Day) + 
  theme(axis.text.x= element_text(angle= 90)) + 
  scale_fill_gradientn(colours=myPalette(100), limits= c(-1,1)*max(abs(pltdd[Treat!="CombinationRibo"]$obs.z)), name="Communication strength in \n resistant relative to sensitive \n tumors") + 
  theme(aspect.ratio = 1) +
  labs(y= "Signal receiver population", x= "Signal sender population")+
  theme(axis.title=element_blank(), axis.text.x=element_blank(), axis.text=element_blank(),strip.text = element_blank(),legend.title=element_blank(),legend.text=element_blank())
ggsave(paste0(paperfile,"BLANK Ribo and Letrozole LETROZOLE Significant differences in cell type communication.png"),height=10,width=10)


ggplot(pltdd[Treat=="CombinationRibo"][adjusted_pval<0.001],
       aes(y= ReceptorPhenoCelltype, x= LigandPhenoCelltype, fill= diffOfResponses)) + 
  theme_classic(base_size= 16) +
  geom_tile() + 
  facet_wrap( ~ Day) + 
  theme(axis.text.x= element_text(angle= 90)) + 
  scale_fill_viridis_c(name= "Communication difference in \n resistant versus sensitive \n tumors",option= "A") + 
  theme(aspect.ratio = 1) +
  labs(y= "Signal receiver population", x= "Signal sender population")
#ggsave(filename = "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Strength communication cell types AllArms/signaling strength significant bootstrapped between cell types under Ribo combination therapy.png",width=12,height=6)



### globe plots
ngroup <- 1
set.seed(123);  lumrg <- data.table(Pair.Name= unique(CCI$Pair.Name) , grpvar= sample(1:ngroup, length(unique(CCI$Pair.Name) ) ,replace= T ) )
CCI[,scaleTransduction2:=exp(scale(log(TransductionMu),center=F)),by=c("Pair.Name")] 

average_muln_scaleTransduction <- data.table(CCI%>%
                                               group_by(Patient.Study.ID,LigandPhenoCelltype,ReceptorPhenoCelltype,Day,dynamic_class3,ARM,Treat) %>%
                                               dplyr::summarise(muln_scaleTransduction=mean(log(1+scaleTransduction))) %>% 
                                               group_by(LigandPhenoCelltype,ReceptorPhenoCelltype,Day,dynamic_class3,Treat) %>%
                                               dplyr::summarise(average_muln_scaleTransduction=mean(muln_scaleTransduction)))
average_muln_scaleTransduction$LigandPhenoCelltype <-factor(average_muln_scaleTransduction$LigandPhenoCelltype, 
                  levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Normal epithelial cells","Cancer cells","B cells","Macrophages","CD8+ T cells","CD4+ T cells"))
average_muln_scaleTransduction$ReceptorPhenoCelltype <-factor(average_muln_scaleTransduction$ReceptorPhenoCelltype, 
                   levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Normal epithelial cells","Cancer cells","B cells","Macrophages","CD8+ T cells","CD4+ T cells") )

initState <- data.table(average_muln_scaleTransduction[Day==0]%>%
                          group_by(LigandPhenoCelltype,ReceptorPhenoCelltype)%>% #Treat, dynamic_class3
                          dplyr::summarise(intiState=mean(average_muln_scaleTransduction)) )

average_muln_scaleTransduction2 <- merge(average_muln_scaleTransduction,initState,by=c("LigandPhenoCelltype", "ReceptorPhenoCelltype")) #,"Treat","dynamic_class3"
average_muln_scaleTransduction2[ , lnfoldchange:= (average_muln_scaleTransduction-intiState) ] 

CCIcondensed <- average_muln_scaleTransduction2[LigandPhenoCelltype!="B cells"][ReceptorPhenoCelltype!="B cells"]  

CCIcondensed$LigandPhenoCelltype <-factor(CCIcondensed$LigandPhenoCelltype, 
                       levels= sort(c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Normal epithelial cells","Cancer cells","Macrophages","CD8+ T cells","CD4+ T cells"))) #"B cells",
CCIcondensed$ReceptorPhenoCelltype <-factor(CCIcondensed$ReceptorPhenoCelltype, 
                       levels= sort(c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Normal epithelial cells","Cancer cells","Macrophages","CD8+ T cells","CD4+ T cells") )) #"B cells",
setnames(CCIcondensed, old=c("LigandPhenoCelltype", "ReceptorPhenoCelltype", "lnfoldchange"), new= c("from", "to", "weight"))

# specifcy color of nodes
node_col <- ggsci::pal_npg("nrc")(2)  
# circle layout
n <- length(unique(CCIcondensed$to)) -1
pts.circle <- t(sapply(1:n, function(r)c(cos(2*r*pi/n),sin(2*r*pi/n))))
NodeList <- data.table(c("Cancer cells", "CD4+ T cells" , "Adipocytes", "Fibroblasts", "Normal epithelial cells", "Pericytes"  , "Endothelial cells" ,"Macrophages", "CD8+ T cells" )  ,
 c(0, pts.circle[,1] )  ,  c(0, pts.circle[,2] ) ) #"B cells",
presloc <- data.frame(NodeList)[order(NodeList$V1),]

tmp <- CCIcondensed[Treat== "CombinationRibo"][Day== 180][order( Treat,from,to)]#[Day== 0]
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


tmp <- CCIcondensed[][Day== 180][order( Treat,from,to)]
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
  geom_node_point( col= rep( rep(node_col, each= nrow(createdLayout) ) ,  2) , size= 20 ) + 
  geom_node_text(aes(label= name),size=8 ) + #,col= "grey80" ) + #, size= 9) +
  facet_edges(Treatmentlab~TumorResponse  ) + 
  theme_void(base_size= 1.5* 14 ) +
  theme(aspect.ratio = 1, legend.position = "none" ) + 
  scale_edge_width(range=c(0.0,2))+
  scale_edge_alpha(range=c(0.0,0.97)) +
  scale_edge_color_gradient2()
p0 + scale_edge_color_gradient2(low="white",high="black")
ggsave(   file="~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Globe plot ggrraph/Globe plot ggrraph1.png",height=10, width = 10)


p0<-ggraph(createdLayout , aes(x= x, y= y , col= dynamic_class3 ) ) + 
  geom_edge_loop(  aes(edge_colour=weight,edge_width= (weight),edge_alpha= (weight)) ,  arrow= arrow(length= unit(3, "mm" ) , type= "open"),
                   start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" ) ) +
  geom_edge_arc( aes(edge_colour=weight,edge_width= (weight),edge_alpha= (weight)  ) , arrow= arrow(length= unit(3, "mm" ) , type= "open" ),
                 start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" )  ,strength=0.2) +
  geom_node_point( col= rep( rep(node_col, each= nrow(createdLayout) ) ,  2) , size= 20 ) + 
  geom_node_text(aes(label= name),size=6 ) + #,col= "grey80" ) + #, size= 9) +
  facet_edges(Treatmentlab~TumorResponse  ) + 
  theme_void(base_size= 1.5* 14 ) +
  theme(aspect.ratio = 1, legend.position = "none" ) + 
  scale_edge_width(range=c(0.0,2))+
  scale_edge_alpha(range=c(0.0,0.97)) +
  scale_edge_color_gradient2()
p0 + scale_edge_color_gradient2(low="white", high="black")

p0<-ggraph(createdLayout , aes(x= x, y= y , col= dynamic_class3 ) ) + 
  geom_edge_loop(  aes(edge_colour=weight,edge_width= (weight),edge_alpha= (weight)) ,  arrow= arrow(length= unit(3, "mm" ) , type= "open"),
                   start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" ) ) +
  geom_edge_arc( aes(edge_colour=weight,edge_width= (weight),edge_alpha= (weight)  ) , arrow= arrow(length= unit(3, "mm" ) , type= "open" ),
                 start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" ) ,strength=0.25 ) +
  geom_node_point( col= rep( rep(node_col, each= nrow(createdLayout) ) ,  2) , size= 20 ) + 
  geom_node_text(aes(label= name),size=6, col="black" ) + #,col= "grey80" ) + #, size= 9) +
  facet_edges(Treatmentlab~TumorResponse  ) + 
  theme_void(base_size= 1.5* 18) +
  theme(aspect.ratio = 1, legend.position = "none" ) + 
  scale_edge_width(range=c(0.0,2))+
  scale_edge_alpha(range=c(0.0,0.97)) +
  scale_edge_color_gradient2()
p0 + scale_edge_color_gradient2(low="white",high="black")+coord_cartesian(clip="off")+theme(aspect.ratio=1,
                                                                                            plot.margin=unit(c(10,30,30,30), "points"),
                                                                                            panel.spacing=unit(2, "lines") ,
                                                                                            strip.text=element_text(size=22))

ggsave(   file="~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Globe plot ggrraph/Globe plot ggrraph 2.png",height=12, width = 10)

paperfile<- "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
paperfile<- "~/Dropbox/FELINE Project (1)/Manuscript  Feline immune communication/Figures Communication Project/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Ribo and Letrozole Globe plot ggrraph.png"),height=10,width=10)

p0<-ggraph(createdLayout , aes(x= x, y= y , col= dynamic_class3 ) ) + 
  geom_edge_loop(  aes(edge_colour=weight,edge_width= (weight),edge_alpha= (weight)) ,  arrow= arrow(length= unit(3, "mm" ) , type= "open"),
                   start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" ) ) +
  geom_edge_arc( aes(edge_colour=weight,edge_width= (weight),edge_alpha= (weight)  ) , arrow= arrow(length= unit(3, "mm" ) , type= "open" ),
                 start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" ) ,strength=0.25 ) +
  geom_node_point( col= rep( rep(node_col, each= nrow(createdLayout) ) ,  2) , size= 20 ) + 
  facet_edges(Treatmentlab~TumorResponse  ) + 
  theme_void(base_size= 1.5* 14 ) +
  theme(aspect.ratio = 1, legend.position = "none" ) + 
  scale_edge_width(range=c(0.0,2))+
  scale_edge_alpha(range=c(0.0,0.97)) +
  scale_edge_color_gradient2() + scale_edge_color_gradient2(low="white",high="black")+coord_cartesian(clip="off")+theme(aspect.ratio=1,
                                                                                            plot.margin=unit(c(10,30,30,30), "points"),
                                                                                            panel.spacing=unit(2, "lines") )

p0 + scale_edge_color_gradient2(low="white",high="black") +
  theme(axis.title=element_blank(),  axis.text=element_blank(),strip.text = element_blank())
ggsave(paste0(paperfile,"BLANK Ribo and Letrozole Globe plot ggrraph.png"),height=10,width=10)


tmp <- CCIcondensed[][Day== 0][order( Treat,from,to)]
tmp$DayLab <- paste0("Day ", tmp$Day )
tmp <- data.table( tmp %>% group_by(Treat) %>% mutate(weight = scale(weight) ) )
tmp[,Treatmentlab:= "Letrozole alone"]
tmp[Treat=="CombinationRibo",Treatmentlab:= "Combination ribociclib"]
tmp[,TumorResponse:="Resistant"]
tmp[dynamic_class3=="Response",TumorResponse:="Sensitive"]
tmp[from=="Normal epithelial cells",from:="Normal \n epithelial \n cells"]
tmp[to=="Normal epithelial cells",to:="Normal \n epithelial \n cells"]
tmp[from=="Endothelial cells",from:="Endothelial \n cells"]
tmp[to=="Endothelial cells",to:="Endothelial \n cells"]
tmp[from=="CD8+ T cells",from:="CD8+ \n T cells"]
tmp[to=="CD8+ T cells",to:="CD8+ \n T cells"]
tmp[from=="CD4+ T cells",from:="CD4+ \n T cells"]
tmp[to=="CD4+ T cells",to:="CD4+ \n T cells"]
tmp[from=="Cancer cells",from:="Cancer \n cells"]
tmp[to=="Cancer cells",to:="Cancer \n cells"]
graphdd <- as_tbl_graph( tmp, directed= T)
createdLayout <- create_layout(graphdd, layout= "star")
createdLayout$x <- presloc$V2
createdLayout$y <- presloc$V3

p0<-ggraph(createdLayout , aes(x= x, y= y , col= dynamic_class3 ) ) + 
  geom_edge_loop(  aes(edge_colour=weight,edge_width= (weight),edge_alpha= (weight)) ,  arrow= arrow(length= unit(3, "mm" ) , type= "open"),
                   start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" ) ) +
  geom_edge_arc( aes(edge_colour=weight,edge_width= (weight),edge_alpha= (weight)  ) , arrow= arrow(length= unit(3, "mm" ) , type= "open" ),
                 start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" )  ) +
  geom_node_point( col= rep( rep(node_col, each= nrow(createdLayout) ) ,  2) , size= 20 ) + 
  geom_node_text(aes(label= name) ) + #,col= "grey80" ) + #, size= 9) +
  facet_edges(Treatmentlab~TumorResponse  ) + 
  theme_void(base_size= 14 ) +
  theme(aspect.ratio = 1, legend.position = "none" ) + 
  scale_edge_width(range=c(0.0,2))+
  scale_edge_alpha(range=c(0.0,0.97)) +
  scale_edge_color_gradient2()
p0 + scale_edge_color_gradient2(low="white",high="black")
ggsave(   file="~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Globe plot ggrraph/Globe plot Day 0 ggrraph.png",height=10, width = 10)


















p1<-ggraph(createdLayout , aes(x= x, y= y , col= dynamic_class3 ) ) + 
  geom_edge_loop(  aes(edge_colour=weight,edge_width= exp(weight),edge_alpha= exp(weight)) ,  arrow= arrow(length= unit(5, "mm" ) , type= "closed"),
                   start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" ) ) +
  geom_edge_arc( aes(edge_colour=weight,edge_width= exp(weight),edge_alpha= exp(weight)  ) , arrow= arrow(length= unit(5, "mm" ) , type= "closed" ),
                 start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" )  ) +
  #geom_node_point( col= rep( rep(node_col, each= nrow(createdLayout) ) , each= 2) , size= 20 ) + 
  geom_node_point( col= rep( rep(node_col, each= nrow(createdLayout) ) ,  2) , size= 20 ) + 
  geom_node_text(aes(label= name) ) + #,col= "grey80" ) + #, size= 9) +
  facet_edges(Treat~dynamic_class3  ) + 
  theme_void(base_size= 12 ) +
  theme(aspect.ratio = 1, legend.position = "none" ) + 
  scale_edge_width(range=c(0.0,3.5))+
  scale_edge_alpha(range=c(0.4,0.97)) +
  scale_edge_color_gradient2()
 # scale_edge_color_viridis(option="B")

p1+ scale_edge_color_gradient()
p1+ scale_edge_color_gradient2(low="white",high="black")




tmp <- CCIcondensed[to=="Macrophages"][Day== 0][order( Treat,from,to)]
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





ggraph(createdLayout , aes(x= x, y= y , col= dynamic_class3 ) ) + 
  geom_edge_loop(edge_col="grey50", aes(edge_width= 0.01*sqrt(exp(weight)),edge_alpha= sqrt(exp(weight))) ,  arrow= arrow(length= unit(4, "mm" ) , type= "open"),
                 start_cap= circle(5, "mm" ),  end_cap= circle(3, "mm" ) ) +
  geom_edge_arc(edge_col="grey50", aes(edge_width= 0.01*sqrt(exp(weight)),edge_alpha= sqrt(exp(weight))  ) , arrow= arrow(length= unit(4, "mm" ) , type= "open" ),
                start_cap= circle(5, "mm" ),  end_cap= circle(3, "mm" )  ) +
  geom_node_point( col= rep( rep(node_col, each= nrow(createdLayout) ) , each= 3) , size= 5 ) + 
  geom_node_text(aes(label= name) ) + #,col= "grey80" ) + #, size= 9) +
  facet_edges(dynamic_class3 ~ DayLab ) + 
  theme_void(base_size= 12 ) +
  theme(aspect.ratio = 1, legend.position = "none" ) 


ggraph(createdLayout , aes(x= x, y= y , col= dynamic_class3 ) ) + 
  geom_edge_loop( aes(edge_colour=weight, edge_width= weight,edge_alpha= weight) ,  arrow= arrow(length= unit(4, "mm" ) , type= "open"),
                 start_cap= circle(5, "mm" ),  end_cap= circle(3, "mm" ) ) +
  geom_edge_arc( aes(edge_colour=weight, edge_width= weight,edge_alpha= weight  ) , arrow= arrow(length= unit(4, "mm" ) , type= "open" ),
                start_cap= circle(5, "mm" ),  end_cap= circle(3, "mm" )  ) +
  geom_node_point( col= rep( rep(node_col, each= nrow(createdLayout) ) , each= 3) , size= 5 ) + 
  geom_node_text(aes(label= name) ) + #,col= "grey80" ) + #, size= 9) +
  facet_edges(dynamic_class3 ~ DayLab ) + 
  theme_void(base_size= 12 ) +
  theme(aspect.ratio = 1, legend.position = "none" ) +scale_edge_color_viridis()+
  scale_edge_width(range=c(0.1,2))+
  scale_edge_alpha(range=c(0.1,0.8))



ggraph(createdLayout , aes(x= x, y= y , col= dynamic_class3 ) ) + 
 # geom_edge_loop(edge_col="grey50", aes(edge_width= ((weight)),edge_alpha= sqrt(exp(weight))) ,  arrow= arrow(length= unit(1, "mm" ) , type= "open"),
  #               start_cap= circle(10, "mm" ),  end_cap= circle(10, "mm" ) ) +
  geom_edge_arc(edge_col="grey50", aes(edge_width= (sqrt(weight)),edge_alpha= sqrt(exp(weight))  ) , arrow= arrow(length= unit(4, "mm" ) , type= "open" ),
                start_cap= circle(5, "mm" ),  end_cap= circle(7, "mm" )  ) +
  #geom_node_point( col= rep( rep(node_col, each= nrow(createdLayout) ) , each= 3) , size= 20 ) + 
  #geom_node_text(aes(label= name) ) + #,col= "grey80" ) + #, size= 9) +
  facet_edges(dynamic_class3 ~ DayLab ) + 
  theme_void(base_size= 12 ) +
  theme(aspect.ratio = 1, legend.position = "none" ) +
  scale_edge_width(range=c(0.1,4))+
  scale_edge_alpha(range=c(0.1,0.8))








ggraph(createdLayout ,aes(x=x,y=y )) + geom_edge_arc2(aes(edge_width=weight)) +geom_node_point() + geom_node_text(aes(label=name),col="blue")+
  facet_edges(dynamic_class3~Day)


createdLayout <- layout_tbl_graph_manual(graphdd,x=presloc$V2,y=presloc$V3)


ggraph::layout_tbl_graph_manual(graphdd,x=presloc$V2,y=presloc$V3)
ggraphdata<- CCIcondensed[RiboTreated== T]
ggraphdata$from <- factor(ggraphdata$from, levels= c("Cancer cells", "CD4+ T cells" ,"Adipocytes","Fibroblasts","Normal epithelial cells","Pericytes"  , "Endothelial cells" ,"Macrophages", "B cells","CD8+ T cells" )  )
ggraphdata$to <- factor(ggraphdata$to, levels= c("Cancer cells", "CD4+ T cells" ,"Adipocytes","Fibroblasts","Normal epithelial cells","Pericytes"  , "Endothelial cells" ,"Macrophages", "B cells","CD8+ T cells" )  )
graphdd <- as_tbl_graph( ggraphdata[order(from,to)], directed= T)
#graphdd <- graph_from_data_frame( CCIcondensed[RiboTreated==T],directed=T)















plt1 <- summaryCCI[LigandPhenoCelltype!="B cells"&ReceptorPhenoCelltype!="B cells"]
plt1$LigandPhenoCelltype <-factor(plt1$LigandPhenoCelltype, 
                                  levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Normal epithelial cells","Cancer cells","B cells","Macrophages","CD8+ T cells","CD4+ T cells"))
plt1$ReceptorPhenoCelltype <-factor(plt1$ReceptorPhenoCelltype, 
                                    levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Normal epithelial cells","Cancer cells","B cells","Macrophages","CD8+ T cells","CD4+ T cells") )


ggplot(plt1, aes(y= meanSig, x= ReceptorPhenoCelltype, col= ReceptorPhenoCelltype,fill=dynamic_class3,group=interaction(dynamic_class3,ReceptorPhenoCelltype)))+geom_violin(col=NA,scale="width",alpha=0.5)+ geom_point(position = position_dodge(width=1))+facet_wrap(~Day)+theme(axis.text.x=element_text(angle=90))
ggplot(plt1, aes(y= meanSig, x= LigandPhenoCelltype, col= LigandPhenoCelltype,fill=dynamic_class3,group=interaction(dynamic_class3,LigandPhenoCelltype)))+geom_violin(scale="width",alpha=0.5)+ geom_point()+facet_wrap(~Day)+theme(axis.text.x=element_text(angle=90))

ggplot(plt1, aes(y= log(exp(meanSig)-1), x= LigandPhenoCelltype, col= LigandPhenoCelltype,fill=LigandPhenoCelltype,group=interaction(LigandPhenoCelltype)))+geom_violin(scale="width",alpha=0.5)+ geom_point()+facet_wrap(~Day)+theme(axis.text.x=element_text(angle=90))
ggplot(plt1, aes(y= log(exp(meanSig)-1), x= ReceptorPhenoCelltype, col= ReceptorPhenoCelltype,fill=ReceptorPhenoCelltype,group=interaction(ReceptorPhenoCelltype)))+geom_violin(col=NA,scale="width",alpha=0.5)+ geom_point()+facet_wrap(~Day)+theme(axis.text.x=element_text(angle=90))

ggplot(plt1, aes(y= meanSig, x= LigandPhenoCelltype, col= LigandPhenoCelltype,fill=LigandPhenoCelltype,group=interaction(LigandPhenoCelltype)))+geom_violin(scale="width",alpha=0.5)+ geom_point()+theme_classic(base_size=16)+theme(axis.text.x=element_text(angle=90))+
  theme(aspect.ratio=1,legend.position="none") +labs(y="Contribution to tumor wide communication (mean)",x= "Cell type") 
ggplot(plt1, aes(y= log(exp(meanSig)-1), x= ReceptorPhenoCelltype, col= ReceptorPhenoCelltype,fill=ReceptorPhenoCelltype,group=interaction(ReceptorPhenoCelltype)))+geom_violin(col=NA,scale="width",alpha=0.5)+ geom_point()+theme_classic(base_size=16)+theme(axis.text.x=element_text(angle=90)) +
  theme(aspect.ratio=1,legend.position="none") +labs(y="Strength of communicaiton received (mean)",x= "Cell type") 

exp( alldata[Day==0]$diffOfResponses )
###

par(mfrow = c(1, 3) )
for(i in c(0, 14, 180) ){
  Y <- alldata[Day==i][LigandPhenoCelltype!="B cells" & ReceptorPhenoCelltype!= "B cells"]
  g <- graph.data.frame(Y %>% dplyr::select(LigandPhenoCelltype, ReceptorPhenoCelltype), directed=TRUE)
  E(g)$weight <- abs( Y$diffOfResponses ) #exp((log(CCI_plot$SumSignal)-centVal)/scalsd)    #is.weighted((g))
  # edge color
  coledge<- rep( "red" , length(E(g)$weight) )
  coledge[ Y$diffOfResponses < 0 ] <- "blue"
  coledge[ Y$adjusted_pval >= 0.05 ] <-"grey"
  E(g)$color <- coledge
  
  # specifcy color of nodes
  resp_cols <- "green" #ggsci::pal_npg("nrc")(2)
  V(g)$color <- resp_cols[1]
  
  # circle layout.
  n = length(unique(Y$ReceptorPhenoCelltype)) -1
  pts.circle <- t( sapply(1:n, function(r)c(cos(2*r*pi/n), sin(2*r*pi/n))) )
  
  #NodeList <- data.table(sort( unique(SignalReceiverRes$ReceptorPhenoCelltype) ), c(4,9.5,5,8,9,4,2,9,1,1.5 ) /5-1 ,c(9,5,5,8,6,1,7,2,4,1.5)  /5-1 )
  NodeList <- data.table(c("Cancer cells", "CD4+ T cells" ,"Adipocytes","Fibroblasts","Normal epithelial cells","Pericytes"  , "Endothelial cells" ,"Macrophages", "CD8+ T cells" )  , c(0, pts.circle[,1] )  ,  c(0, pts.circle[,2] ) )
  presloc <- NodeList[na.omit((match(names(V(g) ), NodeList$V1  ))), ]
  
  plot.igraph(g, rescale = FALSE,#vertex.label=V(g)$name,
              layout = as.matrix(presloc %>% dplyr::select(-V1))  , 
              xlim = c(-1, 1), ylim = c(-1,1),#layout=layout.fruchterman.reingold, 
              vertices= NodeList[V1 %in% presloc$V1],
              edge.label.color = E(g)$color , edge.arrow.size=0.02, edge.arrow.width=0.2,
              edge.weight= 1000*( 1 + E(g)$weight ) )
  
}


#######
#
#all six globes
par(mfrow=c(1,1))
luglode <- data.table( expand.grid(Day=c(0,14,180), dynamic_class3=c("Non-response","Response")) )
centVal <- mean(log( CCI$SumSignal ))
scalsd<-sd(log(CCI$SumSignal))
RiboTreat <- TRUE


for(i in 1:nrow(luglode)){  #median(log(1+scaleTransduction))
  CCI_plot <- CCI[CCI$SumSignal > centVal][dynamic_class3== luglode[i]$dynamic_class3][Day== luglode[i]$Day][RiboTreated== RiboTreat][order(LigandPhenoCelltype, ReceptorPhenoCelltype)]
  g <- graph.data.frame(CCI_plot %>% dplyr::select(-Patient.Study.ID), directed=TRUE)
  E(g)$weight <- exp( (log(CCI_plot$SumSignal)-centVal)/scalsd ) #E(g)$weight <- CCI_plot$scaleTransduction    #is.weighted((g))
  # specifcy color of nodes
  resp_cols <- ggsci::pal_npg("nrc")(2)
  if(luglode[i]$dynamic_class3 == "Non-response"){V(g)$color <- resp_cols[1] }else{if(luglode[i]$dynamic_class3 == "Response"){V(g)$color <- resp_cols[2]}}
  # circle layout.
  n= length(unique(CCI_plot$ReceptorPhenoCelltype)) -1
  pts.circle <- t(sapply(1:n, function(r)c(cos(2*r*pi/n),sin(2*r*pi/n))))
  #NodeList <- data.table(sort( unique(SignalReceiverRes$ReceptorPhenoCelltype) ), c(4,9.5,5,8,9,4,2,9,1,1.5 ) /5-1 ,c(9,5,5,8,6,1,7,2,4,1.5)  /5-1 )
  NodeList <- data.table(c("Cancer cells", "CD4+ T cells" ,"Adipocytes","Fibroblasts","Normal epithelial cells","Pericytes"  , "Endothelial cells" ,"Macrophages", "B cells","CD8+ T cells" )  , c(0, pts.circle[,1] )  ,  c(0, pts.circle[,2] ) )
  presloc <- NodeList[na.omit((match(names(V(g) ), NodeList$V1  ))), ]
  
  # light version
  pdf(file= paste0("/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Globes_Communication globes3/GlobeCommunication",
                   luglode[i]$dynamic_class3,"_", luglode[i]$Day, ".pdf" ), width= 20, height= 20 )
  par(mfrow= c(1, 1))
  # option1 
  #plot.igraph(g,rescale=FALSE,vertex.label=NULL,#vertex.label=V(g)$name,
  #         layout = as.matrix(presloc%>%select(-V1))   , xlim =c(-1,1), ylim =c(-1,1),#layout=layout.fruchterman.reingold, 
  #         vertices=NodeList[V1%in%presloc$V1],
  #         edge.label.color = adjustcolor("black", .9),
  #         edge.color= adjustcolor("black", 1),edge.width= 1e-16*E(g)$weight)
  #option2
  #plot.igraph(g, rescale= FALSE,
  #           vertex.label= V(g)$name,
  #           layout = as.matrix(presloc %>% dplyr::select( -V1 ))   , xlim = c(-1, 1), ylim = c(-1, 1),#layout=layout.fruchterman.reingold, 
  #           vertices= NodeList[V1 %in% presloc$V1],
  #           edge.label.color = adjustcolor("black", .75),
  #           edge.color= adjustcolor("black", .005), edge.width= 0.000000001*E(g)$weight)
  #option3
  plot.igraph(g,rescale=FALSE,vertex.label="",#vertex.label=V(g)$name,
              layout = as.matrix(presloc%>%select(-V1))   , xlim =c(-1,1), ylim =c(-1,1),#layout=layout.fruchterman.reingold,
              vertices= NodeList[V1%in%presloc$V1],
              edge.label.color = adjustcolor("black", .9),
              edge.color= adjustcolor("black", 0.08), edge.width= 1e-36*E(g)$weight)
  dev.off()
}


#######









ggplot(alldata[adjusted_pval<0.05][LigandPhenoCelltype!="B cells"][ReceptorPhenoCelltype!="B cells"],
       aes(y=ReceptorPhenoCelltype,x=LigandPhenoCelltype,fill=obs.z))+ theme_classic()+
  geom_tile()+facet_wrap(~Day)+theme(axis.text.x=element_text(angle=90)) + scale_fill_viridis_c(option="D")


ggplot(alldata[adjusted_pval<0.05][],aes(y=ReceptorPhenoCelltype,x=LigandPhenoCelltype,fill=obs.z))+ theme_classic()+
  geom_tile()+facet_wrap(~Day)+theme(axis.text.x=element_text(angle=90)) + scale_fill_viridis_c(option="A")

ggplot(alldata[adjusted_pval<0.05][Day==0],aes(y=ReceptorPhenoCelltype,x=LigandPhenoCelltype,fill=obs.z))+
  geom_tile()


geom_violin(scale="width",alpha=0.5)+ geom_point()+facet_wrap(~Day)+theme(axis.text.x=element_text(angle=90))

ggplot(alldata,aes(y=diffOfResponses,x=LigandPhenoCelltype,fill=adjusted_pval<0.05,col=adjusted_pval<0.05))+geom_violin(scale="width",alpha=0.5)+ geom_point()+facet_wrap(~Day)+theme(axis.text.x=element_text(angle=90))
ggplot(alldata,aes(y=diffOfResponses,x=ReceptorPhenoCelltype,fill=adjusted_pval<0.05,col=adjusted_pval<0.05))+geom_violin(scale="width",alpha=0.5)+ geom_point()+facet_wrap(~Day)+theme(axis.text.x=element_text(angle=90))

ggplot(summaryCCIdiff,aes(y=diffOfResponses,x=ReceptorPhenoCelltype))+geom_violin(alpha=0.5)+ geom_point()+facet_wrap(~Day)+theme(axis.text.x=element_text(angle=90))
ggplot(summaryCCIdiff,aes(y=diffOfResponses,x=LigandPhenoCelltype))+geom_violin(alpha=0.5)+ geom_point()+facet_wrap(~Day)+theme(axis.text.x=element_text(angle=90))

summaryCCI[,isSenderCancer:=F]
summaryCCI[LigandPhenoCelltype=="Cancer cells",isSenderCancer:=T]
summaryCCI[,isReceiverCancer:=F]
summaryCCI[ReceptorPhenoCelltype=="Cancer cells",isReceiverCancer:=T]
ggplot(summaryCCI[LigandPhenoCelltype!="B cells"&ReceptorPhenoCelltype!="B cells"],aes(y=meanSig,x=ReceptorPhenoCelltype))+geom_violin(alpha=0.5)+ geom_point(aes(col=isSenderCancer))+facet_wrap(~Day)+theme(axis.text.x=element_text(angle=90))
ggplot(summaryCCI[LigandPhenoCelltype!="B cells"&ReceptorPhenoCelltype!="B cells"],aes(y=meanSig,x=LigandPhenoCelltype, col=dynamic_class3))+geom_violin(scale="width",alpha=0.5)+ geom_point()+facet_wrap(~Day)+theme(axis.text.x=element_text(angle=90))
ggplot(summaryCCI[LigandPhenoCelltype=="Cancer cells"],aes(y=meanSig,x=ReceptorPhenoCelltype))+geom_violin(alpha=0.5)+ geom_point(aes(col=isSenderCancer))+facet_wrap(~Day)+theme(axis.text.x=element_text(angle=90))
ggplot(CCI[LigandPhenoCelltype!="B cells"&ReceptorPhenoCelltype!="B cells"],aes(y=log(1+scaleTransduction),x=ReceptorPhenoCelltype))+geom_violin(alpha=0.5)+ geom_point()+facet_wrap(~Day)+theme(axis.text.x=element_text(angle=90))

# non-response - response
ggplot( res1 , aes(obs.z, pval) ) + geom_point()
ggplot( res1 , aes(diffOfResponses, rand.meanDif, col= abs(obs.z) ) ) + geom_point()
ggplot( res1 , aes(diffOfResponses, rand.meanDif, col= -sqrt( adjusted_pval ) ) ) + geom_point()
ggplot( res1 , aes(obs.z, diffOfResponses) ) + geom_point()


ggplot( res1 , aes(obs.z, adjusted_pval,col=LigandPhenoCelltype ) ) + geom_point()+facet_wrap(~Day)
ggplot( res1[adjusted_pval < 0.05] , aes(diffOfResponses, obs.z,col=LigandPhenoCelltype ) ) + geom_point()+facet_wrap(~Day)
ggplot( res1[adjusted_pval < 0.05] , aes(obs.z, ReceptorPhenoCelltype, col=LigandPhenoCelltype ) ) + geom_point()+facet_wrap(~Day)
ggplot( res1 , aes(diffOfResponses, LigandPhenoCelltype , col=ReceptorPhenoCelltype ) ) + geom_point()+facet_wrap(~Day)

res1[adjusted_pval < 0.05][Day==0] %>% dplyr::select(LigandPhenoCelltype, ReceptorPhenoCelltype) %>% table() %>% rowSums()
res1[adjusted_pval < 0.05][Day==0] %>% dplyr::select(LigandPhenoCelltype, ReceptorPhenoCelltype) %>% table() %>% colSums()

res1[adjusted_pval < 0.05][Day==0] %>% dplyr::select(LigandPhenoCelltype, ReceptorPhenoCelltype) %>% table() %>% rowSums()/ length(unique(res1$LigandPhenoCelltype))
res1[adjusted_pval < 0.05][Day==0] %>% dplyr::select(LigandPhenoCelltype, ReceptorPhenoCelltype) %>% table() %>% colSums()/ length(unique(res1$LigandPhenoCelltype))

res1[adjusted_pval < 0.05][Day==14] %>% dplyr::select(LigandPhenoCelltype, ReceptorPhenoCelltype) %>% table() %>% rowSums()/ length(unique(res1$LigandPhenoCelltype))
res1[adjusted_pval < 0.05][Day==14] %>% dplyr::select(LigandPhenoCelltype, ReceptorPhenoCelltype) %>% table() %>% colSums()/ length(unique(res1$LigandPhenoCelltype))

res1[adjusted_pval < 0.05][Day==180] %>% dplyr::select(LigandPhenoCelltype, ReceptorPhenoCelltype) %>% table() %>% rowSums()/ length(unique(res1$LigandPhenoCelltype))
res1[adjusted_pval < 0.05][Day==180] %>% dplyr::select(LigandPhenoCelltype, ReceptorPhenoCelltype) %>% table() %>% colSums()/ length(unique(res1$LigandPhenoCelltype))

res1[adjusted_pval < 0.05][Day==180][ReceptorPhenoCelltype=="CD8+ T cells"]
res1[ adjusted_pval < 0.05 ]

res1



RRR <- "Response"
RRR <- "Non-response"
Day_var <- 180
TTT <- "CombinationRibo"
RiboTreat <- TRUE
# 
CCI_plot <- CCI[dynamic_class3==RRR][Day==Day_var][ ][order(LigandPhenoCelltype,ReceptorPhenoCelltype)][Treat==TTT]

CCI[ Pair.Name== "ADAM10_AXL"][LigandPhenoCelltype== "Adipocytes"][ReceptorPhenoCelltype== "Adipocytes"]

randomizeCCI <- CCI[Treat==TTT]# [ Pair.Name== "ADAM10_AXL"][LigandPhenoCelltype== "Adipocytes"][ReceptorPhenoCelltype== "Adipocytes"] 
randomizeCCI[ ,Rand_dynamic_class3:= sample(dynamic_class3) ,by= c(Pair.Name, LigandPhenoCelltype, Day, ReceptorPhenoCelltype) ]


summaryCCI <- data.table( randomizeCCI[order(LigandPhenoCelltype, ReceptorPhenoCelltype)]%>% 
                            group_by(LigandPhenoCelltype, ReceptorPhenoCelltype,Rand_dynamic_class3,Day) %>% 
                            dplyr::summarise(meanSig= mean(scaleTransduction)) )
spread(summaryCCI, dynamic_class3, meanSig, fill= 0 ) %>% setnames(old= "Non-response", new= "Nonresponse") %>% mutate(diffOfResponses= Nonresponse - Response)



summaryCCI <- data.table( CCI[order(LigandPhenoCelltype, ReceptorPhenoCelltype)][Treat==TTT] %>% 
                            group_by(LigandPhenoCelltype, ReceptorPhenoCelltype,dynamic_class3,Day) %>% 
                            dplyr::summarise(meanSig= mean(scaleTransduction)) )
spread(summaryCCI, dynamic_class3, meanSig, fill= 0 ) %>% setnames(old= "Non-response", new= "Nonresponse") %>% mutate(diffOfResponses= Nonresponse - Response)


# Randomization



# Construct directed graph
g <- graph.data.frame(CCI_plot %>% dplyr::select(LigandPhenoCelltype, ReceptorPhenoCelltype), directed= TRUE) # -Patient.Study.ID), directed= TRUE)
# Add weights to edges of the graph
E(g)$weight <- CCI_plot$scaleTransduction
# Specifcy color of nodes
resp_cols<- ggsci::pal_npg("nrc")(2)
if(RRR=="Non-response"){V(g)$color <- resp_cols[1] }else{if(RRR=="Response"){V(g)$color <- resp_cols[2]}}

# Generate circle layout
n <- length( unique(CCI_plot$ReceptorPhenoCelltype) ) -1

pts.circle <- t(sapply(1:n,function(r)c(cos(2*r*pi/n),sin(2*r*pi/n))))
NodeList <- data.table(c("Cancer cells", "CD4+ T cells" ,"Adipocytes","Fibroblasts","Normal epithelial cells","Pericytes"  , "Endothelial cells" ,"Macrophages", "B cells","CD8+ T cells" )  , c(0, pts.circle[,1] )  ,  c(0, pts.circle[,2] ) )

presloc <- NodeList[na.omit((match(names(V(g) ), NodeList$V1  ))), ]
# Visualize communication with weighted graph
plot.igraph(g, rescale= FALSE,
            layout = as.matrix(presloc %>% select(-V1))   ,
            xlim = c(-1,1), ylim =c(-1,1),
            vertices= NodeList[V1 %in% presloc$V1], edge.label.color = adjustcolor("black", 0.5), edge.color= adjustcolor("black", 0.5), edge.width= 0.1*E(g)$weight )


#####







CCI_plot[,signalCol:="blue"]
CCI_plot[LigandPhenoCelltype=="Cancer cells",signalCol:="red"]
centVal<-mean(log(CCI$TransductionMu))
scalsd<-sd(log(CCI$TransductionMu))


RRR <- "Response"
RRR <- "Non-response"
Day_var <- 180
RiboTreat <- TRUE
#CCI_plot <- CCI[dynamic_class3==RRR][Day==Day_var][Pair.Name%in%summarytable$key_][order(LigandPhenoCelltype,ReceptorPhenoCelltype)]
CCI_plot <- CCI[CCI$SumSignal > centVal][dynamic_class3==RRR][Day==Day_var][RiboTreated==RiboTreat][order(LigandPhenoCelltype,ReceptorPhenoCelltype)]
g <- graph.data.frame(CCI_plot%>%dplyr::select(-Patient.Study.ID), directed=TRUE)
E(g)$weight <- exp((log(CCI_plot$SumSignal)-centVal)/scalsd)
is.weighted((g))

# specifcy color of nodes
resp_cols<- ggsci::pal_npg("nrc")(2)
if(RRR=="Non-response"){V(g)$color <- resp_cols[1] }else{if(RRR=="Response"){V(g)$color <- resp_cols[2]}}


# circle layout.
n= length(unique(CCI_plot$ReceptorPhenoCelltype)) -1
pts.circle <- t(sapply(1:n,function(r)c(cos(2*r*pi/n),sin(2*r*pi/n))))

#NodeList <- data.table(sort( unique(SignalReceiverRes$ReceptorPhenoCelltype) ), c(4,9.5,5,8,9,4,2,9,1,1.5 ) /5-1 ,c(9,5,5,8,6,1,7,2,4,1.5)  /5-1 )
NodeList <- data.table(c("Cancer cells", "CD4+ T cells" ,"Adipocytes","Fibroblasts","Normal epithelial cells","Pericytes"  , "Endothelial cells" ,"Macrophages", "B cells","CD8+ T cells" )  , c(0, pts.circle[,1] )  ,  c(0, pts.circle[,2] ) )

presloc <- NodeList[na.omit((match(names(V(g) ), NodeList$V1  ))), ]

#present_i <- unique( as.vector(unlist(as.matrix(  unique(plotddi[order(from,to)]%>%dplyr::select(from,to)) ) )) )
# plot.igraph(g,rescale=FALSE,#vertex.label=V(g)$name,
#             layout = as.matrix(presloc%>%select(-V1))   , xlim =c(-1,1), ylim =c(-1,1),#layout=layout.fruchterman.reingold, 
#             vertices=NodeList[V1%in%presloc$V1],
#             edge.label.color = adjustcolor("black", .5),
#             edge.color="black",edge.width=0.03*E(g)$weight)
# 
# plot.igraph(g,rescale=FALSE,#vertex.label=V(g)$name,
#             layout = as.matrix(presloc%>%select(-V1))   , xlim =c(-1,1), ylim =c(-1,1),#layout=layout.fruchterman.reingold, 
#             vertices=NodeList[V1%in%presloc$V1],
#             edge.label.color =adjustcolor("black", .5),
#             edge.color=adjustcolor("black", .05),edge.width=0.01*E(g)$weight)
# 



# light version
plot.igraph(g,rescale=FALSE,#vertex.label=V(g)$name,
            layout = as.matrix(presloc%>%dplyr::select(-V1))   , xlim =c(-1,1), ylim =c(-1,1),#layout=layout.fruchterman.reingold, 
            vertices=NodeList[V1%in%presloc$V1],
            edge.label.color =adjustcolor("black", .5),
            edge.color=adjustcolor("black", .005),edge.width=0.000000001*E(g)$weight)



# dark version
plot.igraph(g,rescale=FALSE,#vertex.label=V(g)$name,
            layout = as.matrix(presloc%>%select(-V1))   , xlim =c(-1,1), ylim =c(-1,1),#layout=layout.fruchterman.reingold, 
            vertices=NodeList[V1%in%presloc$V1],
            edge.label.color =adjustcolor("black", .9),
            edge.color=adjustcolor("black", 1),edge.width=1e-16*E(g)$weight)










g <- graph.data.frame(CCI_plot%>%dplyr::select(-Patient.Study.ID), directed=TRUE)
E(g)$weight <- exp((log(CCI_plot$TransductionMu)-centVal)/scalsd)
# # specifcy color of nodes
resp_cols<- ggsci::pal_npg("nrc")(2)
if(RRR=="Non-response"){V(g)$color <- resp_cols[1] }else{if(RRR=="Response"){V(g)$color <- resp_cols[2]}}
# circle layout.
n= length(unique(CCI_plot$ReceptorPhenoCelltype)) -1
pts.circle <- t(sapply(1:n,function(r)c(cos(2*r*pi/n),sin(2*r*pi/n))))
NodeList <- data.table(c("Cancer cells", "CD4+ T cells" ,"Adipocytes","Fibroblasts","Normal epithelial cells","Pericytes"  , "Endothelial cells" ,"Macrophages", "B cells","CD8+ T cells" )  , c(0, pts.circle[,1] )  ,  c(0, pts.circle[,2] ) )

presloc<-NodeList[na.omit((match(names(V(g) ), NodeList$V1  ))), ]
plot.igraph(g,rescale=FALSE,#vertex.label=V(g)$name,
            layout = as.matrix(presloc%>%select(-V1))   , xlim =c(-1,1), ylim =c(-1,1),#layout=layout.fruchterman.reingold, 
            vertices=NodeList[V1%in%presloc$V1],        
            edge.label.color =adjustcolor("black", .5),
            edge.color= adjustcolor(CCI_plot$signalCol, .5),
            # edge.color=adjustcolor("black", .05),
            edge.width=0.1*E(g)$weight)

# 
# 