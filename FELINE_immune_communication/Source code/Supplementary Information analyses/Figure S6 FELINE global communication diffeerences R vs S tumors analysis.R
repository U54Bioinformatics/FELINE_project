rm(list=ls())   
require(rTensor); require(nnTensor); library(abind)
require(data.table); require(dplyr); require(ggplot2); require(tidyr)
require(igraph)
require(ggraph)
require(tidygraph)

### Load Ligand Receptor database list of Ramilowski et al 2015
LRpairsFiltered <- data.table(read.csv( file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/SI data/FigureS4/SourceData_Figure S4 LigandReceptorPairsRamilowski2015.csv"))
LRgenelist <- unique( c(LRpairsFiltered$HPMR.Receptor, LRpairsFiltered$HPMR.Ligand) )

### Load GO database list of GF receptors
#growthFactorReceptors2 <- read.csv(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/SI data/FigureS4/Figure S4 GeneOntologyGrowthFactorReceptors.csv")$x

### Load some signalling data for downstream analysis
savloc<-"/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/Figure3/"
CCI_all <-  read.csv(file=paste0(savloc, "SourceData_Figure3_CellCummunicationIndexTumorWideSingnaling.csv"))
CCI_discovery <- data.table(CCI_all)[Cohort=="Discovery"]
CCI<-CCI_discovery


# slight grouping to denoise
ngroup <- 10
set.seed(123);  lumrg <- data.table(Pair.Name= unique(CCI$Pair.Name) , 
                                    grpvar= sample(1:ngroup, length(unique(CCI$Pair.Name) ) ,replace=T ) )
#grpvar= 1:length(unique(CCI$Pair.Name) ) )
CCIcondensed <- data.table(merge( lumrg, CCI, by ="Pair.Name" ) %>% 
                             group_by(Patient.Study.ID, 
                                      dynamic_class3, Day, RiboTreated, grpvar, LigandPhenoCelltype, ReceptorPhenoCelltype) %>%
                             dplyr::summarise( scaleTransduction= exp( median( log( scaleTransduction) ) )    ))

luglode <- data.table( expand.grid(Day = c(0,  180), dynamic_class3 = c("Non-response", "Response")) )

RiboTreat <- TRUE
pdf(file= paste0("/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Globes_Communication globes3/Ribo_AllGlobeCommunicationFinerPairedD2",
                 ".pdf" ), width= 20, height= 20 )
par(mfrow= c(2, 2))
#set.seed(123);whcprs<-sample(unique(CCI$Pair.Name),200)
for(i in c(1,3,2,4)){  
  CCI_plot <- CCIcondensed[LigandPhenoCelltype!="B cells"][ReceptorPhenoCelltype!="B cells"][dynamic_class3== luglode[i]$dynamic_class3][Day== luglode[i]$Day][RiboTreated== RiboTreat][order(LigandPhenoCelltype, ReceptorPhenoCelltype)]
  g <- graph.data.frame(CCI_plot %>% dplyr::select(LigandPhenoCelltype, ReceptorPhenoCelltype), directed=TRUE)
  E(g)$weight <- CCI_plot$scaleTransduction    
  # specifcy color of nodes
  resp_cols <- ggsci::pal_npg("nrc")(2)
  if(luglode[i]$dynamic_class3 == "Non-response"){V(g)$color <- resp_cols[1] }else{if(luglode[i]$dynamic_class3 == "Response"){V(g)$color <- resp_cols[2]}}
  # circle layout.
  n= length(unique(CCI_plot$ReceptorPhenoCelltype)) -1
  pts.circle <- t(sapply(1:n, function(r)c(cos(2*r*pi/n),sin(2*r*pi/n))))
  NodeList <- data.table(c("Cancer cells", "CD4+ T cells" ,"Adipocytes","Fibroblasts","Normal epithelial cells","Pericytes"  , "Endothelial cells" ,"Macrophages", "CD8+ T cells" )  , c(0, pts.circle[,1] )  ,  c(0, pts.circle[,2] ) )
  presloc <- NodeList[na.omit((match(names(V(g) ), NodeList$V1  ))), ]
  # option1
  plot.igraph(g,rescale=FALSE,vertex.label.cex=3, vertex.size=18,#vertex.label="",
              layout = as.matrix(presloc %>% dplyr::select(-V1))   , xlim = c(-1,1), ylim = c(-1,1), #layout=layout.fruchterman.reingold,
              vertices= NodeList[V1%in%presloc$V1],
              edge.label.color = adjustcolor("black", .9),
              edge.color = adjustcolor("black", 0.03), 
              edge.width = 10*E(g)$weight,
              edge.arrow.size = 0.02, edge.arrow.width= 0.0005)
}
dev.off()

comOut<-CCIcondensed[LigandPhenoCelltype!="B cells"][ReceptorPhenoCelltype!="B cells"]

savloc <-"/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/SI data/FigureS6/"
write.csv(comOut, file=paste0(savloc,"SourceData_FigureS6_GlobalCommunicationDifferencesRibociclibResistantvsSensitiveTumors.csv"))
