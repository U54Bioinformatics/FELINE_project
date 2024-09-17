rm(list=ls())   
require(rTensor); require(nnTensor); library(abind)
require(data.table); require(dplyr); require(ggplot2); require(tidyr)
require(igraph)
require(RColorBrewer)
require(ggraph)
require(tidygraph)
require(scales)

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



### Discovery cohort analysis of communication connectance  
CCI <- CCI_discovery
#ngroup <- 1
#set.seed(123);  lumrg <- data.table(Pair.Name= unique(CCI$Pair.Name) , grpvar= sample(1:ngroup, length(unique(CCI$Pair.Name) ) ,replace= T ) )
CCI[,scaleTransduction2:=exp(scale(log(TransductionMu),center=F)),by=c("Pair.Name")] 

# nested averaging
average_muln_scaleTransduction <- data.table(CCI%>%
                                               group_by(Patient.Study.ID,LigandPhenoCelltype,ReceptorPhenoCelltype,Day,dynamic_class3,ARM,Treat) %>%
                                               dplyr::summarise(muln_scaleTransduction=mean(log(1+scaleTransduction))) %>% 
                                               group_by(LigandPhenoCelltype,ReceptorPhenoCelltype,Day,dynamic_class3,Treat) %>%
                                               dplyr::summarise(average_muln_scaleTransduction=mean(muln_scaleTransduction)))

# reorder cell types
average_muln_scaleTransduction$LigandPhenoCelltype <-factor(average_muln_scaleTransduction$LigandPhenoCelltype, 
                                                            levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Normal epithelial cells","Cancer cells","B cells","Macrophages","CD8+ T cells","CD4+ T cells"))
average_muln_scaleTransduction$ReceptorPhenoCelltype <-factor(average_muln_scaleTransduction$ReceptorPhenoCelltype, 
                                                              levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Normal epithelial cells","Cancer cells","B cells","Macrophages","CD8+ T cells","CD4+ T cells") )

# calculate pre-treatment average 
initState <- data.table(average_muln_scaleTransduction[Day==0]%>%
                          group_by(LigandPhenoCelltype,ReceptorPhenoCelltype)%>% #Treat, dynamic_class3
                          dplyr::summarise(intiState=mean(average_muln_scaleTransduction)) )

# merge with main dataset
average_muln_scaleTransduction2 <- merge(average_muln_scaleTransduction,initState,by=c("LigandPhenoCelltype", "ReceptorPhenoCelltype")) #,"Treat","dynamic_class3"
# calc fold changes in communication
average_muln_scaleTransduction2[ , lnfoldchange:= (average_muln_scaleTransduction-intiState) ] 

# exclude rarely detected B cell population
CCIcondensed <- average_muln_scaleTransduction2[LigandPhenoCelltype!="B cells"][ReceptorPhenoCelltype!="B cells"]  

# Reorder and update names for ploting
CCIcondensed$LigandPhenoCelltype <-factor(CCIcondensed$LigandPhenoCelltype, 
                                          levels= sort(c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Normal epithelial cells","Cancer cells","Macrophages","CD8+ T cells","CD4+ T cells"))) #"B cells",
CCIcondensed$ReceptorPhenoCelltype <-factor(CCIcondensed$ReceptorPhenoCelltype, 
                                            levels= sort(c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Normal epithelial cells","Cancer cells","Macrophages","CD8+ T cells","CD4+ T cells") )) #"B cells",
setnames(CCIcondensed, old=c("LigandPhenoCelltype", "ReceptorPhenoCelltype", "lnfoldchange"), new= c("from", "to", "weight"))

DiscoveryCommunicationNetwork <- data.table(Cohort="Discovery",CCIcondensed)
#save(CCIcondensed,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/temporary/DiscoveryCommunicationNetwork.RData")



### Validation cohort analysis of communication connectance  
CCI <- CCI_validation
CCI[,scaleTransduction2:=exp(scale(log(TransductionMu),center=F)),by=c("Pair.Name")] 

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
ValidationCommunicationNetwork <- data.table(Cohort="Validation",CCIcondensed)

#save(CCIcondensed,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/temporary/ValidationCommunicationNetwork.RData")



### Join datasets together and fold changes
CCIcondensed <- rbind(DiscoveryCommunicationNetwork,ValidationCommunicationNetwork)
CCIcondensed[ , lnfoldchange:= (average_muln_scaleTransduction-intiState) ] 

# specifcy color of nodes
node_col <- ggsci::pal_npg("nrc")(2)  
# circle layout
n <- length(unique(CCIcondensed$to)) -1
pts.circle <- t(sapply(1:n, function(r)c(cos(2*r*pi/n),sin(2*r*pi/n))))
NodeList <- data.table(c("Cancer cells", "CD4+ T cells" , "Adipocytes", "Fibroblasts", "Normal epithelial cells", "Pericytes"  , "Endothelial cells" ,"Macrophages", "CD8+ T cells" )  , c(0, pts.circle[,1] )  ,  c(0, pts.circle[,2] ) ) #"B cells",
presloc <- data.frame(NodeList)[order(NodeList$V1),]

networks <- data.table(CCIcondensed[Day%in%c(0)][order( Treat,from,to)]%>%group_by(Cohort,from,to,dynamic_class3,Treat)%>%
                         dplyr::summarise(lnfoldchange=mean(lnfoldchange)))

#save(networks, file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/Figure3/Discovery and Validation Communication networks pre treatment.RData")
#load( file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/Figure3/Discovery and Validation Communication networks pre treatment.RData")
#write.csv(networks,"/Users/jason/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/Figure3/SourceData_Figure3_Communication network changes_Output.csv")

#### Viualization of output
tmp <- data.table( networks %>% group_by(Treat,
                                         Cohort,dynamic_class3) %>% 
                     mutate(weight = scale(exp(lnfoldchange))^(0.5) ) )
tmp[,CohortLab2 := " "]
tmp[Cohort=="Validation",CohortLab2 := "  "]
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
  geom_edge_loop( aes(edge_colour=weight,edge_width= 5*(weight),edge_alpha= 5*(weight)) ,  arrow= arrow(length= unit(3.5, "mm" ) , type= "open"),
                  start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" ) ) +
  geom_edge_arc( aes(edge_colour=weight,edge_width= 5*(weight),edge_alpha= 5*(weight)  ) , arrow= arrow(length= unit(3.5, "mm" ) , type= "open" ),
                 start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" )  ) +
  geom_node_point( col= rep( rep(node_col, each= nrow(createdLayout) ) ,  4) , size= 20 ) + 
  geom_node_text(aes(label= name),size=7.5 ) + #,col= "grey80" ) + #, size= 9) +
  facet_edges(paste0(Treatmentlab,CohortLab2)~TumorResponse ,ncol=4 ) + 
  theme_void(base_size= 1.5* 22 ) +
  theme(aspect.ratio = 1, legend.position = "none" ) + 
  scale_edge_width(range=c(0.0,2.5))+
  scale_edge_alpha(range=c(0.0,0.97)) +
  scale_edge_color_gradient2()+
  theme(#aspect.ratio=1,plot.margin=unit(c(10,30,30,30), "points"),
    # panel.spacing=unit(2, "lines") ,
    strip.text=element_text(size=22))

paperfile<- "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
p0 + scale_edge_color_gradient2(low="white", high="black")
p0 + scale_edge_color_gradient2(low="white", high="slategrey")










tmp2 <- data.table( networks %>% group_by(Treat,
                                          Cohort,dynamic_class3) %>% 
                      mutate(weight = scale(exp(lnfoldchange))^(0.5) ) )


tmp2 <- data.table(tmp2%>%group_by(Treat,Cohort,dynamic_class3)%>%mutate(
  mu=mean(lnfoldchange,na.rm=T),
  ucl=sd(lnfoldchange,na.rm=T) ))
tmp2[,weightcol:="Low"]
tmp2[lnfoldchange<=(mu+1.2*ucl),weightcol:="Low"]
tmp2[lnfoldchange >(mu+1.2*ucl),weightcol:="High"]

tmp2[,CohortLab2 := " "]
tmp2[Cohort=="Validation",CohortLab2 := "  "]
tmp2[,Treatmentlab:= "Letrozole alone"]
tmp2[Treat=="CombinationRibo",Treatmentlab:= "Combination ribociclib"]
tmp2[,TumorResponse:="Resistant"]
tmp2[dynamic_class3=="Response",TumorResponse:="Sensitive"]
tmp2[from=="Normal epithelial cells",from:="Diploid \n epithelial \n cells"]
tmp2[to=="Normal epithelial cells",to:="Diploid \n epithelial \n cells"]
tmp2[from=="Endothelial cells",from:="Endothelial \n cells"]
tmp2[to=="Endothelial cells",to:="Endothelial \n cells"]
tmp2[from=="CD8+ T cells",from:="CD8+ \n T cells"]
tmp2[to=="CD8+ T cells",to:="CD8+ \n T cells"]
tmp2[from=="CD4+ T cells",from:="CD4+ \n T cells"]
tmp2[to=="CD4+ T cells",to:="CD4+ \n T cells"]
tmp2[from=="Cancer cells",from:="Cancer \n cells"]
tmp2[to=="Cancer cells",to:="Cancer \n cells"]
tmp2[from=="Macrophages",from:="Myeloid \n cells"]
tmp2[to=="Macrophages",to:="Myeloid \n cells"]
graphdd2 <- as_tbl_graph( tmp2, directed= T)
createdLayout2 <- create_layout(graphdd2, layout= "star")
createdLayout2$x <- presloc$V2
createdLayout2$y <- presloc$V3

p1<-ggraph(createdLayout2 , aes(x= x, y= y , col= dynamic_class3 ) ) + 
  geom_edge_loop( aes(edge_colour=weightcol,edge_width= 5*(weight),
                      edge_alpha= 5*(weight)
  ) ,  
  arrow= arrow(length= unit(3.5, "mm" ) , type= "open"),
  start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" ) ) +
  geom_edge_arc( aes(edge_colour=weightcol,edge_width= 5*(weight),
                     edge_alpha= 5*(weight) 
  ) , 
  arrow= arrow(length= unit(3.5, "mm" ) , type= "open" ),
  start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" )  ) +
  geom_node_point( col= rep( rep(node_col, each= nrow(createdLayout) ) ,  4) , size= 20 ) + 
  geom_node_text(aes(label= name),size=7.5 ) + #,col= "grey80" ) + #, size= 9) +
  facet_edges(paste0(Treatmentlab,CohortLab2)~TumorResponse ,ncol=4 ) + 
  theme_void(base_size= 1.5* 22 ) +
  theme(aspect.ratio = 1, legend.position = "none" ) + 
  scale_edge_width(range=c(0.0,2.5))+
  scale_edge_alpha(range=c(0.0,0.97)) +
  scale_edge_color_manual(values=c("black","grey"))+
  theme(#aspect.ratio=1,plot.margin=unit(c(10,30,30,30), "points"),
    # panel.spacing=unit(2, "lines") ,
    strip.text=element_text(size=22))

