rm(list=ls())   
require(data.table); require(dplyr); require(ggplot2); require(tidyr);require(igraph);require(ggraph);require(tidygraph)
### Phenotype data
load(file="~/Dropbox/Cancer_pheno_evo/data/FELINE2/PhenotypesOfAllCellTypesAllArms/PhenotypesOfAllCellTypesAllArms.RData") #allphenotypes,UMAPlocs ,UMAPfiles,umapDImRedloc,umapDImRedfiles,nCellTypes,
# Encode the subtypes that have been analysed using UMAP
tmp <- data.table( allphenotypes %>% group_by(key_) %>% slice(1) ) %>% dplyr::select(Celltype , Celltype_subtype) %>% unique
tmp[ , PhenoCelltype:= Celltype]
tmp[Celltype_subtype %in% c("CD4+ T cells", "Tregs"), PhenoCelltype:= "CD4+ T cells"]
tmp[Celltype_subtype %in% c("CD8+ T cells", "NK cells"), PhenoCelltype:= "CD8+ T cells"]
tmp[Celltype_subtype %in% c("Macrophages", "DC", "Monocytes"), PhenoCelltype:= "Macrophages"]
tmp[Celltype_subtype %in% c("Vas-Endo", "Lym-Endo","Endothelial cells"), PhenoCelltype:= "Endothelial cells"]
tmp[Celltype_subtype %in% c("Cancer cells"), PhenoCelltype:= "Cancer cells"]
tmp[Celltype_subtype %in% c("Normal epithelial cells"), PhenoCelltype:= "Normal epithelial cells"]
allphenotypes <- merge(allphenotypes %>% dplyr::select(-c( paste0("V", 1:5), "nCount_RNA", "nFeature_RNA", "Infercnv_CNA", "Platform","Timepoint","Sample_p_t","prop_change", "file_string", "day_fact" , "Ribo","TreatLab", "Burden_t0" ,"BurdenTracked" ,"Day0", "DayLastBurd", "Dose_lab","TreatCode","TreatCodeOrd","dynamic_class2","rgr_A","rgr_B")), tmp, by= c("Celltype", "Celltype_subtype") )

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
CCI[,scaleTransduction:=scale(TransductionMu,center=F),by=c("Pair.Name")] #"ReceptorPhenoCelltype",
CCI[!is.finite(TransductionMu),scaleTransduction:=0]
CCI[,scalelnTransductionMu:=scale(log(1+TransductionMu)),by=Pair.Name]

### globe plots
ngroup <- 1
set.seed(123);  lumrg <- data.table(Pair.Name= unique(CCI$Pair.Name) , 
                                    grpvar= sample(1:ngroup, length(unique(CCI$Pair.Name) ) ,replace= T ) )
CCI[,scaleTransduction2:=exp(scale(log(TransductionMu),center=F)),by=c("Pair.Name")] 

average_muln_scaleTransduction <- data.table(CCI[LigandPhenoCelltype!="B cells"][ReceptorPhenoCelltype!="B cells"] %>%
                                               group_by(Patient.Study.ID,LigandPhenoCelltype,ReceptorPhenoCelltype,Day,dynamic_class3,Treat) %>%
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
NodeList <- data.table(c("Cancer cells", "CD4+ T cells" , "Adipocytes", "Fibroblasts", "Normal epithelial cells", "Pericytes"  , "Endothelial cells" ,"Macrophages", "CD8+ T cells" )  , c(0, pts.circle[,1] )  ,  c(0, pts.circle[,2] ) ) #"B cells",
presloc <- data.frame(NodeList)[order(NodeList$V1),]

tmp <- CCIcondensed[Treat== "CombinationRibo"][Day== 180][order( Treat,from,to)]#[Day== 0]
tmp$DayLab <- paste0("Day ", tmp$Day )
graphdd <- as_tbl_graph( tmp, directed= T)
tmp$weight <- scale(tmp$weight)
createdLayout <- create_layout(graphdd, layout= "star")
createdLayout$x <- presloc$V2
createdLayout$y <- presloc$V3

visualiseComNet <- function(D){
tmp <- CCIcondensed[][Day== D][order( Treat,from,to)]
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
graphdd <- as_tbl_graph( tmp, directed= T)
createdLayout <- create_layout(graphdd, layout= "star")
createdLayout$x <- presloc$V2
createdLayout$y <- presloc$V3

p0 <- ggraph(createdLayout , aes(x= x, y= y , col= dynamic_class3 ) ) + 
  geom_edge_loop(  aes(edge_colour=weight,edge_width= (weight),edge_alpha= (weight)) ,  arrow= arrow(length= unit(3, "mm" ) , type= "open"),
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
p1<-p0 + scale_edge_color_gradient2(low="white",high="black")

return(p1)
}

Net0   <- visualiseComNet(0)
Net180 <- visualiseComNet(180)

Net0
ggsave(   file="~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Globe plot ggrraph/Day 0 Globe plot ggrraph1.png",height=10, width = 10)

Net180
ggsave(   file="~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Globe plot ggrraph/Day 180 Globe plot ggrraph1.png",height=10, width = 10)





ngroup <- 10
set.seed(123);  lumrg <- data.table(Pair.Name= unique(CCI$Pair.Name) , 
                                    grpvar= sample(1:ngroup, length(unique(CCI$Pair.Name) ) ,replace=T ) )
CCIcondensed <- data.table(merge( lumrg, CCI, by ="Pair.Name" ) %>% 
                             group_by(Patient.Study.ID, 
                                      dynamic_class3, Day, RiboTreated, grpvar, LigandPhenoCelltype, ReceptorPhenoCelltype) %>%
                             dplyr::summarise( scaleTransduction= exp( median( log( scaleTransduction) ) )    ))

luglode <- data.table( expand.grid(Day = c(0,  180), dynamic_class3 = c("Non-response", "Response")) )


RiboTreat <- TRUE
pdf(file= paste0("~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Globes_Communication globes3/Ribo_AllGlobeCommunicationFinerPairedD2",
                 ".pdf" ), width= 20, height= 20 )
par(mfrow= c(2, 2))
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






RiboTreat <- FALSE
pdf(file= paste0("~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Globes_Communication globes3/Letrozole_AllGlobeCommunicationFinerPairedD2",
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






luglode <- data.table( expand.grid( dynamic_class3 = c("Non-response", "Response") , RiboTreated = c("TRUE", "FALSE")) )

Dy <- 0
pdf(file= paste0("~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Globes_Communication globes3/Day0AllGlobeCommunicationFinerPairedD2",
                 ".pdf" ), width= 20, height= 20 )
par(mfrow= c(2, 2))
for(i in c(1,3,2,4)){  
  CCI_plot <- CCIcondensed[LigandPhenoCelltype!="B cells"][ReceptorPhenoCelltype!="B cells"][dynamic_class3== luglode[i]$dynamic_class3][Day== Dy][RiboTreated== luglode[i]$RiboTreated][order(LigandPhenoCelltype, ReceptorPhenoCelltype)]
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






RiboTreat <- FALSE
pdf(file= paste0("~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Globes_Communication globes3/Letrozole_AllGlobeCommunicationFinerPairedD2",
                 ".pdf" ), width= 20, height= 20 )
par(mfrow= c(2, 2))
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
