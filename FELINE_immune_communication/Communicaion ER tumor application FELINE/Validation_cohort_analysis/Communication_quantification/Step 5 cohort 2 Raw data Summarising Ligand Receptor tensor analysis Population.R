
rm(list=ls())   
##devtools::install_github("rikenbit/nnTensor") #install.packages("rTensor")
#install.packages("https://cran.r-project.org/src/contrib/Archive/rTensor/rTensor_1.4.tar.gz", repo=NULL, type="source")
require(rTensor); require(nnTensor); library(abind)
require(data.table); require(dplyr); require(ggplot2); require(tidyr)
require(igraph);require(ggsci)
### Load clinical data
#load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE/Feline clinical input.RData") #save(Clin_resp_dd,Clin_resp_dd_class,Clin_resp_dd_classAdd,patientCode_LU,res_file_46,response_code_dd,response_dd,file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE/Feline clinical input.RData") #gsea_path<-"~/Dropbox/FELINE Project/Data_analysis/scRNA/05_ssGSEA_score/Signature_c2_hallmark/results/"   #"~/Dropbox/FELINE/Data_share/Modeling_Data/pathway_seperated_files/Data_gene_count_per_celltype_model_zinbwave_ssGSEA/FEL001043/"

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
allphenotypes <- merge(allphenotypes %>%dplyr::select(-c( paste0("V", 1:5), "nCount_RNA", "nFeature_RNA",  "Platform","Sample_p_t", "file_string", "day_fact")), tmp, by= c("Celltype", "Celltype_subtype") )

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

#summarytable <- read.csv("/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/TensorAnalysisOutput/Ligand Receptor list fro NTD model of population cellcell communication_AB.csv")

# extract cell cell interactions from communication data.table
CCI <- rbindlist(lapply(1:length(filenamesCCI), function(ii){
  load(file= paste0(savelocCCI, filenamesCCI[ii]))
  cat(ii);
  SumComm <- data.table( Communication %>%
                           group_by(Patient.Study.ID,LigandPhenoCelltype,ReceptorPhenoCelltype,    Day ,Pair.Name,dynamic_class3,ARM)%>%
                           dplyr::summarise(TransductionMu=weighted.mean(Signal,countofvalues),SumSignal=sum(Signal),   Receptor=weighted.mean(Receptor),Ligandtot=sum(Ligand),Ligand_Ntot=sum(Ligand_N))  )[!( is.na(LigandPhenoCelltype)|is.na(ReceptorPhenoCelltype) ) ]
  #dplyr::summarise(TransductionMu=mean(Transduction),  Receptor=mean(Receptor),Ligandtot=sum(Ligand),Ligand_Ntot=sum(Ligand_N))  )[!( is.na(LigandPhenoCelltype)|is.na(ReceptorPhenoCelltype) ) ]
}))
CCI[, RiboTreated:=TRUE]
CCI[ARM=="A", RiboTreated:=FALSE]

#save(CCI,allphenotypes, uu,perIndiv,file ="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/cohort2 Communication output merged/PopulationCommunicationMerged ALL ARMS cohort2.RData" )

load(file ="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/cohort2 Communication output merged/PopulationCommunicationMerged ALL ARMS cohort2.RData" )
centVal <- mean(log( CCI$SumSignal ))
scalsd<-sd(log(CCI$SumSignal))
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








# 
# 
# ### Color by function
# load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/GeneOntology/GeneAnnotation.RData") #my_gene_classProc,my_gene_class,dendprune,my_hclust_allg,receptorFunction, receptorFunctionScore1,
# #centVal <- mean(log(CCI$TransductionMu))
# #scalsd <- sd(log(CCI$TransductionMu))
# RRR <-"Non-response"
# Day_var <- 180
# 
# plotddi <- CCI[dynamic_class3==RRR][Day==Day_var]
# setnames(plotddi,old=c("LigandPhenoCelltype","ReceptorPhenoCelltype"),new=c("from","to"))
# subsetdd <- plotddi[TransductionMu>quantile(TransductionMu,0.8)][order(from,to)]
# subsetdd[, c("Ligand", "Receptor") := tstrsplit(Pair.Name, "_", fixed=TRUE)]
# subsetdd[,lnTransductionMu:=log(TransductionMu)]
# # names of signalling types
# bioproc <- unique(my_gene_classProc$Process)
# 
# # central tendancies to scale data and arrow sizes for plotting
# centVal <- mean(log(CCI$TransductionMu))
# scalsd <- sd(log(CCI$TransductionMu))
# 
# # subsetto plot
# subst<-5000
# # plot all the communication processes.
# par(mfrow=c(3,3));   
# for(i in 1:length(unique(my_gene_classProc$Process))){
#   subsetdd2 <- merge(subsetdd,my_gene_classProc[Process%in%bioproc[i]] ,by="Receptor")
#   # subsetdd2_0 <- merge(subsetdd,my_gene_classProc[Process%in%bioproc[i]] ,by="Receptor")
#   # if(nrow(subsetdd2_0)>subst){
#   #   subsetdd2 <- data.table(subsetdd2_0  %>%slice(sample(1:nrow(subsetdd2_0),subst)))
#   # }else{
#   #   subsetdd2 <- data.table(subsetdd2_0 )
#   # }
#   
#   g <- graph.data.frame(subsetdd2[order(from,to)]%>%dplyr::select(-c(Receptor,Patient.Study.ID)), directed=TRUE)
#   E(g)$weight <- exp((subsetdd2$lnTransductionMu-centVal)/scalsd)   #subsetdd2$val
#   #is.weighted((g))
#   
#   # specifcy color of nodes
#   resp_cols<- ggsci::pal_npg("nrc")(2)
#   if(RRR=="Non-response"){V(g)$color <- resp_cols[1] }else{if(RRR=="Response"){V(g)$color <- resp_cols[2]}}
#   
#   # circle layout
#   n= length(unique(CCI$ReceptorPhenoCelltype)) -1
#   pts.circle <- t(sapply(1:n,function(r)c(cos(2*r*pi/n),sin(2*r*pi/n))))
#   NodeList <- data.table(c("Cancer cells", "CD4+ T cells" ,"Adipocytes","Fibroblasts","Normal epithelial cells","Pericytes"  , "Endothelial cells" ,"Macrophages", "B cells","CD8+ T cells" )  , c(0, pts.circle[,1] )  ,  c(0, pts.circle[,2] ) )
#   presloc <- NodeList[na.omit((match(names(V(g) ), NodeList$V1  ))), ]
#   
#   # colour palette
#   colPaledge <- ggsci::pal_jco("default")( length(my_gene_classProc$Process%>%unique()) )
#   names(colPaledge) <- my_gene_classProc$Process%>%unique()
#   
#   # Plot graph
#   plot.igraph(g,rescale=FALSE,#vertex.label=V(g)$name,
#               layout = as.matrix(presloc%>%select(-V1))   , xlim =c(-1,1), ylim =c(-1,1),#layout=layout.fruchterman.reingold, 
#               vertices=NodeList[V1%in%presloc$V1],
#               edge.label.color =adjustcolor("black", .5),
#               edge.arrow.size=0.1,
#               edge.color=adjustcolor(colPaledge[ as.character(subsetdd2$Process)]  , .05),edge.width=1e-20*E(g)$weight )
# }
# plot.new();plot.new()
# legend("topleft",legend=unique(my_gene_classProc$Process),pch=17,pt.cex=3,col=  ggsci::pal_jco("default")( length(my_gene_classProc$Process%>%unique()) ) )
# 
# 
# 
# 
# 
# 
# 
# 
# RRR<-"Response"
# RRR<-"Non-response"
# Day_var=180
# #CCI_plot <- CCI[dynamic_class3==RRR][Day==Day_var][Pair.Name%in%summarytable$key_][order(LigandPhenoCelltype,ReceptorPhenoCelltype)]
# CCI_plot <- CCI[CCI>centVal][dynamic_class3==RRR][Day==Day_var][][order(LigandPhenoCelltype,ReceptorPhenoCelltype)]
# g <- graph.data.frame(CCI_plot%>%dplyr::select(-Patient.Study.ID), directed=TRUE)
# E(g)$weight <- exp((log(CCI_plot$TransductionMu)-centVal)/scalsd)
# is.weighted((g))
# 
# # specifcy color of nodes
# resp_cols<- ggsci::pal_npg("nrc")(2)
# if(RRR=="Non-response"){V(g)$color <- resp_cols[1] }else{if(RRR=="Response"){V(g)$color <- resp_cols[2]}}
# 
# 
# # circle layout.
# n= length(unique(CCI_plot$ReceptorPhenoCelltype)) -1
# pts.circle <- t(sapply(1:n,function(r)c(cos(2*r*pi/n),sin(2*r*pi/n))))
# 
# #NodeList <- data.table(sort( unique(SignalReceiverRes$ReceptorPhenoCelltype) ), c(4,9.5,5,8,9,4,2,9,1,1.5 ) /5-1 ,c(9,5,5,8,6,1,7,2,4,1.5)  /5-1 )
# NodeList <- data.table(c("Cancer cells", "CD4+ T cells" ,"Adipocytes","Fibroblasts","Normal epithelial cells","Pericytes"  , "Endothelial cells" ,"Macrophages", "B cells","CD8+ T cells" )  , c(0, pts.circle[,1] )  ,  c(0, pts.circle[,2] ) )
# 
# presloc <- NodeList[na.omit((match(names(V(g) ), NodeList$V1  ))), ]
# 
# #present_i <- unique( as.vector(unlist(as.matrix(  unique(plotddi[order(from,to)]%>%dplyr::select(from,to)) ) )) )
# # plot.igraph(g,rescale=FALSE,#vertex.label=V(g)$name,
# #             layout = as.matrix(presloc%>%select(-V1))   , xlim =c(-1,1), ylim =c(-1,1),#layout=layout.fruchterman.reingold, 
# #             vertices=NodeList[V1%in%presloc$V1],
# #             edge.label.color = adjustcolor("black", .5),
# #             edge.color="black",edge.width=0.03*E(g)$weight)
# # 
# # plot.igraph(g,rescale=FALSE,#vertex.label=V(g)$name,
# #             layout = as.matrix(presloc%>%select(-V1))   , xlim =c(-1,1), ylim =c(-1,1),#layout=layout.fruchterman.reingold, 
# #             vertices=NodeList[V1%in%presloc$V1],
# #             edge.label.color =adjustcolor("black", .5),
# #             edge.color=adjustcolor("black", .05),edge.width=0.01*E(g)$weight)
# # 
# 
# 
# 
# # light version
# plot.igraph(g,rescale=FALSE,#vertex.label=V(g)$name,
#             layout = as.matrix(presloc%>%select(-V1))   , xlim =c(-1,1), ylim =c(-1,1),#layout=layout.fruchterman.reingold, 
#             vertices=NodeList[V1%in%presloc$V1],
#             edge.label.color =adjustcolor("black", .5),
#             edge.color=adjustcolor("black", .005),edge.width=0.000000001*E(g)$weight)
# 
# 
# 
# # dark version
# plot.igraph(g,rescale=FALSE,#vertex.label=V(g)$name,
#             layout = as.matrix(presloc%>%select(-V1))   , xlim =c(-1,1), ylim =c(-1,1),#layout=layout.fruchterman.reingold, 
#             vertices=NodeList[V1%in%presloc$V1],
#             edge.label.color =adjustcolor("black", .9),
#             edge.color=adjustcolor("black", 1),edge.width=1e-16*E(g)$weight)
# 
