rm(list=ls())   
library(abind); require(data.table); require(dplyr); require(ggplot2); require(tidyr);require(igraph);require(ggsci)

### Phenotype data
load(file="~/Dropbox/Cancer_pheno_evo/data/FELINE2/PhenotypesOfAllCellTypesAllArms/PhenotypesOfAllCellTypesAllArms.RData") #allphenotypes,UMAPlocs ,UMAPfiles,umapDImRedloc,umapDImRedfiles,nCellTypes,
# Encode the subtypes that have been analysed using UMAP
tmp <- data.table(allphenotypes %>% group_by(key_) %>% slice(1)) %>% dplyr::select(Celltype , Celltype_subtype) %>% unique
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

# extract cell cell interactions from communication data.table
CCI <- rbindlist(lapply(1:length(filenamesCCI), function(ii){
  load(file= paste0(savelocCCI, filenamesCCI[ii]))
  cat(ii);
  SumComm <- data.table( Communication %>%
                           group_by(Patient.Study.ID,LigandPhenoCelltype,ReceptorPhenoCelltype,    Day ,Pair.Name,dynamic_class3,ARM)%>%
                           dplyr::summarise(TransductionMu=weighted.mean(Signal,countofvalues),SumSignal=sum(Signal),   Receptor=weighted.mean(Receptor),Ligandtot=sum(Ligand),Ligand_Ntot=sum(Ligand_N))  )[!( is.na(LigandPhenoCelltype)|is.na(ReceptorPhenoCelltype) ) ]
}))
CCI[, RiboTreated:=TRUE]
CCI[ARM=="A", RiboTreated:=FALSE]

#save(CCI,allphenotypes, uu,perIndiv,file ="~/Dropbox/Cancer_pheno_evo/data/FELINE2/Communication output merged ALL ARMS/PopulationCommunicationMerged ALL ARMS.RData" )


load(file ="~/Dropbox/Cancer_pheno_evo/data/FELINE2/Communication output merged ALL ARMS/PopulationCommunicationMerged ALL ARMS.RData" )
centVal <- mean(log( CCI$SumSignal ))
scalsd<-sd(log(CCI$SumSignal))
RRR <- "Response"
RRR <- "Non-response"
Day_var <- 180
RiboTreat <- TRUE
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
NodeList <- data.table(c("Cancer cells", "CD4+ T cells" ,"Adipocytes","Fibroblasts","Normal epithelial cells","Pericytes"  , "Endothelial cells" ,"Macrophages", "B cells","CD8+ T cells" )  , c(0, pts.circle[,1] )  ,  c(0, pts.circle[,2] ) )

presloc <- NodeList[na.omit((match(names(V(g) ), NodeList$V1  ))), ]


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


