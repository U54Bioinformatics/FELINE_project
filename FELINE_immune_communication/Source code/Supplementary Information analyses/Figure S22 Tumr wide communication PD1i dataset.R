rm(list=ls())   
require(abind); require(data.table); require(dplyr); require(ggplot2); require(tidyr); require(igraph); require(parallel); require("ggalluvial")

savloc<-"/Users/jason/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/SI data/FigureS22/"
# Load cell metadata with fine resolution cell subtype annotation and clinical information
Meta.dd2 <- data.table(read.csv( file=paste0(savloc,"SourceData_FigureS22_PD1iStudyImmuneCellsAnalysed.csv")))
# Read count per million (CPM) gene expression data and join umap locations to this data or the metadata
CPM  <- data.table(read.csv( file=paste0(savloc,"SourceData_FigureS22_PD1iStudyscRNACountsUMAPCoords.csv")))
Meta.dd2[1]
names(CPM)[1]

### Load Ligand Receptor database list of Ramilowski et al 2015
LRpairsFiltered <- data.table(read.csv( file="/Users/jason/Dropbox/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/SI data/FigureS4/SourceData_Figure S4 LigandReceptorPairsRamilowski2015.csv"))
LRgenelist <- unique( c(as.character(LRpairsFiltered$HPMR.Receptor), as.character(LRpairsFiltered$HPMR.Ligand) ))


# Visualize cell type heterogeneity at different resolutions
umap_vis_dd <- merge( CPM %>% dplyr::select(Cell.ID, UMAP1, UMAP2) , Meta.dd2 , by= "Cell.ID")
ggplot(umap_vis_dd, aes(UMAP1, UMAP2, col= Sub_cluster)) + theme_classic() + geom_point(alpha= 0.8, size= 0.5)
#ggplot(umap_vis_dd, aes(UMAP1, UMAP2, col= Cluster)) + theme_classic() + geom_point(alpha= 0.8, size= 0.5) + theme(legend.position= "none")


# Specify  the patient and timepoint codes to select cells from one sample of a tumor
pars <- c(Patient = "HJD33E", TimePoint="C1")


## Subset the data
# Extract subset of cells sampled from one patient: meta and CPM data
WhichCells <- Meta.dd2[Patient.ID %in% pars["Patient"]]
CPMsubset <- CPM[Cell.ID %in% WhichCells$Cell.ID,] # CPMsubset[1:10,1:10]
rm(list="CPM")


### Unique units (uu): clusters of phenotypically similar cells (all cell types) and their umap discretization level
uu <- unique( WhichCells %>% dplyr::select( c("Major_cluster", "Intermediate_Clusters_ML" ,"Sub_cluster","Cluster") ) ) 
ggplot(umap_vis_dd[Cell.ID %in% WhichCells$Cell.ID,], aes(UMAP1, UMAP2, col=Cluster))+ theme_classic() + geom_point(alpha=0.8, size=.5) + theme(legend.position="none")


### Extract Ligand and Receptor genes in the L-R database that are present in the CPM data
# Subset Ligand Receptor genes and umap coordinates from CPM data
LRcpm <- CPMsubset[, c("UMAP1", "UMAP2", "Cell.ID", LRgenelist[ LRgenelist %in% names(CPMsubset) ]), with= FALSE] #LRcpm[1:3,1:4]    # select just the ligand receptor gene expression
# List the genes in the LR cm dataset
LRGene.ID <- LRgenelist[ LRgenelist %in% names(CPMsubset) ] 
rm(list= "CPMsubset")


# Identify which receptor and ligand pairs are both represented in the dataset
# copy the dataset as we will add to it an indicator if each ligand/receptor is present and then retain those cases where both are present
LRpairsFiltered_i <- LRpairsFiltered   
LRpairsFiltered_i[, LigandPresent:= 0 ]
LRpairsFiltered_i[, ReceptorPresent:= 0 ]
LRpairsFiltered_i[LRpairsFiltered_i$HPMR.Ligand %in% LRGene.ID, LigandPresent:= 1 ]
LRpairsFiltered_i[LRpairsFiltered_i$HPMR.Receptor %in% LRGene.ID, ReceptorPresent:= 1 ]
# Selecte the ligand-receptor gene pairs for which we can calculate communication scores because we have both genes
LRpairsFiltered2_i <- LRpairsFiltered_i[ LigandPresent== 1 & ReceptorPresent== 1 ]
LRpairsFiltered2_i[1:10,]

# determine the number of pathways
nrow(LRpairsFiltered2_i)


# Specify timepoint
tau <- pars["TimePoint"]
# Phenotype classifications of all cells in the sample of the tumor at this timepoint
phenotypes_i <- WhichCells[Time.Point == tau]
## Extract clinical data of the patient for merging
clin_i <- phenotypes_i[1, ] %>% dplyr::select(Patient.ID, Time.Point, Responder)

# Merge cell metadata (annotaitons and umap location), sample information and expression of Ligand-Receptor genes
dd2 <-  merge(phenotypes_i, LRcpm, by= "Cell.ID")
# Count the total number of cells in this sample
dd2[ ,samplesize_it:= nrow(dd2) ]   #dd2[,1:40]
# Gather genes into long format and remove wide copy to save memory
dd3 <- data.table( gather( dd2, gene, expression, all_of(LRGene.ID) ) )
rm(list= c("dd2"))


# Visualise specific genes as required
ggplot(dd3[gene=="IFNG"], aes(y= expression, x= Sub_cluster, col= Major_cluster) ) + geom_point()+ylab("IFNG") + theme_classic()+theme(axis.text.x = element_text(angle=90))
ggplot(dd3[gene=="LRP1"], aes(y= expression, x= Sub_cluster, col= Major_cluster) ) + geom_point()+ylab("LRP1") + theme_classic()+theme(axis.text.x = element_text(angle=90))
ggplot(dd3[gene=="CSF1R"], aes(y= expression, x= Sub_cluster, col= Major_cluster) ) + geom_point()+ylab("CSF1R") + theme_classic()+theme(axis.text.x = element_text(angle=90))



### Summarise the average expression of genes in each cell discretization class      
grps <- c("gene", "samplesize_it", "Major_cluster", "Intermediate_Clusters_ML" , "Sub_cluster", "Cluster")   
dd4 <- data.table( dd3 %>% dplyr::group_by_(.dots = grps) %>% dplyr::summarise( expression_bar:= mean(expression), countofvalues = n() ) )
rm(list= "dd3")


## Merge clinical data with average expression data
dd5 <- data.table(clin_i, dd4  )
rm(list= c("dd4"))

# Add the fraction that each celltype contributes to the total tumor sample
dd5[, FracSample:= countofvalues/samplesize_it ]
# Specify that the Cluster column is the key to work with
dd5$key_ <- as.character(dd5$Cluster)

# Composition visualization
ggplot( dd5[gene== dd5$gene[1]] , aes( x= FracSample, y= paste( Patient.ID, Time.Point,sep="_"), fill= Major_cluster)) + geom_bar(stat= "identity") + theme_classic() + coord_flip() + ylab("Sample") +xlab("Fraction")
ggplot( dd5[gene== dd5$gene[1]] , aes( x= FracSample, y= paste( Patient.ID, Time.Point,sep="_"), fill= Sub_cluster)) + geom_bar(stat= "identity") +theme_classic() + coord_flip() + ylab("Sample")+ xlab("Fraction")

#Run communication analysis
## Pre-define look ups for all of the combinations of gene and cell types (discretization class) 
# Two copies of the look up of unique units (uu) where we will modify column names in 2 different ways
lu2 <- uu %>% dplyr::select("Intermediate_Clusters_ML" , "Sub_cluster", "Cluster")
setnames(lu2, old= c("Cluster", "Sub_cluster", "Intermediate_Clusters_ML"), new= c("Receptor", "ReceptorCelltype", "ReceptorPhenoCelltype"))
lu2i <- uu %>% dplyr::select("Intermediate_Clusters_ML" , "Sub_cluster", "Cluster")
setnames(lu2i, old= c("Cluster", "Sub_cluster", "Intermediate_Clusters_ML"), new= c("Ligand", "LigandCelltype", "LigandPhenoCelltype"))
lu2[1:4,]
lu2i[1:4,]

# Perform communication analysis
Communication <- rbindlist(mclapply( 1:nrow(LRpairsFiltered2_i)   , function(p){   
  # Names of ligand and receptor
  LRnms <- as.character(unname( unlist( LRpairsFiltered2_i[p][ , HPMR.Receptor, HPMR.Ligand ] ) ))
  # Select expression of the pair in each cell subtype
  LRpair_it_dd <-  data.table( dd5[gene %in% LRnms ] %>% spread(gene, expression_bar) ,  LRpairsFiltered2_i[p] %>% dplyr::select(HPMR.Receptor, HPMR.Ligand, Pair.Name))
  setnames(LRpair_it_dd , old= LRnms, new= c("Ligand", "Receptor"))
  setcolorder(LRpair_it_dd,c( names(dd5[1] %>% dplyr::select(-c("gene", "expression_bar"))), "Ligand", "Receptor", "HPMR.Receptor", "HPMR.Ligand", "Pair.Name" )  )
  
  # Calculate the expression of the signaller cell type by multiplying single cell average expression by the cell number
  LRpair_it_dd[ , Ligand_N:= Ligand * FracSample]
  
  # Caclulate ligand-receptor signalling between each cell type: outer product matrix -> Signaler on the cols and receiver cell class on the rows
  Rmat <- LRpair_it_dd[, Receptor] %*% t( unlist( LRpair_it_dd[, "Ligand_N"] ) )
  rownames(Rmat) <- colnames(Rmat) <- LRpair_it_dd$key_
  RmatperSignaller <- LRpair_it_dd[, Receptor] %*% t( unlist( LRpair_it_dd[, "Ligand"] ) )
  rownames(RmatperSignaller) <- colnames(RmatperSignaller) <- LRpair_it_dd$key_
  
  # Marginalise signalling matrix to calculate signal transduction to cells of each receiver cell type
  LRpair_it_dd[ , Transduction := rowSums(Rmat) ]
  LRpair_it_dd[ , TransductionperSignaller := rowSums(RmatperSignaller) ]
  
  # Reformat the ligand-receptor signalling matrix into a long dataframe
  Rmatlong <- data.table( gather(as.data.table(Rmat, keep.rownames = T), Ligand, Signal, -1) )[Signal > 0]
  setnames(Rmatlong, old= "rn", new= "Receptor")
  
  # Merge ligand-receptor signalling with cell type information
  Rmatlong2 <- merge( merge(Rmatlong, lu2, by= "Receptor"), lu2i , by= "Ligand")
  # Merge ligand-receptor signalling with transduction data and all clinical information
  setnames(Rmatlong2, old= c("Ligand", "Receptor"),new= c("key_signaller", "key_"))
  output <- data.table( merge( LRpair_it_dd, Rmatlong2, by= "key_", all.x= T )  )
  return( output )
  #cat(p);cat("    ")
},mc.cores= detectCores()-2))
PatientID <- Communication[1]$Patient.ID
cat("Saving output for :        patient "); cat(PatientID); cat("       time     ") ; cat(tau)
#savenm <- paste0("ImmuneCommunicationResults__","PatientID_",PatientID,"__", "Day_",tau,".RData")
#saveloc <- "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/ImmuneCommunicationOutput2/"
#save( PatientID, Communication, tau,  uu, dd5,file=paste0(saveloc,savenm))


# For each sample in the communication output folder, calulate the strength of communciation between each cell type via each communication pathway
filenamesCCI <- list.files(saveloc)
CCI <- rbindlist(lapply(1,#:length(filenamesCCI), 
                        function(ii){
  #load(file= paste0(saveloc, filenamesCCI[ii]))
  #cat(ii);
  SumComm <- data.table( Communication %>%
                           group_by(Patient.ID, LigandPhenoCelltype, ReceptorPhenoCelltype,    Time.Point , Pair.Name, Responder)%>%
                           dplyr::summarise(TransductionMu= weighted.mean(Signal,countofvalues),  Receptor= weighted.mean(Receptor), Ligandtot= sum(Ligand), Ligand_Ntot= sum(Ligand_N))  )[!( is.na(LigandPhenoCelltype)|is.na(ReceptorPhenoCelltype) ) ]
}))
# Scale TME wide communication data (TransductionMu) for each communication pathway to make data comparable across communication types
CCI[, scaleTransduction:= scale(TransductionMu, center=F), by= c("Pair.Name")] 
CCI[!is.finite(TransductionMu), scaleTransduction:= 0]
#save(CCI,WhichCells, uu, perIndiv, file ="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Communication output merged/PopulationCommunicationMerged.RData" )
CCI[1,]

lrpairs1<-c(
  # T cells
  "IL15_IL2RB",
  "IL7_IL7R" ,
  "CCL5_CCR5",
  "CCL28_CCR10",
  "CX3CL1_CX3CR1",# NK cells
  # Monocytes
  "TNF_TNFRSF1A",
  "CSF1_CSF1R",
  # Dendritic cells
  "LTA_LTBR",
  #B cells
  "LTB_CD40",
  # activated platelets
  "TGFB3_TGFBR1"
)

res3<-data.table( CCI %>%spread(ReceptorPhenoCelltype,scaleTransduction,fill=0)%>%gather(ReceptorPhenoCelltype,scaleTransduction,unique(CCI$ReceptorPhenoCelltype)) )[][
  Pair.Name%in%lrpairs1][!ReceptorPhenoCelltype%in%c("Unknown","EB_Cell")]


res3$Pair.Name<- factor(res3$Pair.Name, levels=c(lrpairs1) )
res3$ReceptorPhenoCelltype<- factor(res3$ReceptorPhenoCelltype,rev(c("T_Cell_CD4","T_Cell_CD8","NK_Cell","Monocyte","Dendritic_Cell" ,"B_Cell","Activated_Platelets")))


#ggplot(res3,
#  aes(y=log(1+scaleTransduction),x=ReceptorPhenoCelltype ,col=ReceptorPhenoCelltype ,fill=ReceptorPhenoCelltype) ) + theme_classic()+
#  geom_violin()+ geom_point(pch=21,col="black",size=2.5) +  theme(legend.position="none")+
#  labs(y="Tumor wide LTB_CD40 communication",x="Signal sender")+facet_wrap(~Pair.Name)+coord_flip()

res3$ReceptorPhenoCelltype<- gsub("_"," ",res3$ReceptorPhenoCelltype)
res3$ReceptorPhenoCelltype<- gsub("T Cell CD4","CD4+ T cell",res3$ReceptorPhenoCelltype)
res3$ReceptorPhenoCelltype<- gsub("T Cell CD8","CD8+ T cell",res3$ReceptorPhenoCelltype)
res3$ReceptorPhenoCelltype<- gsub("NK Cell","NK cell",res3$ReceptorPhenoCelltype)
res3$ReceptorPhenoCelltype<- gsub("Dendritic Cell","Dendritic cell",res3$ReceptorPhenoCelltype)
res3$ReceptorPhenoCelltype<- gsub("B Cell","B cell",res3$ReceptorPhenoCelltype)
res3$ReceptorPhenoCelltype<- gsub("Activated Platelets","Activated platelet",res3$ReceptorPhenoCelltype)



res3$Pair.Name<- factor(res3$Pair.Name, levels=rev( c("IL7_IL7R","CCL28_CCR10", "CCL5_CCR5", "CX3CL1_CX3CR1","IL15_IL2RB",
                                                      "TNF_TNFRSF1A", "CSF1_CSF1R",
                                                      "LTA_LTBR",
                                                      "LTB_CD40",
                                                      "TGFB3_TGFBR1") ) )
res3$ReceptorPhenoCelltype<- factor(res3$ReceptorPhenoCelltype,(c("CD4+ T cell","CD8+ T cell","NK cell","Monocyte","Dendritic cell" ,"B cell","Activated platelet")))

ggplot(res3[ReceptorPhenoCelltype!="Activated platelet"],  aes(y=log(1+scaleTransduction),x=Pair.Name ,col=Pair.Name ,fill=Pair.Name) ) + theme_classic(base_size=20)+
  geom_violin(scale="width")+ geom_point(pch=21,col="black",size=3.5) +  theme(legend.position="none")+
  labs(y="Tumor wide communication strength",x="Signal receiver")+facet_wrap(~ReceptorPhenoCelltype)+coord_flip() +
  theme(aspect.ratio=1)

