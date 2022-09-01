rm(list=ls())   
library(abind); require(data.table); require(dplyr); require(ggplot2); require(tidyr); require(igraph)

### Load Ligand Receptor database list of Ramilowski et al 2015
load( "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/LigandReceptor/Filtered_Human-2015-Ramilowski-LR-pairs.RData")
LRgenelist <- unique( c(LRpairsFiltered$HPMR.Receptor, LRpairsFiltered$HPMR.Ligand) )
LRpairsFiltered[1:10, 1:7]

# load cell metadata and clinical information
load( file="/Users/jason/Dropbox/PD1 Analysis/Lance/PD1_combined/PD1_paper_cells_analysed.RData")
load( file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/TWISTER example datasets/PD1_paper_cells_analysed.RData")
Meta.dd2[1:10, 1:10]


# read gene expression data and join umap locations to this data or the metadata
RAW  <- fread(file= "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/TWISTER example datasets/PD1.RawCounts.wUMAPCoords.txt" ) 
#RAW  <- fread(file= "/Users/jason/Dropbox/PD1 Analysis/Lance/PD1_combined/Raw Counts/PD1.RawCounts.wUMAPCoords.txt" ) 
RAW[1:10, 1:10]

# visualize cell type heterogeneity at different resolutions
umap_vis_dd <- merge( RAW %>% dplyr::select(Cell.ID, UMAP1, UMAP2) , Meta.dd2 , by= "Cell.ID")
umap_vis_dd$Intermediate_Clusters_ML2<- umap_vis_dd$Intermediate_Clusters_ML
umap_vis_dd$Intermediate_Clusters_ML2<- gsub("_", " ", umap_vis_dd$Intermediate_Clusters_ML2)
#umap_vis_dd$Intermediate_Clusters_ML2<- gsub("Platelets", "Platelet", umap_vis_dd$Intermediate_Clusters_ML2)
umap_vis_dd$Intermediate_Clusters_ML2<- gsub("T Cell CD4","CD4+ T cell",umap_vis_dd$Intermediate_Clusters_ML2)
umap_vis_dd$Intermediate_Clusters_ML2<- gsub("T Cell CD8","CD8+ T cell",umap_vis_dd$Intermediate_Clusters_ML2)
umap_vis_dd$Intermediate_Clusters_ML2<- gsub("NK Cell","NK cell",umap_vis_dd$Intermediate_Clusters_ML2)
umap_vis_dd$Intermediate_Clusters_ML2<- gsub("Dendritic Cell","Dendritic cell",umap_vis_dd$Intermediate_Clusters_ML2)
umap_vis_dd$Intermediate_Clusters_ML2<- gsub("B Cell","B cell",umap_vis_dd$Intermediate_Clusters_ML2)
umap_vis_dd$Intermediate_Clusters_ML2<- gsub("Activated Platelets","Activated platelet",umap_vis_dd$Intermediate_Clusters_ML2)



ggplot(umap_vis_dd[], aes(UMAP1, UMAP2, col= Cell.class.for.normalisation)) + theme_classic() + geom_point(alpha= 0.8, size= 0.5)
ggplot(umap_vis_dd, aes(UMAP1, UMAP2, col= Sub_cluster)) + theme_classic() + geom_point(alpha= 0.8, size= 0.5)
ggplot(umap_vis_dd, aes(UMAP1, UMAP2, col= Cluster)) + theme_classic() + geom_point(alpha= 0.8, size= 0.5) + theme(legend.position= "none")

ggplot(umap_vis_dd[!Intermediate_Clusters_ML2%in%c("Unknown","RBC","Plasmactoid Dendritic cell","Activated platelet","EB cell")], aes(UMAP1, UMAP2, col= Intermediate_Clusters_ML2)) + theme_classic(base_size=32) + geom_point(alpha= 0.8, size= 0.5)+
  scale_color_discrete(name="Cell type") + guides(color=guide_legend(override.aes=list(size=4.5)))+theme(aspect.ratio=1)

paperfile<- "/Users/jason/Dropbox/FELINE Project (1)/Manuscript  Feline immune communication/Figures Communication Project/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Ribo and Letrozole SI Peripheral blood umap.png"),height=10,width=10)


# Settings to select cells from one sample of a tumor
pars <- c(Patient = "HJD33E", TimePoint="C1")


# Extract subset of cells sampled from one patient: meta and CPM data
WhichCells <- Meta.dd2[Patient.ID %in% pars["Patient"]]
RAWsubset <- RAW[Cell.ID %in% WhichCells$Cell.ID,] # RAWsubset[1:10,1:10]
rm(list="RAW")


### Unique units (uu): clusters of phenotypically similar cells (all cell types) and their umap discretization level
uu <- unique( WhichCells %>% dplyr::select( c("Major_cluster", "Intermediate_Clusters_ML" ,"Sub_cluster","Cluster") ) ) 
ggplot(umap_vis_dd[Cell.ID%in%WhichCells$Cell.ID,], aes(UMAP1, UMAP2, col=Cluster))+ theme_classic() + geom_point(alpha=0.8, size=.5) + theme(legend.position="none")


### Extract Ligand and Receptor genes in the L-R database that are present in the CPM data
# Subset Ligand Receptor genes and umap coordinates from CPM data
LRcpm <- RAWsubset[, c("UMAP1", "UMAP2", "Cell.ID", LRgenelist[ LRgenelist %in% names(RAWsubset) ]), with= FALSE] #LRcpm[1:3,1:4]    # select just the ligand receptor gene expression
# List the genes in the LR cm dataset
LRGene.ID <- LRgenelist[ LRgenelist %in% names(RAWsubset) ] 
rm(list= "RAWsubset")


# Identify which receptor and ligand pairs are both represented in the dataset
# copy the dataset as we will add to it an indicator if each ligand/receptor is present and then retain those cases where both are present
LRpairsFiltered_i <- LRpairsFiltered   
LRpairsFiltered_i[, LigandPresent:= 0 ]
LRpairsFiltered_i[, ReceptorPresent:= 0 ]
LRpairsFiltered_i[LRpairsFiltered_i$HPMR.Ligand %in% LRGene.ID, LigandPresent:= 1 ]
LRpairsFiltered_i[LRpairsFiltered_i$HPMR.Receptor %in% LRGene.ID, ReceptorPresent:= 1 ]
# Selecte the genes for which we can calculate communication scores
LRpairsFiltered2_i <- LRpairsFiltered_i[ LigandPresent== 1 & ReceptorPresent== 1 ]
LRpairsFiltered2_i[1:10,]

# For a specific timepoint and patient sample, merge cell metadata and expression of LR genes
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
ggplot(dd3[gene=="IFNG"], aes(y= expression, x= Sub_cluster, col= Major_cluster) ) + geom_point()+ylab("IFNG")
ggplot(dd3[gene=="LRP1"], aes(y= expression, x= Sub_cluster, col= Major_cluster) ) + geom_point()+ylab("LRP1")
ggplot(dd3[gene=="CSF1R"], aes(y= expression, x= Sub_cluster, col= Major_cluster) ) + geom_point()+ylab("CSF1R")


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

dd5[1:10,]

# Composition visualization
ggplot( dd5[gene==dd5$gene[1]] , aes( x=FracSample, y=paste( Patient.ID, Time.Point,sep="_"), fill=Major_cluster)) + geom_bar(stat="identity") +theme_classic()+coord_flip()+ylab("Sample") +xlab("Fraction")
ggplot( dd5[gene==dd5$gene[1]] , aes( x=FracSample, y=paste( Patient.ID, Time.Point,sep="_"), fill=Sub_cluster)) + geom_bar(stat="identity") +theme_classic()+coord_flip()+ylab("Sample")+xlab("Fraction")
ggplot( dd5[gene==dd5$gene[1]] , aes( x=FracSample, y=paste( Patient.ID, Time.Point,sep="_"), fill=key_)) + geom_bar(stat="identity") +theme_classic()+theme(legend.position="none")+coord_flip()+ylab("Sample")+xlab("Fraction")

plotthis<- dd5[gene==dd5$gene[1]][!Intermediate_Clusters_ML%in%c("Activated_Platelets","EB_Cell","Unknown")]
plotthis$Sub_cluster<- gsub("_"," ",plotthis$Sub_cluster)
ggplot( plotthis , aes( x=FracSample/sum(FracSample), y="", fill=Sub_cluster)) + geom_bar(stat="identity") +theme_classic(base_size = 35)+coord_flip()+ylab("")+xlab("Fraction")+
  scale_fill_discrete(name="Cell type subpopulation")

paperfile<- "/Users/jason/Dropbox/FELINE Project (1)/Manuscript  Feline immune communication/Figures Communication Project/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Ribo and Letrozole SI Peripheral blood composition.png"),height=8,width=10)

###Quantify communication (Ligand-Receptor interaction)  
## Look ups for all of the combinations of gene and cell types (discretization class) 
# Two copies of the look up of unique units (uu) where we will modify column names in 2 different ways
lu2 <- uu %>% dplyr::select("Intermediate_Clusters_ML" , "Sub_cluster", "Cluster")
setnames(lu2, old= c("Cluster", "Sub_cluster", "Intermediate_Clusters_ML"), new= c("Receptor", "ReceptorCelltype", "ReceptorPhenoCelltype"))
lu2i <- uu %>% dplyr::select("Intermediate_Clusters_ML" , "Sub_cluster", "Cluster")
setnames(lu2i, old= c("Cluster", "Sub_cluster", "Intermediate_Clusters_ML"), new= c("Ligand", "LigandCelltype", "LigandPhenoCelltype"))
lu2[1:4,]
lu2i[1:4,]

#run communication analysis
Communication <- rbindlist(mclapply( 1:nrow(LRpairsFiltered2_i)   , function(p){   
  # Names of ligand and receptor
  LRnms <- unname( unlist( LRpairsFiltered2_i[p][ , HPMR.Receptor, HPMR.Ligand ] ) )
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
  output <- data.table( merge( LRpair_it_dd, Rmatlong2, by= "key_", all.x=T )  )
  return( output )
  #cat(p);cat("    ")
},mc.cores= detectCores()-2))
addzeros=FALSE
PatientID <- Communication[1]$Patient.ID
cat("Saving output for :        patient "); cat(PatientID); cat("       time     ") ; cat(tau)
#savenm <- paste0("ImmuneCommunicationResults__","PatientID_",PatientID,"__", "Day_",tau,".RData")
saveloc <- "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/ImmuneCommunicationOutput2/"
#save( PatientID, Communication, tau,  uu, addzeros, dd5,file=paste0(saveloc,savenm))
rm(list=c("dd5", "lu2", "lu2i", "Communication"))

filenamesCCI <- list.files(saveloc)

CCI <- rbindlist(lapply(1:length(filenamesCCI), function(ii){
  load(file= paste0(saveloc, filenamesCCI[ii]))
  cat(ii);
  SumComm <- data.table( Communication %>%
                           group_by(Patient.ID, LigandPhenoCelltype, ReceptorPhenoCelltype,    Time.Point , Pair.Name, Responder)%>%
                           dplyr::summarise(TransductionMu= weighted.mean(Signal,countofvalues),  Receptor= weighted.mean(Receptor), Ligandtot= sum(Ligand), Ligand_Ntot= sum(Ligand_N))  )[!( is.na(LigandPhenoCelltype)|is.na(ReceptorPhenoCelltype) ) ]
}))
# Scale data by communication pathway to make data comparable across communication types
CCI[, scaleTransduction:= scale(TransductionMu, center=F), by=c("Pair.Name")] 
CCI[!is.finite(TransductionMu), scaleTransduction:= 0]

CCI$ReceptorPhenoCelltype<- factor(CCI$ReceptorPhenoCelltype,c("Monocyte","Dendritic_Cell" ,"B_Cell","EB_Cell","T_Cell_CD4","T_Cell_CD8","NK_Cell","Activated_Platelets","Unknown"))

ggplot(CCI[ReceptorPhenoCelltype=="T_Cell_CD8"][ReceptorPhenoCelltype!="Unknown"],
       aes(y=log(1+scaleTransduction),x=LigandPhenoCelltype ,col=LigandPhenoCelltype ,fill=LigandPhenoCelltype) ) + theme_classic()+
  geom_violin(scale ="width")+
  geom_point(pch=21,col="black",size=2.5)+
  labs(y="Tumor wide communication to CD8 T cells",x="Signal sender")


ggplot(CCI[][ReceptorPhenoCelltype!="Unknown"],
       aes(y=log(1+scaleTransduction),x=LigandPhenoCelltype ,col=LigandPhenoCelltype ,fill=LigandPhenoCelltype) ) + theme_classic()+
  geom_violin(scale ="width")+
  geom_point(pch=21,col="black",size=2.5)+
  labs(y="Tumor wide communication to CD8 T cells",x="Signal sender")

ggplot(CCI[][ReceptorPhenoCelltype!="Unknown"],
       aes(y=log(1+scaleTransduction),x=ReceptorPhenoCelltype ,col=ReceptorPhenoCelltype ,fill=ReceptorPhenoCelltype) ) + theme_classic()+
  geom_violin(scale ="width")+
  geom_point(pch=21,col="black",size=2.5)+
  labs(y="Tumor wide communication to CD8 T cells",x="Signal receiver")

# 
# ggplot(CCI[Pair.Name=="IFNG_IFNGR1"][ReceptorPhenoCelltype=="T_Cell_CD8"],
#        aes(y=log(1+scaleTransduction),x=LigandPhenoCelltype ,col=LigandPhenoCelltype ,fill=LigandPhenoCelltype) ) + theme_classic()+
#   geom_point(pch=21,col="black",size=2.5) + theme(legend.position="none")+
#   labs(y="Tumor wide IFNG_IFNGR1 communication to CD8 T cells",x="Signal sender")
# ggplot(CCI[Pair.Name=="CSF1_CSF1R"][ReceptorPhenoCelltype%in%"Monocyte"],
#        aes(y=log(1+scaleTransduction),x=LigandPhenoCelltype ,col=LigandPhenoCelltype ,fill=LigandPhenoCelltype) ) + theme_classic()+
#   geom_point(pch=21,col="black",size=2.5) + theme(legend.position="none")+
#   labs(y="Tumor wide IFNG_IFNGR1 communication to CD8 T cells",x="Signal sender")
# 
# 
# ggplot(CCI[ReceptorPhenoCelltype!="Unknown"][Pair.Name=="IL2_IL2RB"][],
#        aes(y=log(1+scaleTransduction),x=ReceptorPhenoCelltype ,col=ReceptorPhenoCelltype ,fill=ReceptorPhenoCelltype) ) + theme_classic()+
#   geom_point(pch=21,col="black",size=2.5) + theme(legend.position="none")+
#   labs(y="Tumor wide LTB_CD40 communication",x="Signal sender")



res1<-rbindlist( lapply(unique(CCI$Pair.Name), function(i){
  
  out <- tryCatch({
    dat<-CCI[Pair.Name==i][!ReceptorPhenoCelltype%in%c("EB_Cell","Unknown")]
    r<-data.table("Pair.Name"=i, coef(summary( lm( scale( log(1+scaleTransduction))~   -1+ReceptorPhenoCelltype ,data=  dat) )),keep.rownames = T)
    setnames(r,old=c("Std. Error",    "t value",  "Pr(>|t|)"), new=c("Std.Error","tval","pval"))
    r2<-r[pval<0.05]
    r2$n<- nrow(dat)
    return(r2)
  },
  error=function(x){
    data.table("Pair.Name"=i, data.table("rn"=NA, "Estimate"=NA ,'Std.Error'=NA, 'tval'=NA,'pval'=NA,'n'=NA) )
  })
return(out)
})) 

res1<-na.omit(res1)
res1$rn<- gsub("ReceptorPhenoCelltype","",res1$rn)
res1$adjpval<- p.adjust(res1$pval)
res1[adjpval<0.05][order(rn)]$pval%>%hist()
res2<-res1[n>10][adjpval<0.05][order(-adjpval)]

res2[Estimate>0]

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

paperfile<- "/Users/jason/Dropbox/FELINE Project (1)/Manuscript  Feline immune communication/Figures Communication Project/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Ribo and Letrozole SI Peripheral blood communication.png"),height=10,width=10)


sort( unique(CCI$Pair.Name)[!grepl("ITG",unique(CCI$Pair.Name))] )


#save(CCI,WhichCells, uu, perIndiv, file ="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Communication output merged/PopulationCommunicationMerged.RData" )

# select data to plot
CCI_plot <- CCI[
  ReceptorPhenoCelltype%in%c("B_Cell","Dendritic_Cell","Monocyte","NK_Cell","T_Cell_CD4","T_Cell_CD8")][
    LigandPhenoCelltype%in%c("B_Cell","Dendritic_Cell","Monocyte","NK_Cell","T_Cell_CD4","T_Cell_CD8")][
      Time.Point==pars["TimePoint"]][][order(LigandPhenoCelltype,ReceptorPhenoCelltype)]
CCI_plot[ReceptorPhenoCelltype=="B_Cell",ReceptorPhenoCelltype:="B cell"]
CCI_plot[LigandPhenoCelltype=="B_Cell",LigandPhenoCelltype:="B cell"]

CCI_plot[ReceptorPhenoCelltype=="Dendritic_Cell",ReceptorPhenoCelltype:="Dendritic cell"]
CCI_plot[LigandPhenoCelltype=="Dendritic_Cell",LigandPhenoCelltype:="Dendritic cell"]

CCI_plot[ReceptorPhenoCelltype=="NK_Cell",ReceptorPhenoCelltype:="NK cell"]
CCI_plot[LigandPhenoCelltype=="NK_Cell",LigandPhenoCelltype:="NK cell"]

CCI_plot[ReceptorPhenoCelltype=="T_Cell_CD4",ReceptorPhenoCelltype:="CD4+ T cell"]
CCI_plot[LigandPhenoCelltype=="T_Cell_CD4",LigandPhenoCelltype:="CD4+ T cell"]

CCI_plot[ReceptorPhenoCelltype=="T_Cell_CD8",ReceptorPhenoCelltype:="CD8+ T cell"]
CCI_plot[LigandPhenoCelltype=="T_Cell_CD8",LigandPhenoCelltype:="CD8+ T cell"]

# construct directed graph
g <- graph.data.frame(CCI_plot%>%dplyr::select(-Patient.ID), directed=TRUE)
# add weights to edges of the graph 
E(g)$weight <- CCI_plot$scaleTransduction  
# specifcy color of nodes
V(g)$color <- ggsci::pal_npg("nrc")(1) 
V(g)$label.cex<-2
V(g)$label.color<-"orange"

# Visualise communication with weighted graph
set.seed(11)
plot.igraph(g,edge.width=0.1*E(g)$weight,edge.color=adjustcolor("black", .5))
 #filename= Ribo and Letrozole SI Peripheral blood communication.pdf           



,layout=layout.fruchterman.reingold)

,rescale=FALSE,
           # layout = as.matrix(presloc%>%select(-V1))   , xlim =c(-1,1), ylim =c(-1,1),#
           layout=layout.fruchterman.reingold,edge.width=0.1*E(g)$weight)
           , 
            #vertices=NodeList[V1%in%presloc$V1],
            edge.label.color =adjustcolor("black", .5),
            edge.color=adjustcolor("black", .5),edge.width=0.1*E(g)$weight)


# circle layout
n= length(unique(CCI_plot$ReceptorPhenoCelltype)) -1
pts.circle <- t(sapply(1:n,function(r)c(cos(2*r*pi/n),sin(2*r*pi/n))))
NodeList <- data.table( c("Dendritic_Cell", "Monocyte", "T_Cell_CD4", "T_Cell_CD8", "NK_Cell", "B_Cell", "EB_Cell" ,"Activated_Platelets", "Unknown" ) , c(0, pts.circle[,1] )  ,  c(0, pts.circle[,2] ) )
presloc <- NodeList[na.omit((match(names(V(g) ), NodeList$V1  ))), ]

# Visualise communication with weighted graph
plot.igraph(g,rescale=FALSE,
            layout = as.matrix(presloc%>%select(-V1))   , xlim =c(-1,1), ylim =c(-1,1),#layout=layout.fruchterman.reingold, 
            vertices=NodeList[V1%in%presloc$V1],
            edge.label.color =adjustcolor("black", .5),
            edge.color=adjustcolor("black", .5),edge.width=0.1*E(g)$weight)




