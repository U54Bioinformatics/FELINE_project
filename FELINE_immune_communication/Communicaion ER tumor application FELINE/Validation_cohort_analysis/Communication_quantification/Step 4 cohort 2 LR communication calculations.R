rm(list=ls())
require(data.table)
require(dplyr)
require(ggplot2)
require(tidyr)
### Load clinical data
#load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE/Feline clinical input.RData") #save(Clin_resp_dd,Clin_resp_dd_class,Clin_resp_dd_classAdd,patientCode_LU,res_file_46,response_code_dd,response_dd,file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE/Feline clinical input.RData") #gsea_path<-"~/Dropbox/FELINE Project/Data_analysis/scRNA/05_ssGSEA_score/Signature_c2_hallmark/results/"   #"~/Dropbox/FELINE/Data_share/Modeling_Data/pathway_seperated_files/Data_gene_count_per_celltype_model_zinbwave_ssGSEA/FEL001043/"

### Phenotype data
#load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/PhenotypesOfAllCellTypesAllArms/PhenotypesOfAllCellTypesAllArms.RData") #allphenotypes,UMAPlocs ,UMAPfiles,umapDImRedloc,umapDImRedfiles,nCellTypes,
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


### Load counts per million read data for a specific patient's samples and gen Ligand+Receptor genes
#CPMlocs <- "/Users/jason/Dropbox/FELINE Project/Data_analysis/FELINE_data_folder/scRNA_count_CPM/output/"
CPMlocs1 <- "/Users/jason/Dropbox/FELINE Project/FELINE Cohort 2/DNANexusCopy/CPM/HQ/"
CPMlocs2 <- "/Users/jason/Dropbox/FELINE Project/FELINE Cohort 2/DNANexusCopy/CPM/HQLQ/"

#CPMfiles <- list.files(CPMlocs)[ grep("CPM", list.files(CPMlocs) ) ]
CPMfiles1 <- list.files(CPMlocs1)[ grep("CPM", list.files(CPMlocs1) ) ]
CPMfiles2 <- list.files(CPMlocs2)[ grep("CPM", list.files(CPMlocs2) ) ]
n10Xpats <- length(CPMfiles1)


# settings
shouldscale <- FALSE ## should cpm be scaled
saveloc <- "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/CommunicationOutputAllArmsCohort2/"

## Select a patient index to load cpm data
for(i in 1:length(CPMfiles1)){  cat("cell type index    ");cat(i);cat("    loading data    ") #i= 2 #patient index

  ## Load gene expression of ligands and receptors
  cpm_i <- data.table( fread( paste0(CPMlocs1, CPMfiles1[i]) ) )   # load full gene expression
  names(cpm_i)[1]<-"Gene.ID"
  LRpairsFiltered_i <- LRpairsFiltered                      # copy ligand receptor list for this sample 
  LRcpm <- cpm_i[Gene.ID %in% LRgenelist, ] #LRcpm[1:3,1:4]    # select just the ligand receptor gene expression
  rm(list= "cpm_i")
  LRGene.ID0 <- LRcpm$Gene.ID
  
  # Identify which receptor and ligand pairs are both represented in the dataset
  LRpairsFiltered_i[, LigandPresent:= 0 ]
  LRpairsFiltered_i[, ReceptorPresent:= 0 ]
  LRpairsFiltered_i[LRpairsFiltered_i$HPMR.Ligand %in% LRGene.ID0, LigandPresent:= 1 ]
  LRpairsFiltered_i[LRpairsFiltered_i$HPMR.Receptor %in% LRGene.ID0, ReceptorPresent:= 1 ]
  LRpairsFiltered2_i <- LRpairsFiltered_i[ LigandPresent==1 & ReceptorPresent==1 ]
  LRcpm <- LRcpm[Gene.ID %in% unique(c(LRpairsFiltered2_i$HPMR.Ligand, LRpairsFiltered2_i$HPMR.Receptor)), ] #LRcpm[1:3,1:4]    # select just the ligand receptor gene expression
  LRGene.ID <- LRcpm$Gene.ID
  # Transpose the ligand receptor gene expression data
  if(shouldscale==T){
    LRtransposedHQ <- as.data.table( scale( t(LRcpm[Gene.ID%in%LRGene.ID,-1]) , center = FALSE) , keep.rownames = T)
  }else{
    LRtransposedHQ <- as.data.table( t(LRcpm[Gene.ID%in%LRGene.ID,-1])  , keep.rownames = T)
  }
  
  LRtransposedHQ[ is.na(LRtransposedHQ) ] <- 0
  names(LRtransposedHQ) <- c("Cell.ID", unname( unlist( LRcpm[,1] ) ) )
  rm(list="LRcpm")
  
 # LRtransposedHQ[1:10,1:3]
  
  
  cpm_iLQ <- data.table( fread( paste0(CPMlocs2, CPMfiles2[i]) ) )   # load full gene expression
  names(cpm_iLQ)[1]<-"Gene.ID"
  cpm_iLQ <- cpm_iLQ[Gene.ID %in% LRGene.ID, ]
  #rownames(cpm_iLQ) <-cpm_iLQ$Gene.ID
  #cpm_iLQ[1:10,1:3]
  LRtransposedLQ <- as.data.table( t(cpm_iLQ[,-1])  , keep.rownames = T)
  colnames(LRtransposedLQ) <- c("Cell.ID",cpm_iLQ$Gene.ID)
  
  setcolorder(LRtransposedHQ, c("Cell.ID",sort(LRGene.ID) ))
  setcolorder(LRtransposedLQ, c("Cell.ID",sort(LRGene.ID) ))
  
  LRtransposed <- rbind(LRtransposedHQ, LRtransposedLQ)
  save(LRtransposed, file= paste0("/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/cohort2LRIntermediateFiles/tmp",i,"file.RData"))
}  


LRtransposed <-  rbindlist(lapply( 1:length(CPMfiles1) , function(i){
  load( file= paste0("/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/cohort2LRIntermediateFiles/tmp",i,"file.RData"))
  return(LRtransposed)
  }) )
  
lookupSamples <- unique(allphenotypes%>%dplyr::select(Day,orig.ident))

LRpairsFiltered_i <- LRpairsFiltered                      # copy ligand receptor list for this sample 
LRpairsFiltered_i[, LigandPresent:= 0 ]
LRpairsFiltered_i[, ReceptorPresent:= 0 ]
LRpairsFiltered_i[LRpairsFiltered_i$HPMR.Ligand %in% names(LRtransposed)[-1], LigandPresent:= 1 ]
LRpairsFiltered_i[LRpairsFiltered_i$HPMR.Receptor %in% names(LRtransposed)[-1], ReceptorPresent:= 1 ]
LRpairsFiltered2_i <- LRpairsFiltered_i[ LigandPresent==1 & ReceptorPresent==1 ]

# Run analysis on samples from each time point
for(rr in 1:nrow(lookupSamples)){    cat("    sample    ");cat(rr);cat("   ");     #tau =0 # timepoint (0,14,180)
  tau <- lookupSamples[rr]$Day
  samp<- lookupSamples[rr]$orig.ident
  # Load V1:NN with cell ID and other metadata for every cell in patient i at timepoint tau
  phenotypes_i <- allphenotypes[Day == tau][orig.ident == samp]

  # Continue analysis only if some cells sampled at that timepoint
  if(nrow(phenotypes_i) > 0){
      # For a speciic sample day: Add gene expression (y1,y2) to phenotpye and discretized phenotype data 
      dd2 <-  merge(phenotypes_i, LRtransposed, by="Cell.ID")
      dd2[ ,samplesize_it:= nrow(dd2) ]   #dd2[,1:40]
      
      # Gather genes into long format
      dd3 <- data.table( gather( dd2, gene, expression, names( LRtransposed[, -1] ) ) )
      #ggplot(dd3[gene=="EGF"], aes(y=expression,x=Disc_V1,col=gene) ) + geom_point()+facet_wrap(~Celltype) #ggplot(dd3[gene=="EGFR"], aes(y=expression,x=Disc_V1,col=gene) ) + geom_point()+facet_wrap(~Celltype)
      
      ## Extract clinical data of the patient for merging
      clin_i <- dd3[1, ] %>% dplyr::select(Patient.Study.ID, Sample, orig.ident,Day, ARM:dynamic_class3)
      
      ### Summarise the average expression of genes in each cell discretization class      # dd4[,expression_bar:=median(expression),by=c("gene",paste0("Disc_V",1:NN))]
      grps <- c("gene", "samplesize_it", "Celltype", "PhenoCelltype", "key_", paste0("Disc_V", 1:5))  ##"Celltype_subtype",
      dd4 <- data.table( dd3 %>% dplyr::group_by_(.dots = grps) %>% dplyr::summarise(expression_bar:= mean(expression), countofvalues = n()))
      rm(list="dd3")
      LRGene.ID <- names(LRtransposed)[-1]
      ## Look ups of the combinations of gene and cell type combinations that are produced with differing labels of the columns for merging different datasets.
      lu1 <- rbindlist( lapply(1:length(LRGene.ID), function(gg){   data.table(gene= LRGene.ID[gg] , uu)   }) )
      lu2 <- lu1[ gene == lu1$gene[1] ] %>% dplyr::select(Celltype, PhenoCelltype, key_)
      setnames(lu2, old= c("key_", "Celltype", "PhenoCelltype"), new= c("Receptor", "ReceptorCelltype", "ReceptorPhenoCelltype"))
      lu2i <- lu1[ gene == lu1$gene[1] ] %>% dplyr::select(Celltype, PhenoCelltype, key_)
      setnames(lu2i, old= c("key_", "Celltype", "PhenoCelltype"), new= c("Ligand", "LigandCelltype", "LigandPhenoCelltype"))
      
      
      ## Merge clinical data with average expression data and add back in any cell types with zero abundance if needed (usally not)
      addzeros <- FALSE
      if( addzeros== TRUE ){
        dd5 <- data.table(clin_i, merge(lu1, dd4, all= TRUE, fill= 0)  )
        dd5[is.na(expression_bar), expression_bar := 0]
        dd5[is.na(countofvalues), countofvalues := 0]
        dd5[is.na(samplesize_it), samplesize_it := nrow(dd2)]
        dd5[, FracSample := countofvalues/samplesize_it ]
      }else{
        dd5 <- data.table(clin_i, dd4  )
        dd5[,FracSample:=countofvalues/samplesize_it ]
      }
      #dd5[,Treatment:=NULL ]
      dd5[,orig.ident:=NULL ]
      #dd5[,dynamic_class:=NULL ]
      dd5[,Sample:=NULL ]
      
      rm(list=c("dd2","dd4","lu1"))
      cat("    calculating communication    ")
      ###Quantify communication: Ligand-Receptor pair (p) interaction   #p=352
      Communication <- rbindlist(mclapply( 1:nrow(LRpairsFiltered2_i)   , function(p){   #1:nrow(LRpairsFiltered2_i)
        # Names of ligand and receptor
        LRnms <- unname( unlist( LRpairsFiltered2_i[p][ , HPMR.Receptor, HPMR.Ligand ] ) )
        # Select expression of the pair in each cell subtype
        LRpair_it_dd <-  data.table( dd5[gene %in% LRnms ]%>%spread(gene, expression_bar) ,  LRpairsFiltered2_i[p] %>% dplyr::select(HPMR.Receptor, HPMR.Ligand, Pair.Name))
        setnames(LRpair_it_dd ,old= LRnms, new= c("Ligand", "Receptor"))
        setcolorder(LRpair_it_dd,c( names(dd5[1]%>%dplyr::select(-c("gene","expression_bar"))), "Ligand", "Receptor", "HPMR.Receptor", "HPMR.Ligand", "Pair.Name")  )
        
        # Calculate the expression of the signaller cell type by multiplying single cell average expression by the cell number
        LRpair_it_dd[ , Ligand_N:= Ligand*FracSample]
        
        # Caclulate ligand-receptor signalling between each cell type: outer product matrix -> Signaler on the cols and receiver cell class on the rows
        Rmat <- LRpair_it_dd[, Receptor] %*% t( unlist( LRpair_it_dd[, "Ligand_N"] ) )
        rownames(Rmat) <- colnames(Rmat) <- LRpair_it_dd$key_
        RmatperSignaller <- LRpair_it_dd[, Receptor] %*% t( unlist( LRpair_it_dd[, "Ligand"] ) )
        rownames(RmatperSignaller) <- colnames(RmatperSignaller) <- LRpair_it_dd$key_
        
        # Marginalise signalling matrix to calculate signal transduction to each receiver cell type
        LRpair_it_dd[ , Transduction := rowSums(Rmat) ]
        LRpair_it_dd[ , TransductionperSignaller := rowSums(RmatperSignaller) ]
        
        # Refromat the ligand-receptor signalling matrix into a long dataframe
        Rmatlong <- data.table( gather(as.data.table(Rmat, keep.rownames = T), Ligand, Signal, -1) )[Signal>0]
        setnames(Rmatlong, old= "rn", new= "Receptor")
        
        # Merge ligand-receptor signalling with cell type information
        Rmatlong2 <- merge( merge(Rmatlong, lu2, by= "Receptor"), lu2i , by= "Ligand")
        # Merge ligand-receptor signalling with transduction data and all clinical information
        setnames(Rmatlong2, old= c("Ligand", "Receptor"),new= c("key_signaller", "key_"))
        output <- data.table( merge( LRpair_it_dd, Rmatlong2, by= "key_", all.x=T)  )
        output[,Celltype:=NULL]
        output[,PhenoCelltype:=NULL]
        return( output )
        #cat(p);cat("    ")
      },mc.cores= detectCores()-2))
      
      PatientID <- Communication[1]$Patient.Study.ID
      cat("Saving output for :        patient "); cat(PatientID); cat("       time     ") ; cat(tau)
      savenm <- paste0("CommunicationResults__","PatientID_",PatientID,"__", "Day_",tau,".RData")
      #saveloc <- "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/CommunicationOutputAllArms/"
      save( PatientID, Communication, tau, rr, uu, addzeros, file=paste0(saveloc,savenm))
      rm(list=c("dd5","lu2","lu2i","Communication"))
      
    }
  cat("Sample done   ");cat(100*rr/nrow(lookupSamples)); cat("%"); cat("         ")
    
  }
  

# NB) Rmat[1,] == RmatperSignaller[1,] * LRpair_it_dd$FracSample  and so RmatperSignaller[1,]  = Rmat[1,]/LRpair_it_dd$FracSample == 
list.files(saveloc)%>%length
list.files(saveloc)
sapply(1:length(list.files(saveloc)),function(x){
  strsplit(list.files(saveloc)[x],"__")[[1]][2]
  
})%>%unique()
allphenotypes[]$Patient.Study.ID%>%unique()%>%length()
