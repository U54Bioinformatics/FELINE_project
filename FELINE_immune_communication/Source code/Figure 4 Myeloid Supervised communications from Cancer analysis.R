rm(list=ls())   
library(abind)
require(data.table); require(dplyr); require(ggplot2); require(tidyr)
require(igraph)


### Load Ligand Receptor database list of Ramilowski et al 2015
LRpairsFiltered <- data.table(read.csv( file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/SI data/FigureS4/Figure S4 LigandReceptorPairsRamilowski2015.csv"))
LRgenelist <- unique( c(LRpairsFiltered$HPMR.Receptor, LRpairsFiltered$HPMR.Ligand) )

### Load some signalling data for downstream analysis
savloc<-"/Users/jason/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/Figure3/"
CCI_all <-  read.csv(file=paste0(savloc, "SourceData_Figure3_CellCummunicationIndexTumorWideSingnaling.csv"))
CCI_discovery <- data.table(CCI_all)[Cohort=="Discovery"]
CCI_validation <-  data.table(CCI_all)[Cohort=="Validation"]

# analysis for Discovery
CCI <- data.table(CCI_discovery)

# scaling
CCI[,scaleTransduction:=scale(TransductionMu,center=F),by=c("Pair.Name")] #"ReceptorPhenoCelltype",
CCI[!is.finite(TransductionMu),scaleTransduction:=0]
CCI[,scalelnTransductionMu:=scale(log(1+TransductionMu)),by=Pair.Name]

# Select cancer - macrophage communications for communication pathways linked to M2 differentiation.
#corruptsigs0 <- CCI[Day!=180][Treat=="CombinationRibo"][LigandPhenoCelltype=="Cancer cells"][ReceptorPhenoCelltype=="Macrophages"]
corruptsigs0 <- CCI[Day==0][Day!=180][Treat=="CombinationRibo"][LigandPhenoCelltype=="Cancer cells"][ReceptorPhenoCelltype=="Macrophages"]

# create supervised list LR communications
nsm1<- unique(corruptsigs0[Day!=180]$Pair.Name)
extractThese0 <- c(grep("FGF2_",nsm1),#grep("FGF1_",nsm1),
                   grep("IL11_",nsm1),grep("IL12A_",nsm1),
                   grep("CLCF1_",nsm1),grep("MADCAM1_",nsm1),
                   #grep("SEMA3",nsm1),grep("SEMA4",nsm1),
                   grep("PLAU_",nsm1),grep("GDF9_",nsm1))
extractThese1 <- c(grep("_IL4",nsm1),grep("_IL13",nsm1),grep("_IL10",nsm1)         )
extractThese2 <- c(grep("_CSF1R",nsm1),grep("_CSF2RB",nsm1))#,grep("_MERTK",nsm1),grep("GAS_",nsm1) )
extractThese3 <- c(grep("_FGFR1",nsm1),grep("IL6_",nsm1) )
selectList<- nsm1[c(extractThese0,extractThese1,extractThese2,extractThese3)]
selectList<-selectList[selectList!="NCAM1_FGFR1"]
Supervisedcorruptsigs <- corruptsigs0[Pair.Name%in%selectList] 

# Restructure the data
corruptsigstmp <- data.table(Supervisedcorruptsigs%>%dplyr::select(Pair.Name, scalelnTransductionMu, Patient.Study.ID, dynamic_class3, Day)%>% spread(Pair.Name,scalelnTransductionMu,fill=0))[order(dynamic_class3,Patient.Study.ID,Day)]
clustmatResponder<- data.frame( corruptsigstmp[][dynamic_class3=="Response"]%>% dplyr::select(-c(Patient.Study.ID,dynamic_class3,Day))  )
clustmatNonResponder<- data.frame( corruptsigstmp[][dynamic_class3=="Non-response"]%>% dplyr::select(-c(Patient.Study.ID,dynamic_class3,Day))  )
# Reorder samples to place resistant and sensitive tumor samples together
TumorOrder <- c(
  as.character(unique(corruptsigstmp[][dynamic_class3=="Non-response"]$Patient.Study.ID[hclust(dist(clustmatNonResponder))$order]) ),
               as.character(unique(corruptsigstmp[][dynamic_class3=="Response"]$Patient.Study.ID   [hclust(dist(clustmatResponder))$order])))

Supervisedcorruptsigs$Patient.Study.ID <- factor(Supervisedcorruptsigs$Patient.Study.ID , levels=TumorOrder)

corruptsigs <- data.table(Supervisedcorruptsigs%>%dplyr::select(Pair.Name, scalelnTransductionMu, Patient.Study.ID, dynamic_class3, Day)%>% spread(Pair.Name,scalelnTransductionMu,fill=NA))[order(dynamic_class3,Patient.Study.ID,Day)]

# Plot as heatmap
pp <- data.frame( corruptsigs%>% dplyr::select(-c(Patient.Study.ID,dynamic_class3,Day))  )
pp <- pp[,colMeans(!is.na(pp))>0.3]

rownames(pp)<- paste0("samp",1:nrow(pp))
rownam <- data.frame( dynamic_class3= corruptsigs$dynamic_class3 )
rownames(rownam) <- paste0("samp",1:nrow(pp))
wellMeasured <- rowMeans(!is.na(pp))>0.33

## annotations, coloration and labels for plotting
pp<- pp[wellMeasured,]
rownam <- rownam[wellMeasured,]
annotation_row <- data.frame(corruptsigs[Day!=180][wellMeasured,]%>%dplyr::select(dynamic_class3, Day,Patient.Study.ID)) #orig.ident
rownames(annotation_row) <- paste0(annotation_row$Patient.Study.ID, annotation_row$Day);annotation_row$Patient.Study.ID <-NULL;annotation_row$Day <-NULL
names(annotation_row) <- "Tumor response"

annotation_col = list( "Tumor response" = rev(pal_npg("nrc")(2) ))
names(annotation_col[[1]] ) <- c("Response","Non-response")
rownames(pp) <-rownames(annotation_row)
brk <- sum(annotation_row$'Tumor response'=="Non-response")

fillseq <- seq(-2.35,4.45,length=101)
pheatmap::pheatmap((pp), gaps_row=brk, cluster_rows=F, cluster_cols=F, breaks=fillseq, na_col="white", annotation_row = annotation_row, annotation_colors = annotation_col, scale="none",show_rownames=F,annotation_legend=F,border_color=NA)

# Gather results
corruptsigs [,rn:=paste0(Patient.Study.ID, Day)]
discCommsCombinationRibo<- data.table(corruptsigs[rn%in%rownames(pp)])
discCommsCombinationRibo[,rn:=NULL]
discCommsCombinationRibo[,Treat:="CombinationRibo"]
discCommsCombinationRibo[,Cohort:="Discovery"]

## repeat for letrozole arm of discovery cohort
# Select cancer - macrophage communications for communication pathways linked to M2 differentiation.
corruptsigs<-CCI[Day==0][Day!=180][Treat!="CombinationRibo"][LigandPhenoCelltype=="Cancer cells"][ReceptorPhenoCelltype=="Macrophages"]
nsm1<- unique(corruptsigs[Day!=180]$Pair.Name)
supervisedlist<-c("CLCF1_CRLF1","CLCF1_IL6ST","CLCF1_LIFR", "CSF1_CSF1R", "FGF1_FGFR1","FGF17_FGFR1","FGF18_FGFR1","FGF2_CD44","FGF2_FGFR1","FGF2_NRP1","FGF2_SDC1","FGF2_SDC2","FGF2_SDC3","FGF2_SDC4","FGF23_FGFR1","FGF7_FGFR1",
                  "GDF9_ACVR2A", "GDF9_BMPR1A" ,"GDF9_BMPR1B", "GDF9_BMPR2" , "GDF9_TGFBR1","IL11_IL6ST","IL11_IL11RA","IL12A_IL12RB1","IL12A_IL12RB2","IL5_CSF2RB","MADCAM1_CD44", "MADCAM1_ITGA4",
                  "PLAU_IGF2R","PLAU_LRP1" ,"PLAU_LRP2" ,"PLAU_PLAUR")
Supervisedcorruptsigs <- corruptsigs[Pair.Name%in%supervisedlist] 
corruptsigs <- data.table(Supervisedcorruptsigs%>%dplyr::select(Pair.Name, scalelnTransductionMu, Patient.Study.ID,dynamic_class3,Day)%>% spread(Pair.Name,scalelnTransductionMu,fill=NA))[order(dynamic_class3,Patient.Study.ID,Day)]

# Plot as heatmap
pp <- data.frame( corruptsigs%>% dplyr::select(-c(Patient.Study.ID,dynamic_class3,Day))  )
#pp<- pp[,colMeans(!is.na(pp))>0.3]
rownames(pp)<- paste0("samp",1:nrow(pp))
rownam <- data.frame( dynamic_class3= corruptsigs$dynamic_class3 )
rownames(rownam)<- paste0("samp",1:nrow(pp))

annotation_row <-data.frame(corruptsigs[Day!=180]%>%dplyr::select(dynamic_class3, Day,Patient.Study.ID)) #orig.ident
rownames(annotation_row) <- paste0(annotation_row$Patient.Study.ID, annotation_row$Day);annotation_row$Patient.Study.ID <-NULL;annotation_row$Day <-NULL
names(annotation_row) <- "Tumor response"

annotation_col = list( "Tumor response" = rev(pal_npg("nrc")(2) ))
names(annotation_col[[1]] ) <- c("Response","Non-response")
rownames(pp) <-rownames(annotation_row)
brk <- sum(annotation_row$'Tumor response'=="Non-response")
pheatmap::pheatmap((pp),gaps_row=brk,breaks= fillseq,cluster_rows=F,cluster_cols=F,na_col="white",annotation_row = annotation_row, annotation_colors = annotation_col,scale="none",show_rownames=F,annotation_legend=F,border_color=NA)

blankrowannot <- annotation_row
names( blankrowannot)<- c(" ")

names( annotation_col)<- names( blankrowannot)
names(annotation_col[[1]] ) <- c("Response","Non-response")

p1 <- pheatmap::pheatmap((pp), gaps_row=brk, cluster_rows=F, cluster_cols=F, breaks=fillseq, na_col="white", 
                         annotation_row = blankrowannot, annotation_colors = annotation_col, 
                         scale= "none", show_rownames=F, annotation_legend=F, border_color=NA,labels_col=rep("",100),
                         labels_row= rep("",100) ,
                         legend_breaks = -4:2,
                         legend_labels = rep("",7) )

corruptsigs [,rn:=paste0(Patient.Study.ID, Day)]
discCommsLetro<- data.table(corruptsigs[rn%in%rownames(pp)])
discCommsLetro[,rn:=NULL]
discCommsLetro[,Treat:="LetrozoleAlone"]
discCommsLetro[,Cohort:="Discovery"]








##### Matching analysis for Validation cohort
CCI <- data.table(CCI_validation)
corruptsigs<-CCI[Day==0][Day!=180][Treat=="CombinationRibo"][LigandPhenoCelltype=="Cancer cells"][ReceptorPhenoCelltype=="Macrophages"]
nsm1<- unique(corruptsigs[Day!=180]$Pair.Name)
extractThese <- c(grep("CSF1_CSF1R",nsm1) ,grep("FGF2_",nsm1),grep("FGF1_",nsm1),grep("IL11_",nsm1),grep("IL12A_",nsm1),
                  grep("CLCF1_",nsm1),grep("MADCAM1_",nsm1),grep("SEMA3",nsm1),grep("SEMA4",nsm1),grep("PLAU_",nsm1),grep("GDF9_",nsm1))
supervisedlist<-c("CLCF1_CRLF1","CLCF1_IL6ST","CLCF1_LIFR", "CSF1_CSF1R", "FGF1_FGFR1","FGF17_FGFR1","FGF18_FGFR1","FGF2_CD44","FGF2_FGFR1","FGF2_NRP1","FGF2_SDC1","FGF2_SDC2","FGF2_SDC3","FGF2_SDC4","FGF23_FGFR1","FGF7_FGFR1",
                   "GDF9_ACVR2A", "GDF9_BMPR1A" ,"GDF9_BMPR1B", "GDF9_BMPR2" , "GDF9_TGFBR1","IL11_IL6ST","IL11_IL11RA","IL12A_IL12RB1","IL12A_IL12RB2","IL5_CSF2RB","MADCAM1_CD44", "MADCAM1_ITGA4",
                   "PLAU_IGF2R","PLAU_LRP1" ,"PLAU_LRP2" ,"PLAU_PLAUR")
extractThese <- which(nsm1 %in% supervisedlist )
#nsm1[grep("FGF17",nsm1)]

Supervisedcorruptsigs <- corruptsigs[Pair.Name%in%nsm1[extractThese]] 
corruptsigs <- data.table(Supervisedcorruptsigs%>%dplyr::select(Pair.Name, scalelnTransductionMu, Patient.Study.ID,dynamic_class3,Day)%>% spread(Pair.Name,scalelnTransductionMu,fill=NA))[order(dynamic_class3,Day,Patient.Study.ID)]

# Plot as heatmap
pp <- data.frame( corruptsigs%>% dplyr::select(-c(Patient.Study.ID,dynamic_class3,Day))  )
pp<- pp[,colMeans(!is.na(pp))>0.3]
rownames(pp)<- paste0("samp",1:nrow(pp))
rownam <- data.frame( dynamic_class3= corruptsigs$dynamic_class3 )
rownames(rownam) <- paste0("samp",1:nrow(pp))

annotation_row <-data.frame(corruptsigs[Day!=180]%>%dplyr::select(dynamic_class3, Day,Patient.Study.ID)) #orig.ident
rownames(annotation_row) <- paste0(annotation_row$Patient.Study.ID, annotation_row$Day);annotation_row$Patient.Study.ID <-NULL;annotation_row$Day <-NULL
names(annotation_row)[1] <- "Tumor response"

annotation_col = list( "Tumor response" = rev(pal_npg("nrc")(2) ))
names(annotation_col[[1]] ) <- c("Response","Non-response")
rownames(pp) <-rownames(annotation_row)

brk <- sum(annotation_row$'Tumor response'=="Non-response")
colnames(pp) <- gsub("_","-",colnames(pp))


pheatmap::pheatmap((pp),gaps_row=brk,cluster_rows=F,cluster_cols=F,na_col="white",annotation_row = annotation_row, annotation_colors = annotation_col,scale="none",show_rownames=F,annotation_legend=F,border_color=NA)

corruptsigs [,rn:=paste0(Patient.Study.ID, Day)]
validCommsCombinationRibo<- data.table(corruptsigs[rn%in%rownames(pp)])
validCommsCombinationRibo[,rn:=NULL]
validCommsCombinationRibo[,Treat:="CombinationRibo"]
validCommsCombinationRibo[,Cohort:="Validation"]



corruptsigs<-CCI[Day==0][Day!=180][Treat!="CombinationRibo"][LigandPhenoCelltype=="Cancer cells"][ReceptorPhenoCelltype=="Macrophages"]
nsm1<- unique(corruptsigs[Day!=180]$Pair.Name)
extractThese <- c(grep("FGF2_",nsm1),grep("FGF1_",nsm1),grep("IL11_",nsm1),grep("IL12A_",nsm1),
                  grep("CLCF1_",nsm1),grep("MADCAM1_",nsm1),grep("SEMA3",nsm1),grep("SEMA4",nsm1),grep("PLAU_",nsm1),grep("GDF9_",nsm1))
Supervisedcorruptsigs <- corruptsigs[Pair.Name%in%supervisedlist] 
corruptsigs <- data.table(Supervisedcorruptsigs%>%select(Pair.Name, scalelnTransductionMu, Patient.Study.ID,dynamic_class3,Day)%>% spread(Pair.Name,scalelnTransductionMu,fill=NA))[order(dynamic_class3,Patient.Study.ID,Day)]

# Plot as heatmap
pp <- data.frame( corruptsigs%>% select(-c(Patient.Study.ID,dynamic_class3,Day))  )
rownames(pp)<- paste0("samp",1:nrow(pp))
rownam <- data.frame( dynamic_class3= corruptsigs$dynamic_class3 )
rownames(rownam)<- paste0("samp",1:nrow(pp))
annotation_row <-data.frame(corruptsigs[Day!=180]%>%dplyr::select(dynamic_class3, Day,Patient.Study.ID)) #orig.ident
rownames(annotation_row) <- paste0(annotation_row$Patient.Study.ID, annotation_row$Day);annotation_row$Patient.Study.ID <-NULL;annotation_row$Day <-NULL
names(annotation_row) <- "Tumor response"
annotation_col = list( "Tumor response" = rev(pal_npg("nrc")(2) ))
names(annotation_col[[1]] ) <- c("Response","Non-response")
rownames(pp) <-rownames(annotation_row)
brk <- sum(annotation_row$'Tumor response'=="Non-response")
colnames(pp) <- gsub("_","-",colnames(pp))
pheatmap::pheatmap((pp),gaps_row=brk,cluster_rows=F,cluster_cols=F,na_col="white",annotation_row = annotation_row, annotation_colors = annotation_col,scale="none",show_rownames=F,annotation_legend=F,border_color=NA)

corruptsigs [,rn:=paste0(Patient.Study.ID, Day)]
validCommsLetro<- data.table(corruptsigs[rn%in%rownames(pp)])
validCommsLetro[,rn:=NULL]
validCommsLetro[,Treat:="LetrozoleAlone"]
validCommsLetro[,Cohort:="Validation"]




# Gather all results
resultstable<-data.table( rbind(
  discCommsCombinationRibo%>%select(names(validCommsCombinationRibo)),
  discCommsLetro  ,
  validCommsCombinationRibo,
  validCommsLetro) )

savloc<-"/Users/jason/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/Figure4/"
write.csv(resultstable, file= paste0(savloc, "SourceData_Figure4_CancertoMyeloidSupervisedLRComms_Out.csv"))
