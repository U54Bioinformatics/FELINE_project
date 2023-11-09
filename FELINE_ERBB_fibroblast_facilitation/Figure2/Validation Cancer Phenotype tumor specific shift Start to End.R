rm(list=ls())
require(data.table)
require(dplyr)
require(ggplot2)
require(tidyr)
require(ggsci)

### Phenotype data
load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/PhenotypesOfAllCellTypesAllArmsCohort2/PhenotypesOfAllCellTypesAllArmsCohort2.RData") #allphenotypes,UMAPlocs ,UMAPfiles,umapDImRedloc,umapDImRedfiles,nCellTypes,
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

### Load counts per million read data for a specific patient's samples and gen Ligand+Receptor genes
CPMlocs1 <- "/Users/jason/Dropbox/FELINE Project (1)/FELINE Cohort 2/DNANexusCopy/CPM/HQ/"
CPMlocs2 <- "/Users/jason/Dropbox/FELINE Project (1)/FELINE Cohort 2/DNANexusCopy/CPM/HQLQ/"

CPMfiles1 <- list.files(CPMlocs1)[ grep("CPM", list.files(CPMlocs1) ) ]
CPMfiles2 <- list.files(CPMlocs2)[ grep("CPM", list.files(CPMlocs2) ) ]
n10Xpats <- length(CPMfiles1)

discData<-data.table(read.csv(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/DiscoveryCohortDEGPrePost.csv"))
geneList<- as.character(discData$gene%>%unique())

# settings
shouldscale <- FALSE ## should cpm be scaled

x_celldd <- allphenotypes[PhenoCelltype=="Cancer cells"]
#saveloc <- "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/CommunicationOutputAllArms/"


cpm_i1 <- data.table( fread( paste0(CPMlocs1, CPMfiles1[grep("Cancer",CPMfiles1)]) ) )   # load full gene expression
names(cpm_i1)[1]<- "Gene.ID"
transposedHQ <- as.data.table( t(cpm_i1[,-1])  , keep.rownames = T)
colnames(transposedHQ) <- c("Cell.ID", cpm_i1$Gene.ID)
rm(list="cpm_i1")

cpm_i2<- data.table( fread( paste0(CPMlocs2, CPMfiles2[grep("Cancer",CPMfiles2)]) ) )   # load full gene expression
names(cpm_i2)[1]<- "Gene.ID"
transposedLQ <- as.data.table( t(cpm_i2[,-1])  , keep.rownames = T)
colnames(transposedLQ) <- c("Cell.ID", cpm_i2$Gene.ID)
rm(list="cpm_i2")

srtGene.ID0<-sort(intersect(names(transposedHQ)[-1],names(transposedLQ)[-1]))
srtGene.ID <- sort(intersect(srtGene.ID0,geneList))

transposedHQ<- transposedHQ[,c("Cell.ID",srtGene.ID),with=F]
transposedLQ<- transposedLQ[,c("Cell.ID",srtGene.ID),with=F]
transposedHQ[1:10,1:10]
transposedLQ[1:10,1:10]
DEGtransposed <- rbind(transposedHQ, transposedLQ)
rm(list="transposedHQ")
rm(list="transposedLQ")

x_celldd[,Treatment:="letrozole"]
x_celldd[ARM!="A",Treatment:="letrozole + ribo"]
DEGtransposedData<- merge(x_celldd,DEGtransposed, by="Cell.ID")

DEGtransposedData[,hasStart:= any(Day==0), by="Patient.Study.ID"]
DEGtransposedData[,hasEnd:= any(Day==180), by="Patient.Study.ID"]
DEGtransposedData[,hasStartEnd:=0] 
DEGtransposedData[(hasStart==T)&(hasEnd==T),hasStartEnd:=1]

start_enddd <-DEGtransposedData[Day!=14][hasStartEnd==1] 

lu1 <- unique( data.table(start_enddd%>% dplyr::select(Patient.Study.ID,Day, Treatment) ) ) #lu1[Day==0][Treatment=="letrozole"]
lu2 <- unique( data.table(start_enddd%>% dplyr::select(Patient.Study.ID, Treatment) ) )

expressedGenes <-srtGene.ID
#for each patient select cancer cell data (across timepoints but within one tumor) and perform DE analysis using all genes with >1% coverage
allDE <- rbindlist( lapply(1:nrow(lu2), function(i){
  cat("Tumor    ") ;  cat(i) ; cat("\n")
  
  dd_i <- start_enddd[ Patient.Study.ID==lu2[i]$Patient.Study.ID]   
  # Perform pathway analysis 
  pathanalysis0 <- rbindlist(lapply( 1:length(expressedGenes), function(pp){
    cat(pp)
    res <- tryCatch({
      m1 <- lm(formula= paste0("log(1+",expressedGenes[pp],") ~ as.factor(Day)") ,    dd_i )
      data.table(gene=expressedGenes[pp],  coef( summary(m1) ), keep.rownames = T)
    },
    error= function(x){
      return(NULL)
    })
    return(res)  
  } ))
  setnames(pathanalysis0, old= c( "Std. Error", "t value", "Pr(>|t|)" ), new= c( "Std.Error", "tval", "pval" ))
  pathanalysis0$Patient.Study.ID <- dd_i[1]$Patient.Study.ID
  pathanalysis0$ARM <- dd_i[1]$ARM
  pathanalysis0$Treatment <- dd_i[1]$Treatment
  pathanalysis0$dynamic_class <- dd_i[1]$dynamic_class
  pathanalysis0$dynamic_class3 <- dd_i[1]$dynamic_class3
  pathanalysis0$adj.pval <- p.adjust(p= pathanalysis0$pval, method= "fdr" )
  return(pathanalysis0)
}))
save(allDE,  file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Gene expression landscapes/cancer gene expression landscape/Validation Tumor specific shift strt vs end cancer gene expression scrna.RData")

pltthis1 <- allDE[rn!="(Intercept)"][  ] 
pltthis1$gene <- factor(pltthis1$gene , levels=rev( srtGene.ID) )

pltthis1[,dir:="up"]
pltthis1[gene%in%consistentDownreg$gene,dir:="down"]

pltthis1 <- merge(pltthis1, unique( discData%>%dplyr::select(gene,dir) ) , by="gene")
pltthis1

valData <- data.table( pltthis1%>%dplyr::select(-c(dynamic_class) ))
setnames( valData , old=c("rn"), new=c("ModelParamLM"))
valData[Treatment=="letrozole + ribo" ,Treatment:="letrozole + ribociclib"]


DEGdiscData <- data.table(read.csv(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/DiscoveryCohortDEGPrePost.csv"))

getordFun<-function(){
  load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Gene expression landscapes/cancer gene expression landscape/Tumor specific shift strt vs end cancer gene expression scrna.RData")
  load( file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Gene expression landscapes/cancer gene expression landscape/Tumor specific shift strt vs end cancer gene expression scrna summary.RData")
  topN <- 25
  consistentUpreg <- summary_allDE[nsignif> (length(unique(allDE[Treatment=="letrozole"]$Patient.Study.ID)))/2 ][mean_effect>log(1.5)][order(-abs(mean_effect))] [1:topN] 
  consistentUpreg[,dir:="up"]
  consistentDownreg <- summary_allDE[nsignif> (length(unique(allDE[Treatment=="letrozole"]$Patient.Study.ID)))/2 ][mean_effect<log(1.5)][order((mean_effect))] [topN:1] 
  consistentDownreg[,dir:="down"]
  consistentgenesUDwn <- na.omit(rbind(consistentUpreg,consistentDownreg))
  consistentgenesUDwn
}
consistentgenesUDwn <- getordFun()
consistentgenesUDwn$gene
DEGdiscData$gene <- factor(DEGdiscData$gene , levels=rev( consistentgenesUDwn$gene) )


ggplot( DEGdiscData, aes(y=Estimate,x=gene) ) + theme_classic(base_size=26) + geom_point()+
  #jitter(width= 0, height= 0.1) + 
  geom_hline(aes(yintercept= 0),linetype=2) +
  scale_fill_manual(name="Change under ET", labels=c("Downregulation","Upregulation"),values=c("yellow","red"))+
  labs(x= "Gene" , y= "Fold change in cancer cell expression") +facet_wrap(~Treatment)+geom_boxplot(aes(fill=dir,group=gene)) +coord_flip() +
  #scale_fill_aaas(name="Change under ET", labels=c("Downregulation","Upregulation")) +
  theme(aspect.ratio=2)

valData
valData$gene <- factor(valData$gene , levels=rev( consistentgenesUDwn$gene) )

ggplot( valData, aes(y=Estimate,x=gene) ) + theme_classic(base_size=26) + geom_point()+
  #jitter(width= 0, height= 0.1) + 
  geom_hline(aes(yintercept= 0),linetype=2) +
  scale_fill_manual(name="Change under ET", labels=c("Downregulation","Upregulation"),values=c("yellow","red"))+
  labs(x= "Gene" , y= "Fold change in cancer cell expression") +facet_wrap(~Treatment)+geom_boxplot(aes(fill=dir,group=gene)) +coord_flip() +
  #scale_fill_aaas(name="Change under ET", labels=c("Downregulation","Upregulation")) +
  theme(aspect.ratio=2)
ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Validation Single cell most significant changes in cancer cell expression in resistant cells YR.pdf", height=18, width = 18)

#write.csv(valData,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/ValidationCohortDEGPrePost.csv")

#DEGdiscData <- data.table(read.csv(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/Figure2/DiscoveryCohortDEGPrePost.csv"))
#valData <- data.table( read.csv(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/Figure2/ValidationCohortDEGPrePost.csv") )


d_out<-rbind( data.table(cohort="Discovery",DEGdiscData%>%dplyr::select(-X) ),
       data.table(cohort="Validation",valData %>%dplyr::select(-X)  ))

data.table(valData%>%group_by(gene)%>%dplyr::summarise(m=mean(Estimate),m.adj.pval=mean(adj.pval, nr.rm=T),dir=dir[1]))

valData[adj.pval<0.05]$gene%>%unique()
ggplot( d_out, aes(y=Estimate,x=gene,group=interaction(gene,cohort)) ) + theme_classic(base_size=26) + 
  geom_hline(aes(yintercept= 0),linetype=2) +
  scale_fill_manual(name="Change under ET", labels=c("Downregulation","Upregulation"),values=c("yellow","red"))+
  labs(x= "Gene" , y= "Fold change in cancer cell expression") +facet_wrap(~Treatment)+
  geom_boxplot(outlier.colour=NA,position=position_dodge(),aes(fill=dir,group=interaction(gene,cohort))) +
  stat_boxplot(geom = "errorbar",
               width = 0.5)+
  geom_point(aes(shape=cohort), size=1,position=position_dodge(width=1))+
  coord_flip() +
  theme(aspect.ratio=2)+
  scale_shape_discrete(name="Cohort")+
  guides(shape = guide_legend(override.aes = list(size = 2))) 
  
  
ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Discovery and Validation Single cell most significant changes in cancer cell expression in resistant cells YR.pdf", height=18, width = 18)

ggplot( d_out, aes(y=Estimate,x=gene,group=interaction(gene)) ) + theme_classic(base_size=26) + 
  geom_hline(aes(yintercept= 0),linetype=2) +
  scale_fill_manual(name="Change under ET", labels=c("Downregulation","Upregulation"),values=c("yellow","red"))+
  labs(x= "Gene" , y= "Fold change in cancer cell expression") +facet_wrap(~cohort)+
  geom_boxplot(outlier.colour=NA,position=position_dodge(),aes(fill=dir,group=interaction(gene,cohort))) +
  stat_boxplot(geom = "errorbar",
               width = 0.5)+
  geom_point(size=1,position=position_dodge(width=1))+
  coord_flip() +
  theme(aspect.ratio=2)+
  theme(legend.position = "none")

ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Discovery and Validation Single cell most significant changes in cancer cell expression in resistant cells YRCohort.pdf", height=18, width = 18)



test1 <-data.table( coef( summary(lm(Estimate~-1+gene,data=d_out[cohort=="Discovery"])) ) , keep.rownames = T)
test2 <-data.table( coef( summary(lm(Estimate~-1+gene,data=d_out[cohort=="Validation"])) ) , keep.rownames = T)
setnames(test1, old="Pr(>|t|)", new="pval")
setnames(test2, old="Pr(>|t|)", new="pval")

top15<- c( test1[pval<0.05][Estimate>0][order(-Estimate)][1:15,]$rn, test1[pval<0.05][Estimate<0][order(Estimate)]$rn)
top15 <- gsub("gene","",top15)
plotdd<-d_out[gene%in%top15]
plotdd$gene<-factor(plotdd$gene,levels=rev(top15))
plotdd<-plotdd[order(gene)]
ggplot( plotdd, aes(y=Estimate,x=gene,group=interaction(gene)) ) + theme_classic(base_size=26) + 
  geom_hline(aes(yintercept= 0),linetype=2) +
  scale_fill_manual(name="Change under ET", labels=c("Downregulation","Upregulation"),values=c("yellow","red"))+
  labs(x= "Gene" , y= "Fold change in cancer cell expression") +facet_wrap(~cohort)+
  geom_boxplot(outlier.colour=NA,position=position_dodge(),aes(fill=dir,group=interaction(gene,cohort))) +
  stat_boxplot(geom = "errorbar",
               width = 0.5)+
  geom_point(size=1,position=position_dodge(width=1))+
  coord_flip() +
  theme(aspect.ratio=2)+
  theme(legend.position = "none")

ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Discovery and Validation Single cell most significant changes in cancer cell expression in resistant cells YRCohort.pdf", height=10, width = 10)


discStatsdat <-data.table( coef( summary(lm(Estimate~-1+gene,data=d_out[gene%in%top15][cohort=="Discovery"])) ) , keep.rownames = T)
validStatsdat <-data.table( coef( summary(lm(Estimate~-1+gene,data=d_out[gene%in%top15][cohort=="Validation"]))  ) , keep.rownames = T)
d_out[cohort=="Discovery"][gene=="FOS"]
d_out[cohort=="Validation"][gene=="FOS"]
