rm(list=ls())
library(umap)
require(data.table); require(dplyr); require(ggplot2); require(tidyr)
library(car)
require(org.Hs.eg.db)
require(GEOquery )
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE75130", "file=GSE75130_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)

# pre-filter low count genes
# keep genes with at least 2 counts > 10
keep <- rowSums( tbl >= 10 ) >= 2
tbl <- tbl[keep, ]

# log transform raw counts
# instead of raw counts can display vst(as.matrix(tbl)) i.e. variance stabilized counts
dat <- log10(tbl + 1)
dat <- dat[!duplicated(dat), ]
rownames(dat) <-unname(mapIds(
  org.Hs.eg.db,
  keys = rownames(dat),
  column = 'SYMBOL',
  keytype = 'ENTREZID'))
dat <-dat[!is.na(rownames(dat)),] # first remove duplicates

# box-and-whisker plot
par(mar=c(7,4,2,1))
boxplot(dat, boxwex=0.7, notch=T, main="GSE75130", ylab="lg(cnt + 1)", outline=F, las=2)

# meta data
meta <- getGEO("GSE75130", GSEMatrix =F)
meta@gsms$GSM1943688
metadd <-rbindlist(lapply(1:length(meta@gsms),function(ii){
  ch1<-meta@gsms[[ii]]@header$characteristics_ch1
  y<-data.table( sample=names(meta@gsms[ii]),
                 cellLine=  gsub("cell line: ","",ch1[grep("cell line: ",ch1)] ) ,
                 cellType=  gsub("cell type: ","",ch1[grep("cell type: ",ch1)] ) ,
                 cultureCondition= gsub("culture condition: ","",ch1[grep("culture condition: ",ch1)] ) ,
                 immunopathClass= gsub("immunopathological classification: ","",ch1[grep("immunopathological classification: ",ch1)] ))
  return(y)
}))

# UMAP plot (dimensionality reduction)
ump <- umap(t(dat), n_neighbors = 7, random_state = 123)
plot(ump$layout, main="GSE75130 UMAP plot, nbrs =7", xlab="", ylab="", pch=20, cex=1.5)
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)


# DE analysis
dat_t<-data.table(t(dat),keep.rownames = T)
setnames(dat_t, old="rn",new="sample")
genecodes <- names(dat_t)[-1]
dat_de <- merge(metadd[cellType=="Monocyte"], dat_t,by="sample")
dat_de$cultureCondition<-factor(dat_de$cultureCondition,levels=c("single","co-culture"))

deres <- rbindlist(lapply(1:length(genecodes),function(ii){
  cat(ii)
  tryCatch({
    data.table(gene=genecodes[ii],
               coef(summary(lm(paste0(genecodes[ii],"~cultureCondition"),data=dat_de))),
               keep.rownames = T)[rn!="(Intercept)"]
  }  , error = function(cond) {return(NULL)} )
}))
setnames(deres, old=c("Pr(>|t|)","t value","Std. Error"),new=c("pvalue","tvalue","Std.Error"))

deres[,FDR:=p.adjust(pvalue,method="fdr")]
signif <- deres[FDR<0.05][order(-Estimate)]
signif[1:25]
signif[abs(Estimate)>1.5][order(-abs(Estimate))][1:50][order(-Estimate)]
ggplot(signif[abs(Estimate)>1.5], aes(y=FDR, Estimate)) + 
  geom_point()
deres[gene%in%"MSR1"]
deres[gene%in%"CD163"]
deres[gene%in%"CSTB"]
deres[gene%in%"MRC1"]

markersdd<-data.table(read.csv("/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/M1M2like gene list/M1M2like gene list.csv"))
M2MarkerGenes<- c("ARG1","ARG2","IL10","CD163","CD23","CD200R1","PDCD1LG2" ,"CD274",
                "MARCO" ,   "CSF1R" ,   "MRC1" ,    "Il1RA" ,   "Il1R2" ,   "IL4R"  ,   "CCL4"  ,   "CCL13" ,  
                "CCL20" ,   "CCL17" ,   "CCL18",    "CCL22" ,   "CCL24" ,   "LYVE1" ,   "VEGFA" ,   "VEGFB" ,  
                "VEGFC"  ,  "VEGFD" ,   "EGF" ,     "CTSA"  ,   "CTSB" ,    "CSTC"  ,   "CTSD"   ,  "TGFB1" ,  
                "TGFB2"  ,  "TGFB3" ,   "MMP14",    "MMP19",    "MMP9" ,    "CLEC7A" ,  "WNT7B" ,   "FASL"  ,  
                "TNFSF12",  "TNFSF8" ,  "CD276" ,   "VTCN1" ,   "MSR1"  ,   "FN1"   ,   "IRF4"  ,   "CD36" ,   
                "PPARG" ,   "LIPA" ,    "CYP27A",   "DHRS9" ,   "FABP4"  ,  "SPP1"  ,   "FAM20C" ,  "LRP1" ,   
                "NRP1" ,    "MITF" ,    "SPOCD1" ,  "FABP5" ,   "DOCK3" ,   "LPL")
  
ggplot(deres[gene%in% M2MarkerGenes ],
       aes(y=FDR, Estimate,col=FDR<0.05)) + 
  geom_point()

signifmarkers<-deres[gene%in% M2MarkerGenes ][order(-Estimate)]
signifmarkers[,FDR_mark:=p.adjust(pvalue,method="fdr")]
signifmarkers[FDR_mark<0.05]
signifmarkers$gene<- factor(signifmarkers$gene, levels=rev(signifmarkers$gene))
ggplot(signifmarkers[FDR_mark<0.05]  ,
  aes(y=gene, x=rn, fill=Estimate))+ theme_classic(base_size=26)+
  geom_tile()+
  scale_fill_gradient2(name="M2-like macrophage \n marker gene \n fold change in \n cancer co-culture \n (vs monoculture)",low="darkblue",high="darkred", mid="white",midpoint=0)+
  labs(y="Gene", x="")+
  scale_x_discrete(labels="")+
  theme(aspect.ratio=1.5, axis.ticks.x = element_blank())
paperfile<- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"In vitro and Validation Myeloid coculture M2 Markers.png"),height=10,width=10, dpi=320)

resout<-signifmarkers[FDR_mark<0.05]
savloc <- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/SI data/FigureS14/"
write.csv(resout, file=paste0(savloc,"SourceData_FigureS14_InVitroCancerStimulatesM2MyeloidDifferentiation.csv"))
