rm(list=ls())
require(data.table)
require(dplyr)
require(ggplot2)
require(tidyr)
require(ggsci)
require(lmerTest)
require(slingshot)
library(mclust, quietly = TRUE)
#library(TSCAN)
#library(scater)
library(tradeSeq)
require(mgcv)


# myeloid UMAP
Sourceloc<-"/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/Figure4/"
SourceNam<-"SourceData_Figure4_MyeloidGeneUMAP.csv"
u_dat<- data.table(read.csv(file=paste0(Sourceloc,SourceNam)))
# first 15 rows contain metadata
gene_nms<-names(u_dat)[-c(1:15)]

# add label for treatment group
u_dat[,TreatmentLab:="Combination ribociclib"]
u_dat[ARM=="A",TreatmentLab:="Letrozole alone"]


# gene expression data per cell
u_datx <- as.matrix(u_dat%>%dplyr::select(any_of(gene_nms)) )

# Rotate Dimension 1 by 180 degrees to align with differentiation in dimension 2 of UMAP space.
u_dat$UMAP1<- -u_dat$UMAP1

sz <- 1.25


ggplot( u_dat[ARM!="A"] , aes(UMAP1,UMAP2,col=  dynamic_class3)) +  
  geom_point(size=sz)+#scale_color_ucscgb(name="Myleoid cell types")+
  theme_classic()+facet_wrap(TreatmentLab~Cohort)+
  theme(aspect.ratio=1,axis.text=element_blank(),axis.title=element_blank(),legend.position="none")
ggplot( u_dat[ARM!="A"] , aes(UMAP1,UMAP2,col=  dynamic_class3)) +  
  geom_point(size=sz)+#scale_color_ucscgb(name="Myleoid cell types")+
  theme_classic()+facet_wrap(dynamic_class3~Cohort)+
  theme(aspect.ratio=1,axis.text=element_blank(),axis.title=element_blank(),legend.position="none")

ggplot( u_dat[ARM!="A"] , aes(dynamic_class3,UMAP2,col=  dynamic_class3)) +  
  geom_point(size=sz)+geom_boxplot()+#scale_color_ucscgb(name="Myleoid cell types")+
  theme_classic()+facet_wrap(TreatmentLab~Cohort)+
  theme(aspect.ratio=1,axis.text=element_blank(),axis.title=element_blank(),legend.position="none")

# umap coordinates
ucoords<- u_dat %>% dplyr::select(UMAP1,UMAP2)

# GMM clusteing using BIC model selection
set.seed(123)
clmod1 <- Mclust(ucoords, G=5:21, modelNames = "EVV")#max(clmod1$classification)
cl1 <- clmod1$classification
u_dat$classification <-cl1



#u_datCombo$classification <-cl1
p<-ggplot( u_dat , aes(UMAP1,UMAP2)) +  
  geom_point(aes(col=  as.factor(classification)),size=sz)+#scale_color_ucscgb(name="Myleoid cell types")+
  theme_classic()+#facet_wrap(~classification)+
  theme(aspect.ratio=1,axis.text=element_blank(),axis.title=element_blank(),legend.position="none")
p
trajstart<-names(u_dat[Celltype_subtype=="Monocytes"]$classification%>%table()%>%which.max())
#ggplot( u_dat , aes(UMAP1,UMAP2)) + geom_point(aes(col=  (classification)==trajstart),size=sz)#scale_color_ucscgb(name="Myleoid cell types")+
sce <- slingshot(as.matrix(ucoords), clusterLabels = cl1, reducedDim = NULL,start.clus =trajstart)
u_dat$common.pseudo <- averagePseudotime(sce,i=1)  #u_dat$common.pseudo2 <- averagePseudotime(sce,i=2) 
sce_dd <- as.SlingshotDataSet(sce)
curves <- slingCurves(sce, as.df = TRUE)

p + geom_path(data = curves %>% arrange(Order),
              aes(UMAP1,UMAP2,group = Lineage)) 

ggplot(u_dat,aes(UMAP1,UMAP2))+geom_point(aes(col = common.pseudo))+
  geom_path(data = data.table(curves) %>% arrange(Order),
            aes(UMAP1,UMAP2,group = Lineage)) +
  scale_color_viridis_c(option = "B")+
  theme(aspect.ratio=1)

hist(log(u_dat$common.pseudo), breaks=100)
clustB<-Mclust( log(u_dat$common.pseudo),G=4, modelNames="E")#),                G=2:4)
clustB<-Mclust( log(u_dat$common.pseudo),G=2, modelNames="E")#),                G=2:4)
u_dat$common.pseudoclass<- clustB$classification
u_dat$low<-u_dat$common.pseudo<2
undifferentiated <- u_dat[low==T]
landlm <- gam(UMAP2 ~s(UMAP1, k=3) , data=undifferentiated)
undifferentiated$pred<-predict(landlm)

u_dat$preddivide<-predict(landlm,newdata=u_dat)
u_dat [ , isinflam:=F]
u_dat [UMAP2>preddivide , isinflam:=T]
u_dat[, polarization:= common.pseudo *( (2*(UMAP2-preddivide))-1) ]

SourceDataMyeloidPhenotype<-u_dat%>%dplyr::select(c( names(u_dat)[1:15],polarization))
#savloc4<-"/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/Figure4/"
#write.csv(SourceDataMyeloidPhenotype, file=paste0(savloc4,"SourceData_Figure4_Myeloid UMAP Phenotype Polarization_Out.csv"))

SourceDataMyeloidPhenotypeSI<-u_dat%>%dplyr::select(c( names(u_dat)[1:15],"classification","common.pseudo","polarization", #CD36, LIPA, PPARG,CYP27A1, DHRS9,
                                                       "CD36","PPARG","DHRS9","PLA2G7","DOCK3","FABP4","FABP5","LIPA","CYP27A1", "MITF", "KCNE1",
                                                       "CSF2RA", "CIITA","C3","ERBB4","ENOX1","ZEB1","FLT3","WNT5B"
                                                       ))

#savlocS8<-"/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/SI data/FigureS8/"
#write.csv(SourceDataMyeloidPhenotypeSI, file=paste0(savlocS8,"SourceData_Figure4_Myeloid UMAP Phenotype Polarization_Out.csv"))


SourceDataMyeloidPhenotypeSI10<-u_dat%>%dplyr::select(c( names(u_dat)[1:15],
                                                       "CD36","DHRS9","PLA2G7","DOCK3","FABP4","FABP5","LIPA","CYP27A1", "MITF", "KCNE1",
                                                       "CSF2RA", "CIITA","C3","ERBB4","ENOX1","ZEB1","FLT3","WNT5B"
))

#savlocS10<-"/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/SI data/FigureS10/"
#write.csv(SourceDataMyeloidPhenotypeSI10, file=paste0(savlocS10,"SourceData_FigureS10_Myeloid Marker Overlay_Out.csv"))



ggplot(u_dat,aes(UMAP1,UMAP2))+geom_point(aes(col = common.pseudo))+
  geom_path(data = data.table(curves) %>% arrange(Order),
            aes(UMAP1,UMAP2,group = Lineage)) +
  scale_color_viridis_c(option = "B")+
  theme(aspect.ratio=1)

p1<-ggplot(u_dat,aes(UMAP1,UMAP2))+geom_point(aes(col = polarization))+
  geom_path(data = data.table(curves) %>% arrange(Order),
            aes(UMAP1,UMAP2,group = Lineage)) +
  theme_classic(base_size=26)+
  scale_color_gradientn(name="Polarization",
                        colors = rev(c(rep(pal_npg(palette = c("nrc"))(2)[1],2),
                                       'lightgrey', 
                                       rep(pal_npg(palette = c("nrc"))(2)[2],2))),
                        values = scales::rescale(
                          c(min(u_dat$polarization), 2, -0.12, 2, max(u_dat$polarization)), 0:1)
  )+theme(aspect.ratio=1)+
  labs(y="Myeloid phenotype: UMAP 2",x="Myeloid phenotype: UMAP 1")
paperfile<- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Discovery and Validation Ribo and Letrozole Myeloid Polarization.png"),height=12,width=12, dpi=320)

sz <- 1.25
ggplot(u_dat,aes(UMAP1,UMAP2))+geom_point(aes(col = polarization),size=sz)+
  geom_path(data = data.table(curves) %>% arrange(Order),
            aes(UMAP1,UMAP2,group = Lineage),linewidth=sz) +
  theme_classic(base_size=26)+
  scale_color_gradientn(name="Polarization",
                        colors = rev(c(rep(pal_npg(palette = c("nrc"))(2)[1],2),
                                       'lightgrey', 
                                       rep(pal_npg(palette = c("nrc"))(2)[2],2))),
                        values = scales::rescale(
                          c(min(u_dat$polarization), 2, -0.12, 2, max(u_dat$polarization)), 0:1)
  )+  theme(aspect.ratio=1,axis.text=element_blank(),axis.title=element_blank(),legend.position="bottom",legend.title=element_blank())
paperfile<- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"BLANK Discovery and Validation Ribo and Letrozole Myeloid Polarization.png"),height=8,width=8, dpi=320)


ggplot( u_dat[] ,aes(UMAP1,UMAP2)) +    geom_point(aes(col = polarization),size=sz)+
  geom_path(data = data.table(curves) %>% arrange(Order),
            aes(UMAP1,UMAP2,group = Lineage),linewidth=sz) +
  scale_color_gradientn(name="Polarization",
                        colors = rev(c(rep(pal_npg(palette = c("nrc"))(2)[1],2),
                                       'lightgrey', 
                                       rep(pal_npg(palette = c("nrc"))(2)[2],2))),
                        values = scales::rescale(
                          c(min(u_dat$polarization), 2, -0.12, 2, max(u_dat$polarization)), 0:1)
  )+
  theme_classic()+
  theme(aspect.ratio=1,axis.text=element_blank(),axis.title=element_blank(),legend.text=element_blank(),legend.title=element_blank(),
        legend.key.size = unit(1.25, 'cm'))
#ggsave(paste0(paperfile,"BLANK Discovery and Validation Ribo and Letrozole Myeloid Polarization2.png"),height=8,width=8, dpi=320)

ggplot( u_dat[Cohort=="Validation"],aes(UMAP1,UMAP2,col=  Celltype_subtype)) +  
  geom_point(size=sz)+#scale_color_ucscgb(name="Myleoid cell types")+
  geom_point(data=u_dat[Cohort=="Discovery"],size=sz)+
  theme_classic()+
  theme(aspect.ratio=1,axis.text=element_blank(),axis.title=element_blank(),legend.position="none")


ggplot(u_dat,aes(UMAP1,UMAP2))+geom_point(aes(col = as.factor(classification)))+
  geom_path(data = data.table(curves) %>% arrange(Order),
            aes(UMAP1,UMAP2,group = Lineage)) +
  theme_classic(base_size=26)+
  theme(aspect.ratio=1, legend.position="none")+
  labs(y="Myeloid phenotype: UMAP 2",x="Myeloid phenotype: UMAP 1")
paperfile<- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Discovery and Validation Ribo and Letrozole Myeloid Clusters.png"),height=8,width=8, dpi=320)

ggplot(u_dat,aes(UMAP1,UMAP2))+geom_point(aes(col = as.factor(classification)))+
  #geom_path(data = data.table(curves) %>% arrange(Order),
  #          aes(UMAP1,UMAP2,group = Lineage)) +
  theme_classic(base_size=26)+
  theme(aspect.ratio=1, legend.position="none")+
  labs(y="Myeloid phenotype: UMAP 2",x="Myeloid phenotype: UMAP 1")
paperfile<- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Discovery and Validation Ribo and Letrozole Myeloid Clusters NoTrajB.png"),height=10,width=10, dpi=320)


ggplot(u_dat,aes(UMAP1,UMAP2))+geom_point(aes(col = as.factor(classification)))+
  geom_path(data = data.table(curves) %>% arrange(Order),
            aes(UMAP1,UMAP2,group = Lineage),linewidth=sz) +
  theme_classic(base_size=26)+
  theme(aspect.ratio=1, legend.position="none")+
  theme(aspect.ratio=1,axis.text=element_blank(),axis.title=element_blank())
#ggsave(paste0(paperfile,"BLANK Discovery and Validation Ribo and Letrozole Myeloid Clusters.png"),height=8,width=8, dpi=320)



ggplot(u_dat,aes(UMAP1,UMAP2))+geom_point(aes(col = common.pseudo))+
  geom_path(data = data.table(curves) %>% arrange(Order),
            aes(UMAP1,UMAP2,group = Lineage),size=sz) +
  scale_color_viridis_c(name="Pseudotime",option = "B")+
  theme_classic(base_size=26)+
  theme(aspect.ratio=1)+
  labs(y="Myeloid phenotype: UMAP 2",x="Myeloid phenotype: UMAP 1")
paperfile<- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Discovery and Validation Ribo and Letrozole Myeloid Average Pseudotime TrajB.png"),height=12,width=12, dpi=320)

#gam(formula=paste0("polarization ~s(",corgenes$rn[1],",k=4)") , family="negbin", data=u_dat)
#gamtst<-gam(formula=polarization ~s(scale(log(1+CD36)),k=3, bs="ts") ,theta=1, data=u_dat)
#summary(gamtst)

# Assess dynamically changing gene expression with polarization
assoctest<-rbindlist(lapply(1:length(colnames(u_datx)), function(i){
  gammd<-gam(formula=as.formula(paste0("polarization ~s(scale(log(1+",
                                       colnames(u_datx)[i],
                                       ")),k=3, bs='ts')")) , data=u_dat)
  data.table(rn=colnames(u_datx)[i], summary(gammd)$ s.table , r.sq=summary(gammd)$ r.sq)
}))
setnames(assoctest, old=c("p-value","F"), new=c("p.value","F.value"), skip_absent=TRUE)
assoctest[p.value<0.05]
ggplot(assoctest, aes(r.sq,edf))+geom_point()
assoctest[order(-r.sq)][1:10]
#write.csv(assoctest, file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/M1M2like gene list/Macrophage polarization gene associations.csv")
assoctest <- data.table(read.csv(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/M1M2like gene list/Macrophage polarization gene associations.csv") )

# Generate UMPA overlays for key genes
paneloverlay_u_dat<- function(gene,alpha=1,sz=1.25,save_im=F){
  ggplot( u_dat,aes_string("UMAP1","UMAP2",col=  paste0("log(1+",gene, ")"))) +  
    geom_point(size=sz,alpha=alpha)+#scale_color_ucscgb(name="Myleoid cell types")+
    theme_classic()+#facet_wrap(dynamic_class3~Cohort)+
    theme(aspect.ratio=1,axis.text=element_blank(),axis.title=element_blank(),legend.position="none")+
    ggtitle(gene)+ theme(plot.title = element_text(size=26))
  if(save_im==T){
    paperfile<- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
    ggsave(paste0(paperfile,"Discovery and Validation Ribo and Letrozole Myeloid MarkerExpression_",gene,".png"),height=8,width=8, dpi=320)
    
  }
}
paneloverlay_u_dat("CD80",save_im=F)
paneloverlay_u_dat("CD36",save_im=T)
paneloverlay_u_dat("DHRS9",save_im=T)
paneloverlay_u_dat("DOCK3",save_im=T)
paneloverlay_u_dat("FABP4",save_im=T)
paneloverlay_u_dat("PLA2G7",save_im=T)
paneloverlay_u_dat("FABP5",save_im=T)
paneloverlay_u_dat("LIPA",save_im=T)
paneloverlay_u_dat("CYP27A1",save_im=T)
paneloverlay_u_dat("MITF",save_im=T)
paneloverlay_u_dat("KCNE1",save_im=T)

negassoc <- data.table( cor(u_datx,u_dat %>%select(polarization)) , keep.rownames = T)[order(polarization)][1:40]
posassoc <- data.table( cor(u_datx,u_dat %>%select(polarization)) , keep.rownames = T)[order(-polarization)][1:40]

assoctest[rn%in%negassoc$rn][order(-r.sq)][1:10]
paneloverlay_u_dat("CCR5")
paneloverlay_u_dat("CD80")
paneloverlay_u_dat("TNFAIP3") #ok
paneloverlay_u_dat("CCL3")
paneloverlay_u_dat("CCL28")#good
paneloverlay_u_dat("IRF5")
paneloverlay_u_dat("CSF2RA")
paneloverlay_u_dat("TRPS1") #
paneloverlay_u_dat("CD40") # weak signal in DCs
paneloverlay_u_dat("FKBP5") # FKBP5 regulates M1 macrophage polarization
paneloverlay_u_dat("KYNU")#,sz=0.7,alpha=1)# weak signal in DCs
paneloverlay_u_dat("PALD1")#,sz=0.7,alpha=1)# weak signal in DCs
paneloverlay_u_dat("PALD1")#,sz=0.7,alpha=1)# lower right traj
paneloverlay_u_dat("CSF2RA",save_im=T)#,sz=0.7,alpha=1)# lower right traj
paneloverlay_u_dat("FLT3",save_im=T)#,sz=0.7,alpha=1)# lower left DC traj
paneloverlay_u_dat("INPP5D")#,sz=0.7,alpha=1)# weak lower section signal
paneloverlay_u_dat("ZEB1",save_im=T)#cDC1 marker
paneloverlay_u_dat("ABCC4")#Migratory DC marker
paneloverlay_u_dat("ENOX1",save_im=T)# general DC marker
paneloverlay_u_dat("CIITA" ,save_im=T)# M1 marker: Macrophage-specific MHCII expression is regulated by a remote Ciita # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6219740/ 
paneloverlay_u_dat("ERBB4",save_im=T)# M1 marker: Macrophage-specific MHCII expression is regulated by a remote Ciita # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6219740/ 


#monocyte-derived macrophages #https://medicine.uky.edu/sites/default/files/file-uploads/2023-10/Cell-autonomous%20regulation%20of%20complement%20C3%20by%20factor%20H%20limits%20macrophage%20efferocytosis%20and%20exacerbates%20atherosclerosis.pdf
paneloverlay_u_dat("C3",save_im=T)# Central marker :C3 Produced by Macrophages Promotes Renal Fibrosis via IL-17A Secretion# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6204358/
paneloverlay_u_dat("WNT5B",save_im=T) #Migratory DC marker #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7039855/

# Heatmap of changes in gene expression with polarization
topgenes <- c("CD36","DHRS9","PLA2G7","DOCK3","FABP4","FABP5","LIPA","CYP27A1", "MITF", "KCNE1",
              "CSF2RA", "CIITA","C3","ERBB4","ENOX1","ZEB1","FLT3","WNT5B")
pst.ord <- order(u_dat$polarization, na.last = NA)
heatdata <- u_datx[pst.ord, topgenes ]
heatclus <- u_dat$classification[pst.ord]

# color pallete definition and plotting
gg_color_hue <- function(n) {  hcl(h = seq(15, 375, length = n + 1), l = 65, c = 100)[1:n] }
heatmap(log1p(t(heatdata)), Rowv=NA,scale="row",revC=T,Colv = NA, ColSideColors =  gg_color_hue(max(heatclus))[heatclus])

# Larger gene list if needed
#topgenes <- c(assoctest[rn%in%negassoc$rn][order(-r.sq)][1:20]$rn,
#              assoctest[rn%in%posassoc$rn][order(-r.sq)][1:20]$rn)
#pst.ord <- order(u_dat$polarization, na.last = NA)
#heatdata <- u_datx[pst.ord, topgenes ]
#heatclus <- u_dat$classification[pst.ord]
#gg_color_hue <- function(n) {  hcl(h = seq(15, 375, length = n + 1), l = 65, c = 100)[1:n] }
#heatmap(log1p(t(heatdata)), scale="row",revC=T,Colv = NA, ColSideColors =  gg_color_hue(max(heatclus))[heatclus])

# load existing m2/m1 signature list
savloc<- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/SI data/FigureS9/"
markersdd <- data.table(read.csv(paste0(savloc,"SourceData_FigureS9_Myeloid marker gene list.csv")))
tmp <- unique(markersdd )[order(-State,Gene)]

# Calculate GSEA scires
require(GSVA)
geneSets <-list(m2Set=unique(markersdd$Gene)[unique(markersdd$Gene)%in%colnames(u_datx)])
#tmp[,M2geneset:=FALSE]
#tmp[Gene%in%geneSets[[1]],M2geneset:=TRUE]
#write.csv(tmp,file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/M1M2like gene list/Referenced M1M2like gene list.csv")

gsvapar <- gsvaParam(t(u_datx), geneSets, maxDiff=TRUE)
gsva_es <- gsva(gsvapar)
# write.csv(gsva_es,"/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/M1M2like gene list/M1M2like gene gsva M2 sig.csv")
#gsva_es <- read.csv("/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/M1M2like gene list/M1M2like gene gsva M2 sig.csv")[-1]

u_dat$gsva_scores<- unlist(as.vector(gsva_es))

ggplot(u_dat,
       aes(y=polarization,x=gsva_scores,col=UMAP2)) + 
  geom_point(size=0.6)


# gsva enrichment scores
ggplot(u_dat, aes(UMAP1,UMAP2)) + geom_point(aes(col = gsva_scores))+
  geom_path(data = data.table(curves) %>% arrange(Order),
            aes(UMAP1, UMAP2, group = Lineage)) +
  scale_color_viridis_c(name="M2 signature \n (GSVA \nenrichment \n scores)", option = "B")+
  theme_classic() +
  theme(aspect.ratio= 1) +
  labs(y= "Myeloid phenotype: UMAP 2", x= "Myeloid phenotype: UMAP 1")


# PLAGE pathway scores (svd pca) of macrophage markers
pcamod <- prcomp(log(1+u_datx[,geneSets[[1]]]))
plot((pcamod$x[,1:2]))
markerloadings<- data.table(pcamod$rotation[,1,drop=F],keep.rownames=T)[order(-abs(PC1))]
u_dat$PLAGE<- -pcamod$x[,1]
gam(PLAGE~s(polarization,k=3) , data=u_dat)%>%summary()
gam(gsva_scores~s(polarization,k=3) , data=u_dat)%>%summary()
gam(gsva_scores~s(UMAP2) , data=u_dat)%>%summary()
cor(u_dat$polarization,u_dat$PLAGE, method="pearson")
cor(u_dat$gsva_scores,u_dat$PLAGE, method="pearson")
cor(u_dat[ARM!="A"]$dynamic_class3=="Response",u_dat[ARM!="A"]$PLAGE)

cordd<-u_dat%>%select(polarization,gsva_scores,PLAGE,UMAP2)
setnames(cordd,old=c("polarization","gsva_scores","PLAGE","UMAP2" ),new=c("M2-like \n Polarization","M2-like \n GSVA \n score","M2-like \n PLAGE \n score","M2-like \n Differentiation" ))
corrplot::corrplot(cor(cordd, method="pearson"),
                   type="upper",diag = F,addCoefasPercent=T,order="hclust",
                   col=colorRampPalette(c("blue","white","red"))(100),
)$corrPos -> p1

write.csv()
cordd
#savlocS9<-"/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/SI data/FigureS9/"
#write.csv(cordd, file=paste0(savlocS9,"SourceData_FigureS9_Myeloid Phenotype Metrics_Out.csv"))

corrplot::corrplot.mixed(cor(cordd, method="pearson"))
corrplot::corrplot(cor(cordd, method="pearson"))
ggplot(u_dat, aes(UMAP1,UMAP2)) + geom_point(size=2.5,aes(col = PLAGE))+
  geom_path(data = data.table(curves) %>% arrange(Order),
            aes(UMAP1, UMAP2, group = Lineage),linewidth=2) +
  scale_color_viridis_c(name="M2 signature \n (PLAGE \nenrichment \n scores)", option = "B")+
  theme_classic(base_size=26) +
  theme(aspect.ratio= 1) +
  labs(y= "Myeloid phenotype: UMAP 2", x= "Myeloid phenotype: UMAP 1")
paperfile<- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Discovery and Validation Ribo and Letrozole Myeloid PLAGE ovelay.png"),height=10,width=10, dpi=320)

ggplot(u_dat, aes(UMAP1,UMAP2)) + geom_point(size=2.5,aes(col = gsva_scores))+
  geom_path(data = data.table(curves) %>% arrange(Order),
            aes(UMAP1, UMAP2, group = Lineage),linewidth=2) +
  scale_color_viridis_c(name="M2 signature \n (GSVA \nenrichment \n scores)", option = "B")+
  theme_classic(base_size=26) +
  theme(aspect.ratio= 1) +
  labs(y= "Myeloid phenotype: UMAP 2", x= "Myeloid phenotype: UMAP 1")
paperfile<- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Discovery and Validation Ribo and Letrozole Myeloid GSVA ovelay.png"),height=10,width=10, dpi=320)


ggplot(u_dat, aes(y=gsva_scores, x=polarization)) + 
  geom_point(alpha=0.3,aes(col = PLAGE))+
  theme(aspect.ratio=1)+
  scale_color_viridis_c(name="M2-like gene set \n PLAGE \n enrichment score", option = "B")+
  theme_classic() +
  geom_smooth(method="gam",formula= y~s(x,k=4))+ theme(aspect.ratio= 1) +
  labs(x="Polarization (M1/M2-like differentiation)", y="M2-like gene set \n GSVA enrichment score")

ggplot(u_dat, aes(y=PLAGE, x=polarization)) + 
  geom_point(size=2.5,alpha=0.3,aes(col = gsva_scores))+
  theme(aspect.ratio=1)+
  scale_color_viridis_c(name="M2-like \n GSVA score \n gene set \n enrichment ", option = "B")+
  theme_classic(base_size=26) +
  geom_smooth(linewidth=2,col="black",se=F,method="gam",formula= y~s(x,k=3))+ theme(aspect.ratio= 1) +
  #labs(x="Polarization (M1/M2-like pseudotime differentiation)", y="M2-like gene set \n PLAGE enrichment score")+
  labs(x="M2-like Polarization (pseudotime score)", y="M2-like PLAGE score \n (gene set enrichment)")

paperfile<- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Discovery and Validation Ribo and Letrozole Myeloid Comparison of PLAGE GSVA and Polarization.png"),height=11,width=11, dpi=320)


#M2-like \n macrophage \n marker
ggplot(u_dat, aes(y=gsva_scores, x=polarization)) + 
  geom_point(alpha=0.3,aes(col = log(1+MSR1)))+
  theme(aspect.ratio=1)+
  scale_color_viridis_c(name="MSR1 \n expression \n log(1+x)", option = "B")+
  theme_classic() +
  geom_smooth(method="gam",formula= y~s(x,k=4))+ theme(aspect.ratio= 1) +
  labs(x="Polarization (M1/M2-like differentiation)", y="M2-like gene set \n GSVA enrichment score")

ggplot(u_dat, aes(y=PLAGE, x=polarization)) + geom_point(size=2,alpha=0.3,aes(col = log(1+MSR1)))+
  theme(aspect.ratio=1)+
  scale_color_viridis_c(name="MSR1 \n expression \n log(1+x)", option = "B")+
  theme_classic(base_size=26) +
  geom_smooth(method="gam",formula= y~s(x,k=3))+ theme(aspect.ratio= 1) +
  labs(x="M2-like Polarization (pseudotime score)", y="M2-like gene set \n PLAGE enrichment score")
paperfile<- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Discovery and Validation Ribo and Letrozole MSR1 marker M2 Myeloid polarization.png"),height=10,width=10, dpi=320)

ggplot(u_dat, aes(y=PLAGE, x=polarization)) + geom_point(size=2,alpha=0.3,aes(col = log(1+CTSB)))+
  theme(aspect.ratio=1)+
  scale_color_viridis_c(name="CTSB \n expression \n log(1+x)", option = "B")+
  theme_classic(base_size=26) +
  geom_smooth(method="gam",formula= y~s(x,k=3))+ theme(aspect.ratio= 1) +
  labs(x="M2-like Polarization (pseudotime score)", y="M2-like gene set \n PLAGE enrichment score")
paperfile<- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Discovery and Validation Ribo and Letrozole CTSB marker M2 Myeloid polarization.png"),height=10,width=10, dpi=320)

ggplot(u_dat, aes(y=PLAGE, x=polarization)) + geom_point(size=2,alpha=0.3,aes(col = log(1+CSF1R)))+
  theme(aspect.ratio=1)+
  scale_color_viridis_c(name="CSF1R \n expression \n log(1+x)", option = "B")+
  theme_classic(base_size=26) +
  geom_smooth(method="gam",formula= y~s(x,k=3))+ theme(aspect.ratio= 1) +
  labs(x="M2-like Polarization (pseudotime score)", y="M2-like gene set \n PLAGE enrichment score")
paperfile<- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Discovery and Validation Ribo and Letrozole CSF1R marker M2 Myeloid polarization.png"),height=10,width=10, dpi=320)

ggplot(u_dat, aes(y=PLAGE, x=polarization)) + geom_point(size=2,alpha=0.4,aes(col = log(1+CD163)))+
  theme(aspect.ratio=1)+
  scale_color_viridis_c(name="CD163 \n expression \n log(1+x)", option = "B")+
  theme_classic(base_size=26) +
  geom_smooth(method="gam",formula= y~s(x,k=3))+ theme(aspect.ratio= 1) +
  labs(x="M2-like Polarization (pseudotime score)", y="M2-like gene set \n PLAGE enrichment score")
paperfile<- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Discovery and Validation Ribo and Letrozole CD163 marker M2 Myeloid polarization.png"),height=10,width=10, dpi=320)

ggplot(u_dat, aes(y=PLAGE, x=polarization)) + geom_point(size=2,alpha=0.3,aes(col = log(1+MRC1)))+
  theme(aspect.ratio=1)+
  scale_color_viridis_c(name="MRC1 \n expression \n log(1+x)", option = "B")+
  theme_classic(base_size=26) +
  geom_smooth(method="gam",formula= y~s(x,k=3))+ theme(aspect.ratio= 1) +
  labs(x="M2-like Polarization (pseudotime score)", y="M2-like gene set \n PLAGE enrichment score")
paperfile<- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Discovery and Validation Ribo and Letrozole MRC1 marker M2 Myeloid polarization.png"),height=10,width=10, dpi=320)

ggplot(u_dat, aes(y=gsva_scores, x=UMAP2)) + geom_point(size=2,alpha=0.1,aes(col = gsva_scores))+
  scale_color_viridis_c(name="M2 signature \n (GSVA \nenrichment \n scores)", option = "B")+
  theme_classic(base_size=26) +
  geom_smooth(method="lm")+ theme(aspect.ratio= 1) #+

enrichCompar <- data.table( u_dat%>%group_by(Cohort,TreatmentLab,Patient.Study.ID, Day,dynamic_class3)%>%
                              summarise(UMAP2=mean(UMAP2),polarization=mean(polarization), gsva_scores=mean(gsva_scores),PLAGE=mean(PLAGE)) )

lm(gsva_scores~I(dynamic_class3=="Response"), data=enrichCompar[Day<180][TreatmentLab=="Combination ribociclib"]) %>% summary()
lm(PLAGE~I(dynamic_class3=="Response"), data=enrichCompar[Day<180][TreatmentLab=="Combination ribociclib"]) %>% summary()
lm(polarization~I(dynamic_class3=="Response"), data=enrichCompar[Day<180][TreatmentLab=="Combination ribociclib"]) %>% summary()

lm(gsva_scores~I(dynamic_class3=="Response"), data=enrichCompar[Day<180][TreatmentLab!="Combination ribociclib"]) %>% summary()
lm(PLAGE~I(dynamic_class3=="Response"), data=enrichCompar[Day<180][TreatmentLab!="Combination ribociclib"]) %>% summary()
lm(polarization~I(dynamic_class3=="Response"), data=enrichCompar[Day<180][TreatmentLab!="Combination ribociclib"]) %>% summary()



glm(I(dynamic_class3=="Response") ~ gsva_scores, data=enrichCompar[Day<180][TreatmentLab=="Combination ribociclib"], family= "binomial") %>% summary()
glm(I(dynamic_class3=="Response") ~ gsva_scores, data=enrichCompar[Day<180][TreatmentLab!="Combination ribociclib"], family= "binomial") %>% summary()
glm(I(dynamic_class3=="Response") ~ gsva_scores*TreatmentLab, data=enrichCompar[Day<180], family= "binomial") %>% summary()
glm(I(dynamic_class3=="Response") ~ polarization*TreatmentLab, data=enrichCompar[Day<180], family= "binomial") %>% summary()

ploty<- data.table(enrichCompar%>%gather(var,val, UMAP2:PLAGE))[var!="UMAP2"]
ploty[,method:="M2-like \n Polarization"]
ploty[var=="gsva_scores",method:="M2-like \n GSVA score"]
ploty[var=="PLAGE",method:="M2-like \n PLAGE score"]
ggplot(ploty[Day<180], 
       aes(y=val, x=dynamic_class3, fill=dynamic_class3)) + 
  geom_boxplot(col="black",outlier.color = NA)+facet_grid(method~TreatmentLab, scales="free_y")+
  stat_boxplot(geom = "errorbar",
               width = 0.45) +
  geom_jitter(aes(shape=Cohort),height=0,width=0.05,size=1.5*sz)+
  scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  theme_classic(base_size=26) +
  theme(aspect.ratio= 1) +labs(y="Myeloid M2-like phenotype (mean per tumor)", x="Tumor response")
paperfile<- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Discovery and Validation Ribo and Letrozole Boxplot Myeloid differentiation metrics all predict response to ribo not letrozole.png"),height=20,width=13.5, dpi=320)


EarlyMyeloidDiffMetrics <- ploty[Day<180]
#savlocS11<-"/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/SI data/FigureS11/"
#write.csv(EarlyMyeloidDiffMetrics, file=paste0(savlocS11,"SourceData_FigureS11_Early Myeloid Phenotype Differentiation Metrics_Out.csv"))
           


ggplot(enrichCompar, 
       aes(y=PLAGE, x=dynamic_class3, fill=dynamic_class3)) + 
  geom_boxplot(col="black",outlier.color = NA)+facet_grid(.~TreatmentLab)+
  stat_boxplot(geom = "errorbar",
               width = 0.45) +
  geom_jitter(aes(shape=Cohort),height=0,width=0.05)+
  scale_color_npg()+
  theme_classic() +
  theme(aspect.ratio= 1) 

ggplot(enrichCompar[Day<180], 
       aes(y=polarization, x=dynamic_class3, fill=dynamic_class3)) + 
  geom_boxplot(col="black",outlier.color = NA)+facet_grid(.~TreatmentLab)+
  stat_boxplot(geom = "errorbar",
               width = 0.45) +
  geom_jitter(aes(shape=Cohort),height=0,width=0.08)+
  scale_color_npg()+
  theme_classic() +
  theme(aspect.ratio= 1) +
  labs(y="Mean ")


enrichCompar[,success:=0]
enrichCompar[dynamic_class3=="Response",success:=1]
glm1<-glm(success~PLAGE*TreatmentLab, family = "binomial", data=enrichCompar)
glm1<-glm(success~polarization*TreatmentLab, family = "binomial", data=enrichCompar)
enrichCompar$glmPred<- predict(glm1, type="response")

ggplot(enrichCompar, 
       aes(x=polarization, y=1*success)) + 
  geom_boxplot(aes(fill=dynamic_class3),col="black",outlier.color = NA)+
  facet_grid(.~TreatmentLab)+
  stat_boxplot(aes(group=dynamic_class3),geom = "errorbar",
               width = 0.45) +
  geom_jitter(aes(shape=Cohort),height=0,width=0.05)+
  scale_color_npg()+
  theme_classic() +
  geom_smooth(stat="identity",aes(y=glmPred, group=TreatmentLab))+
  theme(aspect.ratio= 1) 


ggplot(enrichCompar,
       aes(y=PLAGE,x=log(1+Day),col=dynamic_class3)) + 
  geom_point(size=0.6)+facet_grid(.~TreatmentLab)+
  geom_smooth(method="gam", formula=y~s(x, k=3))+
  scale_x_continuous(breaks=log(1+c(0,14,180)), labels=c(0,14,180) )+
  labs(x="Day", y="Myeloid M2-like PLAGE score")+
  scale_color_npg(name="Tumor \n response", labels=c("Resistant","Sensitive"))

ggplot(enrichCompar,
       aes(y=gsva_scores,x=log(1+Day),col=dynamic_class3)) + 
  geom_point(size=0.6)+facet_grid(.~TreatmentLab)+
  geom_smooth(method="gam", formula=y~s(x, k=3))+
  scale_x_continuous(breaks=log(1+c(0,14,180)), labels=c(0,14,180) )+
  labs(x="Day", y="Myeloid M2-like GSVA score")+
  scale_color_npg(name="Tumor \n response", labels=c("Resistant","Sensitive"))

getplotdd<-data.table(enrichCompar%>%select(-c(success,glmPred))%>%gather(var,val,UMAP2, polarization, gsva_scores,       PLAGE ))
getplotdd[var=="UMAP2",var:= "M2-like \n Differentiation"]
getplotdd[var=="PLAGE",var:= "M2-like \n PLAGE score"]
getplotdd[var=="gsva_scores",var:= "M2-like \n GSVA score"]
getplotdd[var=="polarization",var:= "M2-like \n Polarization"]
ggplot(getplotdd,aes(y=val,x=log(1+Day),col=dynamic_class3,fill=dynamic_class3)) + 
  geom_point(size=2)+facet_grid(var~TreatmentLab, scales="free_y")+
  theme_classic(base_size=26)+
  geom_smooth(method="gam", formula=y~s(x, k=3))+
  theme(aspect.ratio = 1)+
  scale_x_continuous(breaks=log(1+c(0,14,180)), labels=c(0,14,180) )+
  labs(x="Day", y="Myeloid M2-like phenotype (mean per sample)")+
  scale_color_npg(name="Tumor \n response", labels=c("Resistant","Sensitive"))+
  scale_fill_npg(name="Tumor \n response", labels=c("Resistant","Sensitive"))

paperfile<- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Discovery and Validation Ribo and Letrozole Myeloid M2-like phenotype metrics over time.png"),height=14,width=13, dpi=320)

TemporalMyeloidDiffMetrics <- ploty[Day<180]
#savlocS12<-"/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/SI data/FigureS12/"
#write.csv(TemporalMyeloidDiffMetrics, file=paste0(savlocS12,"SourceData_FigureS12_Temporal Myeloid Phenotype Differentiation Metrics_Out.csv"))


ggplot(enrichCompar,
       aes(y=polarization,x=log(1+Day),col=dynamic_class3)) + 
  geom_point(size=0.6)+facet_grid(.~TreatmentLab)+
  geom_smooth(method="gam", formula=y~s(x, k=3))+
  scale_x_continuous(breaks=log(1+c(0,14,180)), labels=c(0,14,180) )+
  labs(x="Day", y="Myeloid M2-like polarization")+
  scale_color_npg(name="Tumor \n response", labels=c("Resistant","Sensitive"))

ggplot(enrichCompar,
       aes(y=UMAP2,x=log(1+Day),col=dynamic_class3)) + 
  geom_point(size=0.6)+facet_grid(.~TreatmentLab)+
  geom_smooth(method="gam", formula=y~s(x, k=3))+
  scale_x_continuous(breaks=log(1+c(0,14,180)), labels=c(0,14,180) )+
  labs(x="Day", y="Myeloid M2-like differentiation")+
  scale_color_npg(name="Tumor \n response", labels=c("Resistant","Sensitive"))




#### Ealy treatment myeloid differentiation of average myeloid cell vs response
u_dat[,Treatment:="letrozole + ribo"]
u_dat[ARM=="A",Treatment:="letrozole"]

summarizedCohorts<- data.table(u_dat[Day!=180] %>% group_by(Treatment,#Day,
                                                               Patient.Study.ID,dynamic_class3,Cohort,ARM,Celltype) %>%
                                 summarise(n= n(),UMAP1= mean(UMAP1), UMAP2= mean(UMAP2) ,n= n() ) %>%
                                 group_by(Cohort) %>% mutate(cohortmu= mean(UMAP2))%>%
                                 group_by(Treatment,Cohort)%>%mutate(cohortTreatmu= median(UMAP2)))

summarizedCohorts[,Success:=1]
summarizedCohorts[dynamic_class3=="Non-response",Success:=0]
summarizedCohorts[,TumorDifferScore:=UMAP2-cohortmu ]

ggplot(summarizedCohorts, 
       aes(y=UMAP2-cohortmu, x=dynamic_class3,col=Cohort, group=interaction(Cohort,dynamic_class3) )) +
  geom_violin(scale="width")+
  geom_point(aes(shape=Treatment),positio=position_dodge(width=0.9))+
  theme_classic()+theme(aspect.ratio = 1)+facet_wrap(~Treatment)

summarizedCohorts[,Treatmentlab:= "Combination ribociclib"]
summarizedCohorts[Treatment=="letrozole",Treatmentlab:= "Letrozole alone"]


ggplot(summarizedCohorts,  aes(y=TumorDifferScore, x=dynamic_class3,col=dynamic_class3, fill=dynamic_class3, shape=Cohort, group=interaction(dynamic_class3) )) +
  # geom_boxplot(scale="width")+
  # stat_boxplot(scale="errorbar",linetype=1,width=2,col="black")+
  
  geom_boxplot(outlier.colour=NA, position= position_dodge() ,col="black")+
  stat_boxplot(geom="errorbar",position=position_dodge(1.75),width=0.5)+#geom_smooth(method="lm")+
  
  geom_point(col="black",size=2.5)+
  scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"),guide="none")+
  scale_fill_npg(name="Tumor response", labels=c("Resistant", "Sensitive"),guide="none")+
  scale_x_discrete(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  theme_classic(base_size=26)+theme(aspect.ratio = 1)+facet_wrap(~Treatmentlab,ncol=1)+
  labs(y="Myeloid M2 differentiation \n (tumor mean)")+theme(legend.position="bottom")
#ggsave(filename=paste0("/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cohort 2 results/Boxplot Myeloid differentiation predicts response to ribo not letrozole.png"),width=5,height=5)
#paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Ribo and Letrozole Boxplot Myeloid differentiation predicts response to ribo not letrozole.png"),height=10,width=10)

output<-summarizedCohorts%>%select(Treatment,Treatmentlab,Patient.Study.ID,Cohort,dynamic_class3,TumorDifferScore)
#savloc4<-"/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/Figure4/"
#write.csv(output, file=paste0(savloc4,"SourceData_Figure4_Early treatment Myeloid M2like Phenotype_Out.csv"))

summary(lm(TumorDifferScore~dynamic_class3*Treatmentlab,summarizedCohorts))
lm(TumorDifferScore~dynamic_class3*Treatmentlab, summarizedCohorts)%>%summary()
lm(TumorDifferScore~dynamic_class3, summarizedCohorts[Treatmentlab=="Letrozole alone"])%>%summary()
lm(TumorDifferScore~dynamic_class3, summarizedCohorts[Treatmentlab!="Letrozole alone"])%>%summary()






#### Ealy treatment fraction of myeloid cells in  M2-like state vs response

# Load cell annotation data
#data.loc<- "/Users/jason/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/Figure2/"
data.loc<- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/Figure2/"
data.name <-"SourceData_Figure2_CellAnnotations.csv"
dd0 <- data.table(read.csv( file=paste0(data.loc, data.name)))

# count of each cell type's abundance and frequency
countTable0 <- data.table( dd0 %>% group_by(Patient.Study.ID, Day, Cohort, Treatmentlab, dynamic_class3,
                                            CellAnnot2,totalCount) %>% summarise(count=n() , frac=n()/unique(totalCount) ) ) 

countTable0wideM1M2<-data.table(countTable0[CellAnnot2 %in% unique(dd0[CellAnnot1=="Macrophages"]$CellAnnot2) ]%>% select(-count)%>%spread(CellAnnot2,frac,fill=0) )
names(countTable0wideM1M2) <- gsub(" ","_",names(countTable0wideM1M2))
countTable0wideM1M2[,M1M2ratio:= M1_macrophages/(  M1_macrophages + M2_macrophages ) ]
countTable0wideM1M2[,Myeloidtot := ( DC + M1_macrophages + M2_macrophages + Monocytes)  ]

countTable0wideM1M2[,Success:=1]
countTable0wideM1M2[dynamic_class3=="Non-response",Success:=0]
countTable0wideM1M2<-countTable0wideM1M2[is.finite(M1M2ratio)]
countTable0wideM1M2[,x:=M1M2ratio ]
countTable0wideM1M2[,y:=log(Myeloidtot) ]


ggplot(countTable0wideM1M2[Day!=180][(Myeloidtot*totalCount)>20], 
       aes(y= (M1M2ratio), x= dynamic_class3, fill= dynamic_class3 )) + theme_classic(base_size= 28)+
  #geom_violin(aes(group= interaction(dynamic_class3)),alpha=0.6) + 
  geom_boxplot()+
  geom_boxplot(outlier.colour=NA, position= position_dodge() ,col="black")+
  stat_boxplot(geom="errorbar",position=position_dodge(1.75),width=0.5)+#geom_smooth(method="lm")+
  facet_wrap(~ Treatmentlab,ncol=1) + theme(aspect.ratio= 1) +
  #geom_smooth(method="gam", formula=y~s(x,k=3),se=F,size=2)+
  scale_fill_npg(name= "Tumor response", labels= c("Resistant", "Sensitive"))+
  labs(y="Proportion myeloid cells in M1 state",x="Tumor response") + 
  geom_point(aes(shape=Cohort),size=2.5, col="black")+theme(legend.position="none")+
  scale_x_discrete(labels=c("Resistant", "Sensitive"))

#paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Discovery and Validation Ribociclib and Letrozole Early treatment M1 abundance high in response to ribo not letrozole boxplot.png"),height=10,width=10, dpi=320)


ggplot(countTable0wideM1M2[Day!=180], 
       aes(y= (M1M2ratio), x= dynamic_class3, fill= dynamic_class3 )) + theme_classic(base_size= 28)+
  #geom_violin(aes(group= interaction(dynamic_class3)),alpha=0.6) + 
  geom_boxplot()+
  geom_boxplot(outlier.colour=NA, position= position_dodge() ,col="black")+
  stat_boxplot(geom="errorbar",position=position_dodge(1.75),width=0.5)+#geom_smooth(method="lm")+
  facet_wrap(~ Treatmentlab,ncol=1) + theme(aspect.ratio= 1) +
  #geom_smooth(method="gam", formula=y~s(x,k=3),se=F,size=2)+
  scale_fill_npg(name= "Tumor response", labels= c("Resistant", "Sensitive"))+
  labs(y="Proportion myeloid cells in M1 state",x="Tumor response") + 
  geom_point(aes(shape=Cohort),size=2.5, col="black")+theme(legend.position="none")+
  scale_x_discrete(labels=c("Resistant", "Sensitive"))

#paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Discovery and Validation Ribociclib and Letrozole Early treatment M1 abundance high in response to ribo not letrozole boxplotB.png"),height=10,width=10, dpi=320)


glmmodC<-lm( M1M2ratio~ dynamic_class3*Day*Cohort, data=countTable0wideM1M2[Day!=180][Treatmentlab=="Combination ribociclib"])
glmmodL <-lm( M1M2ratio ~ dynamic_class3*Day*Cohort, data=countTable0wideM1M2[Day!=180][Treatmentlab=="Letrozole alone"] )
summary(glmmodC)
summary(glmmodL)

resdat<- countTable0wideM1M2[Day!=180] %>%select(M1M2ratio,dynamic_class3,Treatmentlab,Cohort) 
#savloc4<-"/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/Figure4/"
#write.csv(output, file=paste0(savloc4,"SourceData_Figure4_Early treatment Myeloid M2like Phenotype_Out.csv"))
write.csv(resdat , file=paste0(savloc4,"SourceData_Figure4_Early treatment Myeloid M1like state proportion_Out.csv") )
#save(resdat, file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/Figure4/Discovery and Validation Myeloid M1 proportion early treatment.RData")
load(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/Figure4/Discovery and Validation Myeloid M1 proportion early treatment.RData")
#resdat




