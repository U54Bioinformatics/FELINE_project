rm(list=ls())
library(glmnet);require(dplyr);require(data.table);require(tidyr)
require(ggraph);require(tidygraph)

load( "~/Dropbox/Cancer_pheno_evo/data/FELINE2/LigandReceptor/Filtered_Human-2015-Ramilowski-LR-pairs.RData" )
load(file ="~/Dropbox/Cancer_pheno_evo/data/FELINE2/cohort2 Communication output merged/PopulationCommunicationMerged ALL ARMS cohort2.RData" )#CCI,allphenotypes, uu,perIndiv,
LRgenelist <- unique( data.table(merge(CCI[],LRpairsFiltered,by="Pair.Name"))$HPMR.Receptor) #unique( c(LRpairsFiltered$HPMR.Receptor) ) #, LRpairsFiltered$HPMR.Ligand) )
rm(list= c("LRpairs", "LRpairsFiltered"))

CCI[,Treat:="CombinationRibo"]
CCI[ARM=="A",Treat:="LetrozoleAlone"]
CCI[,scaleTransduction:=scale(TransductionMu,center=F),by=c("Pair.Name")] 
CCI[!is.finite(TransductionMu),scaleTransduction:=0]
CCI[,scalelnTransductionMu:=scale(log(1+TransductionMu)),by=Pair.Name]
CCI[,scaleTransduction2:=exp(scale(log(TransductionMu),center=F)),by=c("Pair.Name")] 

load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/GeneOntology/growthFactorReceptors.RData" )

wchComms<- unique( c( "EGFR","HBEGFR","TGFBR1","TGFBR2","TGFBR3"))

isApath<-sapply( 1:length(wchComms),function(ii){grepl(paste0("_",wchComms[ii]),unique(CCI$Pair.Name) ) } ) 
l1 <- unique(CCI$Pair.Name)[rowSums(isApath)>0 ]

average_muln_scaleTransduction <- data.table(CCI[Pair.Name%in%l1]%>%
                                               group_by(Patient.Study.ID,LigandPhenoCelltype,ReceptorPhenoCelltype,Day,ARM,Treat) %>%
                                               dplyr::summarise(muln_scaleTransduction=mean(log(1+scaleTransduction))) %>% 
                                               group_by(LigandPhenoCelltype, ReceptorPhenoCelltype, Day, Treat) %>%
                                               dplyr::summarise(average_muln_scaleTransduction=mean(muln_scaleTransduction)))
average_muln_scaleTransduction$LigandPhenoCelltype <- factor(average_muln_scaleTransduction$LigandPhenoCelltype, 
                                                             levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Normal epithelial cells","Cancer cells","B cells","Macrophages","CD8+ T cells","CD4+ T cells"))
average_muln_scaleTransduction$ReceptorPhenoCelltype <- factor(average_muln_scaleTransduction$ReceptorPhenoCelltype, 
                                                               levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Normal epithelial cells","Cancer cells","B cells","Macrophages","CD8+ T cells","CD4+ T cells") )

initState <- data.table(average_muln_scaleTransduction[Day==0]%>%
                          group_by(LigandPhenoCelltype,ReceptorPhenoCelltype)%>% #Treat, dynamic_class3
                          dplyr::summarise(intiState=mean(average_muln_scaleTransduction)) )

average_muln_scaleTransduction2 <- merge(average_muln_scaleTransduction,initState,by=c("LigandPhenoCelltype", "ReceptorPhenoCelltype")) #,"Treat","dynamic_class3"
average_muln_scaleTransduction2[ , lnfoldchange:= (average_muln_scaleTransduction-intiState) ]      #average_muln_scaleTransduction2[ , lnfoldchange:= (average_muln_scaleTransduction) ] 


CCIcondensed <- average_muln_scaleTransduction2[LigandPhenoCelltype!="B cells"][ReceptorPhenoCelltype!="B cells"]  
CCIcondensed$LigandPhenoCelltype <-factor(CCIcondensed$LigandPhenoCelltype, 
                                          levels= sort(c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Normal epithelial cells","Cancer cells","Macrophages","CD8+ T cells","CD4+ T cells"))) #"B cells",
CCIcondensed$ReceptorPhenoCelltype <-factor(CCIcondensed$ReceptorPhenoCelltype, 
                                            levels= sort(c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Normal epithelial cells","Cancer cells","Macrophages","CD8+ T cells","CD4+ T cells") )) #"B cells",
setnames(CCIcondensed, old=c("LigandPhenoCelltype", "ReceptorPhenoCelltype", "average_muln_scaleTransduction"), new= c("from", "to", "weight"))

tmp <- CCIcondensed[][][to=="Fibroblasts"][order( Treat,from,to)]
tmp$DayLab <- paste0("Day ", tmp$Day )
tmp <- data.table( tmp %>% group_by(Treat) %>% mutate(weight = scale(weight) ) )
tmp[,Treatmentlab:= "Letrozole alone"]
tmp[Treat=="CombinationRibo",Treatmentlab:= "Combination ribociclib"]

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
tmp[from=="Macrophages",from:="Myeloid \n cells"]
tmp[to=="Macrophages",to:="Myeloid \n cells"]
graphdd <- as_tbl_graph( tmp, directed= T)


node_col <- rev(ggsci::pal_npg("nrc")(3)  )

# circle layout
n <- length(unique(CCIcondensed$to)) -1
pts.circle <- t(sapply(1:n, function(r)c(cos(2*r*pi/n),sin(2*r*pi/n))))
NodeList <- data.table(c("Cancer cells", "CD4+ T cells" , "Adipocytes", "Fibroblasts", "Normal epithelial cells", "Pericytes"  , "Endothelial cells" ,"Macrophages", "CD8+ T cells" )  , c(0, pts.circle[,1] )  ,  c(0, pts.circle[,2] ) ) #"B cells",
presloc <- data.frame(NodeList)[order(NodeList$V1),]

createdLayout <- create_layout(graphdd, layout= "star")
createdLayout$x <- presloc$V2
createdLayout$y <- presloc$V3



p0 <- ggraph(createdLayout[] , aes(x= x, y= y , col= Day ) ) + 
  geom_edge_loop(  aes(edge_colour=weight,edge_width= 150*exp(weight),edge_alpha= (weight)) ,  arrow= arrow(length= unit(10, "mm" ) , type= "closed"),
                   start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" ) ) +
  geom_edge_arc( aes(edge_colour=weight,edge_width= 150*exp(weight),edge_alpha= (weight)  ) , arrow= arrow(length= unit(10, "mm" ) , type= "closed" ),
                 start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" )  ) +
  geom_node_point(col=rep( rep("red",#c("green","orange","red"),
                               each=nrow(createdLayout)),1)  ,alpha=0.8,  size= 20)+
  geom_node_text(aes(label= name),size=14, repel = T,col= pal_npg()(3)[3] ) + 
  theme_void(base_size= 1.5* 14 ) +
  theme(aspect.ratio = 1, legend.position = "none" ) + 
  scale_edge_width(range=c(0.1,4))+
  scale_edge_alpha(range=c(0.1,0.97)) +
  scale_edge_color_gradient2()
p0 + scale_edge_color_gradient2(low="white",high="black")
#ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Validation Fibroblast activating stimuli TGFB EGF network with textRED.pdf",width=10, height=10)
#ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Validation Fibroblast activating stimuli TGFB EGF network with textRED.jpeg",width=10, height=10)


p0 <- ggraph(createdLayout[] , aes(x= x, y= y , col= Day ) ) + 
  geom_edge_loop(  aes(edge_colour=weight,edge_width= 150*exp(weight),edge_alpha= (weight)) ,  arrow= arrow(length= unit(10, "mm" ) , type= "closed"),
                   start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" ) ) +
  geom_edge_arc( aes(edge_colour=weight,edge_width= 150*exp(weight),edge_alpha= (weight)  ) , arrow= arrow(length= unit(10, "mm" ) , type= "closed" ),
                 start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" )  ) +
  geom_node_point(col=rep( rep("red",
                               each=nrow(createdLayout)),1)  ,alpha=0.8,  size= 20)+
  theme_void(base_size= 1.5* 14 ) +
  theme(aspect.ratio = 1, legend.position = "none" ) + 
  scale_edge_width(range=c(0.1,4))+
  scale_edge_alpha(range=c(0.1,0.97)) +
  scale_edge_color_gradient2()
p0 + scale_edge_color_gradient2(low="white",high="black")
#ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/BLANK Validation Fibroblast activating stimuli TGFB EGF network with textRED.jpeg",width=10, height=10)

validData <- tmp
write.csv(validData,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/ValidationCohortFibroblastactivatingstimuliTGFBEGF.csv")


ppp<-data.table(CCI[LigandPhenoCelltype!="B cells"][ReceptorPhenoCelltype=="Fibroblasts"][Pair.Name%in%l1]%>%
                  group_by(Patient.Study.ID,LigandPhenoCelltype,ReceptorPhenoCelltype,Day,ARM,Treat) %>%
                  dplyr::summarise(muln_scaleTransduction=mean(log(1+scaleTransduction))) )

pppp<- ppp%>%spread(Day,muln_scaleTransduction,fill=NA)
setnames(pppp,old=c("0", "180"),new=c("Day0scores","Day180scores"))
pppp[, supplem:=Day180scores-Day0scores]
pppp$LigandPhenoCelltype <- factor(pppp$LigandPhenoCelltype , levels=c("Fibroblasts","Cancer cells"  ,"Adipocytes","Macrophages", "Normal epithelial cells","Endothelial cells",
                                                                       "CD8+ T cells","Pericytes","CD4+ T cells","B cells" ))
ggplot(pppp[LigandPhenoCelltype!="B cells" ]
       , aes( y=supplem,x=LigandPhenoCelltype, group=interaction(LigandPhenoCelltype)))+
  geom_boxplot()+geom_point(position=position_dodge(width=1))

validStatsData <- ppp[Day==180][LigandPhenoCelltype!="B cells"] 
validStatsData$medianscore <- median(validStatsData$muln_scaleTransduction)
m0 <- lm(muln_scaleTransduction~ -1+ LigandPhenoCelltype ,offset=medianscore,data= validStatsData)
summary(m0)
#write.csv(validStatsData,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/ValidationCohort_Stats_FibroblastactivatingstimuliTGFBEGF.csv")

ggplot(validStatsData
       , aes( y=muln_scaleTransduction,x=LigandPhenoCelltype, group=interaction(LigandPhenoCelltype)))+
  geom_boxplot()+geom_point(position=position_dodge(width=1))



