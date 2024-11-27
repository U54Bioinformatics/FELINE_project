rm(list=ls())
library(glmnet);require(dplyr);require(data.table);require(tidyr)
require(ggraph);require(tidygraph);require(lmerTest)

# Define data location
SourceDataLoc<-"/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript Feline ERBB facilitation/Molecular Systems Biology/Revision and submission folder/Source Data/"
Intermediateloc<-paste0(SourceDataLoc,"Figure 3/")

# Load discovery communication data
Discovery_communication <- data.table(read.csv( file= paste0(Intermediateloc,"SourceData_Figure_3A_DiscoveryTumortoCancerComms.csv")))

Discovery_communication[,Treat:="CombinationRibo"]
Discovery_communication[ARM=="A",Treat:="LetrozoleAlone"]
Discovery_communication[,scaleTransduction:=scale(TransductionMu,center=F),by=c("Pair.Name")] #"ReceptorPhenoCelltype",
Discovery_communication[!is.finite(TransductionMu),scaleTransduction:=0]
Discovery_communication[,scalelnTransductionMu:=scale(log(1+TransductionMu)),by=Pair.Name]
Discovery_communication[,scaleTransduction2:=exp(scale(log(TransductionMu),center=F)),by=c("Pair.Name")] 

# List communication pathway receptors
wchComms<- unique( c( "EGFR","ERBB2","ERBB3","ERBB4"))
wchComms <- c("ACVR1","ACVR1B","ACVR1C","ACVR2A","ACVR2B","ADIPOR1","ALK","EDNRB","EGFR","EPHA1","EPHB3","EPHB4","ERBB2","ERBB3",    
              "ERBB4","FGFR1","FGFR4","FGFRL1","GPC1","GRM1","INSR","LEPR","LRP2","LRP6","LTBR","MET","NRP1","NRXN3",
              "NTRK2","NTRK3","OGFR","PLXNA2","PLXNA4","PLXNB1","PLXNB2","PTPRF","PTPRS","PTPRZ1","ROBO4","SDC1","SDC3","SDC4",     
              "TGFBR2","TNFRSF11B","TNFRSF14","TNFRSF18") 

# Create list of LR pairs with the required receptors
isApath <- sapply( 1:length(wchComms),function(ii){grepl(paste0("_",wchComms[ii]),unique(Discovery_communication$Pair.Name) ) } ) 
l1 <- unique(Discovery_communication$Pair.Name)[rowSums(isApath)>0 ]

# Aggregate signaling between cell types to sample level
average_muln_scaleTransduction <- data.table(Discovery_communication[Pair.Name%in%l1]%>%
                                               group_by(Patient.Study.ID,LigandPhenoCelltype,ReceptorPhenoCelltype,Day,ARM,dynamic_class3,Treat) %>%
                                               dplyr::summarise(muln_scaleTransduction=mean(log2(1+scaleTransduction))) %>% 
                                               group_by(LigandPhenoCelltype, ReceptorPhenoCelltype, Day) %>%#, Treat) %>%
                                               dplyr::summarise(average_muln_scaleTransduction=mean(muln_scaleTransduction)))
average_muln_scaleTransduction$LigandPhenoCelltype <- factor(average_muln_scaleTransduction$LigandPhenoCelltype, 
                                                             levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Normal epithelial cells","Cancer cells","B cells","Macrophages","CD8+ T cells","CD4+ T cells"))
average_muln_scaleTransduction$ReceptorPhenoCelltype <- factor(average_muln_scaleTransduction$ReceptorPhenoCelltype, 
                                                               levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Normal epithelial cells","Cancer cells","B cells","Macrophages","CD8+ T cells","CD4+ T cells") )

# Baseline signling
initState <- data.table(average_muln_scaleTransduction[Day==0]%>%
                          group_by(LigandPhenoCelltype,ReceptorPhenoCelltype)%>% #Treat, dynamic_class3
                          dplyr::summarise(intiState=mean(average_muln_scaleTransduction)) )

# Fold changes in communication
average_muln_scaleTransduction2 <- merge(average_muln_scaleTransduction,initState,by=c("LigandPhenoCelltype", "ReceptorPhenoCelltype")) #,"Treat","dynamic_class3"
average_muln_scaleTransduction2[ , lnfoldchange:= (average_muln_scaleTransduction-intiState) ]      #average_muln_scaleTransduction2[ , lnfoldchange:= (average_muln_scaleTransduction) ] 

# Plotting data roughly
ppp<-data.table(Discovery_communication[LigandPhenoCelltype!="B cells"][ReceptorPhenoCelltype=="Cancer cells"][Day%in%c(0,180)][Pair.Name%in%l1]%>%
                  group_by(Patient.Study.ID,LigandPhenoCelltype,ReceptorPhenoCelltype,Day,ARM,Treat) %>%
                  dplyr::summarise(muln_scaleTransduction=mean(log(1+scaleTransduction))))

pppp<- data.table(ppp%>%spread(Day,muln_scaleTransduction,fill=NA))
setnames(pppp,old=c("0", "180"),new=c("Day0scores","Day180scores"))
pppp[, supplem:=Day180scores-Day0scores]
pppp$LigandPhenoCelltype <- factor(pppp$LigandPhenoCelltype , levels=c("Fibroblasts","Cancer cells"  ,"Adipocytes","Macrophages", "Normal epithelial cells","Endothelial cells",
                                                                       "CD8+ T cells","Pericytes","CD4+ T cells","B cells" ))
ggplot(pppp[LigandPhenoCelltype!="B cells" ]#[Treat!="CombinationRibo"]
       , aes( y=supplem,x=LigandPhenoCelltype, group=interaction(LigandPhenoCelltype)))+
  geom_boxplot()+geom_point(position=position_dodge(width=1))#+facet_wrap(~Treat)

#  stats
discStatsData <- ppp[][LigandPhenoCelltype!="B cells"] 
m0 <- lm(muln_scaleTransduction~ LigandPhenoCelltype+LigandPhenoCelltype:Day ,data= discStatsData)
summary(m0)
summary(lm(supplem~-1+ LigandPhenoCelltype ,data= pppp[ReceptorPhenoCelltype=="Cancer cells"]))
#write.csv(discStatsData,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/DiscoveryCohort_Stats_SupplementaryGFsignalstocancercellsnetwork.csv")

# construct weighted graph
CCIcondensed <- average_muln_scaleTransduction2[LigandPhenoCelltype!="B cells"][ReceptorPhenoCelltype!="B cells"]  
CCIcondensed$LigandPhenoCelltype <-factor(CCIcondensed$LigandPhenoCelltype, 
                                          levels= sort(c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Normal epithelial cells","Cancer cells","Macrophages","CD8+ T cells","CD4+ T cells"))) #"B cells",
CCIcondensed$ReceptorPhenoCelltype <-factor(CCIcondensed$ReceptorPhenoCelltype, 
                                            levels= sort(c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Normal epithelial cells","Cancer cells","Macrophages","CD8+ T cells","CD4+ T cells") )) #"B cells",
setnames(CCIcondensed, old=c("LigandPhenoCelltype", "ReceptorPhenoCelltype", "lnfoldchange"), new= c("from", "to", "weight"))

plotthis<-CCIcondensed[Day==180][to=="Cancer cells"]
plotthis$from <- factor(plotthis$from,
                        levels=c("Fibroblasts","Adipocytes","Macrophages","Endothelial cells","Normal epithelial cells","Pericytes","Cancer cells","CD4+ T cells","CD8+ T cells"))
ggplot(plotthis, 
       aes( y= weight,x= from, group= interaction(from,Day) )) + theme_classic(base_size= 26)+
  theme(aspect.ratio = 0.4) +
  geom_boxplot(fill="red") + geom_point(position=position_dodge(width=1))+
  labs(y="Supplementary growth factor signaling \n (log fold change pre vs post treatment)", x="Signal sending cell type")

# visualize communication network
tmp <- data.table(CCIcondensed[][Day==180][to=="Cancer cells"]%>%group_by(from, to)%>%dplyr::summarise(weight= mean(weight),"Treat"="all"))[order( Treat,from,to)]
tmp$DayLab <- paste0("Day ", tmp$Day )
tmp <- data.table( tmp %>% group_by(Treat) %>% mutate(weight = scale(weight) ) )
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
  #geom_node_point( col= rep( rep(node_col, each= nrow(createdLayout) ) , each= 2) , size= 20 ) + 
  geom_node_point(col=rep( rep("red",#c("green","orange","red"),
                               each=nrow(createdLayout)),1)  ,alpha=0.8,  size= 20)+#geom_node_point( col= rep( rep(node_col, each= nrow(createdLayout) ) ,  2) , size= 20 ) + 
  geom_node_text(aes(label= name),size=14, repel = T,col= pal_npg()(3)[3] ) + #,col= "grey80" ) + #, size= 9) +
  #facet_edges(~DayLab  ) + 
  theme_void(base_size= 1.5* 14 ) +
  theme(aspect.ratio = 1, legend.position = "none" ) + 
  scale_edge_width(range=c(0.1,4))+
  scale_edge_alpha(range=c(0.1,0.97)) +
  scale_edge_color_gradient2()
#scale_edge_color_viridis(option="B")
p0 + scale_edge_color_gradient2(low="white",high="black")
#ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Discovery Supplementary GF signals to cancer cells network with textRED.pdf",width=10, height=10)
#ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Discovery Supplementary GF signals to cancer cells network with textRED.jpeg",width=10, height=10)

p0 <- ggraph(createdLayout[] , aes(x= x, y= y , col= Day ) ) + 
  geom_edge_loop(  aes(edge_colour=weight,edge_width= 150*exp(weight),edge_alpha= (weight)) ,  arrow= arrow(length= unit(10, "mm" ) , type= "closed"),
                   start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" ) ) +
  geom_edge_arc( aes(edge_colour=weight,edge_width= 150*exp(weight),edge_alpha= (weight)  ) , arrow= arrow(length= unit(10, "mm" ) , type= "closed" ),
                 start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" )  ) +
  #geom_node_point( col= rep( rep(node_col, each= nrow(createdLayout) ) , each= 2) , size= 20 ) + 
  geom_node_point(col=rep( rep("red",#c("green","orange","red"),
                               each=nrow(createdLayout)),1)  ,alpha=0.8,  size= 20)+#geom_node_point( col= rep( rep(node_col, each= nrow(createdLayout) ) ,  2) , size= 20 ) + 
  # geom_node_text(aes(label= name),size=14, repel = T,col= pal_npg()(3)[3] ) + #,col= "grey80" ) + #, size= 9) +
  #facet_edges(~DayLab  ) + 
  theme_void(base_size= 1.5* 14 ) +
  theme(aspect.ratio = 1, legend.position = "none" ) + 
  scale_edge_width(range=c(0.1,4))+
  scale_edge_alpha(range=c(0.1,0.97)) +
  scale_edge_color_gradient2()
#scale_edge_color_viridis(option="B")
p0 + scale_edge_color_gradient2(low="white",high="black")
#ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/BLANK Discovery Supplementary GF signals to cancer cells network with textRED.jpeg",width=10, height=10)

discData<-tmp
#write.csv(discData,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/DiscoveryCohortSupplementaryGFsignalstocancer.csv")








## Replicate analysis in validation cohort

# load validation cohort data
Validation_communication <- data.table(read.csv( file= paste0(Intermediateloc,"SourceData_Figure_3A_ValidationTumortoCancerComms.csv")))
Validation_communication[,Treat:="CombinationRibo"]
Validation_communication[ARM=="A",Treat:="LetrozoleAlone"]
Validation_communication[,scaleTransduction:=scale(TransductionMu,center=F),by=c("Pair.Name")] #"ReceptorPhenoCelltype",
Validation_communication[!is.finite(TransductionMu),scaleTransduction:=0]
Validation_communication[,scalelnTransductionMu:=scale(log(1+TransductionMu)),by=Pair.Name]
Validation_communication[,scaleTransduction2:=exp(scale(log(TransductionMu),center=F)),by=c("Pair.Name")] 

#wchComms<- unique( c( "EGFR","ERBB2","ERBB3","ERBB4"))
#wchComms <- c("ACVR1","ACVR1B","ACVR1C","ACVR2A","ACVR2B","ADIPOR1","ALK","EDNRB","EGFR","EPHA1","EPHB3","EPHB4","ERBB2","ERBB3",    
#  "ERBB4","FGFR1","FGFR4","FGFRL1","GPC1","GRM1","INSR","LEPR","LRP2","LRP6","LTBR","MET","NRP1","NRXN3",
#  "NTRK2","NTRK3","OGFR","PLXNA2","PLXNA4","PLXNB1","PLXNB2","PTPRF","PTPRS","PTPRZ1","ROBO4","SDC1","SDC3","SDC4",     
#  "TGFBR2","TNFRSF11B","TNFRSF14","TNFRSF18") 

isApath <- sapply( 1:length(wchComms),function(ii){grepl(paste0("_",wchComms[ii]),unique(Validation_communication$Pair.Name) ) } ) 
l1 <- unique(Validation_communication$Pair.Name)[rowSums(isApath)>0 ]

average_muln_scaleTransduction <- data.table(Validation_communication[Pair.Name%in%l1]%>%
                                               group_by(Patient.Study.ID,LigandPhenoCelltype,ReceptorPhenoCelltype,Day,ARM) %>%#,Treat) %>%
                                               dplyr::summarise(muln_scaleTransduction=mean(log2(1+scaleTransduction))) %>% 
                                               group_by(LigandPhenoCelltype, ReceptorPhenoCelltype, Day) %>%#, Treat) %>%
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


ppp<-data.table(Validation_communication[LigandPhenoCelltype!="B cells"][ReceptorPhenoCelltype=="Cancer cells"][Day%in%c(0,180)][Pair.Name%in%l1]%>%
                  group_by(Patient.Study.ID,LigandPhenoCelltype,ReceptorPhenoCelltype,Day,ARM,Treat) %>%
                  dplyr::summarise(muln_scaleTransduction=mean(log(1+scaleTransduction))))

pppp<- data.table(ppp%>%spread(Day,muln_scaleTransduction,fill=NA))
setnames(pppp,old=c("0", "180"),new=c("Day0scores","Day180scores"))
pppp[, supplem:=Day180scores-Day0scores]
pppp$LigandPhenoCelltype <- factor(pppp$LigandPhenoCelltype , levels=c("Fibroblasts","Cancer cells"  ,"Adipocytes","Macrophages", "Normal epithelial cells","Endothelial cells",
                                                                       "CD8+ T cells","Pericytes","CD4+ T cells","B cells" ))
ggplot(pppp[LigandPhenoCelltype!="B cells" ]#[Treat!="CombinationRibo"]
       , aes( y=supplem,x=LigandPhenoCelltype, group=interaction(LigandPhenoCelltype)))+
  geom_boxplot()+geom_point(position=position_dodge(width=1))#+facet_wrap(~Treat)

validStatsData <- ppp[][LigandPhenoCelltype!="B cells"] 
m0 <- lm(muln_scaleTransduction~ LigandPhenoCelltype+LigandPhenoCelltype:Day ,data= validStatsData)
summary(m0)

summary(lm(supplem~-1+ LigandPhenoCelltype ,data= pppp[ReceptorPhenoCelltype=="Cancer cells"]))
#write.csv(validStatsData,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/ValidationCohort_Stats_SupplementaryGFsignalstocancercellsnetwork.csv")

CCIcondensed <- average_muln_scaleTransduction2[LigandPhenoCelltype!="B cells"][ReceptorPhenoCelltype!="B cells"]  
CCIcondensed$LigandPhenoCelltype <-factor(CCIcondensed$LigandPhenoCelltype, 
                                          levels= sort(c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Normal epithelial cells","Cancer cells","Macrophages","CD8+ T cells","CD4+ T cells"))) #"B cells",
CCIcondensed$ReceptorPhenoCelltype <-factor(CCIcondensed$ReceptorPhenoCelltype, 
                                            levels= sort(c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Normal epithelial cells","Cancer cells","Macrophages","CD8+ T cells","CD4+ T cells") )) #"B cells",
setnames(CCIcondensed, old=c("LigandPhenoCelltype", "ReceptorPhenoCelltype", "lnfoldchange"), new= c("from", "to", "weight"))

plotthis<-CCIcondensed[Day==180][to=="Cancer cells"]
plotthis$from <- factor(plotthis$from,
                        levels=c("Fibroblasts","Adipocytes","Macrophages","Endothelial cells","Normal epithelial cells","Pericytes","Cancer cells","CD4+ T cells","CD8+ T cells"))
ggplot(plotthis,
       aes( y= weight,x= from, group= interaction(from,Day) )) + theme_classic(base_size= 26)+
  theme(aspect.ratio = 0.4) +
  geom_boxplot(fill="red") + geom_point(position=position_dodge(width=1))+
  labs(y="Supplementary growth factor signaling \n (log fold change pre vs post treatment)", x="Signal sending cell type")

tmp <- data.table(CCIcondensed[][Day==180][to=="Cancer cells"]%>%group_by(from, to)%>%dplyr::summarise(weight= (mean(weight)),"Treat"="all"))[order( Treat,from,to)]
tmp$DayLab <- paste0("Day ", tmp$Day )
tmp <- data.table( tmp %>% group_by(Treat) %>% mutate(weight = (weight) ) )
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
  #geom_node_point( col= rep( rep(node_col, each= nrow(createdLayout) ) , each= 2) , size= 20 ) + 
  geom_node_point(col=rep( rep("royalblue",#c("green","orange","red"),
                               each=nrow(createdLayout)),1)  ,alpha=0.5,  size= 20)+#geom_node_point( col= rep( rep(node_col, each= nrow(createdLayout) ) ,  2) , size= 20 ) + 
  #geom_node_text(aes(label= name),size=10 ) + #,col= "grey80" ) + #, size= 9) +
  #facet_edges(~DayLab  ) + 
  theme_void(base_size= 1.5* 14 ) +
  theme(aspect.ratio = 1, legend.position = "none" ) + 
  scale_edge_width(range=c(0.1,4))+
  scale_edge_alpha(range=c(0.1,0.97)) +
  scale_edge_color_gradient2()
#scale_edge_color_viridis(option="B")
p0 + scale_edge_color_gradient2(low="white",high="black")
#ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Validation cohort Supplementary GF signals to cancer cells network.pdf",width=12, height=12)




p0 <- ggraph(createdLayout[] , aes(x= x, y= y , col= Day ) ) + 
  geom_edge_loop(  aes(edge_colour=weight,edge_width= 150*exp(weight),edge_alpha= (weight)) ,  arrow= arrow(length= unit(10, "mm" ) , type= "closed"),
                   start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" ) ) +
  geom_edge_arc( aes(edge_colour=weight,edge_width= 150*exp(weight),edge_alpha= (weight)  ) , arrow= arrow(length= unit(10, "mm" ) , type= "closed" ),
                 start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" )  ) +
  #geom_node_point( col= rep( rep(node_col, each= nrow(createdLayout) ) , each= 2) , size= 20 ) + 
  geom_node_point(col=rep( rep("red",#c("green","orange","red"),
                               each=nrow(createdLayout)),1)  ,alpha=0.8,  size= 20)+#geom_node_point( col= rep( rep(node_col, each= nrow(createdLayout) ) ,  2) , size= 20 ) + 
  geom_node_text(aes(label= name),size=14, repel = T,col= pal_npg()(3)[3] ) + #,col= "grey80" ) + #, size= 9) +
  #facet_edges(~DayLab  ) + 
  theme_void(base_size= 1.5* 14 ) +
  theme(aspect.ratio = 1, legend.position = "none" ) + 
  scale_edge_width(range=c(0.1,4))+
  scale_edge_alpha(range=c(0.1,0.97)) +
  scale_edge_color_gradient2()
#scale_edge_color_viridis(option="B")
p0 + scale_edge_color_gradient2(low="white",high="black")
#ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Validation Supplementary GF signals to cancer cells network with textRED.pdf",width=10, height=10)
#ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Validation Supplementary GF signals to cancer cells network with textRED.jpeg",width=10, height=10)

p0 <- ggraph(createdLayout[] , aes(x= x, y= y , col= Day ) ) + 
  geom_edge_loop(  aes(edge_colour=weight,edge_width= 150*exp(weight),edge_alpha= (weight)) ,  arrow= arrow(length= unit(10, "mm" ) , type= "closed"),
                   start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" ) ) +
  geom_edge_arc( aes(edge_colour=weight,edge_width= 150*exp(weight),edge_alpha= (weight)  ) , arrow= arrow(length= unit(10, "mm" ) , type= "closed" ),
                 start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" )  ) +
  #geom_node_point( col= rep( rep(node_col, each= nrow(createdLayout) ) , each= 2) , size= 20 ) + 
  geom_node_point(col=rep( rep("red",#c("green","orange","red"),
                               each=nrow(createdLayout)),1)  ,alpha=0.8,  size= 20)+#geom_node_point( col= rep( rep(node_col, each= nrow(createdLayout) ) ,  2) , size= 20 ) + 
  #geom_node_text(aes(label= name),size=14, repel = T,col= pal_npg()(3)[3] ) + #,col= "grey80" ) + #, size= 9) +
  #facet_edges(~DayLab  ) + 
  theme_void(base_size= 1.5* 14 ) +
  theme(aspect.ratio = 1, legend.position = "none" ) + 
  scale_edge_width(range=c(0.1,4))+
  scale_edge_alpha(range=c(0.1,0.97)) +
  scale_edge_color_gradient2()
#scale_edge_color_viridis(option="B")
p0 + scale_edge_color_gradient2(low="white",high="black")
ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/BLANK Validation Supplementary GF signals to cancer cells network with textRED.jpeg",width=10, height=10)


validData<-tmp
write.csv(validData,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/ValidationCohortSupplementaryGFsignalstocancer.csv")













