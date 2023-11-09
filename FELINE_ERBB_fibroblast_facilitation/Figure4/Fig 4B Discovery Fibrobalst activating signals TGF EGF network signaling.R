rm(list=ls())
library(glmnet);require(dplyr);require(data.table);require(tidyr)
require(ggraph);require(tidygraph)

load( "~/Dropbox/Cancer_pheno_evo/data/FELINE2/LigandReceptor/Filtered_Human-2015-Ramilowski-LR-pairs.RData" )
load(file ="~/Dropbox/Cancer_pheno_evo/data/FELINE2/Communication output merged ALL ARMS/PopulationCommunicationMerged ALL ARMS.RData" )#CCI,allphenotypes, uu,perIndiv,
LRgenelist <- unique( data.table(merge(CCI[],LRpairsFiltered,by="Pair.Name"))$HPMR.Receptor) #unique( c(LRpairsFiltered$HPMR.Receptor) ) #, LRpairsFiltered$HPMR.Ligand) )
rm(list= c("LRpairs", "LRpairsFiltered"))

CCI[,Treat:="CombinationRibo"]
CCI[ARM=="A",Treat:="LetrozoleAlone"]
CCI[,scaleTransduction:=scale(TransductionMu,center=F),by=c("Pair.Name")] #"ReceptorPhenoCelltype",
CCI[!is.finite(TransductionMu),scaleTransduction:=0]
CCI[,scalelnTransductionMu:=scale(log(1+TransductionMu)),by=Pair.Name]
CCI[,scaleTransduction2:=exp(scale(log(TransductionMu),center=F)),by=c("Pair.Name")] 

load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/GeneOntology/growthFactorReceptors.RData" )

wchComms<- unique( c( "TGFBR1","TGFBR2","TGFBR3"))

isApath<-sapply( 1:length(wchComms),function(ii){grepl(paste0("_",wchComms[ii]),unique(CCI$Pair.Name) ) } ) 
l1 <- unique(CCI$Pair.Name)[rowSums(isApath)>0 ]
plotx<-CCI[Pair.Name%in%l1][Day!=14][][LigandPhenoCelltype!="B cells"][ReceptorPhenoCelltype=="Fibroblasts"] 
plotx$LigandPhenoCelltype <- factor(plotx$LigandPhenoCelltype , levels=c("Cancer cells" ,"Fibroblasts","Normal epithelial cells","CD4+ T cells","CD8+ T cells","Macrophages","Endothelial cells","Pericytes",
                                                                         "Adipocytes") )
ggplot( plotx ,
        aes(x=LigandPhenoCelltype,fill=as.factor(Day),col=as.factor(Day),y=log(scaleTransduction)  ,group=interaction(LigandPhenoCelltype,Day) )) + 
  theme_classic(base_size=22)+
  geom_boxplot(alpha=0.6,position=position_dodge(width=1))+
  geom_point(position=position_dodge(width=1))+labs(y= "Fibroblast receipt of \n TGF Beta communications",x="Signal sendercell type")+scale_y_continuous(breaks=log(c(10,1,0.1,0.01,0.001)),labels=c(10,1,0.1,0.01,0.001))+
  theme(aspect.ratio=0.5, axis.text.x = element_text(angle = 90))+
  scale_color_manual(values= rev(ggsci::pal_npg("nrc")(2)  ),name="Time point", labels=c("Pre treatment","Post treatment"))+
  scale_fill_manual(values= rev(ggsci::pal_npg("nrc")(2)  ),name="Time point", labels=c("Pre treatment","Post treatment"))
ggsave("~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Fibroblast TGFBR communication start vs end.png",height=10,width=10)


wchComms<- unique( c( "HBEGFR","EGFR"))

isApath<-sapply( 1:length(wchComms),function(ii){grepl(paste0("_",wchComms[ii]),unique(CCI$Pair.Name) ) } ) 
l1 <- unique(CCI$Pair.Name)[rowSums(isApath)>0 ]
plotx<-CCI[Pair.Name%in%l1][Day!=14][][LigandPhenoCelltype!="B cells"][ReceptorPhenoCelltype=="Fibroblasts"] 
plotx$LigandPhenoCelltype <- factor(plotx$LigandPhenoCelltype , levels=c("Cancer cells" ,"Fibroblasts","Normal epithelial cells","CD4+ T cells","CD8+ T cells","Macrophages","Endothelial cells","Pericytes",
                                                                         "Adipocytes") )
ggplot( plotx ,
        aes(x=LigandPhenoCelltype,fill=as.factor(Day),col=as.factor(Day),y=log(scaleTransduction)  ,group=interaction(LigandPhenoCelltype,Day) )) + 
  theme_classic(base_size=22)+
  geom_boxplot(alpha=0.6,position=position_dodge(width=1))+
  geom_point(position=position_dodge(width=1))+labs(y= "Fibroblast receipt of \n EGF communications",x="Signal sendercell type")+scale_y_continuous(breaks=log(c(10,1,0.1,0.01,0.001)),labels=c(10,1,0.1,0.01,0.001))+
  theme(aspect.ratio=0.5, axis.text.x = element_text(angle = 90))+
  scale_color_manual(values= rev(ggsci::pal_npg("nrc")(2)  ),name="Time point", labels=c("Pre treatment","Post treatment"))+
  scale_fill_manual(values= rev(ggsci::pal_npg("nrc")(2)  ),name="Time point", labels=c("Pre treatment","Post treatment"))
ggsave("~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Fibroblast EGFR communication start vs end.png",height=10,width=10)



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
                          group_by(LigandPhenoCelltype,ReceptorPhenoCelltype)%>% 
                          dplyr::summarise(intiState=mean(average_muln_scaleTransduction)) )

average_muln_scaleTransduction2 <- merge(average_muln_scaleTransduction,initState,by=c("LigandPhenoCelltype", "ReceptorPhenoCelltype")) #,"Treat","dynamic_class3"
average_muln_scaleTransduction2[ , lnfoldchange:= (average_muln_scaleTransduction-intiState) ]      #average_muln_scaleTransduction2[ , lnfoldchange:= (average_muln_scaleTransduction) ] 

ggplot( average_muln_scaleTransduction2[][LigandPhenoCelltype!="B cells"][ReceptorPhenoCelltype!="B cells"] ,
        aes(LigandPhenoCelltype,ReceptorPhenoCelltype,fill=average_muln_scaleTransduction-intiState  )) + geom_tile()+
  facet_wrap(Treat~Day)+scale_fill_viridis_c(name="Log fold change \n in communication \n relative to baseline",option="B") + theme_classic()+
  labs(y="Signal receiver",x="Signal sender")+theme(aspect.ratio=1,axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot( average_muln_scaleTransduction2[][][LigandPhenoCelltype!="B cells"][ReceptorPhenoCelltype=="Fibroblasts"] ,
        aes(y=LigandPhenoCelltype,x=as.factor(Day),fill=average_muln_scaleTransduction  )) + geom_tile()+
  scale_fill_viridis_c(name="Log fold change \n in communication \n relative to baseline",option="B") + theme_classic()+
  labs(y="Signal receiver",x="Day")+theme(aspect.ratio=1,axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



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
  geom_edge_loop(  aes(edge_colour=weight,edge_width= 150*(weight),edge_alpha= (weight)) ,  arrow= arrow(length= unit(20, "mm" ) , type= "closed"),
                   start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" ) ) +
  geom_edge_arc( aes(edge_colour=weight,edge_width= 150*(weight),edge_alpha= (weight)  ) , arrow= arrow(length= unit(20, "mm" ) , type= "closed" ),
                 start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" )  ) +
  geom_node_point(col=rep( rep("forestgreen",
                               each=nrow(createdLayout)),1)  ,alpha=0.5,  size= 20)+#geom_node_point( col= rep( rep(node_col, each= nrow(createdLayout) ) ,  2) , size= 20 ) + 
  geom_node_text(aes(label= name),size=12 ) + 
  theme_void(base_size= 1.5* 14 ) +
  theme(aspect.ratio = 1, legend.position = "none" ) + 
  scale_edge_width(range=c(0.1,4))+
  scale_edge_alpha(range=c(0.1,0.97)) +
  scale_edge_color_gradient2()
p0 + scale_edge_color_gradient2(low="white",high="black")
ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Globe plot ggrraph/Fibroblast Growth Factor Communication globe plot all times.png", height=10, width = 10)


p0 <- ggraph(createdLayout[] , aes(x= x, y= y , col= Day ) ) + 
  geom_edge_loop(  aes(edge_colour=weight,edge_width= 150*(weight),edge_alpha= (weight)) ,  arrow= arrow(length= unit(20, "mm" ) , type= "open"),
                   start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" ) ) +
  geom_edge_arc( aes(edge_colour=weight,edge_width= 150*(weight),edge_alpha= (weight)  ) , arrow= arrow(length= unit(20, "mm" ) , type= "open" ),
                 start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" )  ) +
  geom_node_point(col=rep( rep("forestgreen",
                               each=nrow(createdLayout)),1)  ,alpha=0.5,  size= 20)+
  theme_void(base_size= 1.5* 14 ) +
  theme(aspect.ratio = 1, legend.position = "none" ) + 
  scale_edge_width(range=c(0.1,4))+
  scale_edge_alpha(range=c(0.1,0.97)) +
  scale_edge_color_gradient2()
p0 + scale_edge_color_gradient2(low="white",high="black")
#ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Globe plot ggrraph/BLANK Fibroblast Growth Factor Communication globe plot all times.png", height=10, width = 10)

ggsave( file="~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/BLANK Fibroblast Growth Factor Communication globe plot all times.png", height=10, width = 10)


p0 <- ggraph(createdLayout[] , aes(x= x, y= y , col= Day ) ) + 
  geom_edge_loop(  aes(edge_colour=weight,edge_width= 150*exp(weight),edge_alpha= (weight)) ,  arrow= arrow(length= unit(6, "mm" ) , type= "open"),
                   start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" ) ) +
  geom_edge_arc( aes(edge_colour=weight,edge_width= 150*exp(weight),edge_alpha= (weight)  ) , arrow= arrow(length= unit(6, "mm" ) , type= "open" ),
                 start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" )  ) +
  geom_node_point(col=rep( rep("royalblue",
                               each=nrow(createdLayout)),1)  ,alpha=0.5,  size= 20)+
  theme_void(base_size= 1.5* 14 ) +
  theme(aspect.ratio = 1, legend.position = "none" ) + 
  scale_edge_width(range=c(0.1,4))+
  scale_edge_alpha(range=c(0.1,0.97)) +
  scale_edge_color_gradient2()
p0 + scale_edge_color_gradient2(low="white",high="black")
ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Globe plot ggrraph/Letrozole Growth Factor Communication to cancer cells globe plot Day 180 relative to baseline.png", height=10, width = 10)
ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Globe plot ggrraph/Letrozole Growth Factor Communication globe plot Day 180 relative to baseline.png", height=10, width = 10)



p0 <- ggraph(createdLayout[] , aes(x= x, y= y , col= Day ) ) + 
  geom_edge_loop(  aes(edge_colour=weight,edge_width= 150*exp(weight),edge_alpha= (weight)) ,  arrow= arrow(length= unit(10, "mm" ) , type= "closed"),
                   start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" ) ) +
  geom_edge_arc( aes(edge_colour=weight,edge_width= 150*exp(weight),edge_alpha= (weight)  ) , arrow= arrow(length= unit(10, "mm" ) , type= "closed" ),
                 start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" )  ) +
  geom_node_point(col=rep( rep("red",
                               each=nrow(createdLayout)),1)  ,alpha=0.8,  size= 20)+
  geom_node_text(aes(label= name),size=14, repel = T,col= pal_npg()(3)[3] ) + 
  theme_void(base_size= 1.5* 14 ) +
  theme(aspect.ratio = 1, legend.position = "none" ) + 
  scale_edge_width(range=c(0.1,4))+
  scale_edge_alpha(range=c(0.1,0.97)) +
  scale_edge_color_gradient2()
p0 + scale_edge_color_gradient2(low="white",high="black")
ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Discovery Fibroblast activating stimuli TGFB EGF network with textRED.pdf",width=10, height=10)
ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Discovery Fibroblast activating stimuli TGFB EGF network with textRED.jpeg",width=10, height=10)


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
ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/BLANK Discovery Fibroblast activating stimuli TGFB EGF network with textRED.jpeg",width=10, height=10)


discData<-tmp
write.csv(discData,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/DiscoveryCohortFibroblastactivatingstimuliTGFBEGF.csv")


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

discStatsData <- ppp[Day==180][LigandPhenoCelltype!="B cells"] 
discStatsData$medianscore <- median(discStatsData$muln_scaleTransduction)
m0 <- lm(muln_scaleTransduction~ -1+ LigandPhenoCelltype ,offset=medianscore,data= discStatsData)
summary(m0)
#write.csv(discStatsData,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/DiscoveryCohort_Stats_FibroblastactivatingstimuliTGFBEGF.csv")

ggplot(discStatsData
       , aes( y=muln_scaleTransduction,x=LigandPhenoCelltype, group=interaction(LigandPhenoCelltype)))+
  geom_boxplot()+geom_point(position=position_dodge(width=1))



