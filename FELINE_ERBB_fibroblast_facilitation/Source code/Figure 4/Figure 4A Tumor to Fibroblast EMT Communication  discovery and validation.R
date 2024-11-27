rm(list=ls())
library(glmnet);require(dplyr);require(data.table);require(tidyr)
require(ggraph);require(tidygraph)

# Define data location
SourceDataLoc <- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript Feline ERBB facilitation/Molecular Systems Biology/Revision and submission folder/Source Data/"
Intermediateloc <- paste0(SourceDataLoc,"Figure 4/")

discovery_LRcomms <- data.table( read.csv(file= paste0(Intermediateloc, "SourceData_Figure_4AB_DiscoveryTumorComms.csv" ) ) ) 

# Treatment labels
discovery_LRcomms[,Treat:="CombinationRibo"]
discovery_LRcomms[ARM=="A",Treat:="LetrozoleAlone"]
# Scaling and transformation of signals
discovery_LRcomms[,scaleTransduction:=scale(TransductionMu,center=F),by=c("Pair.Name")] #"ReceptorPhenoCelltype",
discovery_LRcomms[!is.finite(TransductionMu),scaleTransduction:=0]
discovery_LRcomms[,scalelnTransductionMu:=scale(log(1+TransductionMu)),by=Pair.Name]

# Specify receptors of interest: TGFB1-3
wchComms<- unique( c( "TGFBR1","TGFBR2","TGFBR3"))#,"IGF1R","AXL" )

# identify relevant LR pairs
isApath<-sapply( 1:length(wchComms),function(ii){grepl(paste0("_",wchComms[ii]),unique(discovery_LRcomms$Pair.Name) ) } ) 
l1 <- unique(discovery_LRcomms$Pair.Name)[rowSums(isApath)>0 ]


# Raw communication plot
plotx<-discovery_LRcomms[Pair.Name%in%l1][Day!=14][][LigandPhenoCelltype!="B cells"][ReceptorPhenoCelltype=="Fibroblasts"] 
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

# Specify receptors of interest: (HB)EGFR
wchComms<- unique( c( "HBEGFR","EGFR"))#,"IGF1R","AXL" )

# identify relevant LR pairs
isApath<-sapply( 1:length(wchComms),function(ii){grepl(paste0("_",wchComms[ii]),unique(discovery_LRcomms$Pair.Name) ) } ) 
l1 <- unique(discovery_LRcomms$Pair.Name)[rowSums(isApath)>0 ]


plotx<-discovery_LRcomms[Pair.Name%in%l1][Day!=14][][LigandPhenoCelltype!="B cells"][ReceptorPhenoCelltype=="Fibroblasts"] 
plotx$LigandPhenoCelltype <- factor(plotx$LigandPhenoCelltype , levels=c("Cancer cells" ,"Fibroblasts","Normal epithelial cells","CD4+ T cells","CD8+ T cells","Macrophages","Endothelial cells","Pericytes","Adipocytes") )
ggplot( plotx ,
        aes(x=LigandPhenoCelltype,fill=as.factor(Day),col=as.factor(Day),y=log(scaleTransduction)  ,group=interaction(LigandPhenoCelltype,Day) )) + 
  theme_classic(base_size=22)+
  geom_boxplot(alpha=0.6,position=position_dodge(width=1))+
  geom_point(position=position_dodge(width=1))+labs(y= "Fibroblast receipt of \n EGF communications",x="Signal sendercell type")+scale_y_continuous(breaks=log(c(10,1,0.1,0.01,0.001)),labels=c(10,1,0.1,0.01,0.001))+
  theme(aspect.ratio=0.5, axis.text.x = element_text(angle = 90))+
  scale_color_manual(values= rev(ggsci::pal_npg("nrc")(2)  ),name="Time point", labels=c("Pre treatment","Post treatment"))+
  scale_fill_manual(values= rev(ggsci::pal_npg("nrc")(2)  ),name="Time point", labels=c("Pre treatment","Post treatment"))
#ggsave("~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Fibroblast EGFR communication start vs end.png",height=10,width=10)


# Combined analysis
wchComms<- unique( c( "EGFR","HBEGFR","TGFBR1","TGFBR2","TGFBR3"))#,"IGF1R","AXL" )

isApath<-sapply( 1:length(wchComms),function(ii){grepl(paste0("_",wchComms[ii]),unique(discovery_LRcomms$Pair.Name) ) } ) 
l1 <- unique(discovery_LRcomms$Pair.Name)[rowSums(isApath)>0 ]

# summarise data by signaling cell type within tumor sample (patient and day)
average_muln_scaleTransduction <- data.table(discovery_LRcomms[Pair.Name%in%l1]%>%
                                               group_by(Patient.Study.ID,LigandPhenoCelltype,ReceptorPhenoCelltype,Day,ARM,Treat) %>%
                                               dplyr::summarise(muln_scaleTransduction=mean(log(1+scaleTransduction))) %>% 
                                               group_by(LigandPhenoCelltype, ReceptorPhenoCelltype, Day, Treat) %>%
                                               dplyr::summarise(average_muln_scaleTransduction=mean(muln_scaleTransduction)))
average_muln_scaleTransduction$LigandPhenoCelltype <- factor(average_muln_scaleTransduction$LigandPhenoCelltype, 
                                                             levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Normal epithelial cells","Cancer cells","B cells","Macrophages","CD8+ T cells","CD4+ T cells"))
average_muln_scaleTransduction$ReceptorPhenoCelltype <- factor(average_muln_scaleTransduction$ReceptorPhenoCelltype, 
                                                               levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Normal epithelial cells","Cancer cells","B cells","Macrophages","CD8+ T cells","CD4+ T cells") )

# Calc scores relative to patient baseline
initState <- data.table(average_muln_scaleTransduction[Day==0]%>%
                          group_by(LigandPhenoCelltype,ReceptorPhenoCelltype)%>% #Treat, dynamic_class3
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


# Construct communication network
discovery_LRcommscondensed <- average_muln_scaleTransduction2[LigandPhenoCelltype!="B cells"][ReceptorPhenoCelltype!="B cells"]  
discovery_LRcommscondensed$LigandPhenoCelltype <-factor(discovery_LRcommscondensed$LigandPhenoCelltype, 
                                          levels= sort(c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Normal epithelial cells","Cancer cells","Macrophages","CD8+ T cells","CD4+ T cells"))) #"B cells",
discovery_LRcommscondensed$ReceptorPhenoCelltype <-factor(discovery_LRcommscondensed$ReceptorPhenoCelltype, 
                                            levels= sort(c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Normal epithelial cells","Cancer cells","Macrophages","CD8+ T cells","CD4+ T cells") )) #"B cells",
setnames(discovery_LRcommscondensed, old=c("LigandPhenoCelltype", "ReceptorPhenoCelltype", "average_muln_scaleTransduction"), new= c("from", "to", "weight"))
#setnames(discovery_LRcommscondensed, old=c("LigandPhenoCelltype", "ReceptorPhenoCelltype", "lnfoldchange"), new= c("from", "to", "weight"))

tmp <- discovery_LRcommscondensed[][][to=="Fibroblasts"][order( Treat,from,to)]
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
n <- length(unique(discovery_LRcommscondensed$to)) -1
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
  #geom_node_point( col= rep( rep(node_col, each= nrow(createdLayout) ) , each= 2) , size= 20 ) + 
  geom_node_point(col=rep( rep("forestgreen",#c("green","orange","red"),
                               each=nrow(createdLayout)),1)  ,alpha=0.5,  size= 20)+#geom_node_point( col= rep( rep(node_col, each= nrow(createdLayout) ) ,  2) , size= 20 ) + 
  geom_node_text(aes(label= name),size=12 ) + #,col= "grey80" ) + #, size= 9) +
  # facet_edges(~DayLab  ) + 
  theme_void(base_size= 1.5* 14 ) +
  theme(aspect.ratio = 1, legend.position = "none" ) + 
  scale_edge_width(range=c(0.1,4))+
  scale_edge_alpha(range=c(0.1,0.97)) +
  scale_edge_color_gradient2()
#scale_edge_color_viridis(option="B")
p0 + scale_edge_color_gradient2(low="white",high="black")
#ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Globe plot ggrraph/Fibroblast Growth Factor Communication globe plot all times.png", height=10, width = 10)


p0 <- ggraph(createdLayout[] , aes(x= x, y= y , col= Day ) ) + 
  geom_edge_loop(  aes(edge_colour=weight,edge_width= 150*(weight),edge_alpha= (weight)) ,  arrow= arrow(length= unit(20, "mm" ) , type= "open"),
                   start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" ) ) +
  geom_edge_arc( aes(edge_colour=weight,edge_width= 150*(weight),edge_alpha= (weight)  ) , arrow= arrow(length= unit(20, "mm" ) , type= "open" ),
                 start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" )  ) +
  #geom_node_point( col= rep( rep(node_col, each= nrow(createdLayout) ) , each= 2) , size= 20 ) + 
  geom_node_point(col=rep( rep("forestgreen",#c("green","orange","red"),
                               each=nrow(createdLayout)),1)  ,alpha=0.5,  size= 20)+#geom_node_point( col= rep( rep(node_col, each= nrow(createdLayout) ) ,  2) , size= 20 ) + 
  # geom_node_text(aes(label= name),size=12 ) + #,col= "grey80" ) + #, size= 9) +
  # facet_edges(~DayLab  ) + 
  theme_void(base_size= 1.5* 14 ) +
  theme(aspect.ratio = 1, legend.position = "none" ) + 
  scale_edge_width(range=c(0.1,4))+
  scale_edge_alpha(range=c(0.1,0.97)) +
  scale_edge_color_gradient2()
#scale_edge_color_viridis(option="B")
p0 + scale_edge_color_gradient2(low="white",high="black")
#ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Globe plot ggrraph/BLANK Fibroblast Growth Factor Communication globe plot all times.png", height=10, width = 10)
#ggsave( file="~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/BLANK Fibroblast Growth Factor Communication globe plot all times.png", height=10, width = 10)


p0 <- ggraph(createdLayout[] , aes(x= x, y= y , col= Day ) ) + 
  geom_edge_loop(  aes(edge_colour=weight,edge_width= 150*exp(weight),edge_alpha= (weight)) ,  arrow= arrow(length= unit(6, "mm" ) , type= "open"),
                   start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" ) ) +
  geom_edge_arc( aes(edge_colour=weight,edge_width= 150*exp(weight),edge_alpha= (weight)  ) , arrow= arrow(length= unit(6, "mm" ) , type= "open" ),
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
#ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Globe plot ggrraph/Letrozole Growth Factor Communication to cancer cells globe plot Day 180 relative to baseline.png", height=10, width = 10)
#ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Globe plot ggrraph/Letrozole Growth Factor Communication globe plot Day 180 relative to baseline.png", height=10, width = 10)



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
#ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Discovery Fibroblast activating stimuli TGFB EGF network with textRED.pdf",width=10, height=10)
#ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Discovery Fibroblast activating stimuli TGFB EGF network with textRED.jpeg",width=10, height=10)


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
#ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/BLANK Discovery Fibroblast activating stimuli TGFB EGF network with textRED.jpeg",width=10, height=10)


discData<-tmp
#write.csv(discData,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/DiscoveryCohortFibroblastactivatingstimuliTGFBEGF.csv")

# Statistical analysis
ppp<-data.table(discovery_LRcomms[LigandPhenoCelltype!="B cells"][ReceptorPhenoCelltype=="Fibroblasts"][Pair.Name%in%l1]%>%
                  group_by(Patient.Study.ID,LigandPhenoCelltype,ReceptorPhenoCelltype,Day,ARM,Treat) %>%
                  dplyr::summarise(muln_scaleTransduction=mean(log(1+scaleTransduction))) )

pppp<- ppp%>%spread(Day,muln_scaleTransduction,fill=NA)
setnames(pppp,old=c("0", "180"),new=c("Day0scores","Day180scores"))
pppp[, supplem:=Day180scores-Day0scores]
pppp$LigandPhenoCelltype <- factor(pppp$LigandPhenoCelltype , levels=c("Fibroblasts","Cancer cells"  ,"Adipocytes","Macrophages", "Normal epithelial cells","Endothelial cells",
                                                                       "CD8+ T cells","Pericytes","CD4+ T cells","B cells" ))
ggplot(pppp[LigandPhenoCelltype!="B cells" ]#[Treat!="CombinationRibo"]
       , aes( y=supplem,x=LigandPhenoCelltype, group=interaction(LigandPhenoCelltype)))+
  geom_boxplot()+geom_point(position=position_dodge(width=1))#+facet_wrap(~Treat)

discStatsData <- ppp[Day==180][LigandPhenoCelltype!="B cells"] 
discStatsData$medianscore <- median(discStatsData$muln_scaleTransduction)

m0 <- lm(muln_scaleTransduction~ -1+ LigandPhenoCelltype ,offset=medianscore,data= discStatsData)
summary(m0)









## Repeat analysis for validation cohort

rm(list=ls())
SourceDataLoc <- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript Feline ERBB facilitation/Molecular Systems Biology/Revision and submission folder/Source Data/"
Intermediateloc <- paste0(SourceDataLoc,"Figure 4/")

validation_LRcomms <- data.table( read.csv(file= paste0(Intermediateloc, "SourceData_Figure_4AB_ValidationTumorComms.csv" ) ) ) 
validation_LRcomms[,Treat:="CombinationRibo"]
validation_LRcomms[ARM=="A",Treat:="LetrozoleAlone"]
validation_LRcomms[,scaleTransduction:=scale(TransductionMu,center=F),by=c("Pair.Name")] #"ReceptorPhenoCelltype",
validation_LRcomms[!is.finite(TransductionMu),scaleTransduction:=0]
validation_LRcomms[,scalelnTransductionMu:=scale(log(1+TransductionMu)),by=Pair.Name]
validation_LRcomms[,scaleTransduction2:=exp(scale(log(TransductionMu),center=F)),by=c("Pair.Name")] 

# Specify relevant receptors: same as in discovery cohort
wchComms<- unique( c( "EGFR","HBEGFR","TGFBR1","TGFBR2","TGFBR3"))#,"IGF1R","AXL" )

# Identify relevant LR pairs
isApath<-sapply( 1:length(wchComms),function(ii){grepl(paste0("_",wchComms[ii]),unique(validation_LRcomms$Pair.Name) ) } ) 
l1 <- unique(validation_LRcomms$Pair.Name)[rowSums(isApath)>0 ]

# Aggregate by sample
average_muln_scaleTransduction <- data.table(validation_LRcomms[Pair.Name%in%l1]%>%
                                               group_by(Patient.Study.ID,LigandPhenoCelltype,ReceptorPhenoCelltype,Day,ARM,Treat) %>%
                                               dplyr::summarise(muln_scaleTransduction=mean(log(1+scaleTransduction))) %>% 
                                               group_by(LigandPhenoCelltype, ReceptorPhenoCelltype, Day, Treat) %>%
                                               dplyr::summarise(average_muln_scaleTransduction=mean(muln_scaleTransduction)))
average_muln_scaleTransduction$LigandPhenoCelltype <- factor(average_muln_scaleTransduction$LigandPhenoCelltype, 
                                                             levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Normal epithelial cells","Cancer cells","B cells","Macrophages","CD8+ T cells","CD4+ T cells"))
average_muln_scaleTransduction$ReceptorPhenoCelltype <- factor(average_muln_scaleTransduction$ReceptorPhenoCelltype, 
                                                               levels= c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Normal epithelial cells","Cancer cells","B cells","Macrophages","CD8+ T cells","CD4+ T cells") )
# Calculate communication relative to patient baseline
initState <- data.table(average_muln_scaleTransduction[Day==0]%>%
                          group_by(LigandPhenoCelltype,ReceptorPhenoCelltype)%>% #Treat, dynamic_class3
                          dplyr::summarise(intiState=mean(average_muln_scaleTransduction)) )

average_muln_scaleTransduction2 <- merge(average_muln_scaleTransduction,initState,by=c("LigandPhenoCelltype", "ReceptorPhenoCelltype")) #,"Treat","dynamic_class3"
average_muln_scaleTransduction2[ , lnfoldchange:= (average_muln_scaleTransduction-intiState) ]      #average_muln_scaleTransduction2[ , lnfoldchange:= (average_muln_scaleTransduction) ] 

# Construct communication network
validation_LRcommscondensed <- average_muln_scaleTransduction2[LigandPhenoCelltype!="B cells"][ReceptorPhenoCelltype!="B cells"]  
validation_LRcommscondensed$LigandPhenoCelltype <-factor(validation_LRcommscondensed$LigandPhenoCelltype, 
                                          levels= sort(c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Normal epithelial cells","Cancer cells","Macrophages","CD8+ T cells","CD4+ T cells"))) #"B cells",
validation_LRcommscondensed$ReceptorPhenoCelltype <-factor(validation_LRcommscondensed$ReceptorPhenoCelltype, 
                                            levels= sort(c("Adipocytes","Pericytes","Fibroblasts","Endothelial cells","Normal epithelial cells","Cancer cells","Macrophages","CD8+ T cells","CD4+ T cells") )) #"B cells",
setnames(validation_LRcommscondensed, old=c("LigandPhenoCelltype", "ReceptorPhenoCelltype", "average_muln_scaleTransduction"), new= c("from", "to", "weight"))

tmp <- validation_LRcommscondensed[][][to=="Fibroblasts"][order( Treat,from,to)]
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
n <- length(unique(validation_LRcommscondensed$to)) -1
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
#ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Validation Fibroblast activating stimuli TGFB EGF network with textRED.pdf",width=10, height=10)
#ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Validation Fibroblast activating stimuli TGFB EGF network with textRED.jpeg",width=10, height=10)


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
#ggsave( file= "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/BLANK Validation Fibroblast activating stimuli TGFB EGF network with textRED.jpeg",width=10, height=10)

validData <- tmp
write.csv(validData,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/ValidationCohortFibroblastactivatingstimuliTGFBEGF.csv")


ppp<-data.table(validation_LRcomms[LigandPhenoCelltype!="B cells"][ReceptorPhenoCelltype=="Fibroblasts"][Pair.Name%in%l1]%>%
                  group_by(Patient.Study.ID,LigandPhenoCelltype,ReceptorPhenoCelltype,Day,ARM,Treat) %>%
                  dplyr::summarise(muln_scaleTransduction=mean(log(1+scaleTransduction))) )

pppp<- ppp%>%spread(Day,muln_scaleTransduction,fill=NA)
setnames(pppp,old=c("0", "180"),new=c("Day0scores","Day180scores"))
pppp[, supplem:=Day180scores-Day0scores]
pppp$LigandPhenoCelltype <- factor(pppp$LigandPhenoCelltype , levels=c("Fibroblasts","Cancer cells"  ,"Adipocytes","Macrophages", "Normal epithelial cells","Endothelial cells",
                                                                       "CD8+ T cells","Pericytes","CD4+ T cells","B cells" ))
ggplot(pppp[LigandPhenoCelltype!="B cells" ]#[Treat!="CombinationRibo"]
       , aes( y=supplem,x=LigandPhenoCelltype, group=interaction(LigandPhenoCelltype)))+
  geom_boxplot()+geom_point(position=position_dodge(width=1))#+facet_wrap(~Treat)

validStatsData <- ppp[Day==180][LigandPhenoCelltype!="B cells"] 
validStatsData$medianscore <- median(validStatsData$muln_scaleTransduction)
m0 <- lm(muln_scaleTransduction~ -1+ LigandPhenoCelltype ,offset=medianscore,data= validStatsData)
summary(m0)
#write.csv(validStatsData,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/ValidationCohort_Stats_FibroblastactivatingstimuliTGFBEGF.csv")

ggplot(validStatsData
       , aes( y=muln_scaleTransduction,x=LigandPhenoCelltype, group=interaction(LigandPhenoCelltype)))+
  geom_boxplot()+geom_point(position=position_dodge(width=1))#+facet_wrap(~Treat)





