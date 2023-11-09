rm(list=ls())
require(data.table)
require(ggplot2)
require(dplyr)
require(ggraph)
require(tidygraph)
require(tidyr)
load(file="~/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/temporary/DiscoveryCommunicationNetwork.RData")
DiscoveryCommunicationNetwork <- data.table(Cohort="Discovery",CCIcondensed)


load(file="~/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/temporary/ValidationCommunicationNetwork.RData")
ValidationCommunicationNetwork <- data.table(Cohort="Validation",CCIcondensed)

CCIcondensed <- rbind(DiscoveryCommunicationNetwork,ValidationCommunicationNetwork)
CCIcondensed[ , lnfoldchange:= (average_muln_scaleTransduction-intiState) ] 
# specifcy color of nodes
node_col <- ggsci::pal_npg("nrc")(2)  
# circle layout
n <- length(unique(CCIcondensed$to)) -1
pts.circle <- t(sapply(1:n, function(r)c(cos(2*r*pi/n),sin(2*r*pi/n))))
NodeList <- data.table(c("Cancer cells", "CD4+ T cells" , "Adipocytes", "Fibroblasts", "Normal epithelial cells", "Pericytes"  , "Endothelial cells" ,"Macrophages", "CD8+ T cells" )  , c(0, pts.circle[,1] )  ,  c(0, pts.circle[,2] ) ) #"B cells",
presloc <- data.frame(NodeList)[order(NodeList$V1),]


networks <- data.table(CCIcondensed[Day%in%c(0)][order( Treat,from,to)]%>%group_by(Cohort,from,to,dynamic_class3,Treat)%>%
                    dplyr::summarise(lnfoldchange=mean(lnfoldchange)))
#save(networks, file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/Figure3/Discovery and Validation Communication networks pre treatment.RData")
tmp <- data.table( networks %>% group_by(Treat,
                                    Cohort,dynamic_class3) %>% 
                     mutate(weight = scale(exp(lnfoldchange))^(0.5) ) )

tmp[,CohortLab2 := " "]
tmp[Cohort=="Validation",CohortLab2 := "  "]
tmp[,Treatmentlab:= "Letrozole alone"]
tmp[Treat=="CombinationRibo",Treatmentlab:= "Combination ribociclib"]
tmp[,TumorResponse:="Resistant"]
tmp[dynamic_class3=="Response",TumorResponse:="Sensitive"]
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
createdLayout <- create_layout(graphdd, layout= "star")
createdLayout$x <- presloc$V2
createdLayout$y <- presloc$V3
require(scales)
p0<-ggraph(createdLayout , aes(x= x, y= y , col= dynamic_class3 ) ) + 
  geom_edge_loop( aes(edge_colour=weight,edge_width= 5*(weight),edge_alpha= 5*(weight)) ,  arrow= arrow(length= unit(3.5, "mm" ) , type= "open"),
                   start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" ) ) +
  geom_edge_arc( aes(edge_colour=weight,edge_width= 5*(weight),edge_alpha= 5*(weight)  ) , arrow= arrow(length= unit(3.5, "mm" ) , type= "open" ),
                 start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" )  ) +
  geom_node_point( col= rep( rep(node_col, each= nrow(createdLayout) ) ,  4) , size= 20 ) + 
  geom_node_text(aes(label= name),size=7.5 ) + #,col= "grey80" ) + #, size= 9) +
  facet_edges(paste0(Treatmentlab,CohortLab2)~TumorResponse ,ncol=4 ) + 
  theme_void(base_size= 1.5* 22 ) +
  theme(aspect.ratio = 1, legend.position = "none" ) + 
  scale_edge_width(range=c(0.0,2.5))+
  scale_edge_alpha(range=c(0.0,0.97)) +
  scale_edge_color_gradient2()+
  theme(strip.text=element_text(size=22))

paperfile<- "~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
p0 + scale_edge_color_gradient2(low="white", high="black")
#ggsave(paste0(paperfile,"Discovery and Validation Ribo and Letrozole Globe plot ggraph.png"),height=12,width=20, dpi=320)
p0 + scale_edge_color_gradient2(low="white", high="slategrey")
#ggsave(paste0(paperfile,"Discovery and Validation Ribo and Letrozole Globe plot ggraph grey2.png"),height=12,width=20, dpi=320)







tmp2 <- data.table( networks %>% group_by(Treat,
                                         Cohort,dynamic_class3) %>% 
                     mutate(weight = scale(exp(lnfoldchange))^(0.5) ) )

tmp2 <- data.table(tmp2%>%group_by(Treat,Cohort,dynamic_class3)%>%mutate(
  mu=mean(lnfoldchange,na.rm=T),
  ucl=sd(lnfoldchange,na.rm=T) ))
tmp2[,weightcol:="Low"]
tmp2[lnfoldchange<=(mu+1.2*ucl),weightcol:="Low"]
tmp2[lnfoldchange >(mu+1.2*ucl),weightcol:="High"]

tmp2[,CohortLab2 := " "]
tmp2[Cohort=="Validation",CohortLab2 := "  "]
tmp2[,Treatmentlab:= "Letrozole alone"]
tmp2[Treat=="CombinationRibo",Treatmentlab:= "Combination ribociclib"]
tmp2[,TumorResponse:="Resistant"]
tmp2[dynamic_class3=="Response",TumorResponse:="Sensitive"]
tmp2[from=="Normal epithelial cells",from:="Diploid \n epithelial \n cells"]
tmp2[to=="Normal epithelial cells",to:="Diploid \n epithelial \n cells"]
tmp2[from=="Endothelial cells",from:="Endothelial \n cells"]
tmp2[to=="Endothelial cells",to:="Endothelial \n cells"]
tmp2[from=="CD8+ T cells",from:="CD8+ \n T cells"]
tmp2[to=="CD8+ T cells",to:="CD8+ \n T cells"]
tmp2[from=="CD4+ T cells",from:="CD4+ \n T cells"]
tmp2[to=="CD4+ T cells",to:="CD4+ \n T cells"]
tmp2[from=="Cancer cells",from:="Cancer \n cells"]
tmp2[to=="Cancer cells",to:="Cancer \n cells"]
tmp2[from=="Macrophages",from:="Myeloid \n cells"]
tmp2[to=="Macrophages",to:="Myeloid \n cells"]
graphdd2 <- as_tbl_graph( tmp2, directed= T)
createdLayout2 <- create_layout(graphdd2, layout= "star")
createdLayout2$x <- presloc$V2
createdLayout2$y <- presloc$V3

p1<-ggraph(createdLayout2 , aes(x= x, y= y , col= dynamic_class3 ) ) + 
  geom_edge_loop( aes(edge_colour=weightcol,edge_width= 5*(weight),
                      edge_alpha= 5*(weight)
                      ) ,  
                  arrow= arrow(length= unit(3.5, "mm" ) , type= "open"),
                  start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" ) ) +
  geom_edge_arc( aes(edge_colour=weightcol,edge_width= 5*(weight),
                     edge_alpha= 5*(weight) 
                     ) , 
                 arrow= arrow(length= unit(3.5, "mm" ) , type= "open" ),
                 start_cap= circle(8, "mm" ),  end_cap= circle(8, "mm" )  ) +
  geom_node_point( col= rep( rep(node_col, each= nrow(createdLayout) ) ,  4) , size= 20 ) + 
  geom_node_text(aes(label= name),size=7.5 ) + 
  facet_edges(paste0(Treatmentlab,CohortLab2)~TumorResponse ,ncol=4 ) + 
  theme_void(base_size= 1.5* 22 ) +
  theme(aspect.ratio = 1, legend.position = "none" ) + 
  scale_edge_width(range=c(0.0,2.5))+
  scale_edge_alpha(range=c(0.0,0.97)) +
  scale_edge_color_manual(values=c("black","grey"))+
  theme(
    strip.text=element_text(size=22))

paperfile<- "~/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
p1 
#ggsave(paste0(paperfile,"Discovery and Validation Ribo and Letrozole Globe plot ggraph dichotomous colour.png"),height=12,width=20, dpi=320)


