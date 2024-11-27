rm(list=ls())
require(ggplot2)
require(data.table)
require(dplyr)
require(tidyr)
require(ggsci)
require(lmerTest)
require(mgcv) 

# load cocuture growth data for fibroblast cells
SourceDataLoc <- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript Feline ERBB facilitation/Molecular Systems Biology/Revision and submission folder/Source Data/"
Intermediateloc <- paste0(SourceDataLoc,"Figure 5/")
Fibroblast_growth <- data.table( read.csv( file=paste0(Intermediateloc,
                                              "SourceData_Figure_5A_Fibroblast growth facilitated by cancer cells.csv")))

# define color palette
collabs<-(ggsci::pal_aaas()(6)[c(2,3,1)])

#visualize fibroblast growth data
fibroblastplot <-ggplot(Fibroblast_growth,aes(y= (RGRfibro-DMSOaloneRGRfibro),#/abs(DMSOaloneRGRfibro),
                                                                                    x= FibroComLab, 
                                                                                    group= interaction(FibroComLab),col=CompositionLab,fill=CompositionLab))+
  stat_boxplot(geom="errorbar")+
  geom_boxplot(col="black")+theme_classic()+
  geom_point(pch=21,col="black")+
  theme_classic(base_size=26)+theme(aspect.ratio = 0.3)+
  facet_grid( wrapLab~., scales="free_x") +
  scale_x_discrete(labels=c("Alone","CAMA1 (CRR)","CAMA1 (ETR)",
                            "MCF7 (CRR)" ,"MCF7 (ETR)",
                            "T47D (CRR)","T47D (ETR)") )+
  scale_color_manual(name="",labels=c("Fibroblast \n alone","Cancer and \n fibroblasts"),values=collabs[-1])+
  scale_fill_manual(name="",labels=c("Fibroblast \n alone","Cancer and \n fibroblasts"),values=collabs[-1])+
  labs(x="Fibroblast coculture", y="Fibroblast growth rate \n (compared to monoculture)")+
  theme(aspect.ratio=1,legend.position = "bottom")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1) ,
        axis.text.x = element_text(vjust=0.5,angle=90))+
  theme(legend.position = "none")
fibroblastplot
#ggsave(fibroblastplot,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Fibroblast_Experiments/coculture/20230223_Feng Afatinib Fulvestrant CAMA1 MCF7 T47D/Modeling_JG/MUTUALISM FINALIZED part1 fibroblast growth rate change under coculture DMSO and fulv treatmentlong.pdf", width=12, height=12, dpi=320)
#ggsave(fibroblastplot,file="~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/MUTUALISM FINALIZED part1 fibroblast growth rate change under coculture DMSO and fulv treatmentlong.pdf",height=10,width=12,dpi = 320)


##### Perform stats : Fibro facilitated by cancer
summary(lmer( I(RGRfibro-DMSOaloneRGRfibro)~ -1+CompositionLab+(1|CancerCellNameResistance),
              Fibroblast_growth[Fulvestrant_nMLab=="Control"]))

summary(lmer( I(RGRfibro-DMSOaloneRGRfibro)~ -1+CompositionLab+(1|CancerCellNameResistance),
              Fibroblast_growth[Fulvestrant_nMLab=="Fulvestrant: 5 nM"]))




