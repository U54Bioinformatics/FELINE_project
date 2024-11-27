rm(list=ls())
require(ggplot2)
require(data.table)
require(dplyr)
require(tidyr)
require(ggsci)
require(lmerTest)
require(mgcv) 

SourceDataLoc <- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript Feline ERBB facilitation/Molecular Systems Biology/Revision and submission folder/Source Data/"
Intermediateloc <- paste0(SourceDataLoc,"Figure 5/")
Cancer_growth_Fulv <- data.table( read.csv( file=paste0(Intermediateloc,
                                              "SourceData_Figure_5D_Cancer facilitation under increasing fulv.csv"  )))

collabs<-(ggsci::pal_aaas()(6)[c(2,3,1)])

ggplot(Cancer_growth_Fulv, aes(y= (RGR-RGRalone), x= sqrt(fulvestrant_nM), 
                     group= interaction(CompositionLab), col= CompositionLab, fill= CompositionLab))+
  stat_boxplot(geom="errorbar",aes(  group= interaction(CellLine, ResistanceLab,fulvestrant_nM,CompositionLab)))+
  #geom_hline(yintercept = 0, linetype="dashed")+
  geom_boxplot(col="black",aes(  group= interaction(CellLine, ResistanceLab,fulvestrant_nM,CompositionLab)))+
  geom_point(aes(  group= interaction(CellLine, ResistanceLab,fulvestrant_nM,CompositionLab)),pch=21, col="black", position=position_dodge(width=1)) +
  theme_classic(base_size= 26)+
  theme(aspect.ratio=1) + 
  geom_smooth(method="gam", formula=y~s(x,k=3,bs="cr"),se=T)+
  #geom_path(aes(group=RepID))+
  #facet_grid(~fulvestrant_nM,scales="free_x")+
  scale_color_manual(name="",values=collabs[-2],labels=c("Cancer \n alone","Cancer and \n fibroblasts"))+
  scale_fill_manual(name="",values=collabs[-2],labels=c("Cancer \n alone","Cancer and \n fibroblasts"))+
  labs(y="Cancer growth rate \n (compared to treatment monoculture)",x="Fulvestrant (nM)")+
  scale_x_continuous(breaks=sqrt(unique(Cancer_growth_Fulv$fulvestrant_nM)) , labels=(unique(Cancer_growth_Fulv$fulvestrant_nM)) )+
  #theme(aspect.ratio=1,legend.position = "bottom")+
  theme(aspect.ratio=1,legend.position = "none")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))#,axis.text.x = element_text(vjust=0.5,angle=90))#+
#ggsave(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Fibroblast_Experiments/coculture/20230305-13_Feng  Fulvestrant dose CAMA1 MCF7 T47D/Modeling_JG/FACILITATION Cancer growth rate in coculture across fulv doses relative to monoculutre across cancer cell lines.pdf", width=8, height=8,dpi=320)
#ggsave( file="~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/FACILITATION Cancer growth rate in coculture across fulv doses relative to monoculutre across cancer cell lines.pdf",height=8.3,width=7,dpi = 320)
#ggsave( file="~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/FACILITATION Cancer growth rate in coculture across fulv doses relative to monoculutre across cancer cell lines.pdf",height=7,width=7,dpi = 320)


##### Stats: Cancer facilitated by Fibro
summary(lmer( I(RGR-RGRalone)~ sqrt(fulvestrant_nM)*CompositionLab+(1+sqrt(fulvestrant_nM)|CancerCellNameResistance),
              Cancer_growth_Fulv))
####
