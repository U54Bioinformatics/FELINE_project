rm(list=ls())
require(ggplot2)
require(data.table)
require(dplyr)
require(tidyr)
require(ggsci)
require(lmerTest)
require(mgcv) 

# load cocuture growth data for cancer cells
SourceDataLoc <- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript Feline ERBB facilitation/Molecular Systems Biology/Revision and submission folder/Source Data/"
Intermediateloc <- paste0(SourceDataLoc,"Figure 5/")
Cancer_growth <- data.table( read.csv( file=paste0(Intermediateloc,
              "SourceData_Figure_5B_Cancer growth facilitated by fibroblast cells.csv")))

# define color palette
collabs<-(ggsci::pal_aaas()(6)[c(2,3,1)])

#visualize cancer growth data
cancerplot <- ggplot(Cancer_growth,aes(y= (RGR-DMSOaloneRGR),
                                                                         x= interaction(ResistanceLab,CancerCellName,CompositionLab) ,#CompositionLab,#CancerComLab, 
                                                                         group= interaction(ResistanceLab,CancerCellName,CompositionLab),col=CompositionLab,fill=CompositionLab))+
  stat_boxplot(geom="errorbar")+
  geom_boxplot(col="black")+theme_classic()+
  geom_point(position=position_dodge(width=1),pch=21,col="black")+
  theme_classic(base_size=26)+theme(aspect.ratio = 0.3)+
  facet_grid( wrapLab~., scales="free_y") +
  scale_color_manual(name="",labels=c("Cancer \n alone","Cancer and \n fibroblasts"),values=collabs[-2])+
  scale_fill_manual(name="",labels=c("Cancer \n alone","Cancer and \n fibroblasts"),values=collabs[-2])+
  labs(x="Cancer cell line", y="Cancer growth rate \n (compared to monoculture)")+
  scale_x_discrete(labels=rep(c("CAMA1 (CRR)","CAMA1 (ETR)",
                                "MCF7 (CRR)" ,"MCF7 (ETR)",
                                "T47D (CRR)","T47D (ETR)") ,2))+
  theme(aspect.ratio=1,legend.position = "bottom")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text.x = element_text(vjust=0.5,angle=90))+
  theme(legend.position = "none")

cancerplot 
#ggsave(cancerplot,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Fibroblast_Experiments/coculture/20230223_Feng Afatinib Fulvestrant CAMA1 MCF7 T47D/Modeling_JG/MUTUALISM FINALIZED part1 cancer growth rate change under coculture DMSO and fulv treatmentlong.pdf", width=12, height=12, dpi=320)
#ggsave(cancerplot,file="~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/MUTUALISM FINALIZED part1 cancer growth rate change under coculture DMSO and fulv treatmentlong.pdf",height=10,width=12,dpi = 320)


### Cancer facilitation stats
# Linear model prediction of effect of fibroblast coculutre in Cancer growth facilitation: DMSO stats
rbindlist(lapply(unique(Cancer_growth$CancerCellNameResistance)[],function(iii){
  data.table(Fulvestrant_nMLab="Control",
             CancerCellNameResistance=iii,coef(summary(lm( I(RGR-DMSOaloneRGR)~-1+CompositionLab,
                                                           Cancer_growth[CancerCellNameResistance==iii][Fulvestrant_nMLab=="Control"]
              ))),keep.rownames = T)[rn!="CompositionLabCancer alone"]
  
}))[order(CancerCellNameResistance)]

# Linear model prediction of effect of fibroblast coculutre in Cancer growth facilitation: Fulv stats
FresCancer<-rbindlist(lapply(unique(Cancer_growth$CancerCellNameResistance)[],function(iii){
  data.table(Fulvestrant_nMLab="Fulvestrant: 5 nM",
             CancerCellNameResistance=iii,coef(summary(lm( I(RGR-DMSOaloneRGR)~-1+CompositionLab,
                                                           Cancer_growth[CancerCellNameResistance==iii][Fulvestrant_nMLab=="Fulvestrant: 5 nM"]
             ))),keep.rownames = T)[rn!="CompositionLabCancer alone"]
  
}))

##### Liear mixed effects model for each drug treatment : measure cancer facilitated by fibroblasts across cancer cell lines
summary(lmer( I(RGR-DMSOaloneRGR)~ -1+CompositionLab+(1|CancerCellNameResistance),
              plotthis1[Afatinib_uM==0][FibroComLab!="Fibroblasts alone"][][Fulvestrant_nMLab=="Control"]))

summary(lmer( I(RGR-DMSOaloneRGR)~ -1+CompositionLab+(1|CancerCellNameResistance),
              plotthis1[Afatinib_uM==0][FibroComLab!="Fibroblasts alone"][][Fulvestrant_nMLab=="Fulvestrant: 5 nM"]))
####
