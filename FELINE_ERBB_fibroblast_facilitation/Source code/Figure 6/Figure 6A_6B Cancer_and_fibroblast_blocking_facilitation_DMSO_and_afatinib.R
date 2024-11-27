rm(list=ls())
require(ggplot2)
require(data.table)
require(dplyr)
require(tidyr)
require(ggsci)
require(lmerTest)
require(mgcv) 

SourceDataLoc <- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript Feline ERBB facilitation/Molecular Systems Biology/Revision and submission folder/Source Data/"
Intermediateloc <- paste0(SourceDataLoc,"Figure 6/")
list.files(Intermediateloc)

FibroCancer_growth_Afat <- data.table( read.csv( file=paste0(Intermediateloc,
                                                        "SourceData_Figure_6A_6B_Cancer_and_fibroblast_blocking_facilitation_DMSO_and_afatinib.csv"  )))

# define color palette
collabs<-(ggsci::pal_aaas()(6)[c(2,3,1)])



fibroblastplot <-ggplot(FibroCancer_growth_Afat[Fulvestrant_nM==0][FibroComLab!="No fibroblast"],aes(y= (RGRfibro-DMSOaloneRGRfibro),#/abs(DMSOaloneRGRfibro),
                                                                                       x= FibroComLab, 
                                                                                       group= interaction(FibroComLab),col=CompositionLab,fill=CompositionLab))+
  geom_hline(yintercept=0, linetype="dashed")+
  stat_boxplot(geom="errorbar")+
  geom_boxplot(col="black")+theme_classic()+
  geom_point(pch=21,col="black")+
  theme_classic(base_size=26)+theme(aspect.ratio = 0.3)+
  #coord_trans(y="log")+
  facet_grid(~labs2 , scales="free_x") +
  scale_color_manual(name="",values=collabs[-1])+
  scale_fill_manual(name="",values=collabs[-1])+
  scale_x_discrete(labels=c("Alone","CAMA1 (CRR)","CAMA1 (ETR)",
                            "MCF7 (CRR)" ,"MCF7 (ETR)",
                            "T47D (CRR)","T47D (ETR)"))+
  labs(x="Fibroblast coculture", y="Fibroblast growth rate \n (compared to DMSO monoculture)")+
  #scale_x_discrete(labels=c("Cancer \n alone", "Cancer + \n unprimed fibroblasts", "Cancer + \n primed fibroblasts"))+
  theme(aspect.ratio=1,legend.position = "none")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1) ,
        axis.text.x = element_text(vjust=0.5,angle=90))#+


cancerplot <- ggplot(FibroCancer_growth_Afat[Fulvestrant_nM==0][FibroComLab!="Alone"],aes(y= (RGR-DMSOaloneRGR),
                                                                            x= interaction(ResistanceLab,CancerCellName,CompositionLab) ,#CompositionLab,#CancerComLab, 
                                                                            group= interaction(ResistanceLab,CancerCellName,CompositionLab),col=CompositionLab,fill=CompositionLab))+
  geom_hline(yintercept=0, linetype="dashed")+
  stat_boxplot(geom="errorbar")+
  geom_boxplot(col="black")+theme_classic()+
  geom_point(position=position_dodge(width=1),pch=21,col="black")+
  theme_classic(base_size=26)+theme(aspect.ratio = 0.3)+
  #coord_trans(y="log")+
  #facet_grid( wrapLab~., scales="free_y") +
  facet_grid(~labs2 , scales="free_x") +
  scale_color_manual(name="",labels=c("Cancer \n alone","Cancer and \n unprimed fibroblasts"),values=collabs[-2])+
  scale_fill_manual(name="",labels=c("Cancer \n alone","Cancer and \n unprimed fibroblasts"),values=collabs[-2])+
  labs(x="Cancer cell line", y="Cancer growth rate \n (compared to monoculture)")+
  scale_x_discrete(labels=rep(c("CAMA1 (CRR)","CAMA1 (ETR)",
                                "MCF7 (CRR)" ,"MCF7 (ETR)",
                                "T47D (CRR)","T47D (ETR)") ,2))+
  theme(aspect.ratio=1,legend.position = "none")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text.x = element_text(vjust=0.5,angle=90))#+


fibroblastplot
#ggsave(fibroblastplot,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Fibroblast_Experiments/coculture/20230223_Feng Afatinib Fulvestrant CAMA1 MCF7 T47D/Modeling_JG/BLOCK MUTUALISM FINALIZED part2 fibroblast growth rate under coculture afatinib treatment.pdf", width=12, height=12, dpi=320)
#ggsave(fibroblastplot,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Fibroblast_Experiments/coculture/20230223_Feng Afatinib Fulvestrant CAMA1 MCF7 T47D/Modeling_JG/BLOCK MUTUALISM FINALIZED part2 fibroblast growth rate under coculture afatinib treatmentB.pdf", width=10, height=10, dpi=320)
#ggsave(fibroblastplot,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Fibroblast_Experiments/coculture/20230223_Feng Afatinib Fulvestrant CAMA1 MCF7 T47D/Modeling_JG/BLOCK MUTUALISM FINALIZED part2 fibroblast growth rate under coculture afatinib treatmentB.pdf", width=10, height=10, dpi=320)
cancerplot 
#ggsave(cancerplot,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Fibroblast_Experiments/coculture/20230223_Feng Afatinib Fulvestrant CAMA1 MCF7 T47D/Modeling_JG/BLOCK MUTUALISM FINALIZED part2 cancer growth rate under coculture afatinib treatment.pdf", width=12, height=12, dpi=320)
#ggsave(cancerplot,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Fibroblast_Experiments/coculture/20230223_Feng Afatinib Fulvestrant CAMA1 MCF7 T47D/Modeling_JG/BLOCK MUTUALISM FINALIZED part2 cancer growth rate under coculture afatinib treatmentB.pdf", width=10, height=10, dpi=320)


#### Stats
##### Fibro facilitated by cancer
summary(lmer( I(RGRfibro-DMSOaloneRGRfibro)~ -1+CompositionLab+(1|CancerCellNameResistance),
              FibroCancer_growth_Afat[Fulvestrant_nM==0][FibroComLab!="No fibroblast"][Afatinib_uMLab=="Control"]))
#Est=5.7e-3,se=6.4e-4,df=25,t=9.0,p=2.5e-9
summary(lmer( I(RGRfibro-DMSOaloneRGRfibro)~ CompositionLab+(1|CancerCellNameResistance),
              FibroCancer_growth_Afat[Fulvestrant_nM==0][FibroComLab!="No fibroblast"][Afatinib_uMLab=="Afatinib: 2.5 uM"]))

##### Fibro growth predicted by afatinib treatment
summary(lmer( I(RGRfibro-DMSOaloneRGRfibro)~ Afatinib_uMLab+(1|CancerCellNameResistance),
              FibroCancer_growth_Afat[Fulvestrant_nM==0][(CompositionLab=="Cancer and unprimed fibroblast") ]))
#Est=-0.0049,se=0.00075,df=29,t=-6.68,p=2.5e-7


summary(lmer( I(RGRfibro-DMSOaloneRGRfibro)~ CompositionLab+(1|CancerCellNameResistance),
              FibroCancer_growth_Afat[Fulvestrant_nM==0][
                (CompositionLab=="Fibroblasts alone"&Afatinib_uMLab=="Control" ) |
                  (CompositionLab=="Cancer and unprimed fibroblast"&Afatinib_uMLab=="Afatinib: 2.5 uM" ) ]))
summary(lmer( I(RGRfibro-DMSOaloneRGRfibro)~ Afatinib_uMLab+(1|CancerCellNameResistance),
              FibroCancer_growth_Afat[Fulvestrant_nM==0][CompositionLab=="Fibroblasts alone"]))
#Est=-5.0e-3,se=1.5e-3,df=14,t=-3.35,p=0.0047 

#summary(lmer( I(RGRfibro-DMSOaloneRGRfibro)~ -1+ CompositionLab*Afatinib_uMLab+(1|CancerCellNameResistance),
#              plotthis1[Fulvestrant_nM==0][FibroComLab!="No fibroblast"][]))

##### Cancer facilitated by Fibro
summary(lmer( I(RGR-DMSOaloneRGR)~ CompositionLab+(1|CancerCellNameResistance),
              FibroCancer_growth_Afat[Fulvestrant_nM==0][FibroComLab!="Alone"][Afatinib_uMLab=="Control"]))
# Est=5.45e-4,se=1.63e-4,df=29,t=3.33,p=0.00235

summary(lmer( I(RGR-DMSOaloneRGR)~ CompositionLab+(1|CancerCellNameResistance),
              FibroCancer_growth_Afat[Fulvestrant_nM==0][FibroComLab!="Alone"][][Afatinib_uMLab=="Afatinib: 2.5 uM"]))
# Est=4.4e-4,se=9.59e-4,df=29,t=-0.46,p=0.65
####

#out<- plotthis1 [Fulvestrant_nM==0]
# save( out, file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/Figure6/Cancer and fibroblast blocking facilitation DMSO and afatinib.RData")
#write.csv( out, file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/Figure6/Cancer and fibroblast blocking facilitation DMSO and afatinib.csv")









fibroblastplot <-ggplot(plotthis1[Fulvestrant_nM==5][FibroComLab!="No fibroblast"],aes(y= (RGRfibro-fulvaloneRGRfibro),#/abs(DMSOaloneRGRfibro),
                                                                                       x= FibroComLab, 
                                                                                       group= interaction(FibroComLab),col=CompositionLab,fill=CompositionLab))+
  geom_hline(yintercept=0, linetype="dashed")+
  stat_boxplot(geom="errorbar")+
  geom_boxplot(col="black")+theme_classic()+
  geom_point(pch=21,col="black")+
  theme_classic(base_size=26)+theme(aspect.ratio = 0.3)+
  #coord_trans(y="log")+
  facet_grid(~labs2 , scales="free_x") +
  scale_color_manual(name="",values=collabs[-1])+
  scale_fill_manual(name="",values=collabs[-1])+
  labs(x="Fibroblast coculture", y="Fibroblast growth rate \n (relative to DMSO monoculutre)")+
  #scale_x_discrete(labels=c("Cancer \n alone", "Cancer + \n unprimed fibroblasts", "Cancer + \n primed fibroblasts"))+
  theme(aspect.ratio=1,legend.position = "bottom")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1) ,
        axis.text.x = element_text(vjust=0.5,angle=90))#+

cancerplot <- ggplot(plotthis1[Fulvestrant_nM==5][FibroComLab!="Alone"],aes(y= (RGR-fulvaloneRGR),
                                                                            x= interaction(ResistanceLab,CancerCellName), 
                                                                            group= interaction(ResistanceLab,CancerCellName,CompositionLab),col=CompositionLab,fill=CompositionLab))+
  geom_hline(yintercept=0, linetype="dashed")+
  stat_boxplot(geom="errorbar")+
  geom_boxplot(col="black")+theme_classic()+
  geom_point(position=position_dodge(width=1),pch=21,col="black")+
  theme_classic(base_size=26)+theme(aspect.ratio = 0.3)+
  #coord_trans(y="log")+
  facet_grid(~labs2 , scales="free_y") +
  scale_color_manual(name="",values=collabs[-2])+
  scale_fill_manual(name="",values=collabs[-2])+
  labs(x="Cancer cell line", y="Cancer growth rate \n (relative to DMSO monoculutre)")+
  scale_x_discrete(labels=c("Resistant \n CAMA1","Sensitive \n CAMA1",
                            "Resistant \n MCF7" ,"Sensitive \n MCF7",
                            "Resistant \n T47D","Sensitive \n T47D"))+
  theme(aspect.ratio=1,legend.position = "bottom")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),axis.text.x = element_text(vjust=0.5,angle=90))#+

fibroblastplot
#ggsave(fibroblastplot,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Fibroblast_Experiments/20230223_Feng Afatinib Fulvestrant CAMA1 MCF7 T47D/Modeling_JG/BLOCK MUTUALISM part2 fibroblast growth rate under coculture afatinib + fulv treatment.pdf", width=12, height=12, dpi=320)

cancerplot 
#ggsave(cancerplot,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Fibroblast_Experiments/20230223_Feng Afatinib Fulvestrant CAMA1 MCF7 T47D/Modeling_JG/BLOCK MUTUALISM part2 cancer growth rate under coculture afatinib + fulv  treatment.pdf", width=12, height=12, dpi=320)




plotthis1[, cancerIDcode:=paste(ResistanceLab,CancerCellName,sep="_")]


# statsFibro<-rbindlist(lapply( as.character(unique(plotthis1[!FibroComLab%in%c("Alone","No fibroblast")]$FibroComLab)) , function(i){
#   data.table(Celltype="Fibroblast",Coculture=i, coef(summary(lm( RGRfibro-fulvaloneRGRfibro~ Afatinib_uMLab*CompositionLab,
#                                                                  data= plotthis1[][Fulvestrant_nM==5][FibroComLab%in%c("Alone",i)]
#   ))) , keep.rownames = T)
# }))
statsFibro<-rbindlist(lapply( as.character(unique(plotthis1[!FibroComLab%in%c("Alone","No fibroblast")]$FibroComLab)) , function(i){
  data.table(Celltype="Fibroblast",Coculture=i, coef(summary(lm( RGRfibro-DMSOaloneRGRfibro~ Afatinib_uMLab*CompositionLab,
                                                                 data= plotthis1[][Fulvestrant_nM==0][FibroComLab%in%c("Alone",i)]
  ))) , keep.rownames = T)
}))
setnames(statsFibro,old=c("Std. Error","t value","Pr(>|t|)"), new=c("Std.Error","tvalue","pvalue"))
statsFibro[pvalue<0.05][rn=="Afatinib_uMLabAfatinib: 2.5 uM:CompositionLabCancer and unprimed fibroblast"]
statsFibro[pvalue<0.05][rn=="CompositionLabCancer and unprimed fibroblast"]
statsFibro[rn=="Afatinib_uMLabAfatinib: 2.5 uM"] [pvalue<0.05]
statsFibro[rn=="CompositionLabCancer and unprimed fibroblast"][pvalue<0.05]

#Resistant CAMA1:est=0.0049,se=0.0012,t=3.90,p=0.001;Sensitive CAMA1:est=0.0052,se=0.0012,t=4.34,p=0.0004;Resistant MCF7:est=0.0029,se=0.0012,t=2.47,p=0.024;Sensitive MCF7:est=0.0034,se=0.0014,t=2.54,p=0.020;Resistant T47D:est=0.0037,se=0.0012,t=3.08,p=0.006;Sensitive T47D:est=0.0037,se=0.0013,t=2.92,p=0.009


plotthis1[,CancerCellNameFibcomp:=cancerIDcode]
plotthis1[CancerPresent==F,CancerCellNameFibcomp:="None"]
summary(lm( RGRfibro-DMSOaloneRGRfibro~ Afatinib_uMLab,
            data= plotthis1[][Fulvestrant_nM==0][CancerCellNameFibcomp=="None"]))
#Afatinib_uMLabAfatinib: 2.5 uM -0.0054069  0.0009858  -5.485 8.04e-05 ***


statsFibro<-rbindlist(lapply( as.character(unique(plotthis1[!FibroComLab%in%c("Alone","No fibroblast")]$CancerCellNameFibcomp)) , function(i){
  data.table(Celltype="Fibroblast",Coculture=i, coef(summary(lm( RGRfibro-DMSOaloneRGRfibro~ Afatinib_uMLab,
                                                                 data= plotthis1[CancerCellNameFibcomp==i][Fulvestrant_nM==0][!FibroComLab%in%c("No fibroblast")]
  ))) , keep.rownames = T)
}))
setnames(statsFibro,old=c("Std. Error","t value","Pr(>|t|)"), new=c("Std.Error","tvalue","pvalue"))
statsFibro[pvalue<0.05][rn=="Afatinib_uMLabAfatinib: 2.5 uM:CompositionLabCancer and unprimed fibroblast"]
statsFibro[pvalue<0.05][rn=="CompositionLabCancer and unprimed fibroblast"]
statsFibro[pvalue<0.05][rn=="Afatinib_uMLabAfatinib: 2.5 uM"]

#Resistant CAMA1:est=-0.0031,se=0.001,t=-3.01,p=0.040;Sensitive CAMA1:est=-0.0039,se=0.001,t=-5.65,p=0.048;Resistant MCF7:est=-0.0037,se=0.0003,t=-11.89,p=0.0003;Sensitive MCF7:est=-0.0084se=0.002,t=-4.94,p=0.0078;Resistant T47D:est=-0.0014,se=0.0004,t=-3.30,p=0.03;Sensitive T47D:est=-0.0045,se=0.001,t=-3.76,p=0.02









statsCancer <- rbindlist(lapply( as.character(unique(plotthis1[!FibroComLab%in%c("Alone","No fibroblast")]$cancerIDcode)) , function(i){
  data.table(Celltype="Cancer",cancerIDcode=i, coef(summary(lm( I(RGR-DMSOaloneRGR) ~ Afatinib_uMLab*CompositionLab,
                                                                data= plotthis1[][Fulvestrant_nM==0][!FibroComLab%in%c("Alone")][cancerIDcode==i]
                                                                
  ))) , keep.rownames = T)
}))
setnames(statsCancer,old=c("Std. Error","t value","Pr(>|t|)"), new=c("Std.Error","tvalue","pvalue"))
statsCancer[pvalue<0.05][rn!="(Intercept)"]
statsCancer[cancerIDcode=="Sensitive_T47D"]
statsCancer[cancerIDcode=="Resistant_MCF7"]
statsCancer[cancerIDcode=="Resistant_CAMA1"]


statsCancerA <- rbindlist(lapply( as.character(unique(plotthis1[!FibroComLab%in%c("Alone","No fibroblast")]$cancerIDcode)) , function(i){
  KW <-     kruskal.test( I(RGR-DMSOaloneRGR)~ Afatinib_uMLab,    data= plotthis1[][Fulvestrant_nM==0][FibroComLab%in%c("No fibroblast")][cancerIDcode==i])
  est=data.table(coef(summary(lm( I(RGR-DMSOaloneRGR)~ Afatinib_uMLab, data= plotthis1[][Fulvestrant_nM==0][FibroComLab%in%c("No fibroblast")][cancerIDcode==i]))))$Estimate[-1]
  data.table(Celltype="Cancer",cancerIDcode=i, method=KW$method, "Kruskal-Wallis chi-squared"=KW$statistic ,est=est,df=KW$parameter,pvalue=KW$p.value)
}))

statsCancerB <- rbindlist(lapply( as.character(unique(plotthis1[!FibroComLab%in%c("Alone","No fibroblast")]$cancerIDcode)) , function(i){
  KW <-     kruskal.test( I(RGR-DMSOaloneRGR)~ Afatinib_uMLab,    data= plotthis1[][Fulvestrant_nM==0][!FibroComLab%in%c("Alone","No fibroblast")][cancerIDcode==i])
  est=data.table(coef(summary(lm( I(RGR-DMSOaloneRGR)~ Afatinib_uMLab, data= plotthis1[][Fulvestrant_nM==0][FibroComLab%in%c("Alone","No fibroblast")][cancerIDcode==i]))))$Estimate[-1]
  data.table(Celltype="Cancer",cancerIDcode=i, method=KW$method, "Kruskal-Wallis chi-squared"=KW$statistic ,est=est,df=KW$parameter,pvalue=KW$p.value)
}))

statsCancer<- rbind(data.table(Coculture="Cancer and Fibroblasts",statsCancerB),data.table(Coculture="Cancer Alone",statsCancerA))
statsCancer[Coculture=="Cancer and Fibroblasts"][pvalue<0.05]
statsCancer[Coculture=="Cancer and Fibroblasts"][pvalue>0.05]
statsCancer[Coculture=="Cancer Alone"][pvalue<0.05]

# coculture p={Resistant CAMA1=0.04630159,Sensitive CAMA1=0.04953461,Resistant MCF7=0.04630159,Sensitive MCF7=0.04953461,Resistant T47D=0.04630159,Sensitive T47D=0.04953461}
#  p={Resistant CAMA1=0.04953461,Sensitive CAMA1=0.04953461,Resistant MCF7=0.82725935,Sensitive MCF7=0.04953461,Resistant T47D=0.51269076,Sensitive T47D=0.04953461}

statsCancerA <- rbindlist(lapply( as.character(unique(plotthis1[!FibroComLab%in%c("Alone","No fibroblast")]$cancerIDcode)) , function(i){
  KW <-     kruskal.test( I(RGR-fulvaloneRGR)~ Afatinib_uMLab,    data= plotthis1[][Fulvestrant_nM==5][FibroComLab%in%c("No fibroblast")][cancerIDcode==i])
  est=data.table(coef(summary(lm( I(RGR-fulvaloneRGR)~ Afatinib_uMLab, data= plotthis1[][Fulvestrant_nM==5][FibroComLab%in%c("No fibroblast")][cancerIDcode==i]))))$Estimate[-1]
  data.table(Celltype="Cancer",cancerIDcode=i, method=KW$method, "Kruskal-Wallis chi-squared"=KW$statistic ,est=est,df=KW$parameter,pvalue=KW$p.value)
}))

statsCancerB <- rbindlist(lapply( as.character(unique(plotthis1[!FibroComLab%in%c("Alone","No fibroblast")]$cancerIDcode)) , function(i){
  KW <-     kruskal.test( I(RGR-fulvaloneRGR)~ Afatinib_uMLab,    data= plotthis1[][Fulvestrant_nM==5][!FibroComLab%in%c("Alone","No fibroblast")][cancerIDcode==i])
  lmmod=data.table(coef(summary(lm( I(RGR-fulvaloneRGR)~ Afatinib_uMLab, data= plotthis1[][Fulvestrant_nM==5][FibroComLab%in%c("Alone","No fibroblast")][cancerIDcode==i]))))
  est<-lmmod$Estimate[-1]
  unlist( lmmod[2,"Pr(>|t|)"] )
  data.table(Celltype="Cancer",cancerIDcode=i, method=KW$method, "Kruskal-Wallis chi-squared"=KW$statistic ,est=est,df=KW$parameter,pvalue=KW$p.value)
}))
statsCancer<- rbind(data.table(Coculture="Cancer and Fibroblasts",statsCancerB),data.table(Coculture="Cancer Alone",statsCancerA))
#setnames(statsCancer,old=c("Std. Error","t value","Pr(>|t|)"), new=c("Std.Error","tvalue","pvalue"))
statsCancer[pvalue<0.05]
statsCancer[!pvalue<0.05]




statsCancer[cancerIDcode=="Sensitive_T47D"]
statsCancer[rn=="Afatinib_uMLabAfatinib: 2.5 uM"]


summary(lmer( I(RGR-DMSOaloneRGR)~ Afatinib_uMLab+(1|cancerIDcode),
              data= plotthis1[][Fulvestrant_nM==0][!FibroComLab%in%c("Alone","No fibroblast")]
))


summary(lmer( I(RGR)~ Afatinib_uMLab*CompositionLab+(1|cancerIDcode),
              data= plotthis1[][Fulvestrant_nM==0][!FibroComLab%in%c("Alone")]
))


data.table( coef(summary(lm( RGRfibro-DMSOaloneRGRfibro~ Afatinib_uMLab*CompositionLab,
                             data= plotthis1[][Fulvestrant_nM==0][FibroComLab%in%c("Alone","Resistant \n CAMA1")]
))) , keep.rownames = T)

summary(lm( RGRfibro-DMSOaloneRGRfibro~ Afatinib_uMLab*CompositionLab,
            data= plotthis1[][Fulvestrant_nM==0][FibroComLab%in%c("Alone","Sensitive \n CAMA1")]
))

summary(lm( RGRfibro-DMSOaloneRGRfibro~ Afatinib_uMLab*CompositionLab,
            data= plotthis1[][Fulvestrant_nM==0][FibroComLab%in%c("Alone","Resistant \n MCF7")]
))

summary(lm( RGRfibro-DMSOaloneRGRfibro~ Afatinib_uMLab*CompositionLab,
            data= plotthis1[][Fulvestrant_nM==0][FibroComLab%in%c("Alone","Sensitive \n MCF7")]
))

summary(lm( RGRfibro-DMSOaloneRGRfibro~ Afatinib_uMLab*CompositionLab,
            data= plotthis1[][Fulvestrant_nM==0][FibroComLab%in%c("Alone","Resistant \n T47D")]
))

summary(lm( RGRfibro-DMSOaloneRGRfibro~ Afatinib_uMLab*CompositionLab,
            data= plotthis1[CancerCellName%in%c("T47D")][Fulvestrant_nM==0][FibroComLab%in%c("Alone","Sensitive \n T47D")]
))

