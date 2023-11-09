rm(list=ls())
require(ggplot2);require(data.table);require(dplyr)
require(tidyr);require(ggsci);require(lmerTest)

fileLoc<-"/Users/jason/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Fibroblast_Experiments/qPCR/20230516 Feng_ERBB_qpcr_CAMA1 MCF7 T47D 2D plus ribo ful combo/"#Feng_ERBB_qpcr_CAMA1 MCF7 T47D/"
fileNames <- c("2023-05-16 140901_CAMA1 2D 7 DAYS 1FBS RIBO PLUS FULV all genes.csv",
               "2023-05-16 163046_MCF7 2D 7 DAYS 1FBS RIBO PLUS FULV RPLP0 TGFB1 HBEGF NRG1 ERBB1-4.csv",
               "2023-05-17 100653_T47D 2D 7 DAYS 1FBS RIBO PLUS FULV RPLP0 TGFB1 HBEGF NRG1 ERBB1-4.csv")  

dd <- rbindlist(lapply(1:length(fileNames),function(ii){
  return( data.table(  read.csv(paste0(fileLoc,fileNames[ii])) ) )
}))

dd2<-data.table( dd%>%group_by(cell,resistant,treatment,genes,cDNA.number)%>%dplyr::summarise(CT=mean(CT,na.rm=T)) )

dd2[genes=="ERBB4_"][cell=="CAMA1_"][resistant=="sen"][treatment%in%c("DMSO","Fulvestrnt 3nM")]
controlgene<- "RPLP0_"
controltreatment <- "DMSO"
controldd <- dd2[genes==controlgene]
setnames(controldd, old="CT",new="controlCT")
controldd[,genes:=NULL]

dd3 <- merge(dd2[genes!=controlgene], controldd, by=names(controldd)[names(controldd)!="controlCT" ])
dd3[, dCT := CT - controlCT ]
dd3[,CT:=NULL]
dd3[,controlCT:=NULL]

dd3[, dCTcontrolaver:=
      sum( (treatment==controltreatment)*dCT )/ sum(treatment==controltreatment ),
    by= c("genes", "cell", "resistant")
    ]
dd3[, ddCT:=dCT-dCTcontrolaver]
dd3[, mRNAfold:= 2^-ddCT]


collabs<-(ggsci::pal_aaas()(6)[c(2,3,1)])

plotthis<-dd3[resistant!="FulR"][][genes%in%c("ERBB1_", "ERBB2_", "ERBB3_", "ERBB4_")]
plotthis$cell<- gsub("_","",plotthis$cell)
plotthis[resistant=="sen",cell:=paste0(cell," (ETR)")]
plotthis[resistant=="riboR",cell:=paste0(cell," (CRR)")]
plotthis$cell<- factor(plotthis$cell, levels=unique(plotthis$cell) )
plotthis$genes<- gsub("_","",plotthis$genes)
plotthis[genes=="ERBB1", genes:="EGFR"]

ggplot( plotthis,aes(y=log(mRNAfold), x=interaction(cell,treatment),fill=treatment, 
                     group=interaction(cell, resistant,treatment, genes)))+
  theme_classic(base_size= 28)+
  geom_hline(yintercept= log(1), linetype= "dashed")+
  stat_boxplot(geom= "errorbar")+
  geom_boxplot()+
  theme(aspect.ratio=1)+
  facet_wrap(~genes)+
  labs(y="Cancer ERBB receptor expression \n (mRNA fold change relative to DMSO)", x="Cancer cell lineage")+
  geom_point(position= position_jitterdodge(dodge.width= 0.75, jitter.width= 0.25),size=2)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_manual(name="Treatment \n (1 week)",
                    values=c("slategrey", "red", "white","pink"),
                    labels=c("DMSO", "Fulvestrant (3nM)","Ribociclib (0.5uM)","Combination"))+
  theme(legend.position="none")+
  scale_x_discrete(labels=rep(unique(plotthis$cell),4))+
  scale_y_continuous(labels=c(0.0125,0.25,0.5,1,2,4,8,16,32,64), breaks=log(c(0.0125,0.25,0.5,1,2,4,8,16,32,64)))+
  theme(axis.text.x=element_text(size=12))

ggsave(paste0(fileLoc,"Modeling_JG/","FINALIZE In vitro ERRB receptor upregulation under treatment.pdf"),width=10,height = 10,dpi=320)


ggplot( plotthis,aes(y=log(mRNAfold), x=interaction(cell,treatment),fill=treatment, 
                     group=interaction(cell, resistant,treatment, genes)))+
  theme_classic(base_size= 18)+
  geom_hline(yintercept= log(1), linetype= "dashed")+
  stat_boxplot(geom= "errorbar")+
  geom_boxplot()+
  theme(aspect.ratio=1)+
  facet_wrap(~genes)+#,scale="free_y")+
  labs(y="Cancer ERBB receptor expression \n (mRNA fold change relative to DMSO)", x="Cancer cell lineage")+
  geom_point(position= position_jitterdodge(dodge.width= 0.75, jitter.width= 0.25),size=2)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_manual(name="Treatment \n (1 week)",
                    values=c("slategrey", "red", "white","pink"),
                    labels=c("DMSO", "Fulvestrant (3nM)","Ribociclib (0.5uM)","Combination"))+
  theme(legend.position="bottom")+
  scale_x_discrete(labels=rep(unique(plotthis$cell),4))+
  scale_y_continuous(labels=c(0.0125,0.25,0.5,1,2,4,8,16,32,64), breaks=log(c(0.0125,0.25,0.5,1,2,4,8,16,32,64)))+
  theme(axis.text.x=element_text(size=12))
ggsave(paste0(fileLoc,"Modeling_JG/","FINALIZE LEGEND In vitro ERRB receptor upregulation under treatment.pdf"),width=10,height = 10,dpi=320)


ggplot( plotthis,aes(y=log(mRNAfold), x=interaction(cell,treatment),fill=treatment, 
                     group=interaction(cell, resistant,treatment, genes)))+
  theme_classic(base_size= 28)+
  geom_hline(yintercept= log(1), linetype= "dashed")+
  stat_boxplot(geom= "errorbar")+
  geom_boxplot()+
  theme(aspect.ratio=1)+
  facet_wrap(~genes)+
  labs(y="Cancer ERBB receptor expression \n (mRNA fold change relative to DMSO)", x="Cancer cell lineage")+
  geom_point(position= position_jitterdodge(dodge.width= 0.75, jitter.width= 0.25),size=2)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_manual(name="Treatment \n (1 week)",
                    values=c("slategrey", "red", "white","pink"),
                    labels=c("DMSO", "Fulvestrant (3nM)","Ribociclib (0.5uM)","Combination"))+
  theme(legend.position="none")+
  scale_x_discrete(labels=rep(rep("",6),4))+
  scale_y_continuous(labels=c(0.0125,0.25,0.5,1,2,4,8,16,32,64), breaks=log(c(0.0125,0.25,0.5,1,2,4,8,16,32,64)))+
  theme(axis.text.x=element_text(size=12))

ggsave(paste0(fileLoc,"Modeling_JG/","FINALIZE In vitro ERRB receptor upregulation under treatment xaxisBLANK.pdf"),width=10,height = 10,dpi=320)





dd4 <- data.table( dd3%>%group_by(cell,resistant,treatment,genes)%>%
  dplyr::summarise(mRNAfold=mean(mRNAfold)) )

ERBBdd<-dd4[resistant!="FulR"][][genes%in%c("ERBB1_", "ERBB2_", "ERBB3_", "ERBB4_")]
ERBBdd[,Treatment:="DMSO"]
ERBBdd[treatment=="Fulvestrnt 3nM",Treatment:="Endocrine"]
ERBBdd[treatment=="Ribociclib 0.5uM",Treatment:="CDK4/6i"]
ERBBdd[treatment=="Ribociclib+Fulv",Treatment:="Combination"]
ERBBdd$Treatment<- factor(ERBBdd$Treatment, levels=c("DMSO","Endocrine","CDK4/6i","Combination"))

ERBBdd[,Gene:=gsub("_","",genes)]
ERBBdd[genes=="ERBB1_", Gene:="EGFR"]
ERBBdd$Gene <- factor(ERBBdd$Gene, levels=c("EGFR","ERBB2","ERBB3","ERBB4"))


ggplot( ERBBdd,
        aes(y=(mRNAfold),x=Gene, fill=Treatment,group=interaction(Treatment,genes)))+
  theme_classic(base_size=26)+
  geom_hline(yintercept=1, linetype="dashed")+
  stat_boxplot(geom="errorbar")+
  geom_boxplot(outlier.colour=NA)+
  geom_point(position=position_jitterdodge(dodge.width=.75, jitter.width=0.15))+
  labs(y="mRNA fold change \n(relative to DMSO mean)",x="Cancer growth factor receptor")+
  scale_fill_manual(name="Treatment \n (1 week)",values=c("slategrey", "red", "white","pink"))+
  theme(aspect.ratio=1)

ggsave("/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/In vitro ERRB receptor upregulation under treatment.pdf",width=9,height = 9,dpi=320)


ggplot( ERBBdd,
        aes(y=log(mRNAfold),x=Gene, fill=Treatment,group=interaction(Treatment,genes)))+
  theme_classic(base_size=26)+
  geom_hline(yintercept=log(1), linetype="dashed")+
  stat_boxplot(geom="errorbar")+
  geom_boxplot(outlier.colour=NA)+
  geom_point(size=2.5,position=position_jitterdodge(dodge.width=.75, jitter.width=0.15))+
  labs(y="mRNA fold change \n(relative to DMSO mean)",x="Cancer growth factor receptor")+
  scale_fill_manual(name="Treatment \n (1 week)",values=c("slategrey", "red", "white","pink"))+
  theme(aspect.ratio=1)+
  scale_y_continuous(labels=c(0.5,1,2,4,8), breaks=log(c(0.5,1,2,4,8)))
ggsave("/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/In vitro ERRB receptor upregulation under treatmentlogscale.pdf",width=9,height = 9,dpi=320)

m0<-lm(mRNAfold~0+Gene+Treatment*Gene ,data= ERBBdd )
m1<-lm(mRNAfold~Gene+Treatment ,data= ERBBdd )
m2<-lm(mRNAfold~0+Treatment ,data= ERBBdd )

anova(m0,m1)
anova(m2,m1)

summary( m1 )

ERBBstatsdd <- dd3[resistant!="FulR"][][genes%in%c("ERBB1_", "ERBB2_", "ERBB3_", "ERBB4_")]
summary(lmer( mRNAfold~treatment+(1+genes|cell),
  data=ERBBstatsdd
))
#write.csv(ERBBstatsdd, file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/Figure2/in vitro ERBB receptor upreg under treatment Stats data.csv")
#write.csv(ERBBdd, file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/Figure2/in vitro ERBB receptor upreg under treatment data.csv")

ERBBdd[,Resistance:="Resistant"]
ERBBdd[resistant=="sen",Resistance:="Sensitive"]

ggplot( ERBBdd,
        aes(y=log(mRNAfold),x=Gene, fill=Treatment,group=interaction(Treatment,genes)))+
  theme_classic(base_size=26)+
  geom_hline(yintercept=log(1), linetype="dashed")+
  stat_boxplot(geom="errorbar")+
  geom_boxplot()+
  geom_point(aes(shape=Resistance),size=2,position=position_jitterdodge(dodge.width=.75, jitter.width=0.25))+
  labs(y="mRNA fold change \n(relative to DMSO mean)",x="Cancer growth factor receptor")+
  scale_fill_manual(name="Treatment \n (1 week)",values=c("slategrey", "red", "white","pink"))+
  theme(aspect.ratio=1)+
  scale_y_continuous(labels=c(0.5,1,2,4,8), breaks=log(c(0.5,1,2,4,8)))

ggsave("/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/In vitro ERRB receptor upregulation under treatmentlogscaleB.pdf",width=9,height = 9,dpi=320)
ggsave(paste0(fileLoc,"Modeling_JG/","In vitro ERRB receptor upregulation under treatmentlogscaleB.pdf"),width=9,height = 9,dpi=320)





ggplot( dd3[resistant!="FulR"][][genes%in%c("HBEGF", "TGFb1_")],
        aes(y=log(mRNAfold),x=gsub("_","",cell), fill=treatment,group=interaction(treatment,cell)))+
  theme_classic(base_size=26)+
  geom_hline(yintercept=log(1), linetype="dashed")+
  stat_boxplot(geom="errorbar")+
  geom_boxplot()+
  geom_point(position=position_jitterdodge(dodge.width=.75, jitter.width=0.25))+
  labs(y="Cancer ligand mRNA fold change \n(relative to DMSO mean)",x="Cancer cell line")+
  scale_fill_manual(name="Treatment \n (1 week)",values=c("slategrey", "red", "white","pink"))+
  theme(aspect.ratio=1)+
  scale_y_continuous(labels=c(0.5,1,2,4,8), breaks=log(c(0.5,1,2,4,8)))+
  facet_grid(gsub("_","",genes)~.)
ggsave(paste0(fileLoc,"Modeling_JG/","In vitro TGFEGF ligand upregulation under treatmentlogscaleB.pdf"),width=9,height = 9,dpi=320)

TFGHBEGFdd<-dd4[resistant!="FulR"][][genes%in%c("HBEGF", "TGFb1_","NRG1_")]
TFGHBEGFdd<-dd4[resistant!="FulR"][][genes%in%c("HBEGF", "TGFb1_")]
TFGHBEGFdd[,Treatment:="DMSO"]
TFGHBEGFdd[treatment=="Fulvestrnt 3nM",Treatment:="Endocrine"]
TFGHBEGFdd[treatment=="Ribociclib 0.5uM",Treatment:="CDK4/6i"]
TFGHBEGFdd[treatment=="Ribociclib+Fulv",Treatment:="Combination"]
TFGHBEGFdd$Treatment<- factor(TFGHBEGFdd$Treatment, levels=c("DMSO","Endocrine","CDK4/6i","Combination"))

TFGHBEGFdd[,Gene:=gsub("_","",genes)]
TFGHBEGFdd[genes=="TGFb1_", Gene:="TGFB1"]
TFGHBEGFdd$Gene <- factor(TFGHBEGFdd$Gene, levels=c("HBEGF","TGFB1"))


ggplot( TFGHBEGFdd[]%>%dplyr::select(-genes)%>%spread(Gene,mRNAfold),
         aes(x=log(HBEGF), col=(cell),y=log(TGFB1)))+
   geom_point()+
   facet_wrap(~Treatment)

ggplot( dd3[resistant!="FulR"][cell=="T47D_"][genes%in%c("HBEGF", "TGFb1_","NRG1_")],
        aes(y=log(mRNAfold),x=genes, fill=treatment,group=interaction(treatment,genes)))+
  theme_classic(base_size=26)+
  geom_hline(yintercept=log(1), linetype="dashed")+
  stat_boxplot(geom="errorbar")+
  geom_boxplot()+
  geom_point(position=position_jitterdodge(dodge.width=.75, jitter.width=0.25))+
  labs(y="mRNA fold change \n(relative to DMSO mean)",x="Cancer growth factor ligands")+
  scale_fill_manual(name="Treatment \n (1 week)",values=c("slategrey", "red", "white","pink"))+
  theme(aspect.ratio=1)+
  scale_y_continuous(labels=c(0.5,1,2,4,8), breaks=log(c(0.5,1,2,4,8)))



ggplot( dd3[resistant!="FulR"][!genes%in%c("ERBB1_","ERBB2_","ERBB3_","ERBB4_")],
        aes(y=log(mRNAfold),x=gsub("_","",cell), fill=treatment,group=interaction(treatment,cell)))+
  theme_classic(base_size=26)+
  geom_hline(yintercept=log(1), linetype="dashed")+
  stat_boxplot(geom="errorbar")+
  geom_boxplot()+
  geom_point(position=position_jitterdodge(dodge.width=.75, jitter.width=0.25))+
  labs(y="mRNA fold change \n(relative to DMSO mean)",x="Cancer growth factor ligands")+
  scale_fill_manual(name="Treatment \n (1 week)",values=c("slategrey", "red", "white","pink"))+
  theme(aspect.ratio=1)+
  scale_y_continuous(labels=c(0.125,0.25,0.5,1,2,4,8), breaks=log(c(0.125,0.25,0.5,1,2,4,8)))+
  facet_wrap(~gsub("_","",genes))
ggsave(paste0(fileLoc,"Modeling_JG/","In vitro TGF HBEGF NRG1 ligand upregulation under treatmentlogscaleB.pdf"),width=16,height = 12,dpi=320)


ggplot( dd3[resistant!="FulR"][!genes%in%c("ERBB1_","ERBB2_","ERBB3_","ERBB4_")],
        aes(y=log(mRNAfold),x=gsub("_","",genes), fill=treatment,group=interaction(treatment,genes)))+
  theme_classic(base_size=26)+
  geom_hline(yintercept=log(1), linetype="dashed")+
  stat_boxplot(geom="errorbar")+
  geom_boxplot()+
  geom_point(position=position_jitterdodge(dodge.width=.75, jitter.width=0.25))+
  labs(y="mRNA fold change \n(relative to DMSO mean)",x="Cancer growth factor ligands")+
  scale_fill_manual(name="Treatment \n (1 week)",values=c("slategrey", "red", "white","pink"))+
  theme(aspect.ratio=1)+
  scale_y_continuous(labels=c(0.125,0.25,0.5,1,2,4,8), breaks=log(c(0.125,0.25,0.5,1,2,4,8)))+
  facet_grid(resistant~gsub("_","",cell))
ggsave(paste0(fileLoc,"Modeling_JG/","In vitro TGF HBEGF NRG1 ligand upregulation under treatmentlogscale.pdf"),width=16,height = 12,dpi=320)


ggplot( TFGHBEGFdd,
        aes(y=log(mRNAfold),x=Gene, fill=Treatment,group=interaction(Treatment,Gene)))+
  theme_classic(base_size=26)+
  geom_hline(yintercept=log(1), linetype="dashed")+
  stat_boxplot(geom="errorbar")+
  geom_boxplot()+
  geom_point(position=position_jitterdodge(dodge.width=.75, jitter.width=0.25))+
  labs(y="mRNA fold change \n(relative to DMSO mean)",x="Cancer growth factor ligands")+
  scale_fill_manual(name="Treatment \n (1 week)",values=c("slategrey", "red", "white","pink"))+
  theme(aspect.ratio=1)+
  scale_y_continuous(labels=c(0.25,0.5,1,2,4,8), breaks=log(c(0.25,0.5,1,2,4,8)))+
  facet_wrap(~cell)



ggplot( TFGHBEGFdd,
        aes(y=log(mRNAfold),x=Gene, fill=Treatment,group=interaction(Treatment,genes)))+
  theme_classic(base_size=26)+
  geom_hline(yintercept=log(1), linetype="dashed")+
  stat_boxplot(geom="errorbar")+
  geom_boxplot()+
  geom_point(position=position_jitterdodge(dodge.width=.75, jitter.width=0.25))+
  labs(y="mRNA fold change \n(relative to DMSO mean)",x="Cancer growth factor ligands")+
  scale_fill_manual(name="Treatment \n (1 week)",values=c("slategrey", "red", "white","pink"))+
  theme(aspect.ratio=1)+
  scale_y_continuous(labels=c(0.25,0.5,1,2,4,8), breaks=log(c(0.25,0.5,1,2,4,8)))

#ggsave("/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/In vitro TGFB HBEGF receptor and NGR1 ligand upregulation under treatmentlogscale.pdf",width=9,height = 9,dpi=320)

ggplot( TFGHBEGFdd,
        aes(y=(mRNAfold),x=Gene, fill=Treatment,group=interaction(Treatment,genes)))+
  theme_classic(base_size=26)+
  geom_hline(yintercept=(1), linetype="dashed")+
  stat_boxplot(geom="errorbar")+
  geom_boxplot()+
  geom_point(position=position_jitterdodge(dodge.width=.75, jitter.width=0.25))+
  labs(y="mRNA fold change \n(relative to DMSO mean)",x="Cancer growth factor ligands")+
  scale_fill_manual(name="Treatment \n (1 week)",values=c("slategrey", "red", "white","pink"))+
  theme(aspect.ratio=1)#+
  #scale_y_continuous(labels=c(0.25,0.5,1,2,4,8), breaks=(c(0.25,0.5,1,2,4,8)))

#ggsave("/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/In vitro TGFB HBEGF receptor and NGR1 ligand upregulation under treatment.pdf",width=9,height = 9,dpi=320)

m0<-lm(mRNAfold~0+Gene+Treatment*Gene ,data= TFGHBEGFdd )
m0b<-lm(mRNAfold~0+Gene+Treatment:Gene ,data= TFGHBEGFdd )
m1<-lm(mRNAfold~Gene+Treatment ,data= TFGHBEGFdd )
m2<-lm(mRNAfold~0+Treatment ,data= TFGHBEGFdd )

summary(lm(mRNAfold~Treatment ,data= TFGHBEGFdd[Gene=="HBEGF"] ))
summary(lm(mRNAfold~Treatment ,data= TFGHBEGFdd[Gene=="TGFB1"] ))
summary(lm(mRNAfold~Treatment ,data= TFGHBEGFdd[Gene=="NRG1"] ))

anova(m0,m1)
anova(m2,m1)

summary( m0b )

require(lmerTest)
ERBBstatsdd <- dd3[resistant!="FulR"][][genes%in%c("ERBB1_", "ERBB2_", "ERBB3_", "ERBB4_")]
summary(lmer( mRNAfold~treatment+(1+genes|cell),
              data=ERBBstatsdd
))
#write.csv(ERBBstatsdd, file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/Figure2/in vitro ERBB receptor upreg under treatment Stats data.csv")
#write.csv(ERBBdd, file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/Figure2/in vitro ERBB receptor upreg under treatment data.csv")

