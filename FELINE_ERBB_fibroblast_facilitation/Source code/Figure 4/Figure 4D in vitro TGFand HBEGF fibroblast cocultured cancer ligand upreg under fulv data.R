rm(list=ls())
require(ggplot2)
require(data.table)
require(dplyr)
require(tidyr)
require(ggsci)
require(lmerTest)


# Define sourvce data location and read data 
SourceDataLoc <- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript Feline ERBB facilitation/Molecular Systems Biology/Revision and submission folder/Source Data/"
Intermediateloc <- paste0(SourceDataLoc,"Figure 4/")
dd <- data.table( read.csv( file=paste0(Intermediateloc,
                "SourceData_Figure_4D_in vitro TGFand HBEGF fibroblast cocultured cancer ligand upreg under fulv data.csv")))

# aggregate to sample level
dd2<-data.table( dd%>%group_by(cell,resistant,coculture,treatment,genes,cDNA.number)%>%dplyr::summarise(CT=mean(CT,na.rm=T)) )
dd2[,coculture:=gsub("_","",coculture)]
dd2[,cell:=gsub("_","",cell)]

#define controls
controlgene<- "RPLP0_"
controlcoculture <- "No"
controltreatment <- "DMSO"
controldd <- dd2[genes==controlgene ]
setnames(controldd, old="CT",new="controlCT")
controldd[,genes:=NULL]

# merge data with relevant controls
dd3 <- merge(dd2[genes!=controlgene], controldd, by=names(controldd)[names(controldd)!="controlCT" ])
dd3[, dCT := CT - controlCT ]
dd3[,CT:=NULL]
dd3[,controlCT:=NULL]


# Normalize relative to controls
dd3[, dCTcontrolaver:=
      sum( (coculture==controlcoculture & treatment==controltreatment)*dCT ,na.rm=T)/ sum(coculture==controlcoculture& treatment==controltreatment ,na.rm=T),
    by= c("genes", "cell", "resistant")
]

dd3[, dCTcontrolaver:=
      sum( (coculture==controlcoculture )*dCT ,na.rm=T)/ sum(coculture==controlcoculture,na.rm=T),
    by= c("genes", "cell", "resistant","treatment")
]

dd3[, ddCT:=dCT-dCTcontrolaver]
dd3[, mRNAfold:= 2^-ddCT]

# Order cell compositions
dd3$coculture<- factor(dd3$coculture, levels=c(  "No", "MRC5", "MRC5-EMT"))

# Adjust label names
dd3[,genes:=gsub("_", "", genes)]
dd3[,cellresistant:=paste0(cell,resistant)]

# define plot color palette
collabs<-(ggsci::pal_aaas()(6)[c(2,3,1)])

# linear model using coculture status, treatment and cancer resistance state to predict expression change vs controls
summary(lm(log(mRNAfold)~(coculture+
                            treatment+cellresistant)^2 ,data=dd3[genes%in%c("TGFb1")]))

# normalize to monoculture and untreated culture state
dd3[, dCTcontrolaverb:=
      sum( (coculture==controlcoculture&treatment==controltreatment )*dCT ,na.rm=T)/ sum(coculture==controlcoculture&treatment==controltreatment,na.rm=T),
    by= c("genes", "cell", "resistant")
]
dd3[, ddCTb:=dCT-dCTcontrolaverb]
dd3[, mRNAfoldb:= 2^-ddCTb]

# perform Linear mixed effects model to predict log fold changes in mRNA from coculture status (with cancer lineage as random effect)
# repeat for different genes (2) , treatments (2) and coculture cell types (3)
summary(lmer(log(mRNAfold)~coculture+(1|cellresistant)
             ,data=dd3[genes%in%c("TGFb1")][treatment%in%c("DMSO")]))
summary(lmer(log(mRNAfold)~coculture+(1|cellresistant)
             ,data=dd3[genes%in%c("TGFb1")][coculture!="No"][treatment%in%c("DMSO")]))
summary(lmer(log(mRNAfold)~coculture+(1|cellresistant)
             ,data=dd3[genes%in%c("TGFb1")][treatment%in%c("Fulvestrnt 5nM")]))
summary(lmer(log(mRNAfold)~coculture+(1|cellresistant)
             ,data=dd3[genes%in%c("TGFb1")][coculture!="No"][treatment%in%c("Fulvestrnt 5nM")]))

summary(lmer(log(mRNAfold)~coculture+(1|cellresistant)
             ,data=dd3[coculture!="MRC5"][genes%in%c("TGFb1")][treatment%in%c("DMSO")]))

summary(lmer(log(mRNAfold)~coculture+(1|cellresistant)
             ,data=dd3[coculture!="MRC5"][genes%in%c("TGFb1")][treatment%in%c("Fulvestrnt 5nM")]))

summary(lmer(log(mRNAfold)~coculture+(1|cellresistant)
             ,data=dd3[coculture!="MRC5"][genes%in%c("HBEGF")][treatment%in%c("DMSO")]))

summary(lmer(log(mRNAfold)~coculture+(1|cellresistant)
             ,data=dd3[coculture!="MRC5"][genes%in%c("HBEGF")][treatment%in%c("Fulvestrnt 5nM")]))

summary(lmer(log(mRNAfold)~coculture+(1|cellresistant)
             ,data=dd3[coculture!="MRC5-EMT"][genes%in%c("HBEGF")][treatment%in%c("DMSO")]))

summary(lmer(log(mRNAfold)~coculture+(1|cellresistant)
             ,data=dd3[coculture!="MRC5-EMT"][genes%in%c("TGFb1")][treatment%in%c("DMSO")]))

summary(lmer(log(mRNAfold)~coculture+(1|cellresistant)
             ,data=dd3[coculture!="MRC5-EMT"][genes%in%c("HBEGF")][treatment%in%c("Fulvestrnt 5nM")]))

summary(lmer(log(mRNAfold)~coculture+(1|cellresistant)
             ,data=dd3[coculture!="MRC5-EMT"][genes%in%c("TGFb1")][treatment%in%c("Fulvestrnt 5nM")]))

#######
# Visualize data
cancermyCAFsigs <- dd3[genes%in%c("HBEGF","TGFb1")][treatment%in%c("DMSO","Fulvestrnt 5nM")][coculture!="MRC5-EMT"]
cancermyCAFsigs$treatment<- factor(cancermyCAFsigs$treatment, levels=c("DMSO","Fulvestrnt 5nM"))
cancermyCAFsigs$cell2<- cancermyCAFsigs$cell
cancermyCAFsigs[resistant=="riboR",cell2:=paste0(cell2," (CRR)")]
cancermyCAFsigs[resistant=="sen",cell2:=paste0(cell2," (ETR)")]
cancermyCAFsigs[treatment=="Fulvestrnt 5nM",treatment:="Fulvestrant (5nM)"]
ggplot( cancermyCAFsigs,
        aes(y=log(mRNAfold), x=interaction(resistant,cell,coculture),#treatment, 
            fill=coculture, group=interaction(treatment,cell,resistant,coculture, genes)))+
  theme_classic(base_size= 26)+
  facet_grid(genes~treatment,scales="free_x")+
  geom_hline(yintercept= log(1), linetype= "dashed")+
  stat_boxplot(geom= "errorbar")+
  geom_boxplot()+
  geom_point(position= position_jitterdodge(dodge.width= 0.75, jitter.width= 0.25))+
  labs(y="Cancer mRNA fold change \n(relative to treatment monoculture mean)", x="Cancer cell line")+
  scale_fill_manual(name="Coculture \n (1 week)",
                    values=collabs[-2],
                    labels=c("None \n","Unprimed \n fibroblast"))+
  theme(aspect.ratio=1)+
  theme(axis.text.x = element_text(size=rel(0.9),angle = 90, vjust = 0.5, hjust=1))+
  scale_y_continuous(labels=c(0.0125,0.25,0.5,1,2,4,8), breaks=log(c(0.0125,0.25,0.5,1,2,4,8)))+
  scale_x_discrete(labels= rep(unique(cancermyCAFsigs$cell2),2))+
  theme(legend.position="none")

#ggsave(paste0(fileLoc,"Modeling_JG/","FINALIZE In vitro cancer coculture HBEGF TGFB1 ligand expression under coculture across treatment logscaleC.pdf"),width=8,height = 10.5,dpi=320)
#ggsave(paste0(fileLoc,"Modeling_JG/","FINALIZE In vitro cancer coculture HBEGF TGFB1 ligand expression under coculture across treatment logscaleB.pdf"),width=9,height = 11,dpi=320)
#ggsave("~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/FINALIZE In vitro cancer coculture HBEGF TGFB1 ligand expression under coculture across treatment logscaleB.pdf",height=9,width=11,dpi = 320)
#write.csv(cancermyCAFsigs, file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/Figure4/in vitro TGFand HBEGF cocultured cancer ligand upreg under fulv data.csv")