rm(list=ls())
require(ggplot2)
require(data.table)
require(dplyr)
require(tidyr)
require(ggsci)

# Define sourvce data location and read data 
SourceDataLoc <- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript Feline ERBB facilitation/Molecular Systems Biology/Revision and submission folder/Source Data/"
Intermediateloc <- paste0(SourceDataLoc,"Figure 4/")
dd <- data.table( read.csv( file=paste0(Intermediateloc,"SourceData_Figure_4C_in vitro TGFand HBEGF cancer ligand upreg under fulv data.csv")))

# aggregate data to sample level and adjust annotations for names
dd2<-data.table( dd%>%group_by(cell,resistant,culutre.type,treatment,genes,cDNA.number)%>%dplyr::summarise(CT=mean(CT,na.rm=T)) )
dd2[,coculture:=gsub("_","",culutre.type)]
dd2[,cell:=gsub("_","",cell)]
dd2[,genes:=gsub("_","",genes)]

# define controls
controlgene<- "RPLP0"
#controlcoculture <- "No"
controltreatment <- "DMSO"# "control"
controldd <- dd2[genes==controlgene ]
setnames(controldd, old="CT",new="controlCT")
controldd[,genes:=NULL]

# normalize relative to controls
dd3 <- merge(dd2[genes!=controlgene], controldd, by=names(controldd)[names(controldd)!="controlCT" ])
dd3[, dCT := CT - controlCT ]
dd3[,CT:=NULL]
dd3[,controlCT:=NULL]

dd3[, dCTcontrolaver:=
      sum( ( treatment==controltreatment)*dCT ,na.rm=T)/ sum( treatment==controltreatment ,na.rm=T),
    by= c("genes", "cell", "resistant")
]

dd3[, ddCT:=dCT-dCTcontrolaver]
dd3[, mRNAfold:= 2^-ddCT]

dd3$treatment<-factor(dd3$treatment, levels= c("DMSO","Fulvestrnt 3nM","Fulvestrnt 10nM" ,
                                               "afatinib 0.5uM","afatinib 0.5uM+Fulv 3nM","afatinib 0.5uM+Fulv 10nM" )
)

#extract data for EGFR/TGFBR ligands
plotthis<-dd3[resistant!="fulvR"][genes%in%c("HBEGF","TGFb1")][treatment%in%c("DMSO","Fulvestrnt 10nM"  )]

# create treatment labels
plotthis[,Treatment:="Fulvestrant (3nM)"]
plotthis[treatment=="DMSO",Treatment:="DMSO"]

# order cell line presentation
ordlevs<-c("Resistant \n CAMA1","Sensitive \n CAMA1",
           "Resistant \n MCF7","Sensitive \n MCF7",
           "Resistant \n T47D","Sensitive \n T47D")
plotthis2<-merge( plotthis, data.table(unique(plotthis%>%dplyr::select(cell,resistant )),
                                       cellresistant=ordlevs ), by=c("cell","resistant"))

plotthis2$cellresistant<- factor(plotthis2$cellresistant, levels=ordlevs)


#define color palette
collabs<-(ggsci::pal_aaas()(6)[c(2,3,1)])

# Adjust cell lineage names to match other figures
plotthis2$cell2<- plotthis2$cell
plotthis2[resistant=="riboR",cell2:=paste0(cell2," (CRR)")]
plotthis2[resistant=="sen",cell2:=paste0(cell2," (ETR)")]
#plotthis2$cell2<- factor(plotthis2$cell2, levels=c("CAMA1 (ETR)","CAMA1 (CRR)","MCF7 (ETR)","MCF7 (CRR)","T47D (ETR)","T47D (CRR)" ))

# Visualize data
ggplot( plotthis2,aes(y=log(mRNAfold), x=interaction(resistant,cell,Treatment),fill=Treatment, 
                      group=interaction(cell, Treatment, resistant, coculture, genes)))+
  theme_classic(base_size= 28)+
  geom_hline(yintercept= log(1), linetype= "dashed")+
  stat_boxplot(geom= "errorbar")+
  geom_boxplot()+
  theme(aspect.ratio=1)+
  labs(y="Cancer mRNA fold change \n(relative to DMSO mean)", x="Cancer cell line")+
  scale_y_continuous(labels=c(0.0125,0.25,0.5,1,2,4,8), breaks=log(c(0.0125,0.25,0.5,1,2,4,8)))+
  geom_point(position= position_jitterdodge(dodge.width= 0.75, jitter.width= 0.25))+facet_grid(genes~.,scale="free_y")+
  theme(axis.text.x = element_text(size=rel(0.8),angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_manual(name="Treatment",
                    values=collabs[-1],
                    labels=c("DMSO", "Fulvestrant \n (3nM)"))+
  theme(legend.position="none")+
  scale_x_discrete(labels= rep(unique(plotthis2$cell2),2))

#ggsave(paste0(fileLoc,"Modeling_JG/","FINALIZE In vitro cancer coculture HBEGF TGFB1 ligand expression under fulv vs dmso logscaleb.pdf"),width=7.5,height = 9.75,dpi=320)
#ggsave(paste0("~/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/","FINALIZE In vitro cancer coculture HBEGF TGFB1 ligand expression under fulv vs dmso logscaleb.pdf"),width=7.5,height = 9.75,dpi=320)


#write.csv(plotthis2, file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/Figure4/in vitro TGFand HBEGF cancer ligand upreg under fulv data.csv")
# TGFbeta stats
unique(plotthis2$cellresistant)
tres<-rbindlist(lapply(unique(plotthis2$cellresistant),function(i){
  data.table(genes="TGFb1",cellresistant=i,coef(summary(lm(log(mRNAfold)~-1+Treatment,
                                                           plotthis2[genes=="TGFb1"][cellresistant==i]))), keep.rownames = T)
  
}))
setnames(tres, old=c("Std. Error","t value","Pr(>|t|)"), new=c("Std.Error","tvalue","pvalue"))

# HBEGF stats
hres<-rbindlist(lapply(unique(plotthis2$cellresistant),function(i){
  data.table(genes="HBEGF",cellresistant=i,coef(summary(lm(log(mRNAfold)~-1+Treatment,
                                                           plotthis2[genes=="HBEGF"][cellresistant==i]))), keep.rownames = T)
  
}))
setnames(hres, old=c("Std. Error","t value","Pr(>|t|)"), new=c("Std.Error","tvalue","pvalue"))

tres[rn!="TreatmentDMSO"]
hres[rn!="TreatmentDMSO"]


