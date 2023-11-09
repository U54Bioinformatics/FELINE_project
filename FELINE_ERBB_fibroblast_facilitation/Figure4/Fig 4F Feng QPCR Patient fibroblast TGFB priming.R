rm(list=ls())
require(ggplot2)
require(data.table)
require(dplyr)
require(tidyr)
require(ggsci)
fileLoc<-"~/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Fibroblast_Experiments/qPCR/2-2-2023 patient CAF TGFb1 activation/"
fileNames <- c("2023-02-02 093714_PATIENT CAF TGFB1 96HRS ALL GENES.csv" )  

cols <- c("Well","Well.Position", "cell","culutre.type" ,
          "treatment","cDNA.number","genes",
          "CT","EQ.Ct.Mean","EQ.Ct.SE" ) 
dd<- data.table(  read.csv(paste0(fileLoc,fileNames[1])))

dd2<-data.table( dd%>%group_by(cell,culutre.type,treatment,genes,cDNA.number)%>%dplyr::summarise(CT=mean(CT,na.rm=T)) )
dd2[,coculture:=gsub("_","",culutre.type)]
dd2[,cell:=gsub("_","",cell)]
dd2[,genes:=gsub("_","",genes)]

controlgene<- "RPLP0"
controltreatment <- "control"
controldd <- dd2[genes==controlgene ]
setnames(controldd, old="CT",new="controlCT")
controldd[,genes:=NULL]


dd3 <- merge(dd2[genes!=controlgene], controldd, by=names(controldd)[names(controldd)!="controlCT" ])
dd3[, dCT := CT - controlCT ]
dd3[,CT:=NULL]
dd3[,controlCT:=NULL]

dd3[, dCTcontrolaver:=
      sum( ( treatment==controltreatment)*dCT ,na.rm=T)/ sum( treatment==controltreatment ,na.rm=T),
    by= c("genes", "cell")
]

dd3[, ddCT:=dCT-dCTcontrolaver]
dd3[, mRNAfold:= 2^-ddCT]

dd3$treatment<-factor(dd3$treatment, levels= c("control","TGFb1 20ng/mL" ))

plotthis<-dd3[][genes%in%c("HBEGF","NRG1")]
plotthis[,Treatment:="TGFb1 (20ng/mL)"]
plotthis[treatment=="control",Treatment:="DMSO"]


plotthis$genes<- factor(plotthis$genes, levels=rev(unique(plotthis$genes)))
collabs<-(ggsci::pal_aaas()(6)[c(2,3,1)])


ggplot( plotthis,aes(y=log(mRNAfold), x=interaction(cell,Treatment),fill=Treatment, 
                     group=interaction(cell, Treatment, coculture, genes)))+
  theme_classic(base_size= 28)+
  geom_hline(yintercept= log(1), linetype= "dashed")+
  stat_boxplot(geom= "errorbar")+
  geom_boxplot()+
  theme(aspect.ratio=1)+
  labs(y="Fibroblast ERBB ligand mRNA \n (fold change under TGFB1 vs control)", x="Primary patient-derived fibroblast")+
  scale_y_continuous(labels=c(0.0125,0.25,0.5,1,2,4,8,16,32,64), breaks=log(c(0.0125,0.25,0.5,1,2,4,8,16,32,64)))+
  geom_point(position= position_jitterdodge(dodge.width= 0.75, jitter.width= 0.25))+facet_grid(genes~.,scale="free_y")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_manual(name="Treatment",
                    values=collabs[-1],
                    labels=c("DMSO", "TGFB1 (20ng/mL)"))+
  theme(legend.position="none")+
  scale_x_discrete(labels=rep(
    paste0("CAF ",1:length(unique(plotthis$cell))) ,2))
ggsave("~/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/FINALIZED In vitro patient fibroblast HBEGF NRG1 ligand expression under tgfb1 vs dmso logscale.pdf",height=9,width=11,dpi = 320)

for(i in 1:length(unique(plotthis$cell))){
  plotthis[cell==unique(plotthis$cell)[i] ,cell:=paste0("CAF ",i)]
}

write.csv(plotthis,file="~/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/Figure4/Fig4F/qPCR TGFBeta fibroblast priming.csv")


hres<-rbindlist(lapply(unique(plotthis$cell),function(i){
  data.table(genes="HBEGF",cellresistant=i,coef(summary(lm(log(mRNAfold)~-1+Treatment,
                                                           plotthis[genes=="HBEGF"][cell==i]))), keep.rownames = T)
  
}))
setnames(hres, old=c("Std. Error","t value","Pr(>|t|)"), new=c("Std.Error","tvalue","pvalue"))

nres<-rbindlist(lapply(unique(plotthis$cell),function(i){
  data.table(genes="NRG1",cellresistant=i,coef(summary(lm(log(mRNAfold)~-1+Treatment,
                                                           plotthis[genes=="NRG1"][cell==i]))), keep.rownames = T)
  
}))
setnames(nres, old=c("Std. Error","t value","Pr(>|t|)"), new=c("Std.Error","tvalue","pvalue"))

nres[rn!="TreatmentDMSO"]
hres[rn!="TreatmentDMSO"]

