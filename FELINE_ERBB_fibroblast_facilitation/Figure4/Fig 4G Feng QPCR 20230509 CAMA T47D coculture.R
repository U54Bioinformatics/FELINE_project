rm(list=ls())
require(ggplot2)
require(data.table)
require(dplyr)
require(tidyr)
require(ggsci)
require(lmerTest)
require(lmerTest)

fileLoc<-"/Users/jason/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Fibroblast_Experiments/qPCR/20230509 _QPCR_CAMA1 T47D coculture with MRC5 3D sorting/"
fileNames <- c("2023-05-09 CAMA1 coculture MRC5 all genes.csv" ,"2023-05-09 152147_T47D coculture MRC5 all genes.csv" )  

fileLocMCF7<-"/Users/jason/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Fibroblast_Experiments/qPCR/09-05-2023 164758 MCF7 COCULTURE WITH MRC5 3D 72HRS FACS SORTING/"

fileNamesMCF7 <- c("2023-09-05 164758 MCF7 COCULTURE WITH MRC5 3D 72HRS FACS SORTING all genes.csv")

ddMCF7 <- rbindlist(lapply(1:length(fileNamesMCF7),function(ii){
  return( data.table(  read.csv(paste0(fileLocMCF7,fileNamesMCF7[ii])) ) )
}))
setnames(ddMCF7,old=c("gene"),new=c("genes"))

dd0 <- rbindlist(lapply(1:length(fileNames),function(ii){
  return( data.table(  read.csv(paste0(fileLoc,fileNames[ii])) ) )
}))
setnames(dd0,old=c("culutre.type"),new=c("culture.type"))

dd<- rbind(dd0,ddMCF7)
dd2<-data.table( dd%>%group_by(cell,resistant,coculture,treatment,genes,cDNA.number)%>%dplyr::summarise(CT=mean(CT,na.rm=T)) )
dd2[,coculture:=gsub("_","",coculture)]
dd2[,cell:=gsub("_","",cell)]

controlgene<- "RPLP0_"
controlcoculture <- "No"
controltreatment <- "DMSO"
controldd <- dd2[genes==controlgene ]
setnames(controldd, old="CT",new="controlCT")
controldd[,genes:=NULL]

dd3 <- merge(dd2[genes!=controlgene], controldd, by=names(controldd)[names(controldd)!="controlCT" ])
dd3[, dCT := CT - controlCT ]
dd3[,CT:=NULL]
dd3[,controlCT:=NULL]

dd3[, dCTcontrolaver:=
      sum( (coculture==controlcoculture )*dCT ,na.rm=T)/ sum(coculture==controlcoculture,na.rm=T),
    by= c("genes", "cell", "resistant","treatment")
    ]

dd3[, ddCT:=dCT-dCTcontrolaver]
dd3[, mRNAfold:= 2^-ddCT]
dd3$coculture<- factor(dd3$coculture, levels=c(
  "No", "MRC5", "MRC5-EMT"
))

dd3[,genes:=gsub("_","",genes)]
dd3[,cellresistant:=paste0(cell,resistant)]

collabs<-(ggsci::pal_aaas()(6)[c(2,3,1)])


summary(lmer(log(mRNAfold)~coculture+(1|cellresistant)
                             ,data=dd3[genes%in%c("TGFb1")][treatment%in%c("DMSO")]))
summary(lmer(log(mRNAfold)~coculture+(1|cellresistant)
             ,data=dd3[genes%in%c("TGFb1")][coculture!="No"][treatment%in%c("DMSO")]))
summary(lmer(log(mRNAfold)~coculture+(1|cellresistant)
             ,data=dd3[genes%in%c("TGFb1")][treatment%in%c("Fulvestrnt 5nM")]))
summary(lmer(log(mRNAfold)~coculture+(1|cellresistant)
             ,data=dd3[genes%in%c("TGFb1")][coculture!="No"][treatment%in%c("Fulvestrnt 5nM")]))

summary(lm(log(mRNAfold)~coculture  ,data=dd3[cellresistant=="CAMA1riboR"][coculture!="MRC5"][genes%in%c("HBEGF")][treatment%in%c("Fulvestrnt 5nM")]))
summary(lm(log(mRNAfold)~coculture  ,data=dd3[cellresistant=="CAMA1riboR"][coculture!="MRC5"][genes%in%c("TGFb1")][treatment%in%c("Fulvestrnt 5nM")]))
summary(lm(log(mRNAfold)~coculture  ,data=dd3[cellresistant=="CAMA1riboR"][coculture!="MRC5"][genes%in%c("HBEGF")][treatment%in%c("DMSO")]))
summary(lm(log(mRNAfold)~coculture  ,data=dd3[cellresistant=="CAMA1riboR"][coculture!="MRC5"][genes%in%c("TGFb1")][treatment%in%c("DMSO")]))

summary(lmer(log(mRNAfold)~coculture+(1|cellresistant)
             ,data=dd3[coculture!="MRC5-EMT"][genes%in%c("HBEGF")][treatment%in%c("DMSO")]))

summary(lmer(log(mRNAfold)~coculture+(1|cellresistant)
             ,data=dd3[coculture!="MRC5-EMT"][genes%in%c("TGFb1")][treatment%in%c("DMSO")]))

summary(lmer(log(mRNAfold)~coculture+(1|cellresistant)
             ,data=dd3[coculture!="MRC5-EMT"][genes%in%c("HBEGF")][treatment%in%c("Fulvestrnt 5nM")]))

summary(lmer(log(mRNAfold)~coculture+(1|cellresistant)
             ,data=dd3[coculture!="MRC5-EMT"][genes%in%c("TGFb1")][treatment%in%c("Fulvestrnt 5nM")]))


#######
cancermyCAFsigs<-dd3[genes%in%c("HBEGF","TGFb1")][treatment%in%c("DMSO","Fulvestrnt 5nM")][coculture!="MRC5-EMT"]
cancermyCAFsigs$treatment<- factor(cancermyCAFsigs$treatment, levels=c("DMSO","Fulvestrnt 5nM"))
cancermyCAFsigs$cell2<- cancermyCAFsigs$cell
cancermyCAFsigs[resistant=="riboR",cell2:=paste0(cell2," (CCR)")]
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
#ggsave("~/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/FINALIZE In vitro cancer coculture HBEGF TGFB1 ligand expression under coculture across treatment logscaleB.pdf",height=9,width=11,dpi = 320)

cancermyCAFsigs
#write.csv(cancermyCAFsigs, file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/Figure4/in vitro TGFand HBEGF cocultured cancer ligand upreg under fulv data.csv")

