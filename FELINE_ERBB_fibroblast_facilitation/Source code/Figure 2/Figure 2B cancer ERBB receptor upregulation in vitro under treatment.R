rm(list=ls())
require(mgcv);require(data.table);require(dplyr);require(ggplot2);require(tidyr);require(lme4);require(lmerTest);require(parallel)
library(effects);require(umap);require(Rfast);require(ider);library("dendextend");library(ggdendro);require(ggsci);require(viridis);require("Rdimtools")
require(merTools)

SourceDataLoc <- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript Feline ERBB facilitation/Molecular Systems Biology/Revision and submission folder/Source Data/"
Intermediateloc <- paste0(SourceDataLoc,"Figure 2/")
list.files(Intermediateloc)

dd <- data.table(read.csv(paste0(Intermediateloc, "SourceData_Figure_2B_in vitro ERBB receptor upreg under treatment data.csv")))

# Summarise mean CT across replicates
dd2<-data.table( dd %>% group_by(cell,resistant,treatment,genes,cDNA.number)%>%dplyr::summarise(CT=mean(CT,na.rm=T)) )

# define controls
controlgene<- "RPLP0_"
controltreatment <- "DMSO"
controldd <- dd2[genes==controlgene]
setnames(controldd, old="CT",new="controlCT")
controldd[,genes:=NULL]

# narmalize and compare expression to controls
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

ggplot( dd3[resistant=="sen"][cell=="CAMA1_"][genes=="TGFb1_"],
        aes(y=mRNAfold,x=treatment))+geom_point()


ggplot( dd3[resistant!="FulR"][genes%in%c("ERBB1_", "ERBB2_", "ERBB3_", "ERBB4_")],
        aes(y=log(mRNAfold),x=genes, fill=treatment,group=interaction(treatment,genes)))+
  geom_boxplot()+
  geom_point(position=position_dodge(width=.75))+
  facet_grid(resistant~cell)

ggplot( dd3[resistant!="FulR"][][genes%in%c("ERBB1_", "ERBB2_", "ERBB3_", "ERBB4_")],
        aes(y=mRNAfold,x=genes, fill=treatment,group=interaction(treatment,genes)))+
  geom_boxplot()+
  geom_point(position=position_dodge(width=.75))+
  facet_grid(~cell)


ggplot( dd3[resistant!="FulR"][][genes%in%c("ERBB1_", "ERBB2_", "ERBB3_", "ERBB4_")],
        aes(y=mRNAfold,x=genes, fill=treatment,group=interaction(treatment,genes)))+
  geom_boxplot()+
  geom_point(position=position_dodge(width=.75))+
  facet_grid(~resistant)

collabs<-(ggsci::pal_aaas()(6)[c(2,3,1)])

plotthis<-dd3[resistant!="FulR"][][genes%in%c("ERBB1_", "ERBB2_", "ERBB3_", "ERBB4_")]
plotthis$cell<- gsub("_","",plotthis$cell)
plotthis[resistant=="sen",cell:=paste0(cell," (ETR)")]
plotthis[resistant=="riboR",cell:=paste0(cell," (CRR)")]
plotthis$cell<- factor(plotthis$cell, levels=unique(plotthis$cell) )
plotthis$genes<- gsub("_","",plotthis$genes)
plotthis[genes=="ERBB1", genes:="EGFR"]

# 
# ggplot( plotthis,aes(y=log(mRNAfold), x=interaction(cell,treatment),fill=treatment, 
#                      group=interaction(cell, resistant,treatment, genes)))+
#   theme_classic(base_size= 28)+
#   geom_hline(yintercept= log(1), linetype= "dashed")+
#   stat_boxplot(geom= "errorbar")+
#   geom_boxplot()+
#   theme(aspect.ratio=4)+
#  # facet_wrap(~treatment,scales="free_x")+#,scale="free_y")+
#   facet_grid(~treatment,scale="free_x")+
#   labs(y="Cancer ERBB receptor expression \n (mRNA fold change relative to DMSO)", x="Cancer cell lineage")+
#   geom_point(aes(shape=genes),position= position_jitterdodge(dodge.width= 0.75, jitter.width= 0.25),size=2)+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   scale_fill_manual(name="Treatment \n (1 week)",
#                     values=c("slategrey", "red", "white","pink"),
#                     labels=c("DMSO", "Fulvestrant (3nM)","Ribociclib (0.5uM)","Combination"))+
#   theme(legend.position="none")+
#   scale_x_discrete(labels=rep(unique(plotthis$cell),4))+
#   scale_y_continuous(labels=c(0.0125,0.25,0.5,1,2,4,8,16,32,64), breaks=log(c(0.0125,0.25,0.5,1,2,4,8,16,32,64)))+
#   theme(axis.text.x=element_text(size=12))+
#   theme(
#     strip.background = element_blank(),
#     strip.text.x = element_blank()
#   )
# ggsave(paste0(fileLoc,"Modeling_JG/","FINALIZE B In vitro ERRB receptor upregulation under treatment.pdf"),width=10,height = 10,dpi=320)


ggplot( plotthis,aes(y=log(mRNAfold), x=interaction(cell,treatment),fill=treatment, 
                     group=interaction(cell, resistant,treatment, genes)))+
  theme_classic(base_size= 28)+
  geom_hline(yintercept= log(1), linetype= "dashed")+
  stat_boxplot(geom= "errorbar")+
  geom_boxplot()+
  theme(aspect.ratio=1)+
  facet_wrap(~genes)+#,scale="free_y")+
  #   facet_grid(genes~.,scale="free_y")+
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

#ggsave(paste0(fileLoc,"Modeling_JG/","FINALIZE In vitro ERRB receptor upregulation under treatment.pdf"),width=10,height = 10,dpi=320)



ggplot( plotthis[treatment%in%c("DMSO","Fulvestrnt 3nM")],aes(y=log(mRNAfold), x=interaction(cell,treatment),fill=treatment, 
                                                              group=interaction(cell, resistant,treatment, genes)))+
  theme_classic(base_size= 28)+
  geom_hline(yintercept= log(1), linetype= "dashed")+
  stat_boxplot(geom= "errorbar")+
  geom_boxplot()+
  theme(aspect.ratio=1)+
  facet_wrap(~genes)+#,scale="free_y")+
  #   facet_grid(genes~.,scale="free_y")+
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

#ggsave(paste0(fileLoc,"Modeling_JG/","FINALIZE MAIN In vitro ERRB receptor upregulation under treatment.pdf"),width=10,height = 10,dpi=320)



