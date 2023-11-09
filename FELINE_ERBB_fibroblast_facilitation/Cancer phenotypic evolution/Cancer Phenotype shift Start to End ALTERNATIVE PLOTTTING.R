rm(list=ls())
require(mgcv);require(data.table);require(dplyr);require(ggplot2);require(tidyr);require(lme4);require(lmerTest);require(parallel)
library(effects);require(umap)
require(Rfast);require(ider)
library("dendextend");library(ggdendro);require(ggsci);require(viridis)
require("Rdimtools")

getstartend<- function(){
  load( file=paste0("~/Dropbox/Cancer_pheno_evo/data/FELINE2/PhenotypesAllArms/newjointNormal_epithelial_cellsCancer cells.RData"))
  #u_dat[TreatLab=="Letrozole"][1:10,1:30]
  ERBB_ssgsea <- read.csv( file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Gene expression landscapes/cancer gene expression landscape/ERBB ssGSEA pathways  shift during Endocrine therapy in cancer names.csv")$x
  
  u_dattrim <- data.table( u_dat[Day %in% c(0, 180)] %>% dplyr::select(c(Cell.ID:V5, one_of(ERBB_ssgsea))) )
  u_dattrim[,hasDay0:=any(Day==0), by=Patient.Study.ID]
  u_dattrim[,hasDay180:=any(Day==180), by=Patient.Study.ID]
  u_dattrim[,hasStartEnd:=hasDay0*hasDay180]
  u_dattrim[,ncellsSampled:=length(unique(Cell.ID)), by=c("Patient.Study.ID","Celltype","Day")]
  u_dattrim[,ncellsSampledDay0:=sum(Day==0), by=c("Patient.Study.ID","Celltype")]
  u_dattrim[,ncellsSampledDay180:=sum(Day==180), by=c("Patient.Study.ID","Celltype")]
  
  start_enddd <-  u_dattrim[hasStartEnd==1]  #start_enddd[1:10,1:20]
  rm(list="u_dattrim")
  return(start_enddd) }
start_enddd <-  getstartend()

lu1 <- data.table( start_enddd%>%group_by(Patient.Study.ID,Timepoint,Celltype)%>%summarise(ncells=n()) %>%spread(Timepoint,ncells, fill=0 ) )
lu1[,Spres:=1]
lu1[S<5,Spres:=0]
lu1[,Epres:=1]
lu1[E<5,Epres:=0]

filtstep<-sort(lu1[Epres&Spres]$Patient.Study.ID%>%table(), decreasing =T) 
filtstep[filtstep>1]
start_enddd[,Timepoint2:="Sensitive"]
start_enddd[Timepoint=="E",Timepoint2:="Resistant"]
start_enddd$Timepoint2 <- factor(start_enddd$Timepoint2 , levels=c("Sensitive","Resistant"))

ERBB_ssgsea <- read.csv( file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Gene expression landscapes/cancer gene expression landscape/ERBB ssGSEA pathways  shift during Endocrine therapy in cancer names.csv")$x

start_endddERBB <- data.table( start_enddd[Patient.Study.ID%in%names(filtstep[filtstep>1]) ]%>%select(one_of(c("Cell.ID","Celltype","Timepoint","TreatLab","Patient.Study.ID",as.character(ERBB_ssgsea) ) )) %>% gather(var,val,all_of(ERBB_ssgsea)) )
start_endddERBB$var <- factor(start_endddERBB$var , levels=ERBB_ssgsea)
start_endddERBB[, valDay0Patient:= sum( val*(Timepoint=="S"))/sum(Timepoint=="S") , by=c( "Celltype", "Patient.Study.ID", "var")]
start_endddERBB[, normscore:= val-valDay0Patient]

start_endddERBB[,Timepoint2:="Sensitive"]
start_endddERBB[Timepoint=="E",Timepoint2:="Resistant"]
start_endddERBB$Timepoint2 <- factor(start_endddERBB$Timepoint2 , levels=c("Sensitive","Resistant"))
start_endddERBB$var <- factor(start_endddERBB$var , levels=rev(ERBB_ssgsea))

# ggplot( start_endddERBB[var=="NAGASHIMA_EGF_SIGNALING_UP"]  , aes(y=normscore,x= Timepoint2, group=Timepoint2) )+ geom_boxplot()+ geom_jitter(width=0.2,size=0.2,alpha=0.15)+
#   #facet_wrap(~Treatment)+
#   scale_x_discrete(name="Cell type", labels= c("Sensitive","Resistant")) +
#   theme_classic(base_size=18)+theme(aspect.ratio=1)+facet_wrap(TreatLab~Celltype)
#

ggplot( start_endddERBB  , aes(y=normscore,x= Timepoint2 , fill=var, group=interaction(var,Timepoint2)) )+ 
  geom_boxplot()+ #geom_jitter(width=0.2,size=0.2,alpha=0.15)+
  #facet_wrap(~Treatment)+
  scale_x_discrete(name="Cell type", labels= c("Sensitive","Resistant")) +
  theme_classic(base_size=18)+theme(aspect.ratio=1)+facet_wrap(TreatLab~Celltype) +
  labs(y="ssGSEA score \n (relative to mean baseline)")
ggsave(file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Cancer vs Normal EGF response of EGFand ERBB during ET ssGSEA single cell.png", height=12, width=24)





ggplot( start_enddd[Patient.Study.ID%in%names(filtstep[filtstep>1]) ]  , aes(y=NAGASHIMA_EGF_SIGNALING_UP,x= Timepoint2, group=Timepoint2) )+ geom_boxplot()+ geom_jitter(width=0.2,size=0.2,alpha=0.15)+
  #facet_wrap(~Treatment)+
  scale_x_discrete(name="Cell type", labels= c("Sensitive","Resistant")) +
  theme_classic(base_size=18)+theme(aspect.ratio=1)+facet_wrap(Treatment~Celltype)


ggplot( start_enddd[Patient.Study.ID== "001-101"]  , aes(y=NAGASHIMA_EGF_SIGNALING_UP,x= Timepoint, group=Timepoint) )+ geom_boxplot()+ geom_jitter(width=0.2,size=0.2,alpha=0.15)+facet_wrap(~Treatment)+
  scale_x_discrete(name="Cell type", labels= c("Sensitive","Resistant")) +
  theme_classic(base_size=18)+theme(aspect.ratio=1)+facet_wrap(~Celltype)



lu1 <- unique( data.table(start_enddd%>% dplyr::select(Patient.Study.ID,Day, Treatment) ) ) #lu1[Day==0][Treatment=="letrozole"]
lu2 <- unique( data.table(start_enddd%>% dplyr::select(Patient.Study.ID, Treatment) ) )

#for each patient select cancer cell data (across timepoints but within one tumor) and perform DE analysis using all genes with >1% coverage
allDE <- rbindlist( lapply(1:nrow(lu2), function(i){
  cat("Tumor    ") ;  cat(i) ; cat("\n")
  dd_i <- start_enddd[ Patient.Study.ID==lu2[i]$Patient.Study.ID]   #dd_i %>% dim()   #ggplot(dd_i, aes(y=sqrt(JUN), x=Day , group=Day)  ) + geom_boxplot() + geom_point() #m1 <- lm(formula= paste0("sqrt(","JUN",") ~ as.factor(Day)") ,    dd_i )
  
  # Perform pathway analysis 
  pathanalysis0 <- rbindlist(lapply( 1:length(common_ssgsea), function(pp){
    cat(pp)
    res <- tryCatch({
      m1 <- lm(formula= paste0(common_ssgsea[pp]," ~ as.factor(Day)") ,    dd_i )
      data.table(gene=common_ssgsea[pp],  coef( summary(m1) ), keep.rownames = T)
    },
    error= function(x){
      return(NULL)
    })
    return(res)  
  } ))
  setnames(pathanalysis0, old= c( "Std. Error", "t value", "Pr(>|t|)" ), new= c( "Std.Error", "tval", "pval" ))
  pathanalysis0$Patient.Study.ID <- dd_i[1]$Patient.Study.ID
  pathanalysis0$ARM <- dd_i[1]$ARM
  pathanalysis0$Treatment <- dd_i[1]$Treatment
  pathanalysis0$dynamic_class <- dd_i[1]$dynamic_class
  pathanalysis0$dynamic_class3 <- dd_i[1]$dynamic_class3
  pathanalysis0$adj.pval <- p.adjust(p= pathanalysis0$pval, method= "fdr" )#,n=length(unique(resultout$Pair.Name)))
  return(pathanalysis0)
}))

allDE[rn!="(Intercept)"][Treatment=="letrozole"][gene=="NAGASHIMA_EGF_SIGNALING_UP"]
# save output
#save(allDE,  file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Gene expression landscapes/cancer gene expression landscape/Tumor specific shift strt vs end cancer SSGSEA expression scrna.RData")
load( file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Gene expression landscapes/cancer gene expression landscape/Tumor specific shift strt vs end cancer SSGSEA expression scrna.RData")

allDE[rn!="(Intercept)"][Treatment=="letrozole"][Estimate!=0][adj.pval<0.05][order(-Estimate)][1:30]
ggplot(allDE[rn!="(Intercept)"][Treatment=="letrozole"][Estimate!=0][adj.pval<0.05] , aes(x=Estimate, y=log(1e-10+adj.pval) , col=abs(Estimate)>log(1.5) ))+geom_point()

summary_allDE <- data.table(allDE[rn!="(Intercept)"][Treatment=="letrozole"][Estimate!=0][adj.pval<0.05]%>%group_by(gene) %>%
                              summarise(nsignif=length(Patient.Study.ID) , 
                                        mean_effect=sum(Estimate)/ (length( unique(allDE[Treatment=="letrozole"]$Patient.Study.ID) ) ) ))
#save(summary_allDE,  file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Gene expression landscapes/cancer gene expression landscape/Tumor specific shift strt vs end cancer SSGSEA expression scrna summary.RData")
load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Gene expression landscapes/cancer gene expression landscape/Tumor specific shift strt vs end cancer SSGSEA expression scrna summary.RData")

ggplot(summary_allDE , aes(x=mean_effect, y=nsignif ) ) + geom_point()
ggplot(summary_allDE[nsignif>=6][order(-mean_effect)] , aes(x=mean_effect, y=nsignif )) + geom_jitter(width=0,height=0.1)
summary_allDE[nsignif>6] [abs(mean_effect)>0.02][order(-abs(mean_effect))]

ggplot( allDE[rn!="(Intercept)"][Treatment=="letrozole"][Estimate!=0][
  gene%in%summary_allDE[nsignif>4]$gene][
    gene%in%summary_allDE[order(-(mean_effect))] [1:20] $gene   ][order(gene)] , 
  aes(Estimate,gene))+geom_jitter(width=0,height=0.1)

#write.csv(summary_allDE,  file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Gene expression landscapes/cancer gene expression landscape/Tumor shift strt vs end Endocrine therapy cancer SSGSEA expression scrna summary.csv")


ConsistentChangePaths <- summary_allDE[!grep("_DN",gene)][nsignif> (length(unique(allDE[Treatment=="letrozole"]$Patient.Study.ID)))*0.5 ][order(-(mean_effect))]
allDE[rn!="(Intercept)"][gene%in%ConsistentChangePaths$gene]
topN <- 5
consistentUpreg <- summary_allDE[!grep("_DN",gene)][nsignif> (length(unique(allDE[Treatment=="letrozole"]$Patient.Study.ID)))*0.5 ][][order(-(mean_effect))] [1:topN] 
consistentUpreg[,dir:="up"]
consistentDownreg <- summary_allDE[!grep("_DN",gene)][nsignif> (length(unique(allDE[Treatment=="letrozole"]$Patient.Study.ID)))*0.5 ][][order((mean_effect))] [topN:1] 
consistentDownreg[,dir:="down"]
consistentgenesUDwn <- na.omit(rbind(consistentUpreg,consistentDownreg))



pltthis1 <- allDE[rn!="(Intercept)"][  gene%in% consistentgenesUDwn $gene] #[Treatment=="letrozole"]
pltthis1$gene <- factor(pltthis1$gene , levels=rev( consistentgenesUDwn$gene) )
pltthis1[,dir:="up"]
pltthis1[gene%in%consistentDownreg$gene,dir:="down"]

pltthis1[gene=="AMIT_DELAYED_EARLY_GENES"]
ggplot( pltthis1, aes(y=Estimate,x=gene) ) + theme_classic(base_size=18) +geom_boxplot(aes(fill=dir,group=gene)) + geom_point()+
  #jitter(width= 0, height= 0.1) + 
  geom_hline(aes(yintercept= 0),linetype=2) +
  labs(x= "Gene set" , y= "Change in cancer cell expression") +facet_wrap(~Treatment)+coord_flip() +
  scale_fill_aaas(name="Change under ET", labels=c("Downregulation","Upregulation")) +theme(aspect.ratio=2)
ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Single cell most significant changes in cancer cell SSGSEA in resistant cells 5.png", height=15, width = 15)


topN <- 50
consistentUpreg <- summary_allDE[!grep("_DN",gene)][nsignif> (length(unique(allDE[Treatment=="letrozole"]$Patient.Study.ID)))*0.5 ][][order(-(mean_effect))] [1:topN] 
consistentUpreg[,dir:="up"]
consistentDownreg <- summary_allDE[!grep("_DN",gene)][nsignif> (length(unique(allDE[Treatment=="letrozole"]$Patient.Study.ID)))*0.5 ][][order((mean_effect))] [topN:1] 
consistentDownreg[,dir:="down"]
consistentgenesUDwn <- na.omit(rbind(consistentUpreg,consistentDownreg))
# write.csv(consistentgenesUDwn , file= "~/Dropbox/Cancer_pheno_evo/data/FELINE2/Gene expression landscapes/cancer gene expression landscape/SignifPathways Tumor specific shift cancer SSGSEA scrna summary.csv")

pltthis1 <- allDE[rn!="(Intercept)"][  gene%in% consistentgenesUDwn $gene] #[Treatment=="letrozole"]
pltthis1$gene <- factor(pltthis1$gene , levels=rev( consistentgenesUDwn$gene) )
pltthis1[,dir:="up"]
pltthis1[gene%in%consistentDownreg$gene,dir:="down"]

ggplot( pltthis1, aes(y=Estimate,x=gene) ) + theme_classic(base_size=18) + geom_point()+
  #jitter(width= 0, height= 0.1) + 
  
  labs(x= "Gene set" , y= "Change in cancer cell expression") +facet_wrap(~Treatment)+geom_boxplot(aes(fill=dir,group=gene)) +coord_flip() +
  geom_hline(aes(yintercept= 0),linetype=2) +
  scale_fill_aaas(name="Change under ET", labels=c("Downregulation","Upregulation")) +theme(aspect.ratio=3, axis.text.y =element_text(size=4) )#element_blank())
ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Single cell most significant changes in cancer cell SSGSEA in resistant cells 100.png", height=12, width = 12)





allDE[rn!="(Intercept)"][Treatment=="letrozole"]
allDE[, Direction:=1]
allDE[Estimate>0 , Direction:=1]
allDE[Estimate<=0 , Direction:=-1]
allDE[,nDir:=length(unique(Patient.Study.ID))  , by=c("Direction","gene","rn","Treatment")]
allDE[,isSignif:=adj.pval<0.05]
allDE[,nSignifDir:=length(isSignif)  , by=c("Direction","gene","rn","Treatment")]
allDE[,mean_effect:=mean(Estimate, na.rm=T)  , by=c("gene","rn","Treatment")]

summary_allDE2 <- data.table( allDE[rn!="(Intercept)"][Treatment=="letrozole"][order(-nSignifDir)]%>%group_by(gene,Direction,rn,Treatment)%>%summarise(
  mean_effect=mean(Estimate, na.rm=T)  , nSignifDir= mean(nSignifDir) ) )

ConsistentChangePaths <- summary_allDE2[!grep("_DN",gene)][nSignifDir> (length(unique(allDE[Treatment=="letrozole"]$Patient.Study.ID)))*0.5 ][mean_effect>0.02][order(-(mean_effect))]

topN <- 5
consistentUpreg0 <- summary_allDE2[!grep("_DN",gene)][Direction==1]
#consistentUpreg <-consistentUpreg0[order(-nSignifDir,-mean_effect)][1:topN]
consistentUpreg <-consistentUpreg0[nSignifDir> (length(unique(allDE[Treatment=="letrozole"]$Patient.Study.ID)))*0.5 ][order(-mean_effect,-nSignifDir)][1:topN]
consistentUpreg[,dir:="up"]
consistentDownreg0 <- summary_allDE2[!grep("_DN",gene)][Direction==-1]
#consistentDownreg <- consistentDownreg0[order(-nSignifDir,mean_effect)][1:topN]
consistentDownreg <- consistentDownreg0[nSignifDir> (length(unique(allDE[Treatment=="letrozole"]$Patient.Study.ID)))*0.5 ][order(mean_effect,-nSignifDir)][1:topN]
consistentDownreg[,dir:="down"]
consistentgenesUDwn <- na.omit(rbind(consistentUpreg,consistentDownreg))[order(-(mean_effect))]
# 
# topN <- 20
# consistentUpreg <- summary_allDE[Direction==1][!grep("_DN",gene)][nSignifDir> (length(unique(allDE[Treatment=="letrozole"]$Patient.Study.ID)))*0.5 ][][order(-(mean_effect))] [1:topN] 
# consistentUpreg[,dir:="up"]
# consistentDownreg <- summary_allDE[Direction==-1][!grep("_DN",gene)][nSignifDir> (length(unique(allDE[Treatment=="letrozole"]$Patient.Study.ID)))*0.75][][order((mean_effect))] [topN:1] 
# consistentDownreg[,dir:="down"]
# consistentgenesUDwn <- na.omit(rbind(consistentUpreg,consistentDownreg))

pltthis1 <- allDE[rn!="(Intercept)"][  gene%in% consistentgenesUDwn[order(-(mean_effect))] $gene] #[Treatment=="letrozole"]
pltthis1$gene <- factor(pltthis1$gene , levels=rev( consistentgenesUDwn[order(-(mean_effect))]$gene) )
pltthis1[,dir:="up"]
pltthis1[gene%in%consistentDownreg$gene,dir:="down"]

ggplot( pltthis1, aes(y=Estimate,x=gene) ) + theme_classic() + geom_point()+geom_boxplot(aes(fill=dir,group=gene)) +
  #jitter(width= 0, height= 0.1) + 
  geom_hline(aes(yintercept= 0),linetype=2) +
  labs(x= "Gene" , y= "Change in cancer cell gene set enrichment") +facet_wrap(~Treatment)+ coord_flip() +
  scale_fill_aaas(name="Change under ET", labels=c("Downregulation","Upregulation")) +theme(aspect.ratio=2)

ggplot( allDE[rn!="(Intercept)"][grep("WNT", gene)][!grep("_DN", gene)], aes(y=Estimate,x=gene) ) + theme_classic() + geom_point()+
  #jitter(width= 0, height= 0.1) + 
  geom_hline(aes(yintercept= 0),linetype=2) +
  labs(x= "Gene" , y= "Fold change in cancer cell expression") +facet_wrap(~Treatment)+geom_boxplot(aes(group=gene)) +coord_flip() +
  scale_fill_aaas(name="Change under ET", labels=c("Downregulation","Upregulation")) +theme(aspect.ratio=2)

ggplot( allDE[rn!="(Intercept)"][grepl("AKT", gene)][!grep("_DN", gene)], aes(y=Estimate,x=gene) ) + theme_classic() + 
  #jitter(width= 0, height= 0.1) + 
  geom_hline(aes(yintercept= 0),linetype=2) +
  labs(x= "Gene" , y= "Fold change in cancer cell expression") +facet_wrap(~Treatment)+geom_boxplot(aes(group=gene)) +geom_point()+ coord_flip() +
  scale_fill_aaas(name="Change under ET", labels=c("Downregulation","Upregulation")) +theme(aspect.ratio=2)


ggplot( allDE[rn!="(Intercept)"][grepl("ESR1", gene)|grepl("ESTROGEN", gene)][!grep("_DN", gene)], aes(y=Estimate,x=gene) ) + theme_classic() + 
  #jitter(width= 0, height= 0.1) + 
  geom_hline(aes(yintercept= 0),linetype=2) +
  labs(x= "Gene" , y= "Change in cancer cell expression") +facet_wrap(~Treatment)+geom_boxplot(aes(group=gene)) +geom_point()+ coord_flip() +
  scale_fill_aaas(name="Change under ET", labels=c("Downregulation","Upregulation")) +theme(aspect.ratio=2)

plt <- allDE[rn!="(Intercept)"][grepl("_EGF", gene)|grepl("ERBB", gene)][
  gene %in% ConsistentChangePaths$gene
  ][!grepl("_DOWNREGULATION", gene)][!grepl("SHC1_EVENT", gene)][!grep("_DN", gene)][order(-mean_effect)]
ordplt <- data.table(plt[Treatment=="letrozole"]%>%group_by(gene)%>%summarise(mean_effect=mean(mean_effect)))[order(mean_effect)]$gene
plt $gene <- factor(plt $gene , levels= ordplt)
ggplot( plt, aes(y=Estimate, x=gene) ) + theme_classic() + 
  #jitter(width= 0, height= 0.1) + 
  geom_hline(aes(yintercept= 0),linetype=2) +
  labs(x= "Gene" , y= "Change in cancer cell expression") +facet_wrap(~Treatment)+geom_boxplot(aes(group=gene)) +geom_point()+ coord_flip() +
  scale_fill_aaas(name="Change under ET", labels=c("Downregulation","Upregulation")) +theme(aspect.ratio=2)





ggplot( plt, aes(y=Estimate,x=gene) ) + theme_classic() + 
  #jitter(width= 0, height= 0.1) + 
  geom_hline(aes(yintercept= 0),linetype=2) +
  labs(x= "Gene" , y= "Change in cancer cell expression") +geom_boxplot(aes(group=gene)) +geom_point()+ coord_flip() +
  scale_fill_aaas(name="Change under ET", labels=c("Downregulation","Upregulation")) +theme(aspect.ratio=2)


# 
start_enddd$Timepoint <- factor(start_enddd$Timepoint  , levels=c("S","E")) 
start_endddERBB <- data.table( start_enddd%>%select(one_of(c("Cell.ID","Timepoint","TreatLab","Patient.Study.ID",ordplt) ))%>%gather(var,val,all_of(ordplt)) )
start_endddERBB$var <- factor(start_endddERBB$var , levels=ordplt)
start_endddERBB[, valDay0Patient:= sum( val*(Timepoint=="S"))/sum(Timepoint=="S") , by=c( "Patient.Study.ID","var")]
start_endddERBB[, normscore:= val-valDay0Patient]

start_endddERBB[, Meannormscore:= mean(normscore), by=c("Cell.ID")]
start_endddERBB <-start_endddERBB[order(-Meannormscore)]
start_endddERBB$Cell.ID <- factor(start_endddERBB$Cell.ID , levels=start_endddERBB[var==start_endddERBB[1]$var][order(-Meannormscore)]$Cell.ID)

start_endddERBB[,Cell_type:="Sensitive"]
start_endddERBB[Timepoint=="E",Cell_type:="Resistant"]
start_endddERBB$Cell_type <- factor(start_endddERBB$Cell_type , levels=c("Sensitive","Resistant"))


#rand <- start_endddERBB[TreatLab=="Letrozole"][sample(1:579414,1000),]$Cell.ID

ggplot( start_endddERBB#[TreatLab=="Letrozole"] # [Cell.ID%in%rand,] 
        ,aes(y= var, x= Cell.ID, fill= normscore ))+geom_tile() + facet_wrap(TreatLab~Cell_type, scales="free_x")+scale_fill_viridis_c(name="ssGSEA score \n (relative to \n mean baseline)") +
  labs(y="ssGSEA pathway", x="Cell")+ theme(aspect.ratio=2, axis.text.x=element_blank())
ggsave(file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Cancer heatmap EGFand ERBB ssGSEA single cell.png", height=6)

ggplot( start_endddERBB#[TreatLab=="Letrozole"]  [Cell.ID%in%rand,] 
        ,aes(y= var, x= Cell.ID, fill= normscore ))+geom_tile() + facet_wrap(TreatLab~Cell_type, scales="free_x")+
  scale_fill_viridis_c(option="A",name="ssGSEA score \n (relative to \n mean baseline)") +
  labs(y="ssGSEA pathway", x="Cell")+ theme(aspect.ratio=2, axis.text.x=element_blank())
ggsave(file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Cancer heatmap EGFand ERBB ssGSEA single cell col2.png")

ggplot( start_endddERBB[TreatLab=="Letrozole"]  [Cell.ID%in%rand,] 
        ,aes(y= var, x= Cell.ID, fill= normscore ))+geom_tile() + facet_wrap(TreatLab~Cell_type, scales="free_x")+
  scale_fill_viridis_c(option="B",name="ssGSEA score \n (relative to \n mean baseline)") +
  labs(y="ssGSEA pathway", x="Cell")+ theme(aspect.ratio=2, axis.text.x=element_blank())
ggsave(file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Cancer heatmap EGFand ERBB ssGSEA single cell col3.png")



ggplot( start_enddd , aes(y=NAGASHIMA_EGF_SIGNALING_UP,x= Timepoint, group=Timepoint) )+ geom_boxplot()+ geom_jitter(width=0.2,size=0.2,alpha=0.15)+facet_wrap(~Treatment)+
  scale_x_discrete(name="Cell type", labels= c("Sensitive","Resistant")) +
  theme_classic(base_size=18)+theme(aspect.ratio=1)
ggsave(file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Cancer NAGASHIMA EGF.png")

start_enddd$Timepoint <- factor(start_enddd$Timepoint  , levels=c("S","E")) 
ggplot( start_enddd , aes(y=NAGASHIMA_EGF_SIGNALING_UP,x= Timepoint, group=Timepoint) )+ geom_boxplot()+ geom_jitter(width=0.2,size=0.2,alpha=0.15)+facet_wrap(~Treatment)+
  scale_x_discrete(name="Cell type", labels= c("Sensitive","Resistant")) +
  theme_classic(base_size=18)+theme(aspect.ratio=1)
ggsave(file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Cancer NAGASHIMA EGF.png")

ggplot( start_enddd , aes(y=WHITEHURST_PACLITAXEL_SENSITIVITY,x= Timepoint, group=Timepoint) )+ geom_boxplot()+ geom_jitter(width=0.2,size=0.2,alpha=0.15)+facet_wrap(~Treatment)+
  scale_x_discrete(name="Cell type", labels= c("Sensitive","Resistant")) +
  theme_classic(base_size=18)+theme(aspect.ratio=1)


ggplot( start_enddd , aes(y=ZHENG_FOXP3_TARGETS_UP,x= Timepoint, group=Timepoint) )+ geom_boxplot()+ geom_jitter(width=0.2,size=0.2,alpha=0.15)+facet_wrap(~Treatment)+
  scale_x_discrete(name="Cell type", labels= c("Sensitive","Resistant")) +
  theme_classic(base_size=18)+theme(aspect.ratio=1)

ggplot( start_enddd , aes(y=CHEN_PDGF_TARGETS,x= Timepoint, group=Timepoint) )+ geom_boxplot()+ geom_jitter(width=0.2,size=0.2,alpha=0.15)+facet_wrap(~Treatment)+
  scale_x_discrete(name="Cell type", labels= c("Sensitive","Resistant")) +
  theme_classic(base_size=18)+theme(aspect.ratio=1)

ggplot( start_enddd , aes(y=RAY_TUMORIGENESIS_BY_ERBB2_CDC25A_UP,x= Timepoint, group=Timepoint) )+ geom_boxplot()+ geom_jitter(width=0.2,size=0.2,alpha=0.15)+facet_wrap(~Treatment)+
  scale_x_discrete(name="Cell type", labels= c("Sensitive","Resistant")) +
  theme_classic(base_size=18)+theme(aspect.ratio=1)

ggplot( start_enddd , aes(y=RAY_TUMORIGENESIS_BY_ERBB2_CDC25A_UP,x= Timepoint, group=Timepoint) )+ geom_boxplot()+ geom_jitter(width=0.2,size=0.2,alpha=0.15)+facet_wrap(~Treatment)+
  scale_x_discrete(name="Cell type", labels= c("Sensitive","Resistant")) +
  theme_classic(base_size=18)+theme(aspect.ratio=1)


ggplot( start_enddd , aes(y=HALLMARK_ESTROGEN_RESPONSE_EARLY,x= Day, group=Day) )+geom_boxplot()+facet_wrap(~Treatment)
ggplot( start_enddd , aes(y=KEGG_ERBB_SIGNALING_PATHWAY,x= Day, group=Day) )+geom_boxplot( )+facet_wrap(~Treatment)


ggplot( start_enddd , aes(y=ZHENG_FOXP3_TARGETS_UP,x= Day, group=Day) )+geom_boxplot()+facet_wrap(~Treatment)
ggplot( start_enddd , aes(y=CHEN_PDGF_TARGETS,x= Day, group=Day) )+geom_boxplot()+facet_wrap(~Treatment)
ggplot( start_enddd , aes(y=CHEN_PDGF_TARGETS,x= Day, group=Day) )+geom_boxplot()+facet_wrap(~Treatment)
ggplot( start_enddd , aes(y=SMID_BREAST_CANCER_LUMINAL_B_UP,x= Day, group=Day) )+geom_boxplot()+facet_wrap(~Treatment)
ggplot( start_enddd , aes(y=PHONG_TNF_TARGETS_UP,x= Day, group=Day) )+geom_boxplot()+facet_wrap(~Treatment)
ggplot( start_enddd , aes(y=AMIT_SERUM_RESPONSE_20_MCF10A,x= Day, group=Day) )+geom_boxplot()+facet_wrap(~Treatment)
ggplot( start_enddd , aes(y=PID_RAS_PATHWAY,x= Day, group=Day) )+geom_boxplot()+facet_wrap(~Treatment)
ggplot( start_enddd , aes(y=KEGG_ERBB_SIGNALING_PATHWAY,x= Day, group=Day) )+geom_boxplot( )+facet_wrap(~Treatment)
HALLMARK_ESTROGEN_RESPONSE_EARLY
ggplot( start_enddd , aes(y=PID_ERBB4_PATHWAY,x= Day, group=Day) )+geom_boxplot()+facet_wrap(~Treatment)
"PID_ERBB1_DOWNSTREAM_PATHWAY"                                   
[397] "PID_ERBB1_INTERNALIZATION_PATHWAY"                               "PID_ERBB1_RECEPTOR_PROXIMAL_PATHWAY"   
PID_ERBB2_ERBB3_PATHWAY
PID_ERBB4_PATHWAY
ggplot( start_enddd , aes(y=NAGASHIMA_EGF_SIGNALING_UP,x= HALLMARK_ESTROGEN_RESPONSE_EARLY, group=Day) )+geom_point(alpha=0.1)+facet_wrap(~Treatment)

#for each patient select cancer cell data (across timepoints but within one tumor) and perform DE analysis using all genes with >1% coverage
subsetcommon_ssgsea <- common_ssgsea[c(grep("HALLMARK",common_ssgsea),grep("BIOCARTA",common_ssgsea),grep("KEGG",common_ssgsea),grep("PID",common_ssgsea)) ]
HALLMARKDE <- rbindlist( lapply(1:nrow(lu2), function(i){
  cat("Tumor    ") ;  cat(i) ; cat("\n")
  dd_i <- start_enddd[ Patient.Study.ID==lu2[i]$Patient.Study.ID]   #dd_i %>% dim()   #ggplot(dd_i, aes(y=sqrt(JUN), x=Day , group=Day)  ) + geom_boxplot() + geom_point() #m1 <- lm(formula= paste0("sqrt(","JUN",") ~ as.factor(Day)") ,    dd_i )
  
  # Perform pathway analysis 
  pathanalysis0 <- rbindlist(lapply( 1:length(subsetcommon_ssgsea), function(pp){
    cat(pp)
    res <- tryCatch({
      m1 <- lm(formula= paste0(subsetcommon_ssgsea[pp]," ~ as.factor(Day)") ,    dd_i )
      data.table(gene=subsetcommon_ssgsea[pp],  coef( summary(m1) ), keep.rownames = T)
    },
    error= function(x){
      return(NULL)
    })
    return(res)  
  } ))
  setnames(pathanalysis0, old= c( "Std. Error", "t value", "Pr(>|t|)" ), new= c( "Std.Error", "tval", "pval" ))
  pathanalysis0$Patient.Study.ID <- dd_i[1]$Patient.Study.ID
  pathanalysis0$ARM <- dd_i[1]$ARM
  pathanalysis0$Treatment <- dd_i[1]$Treatment
  pathanalysis0$dynamic_class <- dd_i[1]$dynamic_class
  pathanalysis0$dynamic_class3 <- dd_i[1]$dynamic_class3
  pathanalysis0$adj.pval <- p.adjust(p= pathanalysis0$pval, method= "fdr" )#,n=length(unique(resultout$Pair.Name)))
  return(pathanalysis0)
}))
HALLMARKDE[rn!="(Intercept)"][Treatment=="letrozole"]
summary_HALLMARKDE <- data.table(HALLMARKDE[rn!="(Intercept)"][Treatment=="letrozole"][Estimate!=0][adj.pval<0.05]%>%group_by(gene) %>%
                                   summarise(nsignif=length(Patient.Study.ID) , 
                                             mean_effect=sum(Estimate)/ (length( unique(allDE[Treatment=="letrozole"]$Patient.Study.ID) ) ) ))

topN <- 20
consistentUpreg <- summary_HALLMARKDE[nsignif> (length(unique(HALLMARKDE[Treatment=="letrozole"]$Patient.Study.ID)))/2 ][][order(-(mean_effect))] [1:topN] 
consistentUpreg[,dir:="up"]
consistentDownreg <- summary_HALLMARKDE[nsignif> (length(unique(HALLMARKDE[Treatment=="letrozole"]$Patient.Study.ID)))/2 ][][order((mean_effect))] [topN:1] 
consistentDownreg[,dir:="down"]
consistentgenesUDwn <- na.omit(rbind(consistentUpreg,consistentDownreg))

pltthis1 <- HALLMARKDE[rn!="(Intercept)"][  gene%in% consistentgenesUDwn $gene] #[Treatment=="letrozole"]
pltthis1$gene <- factor(pltthis1$gene , levels=rev( consistentgenesUDwn$gene) )
pltthis1[,dir:="up"]
pltthis1[gene%in%consistentDownreg$gene,dir:="down"]

ggplot( pltthis1, aes(y=Estimate,x=gene) ) + theme_classic() + geom_point()+
  #jitter(width= 0, height= 0.1) + 
  geom_hline(aes(yintercept= 0),linetype=2) +
  labs(x= "Gene" , y= "Fold change in cancer cell expression") +facet_wrap(~Treatment)+geom_boxplot(aes(fill=dir,group=gene)) +coord_flip() +
  scale_fill_aaas(name="Change under ET", labels=c("Downregulation","Upregulation")) +theme(aspect.ratio=2)





#for each patient select cancer cell data (across timepoints but within one tumor) and perform DE analysis using all genes with >1% coverage
subsetcommon_ssgsea <- unique(common_ssgsea[c(grep("ERBB",common_ssgsea),grep("ESTROGEN",common_ssgsea),grep("EGF",common_ssgsea),grep("ERK",common_ssgsea),grep("AKT",common_ssgsea),grep("MTOR",common_ssgsea),grep("MAPK",common_ssgsea),grep("INSILIN",common_ssgsea)) ])
HALLMARKDE <- rbindlist( lapply(1:nrow(lu2), function(i){
  cat("Tumor    ") ;  cat(i) ; cat("\n")
  dd_i <- start_enddd[ Patient.Study.ID==lu2[i]$Patient.Study.ID]   #dd_i %>% dim()   #ggplot(dd_i, aes(y=sqrt(JUN), x=Day , group=Day)  ) + geom_boxplot() + geom_point() #m1 <- lm(formula= paste0("sqrt(","JUN",") ~ as.factor(Day)") ,    dd_i )
  
  # Perform pathway analysis 
  pathanalysis0 <- rbindlist(lapply( 1:length(subsetcommon_ssgsea), function(pp){
    cat(pp)
    res <- tryCatch({
      m1 <- lm(formula= paste0(subsetcommon_ssgsea[pp]," ~ as.factor(Day)") ,    dd_i )
      data.table(gene=subsetcommon_ssgsea[pp],  coef( summary(m1) ), keep.rownames = T)
    },
    error= function(x){
      return(NULL)
    })
    return(res)  
  } ))
  setnames(pathanalysis0, old= c( "Std. Error", "t value", "Pr(>|t|)" ), new= c( "Std.Error", "tval", "pval" ))
  pathanalysis0$Patient.Study.ID <- dd_i[1]$Patient.Study.ID
  pathanalysis0$ARM <- dd_i[1]$ARM
  pathanalysis0$Treatment <- dd_i[1]$Treatment
  pathanalysis0$dynamic_class <- dd_i[1]$dynamic_class
  pathanalysis0$dynamic_class3 <- dd_i[1]$dynamic_class3
  pathanalysis0$adj.pval <- p.adjust(p= pathanalysis0$pval, method= "fdr" )#,n=length(unique(resultout$Pair.Name)))
  return(pathanalysis0)
}))
HALLMARKDE[rn!="(Intercept)"][Treatment=="letrozole"]
HALLMARKDE[, Direction:=1]
HALLMARKDE[Estimate>0 , Direction:=1]
HALLMARKDE[Estimate<=0 , Direction:=-1]
HALLMARKDE[,nDir:=length(unique(Patient.Study.ID))  , by=c("Direction","gene","rn","Treatment")]
HALLMARKDE[,isSignif:=adj.pval<0.05]
HALLMARKDE[gene=="YOSHIMURA_MAPK8_TARGETS_DN"]
HALLMARKDE[,nSignifDir:=length(isSignif)  , by=c("Direction","gene","rn","Treatment")]

summary_HALLMARKDE <- data.table( HALLMARKDE[rn!="(Intercept)"][Treatment=="letrozole"][order(-nSignifDir)]%>%group_by(gene,Direction,rn,Treatment)%>%summarise(
  mean_effect=mean(Estimate, na.rm=T)  , nSignifDir= mean(nSignifDir) ) )

summary_HALLMARKDE[Direction==1][order(-nSignifDir,-mean_effect)]
summary_HALLMARKDE[Direction==1][order(-mean_effect,-nSignifDir)]

topN <- 20
consistentUpreg <- summary_HALLMARKDE[Direction==1][!grep("_DN",gene)][nSignifDir> (length(unique(HALLMARKDE[Treatment=="letrozole"]$Patient.Study.ID)))/2 ][][order(-(mean_effect))] [1:topN] 
consistentUpreg[,dir:="up"]
consistentDownreg <- summary_HALLMARKDE[Direction==-1][!grep("_DN",gene)][nSignifDir> (length(unique(HALLMARKDE[Treatment=="letrozole"]$Patient.Study.ID)))/2 ][][order((mean_effect))] [topN:1] 
consistentDownreg[,dir:="down"]
consistentgenesUDwn <- na.omit(rbind(consistentUpreg,consistentDownreg))

pltthis1 <- HALLMARKDE[rn!="(Intercept)"][  gene%in% consistentgenesUDwn $gene] #[Treatment=="letrozole"]
pltthis1$gene <- factor(pltthis1$gene , levels=rev( consistentgenesUDwn$gene) )
pltthis1[,dir:="up"]
pltthis1[gene%in%consistentDownreg$gene,dir:="down"]

ggplot( pltthis1, aes(y=Estimate,x=gene) ) + theme_classic() + geom_point()+
  #jitter(width= 0, height= 0.1) + 
  geom_hline(aes(yintercept= 0),linetype=2) +
  labs(x= "Gene" , y= "Fold change in cancer cell expression") +facet_wrap(~Treatment)+geom_boxplot(aes(fill=dir,group=gene)) +coord_flip() +
  scale_fill_aaas(name="Change under ET", labels=c("Downregulation","Upregulation")) +theme(aspect.ratio=2)


ggplot( start_enddd, aes(y=RAY_TUMORIGENESIS_BY_ERBB2_CDC25A_UP, x=Patient.Study.ID,fill=Day, group=interaction(Patient.Study.ID,Day))) +  geom_boxplot()+theme_classic(base_size=20)+ facet_wrap(~Treatment,scales="free_x")
ggplot( start_enddd, aes(y=AMIT_EGF_RESPONSE_40_HELA, x=Patient.Study.ID,fill=Day, group=interaction(Patient.Study.ID,Day))) +  geom_boxplot()+theme_classic(base_size=20)+ facet_wrap(~Treatment,scales="free_x")
ggplot( start_enddd, aes(y=LANDIS_ERBB2_BREAST_TUMORS_324_UP, x=Patient.Study.ID,fill=Day, group=interaction(Patient.Study.ID,Day))) +  geom_boxplot()+theme_classic(base_size=20)+ facet_wrap(~Treatment,scales="free_x")
ggplot( start_enddd, aes(y=NAGASHIMA_EGF_SIGNALING_UP, x=Patient.Study.ID,fill=Day, group=interaction(Patient.Study.ID,Day))) +  geom_boxplot()+theme_classic(base_size=20)+ facet_wrap(~Treatment,scales="free_x")

HALLMARKDE[rn!="(Intercept)"][Treatment=="letrozole"][gene=="ST_ERK1_ERK2_MAPK_PATHWAY"]



summary_HALLMARKDE <- data.table(HALLMARKDE[rn!="(Intercept)"][Treatment=="letrozole"][Estimate!=0][adj.pval<0.05]%>%group_by(gene) %>%
                                   summarise(nsignif=length(Patient.Study.ID) , 
                                             mean_effect=sum(Estimate)/ (length( unique(allDE[Treatment=="letrozole"]$Patient.Study.ID) ) ) ))

topN <- 20
consistentUpreg <- summary_HALLMARKDE[!grep("_DN",gene)][nsignif> (length(unique(HALLMARKDE[Treatment=="letrozole"]$Patient.Study.ID)))/2 ][][order(-(mean_effect))] [1:topN] 
consistentUpreg[,dir:="up"]
consistentDownreg <- summary_HALLMARKDE[!grep("_DN",gene)][nsignif> (length(unique(HALLMARKDE[Treatment=="letrozole"]$Patient.Study.ID)))/2 ][][order((mean_effect))] [topN:1] 
consistentDownreg[,dir:="down"]
consistentgenesUDwn <- na.omit(rbind(consistentUpreg,consistentDownreg))

pltthis1 <- HALLMARKDE[rn!="(Intercept)"][  gene%in% consistentgenesUDwn $gene] #[Treatment=="letrozole"]
pltthis1$gene <- factor(pltthis1$gene , levels=rev( consistentgenesUDwn$gene) )
pltthis1[,dir:="up"]
pltthis1[gene%in%consistentDownreg$gene,dir:="down"]

ggplot( pltthis1, aes(y=Estimate,x=gene) ) + theme_classic() + geom_point()+
  #jitter(width= 0, height= 0.1) + 
  geom_hline(aes(yintercept= 0),linetype=2) +
  labs(x= "Gene" , y= "Fold change in cancer cell expression") +facet_wrap(~Treatment)+geom_boxplot(aes(fill=dir,group=gene)) +coord_flip() +
  scale_fill_aaas(name="Change under ET", labels=c("Downregulation","Upregulation")) +theme(aspect.ratio=2)


ggplot( start_enddd, aes(y=REACTOME_MAPK_TARGETS_NUCLEAR_EVENTS_MEDIATED_BY_MAP_KINASES, x=Patient.Study.ID,fill=Day, group=interaction(Patient.Study.ID,Day))) +  geom_boxplot()+theme_classic(base_size=20)+ facet_wrap(~Treatment,scales="free_x")

ggplot( start_enddd, aes(y=NAGASHIMA_EGF_SIGNALING_UP, x=Patient.Study.ID,fill=Day, group=interaction(Patient.Study.ID,Day))) +  geom_boxplot()+theme_classic(base_size=20)+ facet_wrap(~Treatment,scales="free_x")
ggplot( start_enddd, aes(y=HALLMARK_ESTROGEN_RESPONSE_EARLY, x=Patient.Study.ID,fill=Day, group=interaction(Patient.Study.ID,Day))) +  geom_boxplot()+theme_classic(base_size=20)+ facet_wrap(~Treatment,scales="free_x")
ggplot( start_enddd, aes(y=ST_ERK1_ERK2_MAPK_PATHWAY, x=Patient.Study.ID,fill=Day, group=interaction(Patient.Study.ID,Day))) +  geom_boxplot()+theme_classic(base_size=20)+ facet_wrap(~Treatment,scales="free_x")

ggplot( start_enddd, aes(y=HALLMARK_ESTROGEN_RESPONSE_EARLY, x=NAGASHIMA_EGF_SIGNALING_UP,col=Day, group=interaction(Patient.Study.ID,Day))) + 
  geom_point()+theme_classic(base_size=20)+ facet_wrap(~Patient.Study.ID,scales="free_x")



ggplot( u_dat[Celltype=="Cancer cells"][abs(V1)<5&abs(V2)<5&abs(V3)<5], aes(y=V3, x=V5, col=  Day)) +  geom_point()+theme_classic(base_size=20)

+scale_color_viridis_d(name="Cell type")+scale_fill_aaas(name="Cell type")+theme_classic(base_size=20)



ggplot( u_dat[abs(V2)<5&abs(V3)<5], aes(y=V2, x=V3, col=  Celltype)) +  geom_point()+scale_color_viridis_d(name="Cell type")+scale_fill_aaas(name="Cell type")+theme_classic(base_size=20)+ labs(y="Immune response phenotype (Umap 2)" , x="Proliferative phenotype (Umap 1)")
ggsave(file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/PhenotypeLandscapes/CancerNormal epithelial Celltype annotation.png")
ggplot( u_dat[abs(V2)<5&abs(V3)<5], aes(y=V2, x=V3, col=  Patient.Study.ID,shape=Infercnv_CNA))+ scale_shape_manual(values=c( "diamond", "circle")) +  geom_point(size=1.5)+theme_classic(base_size=20)+theme(aspect.ratio=1) +theme(legend.position = "none")+ labs(y="Immune response phenotype (Umap 2)" , x="Proliferative phenotype (Umap 1)")
ggsave(file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/PhenotypeLandscapes/CancerNormal epithelial Patient.Study.ID annotation.png")
#ggplot( u_dat[abs(V2)<5&abs(V3)<5],aes(y=V2, x=V3, col=  Patient.Study.ID ))+   geom_point(size=1.5)+theme(aspect.ratio=1)+theme_classic(base_size=20) +theme(legend.position = "none") + labs(y="Immune response phenotype (Umap 2)" , x="Proliferative phenotype (Umap 1)")

ggplot( u_dat[abs(V2)<5&abs(V3)<5], aes(y=V2, x=V3, col=  dynamic_class3)) +  geom_point(size=1.5)+scale_color_npg(name="Tumor response", labels=c("Resistant","Sensitive"))+theme_classic(base_size=20)+theme(aspect.ratio=1)+ labs(y="Immune response phenotype (Umap 2)" , x="Proliferative phenotype (Umap 1)") #+facet_wrap(~Day)
ggsave(file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/PhenotypeLandscapes/CancerNormal epithelial dynamic_class3 annotation.png")
ggplot( u_dat[abs(V2)<5&abs(V3)<5], aes(y=V2, x=V3, col=  dynamic_class3)) +  geom_point(size=1.5)+scale_color_npg(name="Tumor response", labels=c("Resistant","Sensitive"))+theme_classic(base_size=20)+theme(aspect.ratio=1)+ labs(y="Immune response phenotype (Umap 2)" , x="Proliferative phenotype (Umap 1)") +facet_wrap(~Day)
ggsave(file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/PhenotypeLandscapes/CancerNormal epithelial dynamic_class3 annotation by day.png")

ggplot( u_dat[abs(V2)<5&abs(V3)<5][Day==180],aes(y=V2, x=V3,col=  dynamic_class3)) +  geom_point(size=1.5)+scale_color_npg(name="Tumor response", labels=c("Resistant","Sensitive"))+theme_classic(base_size=20)+theme(aspect.ratio=1)+ labs(y="Immune response phenotype (Umap 2)" , x="Proliferative phenotype (Umap 1)") #+facet_wrap(~Day)


ggplot( u_dat[abs(V2)<5&abs(V3)<5],aes(y=V2, x=V3, col=  EINAV_INTERFERON_SIGNATURE_IN_CANCER)) +  geom_point(size=1.5)+scale_color_viridis_c(name="Einav interferon signature \n in cancer \n response (ssGSEA score)",option="B" )+theme_classic(base_size=20)+theme(aspect.ratio=1)+ labs(y="Immune response phenotype (Umap 2)" , x="Proliferative phenotype (Umap 1)")+facet_wrap(dynamic_class3~Day)

ggplot( u_dat[abs(V2)<5&abs(V3)<5],aes(y=V2, x=V3, col=  EINAV_INTERFERON_SIGNATURE_IN_CANCER)) +  geom_point(size=1.5)+scale_color_viridis_c(name="Einav interferon signature \n in cancer \n response (ssGSEA score)",option="B" )+theme_classic(base_size=20)+theme(aspect.ratio=1)+ labs(y="Immune response phenotype (Umap 2)" , x="Proliferative phenotype (Umap 1)")#+facet_wrap(dynamic_class3~Day)
ggsave(file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/PhenotypeLandscapes/CancerNormal epithelial EINAV_INTERFERON_SIGNATURE_IN_CANCER overlay.png")



ggplot( u_dat[abs(V2)<5&abs(V3)<5],aes(y=V2, x=V3, col=  HALLMARK_INTERFERON_GAMMA_RESPONSE))   +  geom_point(size=1.5)+scale_color_viridis_c(name="Hallmark interferon gamma \n response (ssGSEA score)",option="B" )+theme_classic(base_size=20)+theme(aspect.ratio=1)+ labs(y="Immune response phenotype (Umap 2)" , x="Proliferative phenotype (Umap 1)")#+facet_wrap(dynamic_class3~Day)
ggsave(file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/PhenotypeLandscapes/CancerNormal epithelial HALLMARK_INTERFERON_GAMMA_RESPONSE overlay.png")
ggplot( u_dat[abs(V2)<5&abs(V3)<5],aes(y=V2, x=V3, col=  HALLMARK_G2M_CHECKPOINT)) +  geom_point(size=1.5)+scale_color_viridis_c(name="Hallmark G2-M checkpoint \n (ssGSEA score)",begin=0.2,option="A")+theme_classic()+theme(aspect.ratio=1)+ labs(y="Immune response phenotype (Umap 2)" , x="Proliferative phenotype (Umap 1)")#+facet_wrap(dynamic_class3~Day)
ggsave(file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/PhenotypeLandscapes/CancerNormal epithelial HALLMARK_G2M_CHECKPOINT overlay.png")
ggplot( u_dat[abs(V2)<5&abs(V3)<5],aes(y=V2, x=V3, col=  SMID_BREAST_CANCER_LUMINAL_B_UP))   +  geom_point(size=1.5)+scale_color_viridis_c(name="Smid breast cancer \n luminal B \n (ssGSEA score)",option="B" )+theme_classic(base_size=20)+theme(aspect.ratio=1)+ labs(y="Immune response phenotype (Umap 2)" , x="Proliferative phenotype (Umap 1)")#+facet_wrap(dynamic_class3~Day)
ggsave(file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/PhenotypeLandscapes/CancerNormal epithelial SMID_BREAST_CANCER_LUMINAL_B_UP overlay.png")

corVallDT[grep("BASAL",rn)]

ggplot( u_dat[abs(V2)<5&abs(V3)<5],aes(y=V2, x=V3, col=  SMID_BREAST_CANCER_BASAL_UP))   +  geom_point(size=1.5)+scale_color_viridis_c(name="Smid breast cancer \n basal \n (ssGSEA score)",option="B" )+theme_classic(base_size=20)+theme(aspect.ratio=1)+ labs(y="Immune response phenotype (Umap 2)" , x="Proliferative phenotype (Umap 1)")+facet_wrap(dynamic_class3~Day)

corVallDT[grep("JNK",rn)]
ggplot( u_dat[abs(V2)<5&abs(V3)<5],aes(y=V2, x=V3, col=  ST_JNK_MAPK_PATHWAY))   +  geom_point(size=1.5)+scale_color_viridis_c(name="JNK \n (ssGSEA score)",option="B" )+theme_classic(base_size=20)+theme(aspect.ratio=1)+ labs(y="Immune response phenotype (Umap 2)" , x="Proliferative phenotype (Umap 1)")+facet_wrap(dynamic_class3~Day)

corVallDT[grep("ESTROGEN",rn)]
ggplot( u_dat[abs(V2)<5&abs(V3)<5],aes(y=V2, x=V3, col=  ST_JNK_MAPK_PATHWAY))   +  geom_point(size=1.5)+scale_color_viridis_c(name="JNK \n (ssGSEA score)",option="B" )+theme_classic(base_size=20)+theme(aspect.ratio=1)+ labs(y="Immune response phenotype (Umap 2)" , x="Proliferative phenotype (Umap 1)")+facet_wrap(dynamic_class3~Day)


ggplot( u_dat,aes(y=V5,x=V2,group=Day,col=  HALLMARK_KRAS_SIGNALING_DN)) +  geom_point()+
  theme_classic()+facet_wrap(~Day)

ggplot( u_dat,aes(y=V1,x=V3,group=Day,col=  HALLMARK_G2M_CHECKPOINT)) +  geom_point()+
  theme_classic()+facet_wrap(~Day)
ggplot( u_dat,aes(y=V1,x=V3,group=Day,col=  HALLMARK_ESTROGEN_RESPONSE_EARLY)) +  geom_point()+
  scale_color_npg(name="Cell type")+scale_fill_npg(name="Cell type")+theme_classic()+facet_wrap(~Day)

ggplot( u_dat,aes(y=V1,x=V3,group=Day,col=  REACTOME_CELL_CYCLE_CHECKPOINTS)) +  geom_point()+
  scale_color_viridis_c()+theme_classic()


ggplot( u_dat,aes(y=V1,x=V2,group=Day,col=  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION)) +  geom_point()+
  scale_color_npg(name="Cell type")+scale_fill_npg(name="Cell type")+theme_classic()+facet_wrap(~Day)

ggplot( u_dat,aes(y=V1,x=V2,group=Day,col=  HALLMARK_INTERFERON_GAMMA_RESPONSE)) +  geom_point()+
  scale_color_npg(name="Cell type")+scale_fill_npg(name="Cell type")+theme_classic()+facet_wrap(~Day)



ggplot( u_dat,aes(y=V3,x=V4,group=Day,fill=  Celltype,col=  Celltype)) +  geom_point()+scale_color_npg(name="Cell type")+scale_fill_npg(name="Cell type")+theme_classic()+facet_wrap(~Day)
ggplot( u_dat,aes(y=V3,x=V4,group=Day,fill=  Patient.Study.ID,col=  Patient.Study.ID)) +  geom_point(size=0.5,alpha=0.5)+ theme_classic() + 
  theme(legend.position = "none")+facet_wrap(~Celltype)
ggplot( u_dat,aes(y=V3,x=V4,group=Day,fill=  dynamic_class3,col=  dynamic_class3)) +  geom_point()+scale_color_npg(name="Cell type")+scale_fill_npg(name="Cell type")+theme_classic()+facet_wrap(~Day)



### Phenotype data
load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/PhenotypesOfAllCellTypesAllArms/PhenotypesOfAllCellTypesAllArms.RData") #allphenotypes,UMAPlocs ,UMAPfiles,umapDImRedloc,umapDImRedfiles,nCellTypes,
# Encode the subtypes that have been analysed using UMAP
tmp <- data.table(allphenotypes %>% group_by(key_) %>% slice(1)) %>% dplyr::select(Celltype , Celltype_subtype) %>% unique
tmp[ , PhenoCelltype:= Celltype]
tmp[Celltype_subtype %in% c("CD4+ T cells", "Tregs"), PhenoCelltype:= "CD4+ T cells"]
tmp[Celltype_subtype %in% c("CD8+ T cells", "NK cells"), PhenoCelltype:= "CD8+ T cells"]
tmp[Celltype_subtype %in% c("Macrophages", "DC", "Monocytes"), PhenoCelltype:= "Macrophages"]
tmp[Celltype_subtype %in% c("Vas-Endo", "Lym-Endo","Endothelial cells"), PhenoCelltype:= "Endothelial cells"]
tmp[Celltype_subtype %in% c("Cancer cells"), PhenoCelltype:= "Cancer cells"]
tmp[Celltype_subtype %in% c("Normal epithelial cells"), PhenoCelltype:= "Normal epithelial cells"]
allphenotypes <- merge(allphenotypes %>% dplyr::select(-c( paste0("V", 1:5), "nCount_RNA", "nFeature_RNA", "Infercnv_CNA", "Platform","Timepoint","Sample_p_t","prop_change", "file_string", "day_fact" , "Ribo","TreatLab", "Burden_t0" ,"BurdenTracked" ,"Day0", "DayLastBurd", "Dose_lab","TreatCode","TreatCodeOrd","dynamic_class2","rgr_A","rgr_B")), tmp, by= c("Celltype", "Celltype_subtype") )


# 
# # Overlaying signal received on cancer umap
# load(file= "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Cancer gene expression wrangling/CancerSubclassCCI.RData")
# #CCIperReceiverclass
# CancerSubtypeCommCCI <- (CCIperReceiverclass %>%
#                            dplyr::select(Patient.Study.ID, Day, dynamic_class3, ARM, Pair.Name, key_, LigandPhenoCelltype, ReceptorPhenoCelltype,SumSignal,TransductionMu))
# rm(list="CCIperReceiverclass")
# CancerSubtypeCommCCI[, scaleTransduction:= exp(scale(log(SumSignal), center=F)), by=c("Pair.Name")] 
# 
# CancersummaryInteraction <- data.table(CancerSubtypeCommCCI%>%
#                                          group_by(Patient.Study.ID,LigandPhenoCelltype,ReceptorPhenoCelltype,Day,dynamic_class3,ARM,key_) %>%
#                                          dplyr::summarise(muln_scaleTransduction=mean(log(1+scaleTransduction))))
# rm(list="CancerSubtypeCommCCI")
# 
# 
# ##Norm epith overlay
# load(file= "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Cancer gene expression wrangling/Normal epithelial SubclassCCI.RData")
# #CCIperReceiverclass
# NormEpSubtypeCommCCI <- (CCIperReceiverclass %>%
#                            dplyr::select(Patient.Study.ID, Day, dynamic_class3, ARM, Pair.Name, key_, LigandPhenoCelltype, ReceptorPhenoCelltype,SumSignal,TransductionMu))
# rm(list="CCIperReceiverclass")
# NormEpSubtypeCommCCI[, scaleTransduction:= exp(scale(log(SumSignal), center=F)), by=c("Pair.Name")] 
# 
# NormEpsummaryInteraction <- data.table(NormEpSubtypeCommCCI%>%
#                                          group_by(Patient.Study.ID,LigandPhenoCelltype,ReceptorPhenoCelltype,Day,dynamic_class3,ARM,key_) %>%
#                                          dplyr::summarise(muln_scaleTransduction=mean(log(1+scaleTransduction))))
# rm(list="NormEpSubtypeCommCCI")
# 
# 
# 
# # reshape
# CancersummaryInteraction2 <- data.table( rbind(CancersummaryInteraction,NormEpsummaryInteraction) %>%
#                                            dplyr::select(Patient.Study.ID, Day, dynamic_class3, ARM,  key_, LigandPhenoCelltype, ReceptorPhenoCelltype, muln_scaleTransduction )%>%
#                                            spread(LigandPhenoCelltype,muln_scaleTransduction) )
# 
# setnames(CancersummaryInteraction2,old=c("Adipocytes", "B cells", "Cancer cells", "CD4+ T cells", "CD8+ T cells", "Endothelial cells", "Fibroblasts", "Macrophages", "Normal epithelial cells", "Pericytes"),
#          new=c("FromAdipocytes", "FromBcells", "FromCancercells", "FromCD4Tcells", "FromCD8Tcells", "FromEndothelialcells", "FromFibroblasts", "FromMacrophages", "FromNormalepithelialcells", "FromPericytes"))
# 
# CancerNormalsummaryInteractionTotal <- CancersummaryInteraction2
#save(CancerNormalsummaryInteractionTotal,    file= "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Cancer gene expression wrangling/Cancer Normal epithelial summaryInteractionTotal.RData")


load(    file= "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Cancer gene expression wrangling/Cancer Normal epithelial summaryInteractionTotal.RData")
#CancerNormalsummaryInteractionTotal

# allphenotypes%>%dplyr::select("Cell.ID","key_", "Day","ARM", "dynamic_class3")
# u_dat%>%dplyr::select("Cell.ID","key_", "Day","ARM", "dynamic_class3")  

u_datmerg<- merge(u_dat[abs(V2)<5&abs(V3)<5],allphenotypes[Celltype%in%c("Cancer cells","Normal epithelial cells")]%>%dplyr::select("Cell.ID","key_", "Day","ARM", "dynamic_class3") , by="Cell.ID")
rm(list="u_dat")
rm(list="SSGSEA")
rm(list="umap_data_tt")
u_datmerg[,Day:=Day.x]
u_datmerg[,Day.x:=NULL]
u_datmerg[,Day.y:=NULL]
udataComm <- merge( u_datmerg, CancerNormalsummaryInteractionTotal , by=c("key_", "Patient.Study.ID" ,"Day")) 

udataComm[1:10,1:10]

#save(udataComm, file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Cancer gene expression wrangling/cancer and normal Communication on ssgsea landscape.RData")


immuneplotdd <- udataComm[ARM != "A"] 
immuneplotdd[is.na(FromCD8Tcells), FromCD8Tcells:= 0]
#immuneplotdd[is.na(FromFibroblasts), FromFibroblasts:= 0]

# ggplot(immuneplotdd,aes(y=V2, x=V3, col=  log(1+FromFibroblasts) )) +  geom_point(size=1.5)+theme_classic(base_size=20)+theme(aspect.ratio=1)+ labs(y="Immune response phenotype (Umap 2)" , x="Proliferative phenotype (Umap 1)")+
#   scale_color_viridis_c(name="Fibroblasts communication \n to cancer and normal epithelial cells ",option="B" )+facet_wrap(dynamic_class3~Day)
#
# ggplot(immuneplotdd,aes(y=V2, x=V3, col=  log(1+FromMacrophages) )) +  geom_point(size=1.5)+theme_classic(base_size=20)+theme(aspect.ratio=1)+ labs(y="Immune response phenotype (Umap 2)" , x="Proliferative phenotype (Umap 1)")+
#   scale_color_viridis_c(name="Macrophages communication \n to cancer and normal epithelial cells ",option="B" )+facet_wrap(dynamic_class3~Day)


ggplot(immuneplotdd,aes(y=V2, x=V3, col=  log(1+FromCD8Tcells) )) +  geom_point(size=1.5)+theme_classic(base_size=20)+theme(aspect.ratio=1)+ labs(y="Immune response phenotype (Umap 2)" , x="Proliferative phenotype (Umap 1)")+
  scale_color_viridis_c(name="CD8 Tcells communication \n to cancer and normal epithelial cells ",option="B" )+facet_wrap(dynamic_class3~Day)
ggsave(file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/PhenotypeLandscapes/CancerNormal epithelial CD8 Tcells communication overlay.png")



ggplot(immuneplotdd,aes(y=V2, x=V3, col=  log(1+FromCD8Tcells) )) +  geom_point(size=1.5)+theme_classic(base_size=20)+theme(aspect.ratio=1)+ labs(y="Immune response phenotype (Umap 2)" , x="Proliferative phenotype (Umap 1)")+
  scale_color_viridis_c(name="CD8 Tcells communication \n to cancer and normal epithelial cells ",option="B" )+facet_wrap(~dynamic_class3)
ggsave(file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/PhenotypeLandscapes/CancerNormal epithelial CD8 Tcells communication overlay mixed Days.png")


ggplot(immuneplotdd[Day==180],aes(y=V2, x=V3, col=  log(1+FromCD8Tcells) )) +  geom_point(size=1.5)+theme_classic(base_size=20)+theme(aspect.ratio=1)+ labs(y="Immune response phenotype (Umap 2)" , x="Proliferative phenotype (Umap 1)")+
  scale_color_viridis_c(name="CD8 Tcells communication \n to cancer and normal epithelial cells ",option="B" )+facet_wrap(~dynamic_class3)
ggsave(file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/PhenotypeLandscapes/CancerNormal epithelial CD8 Tcells communication overlay Day180.png")

ggplot(immuneplotdd[Celltype=="Cancer cells"][],aes(y=log(1+FromCD8Tcells), x=Day,group=Day,  col=dynamic_class3 )) + 
  geom_violin()+theme_classic(base_size=20)+theme(aspect.ratio=1)+ 
  facet_wrap(~dynamic_class3)

ggplot(immuneplotdd[Celltype=="Cancer cells"],aes(y=HALLMARK_INTERFERON_GAMMA_RESPONSE, x=log(1+Day),group=interaction(dynamic_class3,Day),  col=dynamic_class3,  fill=dynamic_class3  )) + 
  geom_violin()+theme_classic(base_size=20)+theme(aspect.ratio=1)#+ 
#labs(y="Immune response phenotype (Umap 2)" , x="Proliferative phenotype (Umap 1)")+
#facet_wrap(~dynamic_class3)

ggplot(immuneplotdd[Celltype=="Cancer cells"],aes(y=(sqrt(FromCD8Tcells) ), x=log(1+Day),group=interaction(dynamic_class3,Day),  col=dynamic_class3,  fill=dynamic_class3  )) + 
  geom_violin()+theme_classic(base_size=20)+theme(aspect.ratio=1)#+ 

immuneplotdd2 <- immuneplotdd[Celltype=="Cancer cells"][]
immuneplotdd2[ ,HALLMARK_INTERFERON_GAMMA_RESPONSEcat:= round( 100 * HALLMARK_INTERFERON_GAMMA_RESPONSE)/100]


ggplot(data.table(immuneplotdd2 %>% group_by(Patient.Study.ID,Day,dynamic_class3)%>%summarise(lnFromCD8Tcells = median(sqrt(FromCD8Tcells) ),
                                                                                              HALLMARK_INTERFERON_GAMMA_RESPONSE=median(HALLMARK_INTERFERON_GAMMA_RESPONSE),
                                                                                              N=n() )),
       aes(y=HALLMARK_INTERFERON_GAMMA_RESPONSE, x=lnFromCD8Tcells,col=dynamic_class3)) +  
  geom_point(aes(col=dynamic_class3),size=2.5)+theme_classic(base_size=20)+theme(aspect.ratio=1)#+ 
#labs(y="Immune response phenotype (Umap 2)" , x="Proliferative phenotype (Umap 1)") +
#geom_smooth(method="gam")#+
#facet_wrap(~Day)


ggplot(data.table(immuneplotdd2%>%group_by(Patient.Study.ID,Day,dynamic_class3)%>%summarise(lnFromCD8Tcells = median(sqrt(FromCD8Tcells) ),
                                                                                            HALLMARK_INTERFERON_GAMMA_RESPONSE=median(HALLMARK_INTERFERON_GAMMA_RESPONSE),
                                                                                            ST_JNK_MAPK_PATHWAY=median(ST_JNK_MAPK_PATHWAY),
                                                                                            N=n() )),
       aes(x=HALLMARK_INTERFERON_GAMMA_RESPONSE, y=lnFromCD8Tcells,group=dynamic_class3)) +  
  geom_point(aes(col=ST_JNK_MAPK_PATHWAY))+theme_classic(base_size=20)+theme(aspect.ratio=1)+ 
  #labs(y="Immune response phenotype (Umap 2)" , x="Proliferative phenotype (Umap 1)")+
  geom_smooth(method="gam")#+facet_wrap(~Day)

# ggplot(data.table(immuneplotdd2[Day==180]%>%group_by(HALLMARK_INTERFERON_GAMMA_RESPONSEcat,dynamic_class3)%>%summarise(lnFromCD8Tcells = mean(log(1+FromCD8Tcells) ),
#                                                                                    lcl = mean(log(1+FromCD8Tcells) ) - 1.96*sd(log(1+FromCD8Tcells) ),
#                                                                                    ucl = mean(log(1+FromCD8Tcells) ) + 1.96*sd(log(1+FromCD8Tcells) ),
#                                                                                    N=n() ))[N>3],
#        aes(x=HALLMARK_INTERFERON_GAMMA_RESPONSEcat, y=lnFromCD8Tcells,col=dynamic_class3 )) +  geom_point(size=1.5)+theme_classic(base_size=20)+theme(aspect.ratio=1)+ 
#   #labs(y="Immune response phenotype (Umap 2)" , x="Proliferative phenotype (Umap 1)")+
#   geom_smooth(method="gam")+#+facet_wrap(~Day)+
# geom_errorbar(aes(ymax=ucl,ymin=lcl))
#ggsave(file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/PhenotypeLandscapes/CancerNormal epithelial CD8 Tcells communication vs HALLMARK_INTERFERON_GAMMA_RESPONSE.png")

ggplot(immuneplotdd2[Day==180]
       ,aes(y=HALLMARK_INTERFERON_GAMMA_RESPONSE, x=log(1+FromCD8Tcells))) +  
  geom_jitter(aes (col=dynamic_class3),width=0.02,size=1.5)+theme_classic(base_size=20)+theme(aspect.ratio=1)+ 
  #labs(y="Immune response phenotype (Umap 2)" , x="Proliferative phenotype (Umap 1)")+
  geom_smooth(method="lm")#+facet_wrap(~Day)
ggsave(file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/PhenotypeLandscapes/CancerNormal epithelial CD8 Tcells communication vs HALLMARK_INTERFERON_GAMMA_RESPONSE Day 180.png")

ggplot(immuneplotdd2[Day==0]
       ,aes(y=HALLMARK_INTERFERON_GAMMA_RESPONSE, x=log(1+FromCD8Tcells),col=dynamic_class3)) +  
  geom_point(size=1.5)+theme_classic(base_size=20)+theme(aspect.ratio=1)+ 
  #labs(y="Immune response phenotype (Umap 2)" , x="Proliferative phenotype (Umap 1)")+
  geom_smooth(method="lm")+facet_wrap(~Day)
ggsave(file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/PhenotypeLandscapes/CancerNormal epithelial CD8 Tcells communication vs HALLMARK_INTERFERON_GAMMA_RESPONSE Day 0.png")


ggplot(immuneplotdd2[]
       ,aes(y=HALLMARK_INTERFERON_GAMMA_RESPONSE, x=log(1+FromCD8Tcells),col=dynamic_class3)) +  
  geom_point(size=1.5)+theme_classic(base_size=20)+theme(aspect.ratio=1)+ 
  #labs(y="Immune response phenotype (Umap 2)" , x="Proliferative phenotype (Umap 1)")+
  geom_smooth(method="lm")+facet_wrap(~Day)
ggsave(file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/PhenotypeLandscapes/CancerNormal epithelial CD8 Tcells communication vs HALLMARK_INTERFERON_GAMMA_RESPONSE Day all.png")



ggplot( immuneplotdd, aes(V1, V4, col= log(1+FromCD8Tcells)  ) ) +  geom_point(size=0.5) + theme_classic(base_size=18) + theme(aspect.ratio=1)+
  labs(y="Umap B" , x="Umap A")+scale_color_viridis(option="B")+facet_wrap(dynamic_class3~Day)








#ggplot( u_dat[abs(V2)<5&abs(V3)<5],aes(y=V2, x=V3, col=  HALLMARK_ANGIOGENESIS)) +  geom_point(size=1.5)+scale_color_viridis_c(option="D")+theme_classic()+theme(aspect.ratio=1)+ labs(y="Immune response phenotype (Umap 2)" , x="Proliferative phenotype (Umap 1)")#+facet_wrap(dynamic_class3~Day)



ggplot( u_dat,aes(y=EINAV_INTERFERON_SIGNATURE_IN_CANCER,x=HALLMARK_ANGIOGENESIS,group=Day, col=  dynamic_class3)) +  geom_point()+scale_color_npg( )+theme_classic()+facet_wrap(dynamic_class3~Day)


ggplot( u_dat,aes(y=V3,x=V4,group=Day,col=  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION)) +  geom_point()+scale_color_viridis_c()+theme_classic()
ggplot( u_dat,aes(y=V3,x=V4,group=Day,col=  HALLMARK_ESTROGEN_RESPONSE_EARLY)) +  geom_point()+scale_color_viridis_c()+theme_classic()
ggplot( u_dat,aes(y=V3,x=V4,group=Day,col=  HALLMARK_INTERFERON_GAMMA_RESPONSE)) +  geom_point()+scale_color_viridis_c()+theme_classic()
ggplot( u_dat,aes(y=V3,x=V4,group=Day,col=  BOQUEST_STEM_CELL_UP)) +  geom_point()+scale_color_viridis_c()+theme_classic()
ggplot( u_dat,aes(y=V3,x=V4,group=Day,col=  SMID_BREAST_CANCER_BASAL_UP)) +  geom_point()+scale_color_viridis_c()+theme_classic()
ggplot( u_dat,aes(y=V3,x=V4,group=Day,col=  HALLMARK_MYOGENESIS)) +  geom_point()+scale_color_viridis_c()+theme_classic()
ggplot( u_dat,aes(y=V3,x=V4,group=Day,col=  BOQUEST_STEM_CELL_UP)) +  geom_point()+scale_color_viridis_c()+theme_classic()

ggplot( u_dat,aes(y=V3,x=V4,group=Day,col=  WHITFIELD_CELL_CYCLE_G2)) +  geom_point()+
  scale_color_viridis_c()+theme_classic() + facet_wrap(~Day)

ggplot( u_dat,aes(y=V3,x=V4,group=Day,col=  HALLMARK_ESTROGEN_RESPONSE_EARLY)) +  geom_point()+
  scale_color_viridis_c()+theme_classic() + facet_wrap(~Day)

ggplot( u_dat[ARM!="A"][Celltype!="Cancer cells"][V3< -2.5],aes(y=V3,x=V4,group=Day,col=  dynamic_class3)) +  geom_point(size=0.3)+
  scale_color_npg()+theme_classic() + facet_wrap(~Day)





ggplot( u_dat[ARM!="A"][Celltype=="Cancer cells"],aes(y=V1,x=log(1+Day),group=Day,fill=  dynamic_class3,col=  dynamic_class3)) +  geom_violin()+scale_color_npg(name="Tumor response")+scale_fill_npg(name="Tumor response")+theme_classic()+facet_wrap(~dynamic_class3)
ggplot( u_dat[ARM!="A"][Celltype=="Cancer cells"],aes(y=V2,x=log(1+Day),group=Day,fill=  dynamic_class3,col=  dynamic_class3)) +  geom_violin()+scale_color_npg(name="Tumor response")+scale_fill_npg(name="Tumor response")+theme_classic()+facet_wrap(~dynamic_class3)
ggplot( u_dat[ARM!="A"][Celltype=="Cancer cells"],aes(y=V3,x=log(1+Day),group=Day,fill=  dynamic_class3,col=  dynamic_class3)) +  geom_violin()+scale_color_npg(name="Tumor response")+scale_fill_npg(name="Tumor response")+theme_classic()+facet_wrap(~dynamic_class3)
ggplot( u_dat[ARM!="A"][Celltype=="Cancer cells"],aes(y=V4,x=log(1+Day),group=Day,fill=  dynamic_class3,col=  dynamic_class3)) +  geom_violin()+scale_color_npg(name="Tumor response")+scale_fill_npg(name="Tumor response")+theme_classic()+facet_wrap(~dynamic_class3)
ggplot( u_dat[ARM!="A"][Celltype=="Cancer cells"],aes(y=V5,x=log(1+Day),group=Day,fill=  dynamic_class3,col=  dynamic_class3)) +  geom_violin()+scale_color_npg(name="Tumor response")+scale_fill_npg(name="Tumor response")+theme_classic()+facet_wrap(~dynamic_class3)

ggplot( u_dat[ARM!="A"][Celltype=="Cancer cells"]%>%group_by(Day,dynamic_class3,Patient.Study.ID,Celltype)%>%summarise(V1=median(V1)),aes(y=V1,x=log(1+Day),group=Day,fill=  dynamic_class3,col=  dynamic_class3)) +  geom_violin()+scale_color_npg(name="Tumor response")+scale_fill_npg(name="Tumor response")+theme_classic()+facet_wrap(~dynamic_class3)
ggplot( u_dat[ARM!="A"][Celltype=="Cancer cells"]%>%group_by(Day,dynamic_class3,Patient.Study.ID,Celltype)%>%summarise(V2=median(V2)),aes(y=V2,x=log(1+Day),group=Day,fill=  dynamic_class3,col=  dynamic_class3)) +  geom_violin()+scale_color_npg(name="Tumor response")+scale_fill_npg(name="Tumor response")+theme_classic()+facet_wrap(~dynamic_class3)
ggplot( u_dat[ARM!="A"][Celltype=="Cancer cells"]%>%group_by(Day,dynamic_class3,Patient.Study.ID,Celltype)%>%summarise(V3=median(V3)),aes(y=V3,x=log(1+Day),group=Day,fill=  dynamic_class3,col=  dynamic_class3)) +  geom_violin()+scale_color_npg(name="Tumor response")+scale_fill_npg(name="Tumor response")+theme_classic()+facet_wrap(~dynamic_class3)
ggplot( u_dat[ARM!="A"][Celltype=="Cancer cells"]%>%group_by(Day,dynamic_class3,Patient.Study.ID,Celltype)%>%summarise(V4=median(V4)),aes(y=V4,x=log(1+Day),group=Day,fill=  dynamic_class3,col=  dynamic_class3)) +  geom_violin()+scale_color_npg(name="Tumor response")+scale_fill_npg(name="Tumor response")+theme_classic()+facet_wrap(~dynamic_class3)
ggplot( u_dat[ARM!="A"][Celltype=="Cancer cells"]%>%group_by(Day,dynamic_class3,Patient.Study.ID,Celltype)%>%summarise(V5=median(V5)),aes(y=V5,x=log(1+Day),group=Day,fill=  dynamic_class3,col=  dynamic_class3)) +  geom_violin()+scale_color_npg(name="Tumor response")+scale_fill_npg(name="Tumor response")+theme_classic()+facet_wrap(~dynamic_class3)

ggplot( u_dat[ARM!="A"][Celltype=="Cancer cells"]%>%group_by(Day,dynamic_class3,Patient.Study.ID,Celltype)%>%summarise(V5=median(V5)),
        aes(y=V5,x=log(1+Day),group=Day,fill=  RGR,col=  RGR)) +  geom_violin()+scale_color_npg(name="Tumor response")+scale_fill_npg(name="Tumor response")+theme_classic()+facet_wrap(~dynamic_class3)


ggplot( u_dat[ARM!="A"][Celltype=="Cancer cells"][Day==180], aes(y=V5,x=V2, group=Day,fill=  dynamic_class3,col=  V4)) +  geom_point() +facet_wrap(~dynamic_class3)
ggplot( u_dat[ARM!="A"][Celltype=="Cancer cells"], aes(y=V5,x=V2, group=Day,fill=  dynamic_class3,col=  dynamic_class3)) +  geom_point()+scale_color_npg(name="dynamic_class3")+scale_fill_npg(name="dynamic_class3")+theme_classic()+facet_wrap(dynamic_class3~Day)


corVall <- cor(u_dat[ARM!="A"][Celltype=="Cancer cells"]%>%dplyr::select(paste0("V",1:NN)),u_dat[ARM!="A"][Celltype=="Cancer cells"]%>%dplyr::select(pathwaysNamesUsed),method=methodcor)#corV1<-cor(u_dat$V1,full_dd_transp)

corVallDT <- data.table(t(corVall), keep.rownames = T)
corVallDT[(order(-abs(V1)))][1:10]
corVallDT[(order(-abs(V2)))][1:30]
corVallDT[(order(-abs(V3)))][1:20]
corVallDT[(order(-abs(V4)))][1:20]
corVallDT[(order(-abs(V5)))][1:30]

corVallDT[grep("HALLMARK",rn)][(order(-abs(V2)))][1:30]
corVallDT[grep("HALLMARK",rn)][(order(-abs(V3)))][1:30]

corVallDT[grep("HALLMARK",rn)][(order(-abs(V1)))][1:20] #HALLMARK_PI3K_AKT_MTOR_SIGNALING, HALLMARK_MYOGENESIS, HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION ;
corVallDT[grep("HALLMARK",rn)][(order(-abs(V2)))][1:20] #HALLMARK_INTERFERON_GAMMA_RESPONSE, HALLMARK_INTERFERON_ALPHA_RESPONSE, HALLMARK_APOPTOSIS
corVallDT[grep("HALLMARK",rn)][(order(-abs(V3)))][1:20] #HALLMARK_MYOGENESIS, HALLMARK_E2F_TARGETS, HALLMARK_G2M_CHECKPOINT; 
corVallDT[grep("HALLMARK",rn)][(order(-abs(V4)))][1:20] #HALLMARK_HEDGEHOG_SIGNALING,HALLMARK_XENOBIOTIC_METABOLISM, HALLMARK_COAGULATION,HALLMARK_ESTROGEN_RESPONSE_LATE ; 
corVallDT[grep("HALLMARK",rn)][(order(-abs(V5)))][1:20] #HALLMARK_E2F_TARGETS , HALLMARK_TGF_BETA_SIGNALING,  HALLMARK_ESTROGEN_RESPONSE_EARLY


ggplot( u_dat[ARM!="A"][Celltype=="Cancer cells"],
        aes(y=V5,x=V2,group=Day,col=  HALLMARK_TGF_BETA_SIGNALING)) +  
  geom_point()+scale_color_viridis_c()+
  theme_classic()+facet_wrap(dynamic_class3~Day)

ggplot( u_dat[ARM!="A"][Celltype=="Cancer cells"],
        aes(y=V5,x=V2,group=Day,col=  HALLMARK_ESTROGEN_RESPONSE_EARLY)) +  
  geom_point()+scale_color_viridis_c()+
  theme_classic()+facet_wrap(dynamic_class3~Day)


ggplot( u_dat[ARM!="A"][Celltype=="Cancer cells"],
        aes(y=V5,x=V2,group=Day,col=  sqrt(HALLMARK_INTERFERON_GAMMA_RESPONSE- min( HALLMARK_INTERFERON_GAMMA_RESPONSE)) ) ) +  
  geom_point()+scale_color_viridis_c(option="B", name="name")+
  theme_classic()+facet_wrap(dynamic_class3~Day)


ggplot( u_dat[ARM!="A"][Celltype=="Cancer cells"][Day==180],
        aes(y=V5,x=V2,group=Day,col=  (FARMER_BREAST_CANCER_CLUSTER_7- min( FARMER_BREAST_CANCER_CLUSTER_7)) ) ) +  
  geom_point()+scale_color_viridis_c(option="B", name="name")+
  theme_classic()+facet_wrap(dynamic_class3~Day)


ggplot( u_dat[ARM!="A"][Celltype=="Cancer cells"],
        aes(y=V5,x=V2,group=Day,col=  RIGGINS_TAMOXIFEN_RESISTANCE_DN)) +  
  geom_point()+scale_color_viridis_c()+
  theme_classic()+facet_wrap(dynamic_class3~Day)

ggplot( u_dat[ARM!="A"][Celltype=="Cancer cells"],
        aes(y=V5,x=V2,group=Day,col=  HALLMARK_INTERFERON_GAMMA_RESPONSE)) +  
  geom_point()+scale_color_viridis_c()+
  theme_classic()+facet_wrap(dynamic_class3~Day)



ggplot( u_dat[ARM!="A"][Celltype=="Cancer cells"],
        aes(y=V5,x=V2,group=Day,col=  AZARE_NEOPLASTIC_TRANSFORMATION_BY_STAT3_UP)) +  
  geom_point()+scale_color_viridis_c()+
  theme_classic()+facet_wrap(dynamic_class3~Day)



u_dat[ARM!="A"][Celltype=="Cancer cells"] 


ggplot( u_dat[ARM!="A"][Celltype=="Cancer cells"],aes(y=V5,x=V2,group=Day,fill=  SMID_BREAST_CANCER_BASAL_UP,col=  SMID_BREAST_CANCER_BASAL_UP)) +  geom_point()+scale_color_viridis_c(name="SMID_BREAST_CANCER_BASAL_DN")+scale_fill_viridis_c(name="SMID_BREAST_CANCER_BASAL_DN")+theme_classic()+facet_wrap(dynamic_class3~Day)
ggplot( u_dat[ARM!="A"][Celltype=="Cancer cells"],aes(y=V5,x=V2,group=Day,fill=  ST_JNK_MAPK_PATHWAY,col=  ST_JNK_MAPK_PATHWAY)) +  geom_point()+scale_color_viridis_c(name="ST_JNK_MAPK_PATHWAY")+scale_fill_viridis_c(name="ST_JNK_MAPK_PATHWAY")+theme_classic()+facet_wrap(dynamic_class3~Day)
ggplot( u_dat[ARM!="A"][Celltype=="Cancer cells"],aes(y=V5,x=V2,group=Day,fill=  KOYAMA_SEMA3B_TARGETS_DN,col=  KOYAMA_SEMA3B_TARGETS_DN)) +  geom_point()+scale_color_viridis_c(name="KOYAMA_SEMA3B_TARGETS")+scale_fill_viridis_c(name="KOYAMA_SEMA3B_TARGETS_DN")+theme_classic()+facet_wrap(dynamic_class3~Day)

ggplot( u_dat[ARM!="A"][Celltype=="Cancer cells"],aes(y=V5,x=V2,group=Day,fill=  -XU_GH1_AUTOCRINE_TARGETS_DN,col=  -XU_GH1_AUTOCRINE_TARGETS_DN)) +  geom_point()+scale_color_viridis_c(name="XU_GH1_AUTOCRINE_TARGETS_DN")+scale_fill_viridis_c(name="XU_GH1_AUTOCRINE_TARGETS_DN")+theme_classic()+facet_wrap(dynamic_class3~Day)


ggplot( u_dat[ARM!="A"][Celltype=="Cancer cells"],aes(y=V3,x=V4,group=Day,fill=  SMID_BREAST_CANCER_BASAL_DN,col=  SMID_BREAST_CANCER_BASAL_DN)) +  geom_point()+scale_color_viridis_c(name="SMID_BREAST_CANCER_BASAL_DN")+scale_fill_viridis_c(name="SMID_BREAST_CANCER_BASAL_DN")+theme_classic()+facet_wrap(dynamic_class3~Day)
ggplot( u_dat[ARM!="A"][Celltype=="Cancer cells"],aes(y=V3,x=V4,group=Day,fill=  KEGG_MAPK_SIGNALING_PATHWAY,col=  KEGG_MAPK_SIGNALING_PATHWAY)) +  geom_point()+scale_color_viridis_c(name="KEGG_MAPK_SIGNALING_PATHWAY")+scale_fill_viridis_c(name="KEGG_MAPK_SIGNALING_PATHWAY")+theme_classic()+facet_wrap(dynamic_class3~Day)
ggplot( u_dat[ARM!="A"][Celltype=="Cancer cells"],aes(y=V3,x=V4,group=Day,fill=  ST_JNK_MAPK_PATHWAY,col=  ST_JNK_MAPK_PATHWAY)) +  geom_point()+scale_color_viridis_c(name="ST_JNK_MAPK_PATHWAY")+scale_fill_viridis_c(name="ST_JNK_MAPK_PATHWAY")+theme_classic()+facet_wrap(dynamic_class3~Day)

ggplot( u_dat[ARM!="A"][Celltype=="Cancer cells"],aes(y=V3,x=V4,group=Day,fill=  HALLMARK_PI3K_AKT_MTOR_SIGNALING,col=  HALLMARK_PI3K_AKT_MTOR_SIGNALING)) +  geom_point()+scale_color_viridis_c(name="HALLMARK_PI3K_AKT_MTOR_SIGNALING")+scale_fill_viridis_c(name="HALLMARK_PI3K_AKT_MTOR_SIGNALING")+theme_classic()+facet_wrap(dynamic_class3~Day)

ggplot( u_dat[ARM!="A"][Celltype=="Cancer cells"],aes(y=V5,x=V4,group=Day,fill=  XU_GH1_AUTOCRINE_TARGETS_DN,col=  XU_GH1_AUTOCRINE_TARGETS_DN)) +  geom_point()+scale_color_viridis_c(name="XU_GH1_AUTOCRINE_TARGETS_DN")+scale_fill_viridis_c(name="XU_GH1_AUTOCRINE_TARGETS_DN")+theme_classic()+facet_wrap(dynamic_class3~Day)



ggplot( u_dat,aes(y=V4,x=V2,group=Day,fill=  HALLMARK_INTERFERON_GAMMA_RESPONSE,col=  HALLMARK_INTERFERON_GAMMA_RESPONSE)) +  geom_point()+scale_color_viridis_c(name="HALLMARK_INTERFERON_GAMMA_RESPONSE")+scale_fill_viridis_c(name="HALLMARK_INTERFERON_GAMMA_RESPONSE")+theme_classic()+facet_wrap(dynamic_class3~Day)
