rm(list=ls())
require(data.table)
require(dplyr)
require(ggplot2)
require(tidyr)
require(ggsci)
require(umap)
load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Cancer gene expression wrangling/Cancer summary gene expression.RData") #gene_summary
expressedGenes <- gene_summary[][coverage>0.01]$GeneName

load( file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Gene expression landscapes/cancer gene expression landscape/cancer gene expression raw and summarised.RData")
start_enddd <-  u_dattrim[hasStartEnd==1]  
rm(list="u_dattrim")

lu1 <- unique( data.table(start_enddd%>% dplyr::select(Patient.Study.ID,Day, Treatment) ) ) #lu1[Day==0][Treatment=="letrozole"]
lu2 <- unique( data.table(start_enddd%>% dplyr::select(Patient.Study.ID, Treatment) ) )

#for each patient select cancer cell data (across timepoints but within one tumor) and perform DE analysis using all genes with >1% coverage
allDE <- rbindlist( lapply(1:nrow(lu2), function(i){
  cat("Tumor    ") ;  cat(i) ; cat("\n")
  
  dd_i <- start_enddd[ Patient.Study.ID==lu2[i]$Patient.Study.ID]   #dd_i %>% dim()   #ggplot(dd_i, aes(y=sqrt(JUN), x=Day , group=Day)  ) + geom_boxplot() + geom_point() #m1 <- lm(formula= paste0("sqrt(","JUN",") ~ as.factor(Day)") ,    dd_i )
  
  # Perform pathway analysis 
  pathanalysis0 <- rbindlist(lapply( 1:length(expressedGenes), function(pp){
    cat(pp)
    res <- tryCatch({
      m1 <- lm(formula= paste0("log(1+",expressedGenes[pp],") ~ as.factor(Day)") ,    dd_i )
      data.table(gene=expressedGenes[pp],  coef( summary(m1) ), keep.rownames = T)
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
# save output
#save(allDE,  file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Gene expression landscapes/cancer gene expression landscape/Tumor specific shift strt vs end cancer gene expression scrna.RData")
load(  file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Gene expression landscapes/cancer gene expression landscape/Tumor specific shift strt vs end cancer gene expression scrna.RData")
allDE[rn!="(Intercept)"][gene=="MYC"]$Estimate%>%hist()
allDE[rn!="(Intercept)"][Treatment=="letrozole"][Estimate!=0][adj.pval<0.05][order(-Estimate)][1:30]

allDE[rn!="(Intercept)"][Treatment=="letrozole"][abs(Estimate)>log(2)][adj.pval<0.05][order(-Estimate)][1:30]
ggplot(allDE[rn!="(Intercept)"][Treatment=="letrozole"][abs(Estimate)>log(2)][adj.pval<0.05] , aes(x=Estimate, y=log(1e-10+adj.pval)))+geom_point()

allDE[rn!="(Intercept)"][Treatment=="letrozole"][Estimate!=0][adj.pval<0.05][order(-Estimate)][1:30]

ggplot(allDE[rn!="(Intercept)"][Treatment=="letrozole"][Estimate!=0][adj.pval<0.05] , aes(x=Estimate, y=log(1e-10+adj.pval) , col=abs(Estimate)>log(1.5) ))+geom_point()

summary_allDE <- data.table(allDE[rn!="(Intercept)"][Treatment=="letrozole"][Estimate!=0][adj.pval<0.05]%>%group_by(gene) %>%
  summarise(nsignif=length(Patient.Study.ID) , 
            mean_effect=sum(Estimate)/ (length( unique(allDE[Treatment=="letrozole"]$Patient.Study.ID) ) ) ))
#save(summary_allDE,  file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Gene expression landscapes/cancer gene expression landscape/Tumor specific shift strt vs end cancer gene expression scrna summary.RData")
load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Gene expression landscapes/cancer gene expression landscape/Tumor specific shift strt vs end cancer gene expression scrna summary.RData")

summary_allDE2 <- data.table(allDE[rn!="(Intercept)"][Treatment=="letrozole"][abs(Estimate)>log(1.5)][adj.pval<0.05]%>%group_by(gene) %>%
                               summarise(nsignif=length(Patient.Study.ID) , 
                                         mean_effect=sum(Estimate)/ (length( unique(allDE[Treatment=="letrozole"]$Patient.Study.ID) ) ) ))


ggplot(summary_allDE , aes(x=mean_effect, y=nsignif ) ) + geom_point()
ggplot(summary_allDE[nsignif>=6][order(-mean_effect)] , aes(x=mean_effect, y=nsignif )) + geom_jitter(width=0,height=0.1)
summary_allDE[nsignif>6] [abs(mean_effect)>0.5][order(-abs(mean_effect))]
allDE[rn!="(Intercept)"][Treatment=="letrozole"][gene=="CDKN1A"]

ggplot( allDE[rn!="(Intercept)"][Treatment=="letrozole"][Estimate!=0][
  gene%in%summary_allDE[nsignif>4]$gene][
  gene%in%summary_allDE[order(-(mean_effect))] [1:20] $gene   ][order(gene)] , 
        aes(Estimate,gene))+geom_jitter(width=0,height=0.1)

#write.csv(summary_allDE,  file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Gene expression landscapes/cancer gene expression landscape/Tumor shift strt vs end Endocrine therapy cancer gene expression scrna summary.csv")
allDE[Treatment=="letrozole"]$Patient.Study.ID%>%unique()


sort( allDE[rn!="(Intercept)"][Treatment=="letrozole"][adj.pval<0.05][abs(Estimate) > log(2)]$gene%>%table() ,decreasing =T)[1:20]
sort( allDE[rn!="(Intercept)"][Treatment=="letrozole"][adj.pval<0.05][(Estimate) > log(2)]$gene%>%table() ,decreasing =T)[1:20]
sort( allDE[rn!="(Intercept)"][][adj.pval<0.05][(Estimate) > log2(1)]$gene%>%table() ,decreasing =T)[1:20]

# VTCN1
sort( allDE[rn!="(Intercept)"][Treatment=="letrozole"][adj.pval<0.05][(Estimate) < log(2)]$gene%>%table() ,decreasing =T)[1:20]

topN <- 25
consistentUpreg <- summary_allDE[nsignif> (length(unique(allDE[Treatment=="letrozole"]$Patient.Study.ID)))/2 ][mean_effect>log(1.5)][order(-abs(mean_effect))] [1:topN] 
consistentUpreg[,dir:="up"]
consistentDownreg <- summary_allDE[nsignif> (length(unique(allDE[Treatment=="letrozole"]$Patient.Study.ID)))/2 ][mean_effect<log(1.5)][order((mean_effect))] [topN:1] 
consistentDownreg[,dir:="down"]
consistentgenesUDwn <- na.omit(rbind(consistentUpreg,consistentDownreg))

pltthis1 <- allDE[rn!="(Intercept)"][  gene%in% consistentgenesUDwn $gene] #[Treatment=="letrozole"]
pltthis1$gene <- factor(pltthis1$gene , levels=rev( consistentgenesUDwn$gene) )
pltthis1[,dir:="up"]
pltthis1[gene%in%consistentDownreg$gene,dir:="down"]

ggplot( pltthis1, aes(y=Estimate,x=gene) ) + theme_classic() + geom_point()+
  #jitter(width= 0, height= 0.1) + 
  geom_hline(aes(yintercept= 0),linetype=2) +
  labs(x= "Gene" , y= "Fold change in cancer cell expression") +facet_wrap(~Treatment)+geom_boxplot(aes(fill=dir,group=gene)) +coord_flip() +
  scale_fill_aaas(name="Change under ET", labels=c("Downregulation","Upregulation")) +theme(aspect.ratio=2)
ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Single cell most significant changes in cancer cell expression in resistant cells.png", height=8, width = 8)

pltthis1[Treatment=="letrozole + ribo",Treatment:="letrozole + \n ribociclib"]
ggplot( pltthis1, aes(y=Estimate,x=gene) ) + theme_classic(base_size=26) + geom_point()+
  #jitter(width= 0, height= 0.1) + 
  geom_hline(aes(yintercept= 0),linetype=2) +
  scale_fill_manual(name="Change under ET", labels=c("Downregulation","Upregulation"),values=c("yellow","red"))+
  labs(x= "Gene" , y= "Fold change in cancer cell expression") +facet_wrap(~Treatment)+geom_boxplot(aes(fill=dir,group=gene)) +coord_flip() +
  #scale_fill_aaas(name="Change under ET", labels=c("Downregulation","Upregulation")) +
  theme(aspect.ratio=2)
ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/Single cell most significant changes in cancer cell expression in resistant cells YR.pdf", height=18, width = 18)


discData <- data.table( pltthis1%>%dplyr::select(-c(dynamic_class) ))
setnames( discData , old=c("rn"), new=c("ModelParamLM"))
#write.csv(discData,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/DiscoveryCohortDEGPrePost.csv")

#
pltthisAdded <- allDE[rn!="(Intercept)"][  gene%in% c("EGFR","ERBB2","ERBB3","ERBB4","MYC"  )] #[Treatment=="letrozole"]
pltthisAdded$gene <- factor(pltthisAdded$gene , levels=c("EGFR","ERBB2","ERBB3","ERBB4","MYC"))
pltthisAdded[,dir:="up"]
discDataAdded <- data.table( pltthisAdded%>%dplyr::select(-c(dynamic_class) ))
setnames( discDataAdded , old=c("rn"), new=c("ModelParamLM"))
write.csv(discDataAdded,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/DiscoveryCohortDEGaddedERBBPrePost.csv")



ggplot( allDE[rn!="(Intercept)"][Treatment=="letrozole"][#Estimate!=0][
  gene%in% summary_allDE[nsignif> (length(unique(allDE[Treatment=="letrozole"]$Patient.Study.ID)))/2 ][order(-(mean_effect))] [1:topN] $gene]#[
    #gene%in%summary_allDE[order(-(mean_effect))] [1:20] $gene   ][order(gene)]
  , aes(Estimate,gene)) + geom_jitter(width=0,height=0.1)


summary_allDE[order(-(mean_effect))] [1:20]

summary_allDE[order(-abs(mean_effect))] [1:20]
summary_allDE[order(-abs(nsignif))] [1:20]
summary_allDE[mean_effect>0][order(-abs(nsignif))] [1:20]
summary_allDE[nsignif>4]
summary_allDE[order(-(mean_effect))] [1:40]


ggplot( start_enddd [] , aes(y=sqrt(GFRA1), x=Day ,group= Day)  ) + theme_classic() + 
  geom_boxplot(fill="red", alpha=0.2) +   geom_point() +facet_wrap(~Treatment)

ggplot( start_enddd [] , aes(y=sqrt(RERG), x=Day ,group= Day)  ) + theme_classic() + 
  geom_boxplot(fill="red", alpha=0.2) +   geom_point() +facet_wrap(~Treatment)

ggplot( start_enddd [] , aes(y=sqrt(BHLHE40), x=Day ,group= Day)  ) + theme_classic() + 
  geom_boxplot(fill="red", alpha=0.2) +   geom_point() +facet_wrap(~Treatment)





allUMAP <- rbindlist( lapply(1:nrow(lu2), function(i){
  cat("Tumor    ") ;  cat(i) ; cat("\n")
  dd_i <- start_enddd[ Patient.Study.ID==lu2[i]$Patient.Study.ID]   #dd_i %>% dim()   #ggplot(dd_i, aes(y=sqrt(JUN), x=Day , group=Day)  ) + geom_boxplot() + geom_point() #m1 <- lm(formula= paste0("sqrt(","JUN",") ~ as.factor(Day)") ,    dd_i )
  # Perform umap analysis
  set.seed(12345)
  umap_mod <- umap( sqrt(as.matrix(dd_i %>% dplyr::select( one_of(expressedGenes) ) ) ) ,n_components=2 , n_neighbors=10 )
  ures <- cbind(dd_i %>% dplyr::select(Cell.ID,Patient.Study.ID, Day,Treatment,ARM,dynamic_class,dynamic_class3), umap_mod$layout)
  return(ures)
}))
save(allUMAP,  file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Gene expression landscapes/cancer gene expression landscape/Tumor specific UMAPS strt vs end cancer gene expression scrna.RData")

ggplot(allUMAP[Treatment=="letrozole"][dynamic_class3!="Non-response"], aes(V1,V2,col=Day))+geom_point()+facet_wrap(dynamic_class3~Patient.Study.ID,scales="free")
ggplot(allUMAP[Treatment=="letrozole"], aes(V1,V2,col=Day))+geom_point()+facet_wrap(dynamic_class3~Patient.Study.ID,scales="free", ncol=4)
ggplot(allUMAP[Treatment!="letrozole"], aes(V1,V2,col=Day))+geom_point()+facet_wrap(dynamic_class3~Patient.Study.ID,scales="free", ncol=4)

ggplot(allUMAP[Treatment=="letrozole"][Patient.Study.ID=="001-103"], aes(V1,V2,col=Day))+geom_point()+facet_wrap(dynamic_class3~Patient.Study.ID,scales="free")



set.seed(12345)
umap_mod <- umap( sqrt(as.matrix(dd_i %>% dplyr::select( one_of(expressedGenes) ) ) ) ,n_components=2 ,n_neighbors=10 )
dd_i_umap <- cbind(dd_i, umap_mod$layout)
ures<- cbind(dd_i%>% dplyr::select(Cell.ID,Patient.Study.ID, Day,Treatment,ARM,dynamic_class,dynamic_class3), umap_mod$layout)

ggplot( dd_i_umap , aes( V1, V2 , col= as.factor(Day) ) ) + theme_classic(base_size=18) + theme(aspect.ratio=1) + geom_point(size= 1) + #facet_wrap(~TreatmentLb2) +
  labs(y="Umap B",x="Umap A")
ggsave(   file="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cancer GF analyses/UMAP of signif Tumor Mean changes in cancer cell expression in resistant cells.png", height=4, width = 4)

corumap <- data.table( t(cor( dd_i_umap%>%dplyr::select(V1,V2) , sqrt( dd_i_umap%>%dplyr::select(one_of(expressedGenes ))) )) , keep.rownames = T)
corumap[order(-abs(V1))][1:10]
corumap[order(-abs(V2))][1:30]
corumap[order(-(V2))][V2>0.5]
ggplot( summ_wide_umap , aes( x=V1, y= V2 , col= TRIB1 ) ) + theme_classic() + geom_point( size= 4) + labs(y="Umap B",x="Umap A") + scale_color_viridis_c()
ggplot( dd_i_umap , aes( x=V1, y= V2 , col= log(1+ZBTB16) ) ) + theme_classic() + geom_point( size= 4) + labs(y="Umap B",x="Umap A") + scale_color_viridis_c()
ggplot( dd_i_umap , aes( x=V1, y= V2 , col= log(1+B2M) ) ) + theme_classic() + geom_point( size= 4) + labs(y="Umap B",x="Umap A") + scale_color_viridis_c()
ggplot( dd_i_umap , aes( x=V1, y= V2 , col= sqrt(1+ESR1) ) ) + theme_classic() + geom_point( size= 4) + labs(y="Umap B",x="Umap A") + scale_color_viridis_c()


require(glmnet)
dd_i[, response:=0]
dd_i[Day==180, response:=1]

glmnet_binomB <- cv.glmnet(x= scale( sqrt(as.matrix(dd_i %>% dplyr::select( one_of(expressedGenes) ) ) )) , y= dd_i$response,family= "binomial", type.measure="class")
coef_est <- coef(glmnet_binomB,s="lambda.min")
coefnonzero<-data.table( t(t(coef_est [which(coef_est!=0),])) ,keep.rownames = T)
coefvals <- data.table( coefnonzero[order(-abs(V1))] )
#plot(glmnet_binomB)

coefvals[rn!="(Intercept)"][order(-abs(V1))][1:10]
pathanalysis0[rn!="(Intercept)"][Estimate!=0][adj.pval<0.05][order(-abs(Estimate))][1:10]

ggplot(dd_i, aes(y=sqrt(IL1RAPL1), x=Day , group=Day)  ) + geom_boxplot() + geom_point() #m1 <- lm(formula= paste0("sqrt(","JUN",") ~ as.factor(Day)") ,    dd_i )
ggplot(dd_i, aes(y=sqrt(MALAT1), x=Day , group=Day)  ) + geom_boxplot() + geom_point() #m1 <- lm(formula= paste0("sqrt(","JUN",") ~ as.factor(Day)") ,    dd_i )

ggplot(dd_i, aes(y=sqrt(ERRFI1), x=Day , group=Day)  ) + geom_boxplot() + geom_point() #m1 <- lm(formula= paste0("sqrt(","JUN",") ~ as.factor(Day)") ,    dd_i )

lassoIterative <- rbindlist( lapply(1:30,function(ii){
  glmnet_binomB <- cv.glmnet(x= as.matrix(dd_i %>% dplyr::select( one_of(expressedGenes) ) ), y= dd_i$response,family= "binomial",type.measure="deviance")
  #plot(glmnet_binomB)
  coef_est <- coef(glmnet_binomB,s="lambda.min")
  coefnonzero<-data.table( t(t(coef_est [which(coef_est!=0),])) ,keep.rownames = T)
  coefvals <- data.table(repli=ii, coefnonzero[order(-abs(V1))] )
  return(coefvals)
}) )

summarlassoIterative <- data.table( lassoIterative%>%group_by(rn )%>%summarise(ndetections = length(repli) , meancoef=sum(V1)/ length(unique(lassoIterative$repli))))
plot(ndetections ~     meancoef, summarlassoIterative)
summarlassoIterative[ndetections>0.5*max(ndetections)]



pathanalysis0[rn=="as.factor(Day)180"][adj.pval<0.05][order(-abs(Estimate))][1:10] #pathanalysis0[rn=="as.factor(Day)180:TreatmentLb2ribocombo"][adj.pval<0.05][order(-abs(Estimate))]
pathanalysis0[rn=="as.factor(Day)180"][ pval<0.01][order(gene)]
pathanalysis0[rn=="as.factor(Day)180"][ adj.pval<0.05][order(gene)]



