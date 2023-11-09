#i=26
rm(list=ls())
library(glmnet);require(dplyr);require(data.table);require(tidyr)
require(ggplot2);require(pheatmap);require(ggsci)
# Select one pathway to analyze
#ssgseapathways$pathway
#wch_signalTrans <-"HALLMARK_INFLAMMATORY_RESPONSE" # "HALLMARK_INFLAMMATORY_RESPONSE"# "HALLMARK_INTERFERON_GAMMA_RESPONSE"# "ST_JNK_MAPK_PATHWAY"  ## "HALLMARK_TGF_BETA_SIGNALING"  # "HALLMARK_INTERFERON_GAMMA_RESPONSE"     #"PID_BMP_PATHWAY" #"HALLMARK_TGF_BETA_SIGNALING" #"ST_JNK_MAPK_PATHWAY"
wch_signalTrans <-"HALLMARK_INTERFERON_GAMMA_RESPONSE"
#wch_signalTrans <-"KEGG_CELL_CYCLE"
#wch_signalTrans <-"REACTOME_CELL_CYCLE"
#wch_signalTrans <-"REACTOME_CELL_CYCLE_MITOTIC"
dir.create(paste0("/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/connect communication to intracellular transduction/Connect",wch_signalTrans,"/"))
dir.create(paste0("/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/connect communication to intracellular transduction/Connect",wch_signalTrans,"/output_SingleCelldd/"))

## Load general datasets
# Load subclone information
subconeInfo <- data.table( read.csv("/Users/jason/Dropbox/FELINE Project (1)/Manuscript/Sumbission folder/Resubmission version/Acceptance revision/Third publication revisions/Source data files/Bild_SourceData_Fig5.csv") %>%
                             dplyr::select( Cell.ID,	Patient.Study.ID,	ARM,	Day,	Response,	Subclone, CDK6_mu,mapk_axisV1,mapk_axisV2, ST_JNK_MAPK_PATHWAY.ssGSEA))


# Select just the intracellular signal transduction signatures of interest.
ssgseapathways <- data.table( read.csv("/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/SSGSEA list/Cancer SSGSEA pathway list_ab.csv"))[(!is.na(extra))|(pathway==wch_signalTrans),]#[!is.na(fib_ECM)| !is.na(immune)| !is.na(extra)]
# Load ssgsea data and umap landscape of cancer and normal epithelial cells
load(  file= "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Cancer gene expression wrangling/cancer and normal Communication on ssgsea landscape.RData")# udataComm

# names(udataComm)[grep("CELL",names(udataComm))]
# names(udataComm)[grep("CYCLE",names(udataComm))]
# names(udataComm)[grep("REACTOME_DOWNSTREAM_SIGNALING_OF_ACTIVATED_FGFR",names(udataComm))]

ssGSEAdat <- udataComm[Celltype=="Cancer cells"] %>%
  dplyr::select(Cell.ID,orig.ident,dynamic_class3,Patient.Study.ID,Day,orig.ident,key_,Celltype,one_of(c(wch_signalTrans,as.character(ssgseapathways$pathway)) ) )
PatIDunique <- unique(ssGSEAdat$Patient.Study.ID)
rm(list="udataComm")


# Load receptor list
load( "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/LigandReceptor/Filtered_Human-2015-Ramilowski-LR-pairs.RData" )
load(file ="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Communication output merged ALL ARMS/PopulationCommunicationMerged ALL ARMS.RData" )#CCI,allphenotypes, uu,perIndiv,
#unique( data.table(merge(CCI,LRpairsFiltered,by="Pair.Name"))$HPMR.Receptor)
LRgenelist <- unique( data.table(merge(CCI[ARM!="A"],LRpairsFiltered,by="Pair.Name"))$HPMR.Receptor) #unique( c(LRpairsFiltered$HPMR.Receptor) ) #, LRpairsFiltered$HPMR.Ligand) )
rm(list= c("LRpairs", "LRpairsFiltered"))


## Analysis of data subsets
# Select one patient
set.seed(1234)
results <- rbindlist(lapply(1:length(PatIDunique)
                            ,function(i){ #i=20
                              cat(i)
                              whichPatient <- PatIDunique[i]
                              PatIsubconeInfo <- subconeInfo[Patient.Study.ID==whichPatient]
                              PatIssGSEAdat   <- ssGSEAdat[Patient.Study.ID==whichPatient]   #dim(PatIssGSEAdat);dim(PatIsubconeInfo)
                              FELcode <- PatIssGSEAdat[1]$orig.ident
                              
                              # Get total signaling in the TME: For each sample day, extract the total ligand signaling from all cell types for each receptor in that tumor has
                              possDays <- sort(unique(PatIssGSEAdat$Day))
                              Full_LigandsTotal <- rbindlist(lapply(possDays,function(dd){ #dd<-14
                                saveloc <- "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/CommunicationOutputAllArms/"
                                savenm <- paste0("CommunicationResults__", "PatientID_", whichPatient, "__", "Day_", dd, ".RData")
                                load(file=paste0(saveloc, savenm))
                                LigandsTotal <- data.table(
                                  Communication[HPMR.Receptor %in% LRgenelist][ReceptorCelltype== "Cancer cells"] %>% 
                                    group_by(HPMR.Receptor,HPMR.Ligand) %>% dplyr::summarise( y= Transduction[1]/Receptor[1] ) %>% 
                                    group_by(HPMR.Receptor) %>% dplyr::summarise( LigandsTotal= sum(y) )
                                )
                                Full_LigandsTotal <- merge( data.table(HPMR.Receptor= LRgenelist), LigandsTotal, by="HPMR.Receptor", all= T)
                                Full_LigandsTotal[is.na(LigandsTotal), LigandsTotal:=0]
                                Full_LigandsTotal$Day <- dd
                                Full_LigandsTotal$Patient.Study.ID <- whichPatient
                                return(Full_LigandsTotal)
                              }))
                              
                              # Load RNA expression of receptor genes in cells from a patient's tumor samples 
                              CPMlocs <- "/Users/jason/Dropbox/FELINE Project/Data_analysis/FELINE_data_folder/scRNA_count_CPM/output/"
                              cpm_i <- data.table( fread( paste0(CPMlocs,"FEL011046_10x_", FELcode, "_gene_symbols.CPM.txt") ))[Gene.ID %in% LRgenelist]
                              LRGene.ID <- cpm_i$Gene.ID
                              LRtransposed <- as.data.table(  t(cpm_i[,-1])  , keep.rownames = T) #LRtransposed <- as.data.table( scale( t(cpm_i[,-1]) , center = T) , keep.rownames = T)
                              LRtransposed[ is.na(LRtransposed) ] <- 0
                              names(LRtransposed) <- c("Cell.ID", unname( unlist( cpm_i[,1] ) ) )   # hist(log(1+LRtransposed$TNFRSF18)) # LRtransposed[1:10,1:10]
                              
                              
                              # Subset RNA expression matrix by day and scale receptor activity by total ligand communication to that receptor at that time point in the TME
                              PatIcomm <- rbindlist( lapply(possDays, function(dd){
                                # reorder rows of total ligand communication data to match receptor expression matrix
                                tmp <- Full_LigandsTotal[HPMR.Receptor %in% colnames(LRtransposed)[-1] ][Day==dd]
                                setorder(tmp[,.r:=order(colnames(LRtransposed)[-1] )],.r  )[,.r:=NULL]
                                testvec <- tmp$LigandsTotal
                                # multiply recepor score by ligand total communcation
                                R_t <- LRtransposed [Cell.ID %in% PatIssGSEAdat[Day == dd]$Cell.ID ,]
                                res <- data.table(Cell.ID=R_t$Cell.ID, sweep((R_t [,-1]) , MARGIN=2, testvec, "*") )
                                return(res)
                              }) )# PatIcomm[1:10,1:10]
                              
                              ## Merge data ready for analysis
                              # Merge subclone info and gene expression data
                              PatIsubclonereceptor <- merge(PatIsubconeInfo[,!names(PatIsubconeInfo) %in% LRgenelist ,with=FALSE] , PatIcomm, by="Cell.ID")
                              PatIsubclonereceptorR <- merge(PatIsubconeInfo[,!names(PatIsubconeInfo) %in% LRgenelist ,with=FALSE] , LRtransposed, by="Cell.ID")
                              # Merge subclone+gene expression  info and ssgsea data
                              PatIssgseasubclonereceptor <- data.table(merge(PatIssGSEAdat, PatIsubclonereceptor, by=c("Cell.ID", "Patient.Study.ID", "Day")) )
                              PatIsubclonereceptorR <- data.table(merge(PatIssGSEAdat, PatIsubclonereceptorR, by=c("Cell.ID", "Patient.Study.ID", "Day")) )
                              save(PatIssgseasubclonereceptor,file=  paste0("/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/connect communication to intracellular transduction/Connect",wch_signalTrans,"/output_SingleCelldd/",whichPatient,"_",wch_signalTrans,"__AndCommunicationReceptor data.RData"))
                              save(PatIsubclonereceptorR,file=  paste0("/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/connect communication to intracellular transduction/Connect",wch_signalTrans,"/output_SingleCelldd/",whichPatient,"_",wch_signalTrans,"__AndCommunicationReceptor data plot.RData"))
                              
                              SubcloneList <- unique(PatIssgseasubclonereceptor$Subclone)
                              
                              # Select one subclone to analyze
                              receptorsDetected <- rbindlist(lapply(1:length(SubcloneList),function(j){ #j=1
                                cat(paste0("tumor  ", i,"       subclone",j))
                                
                                Subclonalfittingdata <- PatIssgseasubclonereceptor[Subclone == SubcloneList[j]] %>% dplyr::select(wch_signalTrans, Day,one_of(LRgenelist) )
                                #Subclonalfittingdata <- PatIsubclonereceptorR[Subclone == SubcloneList[j]] %>% dplyr::select(wch_signalTrans, Day,one_of(LRgenelist) )
                                #Subclonalfittingdata[1:10,1:10]
                                if(nrow(Subclonalfittingdata)>100){
                                  # Gather Lasso input data, perform 10 fold cross validation of predictions
                                  Y <- unname(unlist(Subclonalfittingdata[,1, with=F]))
                                  X0 <- sqrt( as.matrix(Subclonalfittingdata[,-c(1,2), with=F])  ) 
                                  X <- scale( X0[, colSums( X0 )>0] )
                                  #browser()
                                  mod_cv1 <- cv.glmnet(x= X, y= Y, family= "gaussian", nfolds= 10, intercept= T, alpha= 1)#, lower.limits=0)#
                                  countEffects <- sum(coef(mod_cv1, c(mod_cv1$lambda.min))[-1,1]>0)
                                  #if(countEffects==0){
                                  #  mod_cv <- cv.glmnet(x= X, y= Y, family= "gaussian", nfolds= 10, intercept= F, alpha= 1)#, lower.limits=0)#
                                  #}else{
                                    mod_cv <- mod_cv1
                                  #}
                                  #coef(mod_cv1, c(mod_cv1$lambda.min))
                                  #mod_cv <- cv.glmnet(x= X, y= Y, family= "gaussian", nfolds= 10, intercept= F, alpha= 1)#, lower.limits=0)#
                                  lambda_best <- mod_cv$lambda.min     #lambda_best <- mod_cv$lambda.1se 
                                  #plot(mod_cv)#
                                  #coef(mod_cv, c(mod_cv$lambda.min, mod_cv$lambda.1se))#lasso_model <- glmnet(x= X, y= Y, alpha= 1, lambda= lambda_best )#,lower.limits=0
                                  predictions_train <- predict(mod_cv, s= lambda_best, newx= X, type= "response")
                                  Subclonalfittingdata$prediction<- predictions_train
                                  Subclonalfittingdata$Y<- Y
                                  #mean(sqrt((predictions_train - Y )^2))  #summary( lm (predictions_train ~ Y   )  )$ adj.r.squared
                                  #ggplot(Subclonalfittingdata,aes(Y,prediction))+geom_point()+theme_classic() + geom_smooth(method="gam") + labs(y="Predicted signal transduction",x="Observed signal transduction (points=cells)") #plot(Y, predictions_train)
                                  receptor_predictors <- data.table( beta= as.matrix(coef(mod_cv,c(mod_cv$lambda.min))), keep.rownames = T)[abs(beta.1)!= 0][beta.rn!= "(Intercept)"]
                                  
                                  Subclonalfittingdata$Rsquared <- receptor_predictors$Rsquared <- summary( lm (predictions_train ~ Y   )  )$ adj.r.squared
                                  Subclonalfittingdata$Subclone <- SubcloneList[j]
                                  Subclonalfittingdata$Patient.Study.ID <- whichPatient
                                  Subclonalfittingdata$pathway <- wch_signalTrans
                                  Subclonalfittingdata$dynamic_class3<-unique(PatIssgseasubclonereceptor$dynamic_class3)
                                  Subclonalfittingdata$orig.ident<-unique(PatIssgseasubclonereceptor$orig.ident)
                                  Subclonalfittingdata$Celltype<-unique(PatIssgseasubclonereceptor$Celltype)
                                  Subclonalfittingdata$CDK6_mutation <- unique(PatIssgseasubclonereceptor$CDK6_mu)
                                  # lm fit on selected features
                                  Xdash <- as.matrix(Subclonalfittingdata[,receptor_predictors$beta.rn,with= F],drop=F)
                                  if(ncol(Xdash)>0){
                                    lm0 <- (summary(lm(Y ~  Xdash)) %>% coefficients)   
                                    receptor_predictors$lm.est <- coef(lm(Y ~  Xdash))[-1]    #  lm0[-1,1]
                                    tmp<- !is.na(receptor_predictors$lm.est)
                                    receptor_predictors$pval<-receptor_predictors$tval<-receptor_predictors$se <- NA
                                    receptor_predictors$se[which(tmp)] <-lm0[-1,2]
                                    receptor_predictors$tval[which(tmp)]  <- lm0[-1,3]
                                    receptor_predictors$pval[which(tmp)]  <- lm0[-1,4]
                                    receptor_predictors[pval<0.05][order(pval)]
                                    receptor_predictors[order(-beta.1)]
                                  }else{
                                    lm0 <- (summary(lm(Y ~  1)) %>% coefficients)
                                    receptor_predictors <- data.table(beta.rn=NA,beta.1=NA,Rsquared=NA,lm.est=NA,se=NA,tval=NA,pval=NA)
                                  }
                                  receptor_predictors$Subclone <- SubcloneList[j]
                                  receptor_predictors$Patient.Study.ID <- whichPatient
                                  receptor_predictors$pathway <- wch_signalTrans
                                  receptor_predictors$dynamic_class3 <- unique(PatIssgseasubclonereceptor$dynamic_class3)
                                  receptor_predictors$orig.ident <- unique(PatIssgseasubclonereceptor$orig.ident)
                                  receptor_predictors$Celltype <- unique(PatIssgseasubclonereceptor$Celltype)
                                  receptor_predictors$CDK6_mutation <- unique(PatIssgseasubclonereceptor$CDK6_mu)
                                  
                                  save(Subclonalfittingdata,receptor_predictors,mod_cv,mod_cv,lm0,file=paste0("/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/connect communication to intracellular transduction/Connect",wch_signalTrans,"/",
                                                                                                              wch_signalTrans,"in",whichPatient,SubcloneList[j],".RData"))
                                  return(receptor_predictors)
                                }
                              }))
                              
                              
                              return(receptorsDetected)
                            })
)

#save(results,file=paste0("/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/connect communication to intracellular transduction/Connect",wch_signalTrans,"/","OverallResultsReceptorsAssociatedwith",wch_signalTrans,".RData"))
#
wch_signalTrans <-"HALLMARK_INTERFERON_GAMMA_RESPONSE"
load(file=paste0("/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/connect communication to intracellular transduction/Connect",wch_signalTrans,"/","OverallResultsReceptorsAssociatedwith",wch_signalTrans,".RData"))
require(ggplot2)
nn <- nrow(unique(results[Patient.Study.ID%in%c(unique(allphenotypes[ARM!="A"]$Patient.Study.ID))]%>%dplyr::select(Patient.Study.ID,Subclone)))
#summaryOutput<-data.table( results[beta.1>0]%>%group_by(beta.rn,pathway)%>%summarise(N=n(), prop=n()/nn    ) )[N>2] [order(-N)]
#setnames(summaryOutput, old="beta.rn",new="receptor")
summaryOutput<-data.table( results[beta.1>0][lm.est>0][!grep("ITG",beta.rn)][ Patient.Study.ID%in%c(unique(allphenotypes[ARM!="A"]$Patient.Study.ID)) ]%>%group_by(beta.rn,pathway)%>%summarise(N=n(), prop=n()/nn    ) )[N>=2] [order(-N)]
setnames(summaryOutput, old="beta.rn",new="receptor")
summaryOutput$receptor <- factor(summaryOutput$receptor , levels=summaryOutput$receptor)
ggplot(summaryOutput[N>5][prop>0.1]    ,aes(x=receptor, y=prop))+geom_point(size=2.5) +theme_classic()+labs(y="Proportion subclones with \n IFN immune response associated to receptor \n (Hallmark IFN gamma response)",x="Cancer communication receptor")+theme(aspect.ratio=1, axis.text.x = element_text(angle = 90))
#ggsave("/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/connect communication to intracellular transduction/Finalised/Frequency of subclonal association of receptor and Hallmark IFN gamma.png",height=5,width=5)
# 9*9 plot : Frequency of subclonal association of receptor and WHITFIELD_CELL_CYCLE_G2_M #PID_BMP_PATHWAY

results[, Patient.Study.ID_Subclone:= paste0(orig.ident, Subclone)]
resMat <- results[Patient.Study.ID%in%c(unique(allphenotypes[ARM!="A"]$Patient.Study.ID))][beta.rn%in%summaryOutput[N>5] $receptor[1:200]]%>%dplyr::select(beta.rn,beta.1,Patient.Study.ID_Subclone)%>%spread(beta.rn,beta.1,fill=0)
resMat1 <- as.matrix(resMat[,-1])
resMat1[resMat1<0]<-0
rownames(resMat1) <- resMat[,1]$Patient.Study.ID_Subclone

annotation_row <-unique(results%>%dplyr::select(dynamic_class3, Patient.Study.ID_Subclone)) #orig.ident
rownames(annotation_row) <- annotation_row$Patient.Study.ID_Subclone;annotation_row$Patient.Study.ID_Subclone <-NULL
require(ggsci)
annotation_col = list(
  dynamic_class3 = rev(pal_npg("nrc")(2) ))
names(annotation_col[[1]] ) <- c("Response","Non-response")
pheatmap::pheatmap(sqrt(resMat1),annotation_row = annotation_row,annotation_colors=annotation_col)
#10*15 plot :WHITFIELD_CELL_CYCLE_G2_M continuous receptor drivers

pheatmap::pheatmap(1*(resMat1>0),annotation_row = annotation_row,annotation_colors=annotation_col)
#10*15 plot : WHITFIELD_CELL_CYCLE_G2_M binary receptor drivers





llfiles <- list.files(paste0("/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/connect communication to intracellular transduction/Connect",wch_signalTrans,"/output_SingleCelldd/"))
llfiles <- llfiles[grep("plot.RData",llfiles)]
results2 <- data.table( rbindlist(lapply(1:length(llfiles),function(j){ 
  load(file=        paste0("/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/connect communication to intracellular transduction/Connect",wch_signalTrans,"/output_SingleCelldd/", llfiles[j]) )
  PatIsubclonereceptorR[ ,Patient.Study.ID_Subclone:=paste0(orig.ident, Subclone) ]
  return(PatIsubclonereceptorR)  }) ) %>%
    dplyr::select(ARM,Patient.Study.ID,Day,dynamic_class3,Patient.Study.ID_Subclone,
                 c("IL18R1", summaryOutput[N>5]$receptor,one_of(wch_signalTrans))) )

results2[,Treatment:="Combination ribociclib"]
results2[ARM=="A",Treatment:="Letrozole alone"]

ggplot(results2[][], aes(x= sqrt(IL18R1), y= HALLMARK_INTERFERON_GAMMA_RESPONSE ) ) +  #, col=Patient.Study.ID_Subclone
  geom_point(alpha= 0.5) +
  geom_smooth(method="gam",formula=y ~ s(x,k=3),se=F) + theme_classic() +theme(legend.positio="none")#+facet_wrap(~Day)

ggplot(results2[IL18R1>0][IL15RA>0], aes(x= sqrt(IL18R1), y=  sqrt(IL15RA) ) ) +  geom_point(alpha= 0.5) +
  geom_smooth(method="gam",formula=y ~ s(x,k=3),se=F) + theme_classic() +theme(legend.positio="none") #+facet_wrap(~Day)

ggplot(results2[IL15RA>0][sqrt(IL15RA)<20], aes(x= (IL15RA), y= HALLMARK_INTERFERON_GAMMA_RESPONSE ) ) +  #, col=Patient.Study.ID_Subclone
  geom_point(alpha= 0.25, size= 1.5 ) +
  geom_smooth(method="gam",formula=y ~ s(x,k=3), se=T ) + theme_classic() + theme(legend.positio="none")#+facet_wrap(~Day)

ggplot(results2[IL15RA>0][sqrt(IL15RA)<25], aes(x= sqrt(IL15RA), y= HALLMARK_INTERFERON_GAMMA_RESPONSE ) ) +  #, col=Patient.Study.ID_Subclone
 geom_point(alpha= 0.25, size= 1.5,aes(col=Patient.Study.ID_Subclone) ) +
 geom_smooth(method="gam",formula=y ~ s(x,k=3), se=T ) + theme_classic() + theme(legend.positio="none")#+facet_wrap(~Day)

m1<-lm(HALLMARK_INTERFERON_GAMMA_RESPONSE~sqrt(IL15RA) , data= results2[IL15RA>0][] )
m1<-lm(sqrt(IL15RA)~ dynamic_class3, data= results2[IL15RA>0][Day==180] )
m1<-lm(HALLMARK_INTERFERON_GAMMA_RESPONSE~ dynamic_class3, data= results2[IL15RA>0][Day==180] )
summary(m1)
ggplot(results2[][], aes(x= sqrt(IL15RA), y= HALLMARK_INTERFERON_GAMMA_RESPONSE ) ) +  #, col=Patient.Study.ID_Subclone
  geom_point(alpha= 0.5) +
  geom_smooth(method="gam",formula=y ~ s(x,k=3),se=F) + theme_classic() +theme(legend.positio="none") +facet_wrap(Treatment~Day)


ggplot(results2[IL12RB1>0][], aes(x= sqrt(IL12RB1), y= HALLMARK_INTERFERON_GAMMA_RESPONSE ) ) +  #, col=Patient.Study.ID_Subclone
  geom_point(alpha= 0.5) +
  geom_smooth(method="gam",formula=y ~ s(x,k=3),se=F) + theme_classic() +theme(legend.positio="none")#+facet_wrap(~Day)

ggplot(results2[][], aes(y= sqrt(TLR2), x= dynamic_class3 ) ) +  #, col=Patient.Study.ID_Subclone
  geom_point(alpha= 0.5) + geom_violin(scale="width")
  geom_smooth(method="gam",formula=y ~ s(x,k=3),se=F) + theme_classic() +theme(legend.positio="none")#+facet_wrap(~Day)

#results2 <-data.table(results2%>%group_by(Patient.Study.ID_Subclone,Day) %>% mutate(ITGB8_scale=exp(scale(log(1+ITGB8))),
 #                                                                                   IL15RA_scale=exp(scale(log(1+IL15RA)))
#                                                                                    ))
names(results2) <- gsub(" ",".",names(results2))
results2[,discsqrtIL15RA_b:= 2.25*round(sqrt(IL15RA)/2.25 )]


#save(results2,file=paste0("/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/connect communication to intracellular transduction/Connect",wch_signalTrans,"/Singularcellres_ssgseavsCommunicationReceptordata.RData"))
results2[,responselabel:="Sensitive"]
results2[dynamic_class3=="Non-response",responselabel:="Resistant"]
#summarresults2 <- na.omit( data.table(results2[IL15RA>=0][] %>% group_by(Day,dynamic_class3,discsqrtIL15RA_b)%>%
#                               summarise(mu =mean(HALLMARK_INTERFERON_GAMMA_RESPONSE), 
#                                         ucl=mean(HALLMARK_INTERFERON_GAMMA_RESPONSE) + sd(HALLMARK_INTERFERON_GAMMA_RESPONSE)/sqrt(n()) ,
 #                                        lcl=mean(HALLMARK_INTERFERON_GAMMA_RESPONSE) - sd(HALLMARK_INTERFERON_GAMMA_RESPONSE)/sqrt(n()) ,
 #                                        ncell=n()  ) ) )[discsqrtIL15RA_b>0][ncell>30]
trm<-10
summarresults2 <- na.omit( data.table(results2[IL15RA>=0][] %>% group_by(Treatment,#Day,
                                                                         dynamic_class3,discsqrtIL15RA_b)%>%
                                        summarise(mu =mean(HALLMARK_INTERFERON_GAMMA_RESPONSE), 
                                                  ucl=mean(HALLMARK_INTERFERON_GAMMA_RESPONSE) + sd(HALLMARK_INTERFERON_GAMMA_RESPONSE)/sqrt(n()) ,
                                                  lcl=mean(HALLMARK_INTERFERON_GAMMA_RESPONSE) - sd(HALLMARK_INTERFERON_GAMMA_RESPONSE)/sqrt(n()) ,
                                                  ncell=n()  ) ) )[discsqrtIL15RA_b>0][ncell>trm]



summarresults2[,prop:=ncell/sum(ncell),by=c("Treatment","dynamic_class3")]
summarresults2[,responselabel:="Sensitive"]
summarresults2[dynamic_class3=="Non-response",responselabel:="Resistant"]


ggplot(summarresults2[],aes(x=discsqrtIL15RA_b, y=mu,col= dynamic_class3,fill=dynamic_class3)) + 
  geom_point(aes(size=100*prop),alpha=0.95)+ theme_classic(base_size=16)+#+theme( legend.position = "none")+
  geom_errorbar(aes(ymax=ucl,ymin=lcl)) +
  geom_smooth(data=results2[Treatment=="Combination ribociclib"][dynamic_class3=="Non-response"] [
                discsqrtIL15RA_b<=max( summarresults2[Treatment=="Combination ribociclib"][dynamic_class3=="Non-response"]$discsqrtIL15RA_b ) &
                discsqrtIL15RA_b>=min( summarresults2[Treatment=="Combination ribociclib"][dynamic_class3=="Non-response"]$discsqrtIL15RA_b ) ],
              aes(x=sqrt(IL15RA),y=HALLMARK_INTERFERON_GAMMA_RESPONSE ), method="gam",formula=y~s(x,k=4),alpha=0.3) +
  geom_smooth(data=results2[Treatment=="Combination ribociclib"][dynamic_class3=="Response"] [
    discsqrtIL15RA_b<=max( summarresults2[Treatment=="Combination ribociclib"][dynamic_class3=="Response"]$discsqrtIL15RA_b ) &
      discsqrtIL15RA_b>=min( summarresults2[Treatment=="Combination ribociclib"][dynamic_class3=="Response"]$discsqrtIL15RA_b ) ],
    aes(x=sqrt(IL15RA),y=HALLMARK_INTERFERON_GAMMA_RESPONSE ), method="gam",formula=y~s(x,k=4),alpha=0.3) +
  geom_smooth(data=results2[Treatment!="Combination ribociclib"][dynamic_class3=="Non-response"] [
    discsqrtIL15RA_b<=max( summarresults2[Treatment!="Combination ribociclib"][dynamic_class3=="Non-response"]$discsqrtIL15RA_b ) &
      discsqrtIL15RA_b>=min( summarresults2[Treatment!="Combination ribociclib"][dynamic_class3=="Non-response"]$discsqrtIL15RA_b ) ],
    aes(x=sqrt(IL15RA),y=HALLMARK_INTERFERON_GAMMA_RESPONSE ), method="gam",formula=y~s(x,k=4),alpha=0.3) +
  geom_smooth(data=results2[Treatment!="Combination ribociclib"][dynamic_class3=="Response"] [
    discsqrtIL15RA_b<=max( summarresults2[Treatment!="Combination ribociclib"][dynamic_class3=="Response"]$discsqrtIL15RA_b ) &
      discsqrtIL15RA_b>=min( summarresults2[Treatment!="Combination ribociclib"][dynamic_class3=="Response"]$discsqrtIL15RA_b ) ],
    aes(x=sqrt(IL15RA),y=HALLMARK_INTERFERON_GAMMA_RESPONSE ), method= "gam",formula= y~s(x,k=4), alpha=0.3) +
 labs(y="Cancer interferon response phenotype \n (ssGSEA Hallmark IFN gamma response)", x="IL15RA communication to cancer cells")+
 theme(aspect.ratio=1)+
 scale_color_npg(name="Tumor response", labels= c("Resistant","Sensitive"))+
 scale_fill_npg(name="Tumor response", labels= c("Resistant","Sensitive"))+
 scale_size_continuous(name="% cells") + facet_wrap(~Treatment)+
 scale_x_continuous(breaks=sqrt(c(50,100,200,400,800)) ,labels=c(50,100,200,400,800))
#ggsave("/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/connect communication to intracellular transduction/Finalised/Ribo and Letrozole Hallmark IFN gamma vs IL15RA communication.png",height=10,width=10)

summarresults2 <- na.omit( data.table(results2[IL15RA>=0][] %>% group_by(#Day,
                                                                         discsqrtIL15RA_b)%>%
                                        summarise(mu =mean(HALLMARK_INTERFERON_GAMMA_RESPONSE), 
                                                  ucl=mean(HALLMARK_INTERFERON_GAMMA_RESPONSE) + sd(HALLMARK_INTERFERON_GAMMA_RESPONSE)/sqrt(n()) ,
                                                  lcl=mean(HALLMARK_INTERFERON_GAMMA_RESPONSE) - sd(HALLMARK_INTERFERON_GAMMA_RESPONSE)/sqrt(n()) ,
                                                  ncell=n()  ) ) )[discsqrtIL15RA_b>0][ncell>trm]



summarresults2[,prop:=ncell/sum(ncell)]
results2[,discsqrtIL15RA_b:= 4.*round(sqrt(IL15RA)/4. )]
results2[,discsqrtIL15RA_b:= 3.8*round(sqrt(IL15RA)/3.8 )]


plotthis<-results2[discsqrtIL15RA_b>0]
#save( plotthis, file= "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Macrophage T cell communication/Ribo and Letrozole IL15 activation of IFN response Cohort1.RData" )
load( file= "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/Macrophage T cell communication/Ribo and Letrozole IL15 activation of IFN response Cohort1.RData" )


ggplot(results2[discsqrtIL15RA_b>0],aes(x=discsqrtIL15RA_b, y=HALLMARK_INTERFERON_GAMMA_RESPONSE,
                                        #col=as.factor(discsqrtIL15RA_b),
                                        )) + 
  geom_point(size=2,alpha=0.95)+ theme_classic(base_size = 22)+
  #+theme( legend.position = "none")+
  geom_boxplot(alpha=0.5,aes(fill=as.factor(discsqrtIL15RA_b),group=discsqrtIL15RA_b)) +
  scale_fill_brewer(palette = "Spectral", direction= 1)+
  #scale_color_brewer(palette = "Spectral", direction= 1)+
  theme(legend.position = "none")+  theme(aspect.ratio=1)+
  labs(y="Cancer interferon response phenotype \n (ssGSEA Hallmark IFN gamma response)", x="IL15RA communication to cancer cells")+
  scale_x_continuous(breaks=sqrt(c(50,100,200,400,800)) ,labels=c(50,100,200,400,800)) +
  geom_smooth(method="gam",formula=y~s(x,k=4),alpha=0.3)
  
paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Ribo and Letrozole Hallmark IFN gamma versus Tumor wide IL15.png"),height=10,width=10)



ggplot(results2[discsqrtIL15RA_b>0],aes(x=discsqrtIL15RA_b, y=HALLMARK_INTERFERON_GAMMA_RESPONSE,
                                        #col=as.factor(discsqrtIL15RA_b),
)) + 
  geom_point(size=2,alpha=0.95)+ theme_classic(base_size = 22)+
  #+theme( legend.position = "none")+
  geom_boxplot(alpha=0.5,aes(fill=as.factor(discsqrtIL15RA_b),group=discsqrtIL15RA_b)) +
  scale_fill_brewer(palette = "Spectral", direction= 1)+
  #scale_color_brewer(palette = "Spectral", direction= 1)+
  theme(legend.position = "none")+  theme(aspect.ratio=1)+
  labs(y="Cancer interferon response phenotype \n (ssGSEA Hallmark IFN gamma response)", x="IL15RA communication to cancer cells")+
  scale_x_continuous(breaks=sqrt(c(50,100,200,400,800)) ,labels=c(50,100,200,400,800)) +
  geom_smooth(method="gam",formula=y~s(x,k=4),alpha=0.3)+
  theme(axis.title=element_blank(),  axis.text=element_blank(),strip.text = element_blank())

paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"BLANK Ribo and Letrozole Hallmark IFN gamma versus Tumor wide IL15.png"),height=10,width=10)


ggplot(results2[discsqrtIL15RA_b>0],aes(x=discsqrtIL15RA_b, y=HALLMARK_INTERFERON_GAMMA_RESPONSE)) + 
  geom_point(aes(col=as.factor(discsqrtIL15RA_b)),size=2,alpha=0.95)+ theme_classic(base_size = 32)+
  #+theme( legend.position = "none")+
  geom_boxplot(alpha=0.5,col="black",outlier.colour =NA,aes(fill=as.factor(discsqrtIL15RA_b),group=discsqrtIL15RA_b)) +
  stat_boxplot(geom ='errorbar',aes(group=discsqrtIL15RA_b),position=position_dodge(width=1),col="black") + 
  # geom_boxplot(alpha = 0.6,col="black")+
  
  scale_fill_grey(start=0.8, end=0.2)+
  scale_color_grey(start=0.8, end=0.2)+
  #scale_fill_brewer(palette = "Greys", direction= 1)+
  #scale_color_brewer(palette = "Greys", direction= 1)+
  
  theme(legend.position = "none")+  theme(aspect.ratio=1)+
  labs(y="Cancer interferon response phenotype \n (ssGSEA Hallmark IFN gamma response)", x="IL15RA communication to cancer cells")+
  scale_x_continuous(breaks=sqrt(c(50,100,200,400,800)) ,labels=c(50,100,200,400,800)) +
  geom_smooth(method="gam",formula=y~s(x,k=4),alpha=0.3, size=2,col="darkgreen",fill="darkgreen")

paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Ribo and Letrozole Hallmark IFN gamma versus Tumor wide IL15 grey.png"),height=10,width=10)

ggplot(results2[discsqrtIL15RA_b>0],aes(x=discsqrtIL15RA_b, y=HALLMARK_INTERFERON_GAMMA_RESPONSE)) + 
  geom_point(aes(col=as.factor(discsqrtIL15RA_b)),size=2,alpha=0.95)+ theme_classic(base_size = 22)+
  #+theme( legend.position = "none")+
  geom_boxplot(alpha=0.5,col="black",outlier.colour =NA,aes(fill=as.factor(discsqrtIL15RA_b),group=discsqrtIL15RA_b)) +
  stat_boxplot(geom ='errorbar',aes(group=discsqrtIL15RA_b),position=position_dodge(width=1),col="black") + 
  # geom_boxplot(alpha = 0.6,col="black")+
  
  scale_fill_grey(start=0.8, end=0.2)+
  scale_color_grey(start=0.8, end=0.2)+
  #scale_fill_brewer(palette = "Greys", direction= 1)+
  #scale_color_brewer(palette = "Greys", direction= 1)+
  
  theme(legend.position = "none")+  theme(aspect.ratio=1)+
  labs(y="Cancer interferon response phenotype \n (ssGSEA Hallmark IFN gamma response)", x="IL15RA communication to cancer cells")+
  scale_x_continuous(breaks=sqrt(c(50,100,200,400,800)) ,labels=c(50,100,200,400,800)) +
  geom_smooth(method="gam",formula=y~s(x,k=4),alpha=0.3, size=2,col="darkgreen",fill="darkgreen") +
  theme(axis.title=element_blank(),  axis.text=element_blank(),strip.text = element_blank())
#ggsave(paste0(paperfile,"BLANK Ribo and Letrozole Hallmark IFN gamma versus Tumor wide IL15 grey.png"),height=10,width=10)

results2[,discsqrtIL15RA_b:= 0.2*round(sqrt(IL15RA)/0.2 )]
ggplot(results2[discsqrtIL15RA_b>0][discsqrtIL15RA_b<25],aes(x=discsqrtIL15RA_b, y=HALLMARK_INTERFERON_GAMMA_RESPONSE)) + 
  #geom_point(aes(col=as.factor(discsqrtIL15RA_b)),size=2,alpha=0.95)+
  theme_classic(base_size = 22)+ylim(-0.05,0.3)+
  geom_boxplot(alpha=0.5,col="black",outlier.colour =NA,aes(fill=as.factor(discsqrtIL15RA_b),group=discsqrtIL15RA_b)) +
  stat_boxplot(geom ='errorbar',aes(group=discsqrtIL15RA_b),position=position_dodge(width=1),col="black") + 
  scale_fill_grey(start=0.8, end=0.2)+
  scale_color_grey(start=0.8, end=0.2)+

  theme(legend.position = "none")+  theme(aspect.ratio=1)+
  labs(y="Cancer interferon response phenotype \n (ssGSEA Hallmark IFN gamma response)", x="IL15RA communication to cancer cells")+
  scale_x_continuous(breaks=sqrt(c(50,100,200,400,800)) ,labels=c(50,100,200,400,800)) +
  geom_smooth(method="gam",formula=y~s(x,k=4),alpha=0.3, size=2,col="darkgreen",fill="darkgreen") +
  theme(axis.title=element_blank(),  axis.text=element_blank(),strip.text = element_blank())
ggsave(paste0(paperfile,"BLANK Ribo and Letrozole Phenotype versus communication schematic.png"),height=10,width=10)







ggplot(results2[discsqrtIL15RA_b>0],aes(x=discsqrtIL15RA_b, y=HALLMARK_INTERFERON_GAMMA_RESPONSE,
                                        #col=as.factor(discsqrtIL15RA_b),
)) + 
  geom_point(size=2,alpha=0.95)+ theme_classic(base_size = 22)+
  #+theme( legend.position = "none")+
  geom_boxplot(alpha=0.5,aes(fill=as.factor(discsqrtIL15RA_b),group=discsqrtIL15RA_b)) +
  scale_fill_brewer(palette = "Spectral", direction= 1)+
  #scale_color_brewer(palette = "Spectral", direction= 1)+
  theme(legend.position = "none")+  theme(aspect.ratio=1)+
  labs(y="Cancer interferon response phenotype \n (ssGSEA Hallmark IFN gamma response)", x="IL15RA communication to cancer cells")+
  scale_x_continuous(breaks=sqrt(c(50,100,200,400,800)) ,labels=c(50,100,200,400,800)) +
  geom_smooth(method="gam",formula=y~s(x,k=4),alpha=0.3)+
  theme(axis.title=element_blank(),  axis.text=element_blank(),strip.text = element_blank())

paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"BLANK Ribo and Letrozole Hallmark IFN gamma versus Tumor wide IL15.png"),height=10,width=10)



#ggsave("/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/connect communication to intracellular transduction/Finalised/Ribo and Letrozole Merged result Hallmark IFN gamma vs IL15RA communication.png",height=10,width=10)

summary( lm(HALLMARK_INTERFERON_GAMMA_RESPONSE ~discsqrtIL15RA_b ,data= results2[discsqrtIL15RA_b>0]) )


ggplot(summarresults2[],aes(x=discsqrtIL15RA_b, y=mu)) + 
  geom_point(size=3,alpha=0.95)+ theme_classic(base_size=16)+#+theme( legend.position = "none")+
  geom_errorbar(aes(ymax=ucl,ymin=lcl)) +
  geom_smooth(data=results2[    discsqrtIL15RA_b<=max( summarresults2[]$discsqrtIL15RA_b ) &
                                  discsqrtIL15RA_b>=min(summarresults2[]$discsqrtIL15RA_b )  ],
    aes(x=sqrt(IL15RA),y=HALLMARK_INTERFERON_GAMMA_RESPONSE ), method="gam",formula=y~s(x,k=4),alpha=0.3) +
  labs(y="Cancer interferon response phenotype \n (ssGSEA Hallmark IFN gamma response)", x="IL15RA communication to cancer cells")+
  theme(aspect.ratio=1)+
  scale_color_npg(name="Tumor response", labels= c("Resistant","Sensitive"))+
  scale_fill_npg(name="Tumor response", labels= c("Resistant","Sensitive"))+
  scale_size_continuous(name="% cells") + #facet_wrap(~Treatment)+
  scale_x_continuous(breaks=sqrt(c(50,100,200,400,800)) ,labels=c(50,100,200,400,800))
#ggsave("/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/connect communication to intracellular transduction/Finalised/Ribo and Letrozole merged Hallmark IFN gamma vs IL15RA communication.png",height=10,width=10)


pltord<-data.table(results2[Day==180]%>%group_by(Patient.Study.ID)%>%summarise(m=median(HALLMARK_INTERFERON_GAMMA_RESPONSE)))[order(-m)]$Patient.Study.ID
results2$Patient.Study.ID<- factor(results2$Patient.Study.ID, levels=pltord )
ggplot( results2[Day==180][] ,aes(y=HALLMARK_INTERFERON_GAMMA_RESPONSE,x=Patient.Study.ID,fill=responselabel,col=responselabel))+
  geom_boxplot()+ theme_classic(base_size=20)+facet_wrap(~Treatment,scales="free_x",nrow=2)
 


pltord<-data.table(results2[Day==180]%>%group_by(Patient.Study.ID)%>%summarise(m=median(HALLMARK_INTERFERON_GAMMA_RESPONSE)))[order(-m)]$Patient.Study.ID
results2$Patient.Study.ID<- factor(results2$Patient.Study.ID, levels=pltord )
results2[,Daylab:=paste0("Day ",Day)]
ggplot( results2[Day==180][] ,aes(y=HALLMARK_INTERFERON_GAMMA_RESPONSE,x=responselabel,fill=responselabel,col=responselabel))+
  geom_jitter(size=0.5,width=0.01)+ geom_boxplot(col="black")+ theme_classic(base_size=20)+facet_wrap( Treatment~Daylab,ncol=1)+
  theme_classic(base_size=26)+
  theme(aspect.ratio=1)+
  scale_color_npg(name="Tumor response", labels=c("Resistant","Sensitive"))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant","Sensitive")) +
  labs(y="Cancer  interferon response phenotype \n (ssGSEA Hallmark IFN Gamma response)",
       x="Tumor response")+theme(legend.position="none")
#ggsave("/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/connect communication to intracellular transduction/Finalised/Ribo and Letrozole Hallmark IFN gamma distribution.png",height=10,width=10)

paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Ribo and Letrozole Hallmark IFN gamma distribution.png"),height=10,width=10)

ggplot( results2[Day==180][] ,aes(y=HALLMARK_INTERFERON_GAMMA_RESPONSE,x=responselabel,fill=responselabel,col=responselabel))+
  geom_jitter(size=0.95,width=0.01)+ geom_boxplot(col="black")+ theme_classic(base_size=20)+facet_wrap(~Treatment,scales="free_x", ncol=1)+
  theme_classic(base_size=20)+ theme(aspect.ratio=1)+
  scale_color_npg(name="Tumor response", labels=c("Resistant","Sensitive"))+
  scale_fill_npg(name="Tumor response", labels=c("Resistant","Sensitive")) +
  labs(y="Cancer  interferon response phenotype \n (ssGSEA Hallmark IFN Gamma response)",
       x="Tumor response")+theme(legend.position="none") +
  theme(axis.title=element_blank(),  axis.text=element_blank(),strip.text = element_blank())

paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"BLANK Ribo and Letrozole Hallmark IFN gamma distribution.png"),height=10,width=10)


summary(lm(HALLMARK_INTERFERON_GAMMA_RESPONSE~responselabel*Treatment ,data=results2[Day==180] ))
require(lmerTest)
summary(lmer(HALLMARK_INTERFERON_GAMMA_RESPONSE~responselabel + (1|Patient.Study.ID) ,data=results2[Day==180][Treatment=="Combination ribociclib"] ))
summary(lmer(HALLMARK_INTERFERON_GAMMA_RESPONSE~responselabel + (1|Patient.Study.ID) ,data=results2[Day==180][Treatment=="Combination ribociclib"] ))


load( "/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/LigandReceptor/Filtered_Human-2015-Ramilowski-LR-pairs.RData" )


results2 <- data.table( rbindlist(lapply(1:length(llfiles),function(j){ 
  load(file=        paste0("/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/connect communication to intracellular transduction/Connect",wch_signalTrans,"/output_SingleCelldd/", llfiles[j]) )
  PatIsubclonereceptorR[ ,Patient.Study.ID_Subclone:=paste0(orig.ident, Subclone) ]
  return(PatIsubclonereceptorR)  }) ) 
  ) 

results3<- results2 %>% dplyr::select(one_of(LRpairsFiltered$HPMR.Receptor))
dim(results3)
data.table(results2[IL15RA>0]$Patient.Study.ID_Subclone%>%table())[N>30]
data.table(results2[IL18R1>0]$Patient.Study.ID_Subclone%>%table())[N>30]
results2$Patient.Study.ID_Subclone%>%unique()%>%length()
corIL15<- data.table(t(  cor( sqrt(results3$IL15RA), sqrt(results3), method="spearman" )),keep.rownames = T)
corIL15[order(-abs(V1))][1:20]
corIL15[grep("IL18",rn)]
