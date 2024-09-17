
rm(list=ls())   
##devtools::install_github("rikenbit/nnTensor") #install.packages("rTensor")
#install.packages("https://cran.r-project.org/src/contrib/Archive/rTensor/rTensor_1.4.tar.gz", repo=NULL, type="source")
require(rTensor); require(nnTensor); library(abind)
require(data.table); require(dplyr); require(ggplot2); require(tidyr)
require(igraph);require(ggsci)
require(lmerTest)

savloc<- "/Users/jason/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/SI data/FigureS13/"
Summar_MC<- data.table(read.csv(file=paste0(savloc,"SourceData_FigureS13_CancerCelltoMyeloidIndividualCellCellCommunicationsInput.csv")))
countPairdd <- data.table(Summar_MC %>% group_by(Treat,Pair.Name,Differentiation) %>% summarise(n=n()) %>% group_by(Treat,Pair.Name)%>%
                            summarise(n=min(n))
)
Summar_MC<-Summar_MC[Pair.Name %in% countPairdd[n>10]$Pair.Name]


TrendsResMfromCV2diffRibo <- rbindlist( lapply( unique( Summar_MC$Pair.Name )  , function(pp){
  cat(pp)
  out <- tryCatch({
    mm1 <- lmer( I(log(1 + TransductionMu) ) ~ Differentiation + (1|Disc_V1) + (1|Patient.Study.ID), 
                 data= Summar_MC[Pair.Name==pp][Treat=="CombinationRibo"] )
    if( length(coef(summary(mm1)))==0 ){ stop() }else{
      data.table(Treat="CombinationRibo",Pair.Name= pp, data.table( coef(summary(mm1)) , keep.rownames = T))
      
    }
  },
  error=function(x){
    data.table(Treat="CombinationRibo",'Pair.Name'=pp, data.table("rn"=NA, "Estimate"=NA ,'Std. Error'=NA, df =NA,'t value'=NA,'Pr(>|t|)'=NA) )
  })
  return(out)
}))

TrendsResMfromCV2diffLetro <- rbindlist( lapply( unique( Summar_MC$Pair.Name )  , function(pp){
  cat(pp)
  out <- tryCatch({
    mm1 <- lmer( I(log(1 + TransductionMu) ) ~ Differentiation + (1|Disc_V1) + (1|Patient.Study.ID), 
                 data= Summar_MC[Pair.Name==pp][Treat=="LetrozoleAlone"] )
    if( length(coef(summary(mm1)))==0 ){ stop() }else{
      data.table(Treat="LetrozoleAlone",Pair.Name= pp, data.table( coef(summary(mm1)) , keep.rownames = T))
    }
  },
  error=function(x){
    data.table(Treat="LetrozoleAlone",'Pair.Name'=pp, data.table("rn"=NA, "Estimate"=NA ,'Std. Error'=NA, df =NA,'t value'=NA,'Pr(>|t|)'=NA) )
  })
  return(out)
}))

setnames(TrendsResMfromCV2diffRibo,old=c("Std. Error","t value","Pr(>|t|)"),new=c("Std.Error","tval","pval"))
TrendsResMfromCV2diffRibo[,c("L","R"):= tstrsplit(Pair.Name,"_",fixed=T)]
TrendsResMfromCV2diffRibo$ adjpval<- p.adjust(TrendsResMfromCV2diffRibo$pval)
TrendsResMfromCV2diffRibo[Estimate>0][order(pval)][pval<0.05][rn=="DifferentiationM2"]
TrendsResMfromCV2diffRibo[Estimate<0][order(pval)][pval<0.05][rn=="DifferentiationM2"]
TrendsResMfromCV2diffRibo[][order(pval)][pval<0.05][rn=="DifferentiationM2"][1:20]

setnames(TrendsResMfromCV2diffLetro,old=c("Std. Error","t value","Pr(>|t|)"),new=c("Std.Error","tval","pval"))
TrendsResMfromCV2diffLetro[,c("L","R"):= tstrsplit(Pair.Name,"_",fixed=T)]
TrendsResMfromCV2diffLetro$ adjpval<- p.adjust(TrendsResMfromCV2diffLetro$pval)
TrendsResMfromCV2diffLetro[][order(pval)][adjpval<0.05][rn=="DifferentiationM2"]
TrendsResMfromCV2diffLetro[][order(pval)][pval<0.05][rn=="DifferentiationM2"][1:20]

nn <- 5#25

short.table <- rbind(data.table(ChangeInM2 = "up", TrendsResMfromCV2diffRibo[Estimate>0][order(pval)][pval<0.05][rn=="DifferentiationM2"][1:nn]),
      data.table(ChangeInM2 = "up",TrendsResMfromCV2diffLetro[Estimate>0][order(pval)][pval<0.05][rn=="DifferentiationM2"][1:nn]),
      data.table(ChangeInM2 = "down",TrendsResMfromCV2diffRibo[Estimate<0][order(pval)][pval<0.05][rn=="DifferentiationM2"][1:nn]),
      data.table(ChangeInM2 = "down",TrendsResMfromCV2diffLetro[Estimate<0][order(pval)][pval<0.05][rn=="DifferentiationM2"][1:nn]))[order(Estimate)]

listplotpairs <- unique( short.table$Pair.Name )
listplotpairs <- unique( c(
 as.character( TrendsResMfromCV2diffRibo[Estimate>0][order(pval)][pval<0.05][rn=="DifferentiationM2"][1:nn]$Pair.Name),
 as.character( TrendsResMfromCV2diffLetro[Estimate>0][order(pval)][pval<0.05][rn=="DifferentiationM2"][1:nn]$Pair.Name),
 as.character( TrendsResMfromCV2diffRibo[Estimate<0][order(pval)][pval<0.05][rn=="DifferentiationM2"][1:nn]$Pair.Name),
 as.character( TrendsResMfromCV2diffLetro[Estimate<0][order(pval)][pval<0.05][rn=="DifferentiationM2"][1:nn]$Pair.Name)))

ddplot <- Summar_MC[Pair.Name%in% listplotpairs  ][Pair.Name != "VTN_ITGB6"]
ddplot[,scalelntransduction := scale( log(1 + TransductionMu) ) , by= Pair.Name ]
ddplot$Pair.Name <- factor( ddplot$Pair.Name, levels = unique(listplotpairs))
ddplot[,Treatment:= "Combination ribociclib"]
ddplot[Treat=="LetrozoleAlone",Treatment:= "Letrozole alone"]

ggplot(ddplot, aes(x=Pair.Name, y= scalelntransduction ,col= Differentiation,fill= Differentiation  )) + 
  theme_classic(base_size=18)+
  geom_point(position=position_dodge(width=0.9))+
  geom_violin(aes(group=interaction(Differentiation,Pair.Name)),position=position_dodge(width=0.9),col=NA,scale="width",alpha=0.65)+
  labs(y="Cancer to Myeloid communication",x="Communication pathway"  )+
  scale_color_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
  scale_fill_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
  facet_wrap(~Treatment,scale="free_y",nrow=4)+theme(axis.text.x = element_text(angle=90),legend.position = "bottom")
ggsave(filename=paste0("/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cohort 2 results/Cancer communication with M2 vs M1 Macrophages.png"),width=8,height=8)
 
ggplot(ddplot[Pair.Name%in%c("ADAM10_AXL","CSF1_CSF1R","TGFB1_TGFBR1","ZP3_MERTK","C3_ITGAX","SPP1_CD44","FN1_SDC2")], 
       aes(x=as.character(Pair.Name), y= scalelntransduction ,col= Differentiation,fill= Differentiation  )) + 
  theme_classic(base_size=18)+
  geom_point(position=position_dodge(width=0.9))+
  geom_violin(aes(group=interaction(Differentiation,Pair.Name)),position=position_dodge(width=0.9),col=NA,scale="width",alpha=0.65)+
  labs(y="Cancer to Myeloid communication",x="Communication pathway"  )+
  scale_color_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
  scale_fill_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
  facet_wrap(~Treatment,scale="free_y",nrow=4)+
  theme(aspect.ratio=0.5,axis.text.x = element_text(angle=90),legend.position = "bottom")
ggsave(filename=paste0("/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cohort 2 results/Trimmed Cancer communication with M2 vs M1 Macrophages.png"),width=8,height=8)

ggplot(ddplot, aes(x=Pair.Name, y= scalelntransduction ,col= Differentiation,fill= Differentiation  )) + 
  theme_classic(base_size=16)+
  geom_point(position=position_dodge(width=0.8))+
  geom_violin(aes(group=interaction(Differentiation,Pair.Name)),position=position_dodge(width=0.8),color=NA,scale="width",alpha=0.7)+
  labs(y="Cancer to Myeloid communications",x="Communication pathway"  )+
  scale_color_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
  scale_fill_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
  theme(aspect.ratio=0.5,axis.text.x = element_text(angle=90),legend.position = "bottom")
ggsave(filename=paste0("/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cohort 2 results/Cancer communication with M2 vs M1 Macrophages combined ribo and letrozole alone.png"),width=8,height=8)

savloc<-"/Users/jason/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/SI data/FigureS13/"
write.csv(ddplot,file=paste0(savloc,"SourceData_FigureS13_Myeloid subpopulationsTargettedbyCancerSignals.csv"))

datatoorder <- data.table( ddplot%>%group_by(Pair.Name,Treat,Differentiation)%>%summarise(mu=median(scalelntransduction))%>%spread(Differentiation,mu))
datatoorder[,diffM2_1:=M2 - M1]
datatoorder_wide <- datatoorder%>%select(Pair.Name,Treat,diffM2_1)%>%spread(Treat,diffM2_1)
ggplot( datatoorder_wide , aes(y= CombinationRibo, x=LetrozoleAlone)) + geom_point(size=2.5) + 
  geom_smooth(method="lm") + theme_classic(base_size=18) + theme(aspect.ratio = 1) +
  theme(aspect.ratio = 1) +
  labs(y="Cancer to Myeloid communication \n in ribociclib treatment arm", 
       x="Cancer to Myeloid communication \n in letrozole treatment arm" )
ggsave(filename= paste0("/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Cohort 2 results/Cancer communication with M2 vs M1 Macrophages comparing combined ribo and letrozole alone.png"), width= 6 , height= 6)




unique(datatoorder[order(-diffM2_1)]$Pair.Name)

ddplot$Pair.Name <- factor(ddplot$Pair.Name, levels= unique(datatoorder[Treat=="CombinationRibo"][order(-diffM2_1)]$Pair.Name) )
ddplot[,Treatment:="Combination ribociclib"]
ddplot[Treat=="LetrozoleAlone",Treatment:="Letrozole alone"]

ggplot(ddplot, aes(x=Pair.Name, y= scalelntransduction ,col= Differentiation,fill= Differentiation  )) + 
  theme_classic()+
  geom_point(position=position_dodge(width=0.8))+
  geom_violin(aes(group=interaction(Differentiation,Pair.Name)),col=NA,position=position_dodge(width=0.8),scale="width",alpha=0.6)+
  labs(y="Cancer to Myeloid communications",x="Communication pathway"  )+
  scale_color_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
  scale_fill_viridis(name="Myeloid phenotype",labels=c("M1","M2"), discrete=TRUE,option="B",begin = 0.1, end = 0.8) +
  facet_wrap(~Treatment,scale="free_y",nrow=4) +theme(axis.text.x = element_text(angle=90),legend.position = "bottom")



