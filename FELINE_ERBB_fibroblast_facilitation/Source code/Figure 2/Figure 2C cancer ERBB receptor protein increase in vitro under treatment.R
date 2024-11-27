rm(list=ls())
require(data.table)
require(tidyr)
require(ggplot2)
require(dplyr)
require(lme4); require(lmerTest)
require(ggsci)

# Define data location
SourceDataLoc <- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript Feline ERBB facilitation/Molecular Systems Biology/Revision and submission folder/Source Data/"
Intermediateloc <- paste0(SourceDataLoc,"Figure 2/")

# laod data
dd<-data.table(read.csv(file= paste0(Intermediateloc,"SourceData_Figure_2C_Western ERBB quantification Fulv and Afatinib.csv")))

# Quantify baseline conditions relative to b-actin control
dd[,baseline_condition:=(fulv.nM==0& afatinib.uM==0)]
dd[,baseline_intensity:= sum(baseline_condition*Relative_intesity_vsbactin), by=c("Cell_type","Cell_line","Replicate","Protein")]

# calc log2FC
dd[,log2Relative_intesity :=log2(Relative_intesity_vsbactin)]
dd[,log2FC:=log2(Relative_intesity_vsbactin)-log2(baseline_intensity)]

#initial plot of data distribution
dd$ Relative_intesity_vsbactin%>%log2()%>%hist()
dd$ log2FC%>%hist()

# annotate treatments and lineages (cell lines that are resistant or sensitive to drug)
dd[, Lineage:=interaction(Cell_line, Cell_type)]
dd[,TreatmentLab:=paste(fulv.nM,"nM Fulv +",afatinib.uM,"uM Afat")]
dd$TreatmentLab <- factor(dd$TreatmentLab, levels=dd$TreatmentLab%>%unique())

# Across cell line nested random effects of Fulv + Afatinib treatment  (Nested random regression intercept, fulv slope and afatinib effect for each replicate)
# example
tmp <- dd[Protein=="EGFR"]
me2i <- lmer(log2Relative_intesity ~ sqrt(fulv.nM) + afatinib.uM +
               (1 + fulv.nM + afatinib.uM | Lineage / Replicate ), data= tmp)
summary(me2i)
tmp$preds<-predict(me2i)
#me2<-lmer(log2(Relative_intesity_vsbactin) ~ sqrt(fulv.nM) + afatinib.uM +#Lineage*Replicate+ 
#          (1+sqrt(fulv.nM)+afatinib.uM|Cell_line/Cell_type/Replicate), data=tmp)
#me2a<-lmer(log2(Relative_intesity_vsbactin) ~ sqrt(fulv.nM) + afatinib.uM +  (1|Lineage) +  
#             (1+fulv.nM+afatinib.uM|Cell_line/Cell_type/Replicate), data=tmp)
#me2b <- lmer(log2FC ~ -1+ sqrt(fulv.nM) + afatinib.uM +  (0+(sqrt(fulv.nM))+afatinib.uM|Cell_line/Cell_type/Replicate), data=tmp)
#summary(me2)
#summary(me2a)
#summary(me2b)
#tmp$preds<-predict(me2b)
#tmp$preds_y<-predict(me2)
#tmp$preds_ya<-predict(me2a)
# ggplot(tmp, aes(y=log2FC   , col= afatinib.uM, x= sqrt(fulv.nM) )) + geom_point()+
#   facet_grid(Lineage~Replicate) + theme_classic()+theme(aspect.ratio=1)+
#   labs(y="Intensity (log2(FC))", x="Fulvestrant (sqrt(nM))")+#geom_smooth()+
#   geom_point(aes(y=preds) ,shape=2) +
#   geom_line(aes(y=preds) ,shape=2) +
#   scale_colour_continuous(name="Afatinib")#,col="black")


ggplot(tmp, aes(y=log2Relative_intesity   , col= afatinib.uM, x= sqrt(fulv.nM) )) + geom_point()+
  facet_grid(Lineage~Replicate) + theme_classic()+theme(aspect.ratio=1)+
  labs(y="Intensity (log2(FC))", x="Fulvestrant (sqrt(nM))")+#geom_smooth()+
  #geom_point(aes(y=preds) ,shape=2) +
  geom_line(aes(y=preds) ,shape=2) +
  scale_colour_continuous(name="Afatinib")#,col="black")

ggplot(tmp, aes(y=log2Relative_intesity-log2(baseline_intensity)   , col= afatinib.uM, x= sqrt(fulv.nM) )) +
  geom_boxplot(aes(y=preds-log2(baseline_intensity), fill= afatinib.uM, group=TreatmentLab) ,outlier.color=NA,alpha=0.6) +
  geom_point(aes(group=afatinib.uM), position=position_dodge(width=0.5))+
  #facet_grid(Lineage~Replicate) +
  theme_classic()+theme(aspect.ratio=1)+
  labs(y="Intensity (log2(FC))", x="Fulvestrant (sqrt(nM))")+#geom_smooth()
  scale_fill_continuous(name="Afatinib")+
  scale_colour_continuous(name="Afatinib")#,col="black")


# Predict using model, show fit to data and extract and vizualize significance
# Statistics
out<- rbindlist( lapply( unique(dd$Protein), function(i){
  tmp <- dd[Protein==i]
  me2i <- lmer(log2Relative_intesity ~ sqrt(fulv.nM) + afatinib.uM +
                 (1 + fulv.nM + afatinib.uM | Lineage / Replicate ), data= tmp )
  stats<- data.table( coef(summary(me2i)) , keep.rownames = T)
  setnames(stats, old=c("Std. Error", "t value", "Pr(>|t|)"), new=c("Std.Error", "tvalue", "pvalue"))
  stats[,Protein:=i ]
}) )

out[rn=="sqrt(fulv.nM)"]
out[rn=="afatinib.uM"]


# Fits
outpred <- rbindlist( lapply( unique(dd$Protein), function(i){
  dat <- dd[Protein==i]
  me2i <- lmer(log2Relative_intesity ~ sqrt(fulv.nM) + afatinib.uM +
                 (1 + fulv.nM + afatinib.uM | Lineage / Replicate ), data= dat)
  dat$preds<-predict(me2i)
  return(dat)
}) )

#Plot Fits
outpred[,afatinibLab:="0.5uM"]
outpred[afatinib.uM==0,afatinibLab:="DMSO"]
outpred$afatinibLab <- factor(outpred$afatinibLab, levels=c("DMSO","0.5uM"))
ggplot(outpred, aes(y=log2Relative_intesity-log2(baseline_intensity)   , col= afatinibLab, x= sqrt(fulv.nM) )) +
  geom_boxplot(aes(y=preds-log2(baseline_intensity), fill= afatinibLab, group=TreatmentLab) ,outlier.color=NA,alpha=0.6) +
  geom_point(aes(group=afatinibLab), position=position_dodge(width=0.5))+
  facet_wrap(~Protein, scales="free", ncol=4) +
  theme_classic(base_size=16)+theme(aspect.ratio=1)+
  labs(y="Protein level \n Relative intensity (log2(FC vs control))", x="Fulvestrant (sqrt(nM))")+#geom_smooth()
  scale_fill_discrete(name="Afatinib")+
  scale_colour_discrete(name="Afatinib")#,col="black")
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript Feline ERBB facilitation/Supplementary Figures/Western Quantification boxplots.pdf", height=6,width=10)

#Plot example fits
outpred_eg <- rbindlist( lapply( "EGFR", function(i){
  dat <- dd[Cell_line=="T47D"][Cell_type=="ETR"][Protein==i]
  me2i <- lmer(log2Relative_intesity ~ sqrt(fulv.nM) + afatinib.uM +
                 (1 + fulv.nM + afatinib.uM |  Replicate ), data= dat)
  dat$preds<-predict(me2i)
  return(dat)
}) )

#outpred_eg<-outpred[Cell_line=="T47D"][Cell_type=="ETR"][Protein=="EGFR"]
outpred_eg[,afatinibLab:="0.5uM"]
outpred_eg[afatinib.uM==0,afatinibLab:="DMSO"]
outpred_eg$afatinibLab <- factor(outpred_eg$afatinibLab, levels=c("DMSO","0.5uM"))
outpred_eg[ ,ReplicateLab:= "C"]
outpred_eg[Replicate=="1005-B"  ,ReplicateLab:= "B"]
outpred_eg[Replicate=="original" ,ReplicateLab:= "A"]

ggplot(outpred_eg, aes(y=log2Relative_intesity-log2(baseline_intensity)   , col= afatinibLab, x= sqrt(fulv.nM) )) +
  geom_line(data=outpred_eg[afatinib.uM ==0], aes(y=preds-log2(baseline_intensity), fill= afatinibLab, group=ReplicateLab) )+
  geom_boxplot(data=outpred_eg[afatinib.uM !=0],aes(y=preds-log2(baseline_intensity), fill= afatinibLab, group=TreatmentLab) ,outlier.color=NA,alpha=0.6,outliers =F) +
  geom_point(aes(shape=ReplicateLab,group=afatinibLab), position=position_dodge(width=0.5), size=2)+
  facet_wrap(paste0(Cell_line," (",Cell_type,")")~Protein, scales="free", ncol=4) +
  theme_classic(base_size=16)+theme(aspect.ratio=1)+
  labs(y="Protein level \n Relative intensity (log2(FC vs control))", x="Fulvestrant (sqrt(nM))")+#geom_smooth()
  scale_shape_discrete(name="Replicate")+
  scale_fill_discrete(name="Afatinib")+
  scale_colour_discrete(name="Afatinib")#,col="black")
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript Feline ERBB facilitation/Supplementary Figures/Western Quantification Example.pdf", height=6,width=6)



# Plot effects and significance
ggplot( out[rn!="(Intercept)"] , aes(x=Estimate,y=-log10(pvalue), col=rn))+
  geom_point()+
  theme_classic(base_size=16)+theme(aspect.ratio=1)+
  geom_hline(yintercept = -log10(0.05), linetype="dashed", col="grey")+
  scale_color_jco(name="Treatment", labels=c("Afatinib", "Fulvestrant"))+
  geom_text(aes(label=Protein ),nudge_x =0.2,nudge_y =-0.07, show.legend=F)+
  labs(y="Significance \n (-log10(pvalue))", x="Treatment effect size \n (lmer estimate)")+
  lims(x=c(-6,2))#geom_smooth()
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript Feline ERBB facilitation/Supplementary Figures/Western Quantification volcanoplot.pdf", height=7,width=7)

