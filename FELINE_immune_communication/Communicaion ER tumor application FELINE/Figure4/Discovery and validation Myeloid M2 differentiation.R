load(  file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/Figure4/Discovery and Validation myeloid differentiation.RData")
summarizedCohorts

ggplot(summarizedCohorts,  aes(y=TumorDifferScore, x=dynamic_class3,col=dynamic_class3, fill=dynamic_class3, shape=Cohort, group=interaction(dynamic_class3) )) +
  geom_boxplot(outlier.colour=NA, position= position_dodge() ,col="black")+
  stat_boxplot(geom="errorbar",position=position_dodge(1.75),width=0.5)+
  geom_point(col="black",size=2.5)+
  scale_color_npg(name="Tumor response", labels=c("Resistant", "Sensitive"),guide="none")+
  scale_fill_npg(name="Tumor response", labels=c("Resistant", "Sensitive"),guide="none")+
  scale_x_discrete(name="Tumor response", labels=c("Resistant", "Sensitive"))+
  theme_classic(base_size=28)+theme(aspect.ratio = 1)+facet_wrap(~Treatmentlab,ncol=1)+
  labs(y="Myeloid M2 differentiation \n (tumor mean)")+theme(legend.position="none")
paperfile<- "/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Discovery and Validation Ribo and Letrozole Boxplot Myeloid differentiation predicts response to ribo not letrozole.png"),height=10,width=10, dpi=320)

summary(lm(TumorDifferScore~dynamic_class3+Cohort,summarizedCohorts[Treatmentlab=="Combination ribociclib"]))
summary(lm(TumorDifferScore~dynamic_class3+Cohort,summarizedCohorts[Treatmentlab=="Letrozole alone"]))
