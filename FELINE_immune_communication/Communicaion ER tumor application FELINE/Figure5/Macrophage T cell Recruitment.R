rm(list=ls())


load(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/temporary/DiscoveryMyeloidTcellRecruitment.RData")
Recruitmentdiscdat <- data.table(Cohort="Discovery", discdat )

load(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/temporary/DiscoveryMyeloidTcellActivation.RData")
Activationdiscdat <- data.table(Cohort="Discovery", discdat )


load(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/temporary/ValidationMyeloidTcellRecruitment.RData")
Recruitmentvaliddat <- data.table(Cohort="Validation", validdat )

load(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/temporary/ValidationMyeloidTcellActivation.RData")
Activationvaliddat <- data.table(Cohort="Validation", validdat )

Recruitmentalldat <- rbind(
  Recruitmentdiscdat%>%dplyr::mutate(Communication=y)%>%dplyr::select(Function,Cohort,Patient.Study.ID, Day, dynamic_class3, ARM,    Pair.Name,Pair.Name2,key_,LigandPhenoCelltype,ReceptorPhenoCelltype,Communication,Treatmentlab,Tumorresponse),
  Recruitmentvaliddat%>%dplyr::mutate(Communication=deltay)%>%dplyr::select(Function,Cohort,Patient.Study.ID, Day, dynamic_class3, ARM,    Pair.Name,key_,Pair.Name2,LigandPhenoCelltype,ReceptorPhenoCelltype,Communication,Treatmentlab,Tumorresponse)
)[!is.na(Communication)]
Recruitmentalldat[, Daylab:=paste0("Day ",Day)]
save(Recruitmentalldat, file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/Figure5/DiscoveryandValidationMyeloidTcellRecruitment.RData")


ggplot( Recruitmentalldat , aes(y= Communication, Tumorresponse,group=interaction(Tumorresponse,Pair.Name2),col= Pair.Name2,fill= Pair.Name2)) +
  theme_classic(base_size=26) +
  stat_boxplot(geom ='errorbar',position=position_dodge(width=1),col="black") + 
  geom_boxplot(alpha = 0.6,col="black",position=position_dodge(width=1))+
  theme(aspect.ratio=1)+
  geom_point(position=position_dodge(width=1))+ facet_grid(Treatmentlab~Cohort)+
  geom_vline(xintercept = 1.5,linetype="dashed") +
  scale_colour_brewer(name="Communication \n pathway", palette = "YlGn") +
  scale_fill_brewer(name="Communication \n pathway", palette = "YlGn") +
  labs(y="Myeloid to T cell \n recruitment communications post treatment",x="Tumor response")
paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"Discovery and Validation Ribo and letrozole Myeloid to T cell post treatment Recruitment communications.png"),height=11,width=11.75, dpi=320)

ggplot( Recruitmentalldat , aes(y= Communication, x=Tumorresponse, group=interaction(Tumorresponse,Pair.Name2), col= Pair.Name2, fill= Pair.Name2)) +
  theme_classic(base_size=28) +
  stat_boxplot(geom ='errorbar',position=position_dodge(width=1),col="black") + 
  geom_boxplot(alpha = 0.6,col="black",position=position_dodge(width=1))+
  theme(aspect.ratio=1)+
  geom_point(position=position_dodge(width=1))+ facet_wrap(Treatmentlab~Cohort, ncol=2)+
  geom_vline(xintercept = 1.5,linetype="dashed") +
  scale_colour_brewer(name="Communication \n pathway", palette = "YlGn") +
  scale_fill_brewer(name="Communication \n pathway", palette = "YlGn") +
  labs(y="Myeloid to T cell \n recruitment communication post treatment",x="Tumor response") +
  theme(axis.title=element_blank(),  axis.text=element_blank(),strip.text = element_blank(),legend.text=element_blank(),legend.title=element_blank())

paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"BLANK Discovery and Validation Ribo and letrozole Myeloid to T cell post treatment Recruitment communications.png"),height=10,width=10, dpi=320)

Recruitmentalldat[,Cohort2:=" "]
Recruitmentalldat[Cohort=="Validation",Cohort2:="  "]
ggplot( Recruitmentalldat , aes(y= Communication, Tumorresponse,group=interaction(Tumorresponse,Pair.Name2),col= Pair.Name2,fill= Pair.Name2)) +
  theme_classic(base_size=26) +
  stat_boxplot(geom ='errorbar',position=position_dodge(width=1),col="black") + 
  geom_boxplot(alpha = 0.6,col="black",position=position_dodge(width=1))+
  theme(aspect.ratio=1)+
  geom_point(position=position_dodge(width=1))+ facet_wrap(~paste0(Treatmentlab,Cohort2))+
  geom_vline(xintercept = 1.5,linetype="dashed") +
  scale_colour_brewer(name="Communication \n pathway", palette = "YlGn") +
  scale_fill_brewer(name="Communication \n pathway", palette = "YlGn") +
  labs(y="Myeloid to T cell \n recruitment communications post treatment",x="Tumor response")
paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
ggsave(paste0(paperfile,"Discovery and Validation Ribo and letrozole Myeloid to T cell post treatment Recruitment communicationsB.png"),height=12,width=12, dpi=320)
