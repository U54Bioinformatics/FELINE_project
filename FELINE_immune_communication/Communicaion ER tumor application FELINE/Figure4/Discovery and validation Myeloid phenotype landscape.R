rm(list=ls())
require(data.table)
require(dplyr)
require(ggplot2)
require(tidyr)
require(ggsci)
require(lmerTest)

load( file=paste0("~/Dropbox/Cancer_pheno_evo/data/FELINE2/PhenotypesAllArmsCohort2/CPMPhenotpyeLandscape_C2ProjectC2RevisednewMacrophagesDC.RData"))#C1umap,u_dat,DAY,cell_types_all,ARMS,Subtype,
Val_u_dat <- u_dat
rm(list=c("ARMS","C1umap","cell_types_all", "DAY","Subtype","u_dat"  ))
## Load saved data
load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/FELINE2/UMAP genes macrophage only ALLARMS/UMAP genes macrophage only ALLARMS.RData")
Disc_u_dat <- u_dat
rm(list=c("fulldd","gene_summary", "inflammationList","macrocelldd","u_dat","umap_in" ,"umap_in2" ,"umap_mod2", "umap_mod5"  ))

Disc_u_dat[,Cohort:="Discovery"]
Val_u_dat[,Cohort:="Validation"]

joint_dd <- intersect( names(Disc_u_dat),names(Val_u_dat) )[c(8695:8699,1:8694 )]

orderedrows<-c(joint_dd[1:15],sort(joint_dd[-c(1:15)]) )

myeloidUMAP <- rbind( data.table(Disc_u_dat %>%select( orderedrows )) , 
                          data.table(Val_u_dat %>%select( orderedrows )) )
myeloidUMAP[,dynamic_class:=NULL]
setnames(myeloidUMAP,old=paste0("V",1:4),new=paste0("UMAP",1:4))

#write.csv(myeloidUMAP,file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/Figure4/SourceData_Figure4_MyeloidGeneUMAP.csv")

sz <- 1.25
ggplot( Disc_u_dat[!(Celltype_subtype!="DC"&V1>2&V2< -1)],aes(-V1,V2,col=  Celltype_subtype)) +  
  geom_point(size=sz)+
  theme_classic()+
  theme(aspect.ratio=1,axis.text=element_blank(),axis.title=element_blank(),legend.position="none")

u_dat <- rbind(data.table(Cohort="Discovery", Disc_u_dat%>%dplyr::select(Cell.ID,Celltype,Patient.Study.ID,Day,ARM,Celltype_subtype,dynamic_class3,V1,V2,
                                                                         CD36,CYP27A1,DHRS9,LIPA,PPARG)),
               data.table(Cohort="Validation", Val_u_dat%>%dplyr::select(Cell.ID,Celltype,Patient.Study.ID,Day,ARM,Celltype_subtype,dynamic_class3,V1,V2,
                                                                         CD36,CYP27A1,DHRS9,LIPA,PPARG)))

u_dat[,TreatmentLab:="Combination ribociclib"]
u_dat[ARM=="A",TreatmentLab:="Letrozole alone"]

#save(u_dat, file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/temporary/Discovery and Validation Myeloid phenotype landscape.RData")
load(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/Figure4/Discovery and Validation Myeloid phenotype landscape.RData")
#u_dat
u_dat[Patient.Study.ID=="001-171_853686",Patient.Study.ID:="001-171"]
# visualize
sz <- 1.25

ggplot( u_dat[Cohort=="Validation"],aes(-V1,V2,col=  Celltype_subtype)) +  
  geom_point(size=sz)+#scale_color_ucscgb(name="Myleoid cell types")+
  geom_point(data=u_dat[Cohort=="Discovery"],size=sz)+
  theme_classic()+
  theme(aspect.ratio=1,axis.text=element_blank(),axis.title=element_blank(),legend.position="none")
  
paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
#ggsave(paste0(paperfile,"BLANK Discovery and Validation Ribo and Letrozole myeloid classification UMAP.png"),width=8,height=8, dpi=320)

cor( log(1+u_dat$CD36) ,u_dat$V2 )
cor( log(1+u_dat$CYP27A1) ,u_dat$V2 )
cor( log(1+u_dat$DHRS9) ,u_dat$V2 )
cor( log(1+u_dat$LIPA) ,u_dat$V2 )
cor( log(1+u_dat$PPARG) ,u_dat$V2 )

# CD36 : M2 marker of profibrotic macrophages that is upregulated by CSF1 stimulation and functions as a receptor of apoptotic cells promoting their removal and reducing inflammation and subsequent inflammatory signalling
ggplot( u_dat[] ,aes(-V1,V2,col=  log(1+CD36))) +  
  geom_point(size=sz)+scale_color_viridis_c(name="Myleoid cell types")+
  theme_classic()+
  theme(aspect.ratio=1,axis.text=element_blank(),axis.title=element_blank(),legend.text=element_blank(),legend.title=element_blank(),
        legend.key.size = unit(1.25, 'cm'))
#ggsave(paste0(paperfile,"BLANK Discovery and Validation Ribo and Letrozole myeloid M2 differentiation CD36 UMAP.png"),width=8,height=8, dpi=320)


# CYP27A1 : M2 marker of GF production by macrophages (27-Hydroxycholesterol (27-HC) is a cholesterol metabolite that promotes cancer growth): see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6834472/
ggplot( u_dat[] ,aes(-V1,V2,col=  log(1+CYP27A1))) +  
  geom_point(size=sz)+scale_color_viridis_c(name="Myleoid cell types")+
  theme_classic()+
  theme(aspect.ratio=1,axis.text=element_blank(),axis.title=element_blank(),legend.text=element_blank(),legend.title=element_blank(),
        legend.key.size = unit(1.25, 'cm'))
#ggsave(paste0(paperfile,"BLANK Discovery and Validation Ribo and Letrozole myeloid M2 differentiation CYP27A1 UMAP.png"),width=8,height=8, dpi=320)

# DHRS9 : M2b Mreg macrophage marker : see https://pubmed.ncbi.nlm.nih.gov/28594751/
ggplot( u_dat[] ,aes(-V1,V2,col=  log(1+DHRS9))) +  
  geom_point(size=sz)+scale_color_viridis_c(name="Myleoid cell types")+
  theme_classic()+
  theme(aspect.ratio=1,axis.text=element_blank(),axis.title=element_blank(),legend.text=element_blank(),legend.title=element_blank(),
        legend.key.size = unit(1.25, 'cm'))
#ggsave(paste0(paperfile,"BLANK Discovery and Validation Ribo and Letrozole myeloid M2 differentiation DHRS9 UMAP.png"),width=8,height=8, dpi=320)

# LIPA : M2 marker of fatty acid oxidation that distinctively supports M2 metabolism 
ggplot( u_dat[] ,aes(-V1,V2,col=  log(1+LIPA))) +  
  geom_point(size=sz)+scale_color_viridis_c(name="Myleoid cell types")+
  theme_classic()+
  theme(aspect.ratio=1,axis.text=element_blank(),axis.title=element_blank(),legend.text=element_blank(),legend.title=element_blank(),
        legend.key.size = unit(1.25, 'cm'))
#ggsave(paste0(paperfile,"BLANK Discovery and Validation Ribo and Letrozole myeloid M2 differentiation LIPA UMAP.png"),width=8,height=8, dpi=320)

# PPARG : M2 marker of regulated inflamation
ggplot( u_dat[ ] ,aes(-V1,V2,col=  log(1+PPARG))) +  
  geom_point(size=sz)+scale_color_viridis_c(name="Myleoid cell types")+
  theme_classic()+
  theme(aspect.ratio=1,axis.text=element_blank(),axis.title=element_blank(),legend.text=element_blank(),legend.title=element_blank(),
        legend.key.size = unit(1.25, 'cm'))
#ggsave(paste0(paperfile,"BLANK Discovery and Validation Ribo and Letrozole myeloid M2 differentiation PPARG UMAP.png"),width=8,height=8, dpi=320)

ggplot( u_dat[ ][],aes(-V1,V2,col=  dynamic_class3)) +  
  geom_point(size=sz)+ scale_color_npg(name="Tumor response", labels=c("Resistant","Sensitive"))+
  #labs(x="UMAP 1",y="UMAP 2")+
  labs(x="", y="")+
  theme_classic(base_size = 26)+
  theme(aspect.ratio=1,legend.position="none",axis.text=element_blank())+ facet_grid(Cohort~TreatmentLab)
#ggsave(paste0(paperfile,"BLANK Discovery and Validation Ribo and Letrozole TumorResponse myeloid UMAP.png"),width=8,height=8, dpi=320)

ggplot( u_dat[ ][],aes(-V1,V2,col=  dynamic_class3,shape=Cohort)) +  
  geom_point(size=sz)+ scale_color_npg(name="Tumor response", labels=c("Resistant","Sensitive"))+
  #labs(x="UMAP 1",y="UMAP 2")+
  labs(x="", y="")+
  theme_classic(base_size = 26)+
  theme(aspect.ratio=1,legend.position="none",axis.text=element_blank())+ facet_wrap(~TreatmentLab, ncol=1)


require(lmerTest)
M1 <- lmer(V2~dynamic_class3+Cohort+(1|Patient.Study.ID) , data=u_dat[TreatmentLab=="Combination ribociclib"])
M1%>%summary()

M1 <- lmer(V2~dynamic_class3+Cohort+(1|Patient.Study.ID) , data=u_dat[TreatmentLab=="Letrozole alone"])
M1%>%summary()

M1 <- lmer(V2~dynamic_class3+(1|Patient.Study.ID) , data=u_dat[Day!=180][ARM!="A"][!(Celltype_subtype!="DC"&V1>2&V2< -1)])
summary(M1)
M1b <- lmer(V2~dynamic_class3+(1|Patient.Study.ID) , data=u_dat[Day!=180][ARM=="A"][!(Celltype_subtype!="DC"&V1>2&V2< -1)])
summary(M1b)











# visualize
sz <- 1.25
ggplot( u_dat[Day!=180][ARM!="A"][!(Celltype_subtype!="DC"&V1>2&V2< -1)],aes(-V1,V2,col=  dynamic_class3)) +  
  geom_point(size=2.5)+scale_color_npg(name="Tumor response")+
  theme_classic()+
  theme(aspect.ratio=1,axis.text=element_blank(),axis.title=element_blank(),legend.position="none")
ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Macrophage phenotypes AllArms/BLANK Macrophage Ribociclib TumorResponse UMAP AllArms Cohort1 pts1_5.png")
paperfile<- "/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Paper figures Cancer immune communication/"
ggsave(paste0(paperfile,"BLANK Ribo and Letrozole RIBO Early treatment myeloid phenotype UMAP.png"))

ggplot( u_dat[Day==0][ARM=="A"][!(Celltype_subtype!="DC"&V1>2&V2< -1)],aes(-V1,V2,col=  dynamic_class3)) +  
  geom_point(size=2.5)+scale_color_npg(name="Tumor response")+
  theme_classic()+
  theme(aspect.ratio=1,axis.text=element_blank(),axis.title=element_blank(),legend.position="none")
ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Macrophage phenotypes AllArms/BLANK Macrophage Letrozole TumorResponse UMAP AllArms Cohort1 pts1_5.png")
ggsave(paste0(paperfile,"BLANK Ribo and Letrozole LETROZOLE Early treatment myeloid phenotype UMAP.png"))

ggplot( u_dat[!(Celltype_subtype!="DC"&V1>2&V2< -1)],aes(-V1,V2,col=  Celltype_subtype)) +  
  geom_point(size=sz)+#scale_color_ucscgb(name="Myleoid cell types")+
  theme_classic()+
  theme(aspect.ratio=1,axis.text=element_blank(),axis.title=element_blank(),legend.position="none")
ggsave(paste0(paperfile,"BLANK Ribo and Letrozole LETROZOLE myeloid classification UMAP.png"))

# CD36 : M2 marker of profibrotic macrophages that is upregulated by CSF1 stimulation and functions as a receptor of apoptotic cells promoting their removal and reducing inflammation and subsequent inflammatory signalling
ggplot( u_dat[!(Celltype_subtype!="DC"&V1>2&V2< -1)] ,aes(-V1,V2,col=  log(1+CD36))) +  
  geom_point(size=sz)+scale_color_viridis_c(name="Myleoid cell types")+
  theme_classic()+
  theme(aspect.ratio=1,axis.text=element_blank(),axis.title=element_blank(),legend.text=element_blank(),legend.title=element_blank(),
        legend.key.size = unit(1.25, 'cm'))
ggsave(paste0(paperfile,"BLANK Ribo and Letrozole LETROZOLE myeloid M2 differentiation UMAP.png"))


corgenes[order(-V2)][1:20]

# CD36 : M2 marker of profibrotic macrophages that is upregulated by CSF1 stimulation and functions as a receptor of apoptotic cells promoting their removal and reducing inflammation and subsequent inflammatory signalling
ggplot( u_dat[!(Celltype_subtype!="DC"&V1>2&V2< -1)] ,aes(-V1,V2,col=  log(1+CD36))) +  
  geom_point(size=sz)+scale_color_viridis_c(name="Myleoid cell types")+
  theme_classic()+
  theme(aspect.ratio=1,axis.text=element_blank(),axis.title=element_blank(),legend.text=element_blank(),legend.title=element_blank(),
        legend.key.size = unit(1.25, 'cm'))
ggsave(paste0(paperfile,"BLANK Ribo and Letrozole myeloid M2 differentiation CD36 UMAP.png"))


# CYP27A1 : M2 marker of GF production by macrophages (27-Hydroxycholesterol (27-HC) is a cholesterol metabolite that promotes cancer growth): see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6834472/
ggplot( u_dat[!(Celltype_subtype!="DC"&V1>2&V2< -1)] ,aes(-V1,V2,col=  log(1+CYP27A1))) +  
  geom_point(size=sz)+scale_color_viridis_c(name="Myleoid cell types")+
  theme_classic()+
  theme(aspect.ratio=1,axis.text=element_blank(),axis.title=element_blank(),legend.text=element_blank(),legend.title=element_blank(),
        legend.key.size = unit(1.25, 'cm'))
ggsave(paste0(paperfile,"BLANK Ribo and Letrozole myeloid M2 differentiation CYP27A1 UMAP.png"))

# DHRS9 : M2b Mreg macrophage marker : see https://pubmed.ncbi.nlm.nih.gov/28594751/
ggplot( u_dat[!(Celltype_subtype!="DC"&V1>2&V2< -1)] ,aes(-V1,V2,col=  log(1+DHRS9))) +  
  geom_point(size=sz)+scale_color_viridis_c(name="Myleoid cell types")+
  theme_classic()+
  theme(aspect.ratio=1,axis.text=element_blank(),axis.title=element_blank(),legend.text=element_blank(),legend.title=element_blank(),
        legend.key.size = unit(1.25, 'cm'))
ggsave(paste0(paperfile,"BLANK Ribo and Letrozole myeloid M2 differentiation DHRS9 UMAP.png"))

# LIPA : M2 marker of fatty acid oxidation that distinctively supports M2 metabolism 
ggplot( u_dat[!(Celltype_subtype!="DC"&V1>2&V2< -1)] ,aes(-V1,V2,col=  log(1+LIPA))) +  
  geom_point(size=sz)+scale_color_viridis_c(name="Myleoid cell types")+
  theme_classic()+
  theme(aspect.ratio=1,axis.text=element_blank(),axis.title=element_blank(),legend.text=element_blank(),legend.title=element_blank(),
        legend.key.size = unit(1.25, 'cm'))
ggsave(paste0(paperfile,"BLANK Ribo and Letrozole myeloid M2 differentiation LIPA UMAP.png"))

# PPARG : M2 marker of regulated inflamation
ggplot( u_dat[!(Celltype_subtype!="DC"&V1>2&V2< -1)] ,aes(-V1,V2,col=  log(1+PPARG))) +  
  geom_point(size=sz)+scale_color_viridis_c(name="Myleoid cell types")+
  theme_classic()+
  theme(aspect.ratio=1,axis.text=element_blank(),axis.title=element_blank(),legend.text=element_blank(),legend.title=element_blank(),
        legend.key.size = unit(1.25, 'cm'))
ggsave(paste0(paperfile,"BLANK Ribo and Letrozole myeloid M2 differentiation PPARG UMAP.png"))





# CD36 : M2 marker of profibrotic macrophages that is upregulated by CSF1 stimulation and functions as a receptor of apoptotic cells promoting their removal and reducing inflammation and subsequent inflammatory signalling
ggplot( u_dat[!(Celltype_subtype!="DC"&V1>2&V2< -1)] ,aes(-V1,V2,col=  log(1+CD36))) +  
  geom_point(size=sz)+scale_color_viridis_c(name="CD36 \n Expression", breaks=log(1+ c(0,200,20000)), labels=c(0,200,20000))+
  theme_classic(base_size=20)+labs(y="UMAP 2", x="UMAP 1")+
  theme(aspect.ratio=1,#axis.text=element_blank(),axis.title=element_blank(),legend.text=element_blank(),legend.title=element_blank(),
        legend.key.size = unit(1.25, 'cm'))
ggsave(paste0(paperfile,"Ribo and Letrozole myeloid M2 differentiation CD36 UMAP.png"))


# CYP27A1 : M2 marker of GF production by macrophages (27-Hydroxycholesterol (27-HC) is a cholesterol metabolite that promotes cancer growth): see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6834472/
ggplot( u_dat[!(Celltype_subtype!="DC"&V1>2&V2< -1)] ,aes(-V1,V2,col=  log(1+CYP27A1))) +  
  geom_point(size=sz)+scale_color_viridis_c(name="CYP27A1 \n Expression", breaks=log(1+ c(0,10,100,1000)), labels=c(0,10,100,1000))+
  theme_classic(base_size=20)+labs(y="UMAP 2", x="UMAP 1")+
  theme(aspect.ratio=1,#axis.text=element_blank(),axis.title=element_blank(),legend.text=element_blank(),legend.title=element_blank(),
        legend.key.size = unit(1.25, 'cm'))
ggsave(paste0(paperfile,"Ribo and Letrozole myeloid M2 differentiation CYP27A1 UMAP.png"))

# DHRS9 : M2b Mreg macrophage marker : see https://pubmed.ncbi.nlm.nih.gov/28594751/
ggplot( u_dat[!(Celltype_subtype!="DC"&V1>2&V2< -1)] ,aes(-V1,V2,col=  log(1+DHRS9))) +  
  geom_point(size=sz)+scale_color_viridis_c(name="DHRS9 \n Expression", breaks=log(1+ c(0,10,100,1000)), labels=c(0,10,100,1000))+
  theme_classic(base_size=20)+labs(y="UMAP 2", x="UMAP 1")+
  theme(aspect.ratio=1,#axis.text=element_blank(),axis.title=element_blank(),legend.text=element_blank(),legend.title=element_blank(),
        legend.key.size = unit(1.25, 'cm'))
ggsave(paste0(paperfile,"Ribo and Letrozole myeloid M2 differentiation DHRS9 UMAP.png"))

# LIPA : M2 marker of fatty acid oxidation that distinctively supports M2 metabolism 
ggplot( u_dat[!(Celltype_subtype!="DC"&V1>2&V2< -1)] ,aes(-V1,V2,col=  log(1+LIPA))) +  
  geom_point(size=sz)+scale_color_viridis_c(name="LIPA \n Expression", breaks=log(1+ c(0,10,100,1000)), labels=c(0,10,100,1000))+
  theme_classic(base_size=20)+labs(y="UMAP 2", x="UMAP 1")+
  theme(aspect.ratio=1,#axis.text=element_blank(),axis.title=element_blank(),legend.text=element_blank(),legend.title=element_blank(),
        legend.key.size = unit(1.25, 'cm'))
ggsave(paste0(paperfile,"Ribo and Letrozole myeloid M2 differentiation LIPA UMAP.png"))

# PPARG : M2 marker of regulated inflamation
ggplot( u_dat[!(Celltype_subtype!="DC"&V1>2&V2< -1)] ,aes(-V1,V2,col=  log(1+PPARG))) +  
  geom_point(size=sz)+scale_color_viridis_c(name="PPARG \n Expression", breaks=log(1+ c(0,10,100,1000)), labels=c(0,10,100,1000))+
  theme_classic(base_size=20)+labs(y="UMAP 2", x="UMAP 1")+
  theme(aspect.ratio=1,#axis.text=element_blank(),axis.title=element_blank(),legend.text=element_blank(),legend.title=element_blank(),
        legend.key.size = unit(1.25, 'cm'))
ggsave(paste0(paperfile,"Ribo and Letrozole myeloid M2 differentiation PPARG UMAP.png"))



ggplot( u_dat[ARM!="A"][!(Celltype_subtype!="DC"&V1>2&V2< -1)],aes(-V1,V2,col=  dynamic_class3)) +  
  geom_point(size=1.5)+scale_color_npg(name="Tumor response",labels=c("Resistant","Sensitive"))+labs(x="UMAP 1",y="UMAP 2")+
  theme_classic()+
  theme(aspect.ratio=1)
ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Macrophage phenotypes AllArms/Macrophage Ribociclib TumorResponse UMAP AllArms Cohort1 pts1_5.png")

ggplot( u_dat[ARM=="A"][!(Celltype_subtype!="DC"&V1>2&V2< -1)],aes(-V1,V2,col=  dynamic_class3)) +  
  geom_point(size=1.5)+scale_color_npg(name="Tumor response",labels=c("Resistant","Sensitive"))+labs(x="UMAP 1",y="UMAP 2")+
  theme_classic()+
  theme(aspect.ratio=1)
ggsave(filename="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/FELINE 2/Macrophage phenotypes AllArms/Macrophage Letrozole TumorResponse UMAP AllArms Cohort1 pts1_5.png")

