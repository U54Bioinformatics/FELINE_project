rm(list=ls()); require(data.table);require(dplyr);require(ggplot2);require(tidyr);
library(colorspace)
require(segmented)
require("ggsci")

dd <- data.table( read.csv(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Tcell_Experiments/Exp49/RawDataQC_PC/Exp49_Plate136-137_TcellParsed_V2.csv"))
dd[,RepID:=1:nrow(dd)]
dd[,Rep:=1:length(RepID), by=c("Plate","Tcell_Number_Plated","Ribociclib_uM","IL15_ngmL")]

# First we'll continue to work with florescence measured total intensity; filtering ATP conten
dat2 <- data.table(dd %>% dplyr::select(c("RepID","Rep","Tcell_Number_Plated":"IL15_ngmL","Area_4hr":"Area_163hr","ATP_4hr":"ATP_163hr")))


# Convert to Long format 
dat_long <- data.table(gather(dat2, MeasureTime, Value, 'Area_4hr':'ATP_163hr'))
dat_long[, c("Measurement", "Time") := tstrsplit(MeasureTime, "_", fixed=TRUE)]
dat_long[, Hour := as.numeric( gsub("[^0-9]","",Time) )]
dat_long[, MeasureTime := NULL]
dat_long[, Time := NULL]

dat_wide <- data.table( spread(dat_long,Measurement,Value) )
dat_wide[,Treatlab:="Control"]
dat_wide[Ribociclib_uM!=0,Treatlab:="Ribociclib"]

dat_wide[,IL15Level:="IL-15: 1 ngmL"]
dat_wide[as.numeric(as.character(IL15_ngmL))<1,IL15Level:="IL-15: 0 ngmL"]
dat_wide[as.numeric(as.character(IL15_ngmL))>1,IL15Level:="IL-15: 5 ngmL"]

ggplot(dat_wide[Hour==4], aes(Area,Tcell_Number_Plated, col=as.factor(Hour) )) + geom_point()#+facet_wrap(~Hour)

ggplot(dat_wide[!is.na(ATP)], aes(ATP,Area, col=as.factor(Hour) )) + geom_point()+ theme_classic()+
  scale_color_discrete(name="Hour")+theme(aspect.ratio = 1)+
  coord_trans(y="sqrt", x="sqrt")+geom_smooth(method="lm", formula=y~-1+x)#+facet_wrap(~Hour)

ggplot(dat_wide[Hour==4], aes(ATP, Tcell_Number_Plated )) + geom_point()+ theme_classic()+
  coord_trans(y="sqrt", x="sqrt")+geom_smooth(method="lm", formula=y~-1+x)#+facet_wrap(~Hour)

ggplot(dat_wide[Tcell_Number_Plated== 5000], aes(y=Area, x=Hour, col=as.factor(IL15_ngmL) , group=IL15_ngmL)) +geom_smooth(method="gam", formula=y~s(x,k=6), se=F)+
  geom_point(size=1.5,aes(group=RepID))+ theme_classic()+
  facet_grid(~Ribociclib_uM)

mu<- dat_wide[Tcell_Number_Plated== 5000][Hour==min(Hour)][is.finite(ATP)]$ATP%>%mean()
muArea<- dat_wide[Tcell_Number_Plated== 5000][Hour==min(Hour)][is.finite(Area)]$Area%>%mean()

ATPdat <- dat_wide[Tcell_Number_Plated== 5000][Hour==max(Hour)]
ATPdat[,lnFC:=log2(ATP/mu)]
ATPdat[,FC:=(ATP/mu)]

Areadat <- dat_wide[Tcell_Number_Plated== 5000][Hour==max(Hour)]
Areadat[,lnFCArea:=log2(Area/muArea)]
Areadat[,FCArea:=(Area/muArea)]
#save(ATPdat, file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/Figure6/TcellATPIL15.RData")
#save(Areadat, file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/ImmuneFelinePaperSourceData/Figure6/TcellAreaIL15.RData")
require(mgcv)
g1<- gam(lnFC~ s(sqrt(IL15_ngmL), k=4) +s(I( sqrt(IL15_ngmL)*(Treatlab=="Ribociclib")), k=4)  , data=ATPdat)
summary(g1)

g2<- gam(lnFCArea~ s(sqrt(IL15_ngmL), k=4) +s(I( sqrt(IL15_ngmL)*(Treatlab=="Ribociclib")), k=4)  , data=Areadat)
summary(g2)

summary(lm(log(ATP)~Treatlab,data=ATPdat[IL15_ngmL==0][!is.na(ATP)]))
summary(lm(log(Area)~Treatlab,data=ATPdat[IL15_ngmL==0][!is.na(Area)]))


summary(lm(log(ATP)~IL15_ngmL+Treatlab,data=ATPdat[][!is.na(ATP)]))
summary(lm(log(Area)~IL15_ngmL+Treatlab,data=ATPdat[][!is.na(Area)]))

ggplot(ATPdat, aes(y=lnFC, x=sqrt(IL15_ngmL), col=Treatlab , fill=Treatlab ,group=Ribociclib_uM)) +
  geom_smooth(method="gam", formula=y~s(x,k=5), se=T)+
  geom_point(size=2.5,aes(group=RepID))+ theme_classic(base_size=26)+theme(aspect.ratio = 1)+
  coord_trans(y=scales::trans_new("2tox", function(x) 2^x, "log2"),
              x=scales::trans_new("square", function(x) x^2, "sqrt"))+
  scale_x_continuous(name="IL-15 (ng/mL)",breaks=sqrt(c(0,1,5,10)),labels=c(0,1,5,10))+
  scale_y_continuous(name="ATP fold change \n (relative to baseline)",breaks=log2(c(4,8,16)),labels=c(4,8,16))+
  scale_color_jco(name="Treatment")+
  scale_fill_jco(name="Treatment")
#ggsave(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Tcell_Experiments/Exp49/Modelling_JG/T cell alone ATP fold change under ribo and IL15.png", width=8,height=8)

ggplot(dat_wide[Tcell_Number_Plated== 5000][Hour==max(Hour)], aes(y=log2(ATP/mu), x=sqrt(IL15_ngmL), col=Treatlab , fill=Treatlab ,group=Ribociclib_uM)) +
  geom_smooth(method="gam", formula=y~s(x,k=5), se=T)+
  geom_point(size=4.5,aes(group=RepID))+ theme_classic(base_size=26)+theme(aspect.ratio = 1)+
  coord_trans(y="exp",x=scales::trans_new("square", function(x) x^2, "sqrt"))+
  scale_x_continuous(name="IL-15 (ng/mL)",breaks=sqrt(c(0,1,5,10)),labels=c(0,1,5,10))+
  scale_y_continuous(name="ATP fold change \n (relative to baseline)",breaks=log2(c(4,8,16)),labels=c(4,8,16))+
  scale_color_jco(name="Treatment")+
  scale_fill_jco(name="Treatment")+
  theme(legend.position="none",
      axis.title=element_blank(),
      axis.text=element_blank())
#ggsave(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Tcell_Experiments/Exp49/Modelling_JG/BLANK T cell alone ATP fold change under ribo and IL15.png", width=8,height=8)


### nice data 
ggplot(dat_wide[Tcell_Number_Plated== 5000], aes(y=log(Area), x=Hour, col=sqrt(IL15_ngmL) ,fill=sqrt(IL15_ngmL) , group=IL15_ngmL)) +
  geom_smooth(method="gam", formula=y~s(x,k=6), se=T)+
  geom_point(size=1.5,aes(group=RepID))+ theme_classic(base_size = 26)+theme(aspect.ratio = 1)+
  facet_grid(~Treatlab)+labs(y="T cell spheroid size (log)")+
  scale_fill_viridis_c(option="B",end=0.9,name= "IL15 (ng/mL)",breaks=sqrt(c(0,1,5,10)),labels=c(0,1,5,10))+
  scale_color_viridis_c(option="B",end=0.9,name= "IL15 (ng/mL)",breaks=sqrt(c(0,1,5,10)),labels=c(0,1,5,10))
#ggsave(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Tcell_Experiments/Exp49/Modelling_JG/T cell alone area trajectories under ribo and IL15 B.png", width=10,height=6)

ggplot(dat_wide[Tcell_Number_Plated== 5000][Hour==163][is.finite(ATP)], aes(y=log10(ATP), x=sqrt(IL15_ngmL) , col=Treatlab,fill=Treatlab,group=Treatlab)) +
  geom_smooth(method="gam", formula=y~s(x,k=4), se=T)+
  geom_point(size=1.5,aes(group=RepID))+ 
  theme_classic(base_size = 26)+theme(aspect.ratio = 1)+
  scale_y_continuous(breaks=log10(c(20000,40000,80000,160000,320000,640000,1280000,2560000)),labels=c(20000,40000,80000,160000,320000,640000,1280000,2560000))+
  scale_color_jco(name="Treatent")+scale_fill_jco(name="Treatent")







ggplot(dat_wide[Tcell_Number_Plated== 5000][IL15_ngmL%in%c(0,1,5)], aes(y=log(Area), x=Hour, col=Treatlab , group=Treatlab)) +
  geom_smooth(method="gam", formula=y~s(x,k=6), se=F)+
  geom_point(size=1.5,aes(group=RepID))+ theme_classic(base_size=26)+
  facet_grid(~IL15Level)+
  scale_color_jco(name="Treatment")+
  scale_y_continuous(breaks=log(c(50000,500000,5000000)), labels=(c(50000,500000,5000000)))+
  labs(y="T cell spheroid area", x="Time (hours)")+ theme(aspect.ratio=1)


ggplot(dat_wide[Tcell_Number_Plated== 5000][IL15_ngmL%in%c(0,1,5)], aes(y=log(Area), shape=IL15Level,linetype=IL15Level,x=Hour, col=Treatlab , group=interaction(IL15Level,Treatlab))) +
  geom_smooth(method="gam", formula=y~s(x,k=6), se=F)+
  geom_point(size=1.5,aes(group=RepID))+ theme_classic(base_size=26)+
  scale_color_jco(name="Treatment")+
  scale_y_continuous(breaks=log(c(50000,500000,5000000)), labels=(c(50000,500000,5000000)))+
  labs(y="T cell spheroid area", x="Time (hours)")+ theme(aspect.ratio=1) +
  scale_linetype_manual(values=rev(c("solid", "dashed", "dotted")))


ggplot(dat_wide[Tcell_Number_Plated== 5000][IL15_ngmL%in%c(0,1,5)], aes(y=(Area), x=Hour, col=Treatlab , group=Treatlab)) +
  geom_smooth(method="gam", formula=y~s(x,k=6), se=F)+
  geom_point(size=1.5,aes(group=RepID))+ theme_classic()+
  facet_grid(~IL15Level)+
  scale_color_jco(name="Treatment")+
  scale_y_continuous(breaks=(c(100000,200000,400000,800000,1600000,3200000)))+
  labs(y="T cell spheroid area", x="Time (hours)") + theme(aspect.ratio=1)


# missing data for RepID 38
inputdata <- dat_wide[Tcell_Number_Plated== 5000][!is.na(Area)][RepID!=38]
inputdata[,ATP:=NULL]
inputdata[,Tcell_Number_Plated:=NULL]
inputdata[Ribociclib_uM==1][IL15_ngmL==0]
inputdata <- inputdata[order(RepID,Hour)]
inputdata[,Treatlab:="Control"]
inputdata[Ribociclib_uM!=0,Treatlab:="Ribociclib"]

#assume area propotional to count
ggplot(inputdata[Hour<max(Hour)][], aes(y=Area, x=Hour, col=as.factor(IL15_ngmL) , group=IL15_ngmL)) +
  #geom_smooth(method="gam", formula=y~s(x,k=6), se=F)+
  geom_path(size=1.5,aes(group=RepID))+ theme_classic()+
  facet_wrap(~Ribociclib_uM)

ggplot(inputdata[], aes(y=log(Area), x=Hour, col=sqrt(IL15_ngmL) , group=IL15_ngmL)) +
  geom_smooth(method="gam", formula=y~s(x,k=6), se=F)+
  geom_point(size=1.5,aes(group=RepID))+ theme_classic(base_size=26)+theme(aspect.ratio = 1)+
  facet_grid(~Treatlab)+theme(aspect.ratio = 1)+
  scale_color_viridis_c(name="IL-15 ngmL", option="B",end=0.9,
                        labels=c(0,1,5,10), breaks=sqrt(c(0,1,5,10)))+
  labs(y="T cell area", x="Time (hours)")+
  scale_y_continuous(labels=c(50000,500000,5000000), breaks=log(c(50000,500000,5000000)))+
  guides(colour=guide_colourbar(barheight=15))
ggsave(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Tcell_Experiments/Exp49/Modelling_JG/Images/T cell growth trajectories.png", width=10, height=10 )

mu2 <- dat_wide[Tcell_Number_Plated== 5000][Hour==min(Hour)][is.finite(Area)]$Area%>%mean()


ggplot(inputdata[Hour==max(Hour)], aes(y=log2(Area/mu2), x=sqrt(IL15_ngmL), col=Treatlab, fill=Treatlab )) +
  geom_boxplot(col="black", aes(group=interaction(IL15_ngmL,Treatlab) ))+
  stat_boxplot(geom="errorbar", col="black",aes(group=interaction(IL15_ngmL,Treatlab) ))+
  geom_point(size=1.5,aes(group=RepID))+ 
  theme_classic(base_size=26)+theme(aspect.ratio = 1)+
  theme(aspect.ratio = 1)+
  scale_color_jco(name="Treatmemt")+  scale_fill_jco(name="Treatmemt")+
  labs(y="T cell area \n (relative to baseline)", x="IL-15 ngmL")+
  scale_y_continuous(labels=c(2^(2:6)), breaks=(c(2:6)))+
  scale_x_continuous(breaks=sqrt(c(0,1,5,10)), labels=c(0, 1,5,10))+
  coord_trans(x=scales::trans_new("square", function(x) x^2, "sqrt")) 


ggplot(inputdata[Hour==max(Hour)], aes(y=log2(Area/mu2), x=sqrt(IL15_ngmL), col=Treatlab, fill=Treatlab )) +
  geom_smooth(method="gam", formula=y~s(x,k=6), se=T)+
  geom_point(size=3.5,aes(group=RepID))+ 
  theme_classic(base_size=26)+theme(aspect.ratio = 1)+
  theme(aspect.ratio = 1)+
  scale_color_jco(name="Treatmemt")+  scale_fill_jco(name="Treatmemt")+
  labs(y="T cell area \n (relative to baseline)", x="IL-15 ngmL")+
  scale_y_continuous(labels=c(2^(2:6)), breaks=(c(2:6)))+
  scale_x_continuous(breaks=sqrt(c(0,1,5,10)), labels=c(0, 1,5,10))+
  coord_trans(x=scales::trans_new("square", function(x) x^2, "sqrt")) 

ggsave(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Tcell_Experiments/Exp49/Modelling_JG/Final T cell area2.png", width=8, height=8 )


ggplot(inputdata[Hour==max(Hour)], aes(y=log2(Area/mu2), x=sqrt(IL15_ngmL), col=Treatlab, fill=Treatlab )) +
  geom_smooth(method="gam", formula=y~s(x,k=6), se=T)+
  geom_point(size=4.5,aes(group=RepID))+ 
  theme_classic(base_size=26)+theme(aspect.ratio = 1)+
  theme(aspect.ratio = 1)+
  scale_color_jco(name="Treatmemt")+  scale_fill_jco(name="Treatmemt")+
  labs(y="T cell area", x="IL-15 ngmL")+
  scale_y_continuous(labels=c(2^(2:6)), breaks=(c(2:6)))+
  scale_x_continuous(breaks=sqrt(c(0,1,5,10)), labels=c(0, 1,5,10))+
  coord_trans(x=scales::trans_new("square", function(x) x^2, "sqrt")) +
  theme(legend.position="none",
        axis.title=element_blank(),
        axis.text=element_blank())
ggsave(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Tcell_Experiments/Exp49/Modelling_JG/BLANK Final T cell area2.png", width=8,height=8)






