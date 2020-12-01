mutat_dd2<-data.table(fread("/Users/jason/Downloads/FELINE_gene_cnv_barplot.patient.STable_allCNV.three_IDs.txt") %>% dplyr::select(Patient.Study.ID,response,variable,value))
setnames( mutat_dd2 ,old=c("variable","value"),new=c("gene_dna","mutationtype"))

mutat_dd2[,Success:=0]
mutat_dd2[response=="Responder",Success:=1]

mutat_dd2$mutationtype <- factor(mutat_dd2$mutationtype, levels = c("gain", "N","loss"))

mutat_dd2$mutationtype%>%unique
mutat_dd2$gene_dna%>%unique
mutat_dd2[, mutation_status:="Normal"]
mutat_dd2[gene_dna=="AKT1"&mutationtype=="gain", mutation_status:="Mutated"]
mutat_dd2[gene_dna=="AKT3"&mutationtype=="gain", mutation_status:="Mutated"]
mutat_dd2[gene_dna=="CCND1"&mutationtype=="gain", mutation_status:="Mutated"]
mutat_dd2[gene_dna=="CCNE1"&mutationtype=="gain", mutation_status:="Mutated"]
mutat_dd2[gene_dna=="CDK6"&mutationtype=="gain", mutation_status:="Mutated"]
mutat_dd2[gene_dna=="ERBB4"&mutationtype=="gain", mutation_status:="Mutated"]
mutat_dd2[gene_dna=="ESR1"&mutationtype=="loss", mutation_status:="Mutated"]
mutat_dd2[gene_dna=="FGFR1"&mutationtype=="loss", mutation_status:="Mutated"]

mutat_dd2[mutationtype!="N", mutation_status:="Mutated"]
mutat_dd2[gene_dna=="AKT1"&mutationtype=="loss", mutation_status:="Normal"]
mutat_dd2[gene_dna=="AKT1"&mutationtype=="loss", mutation_status:="Normal"]

ggplot(mutat_dd2,aes(y=Success, x=mutation_status))+facet_wrap(~gene_dna)+geom_point()
ggplot(mutat_dd2,aes(y=Success, x=mutationtype))+facet_wrap(~gene_dna)+geom_jitter(width=0, height=0.1,alpha=0.5)


logistRB1<-glm(Success~mutationtype,family="binomial", data= mutat_dd2[gene_dna=="RB1"])
summary(logistRB1)
logistTP53<-glm(Success~mutationtype,family="binomial", data= mutat_dd2[gene_dna=="TP53"])

summary(logistTP53)


logistAKT3<-glm(Success~mutationtype,family="binomial", data= mutat_dd2[gene_dna=="AKT3"])
summary(logistAKT3)
logistAKT1<-glm(Success~mutationtype,family="binomial", data= mutat_dd2[gene_dna=="AKT1"])
summary(logistAKT1)
logistESR1<-glm(Success~mutationtype,family="binomial", data= mutat_dd2[gene_dna=="ESR1"])
summary(logistESR1)

logistPIK3CA<-glm(Success~mutationtype,family="binomial", data= mutat_dd2[gene_dna=="PIK3CA"])
summary(logistPIK3CA)
