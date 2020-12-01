require(data.table);require(dplyr);require(ggplot2);require(tidyr);library(biomaRt)

CycleBaseB<-fread("/Users/jason/Downloads/human_periodic.tsv")
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_listB <- data.table(getBM(filters= "ensembl_peptide_id", attributes= c("ensembl_peptide_id","hgnc_symbol"),values=CycleBaseB$gene,mart= mart))
setnames(G_listB,old=c("ensembl_peptide_id","hgnc_symbol"),new=c("gene","Gene_code"))
CycleBase2B <- merge(unique(CycleBaseB%>%dplyr::select(gene)),G_listB,by="gene")
gene.listB <- CycleBase2B$Gene_code%>%unique
write.csv(gene.listB,file="Cyclebase31_Cell cycle genes_human_periodic.csv")

# CycleBase<-fread("/Users/jason/Downloads/human_experiments_phenotypes.tsv")
# G_list <- data.table(getBM(filters= "ensembl_peptide_id", attributes= c("ensembl_peptide_id","hgnc_symbol"),values=CycleBase$gene,mart= mart))
# setnames(G_list,old=c("ensembl_peptide_id","hgnc_symbol"),new=c("gene","Gene_code"))
# CycleBase2 <- merge(unique(CycleBase%>%dplyr::select(gene)),G_list,by="gene")
# gene.list <- CycleBase2$Gene_code%>%unique