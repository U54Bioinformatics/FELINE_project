
library(tidyverse)
library(RColorBrewer)
library(reshape2)
library(cowplot)

# gene list
gene_list_gain <- c("AKT1", "AKT3", "CCND1", "CCNE2", "CDK6", "ERBB4", "FGFR1", "FGFR2", "MYC", "PIK3CA")
gene_list_gain_label <- paste0(gene_list_gain, "_gain")
gene_list_loss <- c("ESR1", "RB1", "TP53")
gene_list_loss_label <- paste0(gene_list_loss, "_loss")
gene_list <- c(gene_list_gain, gene_list_loss)
gene_list_label <- c(gene_list_gain_label, gene_list_loss_label)

# cnv call
x = read.table("FELINE_FACETS.gene_cnv_call.short_list.STable_CNVcall.txt", header=T, sep="\t")
x$patient <- substr(x$sample, 1, 3)
x %>% select(c("patient", "response", gene_list)) -> x_NR
x_NR
melt(x_NR, id.vars=c("patient", "response")) -> x_NR_m
x_NR_m
## set gain and loss status
x_NR_m$value[x_NR_m$value=="Gain"] = "gain"
x_NR_m$value[x_NR_m$value=="LOH"] = "loss"
x_NR_m$value[x_NR_m$value=="Loss"] = "loss"
## create label
x_NR_m$cnv = paste(x_NR_m$variable, x_NR_m$value, sep="_")
x_NR_m %>% filter(cnv %in% gene_list_label) -> x_NR_m_data
x_NR_m_data <- distinct(x_NR_m_data)
x_NR_m_data
x_NR_m_data_table <- data.frame(table(x_NR_m_data$response, x_NR_m_data$cnv))
names(x_NR_m_data_table) <- c("Response", "CNV", "Num")
x_NR_m_data_table$Label <- paste0(x_NR_m_data_table$CNV, " (", as.integer(100*(x_NR_m_data_table$Num/sum(x_NR_m_data_table$Num))), "%)")
x_NR_m_data_table

# barplot
pdf("FELINE_gene_cnv_barplot.patient.pdf", width=10, height=7)
fontsize=20
p = ggplot(x_NR_m_data_table, aes(x=factor(CNV, levels=gene_list_label), y=Num, fill=Response)) + 
    geom_bar(stat="identity", position=position_dodge()) +
    theme_bw()+theme(axis.line = element_line(colour = "black"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_rect(colour='white'),
                 panel.background = element_blank()) +
    theme(axis.text=element_text(size=fontsize, color='black'),
      axis.text.x =element_text(size=fontsize, angle = 35, color='black', hjust = 0.8, vjust = 0.8),
      axis.title.x =element_text(size=fontsize, face="bold"),
      axis.title.y =element_text(size=fontsize, face="bold"),
      strip.text = element_text(size=fontsize, face="bold")) +
    theme(legend.text = element_text(size=fontsize),
      legend.title = element_blank(),
      legend.position = "top") +
    theme(plot.title=element_text(size=fontsize-4,face="bold"), axis.text=element_text(size=fontsize, face="bold"))
p
dev.off()




