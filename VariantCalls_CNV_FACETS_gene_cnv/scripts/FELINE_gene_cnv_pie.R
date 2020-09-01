
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

# non-responder
x %>% filter(response %in% c("Non-responder")) %>% select(c("response", gene_list)) -> x_NR
melt(x_NR, id.vars=c("response")) -> x_NR_m
## set gain and loss status
x_NR_m$value[x_NR_m$value=="Gain"] = "gain"
x_NR_m$value[x_NR_m$value=="LOH"] = "loss"
x_NR_m$value[x_NR_m$value=="Loss"] = "loss"
## create label
x_NR_m$cnv = paste(x_NR_m$variable, x_NR_m$value, sep="_")
x_NR_m %>% filter(cnv %in% gene_list_label) -> x_NR_m_data
x_NR_m_data_table <- data.frame(table(x_NR_m_data$cnv))
names(x_NR_m_data_table) <- c("CNV", "Num")
x_NR_m_data_table$Label <- paste0(x_NR_m_data_table$CNV, " (", as.integer(100*(x_NR_m_data_table$Num/sum(x_NR_m_data_table$Num))), "%)")
x_NR_m_data_table

# responder
x %>% filter(response %in% c("Responder")) %>% select(c("response", gene_list)) -> x_R
melt(x_R, id.vars=c("response")) -> x_R_m
## set gain and loss status
x_R_m$value[x_R_m$value=="Gain"] = "gain"
x_R_m$value[x_R_m$value=="LOH"] = "loss"
x_R_m$value[x_R_m$value=="Loss"] = "loss"
## create label
x_R_m$cnv = paste(x_R_m$variable, x_R_m$value, sep="_")
x_R_m %>% filter(cnv %in% gene_list_label) -> x_R_m_data
x_R_m_data_table <- data.frame(table(x_R_m_data$cnv))
names(x_R_m_data_table) <- c("CNV", "Num")
x_R_m_data_table$Label <- paste0(x_R_m_data_table$CNV, " (", as.integer(100*(x_R_m_data_table$Num/sum(x_R_m_data_table$Num))), "%)")
x_R_m_data_table


##pei chart
blank_theme <- theme_minimal()+
  theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.border = element_blank(),
  panel.grid=element_blank(),
  axis.ticks = element_blank(),
  plot.title=element_text(size=14, face="bold")
  )

# color
mycolors = colorRampPalette(brewer.pal(8, "Set2"))(length(x_NR_m_data_table$CNV))
x_NR_m_data_table$colors = mycolors

# non-responder label
x_NR_m_data_table <- x_NR_m_data_table %>% 
  mutate(end = 2 * pi * cumsum(Num)/sum(Num),
         start = lag(end, default = 0),
         middle = 0.5 * (start + end),
         hjust = ifelse(middle > pi, 1, 0),
         vjust = ifelse(middle < pi/2 | middle > 3 * pi/2, 0, 1))
x_NR_m_data_table

# Responder label
x_R_m_data_table <- x_R_m_data_table %>%
  mutate(end = 2 * pi * cumsum(Num)/sum(Num),
         start = lag(end, default = 0),
         middle = 0.5 * (start + end),
         hjust = ifelse(middle > pi, 1, 0),
         vjust = ifelse(middle < pi/2 | middle > 3 * pi/2, 0, 1))
x_R_m_data_table


pdf("FELINE_gene_cnv_pie.pdf", width=8, height=7)
# non-responder pie chart
pie_NR <- ggplot(x_NR_m_data_table, aes(x=1, y = Num, fill = CNV)) + geom_bar(stat = "identity")
pie_NR <- pie_NR + coord_polar(theta='y') + blank_theme + scale_fill_manual(values=x_NR_m_data_table$color)
pie_NR <- pie_NR + scale_y_continuous(breaks=cumsum(x_NR_m_data_table$Num) - x_NR_m_data_table$Num / 2, labels= x_NR_m_data_table$Label)
cumsum(x_NR_m_data_table$Num) - x_NR_m_data_table$Num / 2
x_NR_m_data_table$Label
pie_NR_nolegend <- pie_NR + theme(legend.position = "none") 
pie_NR
# responder pie chart
pie_R <- ggplot(x_R_m_data_table, aes(x=1, y = Num, fill = CNV)) + geom_bar(stat = "identity", orientation="y")
pie_R <- pie_R + coord_polar(theta='y', direction = 1) + blank_theme + scale_fill_manual(values=x_NR_m_data_table$color)
pie_R <- pie_R + scale_y_continuous(breaks=cumsum(x_R_m_data_table$Num) - x_R_m_data_table$Num / 2, labels= x_R_m_data_table$Label)
cumsum(x_R_m_data_table$Num) - x_R_m_data_table$Num / 2
x_R_m_data_table$Label
pie_R_nolegend <- pie_R + theme(legend.position = "none")
pie_R
# merge plot
plot_grid(pie_NR_nolegend, pie_R_nolegend, nrow=1, ncol=2, labels=c("NR", "R"), label_size=20)
dev.off()





