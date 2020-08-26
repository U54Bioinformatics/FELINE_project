library("data.table")
library("ggplot2")
library("circlize")
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

prefix=args[1]

clinic = read.table("FELINE_patient_1_46.clinical_data.txt", sep="\t", header=T)
# reformat patient name
clinic$Sample = substring(clinic$Sample, 1, 6)
clinic$Sample = gsub('FEL0', 'P', clinic$Sample)
clinic %>% select(Sample, ARM, Response) %>% arrange(Response, ARM, Sample) -> clinic_ordered

print(clinic_ordered)

data <- fread(paste0(prefix, ".data.txt"), sep="\t", header=T)
data_m <- melt(data, id.var="Sample")
names(data_m) <- c("ID", "Subtype", "Proportion")
em_gsea_h.plot <- data_m
#em_gsea_h.plot$NES_abs <- abs(em_gsea_h.plot$NES)
#em_gsea_h.plot$NES_change[with(em_gsea_h.plot, NES > 0)] <- "Positive"
#em_gsea_h.plot$NES_change[with(em_gsea_h.plot, NES < 0)] <- "Negative"
head(em_gsea_h.plot)
em_gsea_h.plot$GS <- em_gsea_h.plot$ID
head(em_gsea_h.plot)
# reformat patient name
em_gsea_h.plot$GS = substring(em_gsea_h.plot$GS, 1, 6)
em_gsea_h.plot$GS = gsub('FEL0', 'P', em_gsea_h.plot$GS)
head(em_gsea_h.plot)

# rank by sample names
#my_level <- rev(clinic$Sample)
# rank by response then arm then sample names
my_level <- clinic_ordered$Sample
subtype_order = c("S_Basal", "S_LumA", "S_LumB", "S_Her2+", "M_Basal", "M_LumA", "M_LumB", "M_Her2+", "E_Basal", "E_LumA", "E_LumB", "E_Her2+")
print(my_level)
em_gsea_h.plot$Subtype <- factor(em_gsea_h.plot$Subtype, levels = subtype_order)
em_gsea_h.plot$GS <- factor(em_gsea_h.plot$GS, levels = my_level)
head(em_gsea_h.plot)

fontsize=12
pdf(paste0(prefix, ".heatmap.pdf"), height=10, width=8)
p <- ggplot(em_gsea_h.plot, mapping=aes(x=Subtype, y=GS)) +
#geom_point(shape=19) + scale_size_continuous(range = c(0,5)) +
geom_tile(aes(fill = Proportion), color='white') + 
#scale_fill_gradient2(low = "steelblue", high = "tomato") +
scale_fill_gradient2(low = "blue", mid="white", high = "red", na.value='gray80') + 
labs(title=prefix, x="", y = "") +
#scale_color_manual(values=c("blue", "red")) +
theme_classic() +
theme_bw()+theme(axis.line = element_line(colour = "black"),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
panel.border = element_rect(colour = "black", fill=NA, size=1),
panel.background = element_blank()) +
theme(axis.text=element_text(size=fontsize, color='black'),
axis.text.x =element_text(size=fontsize, angle = 65, color='black', hjust = 0, vjust = 1),
axis.title.x =element_blank(),
axis.title.y =element_text(size=fontsize, face="plain"),
axis.title   =element_text(size=fontsize, face="plain"),
strip.text = element_text(face = 'plain', size=12)) +
theme(legend.text = element_text(size=fontsize),
legend.title = element_blank(),
legend.position = "bottom") +
theme(plot.title=element_text(size=fontsize, face="plain", hjust = -0.7, vjust = 0), axis.text=element_text(size=fontsize, face="bold")) +
theme(plot.margin =  margin(t = 1.5, r = 4, b = 1.5, l = 1, unit = "cm")) +
scale_x_discrete(position = "top")
#guides(colour = guide_legend(override.aes = list(size=4), ncol=c(2)))
print(p)
dev.off()

