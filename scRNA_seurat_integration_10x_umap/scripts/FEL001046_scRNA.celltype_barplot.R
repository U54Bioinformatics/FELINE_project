library(data.table)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(cowplot)

barplot_png <- function(filename, title, df, x_lab, y_lab, width=10, height=15){
    pdf(filename, width = width, height = height)
    font_size = 14
    x_level = c("Cancer cells", "Normal epithelial cells", "Stromal cells", "Immune cells")
    p <- ggplot(df, aes(fill=factor(Celltype, levels=x_level), y=Freq, label=sprintf("%0.2f", round(Freq, digits = 2)), x=Patient)) + geom_bar(position="fill", stat="identity") +
    labs(title=title, x=x_lab, y=y_lab) + theme(plot.title=element_text(size=font_size,face="bold"), axis.text=element_text(size=font_size, face="bold"), legend.position = "top", axis.title=element_text(size=font_size,face="bold"), axis.text.x =element_text(size=font_size, angle = 45, color='black', hjust = 1), legend.text=element_text(size=font_size, face="bold"), legend.title=element_blank()) + theme(panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + guides(col = guide_legend(ncol=2), shape = guide_legend(override.aes = list(size = 8), ncol=2)) + theme(plot.margin =  margin(t = 1.5, r = 1.5, b = 1.5, l = 1.5, unit = "cm")) + scale_fill_manual(values=c("Cancer cells" = "#E7298A", "Normal epithelial cells"="#66A61E", "Stromal cells"="#006633", "Immune cells"="light sea green"))

#    facet_grid( ~ Response, labeller = labeller(GeneSet = label_y)) +
#    theme(strip.text = element_text(size = fontsize-8, colour = "black", angle = 0),
#    strip.text.y = element_text(angle = 0),
#    strip.background = element_rect(colour="gray", fill="white"))

    print(p)
    dev.off()
}


x = fread("FEL011046.cell_metadata.txt")
x %>% filter(Platform %in% "10x") -> x
data.frame(table(x$orig.ident, x$Celltype))
data_df = data.frame(table(x$orig.ident, x$Celltype))
names(data_df) <- c("Patient", "Celltype", "Freq")

filename <- "FEL001046_scRNA.celltype.barplot.pdf"
barplot_png(filename, "", data_df, "Patient", "Proportion", width=10, height=4)



