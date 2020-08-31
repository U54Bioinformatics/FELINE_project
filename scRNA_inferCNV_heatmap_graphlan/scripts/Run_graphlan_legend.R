library(RColorBrewer)


clone.colors <- c("#E86A10","#56A4AA", "#3A78C4","#F1AB00","#2F992D", '#8d4891', '#fe9536', '#d7352e')
subclones    <- c(1,2,3,4,5,6,7)

pdf("FELINE_scRNA_tree_legend.pdf", width=7, height=7)
## plot legend
plot(c(5), c(5), xlim=c(0, 10), ylim=c(0, 10))
### timepoint
legend(x=1, y=7, legend=c("Day 0", "Day 14", "Day 180"), border = NULL, title="Timepoint",
   pch=15, col=brewer.pal(n = 3, name = 'PuRd'), pt.cex=2, bty="n", ncol=3)
### subclone
legend(x=1, y=5, legend=subclones[1:length(subclones)-1], border = NULL, title="Subclone",
   pch=15, col=clone.colors[1:length(subclones)-1], pt.cex=2, bty="n", ncol=3)
### gene expression
lgd_ = rep(NA, 3)
lgd_[c(1,6,11)] = c(-2,0,2)
legend(x=4, y=5, legend=rev(lgd_), fill=colorRampPalette(c("red", "white", "blue"))(11),
      border = NULL, title="Normalized gene expression", y.intersp = 0.5, cex = 1.5, text.font = 1, bty="n")

dev.off()
