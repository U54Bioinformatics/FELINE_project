savlocFigS2<-"/Users/jason/Jason Griffiths Dropbox/jason griffiths/FELINE Project (1)/Manuscript  Feline immune communication/Nature communications submission docs/Revision and submission folder/Source Data/SI data/FigureS2/"
LogitCompData <- (read.csv(file=paste0(savlocFigS2,"SourceData_FigureS2_LogitTumorCompositionAnalysis_Output.csv")))


distmat <- vegdist(LogitCompData%>%select(-X), method = "euclidean")
pheatmap::pheatmap(as.matrix(distmat),#annotation_col=rowannot,
                   annotation_row=rowannot,
                   annotation_colors = annotation_col,
                   border_color=NA,
                   cluster_rows = T,cluster_cols = T,
                   labels_row=rep("",500),
                   labels_col=rep("",500),
                   cutree_rows = 5,cutree_cols = 5)


