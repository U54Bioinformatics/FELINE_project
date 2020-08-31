library(HelloRanges)
library(igraph)
library(RColorBrewer)
library(dplyr)

#outdir="./outputRNAT_500"
args <- commandArgs(trailingOnly = TRUE)
patient=args[1]
subclone_dir=args[2]
outdir=args[3]
expr_rds=args[4]

## filenames
patient_orig=substr(patient,0,6)
patient_subclone=paste0(subclone_dir, "/", patient, "_HMM_infer.subclone.anno.txt")


## random seed to have same plot
set.seed(123)

## load meta to set colors
meta = read.table(patient_subclone, sep="\t", header=T)
meta <- meta[,c("Cell.ID", "Sample", "HMM_infer")]
names(meta) <- c("id", "Sample", "Subclone")
root_row =c("root", "root", "root")
names(root_row) <- c("id", "Sample", "Subclone")
root_row <- data.frame(t(root_row))
meta <- rbind(meta, root_row)

## load gene expression (zinbwave)
#count = read.table("FEL015P105_HMM_infer.subclone.anno.txt", sep="\t", header=T)
genes = c("ESR1", "FGFR2", "ERBB4", "CDK6", "FOS", "JUNB", "RORA", "MAPKAP1", "MAP3K5", "MAPK10", "MAP3K1", "MAP3K14", "MAP4K3", "MAP2K5", "MAP3K20", "MAPK3", "MAPK1")
count = readRDS(expr_rds)
count = count[,c("Cell.ID", genes)]
names(count) <- c("id", genes)
count$Sample = do.call(rbind, strsplit(count$id, "_"))[,1]
print(head(count)) 
print(tail(count))
count %>% filter(Sample %in% c(patient_orig)) -> count
print(head(count))
print(tail(count))

root_row = c("root", rep(0,ncol(count)-1))
names(root_row) <- c("id", genes)
count <- rbind(count, root_row)
print(head(count)) 
print(tail(count))

## load tree from MEDALT
treeName=paste0(outdir, "/", patient, ".CNV.tree.txt")
celltree=read.csv(treeName,sep="\t")
nodes=data.frame(id=union(as.character(celltree[,1]),as.character(celltree[,2])),size=5)
nodes_raw = nodes

## set subclone color
nodes = nodes_raw 
nodes = merge(nodes, meta, on="id")
clone.colors <- c("#E86A10","#56A4AA", "#3A78C4","#F1AB00","#2F992D", '#8d4891', '#fe9536', '#d7352e')
subclones <- sort(unique(nodes$Subclone))
print(subclones)
for (i in 1:length(subclones)){
    print(i)
    print(subclones[i])
    print(clone.colors[i])
    if (subclones[i] == "root"){
        # set root as black
        nodes$color_subclone[nodes$Subclone=="root"] = "black"
    }else{
        # set subclone color using clone.colors
        nodes$color_subclone[nodes$Subclone==subclones[i]] = clone.colors[i]
    }
}
#print(nodes)

## set sample color
#nodes = nodes_raw
#nodes= merge(nodes, meta, on="id")
nodes$color_sample[nodes$Sample==paste0(patient, "_S")] = brewer.pal(n = 3, name = 'Dark2')[1]
nodes$color_sample[nodes$Sample==paste0(patient, "_M")] = brewer.pal(n = 3, name = 'Dark2')[2]
nodes$color_sample[nodes$Sample==paste0(patient, "_E")] = brewer.pal(n = 3, name = 'Dark2')[3]
nodes$color_sample[nodes$Sample=="root"] = "black"
#print(nodes)
nodes$Sample <- NULL
nodes$Subclone <- NULL

################################################################
## function to map color to values
map2color<-function(x,pal,limits=NULL){
    if(is.null(limits)) limits=range(x)
    pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}
library(circlize)
heat_cols =colorRamp2(c(-4, 0, 4), c("slateblue", "white", "tomato"))
################################################################
## set gene color
print("nodes dimension before merging with count")
print(dim(nodes))
print(head(nodes))
print(tail(nodes))
nodes = merge(nodes, count, on="id")
print("nodes dimension after merging with count")
print(dim(nodes))
print(head(nodes))
print(tail(nodes))
#mypal = brewer.pal(n = 8, name = "RdBu")
#nodes$color_ESR1 = map2color(nodes$ESR1, mypal,c(-10, 10))
#nodes$color_ESR1    = heat_cols(as.numeric(as.vector(nodes$ESR1)))
#nodes$color_FGFR2   = heat_cols(as.numeric(as.vector(nodes$FGFR2)))
#nodes$color_MAP3K8  = heat_cols(as.numeric(as.vector(nodes$MAP3K8)))
#nodes$color_MAP3K20 = heat_cols(as.numeric(as.vector(nodes$MAP3K20)))
#print("nodes dimension after merging with count, color sets")
#print(dim(nodes))
#print(head(nodes))
#print(tail(nodes))

## plot sc tree with sample colors
#nodes$color="lightblue"
#nodes$color[nodes$id==setdiff(as.character(celltree[,1]),as.character(celltree[,2]))]="black"
net <- graph_from_data_frame(d=celltree, vertices=nodes, directed=T)
pdf(file=paste0(outdir, "/", patient, ".singlecell.tree.sample_color.pdf"), width = 9, height = 6, useDingbats = F)
par(mfrow=c(2, 3), mar=c(1,1,5,1))
## set coords so we can plot different annotation with same tree topology
#coords <- layout_with_fr(net)
coords <- layout_in_circle(net)

## plot sample annotation
plot(net, vertex.frame.color=NA,vertex.color=nodes$color_sample,edge.arrow.size=.2,vertex.label=NA, layout=coords)
legend(x=-1, y=1.5, legend=c("Day 0", "Day 14", "Day 180"), border = NULL, title=patient,
   pch=21, pt.bg=brewer.pal(n = 3, name = 'Dark2'), pt.cex=2, bty="n", ncol=3)

## plot subclone annotation
plot(net, vertex.frame.color=NA,vertex.color=nodes$color_subclone,edge.arrow.size=.2,vertex.label=NA, layout=coords)
#color_set = brewer.pal(n = 8, name = "RdBu")
#color_set = "green"
#plot(net, vertex.frame.color=NA,vertex.color=color_set,edge.arrow.size=.0,vertex.label=NA, layout=coords)
legend(x=-1, y=1.5, legend=subclones[1:length(subclones)-1], border = NULL, title="Subclone",
   pch=21, pt.bg=clone.colors[1:length(subclones)-1], pt.cex=2, bty="n", ncol=3)

## plot gene expression
# legend text position
lgd_ = rep(NA, 11)
lgd_[c(1,6,11)] = c(-4,0,4)
for (gene in genes){
    # color 
    nodes$color_gene    = heat_cols(as.numeric(as.vector(nodes[[gene]])))
    # ESR1
    plot(net, vertex.frame.color=NA,vertex.color=nodes$color_gene,edge.arrow.size=.2,vertex.label=NA,layout=coords)
    legend(x=-1.5, y=1.5, legend=lgd_, fill=colorRampPalette(c("slateblue", "white", "tomato"))(11),
      border = NULL, title=gene, y.intersp = 0.5, cex = 1.5, text.font = 1, bty="n")
}
# other genes
#plot(net, vertex.frame.color=NA,vertex.color=nodes$color_FGFR2,edge.arrow.size=.0,vertex.label=NA, layout=coords)
#title="FGFR2"
#legend(x=-1.5, y=1.5, legend=lgd_, fill=colorRampPalette(c("slateblue", "white", "tomato"))(11),
#      border = NULL, title=title, y.intersp = 0.5, cex = 1.5, text.font = 1, bty="n")
#plot(net, vertex.frame.color=NA,vertex.color=nodes$color_MAP3K8,edge.arrow.size=.0,vertex.label=NA, layout=coords)
#title="MAP3K8"
#legend(x=-1.5, y=1.5, legend=lgd_, fill=colorRampPalette(c("slateblue", "white", "tomato"))(11),
#      border = NULL, title=title, y.intersp = 0.5, cex = 1.5, text.font = 1, bty="n")
#plot(net, vertex.frame.color=NA,vertex.color=nodes$color_MAP3K20,edge.arrow.size=.0,vertex.label=NA, layout=coords)
#title="MAP3K20"
#legend(x=-1.5, y=1.5, legend=lgd_, fill=colorRampPalette(c("slateblue", "white", "tomato"))(11),
#      border = NULL, title=title, y.intersp = 0.5, cex = 1.5, text.font = 1, bty="n")
dev.off()


#write.table(nodes,"test.node.txt",quote=F,sep="\t", row.names=F)
#write.table(meta,"test.meta.txt",quote=F,sep="\t", row.names=F)
