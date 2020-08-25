library(clonevol)
library(fishplot)
library(tidyverse)
library(reshape2)

args <- commandArgs(trailingOnly = TRUE)
prefix=args[1]
patient=args[2]

#prefix="FEL020_pyclone_analysis_FACETS_20_5_0.05_1_10000_none.clonevol_ccf"
load(paste0(prefix, '.clonevol_figure4_fishplot.RData'))
#load("FEL020_pyclone_analysis_FACETS_20_5_0.05_1_10000_none.clonevol_ccf.clonevol_figure4_fishplot.RData")
data_m = fishplot_input$cell.fractions[[1]]
data_p = fishplot_input$parents[[1]]

# reformat data
data_m = data.frame(data_m)
data_m$Cluster = rownames(data_m)
data_m$Parent = data_p
data_mm=melt(data=data_m, id.vars=c("Cluster", "Parent"), measure.vars=sample.names)
## patient id
data_mm$Patient = substring(data_mm$variable, 1, 6)
## timepoint
data_mm$Timepoint[str_detect(data_mm$variable, "_S")] = "0"
data_mm$Timepoint[str_detect(data_mm$variable, "_M")] = "14"
data_mm$Timepoint[str_detect(data_mm$variable, "_E")] = "180"
##reorder
data_mm = data_mm[c("Cluster", "Patient", "Timepoint", "variable", "value", "Parent")]
names(data_mm) = c("Cluster", "Patient", "Timepoint", "Sample", "Fraction", "Parent")
##write into file
ofile = paste0(prefix, ".cluster_frequency.txt")
write.table(data_mm, ofile, sep="\t", row.names=F, quote=F)

