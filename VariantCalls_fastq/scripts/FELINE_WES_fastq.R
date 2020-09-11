library(tidyverse)

# COH068
x = read.table("feline_47.sample.txt", header=TRUE, sep="\t")
x %>% select(Filename, Pair, Sample) -> x_filtered
x_filtered$COHID = rep("COH068", length(x_filtered$Filename)) 

# COH077
y = read.table("feline_7patients.txt", header=TRUE, sep="\t")
y %>% select(Filename, Pair, Sample) -> y_filtered
y_filtered$COHID = rep("COH077", length(y_filtered$Filename))

# COH074 
z = read.table("feline_3patients.txt", header=TRUE, sep="\t")
z %>% select(Filename, Pair, Sample) -> z_filtered
z_filtered$COHID = rep("COH077", length(z_filtered$Filename))

data = rbind(x_filtered, y_filtered, z_filtered)
#data %>% select(Sample, Pair, COHID, Filename) %>% arrange(Sample, Filename, Pair, COHID) -> data_sorted
data %>% select(Sample, Pair, COHID, Filename) %>% filter(!Sample %in% c("FEL040_S", "FEL040_M", "FEL040_E")) -> data_sorted
write.table(data_sorted, "FELINE_WES_fastq.list", quote=F, row.names=F, sep="\t")


