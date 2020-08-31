library(tidyverse)
library(RColorBrewer)

args <- commandArgs(trailingOnly = TRUE)
patient=args[1]

# cell meta and gene expression data
data = readRDS("FEL011046_data_gene_pathway.RDS")
cells = read.table(paste0("./MEDALT_cnv/", patient, ".HMMi6.final.CNV.names.txt"))
data %>% select(Cell.ID, Day, Subclone, ESR1, FGFR2, ERBB4, RORA, CDK6) -> data_sub
#data_sub %>% filter (Cell.ID %in% cells$V1) -> data_sub
data_sub %>% select (Cell.ID, Day, Subclone) -> data_anno
data_sub %>% select (ESR1, FGFR2, ERBB4, RORA, CDK6) -> data_expr
data_expr = scale(as.matrix(data_expr), scale=TRUE, center=TRUE)
data_expr = data.frame(data_expr)
head(data_anno)
head(data_expr)
data_expr$Cell.ID = data_anno$Cell.ID
data_expr$Day     = data_anno$Day
data_expr$Subclone= data_anno$Subclone
data_sub = data_expr
# filter cells to analyze
data_sub %>% filter (Cell.ID %in% cells$V1) -> data_sub
#
print("final dataframe")
head(data_sub)

library(circlize)
#heat_cols   =  colorRamp2(c(-2, 0, 2), c("slateblue", "white", "tomato"))
heat_cols   =  colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
color_esr1  =  heat_cols(as.numeric(as.vector(data_sub$ESR1)))
color_fgfr2 =  heat_cols(as.numeric(as.vector(data_sub$FGFR2)))
color_erbb4 =  heat_cols(as.numeric(as.vector(data_sub$ERBB4)))
color_rora  =  heat_cols(as.numeric(as.vector(data_sub$RORA)))
color_cdk6  =  heat_cols(as.numeric(as.vector(data_sub$CDK6)))

## clade anno
clade_name = c()
clade_type = c()
clade_color = c()

## ring annot
day_clade_name  = c()
day_ring_option = c()
day_ring_level  = c()
day_ring_option_value = c()
#col: Cell.ID, Day, Subclone, ESR1, FGFR2, ERBB4, RORA, CDK6
for (i in 1:nrow(data_sub)) {
    cell = data_sub[i, "Cell.ID"]
    day  = data_sub[i, "Day"]
    subclone = data_sub[i, "Subclone"]
    esr1     = data_sub[i, "ESR1"]
    ## clade size
    clade_name  = append(clade_name, cell)
    clade_type  = append(clade_type, "clade_marker_size")
    clade_color = append(clade_color, 0.1)
    ## clade color 
    clade_name  = append(clade_name, cell)
    clade_type  = append(clade_type, "clade_marker_color")
    clade_color = append(clade_color, "#808080")    

    ring_height_value = 3
    ## timepoint
    color_timepoint = "black"
    if(day == 0) {
        color_timepoint = brewer.pal(n = 3, name = 'PuRd')[1]
    }else if (day == 14) {
        color_timepoint = brewer.pal(n = 3, name = 'PuRd')[2]
    }else if (day == 180) {
        color_timepoint = brewer.pal(n = 3, name = 'PuRd')[3]
    }
    ### ring_height
    day_clade_name  = append(day_clade_name, cell)
    day_ring_option = append(day_ring_option, "ring_height")
    day_ring_level  = append(day_ring_level, 1)
    day_ring_option_value = append(day_ring_option_value, ring_height_value+0.3)
    # ring_color
    day_clade_name  = append(day_clade_name, cell)
    day_ring_option = append(day_ring_option, "ring_color")
    day_ring_level  = append(day_ring_level, 1)
    day_ring_option_value = append(day_ring_option_value, color_timepoint)
    ## subclone
    clone.colors   = c("#E86A10","#56A4AA", "#3A78C4","#F1AB00","#2F992D", '#8d4891', '#fe9536', '#d7352e')
    color_subclone = "black"
    subclone_n     = as.numeric(gsub("Cluster", "", subclone))
    color_subclone = clone.colors[subclone_n] 
    ### ring_height
    day_clade_name  = append(day_clade_name, cell)
    day_ring_option = append(day_ring_option, "ring_height")
    day_ring_level  = append(day_ring_level, 2)
    day_ring_option_value = append(day_ring_option_value, ring_height_value)
    # ring_color
    day_clade_name  = append(day_clade_name, cell)
    day_ring_option = append(day_ring_option, "ring_color")
    day_ring_level  = append(day_ring_level, 2)
    day_ring_option_value = append(day_ring_option_value, color_subclone)
    ## ESR1
    day_clade_name  = append(day_clade_name, cell)
    day_ring_option = append(day_ring_option, "ring_color")
    day_ring_level  = append(day_ring_level, 3)
    day_ring_option_value = append(day_ring_option_value, color_esr1[i])
    ### ring_height
    day_clade_name  = append(day_clade_name, cell)
    day_ring_option = append(day_ring_option, "ring_height")
    day_ring_level  = append(day_ring_level, 3)
    day_ring_option_value = append(day_ring_option_value, ring_height_value)
    ## CDK6
    day_clade_name  = append(day_clade_name, cell)
    day_ring_option = append(day_ring_option, "ring_color")
    day_ring_level  = append(day_ring_level, 4)
    day_ring_option_value = append(day_ring_option_value, color_cdk6[i])
    ### ring_height
    day_clade_name  = append(day_clade_name, cell)
    day_ring_option = append(day_ring_option, "ring_height")
    day_ring_level  = append(day_ring_level, 5)
    day_ring_option_value = append(day_ring_option_value, ring_height_value)
    ## FGFR2
    day_clade_name  = append(day_clade_name, cell)
    day_ring_option = append(day_ring_option, "ring_color")
    day_ring_level  = append(day_ring_level, 5)
    day_ring_option_value = append(day_ring_option_value, color_fgfr2[i])
    ### ring_height
    day_clade_name  = append(day_clade_name, cell)
    day_ring_option = append(day_ring_option, "ring_height")
    day_ring_level  = append(day_ring_level, 4)
    day_ring_option_value = append(day_ring_option_value, ring_height_value)
    ## ERBB4
    day_clade_name  = append(day_clade_name, cell)
    day_ring_option = append(day_ring_option, "ring_color")
    day_ring_level  = append(day_ring_level, 6)
    day_ring_option_value = append(day_ring_option_value, color_erbb4[i])
    ### ring_height
    day_clade_name  = append(day_clade_name, cell)
    day_ring_option = append(day_ring_option, "ring_height")
    day_ring_level  = append(day_ring_level, 6)
    day_ring_option_value = append(day_ring_option_value, ring_height_value)
    ## RORA
    day_clade_name  = append(day_clade_name, cell)
    day_ring_option = append(day_ring_option, "ring_color")
    day_ring_level  = append(day_ring_level, 7)
    day_ring_option_value = append(day_ring_option_value, color_rora[i])
    ### ring_height
    day_clade_name  = append(day_clade_name, cell)
    day_ring_option = append(day_ring_option, "ring_height")
    day_ring_level  = append(day_ring_level, 7)
    day_ring_option_value = append(day_ring_option_value, ring_height_value)
}


# output annot data frame
annot.df = data.frame(day_clade_name, day_ring_option, day_ring_level, day_ring_option_value)
clade.df = data.frame(clade_name, clade_type, clade_color)
write.table(clade.df, paste0(patient, ".annot.txt"), sep="\t", quote=F, row.names=F, col.names=F)
write.table(annot.df, paste0(patient, ".annot.txt"), sep="\t", quote=F, row.names=F, col.names=F, append = TRUE)

