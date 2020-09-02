library(ggplot2)
library(cowplot)
library(magick)

p="FEL015"
cnv="sequenza"
tree0  = ggdraw() + draw_image(paste0(p, "_pyclone_analysis_",cnv ,"_20_5_0.05_1_10000_none.clonevol_ccf.clonevol_figure4_fishplot.png"))
p="FEL020"
cnv="FACETS"
tree1  = ggdraw() + draw_image(paste0(p, "_pyclone_analysis_",cnv ,"_20_5_0.05_1_10000_none.clonevol_ccf.clonevol_figure4_fishplot.png"))
p="FEL041"
cnv="FACETS"
tree2  = ggdraw() + draw_image(paste0(p, "_pyclone_analysis_",cnv ,"_20_5_0.05_1_10000_none.clonevol_ccf.clonevol_figure4_fishplot.png"))
p="FEL043"
cnv="FACETS"
tree3  = ggdraw() + draw_image(paste0(p, "_pyclone_analysis_",cnv ,"_20_5_0.05_1_10000_none.clonevol_ccf.clonevol_figure4_fishplot.png"))
p="FEL044"
cnv="sequenza"
tree4  = ggdraw() + draw_image(paste0(p, "_pyclone_analysis_",cnv ,"_20_5_0.05_1_10000_none.clonevol_ccf.clonevol_figure4_fishplot.png"))
p="FEL019"
cnv="FACETS"
tree5  = ggdraw() + draw_image(paste0(p, "_pyclone_analysis_",cnv ,"_20_5_0.05_1_10000_none.clonevol_ccf.clonevol_figure4_fishplot.png"))
p="FEL034"
cnv="FACETS"
tree6  = ggdraw() + draw_image(paste0(p, "_pyclone_analysis_",cnv ,"_20_5_0.05_1_10000_none.clonevol_ccf.clonevol_figure4_fishplot.png"))
p="FEL013"
cnv="FACETS"
tree7  = ggdraw() + draw_image(paste0(p, "_pyclone_analysis_",cnv ,"_20_5_0.05_1_10000_none.clonevol_ccf.clonevol_figure4_fishplot.png"))
p="FEL026"
cnv="sequenza"
tree8  = ggdraw() + draw_image(paste0(p, "_pyclone_analysis_",cnv ,"_20_5_0.05_1_10000_none.clonevol_ccf.clonevol_figure4_fishplot.png"))
p="FEL035"
cnv="sequenza"
tree9  = ggdraw() + draw_image(paste0(p, "_pyclone_analysis_",cnv ,"_20_5_0.05_1_10000_none.clonevol_ccf.clonevol_figure4_fishplot.png"))
p="FEL036"
cnv="FACETS"
tree10  = ggdraw() + draw_image(paste0(p, "_pyclone_analysis_",cnv ,"_20_5_0.05_1_10000_none.clonevol_ccf.clonevol_figure4_fishplot.png"))
p="FEL039"
cnv="sequenza"
tree11  = ggdraw() + draw_image(paste0(p, "_pyclone_analysis_",cnv ,"_20_5_0.05_1_10000_none.clonevol_ccf.clonevol_figure4_fishplot.png"))
p="FEL014"
cnv="sequenza"
tree12  = ggdraw() + draw_image(paste0(p, "_pyclone_analysis_",cnv ,"_20_5_0.05_1_10000_none.clonevol_ccf.clonevol_figure4_fishplot.png"))
p="FEL023"
cnv="sequenza"
tree14  = ggdraw() + draw_image(paste0(p, "_pyclone_analysis_",cnv ,"_20_5_0.05_1_10000_none.clonevol_ccf.clonevol_figure4_fishplot.png"))
p="FEL028"
cnv="FACETS"
tree15  = ggdraw() + draw_image(paste0(p, "_pyclone_analysis_",cnv ,"_20_5_0.05_1_10000_none.clonevol_ccf.clonevol_figure4_fishplot.png"))
p="FEL045"
cnv="sequenza"
tree17  = ggdraw() + draw_image(paste0(p, "_pyclone_analysis_",cnv ,"_20_5_0.05_1_10000_none.clonevol_ccf.clonevol_figure4_fishplot.png"))
p="FEL022"
cnv="FACETS"
tree18  = ggdraw() + draw_image(paste0(p, "_pyclone_analysis_",cnv ,"_20_5_0.05_1_10000_none.clonevol_ccf.clonevol_figure4_fishplot.png"))
p="FEL025"
cnv="FACETS"
tree19  = ggdraw() + draw_image(paste0(p, "_pyclone_analysis_",cnv ,"_20_5_0.05_1_10000_none.clonevol_ccf.clonevol_figure4_fishplot.png"))
p="FEL027"
cnv="FACETS"
tree20  = ggdraw() + draw_image(paste0(p, "_pyclone_analysis_",cnv ,"_20_5_0.05_1_10000_none.clonevol_ccf.clonevol_figure4_fishplot.png"))
p="FEL037"
cnv="sequenza"
tree22  = ggdraw() + draw_image(paste0(p, "_pyclone_analysis_",cnv ,"_20_5_0.05_1_10000_none.clonevol_ccf.clonevol_figure4_fishplot.png"))
p="FEL046"
cnv="sequenza"
tree23  = ggdraw() + draw_image(paste0(p, "_pyclone_analysis_",cnv ,"_20_5_0.05_1_10000_none.clonevol_ccf.clonevol_figure4_fishplot.png"))
p="FEL031"
cnv="sequenza"
tree25  = ggdraw() + draw_image(paste0(p, "_pyclone_analysis_",cnv ,"_20_5_0.05_1_10000_none.clonevol_ccf.clonevol_figure4_fishplot.png"))

pdf("FELINE_fishplot_by_arm_response.pdf", height=22, width=16)
plot_grid(tree0,tree1,tree2,tree3,tree4,tree5,tree6,tree7,tree8,tree9,tree10,tree11,tree12,tree14,tree15,tree17,tree18,tree19,tree20,tree22,tree23,tree25, nrow=6, ncol=4, labels=c("P15 (Arm A, NR)","P20 (Arm A, NR)","P41 (Arm A, NR)","P43 (Arm A, NR)","P44 (Arm A, NR)","P19 (Arm B, NR)","P34 (Arm B, NR)","P13 (Arm C, NR)","P26 (Arm C, NR)","P35 (Arm C, NR)","P36 (Arm C, NR)","P39 (Arm C, NR)","P14 (Arm A, R)","P23 (Arm A, R)","P28 (Arm A, R)","P45 (Arm A, R)","P22 (Arm B, R)","P25 (Arm B, R)","P27 (Arm B, R)","P37 (Arm B, R)","P46 (Arm B, R)","P31 (Arm C, R)"), label_size=20)
dev.off()

