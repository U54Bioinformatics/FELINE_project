library(ggplot2)
library(cowplot)
library(magick)

p="FEL014P103"
tree12  = ggdraw() + draw_image(paste0("./MEDALT_tree_png/", p, ".tree.png"))
p="FEL016P107"
tree13  = ggdraw() + draw_image(paste0("./MEDALT_tree_png/", p, ".tree.png"))
p="FEL023P121"
tree14  = ggdraw() + draw_image(paste0("./MEDALT_tree_png/", p, ".tree.png"))
p="FEL028P131"
tree15  = ggdraw() + draw_image(paste0("./MEDALT_tree_png/", p, ".tree.png"))
p="FEL029P134"
tree16  = ggdraw() + draw_image(paste0("./MEDALT_tree_png/", p, ".tree.png"))
p="FEL045P122"
tree17  = ggdraw() + draw_image(paste0("./MEDALT_tree_png/", p, ".tree.png"))
p="FEL022P124"
tree18  = ggdraw() + draw_image(paste0("./MEDALT_tree_png/", p, ".tree.png"))
p="FEL025P142"
tree19  = ggdraw() + draw_image(paste0("./MEDALT_tree_png/", p, ".tree.png"))
p="FEL027P125"
tree20  = ggdraw() + draw_image(paste0("./MEDALT_tree_png/", p, ".tree.png"))
p="FEL021P118"
tree21  = ggdraw() + draw_image(paste0("./MEDALT_tree_png/", p, ".tree.png"))
p="FEL031P143"
tree22  = ggdraw() + draw_image(paste0("./MEDALT_tree_png/", p, ".tree.png"))
p="FEL040P701"
tree23  = ggdraw() + draw_image(paste0("./MEDALT_tree_png/", p, ".tree.png"))

pdf("FELINE_scRNA_tree_24patient.Responder.pdf", height=22, width=10)
plot_grid(tree12,tree13,tree14,tree15,tree16,tree17,tree18,tree19,tree20,tree21,tree22,tree23, nrow=6, ncol=2)
dev.off()

