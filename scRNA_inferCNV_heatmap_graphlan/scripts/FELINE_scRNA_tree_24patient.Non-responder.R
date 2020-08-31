library(ggplot2)
library(cowplot)
library(magick)

p="FEL015P105"
tree0  = ggdraw() + draw_image(paste0("./MEDALT_tree_png/", p, ".tree.png"))
p="FEL020P115"
tree1  = ggdraw() + draw_image(paste0("./MEDALT_tree_png/", p, ".tree.png"))
p="FEL041P137"
tree2  = ggdraw() + draw_image(paste0("./MEDALT_tree_png/", p, ".tree.png"))
p="FEL043P202"
tree3  = ggdraw() + draw_image(paste0("./MEDALT_tree_png/", p, ".tree.png"))
p="FEL044P116"
tree4  = ggdraw() + draw_image(paste0("./MEDALT_tree_png/", p, ".tree.png"))
p="FEL019P113"
tree5  = ggdraw() + draw_image(paste0("./MEDALT_tree_png/", p, ".tree.png"))
p="FEL034P601"
tree6  = ggdraw() + draw_image(paste0("./MEDALT_tree_png/", p, ".tree.png"))
p="FEL013P102"
tree7  = ggdraw() + draw_image(paste0("./MEDALT_tree_png/", p, ".tree.png"))
p="FEL026P140"
tree8  = ggdraw() + draw_image(paste0("./MEDALT_tree_png/", p, ".tree.png"))
p="FEL035P119"
tree9  = ggdraw() + draw_image(paste0("./MEDALT_tree_png/", p, ".tree.png"))
p="FEL036P132"
tree10  = ggdraw() + draw_image(paste0("./MEDALT_tree_png/", p, ".tree.png"))
p="FEL039P201"
tree11  = ggdraw() + draw_image(paste0("./MEDALT_tree_png/", p, ".tree.png"))

pdf("FELINE_scRNA_tree_24patient.Non-responder.pdf", height=22, width=10)
plot_grid(tree0,tree1,tree2,tree3,tree4,tree5,tree6,tree7,tree8,tree9,tree10,tree11, nrow=6, ncol=2)
dev.off()

