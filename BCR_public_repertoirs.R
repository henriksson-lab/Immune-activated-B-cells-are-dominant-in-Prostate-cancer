
# Load the package and data
library(immunarch) 
library(RColorBrewer) 

immdata_l <- repLoad("/corgi/martin/R/light_chains/")
immdata_h <- repLoad("/corgi/martin/R/heavy_chains/")

# Checking the number of unique clonotypes
exp_vol_l <- repExplore(immdata_l$data, .method = "volume")
exp_vol_h <- repExplore(immdata_h$data, .method = "volume")

p1 <- vis(exp_vol_l, .title = "")
p2 <- vis(exp_vol_h, .title = "")

# Checking the diversity index D50
div_d50_l <- repDiversity(immdata_l$data, "d50")
div_d50_h <- repDiversity(immdata_h$data, "d50")

p3 <- vis(div_d50_l)
p4 <- vis(div_d50_h)

# Repertoire overlap for light (IgK, IgL) and heavy chains (IgM, IgG) - calculated using "public" method

imm_ov1_l <- repOverlap(immdata_l$data, .method = "public", .verbose = F)
imm_ov1_h <- repOverlap(immdata_h$data, .method = "public", .verbose = F)

p5 <- vis(imm_ov1_l, "heatmap2", .title = "", text.size = 10, .color = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100))
ggsave(filename = "FIGURE_overlap_light.png", plot = p5, width = 4.8, height = 4.5, dpi = 600)

p6 <- vis(imm_ov1_h, "heatmap2", .title = "", text.size = 10, .color = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100))
ggsave(filename = "FIGURE_overlap_heavy.png", plot = p6, width = 4.8, height = 4.5, dpi = 600)