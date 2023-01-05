### Basic statistics for patients 1, 2, and 3

# Load the package and data
library(immunarch)  

immdata_all <- repLoad("/corgi/martin/BCR_Viqar/BCR_v1/mixcr_outputdata/")
immdata_p1 <- repLoad("/corgi/martin/R/patient1/")
immdata_p2 <- repLoad("/corgi/martin/R/patient2/")
immdata_p3 <- repLoad("/corgi/martin/R/patient3/")

# Number of unique clonotypes
exp_vol_1 <- repExplore(immdata_p1$data, .method = "volume")
exp_vol_2 <- repExplore(immdata_p2$data, .method = "volume")
exp_vol_3 <- repExplore(immdata_p3$data, .method = "volume")

p1 <- vis(exp_vol_1, .title = "")
ggsave(filename = "FIGURE_clonotype_p1.png", plot = p1, width = 4.8, height = 4.5, dpi = 600)

p2 <- vis(exp_vol_2, .title = "")
ggsave(filename = "FIGURE_clonotype_p2.png", plot = p2, width = 4.8, height = 4.5, dpi = 600)

p3 <- vis(exp_vol_3, .title = "")
ggsave(filename = "FIGURE_clonotype_p3.png", plot = p3, width = 4.8, height = 4.5, dpi = 600)

# Clonotype abundancy
imm_hom_1 <- repClonality(immdata_p1$data,
  .method = "homeo",
  .clone.types = c(Small = .0001, Medium = .001, Large = .01, Hyperexpanded = 1)
)

imm_hom_2 <- repClonality(immdata_p2$data,
  .method = "homeo",
  .clone.types = c(Small = .0001, Medium = .001, Large = .01, Hyperexpanded = 1)
)

imm_hom_3 <- repClonality(immdata_p3$data,
  .method = "homeo",
  .clone.types = c(Small = .0001, Medium = .001, Large = .01, Hyperexpanded = 1)
)

p4 <- vis(imm_hom_1)
ggsave(filename = "FIGURE_abundancy_p1.png", plot = p4, width = 6, height = 5, dpi = 600)

p5 <- vis(imm_hom_2)
ggsave(filename = "FIGURE_abundancy_p2.png", plot = p5, width = 6, height = 5, dpi = 600)

p6 <- vis(imm_hom_3)
ggsave(filename = "FIGURE_abundancy_p3.png", plot = p6, width = 6, height = 5, dpi = 600)

# Calculating the diversity index D50
div_d50_1 <- repDiversity(immdata_p1$data, "d50")
div_d50_2 <- repDiversity(immdata_p2$data, "d50")
div_d50_3 <- repDiversity(immdata_p3$data, "d50")

p7 <- vis(div_d50_1)
ggsave(filename = "FIGURE_div50_p1.png", plot = p7, width = 4.8, height = 4.5, dpi = 600)

p8 <- vis(div_d50_2)
ggsave(filename = "FIGURE_div50_p2.png", plot = p8, width = 4.8, height = 4.5, dpi = 600)

p9 <- vis(div_d50_3)
ggsave(filename = "FIGURE_div50_p3.png", plot = p9, width = 4.8, height = 4.5, dpi = 600)

# generating Table S2 (public repertoire table using CDR3 aminoacid sequences and V alleles from merged from all three patients)
pr_aav <- pubRep(immdata_all$data, "aa+v", .verbose = F)
write.table(pr_aav, file="/corgi/martin/R/pr_aav.txt")
