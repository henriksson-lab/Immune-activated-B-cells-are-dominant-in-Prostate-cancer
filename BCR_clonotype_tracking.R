### Clonotype tracking

# Load the package and data
library(immunarch)
library(pals)

immdata_K <- repLoad("/corgi/martin/R/IGK/")
immdata_L <- repLoad("/corgi/martin/R/IGL/")
immdata_K1 <- repLoad("/corgi/martin/R/p1_IGK/")
immdata_L1 <- repLoad("/corgi/martin/R/p1_IGL/")

## Tracking the top 5 most abundant IgK and IgL clonotypes per sample in all three patients 
# IgK
target_K <- c("CQQYYSTPRTF", "CQQYYSTPWTF", "CQQYNSYPYTF", "CQQSYSAPPTF", "CQQTYSTPRTF", "CQQYNNWPPWTF", "CQQYYSTPLTF", "CQQYYSTLSLTF", "CQQSYSTPITF", "CQQYNYWPRTF", "CQQRSNWPLTF", "CQQSYSTPQTF", "CQQYNSYPWTF", "CMQALQTPRTF", "CQQYGSSPITF", "CLQYNTYPQYTF")
tc_K <- trackClonotypes(immdata_K$data, target_K, .col = "aa")

p1 <- vis(tc_K) + scale_fill_manual(values=as.vector(cols25(20))) + theme(plot.title = element_blank(), legend.key.size = unit(0.3, 'cm'), legend.key.height = unit(0.3, 'cm'), legend.key.width = unit(0.3, 'cm'), legend.text = element_text(size=7))
ggsave(filename = "FIGURE_IGK_track.png", plot = p1, width = 5, height = 4.5, dpi = 600)

# IgL
target_L <- c("CQTWGTGILVF", "CCSYAGSYTFEWVF", "CQSADNIGPYVVF", "CAAWDDSLSGVF", "CQSYDSSLSGSAVF", "CCSYAGSSTWVF", "CQSADSSGTVVF", "CQTWGTGIQVF", "CCSYAGSSTYVF", "CSSYTSSSTWVF", "CSSYAGSIVVF", "CQSYDSSLSGVVF", "CYSTDSSGTRVF", "CATWDDSLSGVVF", "CQSADSSGAVIF")
tc_L <- trackClonotypes(immdata_L$data, target_L, .col = "aa")

p2 <- vis(tc_L) + scale_fill_manual(values=as.vector(stepped(15))) + theme(plot.title = element_blank(), legend.key.size = unit(0.3, 'cm'), legend.key.height = unit(0.3, 'cm'), legend.key.width = unit(0.3, 'cm'), legend.text = element_text(size=7))
ggsave(filename = "FIGURE_IGL_track.png", plot = p2, width = 5, height = 4.5, dpi = 600)

## Tracking the top 5 most abundant IgK and IgL clonotypes per sample in SN, N-SN and P samples of patient 1
# IgK
target_K1 <- c("CQQYYSTPWTF", "CQQYNSYPYTF", "CQQSYSAPPTF", "CQQTYSTPRTF", "CQQYNNWPPWTF", "CQQYYSTPLTF", "CQQYYSTLSLTF", "CQQSYSTPITF", "CQQYNYWPRTF", "CQQRSNWPLTF", "CQQSYSTPQTF", "CQQYNSYPWTF", "CMQALQTPRTF", "CQQYGSSPITF", "CLQYNTYPQYTF")
tc_K1 <- trackClonotypes(immdata_K1$data, target_K1, .col = "aa")

p3 <- vis(tc_K1) + scale_fill_manual(values=as.vector(cols25(15))) + theme(plot.title = element_blank(), legend.key.size = unit(0.3, 'cm'), legend.key.height = unit(0.3, 'cm'), legend.key.width = unit(0.3, 'cm'), legend.text = element_text(size=7))
ggsave(filename = "FIGURE_p1_IGK_track.png", plot = p3, width = 5, height = 4.5, dpi = 600)

# IgL
target_L1 <- c("CQTWGTGILVF", "CCSYAGSYTFEWVF", "CQSADNIGPYVVF", "CAAWDDSLSGVF", "CQSYDSSLSGSAVF", "CCSYAGSSTWVF", "CQSADSSGTVVF", "CQTWGTGIQVF", "CCSYAGSSTYVF", "CSSYTSSSTWVF", "CSSYAGSIVVF", "CQSYDSSLSGVVF", "CYSTDSSGTRVF", "CATWDDSLSGVVF", "CQSADSSGAVIF")
tc_L1 <- trackClonotypes(immdata_L1$data, target_L1, .col = "aa")

p4 <- vis(tc_L1) + scale_fill_manual(values=as.vector(cols25(15))) + theme(plot.title = element_blank(), legend.key.size = unit(0.3, 'cm'), legend.key.height = unit(0.3, 'cm'), legend.key.width = unit(0.3, 'cm'), legend.text = element_text(size=7))
ggsave(filename = "FIGURE_p1_IGL_track.png", plot = p4, width = 5, height = 4.5, dpi = 600)