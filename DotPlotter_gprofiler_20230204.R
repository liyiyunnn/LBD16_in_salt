# Author: Jasper Lamers (adapted from Maria Chiara)
# Affiliation: Wageningen University & Research

setwd("/Users/liyiyun/Library/CloudStorage/OneDrive-WageningenUniversity&Research/VICI-ENW projects/Data analysis/YanxiaRNAseq2020/20220517_RNAseq_DEA_1timepoint_2treatment_2genotype")

library(ggplot2)
library(dplyr)
library(pheatmap)
library(ggrepel)

#IF
rm(list = ls())
gProfiler <- read.csv("/Users/liyiyun/Library/CloudStorage/OneDrive-WageningenUniversity&Research/VICI-ENW projects/Data analysis/YanxiaRNAseq2020/20220517_RNAseq_DEA_1timepoint_2treatment_2genotype/saltspecific_6h_go.csv")
gProfiler$GeneRatio <- gProfiler$intersection_size/gProfiler$term_size
gProfiler$log_Size <- log10(gProfiler$term_size)
gProfiler$negative_log10_of_adjusted_p_value <- -log10(gProfiler$p_value)
gProfiler$check_show <- c(rep(1,4),rep(0,121),rep(1,4),rep(0,11),rep(1,1),rep(0,18))

gProfiler <- gProfiler %>%  
  mutate(wrap_term_name = stringr::str_wrap(term_name, width = 30))

DP <- ggplot(gProfiler, aes(x=log_Size, y=negative_log10_of_adjusted_p_value)) +
    geom_point(aes(size = GeneRatio)) +
    theme_bw(base_size = 8) +
    #scale_colour_gradient(name="p-value",low="#56b4e9", high ="#000000",limits = c(0,0.05)) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),panel.grid.minor = element_blank(), text = element_text(size =  16)) + scale_size(name="Gene ratio", range=c(0.5,8), breaks = seq(0.1, 0.3, 0.15), limits = c(0, 1)) +
  scale_y_continuous(breaks=seq(0, 14, 2), limits=c(0, 14)) + scale_x_continuous(breaks=seq(1, 5, 1), limits=c(0, 5), labels=c("10","100","1000","10000","100000"),guide = guide_axis(angle = 45)) +
    ylab(bquote(Adjusted~p-value~(-log[10]))) + xlab(bquote(Term~size~(log[10]))) + guides(color=guide_colourbar(title.vjust=3)) + facet_grid(. ~ source) + #, scales = "free_x", space = "free_x") + 
  geom_label_repel(aes(size=0.2,label = ifelse(check_show==1,as.character(wrap_term_name),'')),box.padding = 0.35, point.padding = 0.35,label.padding=0.12,label.size=NA,fill = alpha(c("white"),0),hjust = 0.3,vjust = 0.45,segment.color="black",segment.size=0.2,segment.alpha=0.8,nudge_x = -0.5,nudge_y = 0.1,force = 20,lineheight = 0.75) +
  coord_cartesian(ylim = c(0, 14), xlim = c(1, 5.5), expand = FALSE) +
  theme(plot.margin = unit(c(0.05,0.2,0.05,0.05), "cm"), legend.key.height= unit(0.1, 'cm'),
        legend.key.width= unit(0.1, 'cm'),legend.title = element_text(size=6),legend.text = element_text(size=6))


ggsave(filename = "salt-specific_6h_DotPlot_new_202309.tiff", plot = DP, dpi=400, width =  18, height = 16 , limitsize = FALSE,units = "cm")  



   # DP <- GOenriched %>% ggplot(aes(x = cluster , y = reorder(code, new_order))) +
#   geom_point(aes(size = GeneRatio, color = AverageL2FC)) +
#   theme_bw(base_size = 7) +
#   scale_colour_gradient2(name="Log2FoldChange",low="#0072b2", high ="#d55e00",midpoint=0, mid="#ffffff", limits=c(-0.7,0.7),labels=c(-0.7,-0.35,0,0.35,0.7),breaks=c(-0.7,-0.35,0,0.35,0.7)) + scale_y_discrete(labels=~empty_DF_All[match(.x, empty_DF_All$code),"term"]) +
#   theme(legend.key.size = unit(0.35, 'cm'), axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.minor = element_blank(), text = element_text(size = 7)) + scale_size(name="Gene ratio", range=c(0.1,3.5), breaks = seq(0.1, 1.1, 0.3), limits = c(0, 1)) +
#   ylab(NULL) + xlab("Treatment or mutant") + guides(color=guide_colourbar(title.vjust=3))+ scale_x_discrete(limits=c("CaCl2","cop1-4","etr1-1","moca1","snrk2.2/2.3"), drop=F) + facet_grid(measure ~ ., scales = "free_y", space = "free_y")
# 
# ggsave(filename = "All_v2.pdf", plot = DP, dpi=400, width = 11, height = 15, unit="cm", limitsize = FALSE)
# DP <- NULL