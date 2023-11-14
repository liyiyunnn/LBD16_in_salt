# Author: Jasper Lamers 
# Affiliation: Wageningen University & Research

library(ggplot2)
library(dplyr)
library(pheatmap)
library(ggrepel)

#IF
gProfiler <- read.csv("saltspecific_6h_go.csv")
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

