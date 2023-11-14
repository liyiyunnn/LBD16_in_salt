### heatmap for cell wall related genes

##convert gene name
symconvert <- function(genelist){
  anno <- bitr(geneID = genelist, fromType = "TAIR", toType = "SYMBOL" , OrgDb = "org.At.tair.db", drop = F)
  anno <- anno[!duplicated(anno$TAIR), ]
  for (i in c(1:dim(anno)[1])) {
    if (is.na(anno[i,2])) {
      anno[i,2] <- anno[i,1]
    }
  }
  colnames(anno) <- c("GeneID", "SYMBOL")
  return(anno)
}

# read the files
setwd("/Users/liyiyun/Library/CloudStorage/OneDrive-WageningenUniversity&Research/VICI-ENW projects/Data analysis/YanxiaRNAseq2020/20220517_RNAseq_DEA_1timepoint_2treatment_2genotype/cell wall related")
go <- read.csv("/Users/liyiyun/Library/CloudStorage/OneDrive-WageningenUniversity&Research/VICI-ENW projects/Data analysis/YanxiaRNAseq2020/20220517_RNAseq_DEA_1timepoint_2treatment_2genotype/saltspecific_6h_go.csv")
#data is dds normed
data <- read.csv("/Users/liyiyun/Library/CloudStorage/OneDrive-WageningenUniversity&Research/VICI-ENW projects/Data analysis/YanxiaRNAseq2020/20220517_RNAseq_DEA_1timepoint_2treatment_2genotype/Genotype_comparison/dds_normalized.csv")
rownames(data) <- sub("gene:","",data$X)
data <- data[,-1]

pe <- go[grep("pectin", go$term_name),17]
write.csv(go[grep("pectin", go$term_name),],"pectin-related_GO.csv")
pe.gene <- paste(pe[1],pe[2],sep=",")
pe.gene <- unlist(strsplit(pe.gene,","))
pe.gene <- pe.gene[!duplicated(pe.gene)]


xy <- go[grep("xyloglucan", go$term_name),17]
write.csv(go[grep("xyloglucan", go$term_name),],"xyloglucan-related_GO.csv")
xy.gene <- paste(xy[1],xy[2],sep=",")
xy.gene <- unlist(strsplit(xy.gene,","))
xy.gene <- xy.gene[!duplicated(xy.gene)]

cw <- go[grep("cell wall", go$term_name),17]
write.csv(go[grep("cell wall", go$term_name),],"cellwall-related_GO.csv")
cw.gene <- paste(cw[1],cw[2],cw[3],cw[4],cw[5],cw[6],cw[7],cw[8],cw[9],cw[10],cw[11],cw[12],sep=",")
cw.gene <- unlist(strsplit(cw.gene,","))
cw.gene <- cw.gene[!duplicated(cw.gene)]
cw.spec <- setdiff(cw.gene,xy.gene)
cw.spec <- setdiff(cw.spec,pe.gene)

setdiff(pe.gene,cw.gene)
setdiff(xy.gene,cw.gene)


library(pheatmap)

anno.col <- data.frame(GeneID = c(pe.gene, xy.gene, cw.spec),
                       GO = c(rep("pectin-related",length(pe.gene)),rep("xyloglucan-related",length(xy.gene)),rep("ther cell wall-related", length(cw.spec))))
# anno <- symconvert(anno.col$GeneID)
# anno.col <- merge(anno.col,anno,by = "GeneID")
rownames(anno.col) <- anno.col$GeneID
anno.col <- data.frame(GO = c(rep("pectin-related",length(pe.gene)),rep("xyloglucan-related",length(xy.gene)),rep("Other cell wall-related", length(cw.spec))))
rownames(anno.col) <- c(pe.gene, xy.gene, cw.spec)





plot.data <- data[anno.col$GeneID,]
mean.data <- data.frame(Control_Col0 = rowMeans(plot.data[,c(1,2,3)]),
                        Salt_Col0 = rowMeans(plot.data[,c(4,5,6)]),
                        Control_lbd16 = rowMeans(plot.data[,c(7,8,9)]),
                        Salt_lbd16 = rowMeans(plot.data[,c(10,11,12)]))
anno <- anno.col
anno$gene <- rownames(anno)
mean <- mean.data
mean$gene <- rownames(mean)
new <- merge(anno,mean,  by = "gene")
new$GO <- factor(new$GO,levels = c("pectin-related","xyloglucan-related","other cell wall-related"))
new <- new[order(new$GO),]
rownames(new) <- new$gene
new <- new[,-c(1,2)]


pheatmap(t(new), 
         cluster_rows = F,
         cluster_cols = F,
         color = colorRampPalette(c("navy","yellow"))(10),
         show_rownames = T,
         show_colnames = T,
         annotation_col = anno.col,
         scale = "column",
         cluster_row = FALSE,
         gaps_col=c(length(pe.gene),sum(length(pe.gene),length(xy.gene))),
         gaps_row = c(1,2,3),
         cellheight=7,cellwidth=7,
         fontsize_col=6,
         fontsize_row = 10,
         border_color = NA
         # annotation_colors=ann_colors,
         # filename="heatmap_cellwall.png"
         )



### using symbol

anno.col <- data.frame(GeneID = c(pe.gene, xy.gene, cw.spec),
                       GO = c(rep("Pectin-related",length(pe.gene)),rep("Xyloglucan-related",length(xy.gene)),rep("Other cell wall-related", length(cw.spec))))
anno <- symconvert(anno.col$GeneID)
anno.col <- merge(anno,anno.col,by = "GeneID")
rownames(anno.col) <- anno.col$GeneID






plot.data <- data[anno.col$GeneID,]
mean.data <- data.frame(Control_Col0 = rowMeans(plot.data[,c(1,2,3)]),
                        Salt_Col0 = rowMeans(plot.data[,c(4,5,6)]),
                        Control_lbd16 = rowMeans(plot.data[,c(7,8,9)]),
                        Salt_lbd16 = rowMeans(plot.data[,c(10,11,12)]))
anno <- anno.col
anno$gene <- rownames(anno)
mean <- mean.data
mean$gene <- rownames(mean)
new <- merge(anno,mean,  by = "gene")
new$GO <- factor(new$GO,levels = c("Pectin-related","Xyloglucan-related","Other cell wall-related"))
new <- new[order(new$GO),]
new.anno.col <- data.frame(GO = new$GO)
rownames(new.anno.col) <- new$SYMBOL
rownames(new) <- new$SYMBOL
new <- new[,-c(1,2,3,4)]


pheatmap(t(new), 
         cluster_rows = F,
         cluster_cols = F,
         color = colorRampPalette(c("navy","yellow"))(10),
         show_rownames = T,
         show_colnames = T,
         annotation_col = new.anno.col,
         scale = "column",
         cluster_row = FALSE,
         gaps_col=c(length(pe.gene),sum(length(pe.gene),length(xy.gene))),
         gaps_row = c(1,2,3),
         cellheight=7,cellwidth=7,
         fontsize_col=6,
         fontsize_row = 10,
         border_color = NA
         # annotation_colors=ann_colors,
         # filename="heatmap_cellwall.png"
)


#####seperate one
anno.col <- data.frame(GeneID = c(pe.gene, xy.gene),
                       GO = c(rep("Pectin-related",length(pe.gene)),rep("Xyloglucan-related",length(xy.gene))))
anno <- symconvert(anno.col$GeneID)
anno.col <- merge(anno,anno.col,by = "GeneID")
rownames(anno.col) <- anno.col$GeneID






plot.data <- data[anno.col$GeneID,]
mean.data <- data.frame(Control_Col0 = rowMeans(plot.data[,c(1,2,3)]),
                        Salt_Col0 = rowMeans(plot.data[,c(4,5,6)]),
                        Control_lbd16 = rowMeans(plot.data[,c(7,8,9)]),
                        Salt_lbd16 = rowMeans(plot.data[,c(10,11,12)]))
anno <- anno.col
anno$gene <- rownames(anno)
mean <- mean.data
mean$gene <- rownames(mean)
new <- merge(anno,mean,  by = "gene")
new$GO <- factor(new$GO,levels = c("Pectin-related","Xyloglucan-related"))
new <- new[order(new$GO),]
new.anno.col <- data.frame(GO = new$GO)
rownames(new.anno.col) <- new$SYMBOL
rownames(new) <- new$SYMBOL
new <- new[,-c(1,2,3,4)]


pheatmap(t(new), 
         cluster_rows = F,
         cluster_cols = F,
         color = colorRampPalette(c("navy","yellow"))(10),
         show_rownames = T,
         show_colnames = T,
         annotation_col = new.anno.col,
         scale = "column",
         cluster_row = FALSE,
         gaps_col=c(length(pe.gene),sum(length(pe.gene),length(xy.gene))),
         # gaps_row = c(1,2,3),
         cellheight=7,cellwidth=7,
         fontsize_col=6,
         fontsize_row = 8,
         border_color = NA
         # annotation_colors=ann_colors,
         # filename="heatmap_cellwall.png"
)


anno.col <- data.frame(GeneID = c(cw.spec),
                       GO = c(rep("Other cell wall-related", length(cw.spec))))
anno <- symconvert(anno.col$GeneID)
anno.col <- merge(anno,anno.col,by = "GeneID")
rownames(anno.col) <- anno.col$GeneID






plot.data <- data[anno.col$GeneID,]
mean.data <- data.frame(Control_Col0 = rowMeans(plot.data[,c(1,2,3)]),
                        Salt_Col0 = rowMeans(plot.data[,c(4,5,6)]),
                        Control_lbd16 = rowMeans(plot.data[,c(7,8,9)]),
                        Salt_lbd16 = rowMeans(plot.data[,c(10,11,12)]))
anno <- anno.col
anno$gene <- rownames(anno)
mean <- mean.data
mean$gene <- rownames(mean)
new <- merge(anno,mean,  by = "gene")
new$GO <- factor(new$GO,levels = c("Other cell wall-related"))
new <- new[order(new$GO),]
new.anno.col <- data.frame(GO = new$GO)
rownames(new.anno.col) <- new$SYMBOL
rownames(new) <- new$SYMBOL
new <- new[,-c(1,2,3,4)]


pheatmap(t(new), 
         cluster_rows = F,
         cluster_cols = F,
         color = colorRampPalette(c("navy","yellow"))(10),
         show_rownames = T,
         show_colnames = T,
         annotation_col = new.anno.col,
         scale = "column",
         cluster_row = FALSE,


         cellheight=7,cellwidth=7,
         fontsize_col=6,
         fontsize_row = 8,
         border_color = NA
         # annotation_colors=ann_colors,
         # filename="heatmap_cellwall.png"
)
