####BiocManager::install("DESeq2")
#BiocManager::install("BSgenome.Athaliana.TAIR.TAIR9", version = "3.8")
#BiocManager::install("tximport")
#BiocManager::install("SummarizedExperiment")


library(stringr)
library("clusterProfiler")
setwd("/Users/liyiyun/Desktop/lbd16-1/20220517_RNAseq_DEA_1timepoint_2treatment_2genotype")

#Find Salmon files and add file name to a new csv
setup <- read.csv("/Users/liyiyun/Desktop/lbd16-1/0_Input/RNAseq.csv", header = TRUE)
setup$filename <- list.files("/Users/liyiyun/Desktop/lbd16-1/0_Input/Salmon", recursive = TRUE, include.dirs = T, pattern = ".sf")
setup$Experiment <- ifelse(grepl("W", setup$Sample), "96", "Root")
write.csv(setup, "RNAseq2.csv")#, row.names=FALSE)

#Import all geneIDs for transcriptIDs and remove all n/a
if (!exists("gff")){
  library(tibble)
  library(rtracklayer)
  gff_file = file.path("/Users/liyiyun/Desktop/lbd16-1/Arabidopsis_thaliana.TAIR10.42.gff3")
  gff = import(gff_file)
}



comparison_summary <- setup[(setup$Experiment == "Root" &  setup$Timepoint == "6hrs" &
                               (setup$Genotype == "Col-0" | setup$Genotype == "lbd16-1") &
                               (setup$Treatment == "Control" | setup$Treatment == "NaCl (130mM)" )) , ]

sampleNames <- comparison_summary$filename
sampleGroup1 <- comparison_summary$Genotype
sampleGroup2 <- comparison_summary$Treatment
sampleGroup3 <- comparison_summary$Timepoint
sampleGroup4 <- comparison_summary$Sample
sampleTable <- data.frame(sampleName = sampleNames, Genotype = sampleGroup1, Treatment = sampleGroup2, Timepoint = sampleGroup3, Experiment = sampleGroup4)
sampleTable$Treatment <- sub(")","",sampleTable$Treatment)
sampleTable$Treatment <- sub("\\(","_",sampleTable$Treatment)
sampleTable$Treatment <- sub(" ","",sampleTable$Treatment)

sampleTable <- cbind(sampleTable, Condition = paste(sampleTable$Treatment, sampleTable$Timepoint, sep="_"))



#Find all samples
samples <- read.table("RNAseq2.csv", header=TRUE, sep= ",")

#Full name of file and checks whether it exists
files <- file.path("/Users/liyiyun/Desktop/lbd16-1/0_Input/Salmon", comparison_summary$filename)
files
names(files) <- paste(comparison_summary$Treatment, comparison_summary$Genotype,  comparison_summary$Timepoint, sep="_")
all(file.exists(files))

#Write GeneIDs at TRanscripts levels (sum all splice variants)
tx2gene <- tibble(txid = gff$transcript_id, geneid = as.character(gff$Parent))# ) %>% na.omit()
tx2gene2 <- na.omit(tx2gene)

#Make summary files and import tx (from Salmon)
library(tximport)
library(SummarizedExperiment)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene2)#, txout=TRUE)
txi_counts = txi[2]

#DESeq2
library(DESeq2)
ddsTxi <- DESeqDataSetFromTximport(txi, sampleTable, design = ~ Genotype + Condition + Genotype:Condition)
ddsTxi$Genotype <- relevel(ddsTxi$Genotype, "Col-0")
ddsTxi$Condition <- relevel(ddsTxi$Condition, "Control_6hrs")
dds <- dds[rowSums(counts(dds)) > 10,]
dds <- DESeq(ddsTxi)
Folder_name = "Genotype_comparison"
dir.create(Folder_name, showWarnings = FALSE)
write.csv(txi_counts, file=paste(Folder_name,"/","counts_sum_before_DESeq2.csv", sep=""))
write.csv(sampleTable, file=paste(Folder_name,"/","SampleTable.csv", sep=""))

#topGene <- rownames(resSig)[which.max(resSig$log2FoldChange)]
#data <- plotCounts(dds, gene=topGene, intgroup=c("Treatment"), returnData=TRUE)
#ggplot(data, aes(x=Treatment, y=count, fill=Treatment)) + ggtitle(topGene) + scale_y_log10() + geom_dotplot(binaxis="y", stackdir="center")

#Dispention graph
jpeg(paste(Folder_name,"/",Folder_name,"_Dispention.jpeg", sep=""),quality=100,width = 500, height = 500, pointsize=10)
plotDispEsts(dds)
dev.off()

# Pair-wise DEGs on total normalized dds
for(TR in resultsNames(dds)[-1]) {
  #Create folder for results 
  TR_subfolder <- paste(Folder_name,TR,sep="/" )
  dir.create(TR_subfolder, showWarnings = FALSE)
  
  #Annotate results file with short name AT and ENTREZID
  res <- results(dds, list( c("Genotype_lbd16.1_vs_Col.0",TR)))
  library("AnnotationDbi")
  library("org.At.tair.db")
  convertIDs <- function( ids, from, to, db, ifMultiple=c("putNA", "useFirst")) {
    stopifnot( inherits( db, "AnnotationDb" ) )
    ifMultiple <- match.arg( ifMultiple )
    suppressWarnings( selRes <- AnnotationDbi::select(
      db, keys=ids, keytype=from, columns=c(from,to) ) )
    if ( ifMultiple == "putNA" ) {
      duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ]
      selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ]
    }
    return( selRes[ match( ids, selRes[,1] ), 2 ] )
  }
  
  
  res$GeneID = gsub("gene:", "\\1", row.names(res))
  res$hgnc_symbol <- convertIDs(res$GeneID, "TAIR", "SYMBOL", org.At.tair.db)
  res$entrezgene <- convertIDs(res$GeneID, "TAIR", "ENTREZID", org.At.tair.db)
  
  
  #Write summary Results
  sink(paste(TR_subfolder,"/",TR,"_DESeq2_result_summary.txt", sep=""), type=c("output", "message"), append = FALSE)
  summary(results(dds, name=TR, alpha=0.01))
  resSig <- subset(res, padj < 0.01)
  paste("DEG = ", nrow(resSig), " Genes (p<0.01)", sep="")
  sink(NULL)
  write.csv(res, paste(TR_subfolder,"/",TR,"_DEGs_DESeq2.csv", sep=""))
  
  res.table <-res[complete.cases(res[,c(2,6)]),]
  #  res.table <- res.table[abs(res.table$log2FoldChange)>=1,]
  res.table <- res.table[res.table$padj<0.05,]
  write.csv(res.table, paste(TR_subfolder,"/",TR,"_signDEGs_DESeq2.csv", sep=""))
  
  up.deg <- res.table[res.table$log2FoldChange>0,]
  write.csv(up.deg, paste(TR_subfolder,"/",TR,"_signDEGs_DESeq2_upregulated.csv", sep=""))
  
  down.deg <- res.table[res.table$log2FoldChange<0,]
  write.csv(down.deg, paste(TR_subfolder,"/",TR,"_signDEGs_DESeq2_downregulated.csv", sep=""))
  
  
  if (dim(down.deg)[1] != 0 & dim(up.deg)[1]!=0){
    #GO enrichment
    library("clusterProfiler")
    go.deg <- enrichGO(res.table$GeneID , OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)
    go.deg <- clusterProfiler::simplify(go.deg, cutoff=0.7, by='p.adjust', select_fun = min)
    write.csv(go.deg, paste(TR_subfolder,"/",TR,"_signDEGs_DESeq2_GO.csv", sep=""))
    png(filename = paste(TR_subfolder,"/",TR,"_signDEGs_DESeq2_GO.png", sep=""), width = 600, height = 1000, pointsize = 15)
    barplot(go.deg, showCategory=60,title=paste(TR_subfolder,"/",TR,"_signDEGs_DESeq2", sep=""))
    dev.off()
    
    go.updeg <- enrichGO(up.deg$GeneID , OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)
    go.updeg <- clusterProfiler::simplify(go.updeg, cutoff=0.7, by='p.adjust', select_fun = min)
    write.csv(go.updeg, paste(TR_subfolder,"/",TR,"_signDEGs_DESeq2_upregulated_GO.csv", sep=""))
    png(filename = paste(TR_subfolder,"/",TR,"_signDEGs_DESeq2_upregulated_GO.png", sep=""), width = 600, height = 1000, pointsize = 15)
    barplot(go.updeg, showCategory=60,title=paste(TR_subfolder,"/",TR,"_signDEGs_DESeq2_upregulated", sep=""))
    dev.off()
    
    go.downdeg <- enrichGO(down.deg$GeneID , OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)
    go.downdeg <- clusterProfiler::simplify(go.downdeg, cutoff=0.7, by='p.adjust', select_fun = min)
    write.csv(go.downdeg, paste(TR_subfolder,"/",TR,"_signDEGs_DESeq2_downregulated_GO.csv", sep=""))
    png(filename = paste(TR_subfolder,"/",TR,"_signDEGs_DESeq2_downregulated_GO.png", sep=""), width = 600, height = 1000, pointsize = 15)
    barplot(go.downdeg, showCategory=60,title=paste(TR_subfolder,"/",TR,"_signDEGs_DESeq2_downregulated", sep=""))
    dev.off()
    
    
    
    go.deg <- enrichGO(res.table$GeneID , OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP", pAdjustMethod = "BH", qvalueCutoff = 0.05,readable      = T)
    go.deg <- clusterProfiler::simplify(go.deg, cutoff=0.7, by='p.adjust', select_fun = min)
    write.csv(go.deg, paste(TR_subfolder,"/",TR,"_signDEGs_DESeq2_GO_readable.csv", sep=""))
    
    go.updeg <- enrichGO(up.deg$GeneID , OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP", pAdjustMethod = "BH", qvalueCutoff = 0.05,readable      = T)
    go.updeg <- clusterProfiler::simplify(go.updeg, cutoff=0.7, by='p.adjust', select_fun = min)
    write.csv(go.updeg, paste(TR_subfolder,"/",TR,"_signDEGs_DESeq2_upregulated_GO_readable.csv", sep=""))
    
    go.downdeg <- enrichGO(down.deg$GeneID , OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP", pAdjustMethod = "BH", qvalueCutoff = 0.05,readable      = T)
    go.downdeg <- clusterProfiler::simplify(go.downdeg, cutoff=0.7, by='p.adjust', select_fun = min)
    write.csv(go.downdeg, paste(TR_subfolder,"/",TR,"_signDEGs_DESeq2_downregulated_GO_readable.csv", sep=""))
    
    # go.deg <- enrichGO(res.table$GeneID , OrgDb = "org.At.tair.db", keyType="TAIR", ont="ALL", pAdjustMethod = "BH", qvalueCutoff = 0.05)
    # write.csv(go.deg, paste(TR_subfolder,"/",TR,"_signDEGs_DESeq2_GO_ALL.csv", sep=""))
    # 
    # go.updeg <- enrichGO(up.deg$GeneID , OrgDb = "org.At.tair.db", keyType="TAIR", ont="ALL", pAdjustMethod = "BH", qvalueCutoff = 0.05)
    # write.csv(go.updeg, paste(TR_subfolder,"/",TR,"_signDEGs_DESeq2_upregulated_GO_ALL.csv", sep=""))
    # 
    # go.downdeg <- enrichGO(down.deg$GeneID, OrgDb = "org.At.tair.db", keyType="TAIR", ont="ALL", pAdjustMethod = "BH", qvalueCutoff = 0.05)
    # write.csv(go.downdeg, paste(TR_subfolder,"/",TR,"_signDEGs_DESeq2_downregulated_GO_ALL.csv", sep=""))
  }
  
  #Volcano plot
  res_LFC <- lfcShrink(dds, coef=TR, type="apeglm" )
  write.csv(res_LFC, paste(TR_subfolder,"/",TR,"_resLFC.csv", sep=""))
  
  jpeg(paste(TR_subfolder,"/",TR,"_MA.jpeg", sep=""),quality=100)
  plotMA(res_LFC, ylim=c(-8,8), alpha=0.01, main=TR)
  dev.off()
  
  #Histogram
  jpeg(paste(TR_subfolder,"/",TR,"_HIST.jpeg", sep=""),quality=100,width = 500, height = 500, pointsize=10)
  hist(res$pvalue[res$baseMean > 1], breaks=20, col="grey50", border="white")
  dev.off()
  
  
}


#Write output
write.csv(assay(dds), paste(Folder_name,"/","dds.csv", sep=""))

#Write log_output
rld <- rlog(dds)
write.csv(assay(rld), paste(Folder_name,"/","dds_log.csv", sep=""))

#Dividing counts by size_factor
dds_SF <- estimateSizeFactors(dds)
dds_normalized <- counts(dds_SF, normalized=TRUE)
write.csv(dds_normalized, paste(Folder_name,"/","dds_normalized.csv", sep=""))

#Make heatmap
library( "gplots" )
library( "RColorBrewer" )

sampleDists <- dist( t( assay(rld) ) )
sampleDistMatrix <- as.matrix( sampleDists )

File_info <- data.frame(do.call('rbind', strsplit(as.character(rld$sampleName),'_', fixed=TRUE)))
File_info$Total <- paste(File_info$X2, File_info$X3, sep="-")

rownames(sampleDistMatrix) <- paste( rld$Genotype, rld$Condition,File_info$Total, sep=" - " )
colnames(sampleDistMatrix) <- paste( rld$Treatment, File_info$Total, sep=" - " )
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

jpeg(paste(Folder_name,"/",Folder_name,"_Heatmap.jpeg", sep=""),quality=100,width = 2000, height = 2000, pointsize=10)
heatmap.2(sampleDistMatrix, col=colours, key=TRUE, symkey=FALSE, density.info="none",cexRow=1.6,cexCol=1.6,margins=c(30,30),trace="none",srtCol=45)
graphics.off()


library(ggplot2)
pcaData <- plotPCA(rld, intgroup = "Genotype", returnData=TRUE, ntop = 27586)
percentVar <- round(100 * attr(pcaData, "percentVar"))

max_no <- ceiling(max(rbind(abs(pcaData$PC1),abs(pcaData$PC2))) * 1.1)

ggplot(pcaData, aes(PC1, PC2, color= Treatment, shape = Genotype)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_color_manual(values=c("#cc79a7", "#e69f00", "#56b4e9", "#009e73", "#f0e442", "#0072b2", "#d55e00")) +
  coord_fixed() + ylim(c(-max_no,max_no)) + xlim(c(-max_no,max_no))+
  theme_bw()#+
  #geom_text(aes(label= rld@colData$Experiment), vjust = 0, color = "black")
ggsave(paste(Folder_name,"/",Folder_name,"_PCA_all.jpeg", sep=""), plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 15, height = 15, 
       dpi = 500, limitsize = TRUE)

save(rld,dds, ddsTxi, sampleTable,
     file = "run_20220426.RData")

#summary
first_file_name <- list.files("Genotype_comparison")[c(8:19)]       
dir <- paste("./Genotype_comparison/",first_file_name,sep = "")
n <- length(dir)     

n_sub <- rep(0,n)   
n_sub <- as.data.frame(n_sub)    
n_sub <- t(n_sub)     

data <- read.csv("./Genotype_comparison/Genotype_ds2.3_vs_Col.0/Genotype_ds2.3_vs_Col.0_signDEGs_DESeq2_downregulated.csv")
for (i in 1:n) {   
  b=list.files(dir[i])
  b<- b[grep("_signDEGs_DESeq2",b)]
  n_sub[i]=length(b)    
  for (j in 1:n_sub[i]) {
    file=paste(dir[i],"/",b[j],sep = "")
    new_data <- read.csv(file)
    print(sub(dir[i],"",file))
    print(dim(new_data))
    # new_data <- new_data[-1,]   
    # names(new_data) <- NULL   
  }
}

##salt specific
salt.6 <-read.csv("Genotype_comparison/Genotypelbd16.1.ConditionNaCl_130mM_6hrs/Genotypelbd16.1.ConditionNaCl_130mM_6hrs_signDEGs_DESeq2.csv",row.names = 1)
salt.24 <- read.csv("Genotype_comparison/Genotypelbd16.1.ConditionNaCl_130mM_24hrs/Genotypelbd16.1.ConditionNaCl_130mM_24hrs_signDEGs_DESeq2.csv",row.names = 1)
con.6 <- read.csv("Genotype_comparison/Genotype_lbd16.1_vs_Col.0/Genotype_lbd16.1_vs_Col.0_signDEGs_DESeq2.csv",row.names = 1)
con.24 <- read.csv("Genotype_comparison/Genotypelbd16.1.ConditionControl_24hrs/Genotypelbd16.1.ConditionControl_24hrs_signDEGs_DESeq2.csv",row.names = 1)

salt.spec.6 <- setdiff(rownames(salt.6),rownames(con.6))
salt.spec.6 <- salt.6[salt.spec.6,]
dim(salt.spec.6)
write.csv(salt.spec.6,"saltspecific_6h.csv")

salt.spec.6.up <- salt.spec.6[salt.spec.6$log2FoldChange>0,]
dim(salt.spec.6.up)
write.csv(salt.spec.6.up,"saltspecific_6h_upregulated.csv")
salt.spec.6.down <- salt.spec.6[salt.spec.6$log2FoldChange<0,]
dim(salt.spec.6.down)
write.csv(salt.spec.6.down,"saltspecific_6h_downregulated.csv")

go.deg <- enrichGO(salt.spec.6$GeneID , OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP", pAdjustMethod = "BH", qvalueCutoff = 0.05,readable      = T)
go.deg <- clusterProfiler::simplify(go.deg, cutoff=0.7, by='p.adjust', select_fun = min)
write.csv(go.deg, "6h_signDEGs_DESeq2_saltspecific_GO_readable.csv")

go.updeg <- enrichGO(salt.spec.6.up$GeneID , OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP", pAdjustMethod = "BH", qvalueCutoff = 0.05,readable      = T)
go.updeg <- clusterProfiler::simplify(go.updeg, cutoff=0.7, by='p.adjust', select_fun = min)
write.csv(go.updeg, "6h_signDEGs_DESeq2_saltspecific_upregulated_GO_readable.csv")

go.downdeg <- enrichGO(salt.spec.6.down$GeneID , OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP", pAdjustMethod = "BH", qvalueCutoff = 0.05,readable      = T)
go.downdeg <- clusterProfiler::simplify(go.downdeg, cutoff=0.7, by='p.adjust', select_fun = min)
write.csv(go.downdeg,"6h_signDEGs_DESeq2_saltspecific_downregulated_GO_readable.csv")



go.deg <- enrichGO(salt.spec.6$GeneID , OrgDb = "org.At.tair.db", keyType="TAIR", ont="ALL", pAdjustMethod = "BH", qvalueCutoff = 0.05,readable      = T)
write.csv(go.deg, "6h_signDEGs_DESeq2_saltspecific_GO_ALL_readableL.csv")

go.updeg <- enrichGO(salt.spec.6.up$GeneID , OrgDb = "org.At.tair.db", keyType="TAIR", ont="ALL", pAdjustMethod = "BH", qvalueCutoff = 0.05,readable      = T)
write.csv(go.updeg, "6h_signDEGs_DESeq2_saltspecific_upregulated_GO_AL_readable.csv")


go.downdeg <- enrichGO(salt.spec.6.down$GeneID , OrgDb = "org.At.tair.db", keyType="TAIR", ont="ALL", pAdjustMethod = "BH", qvalueCutoff = 0.05,readable      = T)
write.csv(go.downdeg,"6h_signDEGs_DESeq2_saltspecific_downregulated_GO_AL_readable.csv")
png(filename = "6h_signDEGs_DESeq2_saltspecific_downregulated_GO_readable.png", width = 600, height = 800, pointsize = 15)
barplot(go.downdeg, showCategory=60,title="6h_signDEGs_DESeq2_saltspecific_downregulated_GO")
dev.off()


salt.spec.24 <- setdiff(rownames(salt.24),rownames(con.24))
salt.spec.24 <- salt.24[salt.spec.24,]
dim(salt.spec.24)
write.csv(salt.spec.24,"saltspecific_24h.csv")

salt.spec.24.up <- salt.spec.24[salt.spec.24$log2FoldChange>0,]
dim(salt.spec.24.up)
salt.spec.24.down <- salt.spec.24[salt.spec.24$log2FoldChange<0,]
dim(salt.spec.24.down)

go.deg <- enrichGO(salt.spec.24$GeneID , OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP", pAdjustMethod = "BH", qvalueCutoff = 0.05,readable      = T)
go.deg <- clusterProfiler::simplify(go.deg, cutoff=0.7, by='p.adjust', select_fun = min)
write.csv(go.deg, "24h_signDEGs_DESeq2_saltspecific_GO_readable.csv")

go.updeg <- enrichGO(salt.spec.24.up$GeneID , OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP", pAdjustMethod = "BH", qvalueCutoff = 0.05,readable      = T)
go.updeg <- clusterProfiler::simplify(go.updeg, cutoff=0.7, by='p.adjust', select_fun = min)
write.csv(go.updeg, "24h_signDEGs_DESeq2_saltspecific_upregulated_GO_readable.csv")

go.downdeg <- enrichGO(salt.spec.24.down$GeneID , OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP", pAdjustMethod = "BH", qvalueCutoff = 0.05,readable      = T)
go.downdeg <- clusterProfiler::simplify(go.downdeg, cutoff=0.7, by='p.adjust', select_fun = min)
write.csv(go.downdeg,"24h_signDEGs_DESeq2_saltspecific_downregulated_GO_readable.csv")

go.deg <- enrichGO(salt.spec.24$GeneID , OrgDb = "org.At.tair.db", keyType="TAIR", ont="ALL", pAdjustMethod = "BH", qvalueCutoff = 0.05,readable      = T)
write.csv(go.deg, "24h_signDEGs_DESeq2_saltspecific_GO_ALL_readableL.csv")

go.updeg <- enrichGO(salt.spec.24.up$GeneID , OrgDb = "org.At.tair.db", keyType="TAIR", ont="ALL", pAdjustMethod = "BH", qvalueCutoff = 0.05,readable      = T)
write.csv(go.updeg, "24h_signDEGs_DESeq2_saltspecific_upregulated_GO_AL_readable.csv")

go.downdeg <- enrichGO(salt.spec.24.down$GeneID , OrgDb = "org.At.tair.db", keyType="TAIR", ont="ALL", pAdjustMethod = "BH", qvalueCutoff = 0.05,readable      = T)
write.csv(go.downdeg,"24h_signDEGs_DESeq2_saltspecific_downregulated_GO_AL_readable.csv")


