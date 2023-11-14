# Author: Yiyun Li
# Affiliation: Wageningen University & Research

# Differential expression analysis
# remove probes with multiple geneID
# use affy_ATH1_array_elements-2010-12-20.txt table for annotation
# use decideTests() to get DEGs list
# did not remove any chips
# in this script the quality test is not included

###install bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.11")

BiocManager::install(c("affy", "simpleaffy", "limma", "GEOquery", "affyQCReport"))


library(affy)
library(simpleaffy)
library(limma)

###### Analysis of GSE46205
setwd("/GSE46205/GSE46205_RAW")

###read the .cel file
#path of .cel file
path =  "/GSE46205/GSE46205_RAW"
# import .cel files containing raw probe-level data
data.raw =ReadAffy(celfile.path=path)
QCReport(data.raw,file="QCreport.pdf")


annotation <- read.table(file = "annotation table.txt", sep = "\t", header = T)
head(annotation)

###Preprocessing
data.gcrma <- gcrma(data.raw)

## transfer to expression set to get expression matix
expset <- exprs(data.gcrma)
any(is.na(expset)) # false
write.table(expset, "expset_gcrma_matrix.txt", quote = F, sep = "\t", row.names = T, col.names = T)


#remove filtered probes
##annotation 
probe <- rownames(expset) 
table(rownames(expset) %in% annotation$array_element_name)

expset = expset[rownames(expset) %in% annotation$array_element_name,]
annotation = annotation[match(rownames(expset), annotation$array_element_name),]


tmp = by(expset,
         annotation$locus,
         function(x) rownames(x)[which.max(rowMeans(x))])
probe = as.character(tmp)
expset = expset[rownames(expset) %in% probe,] 
rownames(expset)= annotation[match(rownames(expset),annotation$array_element_name),2]


##grouping
pData(data.gcrma)$treatment <- grouplist<- gl(8, 3, length = length(data.raw), labels = c("MS_1hr", "MS_48hr", 
                                                                                          "NaCl_1hr", "NaCl_3hr", "NaCl_8hr", "NaCl_20hr", "NaCl_32hr", "NaCl_48hr"))
##functions of differential differential expression analysis means
#subset of expression data
expsub <- function(x, y){
  expset.sub = as.data.frame(expset)[,seq(x,y)] ## only keep the epi data
  return(expset.sub)
} 

#subset of normalized data
datasub <- function(x,y){
  data.gcrma.sub <- data.gcrma[,sampleNames(data.raw)[x:y]]
  return(data.gcrma.sub)
}


data.fit <- function(data.sub, exp.sub){
  #factorize both grouping variables:
  groups = pData(data.sub)$treatment
  f = factor(groups,levels=c("MS_1hr", "MS_48hr","NaCl_1hr", "NaCl_3hr", "NaCl_8hr", "NaCl_20hr", "NaCl_32hr", "NaCl_48hr"))
  #create a design matrix
  design = model.matrix(~ 0 + f, data.sub) 
  colnames(design) = c("MS_1hr", "MS_48hr","NaCl_1hr", "NaCl_3hr", "NaCl_8hr", "NaCl_20hr", "NaCl_32hr", "NaCl_48hr")
  data.fit = lmFit(exp.sub,design)
  ##define a contrast matrix defining the contrasts
  contrast.matrix = makeContrasts(NaCl_1hr-MS_1hr, NaCl_3hr-MS_1hr, NaCl_8hr-MS_1hr, NaCl_20hr-MS_1hr, NaCl_32hr-MS_1hr, NaCl_48hr-MS_1hr, levels = design)
  data.fit.con = contrasts.fit(data.fit,contrast.matrix)
  data.fit.eb = eBayes(data.fit.con)
  
  print(summary(decideTests(data.fit.eb)))
  return(data.fit.eb)
}

deg.all <- function(data.fit.eb, layer){
  DEG <- topTable(data.fit.eb, adjust.method="fdr", coef=1:6, n=Inf, sort.by = "B", lfc = 1, p.value = 0.001)
  write.table(DEG, file = paste("(eachtimepoint-MS)DEG-",layer, "-topTable.txt",sep = ""), row.names = T, col.names = T, sep = "\t", na = "NA")
  return(DEG)
}

usedec <- function(data.fit.eb,layer){
  deg.dec <- decideTests(data.fit.eb, method = "separate", p.value = 0.001, lfc = 1)
  de.common <- which(deg.dec[,1]!="0" | deg.dec[,2]!="0" |deg.dec[,3]!="0"|deg.dec[,4]!="0"|deg.dec[,5]!="0"|deg.dec[,6]!="0")
  deg.dec <- cbind(rownames(deg.dec),deg.dec)
  write.table(deg.dec, file = paste("(eachtimepoint-MS)DEG-",layer, "-dec.txt",sep = ""), row.names = F, col.names = T, sep = "\t", na = "NA")
  deg.listdec <- rownames(deg.dec)[de.common]
  return(deg.listdec)
}

###### Analysis of E-MTAB-2565
setwd("/LRdatasets/LRtimedata/")

###read the .cel file
#path of .cel file
path =  "/LRdatasets/LRtimedata/E-MTAB-2565.raw.1"
# import .cel files containing raw probe-level data
data.raw =ReadAffy(celfile.path=path)

###Look at the data
sampleNames(data.raw)
samplename <- sampleNames(data.raw)
samplename <- sub(pattern = "^Voss_1-", "", samplename)
samplename <- sub(pattern = "^[0-9]*_", "", samplename)
samplename <- sub(pattern = "after-lateral-root-initiation_", "", samplename)
samplename <- sub(pattern = "_ATH1.CEL$", "", samplename)
sampleNames(data.raw) <- samplename


#looking at the sample annotation
pData(data.raw)
pData(data.raw)$timepoint <- sub(pattern  = "_Rep[1-9]$", "", samplename)

data.gcrma <- gcrma(data.raw)
expset <- exprs(data.gcrma)
# expset.n <- arrange(expset, colnames(expset))
# expset <- expset[, order(colnames(expset))]

group <- timepoint <- sub(pattern  = "_Rep[1-9]$", "", colnames(expset))


annotation <- read.table(file = "/Users/liyiyun/Desktop/GSE46205/GSE46205_RAW/annotation table.txt", sep = "\t", header = T)

## transfer to expression set to get expression matix
any(is.na(expset)) 
write.table(expset, "expset_gcrma_matrix-LRtime-probe.txt", quote = F, sep = "\t", row.names = T, col.names = T)


#remove filtered probes
##annotation 
probe <- rownames(expset)
table(rownames(expset) %in% annotation$array_element_name)

expset = expset[rownames(expset) %in% annotation$array_element_name,]
annotation = annotation[match(rownames(expset), annotation$array_element_name),]

tmp = by(expset,
         annotation$locus,
         function(x) rownames(x)[which.max(rowMeans(x))])
probe = as.character(tmp)
expset = expset[rownames(expset) %in% probe,] # 过滤有多个探针的基因
rownames(expset)= annotation[match(rownames(expset),annotation$array_element_name),2]
write.table(expset, "/LRdatasets/expset_gcrma_matrix-LRtime-AGI.txt", quote = F, sep = "\t", row.names = T, col.names = T)




#factorize both grouping variables:
groups = pData(data.gcrma)$timepoint
f = factor(groups,levels=c(unique(groups)))
#create a design matrix
design = model.matrix(~ 0 + f, data.gcrma) 
colnames(design) = c(paste("LR",unique(groups),sep = ""))
data.fit = lmFit(expset,design)
##define a contrast matrix defining the contrasts
# contrast.matrix = makeContrasts(LR12h-LR9h,  levels = design)
contrast.matrix = makeContrasts(LR6h-LR0h, LR9h-LR6h, LR12h-LR9h, LR15h-LR12h, LR18h-LR15h, LR21h-LR18h, LR24h-LR21h, 
                                LR27h-LR24h, LR30h-LR27h, LR33h-LR30h,LR36h-LR33h, LR39h-LR36h, LR42h-LR39h, LR45h-LR42h,
                                LR48h-LR45h, LR51h-LR48h, LR54h-LR51h, levels = design)
data.fit.con = contrasts.fit(data.fit,contrast.matrix)
data.fit.eb = eBayes(data.fit.con)


deg.top.1 <- topTable(data.fit.eb, adjust.method="fdr", coef=1:17, n=Inf, sort.by = "B",  p.value = 0.05, lfc = 1)
write.table(deg.top.1, file = "/Users/liyiyun/Desktop/LRdatasets/DEG/DEGlist-LRtime-topTable-eachtoearlier.txt", row.names = T, col.names = T, sep = "\t", na = "NA")


#create a design matrix
design = model.matrix(~ 0 + f, data.gcrma) 
colnames(design) = c(paste("LR",unique(groups),sep = ""))
data.fit = lmFit(expset,design)

##define a contrast matrix defining the contrasts
# contrast.matrix = makeContrasts(LR12h-LR9h,  levels = design)
contrast.matrix = makeContrasts(LR6h-LR0h, LR9h-LR0h, LR12h-LR0h, LR15h-LR0h, LR18h-LR0h, LR21h-LR0h, LR24h-LR0h, 
                                LR27h-LR0h, LR30h-LR0h, LR33h-LR0h,LR36h-LR0h, LR39h-LR0h, LR42h-LR0h, LR45h-LR0h,
                                LR48h-LR0h, LR51h-LR0h, LR54h-LR0h, levels = design)
data.fit.con = contrasts.fit(data.fit,contrast.matrix)
data.fit.eb = eBayes(data.fit.con)

print(summary(decideTests(data.fit.eb)))

deg.top <- topTable(data.fit.eb, adjust.method="fdr", coef=1:17, n=Inf, sort.by = "B",  p.value = 0.05, lfc = 1)
write.table(deg.top, file = "/Users/liyiyun/Desktop/LRdatasets/DEG/DEGlist-LRtime-topTable-to0h.txt", row.names = T, col.names = T, sep = "\t", na = "NA")
