#### Quality control 
##PCA

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.17")
BiocManager::install(c("affy", "simpleaffy", "limma", "affyPLM"))

library(affy)
library(simpleaffy)
library(limma)
library(affyPLM)
setwd("/GSE46205/GSE46205_RAW")

###read the .cel file
#path of .cel file
path =  "/GSE46205/GSE46205_RAW"
# import .cel files containing raw probe-level data
data.raw =ReadAffy(celfile.path=path)


###Look at the data
sampleNames(data.raw)
samplename <- sampleNames(data.raw)
samplename <- sub(pattern = "^GSM[0-9]{7}_", "", samplename)
samplename <- sub(pattern = ".CEL.gz$", "", samplename)
sampleNames(data.raw) <- samplename

data.gcrma <- gcrma(data.raw)
expset <- exprs(data.gcrma)

layer <- celllayer <- gl(4, 24,length = length(data.raw), labels = c('EPI', 'COR', 'END', 'STE'))
group <- gl(8, 3, length = length(data.raw), labels = c("MS_1hr", "MS_48hr", 
                                                        "NaCl_1hr", "NaCl_3hr", "NaCl_8hr", "NaCl_20hr", "NaCl_32hr", "NaCl_48hr"))



#looking at the sample annotation
ph = data.raw@phenoData 
pData(data.raw)
pData(data.raw)$layer <- celllayer <- gl(4, 24,length = length(data.raw), labels = c('EPI', 'COR', 'END', 'STE'))


# Perfect-match probes
pm.data <- pm(data.raw)
head(pm.data)
# Mis-match probes
mm.data <- mm(data.raw)
head(mm.data)




###Quality test
##create microarray pictures
#print the raw intensities of each array
data.len = length(data.raw)
for (i in 1:data.len)
{
  name = paste("image",i,".jpg",sep="")
  jpeg(name)
  image(data.raw[,i],main=ph@data$sample[i])
  dev.off()
}

#create histograms of microarray data
for (i in 1:data.len)
{
  name = paste("histogram",i,".jpg",sep="")
  jpeg(name)
  hist(data.raw[,i],lwd=2,which='pm',ylab='Density',xlab='Log2 intensities',main=ph@data$sample[i])
  dev.off()
}
name = "hisplot.jpg"
jpeg(name)
cols <- rainbow(data.len * 1.2)
hist(data.raw, lty = 1:3, col = cols)
dev.off()

#create boxplots of microarray data
name = "boxplot.jpg"
jpeg(name)
boxplot(data.raw,which='pm',col='red',names=ph@data$sample) 
dev.off()

#create MA plot of microarray data
for (i in 1:data.len)
{
  name = paste("MAplot",i,".jpg",sep="")
  jpeg(name)
  MAplot(data.raw,which=i)
  dev.off()
}
#create MA plot of normalized microarray data
for (i in 1:data.len)
{
  name = paste("MAplot-norm",i,".jpg",sep="")
  jpeg(name)
  MAplot(data.gcrma,which=i)
  dev.off()
}

##Calculate quality measures
data.qc = qc(data.raw)
plot(data.qc)

## transfer to expression set to get expression matix
expset <- exprs(data.gcrma)

#create boxplots of microarray data
name = "boxplot-norm-gcrma.jpg"
jpeg(name)
boxplot(expset,which='pm',col='red',names=ph@data$sample) 
dev.off()

cor <- cor(expset,method = "pearson")
cor <- cbind(rownames(cor), cor)
write.table(cor, file = "correlation of biological replicates.txt", row.names = F, col.names = T, sep = "\t")


expset.pca <- scale(expset,center = T, scale = T)
pca <- prcomp(t(expset.pca))
summary(pca)
pca.t <- cbind(rownames(pca$x),pca$x)
write.table(pca.t, file = "PCA result in QC.txt", sep = "\t", col.names = T, row.names = F)

group2<-data.frame(group) 
layer2<-data.frame(layer)
pca_result<-as.data.frame(pca$x) 
pca_result<-cbind(pca_result,group2, layer2) 

percentage<-round(pca$sdev / sum(pca$sdev) * 100,2)
percentage<-paste(colnames(pca_result),"(", paste(as.character(percentage), "%", ")", sep=""))


p<-ggplot(pca_result,aes(x=pca_result[,1],y=pca_result[,2],color=group,shape = layer))+
  geom_point(size=3)#+
  #geom_text_repel(aes(pca_result[,1], pca_result[,2], label = rownames(ph@data)))
p<-p+theme(legend.title =element_blank()) +
  xlab(percentage[1]) +
  ylab(percentage[2]) 
p 


## Dendrogram
# define nodePar
nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), 
                cex = 0.7, col = "blue")
# cluster
hc=hclust(dist(t(expset)))
hc = hclust(as.dist(1-cor(expset)))
par(mar=c(5,5,5,10)) 

# plot
plot(as.dendrogram(hc), nodePar = nodePar,  horiz = TRUE)
