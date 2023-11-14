# Author: Yiyun Li
# Affiliation: Wageningen University & Research


library(clusterProfiler)
argv <- commandArgs(T)

if (length(argv) != 2){
  cat("Please enter the input and output file!")
  quit('no')
}

input=argv[1]
output=argv[2]

network<-read.table(input, header = F)
id.regulator <- bitr(geneID = network[,1], fromType = "TAIR", toType = "SYMBOL" , OrgDb = "org.At.tair.db", drop = F)
id.regulator <- id.regulator[!duplicated(id.regulator$TAIR), ]
id.target <- bitr(geneID = network[,2], fromType = "TAIR", toType = "SYMBOL" , OrgDb = "org.At.tair.db", drop = F)
id.target <- id.target[!duplicated(id.target$TAIR), ]

for (i in c(1:length(network[,1]))){
  id <- paste(id.regulator[match(network[i,1],id.regulator$TAIR),2],network[i,1],sep= "_")
  network[i,1] <- id
}

for (i in c(1:length(network[,2]))){
  id <- paste(id.target[match(network[i,2],id.target$TAIR),2],network[i,2],sep= "_")
  network[i,2] <- id
}

write.table(network,output,row.names = F,sep = "\t", quote = F,col.names = F)

cat("\n\nTask Finish!\n\n")


