data <- read.table("/Users/liyiyun/Downloads/PlantPAN3.0_LBD16_promoter_analysis_15022022.txt",header = T,sep = "\t")
gene <- data$TF.ID.or.Motif.name
re <- NULL
for (i in c(1:length(gene))) {
  a <- gene[i]
  b <- strsplit(gene[i],";")
  re <- append(re, b)
}
re <- unlist(re)
list <- unique(re)

salt.re <- read.table("/Users/liyiyun/Library/CloudStorage/OneDrive-WageningenUniversity&Research/VICI-ENW projects/Data analysis/GRN/LBD_100run/Salt_specific_network_ste_layer.txt")
salt.re <- sub(".*_","",salt.re$V1)

intersect(list,salt.re)


union <- read.csv("/Users/liyiyun/Library/CloudStorage/OneDrive-WageningenUniversity&Research/VICI-ENW projects/Data analysis/GRN/LBD_100run/Salt_specific_network_all_layer_and_ste_union.csv")
union <- sub(".*_","",union$regulator)
intersect(union,list)



salt <- sub(" - ASL18_AT2G42430","",network.salt.ste$Interaction)
salt <- sub(".*_","",salt)
intersect(salt,list)


con <-

enn1<-venn.diagram(list(Control = network.con.all, Salt = network.salt.all, PlantPAN3.0 = list),fill=c("FEFEBF","D3EFF8","F2F5FA"),"/Users/liyiyun/Desktop/network-grn/LBD_100run/Venn_network_aste_PlantPAN.tiff",
                     +                     width = 800,height = 800,res = 200, scaled =T, cat.dist = 0.05, margin = 0.1)
