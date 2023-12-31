# Author: Yiyun Li
# Affiliation: Wageningen University & Research

# Calculate Freq. of predicted interaction in 100 runs
argv <- commandArgs(T)

if (length(argv) != 3){
  cat("Please enter the input and output file!")
  quit('no')
}

path=argv[1]
output=argv[2]

library(tidyr)

for (i in 1:100) {
  path = paste("/network-grn/run_LBD",i,".txt",sep="")
  var = paste("run",i,sep = "")
  run = as.data.frame(read.table(path, sep = "\t",header = T)[1:20,])
  run = unite(run, "Interaction", regulator, target, sep = " - ", remove = T)
  assign(var, run)
}

run.all <- rbind(run1, run2, run3, run4, run5, run6, run7, run8, run9, run10,
                 run11, run12, run13, run14, run15, run16, run17, run18, run19, run20,
                 run21, run22, run23, run24, run25, run26, run27, run28, run29, run30,
                 run31, run32, run33, run34, run35, run36, run37, run38, run39, run40,
                 run41, run42, run43, run44, run45, run46, run47, run48, run49, run50,
                 run51, run52, run53, run54, run55, run56, run57, run58, run59, run60,
                 run61, run62, run63, run64, run65, run66, run67, run68, run69, run70,
                 run71, run72, run73, run74, run75, run76, run77, run78, run79, run80,
                 run81, run82, run83, run84, run85, run86, run87, run88, run89, run90,
                 run91, run92, run93, run94, run95, run96, run97, run98, run99, run100)
re <- as.data.frame(table(run.all))
re <- re[order(-re[,2]),]
colnames(re) <- c("Interaction","Times")
re <- cbind(re, "Freq" = paste(re$Times/100*100,"%",sep = ''))

write.table(re, "/network-grn/LBD_result.txt", col.names = T, row.names = F, sep = "\t", quote = F)
