con <- t(read.table("/Users/liyiyun/Library/CloudStorage/OneDrive-WageningenUniversity&Research/VICI-ENW projects/Data analysis/GRN/input_94/expset_con_alllayer.txt",sep = "\t"))
salt <-t(read.table("/Users/liyiyun/Library/CloudStorage/OneDrive-WageningenUniversity&Research/VICI-ENW projects/Data analysis/GRN/input_94/expset_salt_alllayer.txt",sep = "\t"))
rownames(con) == rownames(salt)
input = cbind(con,salt)
write.table(input, "/Users/liyiyun/Library/CloudStorage/OneDrive-WageningenUniversity&Research/VICI-ENW projects/Data analysis/GRN/input_94/input_exp.txt",sep = "\t",quote = F)

exp <- read.table("/Users/liyiyun/Desktop/network-grn/expset_salt.txt",sep = "\t")
exp.in <- exp[rownames(con),]
write.table(exp.in, "/Users/liyiyun/Library/CloudStorage/OneDrive-WageningenUniversity&Research/VICI-ENW projects/Data analysis/GRN/input_94/input_94_expression.txt",sep = "\t",quote = F)
