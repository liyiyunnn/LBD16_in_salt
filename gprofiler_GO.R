### salt specific 20220517
library(gprofiler2)
setwd("/Users/liyiyun/Desktop/lbd16-1/20220517_RNAseq_DEA_1timepoint_2treatment_2genotype/gprofiler")

salt <- read.csv("/Users/liyiyun/Library/CloudStorage/OneDrive-WageningenUniversity&Research/VICI-ENW projects/Data analysis/YanxiaRNAseq2020/20220517_RNAseq_DEA_1timepoint_2treatment_2genotype/Genotype_comparison/Genotypelbd16.1.ConditionNaCl_130mM_6hrs/Genotypelbd16.1.ConditionNaCl_130mM_6hrs_signDEGs_DESeq2.csv")
con <- read.csv("/Users/liyiyun/Library/CloudStorage/OneDrive-WageningenUniversity&Research/VICI-ENW projects/Data analysis/YanxiaRNAseq2020/20220517_RNAseq_DEA_1timepoint_2treatment_2genotype/Genotype_comparison/Genotype_lbd16.1_vs_Col.0/Genotype_lbd16.1_vs_Col.0_signDEGs_DESeq2.csv")

salt.spec <- setdiff(salt$GeneID,con$GeneID)


rownames(salt) <- salt$GeneID
salt.spec.exp <- salt[salt.spec,]
write.csv(salt.spec.exp,"salt-specific_DEGs_6hr.csv")

up.deg <- salt.spec.exp[salt.spec.exp$log2FoldChange>0,]
write.csv(up.deg, "salt-specific_upregulated_DEGs_6hr.csv")

down.deg <- salt.spec.exp[salt.spec.exp$log2FoldChange<0,]
write.csv(down.deg, "salt-specific_downregulated_DEGs_6hr.csv")

go.salt.spec=gost(query = salt.spec, 
                    organism = "athaliana", ordered_query = FALSE, 
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = TRUE, 
                    user_threshold = 0.05, correction_method = "fdr", 
                    domain_scope = "annotated", sources = "GO")
head(go.salt.spec$result)
res <- apply(go.salt.spec$result,2,as.character)
write.csv(res,"saltspecific_6h_go.csv")
gostplot(go.salt.spec, capped =F, interactive = T) # interactive view
p<-gostplot(go.salt.spec, capped =F, interactive = F)
pp<- publish_gostplot(p, highlight_terms = c("GO:0071554","GO:0045229","GO:0071555","GO:0005618","GO:0071944","GO:0030312","GO:0016762" ), 
                      width = NA, height = NA, filename = NULL)
ppp<- pp+theme(text = element_text(size = 50)) 
ggsave( "saltspecific_6h_go_font.tiff", dpi=400) 

go.salt.spec.up=gost(query = up.deg$GeneID, 
                  organism = "athaliana", ordered_query = FALSE, 
                  multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                  measure_underrepresentation = FALSE, evcodes = TRUE, 
                  user_threshold = 0.05, correction_method = "fdr", 
                  domain_scope = "annotated", sources = "GO")
head(go.salt.spec.up$result)
res <- apply(go.salt.spec.up$result,2,as.character)
write.csv(res,"saltspecific_6h_go_upregulated.csv")
gostplot(go.salt.spec.up, capped =F, interactive = F)
p<-gostplot(go.salt.spec.up, capped =F, interactive = F)
pp<- publish_gostplot(p, highlight_terms = c("GO:0019758","GO:0019761","GO:0016144","GO:1901659","GO:0044272","GO:0009725","GO:0009507","GO:0043531" ), 
                      width = NA, height = NA, filename = NULL )
ggsave("saltspecific_6h_go_upregulated_withtable.tiff", dpi=300) 

go.salt.spec.down=gost(query = down.deg$GeneID, 
                     organism = "athaliana", ordered_query = FALSE, 
                     multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                     measure_underrepresentation = FALSE, evcodes = TRUE, 
                     user_threshold = 0.05, correction_method = "fdr", 
                     domain_scope = "annotated", sources = "GO")
head(go.salt.spec.down$result)
res <- apply(go.salt.spec.down$result,2,as.character)
write.csv(res,"saltspecific_6h_go_downregulated.csv")
gostplot(go.salt.spec.down, capped =F, interactive = T)
p<-gostplot(go.salt.spec.down, capped =F, interactive = F)
pp<- publish_gostplot(p, highlight_terms = c("GO:0071554","GO:0071555","GO:0045229","GO:0010054","GO:0010053","GO:0005618","GO:0016762" ), 
                      width = NA, height = NA, filename = NULL )
ggsave("saltspecific_6h_go_downregulated_withtable.tiff", dpi=300) 

