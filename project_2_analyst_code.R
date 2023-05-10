# Project 2 --> Analyst role --> Identifying Differentially Expressed Genes Associated with Myocyte Differentiation
cuff_df <- read.delim('/projectnb/bf528/users/group_5/project_2/programmer/cuffdiff_out/gene_exp.diff')
cuff_df_1 <- cuff_df[order(cuff_df$q_value), ] #ordering based on smallest q_values on the top

#Produce a table of the top ten differentially expressed genes, with their names, 
#FPKM values, log fold change, p-value, and q-value in your report.
write.table(cuff_df_1, file = "cuff_df_1.txt", row.names = FALSE, col.names = FALSE)

#Produce a histogram of the log2.foldchange column for all genes with the hist function.
hist(cuff_df_1$log2.fold_change., main = "Distribution of log2fc across all the genes", 
     xlab = "log2foldchange", col = "indianred", breaks = 30, xlim = c(-10,10), ylim = c(0,20000))

#Create a new data frame that contains only the genes where the last column, named significant, is equal to yes. 
sig_yes_df <- subset(cuff_df_1, significant == "yes")

#Subsetting the top 1000 significant genes
thousand_sig_genes_subset <- sig_yes_df[1:1000, ]
write.table(thousand_sig_genes_subset, file = "thousand_sig_genes_subset.txt", row.names = FALSE, col.names = FALSE)

#Create a second histogram of the log2 fold change values only for significant genes.
hist(sig_yes_df$log2.fold_change., main = "Distribution of log2fc across significant genes", 
     xlab = "log2foldchange", col = "darkorchid1", breaks = 30, xlim = c(-10,10))

#Absolute value of log2fold change should be greater than 1 for significant genes 
up_genes <- subset(sig_yes_df, log2.fold_change. > 0)
down_genes <- subset(sig_yes_df, log2.fold_change. < 0)

#write the up and down regulated genes as csv
up_genenames_only <- up_genes$gene
down_genenames_only <- down_genes$gene
write.csv(up_genenames_only, file = "upregulated_genes.csv", row.names = FALSE)
write.csv(down_genenames_only, file = "downregulated_genes.csv", row.names = FALSE)

