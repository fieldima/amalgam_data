
"R get_expression_matrix.r tree.nwk transcriptome_prefix"

library(ape)

#tree_file = "/Users/kf/Dropbox/kfdata/02_Data/04_Convergence_Duplication/20161114_OTU_reduction_for_6-organ_dataset/id814631/proteintree.814631.geneid.Nmasked.maxalign.degap.raxml.notung.treefix.phyml.dated.nwk"
#tc_path_prefix = "/Users/kf/Dropbox/kfdata/02_Data/04_Convergence_Duplication/20160713_Multisp_transcriptome_comparison_on_R/12-sp_6-tissues_tpm_iDEGES-N_BP-N_SVA-Y_log2-Y_pearson/mean_"


args = commandArgs(trailingOnly=TRUE)
tree_file = args[1]
tc_dir = args[2]

tree = read.tree(tree_file)

gene_ids = c()
for (gene in tree$tip.label) {
    pos_underbar = gregexpr("_", gene)[[1]]
    gene_id = substring(gene, pos_underbar[length(pos_underbar)]+1, nchar(gene))
    gene_ids = c(gene_ids, gene_id)
}
#print(paste(gene_ids, collapse=" "))

files = list.files(tc_dir)
tc_files = files
tc_table = NULL
for (tc_file in tc_files) {
    tc_sp = read.table(paste0(tc_dir, tc_file), header=TRUE, stringsAsFactors=FALSE)
    tc_sp_subset = tc_sp[rownames(tc_sp) %in% gene_ids, ]
    tc_table = rbind(tc_table, tc_sp_subset)
}
new_rownames = c()
for (gene_id in rownames(tc_table)) {
    new_rownames = c(new_rownames, tree$tip.label[grep(gene_id, tree$tip.label)])
}
rownames(tc_table) = new_rownames
tc_table = tc_table[tree$tip.label,]

if (nrow(tc_table)==length(tree$tip.label)) {
    print("Gene IDs in the tree and transcriptome are consistent.")
} else {
    print("Gene IDs in the tree and transcriptome are NOT consistent.")    
}

write.table(tc_table, file="expression_matrix.tsv", sep="\t", quote=FALSE, row.names=TRUE)
