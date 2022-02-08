
library(ape)
library(phytools)

if (length(commandArgs(trailingOnly=TRUE))==1) {
    mode = "debug"    
} else {
    mode = "batch"
}
print(paste('mode:', mode))

if (mode=="debug") {
    og = 'OG0000526'
    dir_tmp = '/Users/kf/Dropbox/kfdata/02_Data/my_db/Ensembl/release-91/orthogroup/tmp/'
    dirs = list.files(dir_tmp)
    dir_work = paste0(dir_tmp, dirs[grep(og, dirs)], '/')
    files = list.files(dir_work)
    tree_file = paste0(dir_work, og, '.dated.nwk')
    gff3_dir = "/Users/kf/Dropbox/kfdata/02_Data/my_db/Ensembl/release-91/gff3/"
    setwd(dir_work)
} else if (mode=="batch") {
    args = commandArgs(trailingOnly=TRUE)
    tree_file = args[1]
    gff3_dir = args[2]
    work_dir = getwd()
}

if (! grepl('/$', gff3_dir)) {
    gff3_dir = paste0(gff3_dir, '/')
}

print(paste('tree_file =', tree_file))
print(paste('gff3_dir =', gff3_dir))

my_time = proc.time()

tree = read.newick(tree_file)
options(stringsAsFactors=FALSE)
genus = as.character(data.frame(strsplit(tree$tip.label, '_'))[1,])
species = as.character(data.frame(strsplit(tree$tip.label, '_'))[2,])
spp = unique(paste0(genus, '_', species))
cat('# species in the gene tree =', length(spp), '\n')

gff3_files = list.files(gff3_dir)
gff3_files = gff3_files[grep('.gz$', gff3_files)]
gff3_files = gff3_files[sub('\\..*', '', gff3_files) %in% spp]
cat('# gff3 files to be processed =', length(gff3_files), '\n')
trait = data.frame(matrix(NA, length(tree$tip.label), 3))
colnames(trait) = c('leaf', 'num_intron', 'chromosome')
i=1
records = list()
unique_seqids = list()
for (gff3 in gff3_files) {
    cat('working on:', gff3, '\n')
    sci_name = sub('\\..*', '', gff3)
    df = read.table(paste0(gff3_dir, gff3), sep='\t', quote='', header=FALSE, comment.char='#', stringsAsFactors=FALSE)
    colnames(df) = c('seqid','source','type','start','end','score','strand','phase','attributes')
    unique_seqids[[gff3]] = unique(df$seqid)
    leaf_names = tree$tip.label
    leaf_names = leaf_names[startsWith(leaf_names, sci_name)]
    gene_ids = sub(".*_", "", leaf_names)
    search_phrase = paste(paste0('Parent=gene:',gene_ids, ';'), collapse='|')
    cat('search_phrase:', search_phrase, '\n')
    attr_text_all = df$attributes[grep(search_phrase, df$attributes)]
    cat(sci_name, ': # gene in the gene tree =', length(leaf_names), '\n')
    for (leaf_name in leaf_names) {
        if (! leaf_name %in% names(records)) {
            gene_id = sub(".*_", "", leaf_name)
            attr_text = attr_text_all[grep(gene_id, attr_text_all)]
            if (length(attr_text)>0) {
                records[[leaf_name]] = list()
                cat('processing:', leaf_name, '\n')
                transcript_ids = attr_text
                transcript_ids = sub(".*ID=transcript:", "", transcript_ids)
                transcript_ids = sub(";.*", "", transcript_ids)
                max_intron_num = 0
                for (transcript_id in transcript_ids) {
                    df2 = df[grep(paste0('Parent=transcript:',transcript_id, ';'), df$attributes), ]
                    df2 = df2[(df2$type=='exon'), ]
                    max_intron_num = max(max_intron_num, nrow(df2)-1)
                    records[[leaf_name]][[transcript_id]] = df2
                }
                trait[i,'leaf'] = leaf_name
                trait[i,'num_intron'] = max_intron_num
                trait[i,'chromosome'] = df2$seqid[1]
                i=i+1
            }
        }
    }
}
leaf_order = tree$tip.label
if (mode=='debug') {
    trait = trait[!apply(trait, 1, function(x){all(is.na(x))}),]
    leaf_order = leaf_order[(leaf_order %in% trait$leaf)]
}
trait = na.omit(trait)
rownames(trait) = trait$leaf
trait = trait[leaf_order,]
rownames(trait) = trait$leaf
trait$leaf = NULL
cat("Time elapsed for the trait value extraction from gff3 files\n")
print(proc.time()-my_time)

is_output_num_correct = (length(tree$tip.label)==nrow(trait))
cat('# leaves =', length(tree$tip.label), '# trait row =', nrow(trait), '\n')
if (any(mode=='debug',is_output_num_correct)) {
    trait_intron = trait
    trait_intron$chromosome = NULL
    trait_chromosome = trait
    trait_chromosome$num_intron = NULL
    write.table(trait_intron, file="intron_matrix.tsv", sep="\t", quote=FALSE, row.names=TRUE)
    write.table(trait_chromosome, file="chromosome_matrix.tsv", sep="\t", quote=FALSE, row.names=TRUE)
    cat('Done!\n')
} else {
    cat('Output number is inconsistent. Leaf names not found:\n')
    print(tree$tip.label[(!tree$tip.label %in% rownames(trait))])
    cat('Exiting without generating output files.\n')
}


if (mode=='debug') {
    trait_intron
}

if (mode=='debug') {
    trait_chromosome
}
