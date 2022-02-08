
mode = ifelse(length(commandArgs(trailingOnly=TRUE))==1, 'debug', 'batch')

if (mode=='debug') {
    #devtools::install_github("kfuku52/rkftools", dep=FALSE)
    #remove.packages('rkftools')
    devtools::install_local(path='/Users/kf/Dropbox/kfdata/02_Data/my_projects/rkftools', reload=TRUE, quick=FALSE, local=TRUE, dep=FALSE)
    #devtools::install_git("git://github.com/kfuku52/rkftools.git", branch = "master", dep=FALSE)
}

library(rkftools)
library(phytools)
library(parallel)
options(stringsAsFactors=FALSE)

if (mode=="debug") {    
    og = 'OG0000999'
    dir_tmp = '/Users/kef74yk/Dropbox_w/db/Ensembl/release-91/orthogroup/tmp/'
    dirs = list.files(dir_tmp)
    dir_work = paste0(dir_tmp, dirs[grep(og, dirs)], '/')
    files = list.files(dir_work)
    setwd(dir_work)
    args = c()
    args = c(args, paste0('--tree_file=', dir_work, og, '.dated.nwk'))
    args = c(args, paste0('--trait_file=', dir_work, og, '.chromosome.tsv'))
    args = c(args, paste0('--chr_loc_transition_rate=', '0.001'))
    args = c(args, paste0('--omit_species_prefix=', 'Astyanax_mexicanus|Danio_rerio|Gadus_morhua|Oryzias_latipes|Oreochromis_niloticus|Xenopus_tropicalis|Gallus_gallus|Anolis_carolinensis|Ornithorhynchus_anatinus|Chinchilla_lanigera')) # Excluded because their sex chromosomes are non-homologous to mammals
    args = c(args, paste0('--omit_chromosome_prefix=', 'Un'))
    args = c(args, paste0('--max_chromosome_nchar=', '4'))
    args = c(args, paste0('--nrep=', '3'))
    args = c(args, paste0('--nslots=', '3'))
} else if (mode=="batch") {
    args = commandArgs(trailingOnly=TRUE)
}

cat('arguments:\n')
args = rkftools::get_parsed_args(args, print=TRUE)

tree = read.newick(args[['tree_file']])
trait = read.table(args[['trait_file']], header=TRUE, row.names=1, sep="\t", stringsAsFactors=FALSE)

tree = rkftools::pad_short_edges(tree, threshold=1e-3, external_only=FALSE)
tree = rkftools::force_ultrametric(tree, stop_if_larger_change=0.01)

my_time = proc.time()
tree2 = tree
x = trait$chromosome
names(x) = tree2$tip.label
char_states = c('A','X','Y','MT')
chr_cols = c('white', 'green', 'magenta', 'black')
names(chr_cols) = char_states
num_char = length(char_states)

# Prepare Q matrix
if (args[['chr_loc_transition_rate']]=='empirical') {
    fit = phytools::fitMk(tree=tree2, x=chr_xmatrix, model="ARD")
    Q = matrix(NA, length(fit$states), length(fit$states))
    Q[] = c(0, fit$rates)[fit$index.matrix+1]
    diag(Q) = -0
    diag(Q) = -rowSums(Q)
    colnames(Q) = rownames(Q) = fit$states
} else {
    args[['chr_loc_transition_rate']] = as.numeric(args[['chr_loc_transition_rate']])
    Q = matrix(rep(args[['chr_loc_transition_rate']], num_char^2), num_char, num_char, byrow=TRUE)
    colnames(Q) = char_states
    rownames(Q) = char_states
    for (i in 1:num_char) {
        Q[i,i] = -args[['chr_loc_transition_rate']]
    }
}
print(Q)

# Prepare trait input
x[!((x=='Y')|(x=='X')|(x=='MT'))] = 'A'
chr_xmatrix = to.matrix(x,sort(unique(x)))
original_cols = colnames(chr_xmatrix)
for (chr in char_states) {
    if (! chr %in% colnames(chr_xmatrix)) {
        xdf = data.frame(chr_xmatrix)
        xdf[chr] = rep(0, nrow(chr_xmatrix))
        chr_xmatrix = as.matrix(xdf)
    }
}

# Omit species prefix
if (!is.blank(args[['omit_species_prefix']])) {
    omit_species_prefix = strsplit(args[['omit_species_prefix']], '\\|')[[1]]
    cat('\nOmitted leaves:\n')
    for (pref in omit_species_prefix) {
        ind = grep(pref, rownames(chr_xmatrix))
        cat('omit_species_prefix:', pref, ': number of omitted leaves =', length(ind), '\n')
        if (length(ind)>0) {
            chr_xmatrix[ind,original_cols] = matrix(1/length(original_cols), length(ind), length(original_cols))            
        }
    }    
}

# Omit chromosome prefix
if (!is.blank(args[['omit_chromosome_prefix']])) {
    omit_chromosome_prefix = strsplit(args[['omit_chromosome_prefix']], '\\|')[[1]]
    conditions = FALSE
    for (pref in omit_chromosome_prefix) {
        conditions = ((conditions)|(grepl(pref, trait$chromosome)))    
    }
    is_omitted_already = apply(chr_xmatrix, 1, function(x){!all(x%in%c(0,1))})
    conditions = (conditions)&(!is_omitted_already)
    chr_xmatrix[conditions,original_cols] = matrix(1/length(original_cols), sum(conditions), length(original_cols))
    unmapped_leaves = rownames(chr_xmatrix)[conditions]
    unmapped_annotations = trait[conditions,]
    cat('\nOmitted leaves:\n')
    for (i in 1:length(unmapped_leaves)) {
        cat('omit_chromosome_prefix:', unmapped_leaves[i], 'on', unmapped_annotations[i], '\n')
    }
}

# Max chromosome nchar
if (!is.blank(args[['max_chromosome_nchar']])) {
    conditions = (nchar(as.character(trait$chromosome))>args[['max_chromosome_nchar']])
    is_omitted_already = apply(chr_xmatrix, 1, function(x){!all(x%in%c(0,1))})
    conditions = (conditions)&(!is_omitted_already)
    if (any(conditions)) {
        chr_xmatrix[conditions,original_cols] = matrix(1/length(original_cols), sum(conditions), length(original_cols))
        unmapped_leaves = rownames(chr_xmatrix)[conditions]
        unmapped_annotations = trait[conditions,]
        cat('\nOmitted leaves:\n')
        for (i in 1:length(unmapped_leaves)) {
            cat('max_chromosome_nchar:', unmapped_leaves[i], 'on', unmapped_annotations[i], '\n')
        }
    }
}

cat("\nTime elapsed for trait matrix preparation\n")
print(proc.time()-my_time)

my_time = proc.time()
chr_mtrees = mclapply(1:args[['nslots']], function(n, tree, x, fixedQ) {make.simmap(tree, x, Q=fixedQ, nsim=ceiling(args[['nrep']]/args[['nslots']]))}, 
                         tree=tree2, x=chr_xmatrix, fixedQ=Q, mc.cores=args[['nslots']])
chr_mtrees = do.call(c, chr_mtrees)
if (!("multiSimmap" %in% class(chr_mtrees))) {
    class(chr_mtrees) = c("multiSimmap",class(chr_mtrees))
}
chr_summary = describe.simmap(chr_mtrees, plot=FALSE)
cat("Time elapsed for the stochastic character mapping: chromosome\n")
print(proc.time()-my_time)

cols = c('leaf', 'chromosome', 'A', 'X', 'Y', 'MT')
df_out = data.frame(matrix(NA, 0, length(cols)))
colnames(df_out) = cols
node_labels = c(tree$tip.label, tree$node.label)
df_out[1:length(node_labels),'leaf'] = node_labels
df_out[1:nrow(trait),'chromosome'] = trait$chromosome
if (any(apply(chr_xmatrix, 2, sum)==nrow(chr_xmatrix))) {
    for (col in colnames(chr_xmatrix)) {
        df_out[,col] = ifelse(sum(chr_xmatrix[,col])==nrow(chr_xmatrix), 1, 0)
    }
} else if (nrow(chr_summary$tips)!=nrow(chr_xmatrix)) {
    for (col in colnames(chr_xmatrix)) {
        num_zero = sum(chr_xmatrix[,col]==0)
        num_one = sum(chr_xmatrix[,col]==1)
        anc_state = ifelse(num_one>=num_zero, 1, 0)
        df_out[1:nrow(chr_xmatrix),col] = chr_xmatrix[,col]
        df_out[is.na(df_out[col]),col] = anc_state
    }
} else {
    df_out[,colnames(chr_summary$tips)] = rbind(chr_summary$tips, chr_summary$ace)
}
write.table(df_out, file="chromosome_location_evolution_summary.tsv", sep="\t", quote=FALSE, row.names=FALSE)
if (mode=='debug') {
    df_out
}

fsize=0.35
my_plot = function() {
    plot(chr_summary, colors=chr_cols, fsize=fsize, ftype="reg")
    add.simmap.legend(colors=chr_cols, prompt=FALSE, x=0.9*par()$usr[1], y=5,fsize=fsize)
}

if (any(apply(chr_xmatrix, 2, sum)==nrow(chr_xmatrix))) {
    cat("Tree plotting was skipped because there is no change in the character state.", '\n')
} else {
    cat("Plotting the phylogenetic tree.", '\n')
    pdf("chromosome_location_evolution_plot.pdf", width=7.2, height=length(tree$tip.label)/8)
    my_plot()
    dev.off()
    if (mode=='debug') {
        my_plot()
    }
}

cat('scm chromosome location evolution completed!\n')
