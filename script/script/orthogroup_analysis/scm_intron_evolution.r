
mode = ifelse(length(commandArgs(trailingOnly=TRUE))==1, 'debug', 'batch')

if (mode=='debug') {
    #devtools::install_github("kfuku52/rkftools", dep=FALSE)
    #remove.packages('rkftools')
    devtools::install_local(path='/Users/kef74yk/Dropbox_w/repos/rkftools', reload=TRUE, quick=FALSE, local=TRUE, dep=FALSE)
    options(warn=-1)
}

if (FALSE) {
    library(devtools)
    library(httr)
    set_config(config(ssl_verifypeer = 0L))
    options(repos=structure(c(CRAN="http://cran.rstudio.com/")))
    remove.packages("phytools")
    install.packages(c('animation', 'ape', 'maps', 'nlme', 'phangorn', 'plotrix', 'scatterplot3d'))
    devtools::install_github("liamrevell/phytools", dependency=FALSE)
}
library(phytools)
library(parallel)

if (mode=="debug") {    
    og = 'OG0000999'
    dir_tmp = '/Users/kef74yk/Dropbox_w/db/Ensembl/release-91/orthogroup/tmp/'
    dirs = list.files(dir_tmp)
    dir_work = paste0(dir_tmp, dirs[grep(og, dirs)], '/')
    files = list.files(dir_work)
    setwd(dir_work)
    args = c()
    args = c(args, paste0('--tree_file=', dir_work, og, '.dated.nwk'))
    args = c(args, paste0('--trait_file=', dir_work, og, '.intron.tsv'))
    args = c(args, paste0('--intron_gain_rate=', '0.001'))
    args = c(args, paste0('--retrotransposition_rate=', '0.001'))
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

if (mode=='debug') {
    plot(tree)
    trait
}

# Intron
my_time = proc.time()
tree2 = tree
#tree2$edge.length = tree2$edge.length/max(nodeHeights(tree2)[,2])*scale
tree2$tip.label = paste(trait$num_intron, tree2$tip.label)
introns = trait$num_intron
names(introns) = tree2$tip.label
x = ifelse(introns>0, 'intron_present', 'intron_absent')
names(x) = names(introns)
intron_xmatrix = to.matrix(x,sort(unique(x)))
for (chr in c('intron_present', 'intron_absent')) {
    if (! chr %in% colnames(intron_xmatrix)) {
        xdf = data.frame(intron_xmatrix)
        xdf[chr] = rep(0, nrow(intron_xmatrix))
        intron_xmatrix = as.matrix(xdf)
    }
}
#cols = setNames(palette()[1:length(unique(x))],sort(unique(x)))
intron_cols = c('white', 'black')
names(intron_cols) = c('intron_absent', 'intron_present')

Q = matrix(
    c(
        -args[['intron_gain_rate']], 
        args[['intron_gain_rate']], 
        args[['retrotransposition_rate']], 
        -args[['retrotransposition_rate']]
    )
    , 2, 2, byrow=TRUE)

colnames(Q) = c('intron_absent', 'intron_present')
rownames(Q) = c('intron_absent', 'intron_present')
print(Q)


intron_mtrees = mclapply(1:args[['nslots']], function(n, tree, x, fixedQ) {make.simmap(tree, x, Q=fixedQ, nsim=ceiling(args[['nrep']]/args[['nslots']]))}, 
                         tree=tree2, x=intron_xmatrix, fixedQ=Q, mc.cores=args[['nslots']])
intron_mtrees = do.call(c, intron_mtrees)
if (!("multiSimmap" %in% class(intron_mtrees))) {
    class(intron_mtrees) = c("multiSimmap",class(intron_mtrees))
}
#intron_mtrees = make.simmap(tree2, intron_xmatrix, model='ARD', nsim=args[['nrep']], pi='equal', Q=Q) #Q='empirical', use.empirical=TRUE)
intron_summary = summary(intron_mtrees, plot=FALSE)
print("Time elapsed for the stochastic character mapping: intron")
print(proc.time()-my_time)

cols = c('leaf', 'num_intron', 'intron_present', 'intron_absent')
df_out = data.frame(matrix(NA, 0, length(cols)))
colnames(df_out) = cols
node_labels = c(tree$tip.label, tree$node.label)
df_out[1:length(node_labels),'leaf'] = node_labels
df_out[1:nrow(trait),'num_intron'] = trait$num_intron
if (any(apply(intron_xmatrix, 2, sum)==nrow(intron_xmatrix))) {
    for (col in colnames(intron_xmatrix)) {
        df_out[,col] = ifelse(sum(intron_xmatrix[,col])==nrow(intron_xmatrix), 1, 0)
    }
} else if (nrow(intron_summary$tips)!=nrow(intron_xmatrix)) {
    for (col in colnames(intron_xmatrix)) {
        num_zero = sum(intron_xmatrix[,col]==0)
        num_one = sum(intron_xmatrix[,col]==1)
        anc_state = ifelse(num_one>=num_zero, 1, 0)
        df_out[1:nrow(intron_xmatrix),col] = intron_xmatrix[,col]
        df_out[is.na(df_out[col]),col] = anc_state
    }
} else {
    df_out[,colnames(intron_summary$tips)] = rbind(intron_summary$tips, intron_summary$ace)
}
write.table(df_out, file="intron_evolution_summary.tsv", sep="\t", quote=FALSE, row.names=FALSE)
if (mode=='debug') {
    df_out
}

fsize=0.35
my_plot = function() {
    plot(intron_summary, colors=intron_cols, fsize=fsize, ftype="reg")
    add.simmap.legend(colors=intron_cols, prompt=FALSE, x=0.9*par()$usr[1], y=5,fsize=fsize)
}

if (any(apply(intron_xmatrix, 2, sum)==nrow(intron_xmatrix))) {
    print("Tree plotting was skipped because there is no change in the character state")
} else {
    print("Plotting the phylogenetic tree.")
    pdf("intron_evolution_plot.pdf", width=7.2, height=length(tree$tip.label)/8)
    my_plot()
    dev.off()
    if (mode=='debug') {
        my_plot()
    }
}

cat('scm intron evolution completed!', '\n')


