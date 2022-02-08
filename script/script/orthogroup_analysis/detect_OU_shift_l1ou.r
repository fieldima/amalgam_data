
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
    install.packages(c("igraph", "phylolm", "lars", "grplasso", "magic", "genlasso", "Rcpp"))
    install.packages("knitr")
    remove.packages("l1ou")
    devtools::install_github("khabbazian/l1ou", dep=FALSE)
    #devtools::install_git("git://github.com/kfuku52/rkftools.git", branch = "master", dep=FALSE)
    
    # cd
    # git clone git://github.com/kfuku52/rkftools.git
    # R
    # devtools::install_local(path='/lustre1/home/kfuku/rkftools', reload=TRUE, quick=FALSE, local=TRUE, dep=FALSE)
    # q()
}

library(l1ou)
library(rkftools)
options(stringsAsFactors=FALSE)

# one bootstrap replicate needs approx. 1/2 time of the original fit estimation at nCores=1.
# OUfixedRoot is preferrable, https://github.com/khabbazian/l1ou/pull/5
# fixed.alpha =FALSE might take too long time.

if (mode=="debug") {
    #og = 'OG0002332' # PGK
    og = 'OG0001999'
    #dir_tmp = '/Users/kef74yk/Dropbox_w/db/Ensembl/release-91/orthogroup/tmp/'
    dir_tmp = '/Users/kef74yk/Dropbox_w/db/gastrodia/'
    dirs = list.files(dir_tmp)
    dir_work = paste0(dir_tmp, dirs[grep(og, dirs)], '/')
    files = list.files(dir_work)
    setwd(dir_work)
    args = c()
    #args = c(args, paste0('--tree_file=', dir_work, og, '.dated.nwk'))
    args = c(args, paste0('--tree_file=', dir_work, og, '.dated.pruned.nwk'))
    args = c(args, paste0('--trait_file=', dir_work, og, '.expression.tsv'))
    args = c(args, paste0('--criterion=', 'AICc'))
    args = c(args, paste0('--nslots=', 8))
    args = c(args, paste0('--nbootstrap=', 0))
    #args = c(args, paste0('--fit_ind_file=', dir_work, og, '.l1ou.RData'))
    args = c(args, paste0('--fit_ind_file=', ''))
    args = c(args, paste0('--fit_conv_file=', ''))
    args = c(args, paste0('--alpha_upper=', 'PhylogeneticEM'))
    args = c(args, paste0('--detect_convergence=', 0))
    args = c(args, paste0('--similarity_method=', 'pearson'))
    args = c(args, paste0('--similarity_threshold=', 0.99)) # 0 for no check
    args = c(args, paste0('--ceil_negative=', 0))
    args = c(args, paste0('--max_nshift=', 0)) # 0 for automatic
} else if (mode=="batch") {
    args = commandArgs(trailingOnly=TRUE)
}

cat('arguments:\n')
args = rkftools::get_parsed_args(args, print=TRUE)

# this function was copied from the package PhylogeneticEM
find_grid_alpha <- function(phy, alpha = NULL,
                            nbr_alpha = 10,
                            factor_up_alpha = 2,
                            factor_down_alpha = 3,
                            quantile_low_distance = 0.0001,
                            log_transform = TRUE, ...){
  if (!is.null(alpha)) return(alpha)
  dtips <- cophenetic(phy)
  d_min <- quantile(dtips[dtips > 0], quantile_low_distance)
  h_tree <- node.depth.edgelength(phy)[1]
  alpha_min <- 1 / (factor_down_alpha * h_tree)
  alpha_max <- factor_up_alpha / (2 * d_min)
  alpha_max_machine <- log(.Machine$double.xmax^0.98)/(2*h_tree)
  if (alpha_max > alpha_max_machine){
    warning("The chosen alpha_max was above the machine precision. Taking alpha_max as the largest possible on this machine.")
    alpha_max <- alpha_max_machine
  }
  if (log_transform){
    return(c(0, exp(seq(log(alpha_min), log(alpha_max), length.out = nbr_alpha))))
  } else {
    return(c(0, seq(alpha_min, alpha_max, length.out = nbr_alpha)))
  }
}

my_time = proc.time()
tree = read.tree(args[['tree_file']])
tree_original = tree
trait_table = read.table(args[['trait_file']], header=TRUE, row.names=1, sep="\t")

if ((args[['similarity_method']]=='complementarity')|(args[['ceil_negative']])) {
    if (args[['similarity_method']]=='complementarity') {
        cat('similarity_method is set to complementary: forced to ceil negative values.\n')
    }
    cat('Minimum trait value =', min(trait_table), ', Negative values are ceiled to 0.\n')
    if (min(trait_table) < 0) {
        trait_table[trait_table < 0] = 0
    }
}

tree = rkftools::pad_short_edges(tree, threshold=1e-3, external_only=FALSE)
tree = rkftools::force_ultrametric(tree, stop_if_larger_change=0.01)

original_trait_table = trait_table
original_traits = colnames(original_trait_table)

out = remove_invariant_traits(trait_table, small_dif=1e-3)
trait_table = out[['trait_table']]
removed_traits = out[['removed_traits']]

if (args[['similarity_threshold']]) {
    cat('Start cluster collapsing by similarity thresholding ( cutoff =', args[['similarity_threshold']], ')\n')
    high_comple_clades = rkftools::get_high_similarity_clades(tree, trait_table, args[['similarity_method']], 
                            args[['similarity_threshold']], verbose=FALSE, num_test=0)
    out = rkftools::collapse_clades(tree, trait_table, high_comple_clades)
    tree = out[['tree']]
    trait_table = out[['trait']]
} else {
    cat('Skipped cluster collapsing.', '\n')
}

trait_matrix = as.matrix(trait_table)
if (mode=="debug") {
    trait_matrix
}

cat("Time elapsed for data preparation:\n")
print(proc.time()-my_time)

my_time = proc.time()
fit_ind=NULL
fit_conv=NULL
adj_data=NULL
adj_data = l1ou::adjust_data(tree=tree, Y=trait_matrix, normalize=FALSE)

alpha.lower = NA
if (args[['alpha_upper']]=='Inf') {
    alpha.upper = Inf
} else if (args[['alpha_upper']]=='l1ou') {
    alpha.upper = l1ou:::alpha_upper_bound(adj_data$tree)
} else if (args[['alpha_upper']]=='PhylogeneticEM') {
    alphas = find_grid_alpha(phy=tree, nbr_alpha=10)
    alpha.upper = alphas[11]
    # alpha.lower = alphas[1] # l1ou fails if alpha.lower is set to 0.
} else {
    alpha.upper = as.numeric(args[['alpha_upper']])
}
cat('default alpha upper bound in l1ou =',  l1ou:::alpha_upper_bound(adj_data$tree), '\n')
cat('the alpha upper bound actually used =', alpha.upper, '\n')

args[['max_nshift']] = as.integer(args[['max_nshift']])
num_leaf = length(tree$tip.label)
if (args[['max_nshift']]==0) {
    max_nshift = as.integer(max(min(100, num_leaf/2), floor(num_leaf^0.5)))
} else {
    max_nshift = args[['max_nshift']]
}
cat('maximum number of shifts =', max_nshift, '\n')

ind_flag = FALSE
if ((!is.null(args[['fit_ind_file']]))&(!is.na(args[['fit_ind_file']]))) {
    if (file.exists(args[['fit_ind_file']])) {
        if (file.info(args[['fit_ind_file']])$size) {
            load(args[['fit_ind_file']])
            if(class(fit_ind)!="try-error") {
                cat(args[['fit_ind_file']], 'was found. Loading the shift configurations.\n')
                ind_flag = TRUE
            } else {
                cat(args[['fit_ind_file']], 'is malformed. Printing...\n')
                print(fit_ind)
                rm(fit_ind)
            }
        } else {
            cat(args[['fit_ind_file']], 'has size zero.\n')
        }
    }
}
if (!ind_flag) {
    cat('fit_ind.RData was not found. Starting shift detection.\n')
    fit_ind = try(l1ou::estimate_shift_configuration(
        tree=adj_data$tree, 
        Y=adj_data$Y, 
        max.nShifts=max_nshift,
        criterion=args[['criterion']], 
        root.model="OUfixedRoot", 
        nCores=args[['nslots']], 
        rescale=FALSE, 
        alpha.upper=alpha.upper, 
        alpha.lower=alpha.lower,
        quietly=FALSE))
    save(fit_ind, file='fit_ind.RData')
    cat("Time elapsed for shift detection:\n")
}
print(proc.time()-my_time)

conv_flag = FALSE
if ((!is.null(args[['fit_conv_file']]))&(!is.na(args[['fit_conv_file']]))) {
    if (file.exists(args[['fit_conv_file']])) {
        if (file.info(args[['fit_conv_file']])$size) {
            load(args[['fit_conv_file']])
            if(class(fit_conv)!="try-error") {
                cat(args[['fit_conv_file']], 'was found. Loading the shift configurations.\n')
                conv_flag = TRUE
            } else {
                cat(args[['fit_conv_file']], 'is malformed. Printing...\n')
                print(fit_conv)
                rm(fit_conv)
            }
        } else {
            cat(args[['fit_conv_file']], 'has size zero.\n')
        }
    }
}
if ((!conv_flag)&(args[['detect_convergence']])) {
    cat('Started convergence analysis.\n')
    if (class(fit_ind)=="l1ou") {
        if (fit_ind$nShifts < 2) {
            fit_conv = fit_ind
        } else {
            fit_conv = try(l1ou::estimate_convergent_regimes(
                fit_ind, 
                criterion=args[['criterion']], 
                method="backward", 
                fixed.alpha=TRUE, 
                nCores=args[['nslots']]))
        }
    } else {
        warning("estimate_shift_configuration() failed.")
        warning(attr(fit_ind, "condition"))
    }
    save(fit_conv, file='fit_conv.RData')
    cat("Time elapsed for shift detection + convergence detection:\n")
    print(fit_conv)
    print(proc.time()-my_time)
} else {
    cat('Skipped convergence analysis.\n')
    fit_conv = fit_ind
}


if (class(fit_conv)=="l1ou") {
    if (args[['nbootstrap']] > 0) {
        cat("Start bootstrapping,", args[['nbootstrap']], "replicates,", args[['nslots']], "cores.\n")
        bootstrap_result = l1ou::l1ou_bootstrap_support(
            fit_conv, 
            nItrs=args[['nbootstrap']], 
            multicore=(args[['nslots']]>1), 
            nCores=args[['nslots']])
        bootstrap_table = rkftools::get_bootstrap_table(fit_conv, bootstrap_result, mode='l1ou')
        write.table(bootstrap_table, file="l1ou_bootstrap.tsv", sep="\t", quote=FALSE, row.names=FALSE)
        cat("Time elapsed for shift detection + convergence detection + bootstrapping:\n")
        print(proc.time()-my_time)
    } 
    cat('Start preparing output tables.\n')
    tree_table = rkftools::get_tree_table(fit_conv, mode='l1ou')
    regime_table = rkftools::get_regime_table(fit_conv, mode='l1ou')
    leaf_table = rkftools::get_leaf_table(fit_conv, mode='l1ou')
    for (mt in removed_traits) {
        regime_table[mt] = NA
        leaf_table[mt] = 0 
    }
    regime_table = regime_table[,c(colnames(regime_table)[!colnames(regime_table) %in% original_traits], original_traits)]
    leaf_table = leaf_table[,c(colnames(leaf_table)[!colnames(leaf_table) %in% original_traits], original_traits)]
    
    if (args[['similarity_threshold']]) {
        cat('Starting the recovery of collapsed nodes.\n')
        node_num_mapping = rkftools::map_node_num(tree_original, tree, out[['collapse_leaf_names']], verbose=FALSE)
        tree_table = rkftools::tree_table_collapse2original(tree_table, tree_original)
        regime_table = rkftools::regime_table_collapse2original(regime_table, tree_original, tree, node_num_mapping)
        leaf_table = rkftools::leaf_table_collapse2original(leaf_table, tree_original, tree, node_num_mapping)
    }
    
    write.table(tree_table, file="l1ou_tree.tsv", sep="\t", quote=FALSE, row.names=FALSE) 
    write.table(regime_table, file="l1ou_regime.tsv", sep="\t", quote=FALSE, row.names=FALSE) 
    write.table(leaf_table, file="l1ou_leaf.tsv", sep="\t", quote=FALSE, row.names=FALSE)
    pdf("l1ou_plot.pdf", height=length(tree$tip.label)/5+1, width=ncol(adj_data$Y)/2+5)
    par(mar=c(2,2,2,2))
    try(plot(fit_conv, edge.shift.ann=FALSE))
    graphics.off()
} else {
    warning("estimate_convergent_regimes() failed.")
    warning(attr(fit_conv, "condition"))
}

cat('l1ou completed!\n')

if (mode=='debug') {
    edge_labels = 1:(length(adj_data$tree$tip.label)+length(adj_data$tree$node.label))
    for (branch_index in fit_conv$shift.configuration) {
        shift_node = adj_data$tree$edge[branch_index,2]
        leaf_names = get_tip_labels(adj_data$tree, shift_node)
        leaf_nums = get_node_num_by_name(adj_data$tree, leaf_names)
        mu = fit_conv$mu[leaf_nums,]
        cat(shift_node, 'N leaf_name', length(leaf_names), 'N leaf_num', length(leaf_nums), '\n')
        cat('leaf_names', leaf_names, '\n')
        cat('nrow mu', nrow(mu), '\n')
        print(mu)
        if (!is.null(dim(mu))) {
            mu = mu[1,]
        }
        mu = round(mu, digits=2)
        mu = paste0('\n', paste(mu[1:3], collapse=','), '\n', paste(mu[4:6], collapse=','))
        print(mu)
        #edge_labels[branch_index] = paste0(edge_labels[shift_index], '\n', mu)
    }
    pdf("l1ou_plot.debug.pdf", height=length(adj_data$tree$tip.label)/5+1, width=ncol(adj_data$Y)+5)
    par(mar=c(2,2,2,2))
    try(plot(fit_conv, edge.shift.ann=FALSE, edge.label=edge_labels, edge.label.ann=TRUE, plot.bar=TRUE, cex=0.3))
    graphics.off()
}


