
if (FALSE) {
    library(devtools)
    library(httr)
    set_config(config(ssl_verifypeer = 0L))
    options(repos=structure(c(CRAN="http://cran.rstudio.com/")))
    #remove.packages("phytools")
    #install_github("liamrevell/phytools", dep=FALSE)
    devtools::install_github("kfuku52/rkftools", dep=FALSE)
    #devtools::install_git("git://github.com/kfuku52/rkftools.git", branch = "master")
}

library(ape)
library(phytools)
library(doParallel)
library(rkftools)

options(expressions=20000)

if (length(commandArgs(trailingOnly=TRUE))==1) {
    mode = "debug"    
} else {
    mode = "batch"
}

if (mode=="debug") {
    og = 'OG0000000'
    dir_tmp = '/Users/kf/Dropbox/kfdata/02_Data/my_db/Ensembl/release-91/orthogroup/tmp/'
    dirs = list.files(dir_tmp)
    dir_work = paste0(dir_tmp, dirs[grep(og, dirs)], '/')
    dir_notung = paste0(dir_work, og, '.notung.root/')
    infile = paste0(dir_work, og, '.iqtree.nwk')
    outfile = paste0(dir_work, og, '.iqtree.root.nwk')
    nslots = 4
    
    #dir_work = '/Users/kf/Dropbox/kfdata/02_Data/04_Convergence_Duplication/20171114_mitochondria_sequence_retrieval/tree_search/out.mt_concat.rnas.aln.ta.fasta/'
    #infile=paste0(dir_work, 'mt_concat.rnas.aln.ta.iqtree.ghost.contree.nwk')
    #outfile=paste0(dir_work, 'mt_concat.rnas.aln.ta.iqtree.ghost.contree.r.nwk')
} else if (mode=="batch") {
    args = commandArgs(trailingOnly=TRUE)
    dir_notung = args[1]
    infile = args[2]
    outfile = args[3]
    nslots = as.integer(args[4])
}

get_root_position_dependent_species_overlap_scores = function(phy, nslots) {
    species_overlap_scores = vector(mode="numeric", nrow(phy$edge))
    exp_funs = c('reroot', "get_species_overlap_score", "get_duplication_confidence_score",
                "get_children_num", "get_tip_labels")
    num_parallel = min(nslots, nrow(phy$edge))
    cluster = makeCluster(num_parallel, 'PSOCK', outfile='')
    registerDoParallel(cluster)
    so_score = foreach (i = 1:nrow(phy$edge), .combine=c, .export=exp_funs) %dopar% {
        rt = reroot(tree=phy, node.number=phy$edge[i,2])
        so_score = get_species_overlap_score(phy=rt, dc_cutoff=0)
        names(so_score) = i
        so_score
    }
    stopCluster(cluster)
    species_overlap_scores = so_score[order(as.integer(names(so_score)))]
    return(species_overlap_scores)
}

input_tree = read.tree(infile)
unrooted_tree = unroot(input_tree)
midpoint_tree = phytools::midpoint.root(unrooted_tree)
cat('# leaves:', length(unrooted_tree$tip.label), '\n')

# MAD rooting
start = proc.time()
unrooted_newick = write.tree(unrooted_tree)
gc()
res = try(MAD(unrooted_newick=unrooted_newick, output_mode='custom'))
if(class(res)=="try-error") {
    cat(res, '\n')
    cat('MAD cannot be completed correctly. Proceeding.\n')
    res = vector(mode='list', 7)
}
mad_tree = read.tree(text=res[[1]])
end = proc.time()
cat('elapsed time for MAD rooting:', (end -start)[3], 'sec\n')

# species overlap search
start = proc.time()
check_species_overlap_score = FALSE
if (check_species_overlap_score) {
    species_overlap_scores = get_root_position_dependent_species_overlap_scores(phy=unrooted_tree, nslots=nslots)
    cat('# root positions with the minimum species overlap score:', sum(species_overlap_scores==min(species_overlap_scores)), '/', length(species_overlap_scores), '\n')
    ind_min_so = c(1:length(species_overlap_scores))[species_overlap_scores==min(species_overlap_scores)]
    cat('root positions with the minimum species overlap score:', ind_min_so, '\n')
    cat('species overlap score: minimum:', min(species_overlap_scores), '\n')
    cat('species overlap score: MAD:', species_overlap_scores[res[[4]]], '\n')
    midpoint_so_score = get_species_overlap_score(midpoint_tree)
    cat('species overlap score: midpoint:', midpoint_so_score, '\n')
}
if (is.null(res[[5]])) {
    ind_top_ten_mad = NA
} else {
    ind_top_ten_mad = order(res[[5]])[1:10]
    ind_rho_peak = c(1:length(res[[7]]))[(res[[7]]!=0)&(res[[7]]!=1)]
    cat('root positions with rho peak:', ind_rho_peak, '\n')
    cat('top 10 MAD positions:', ind_top_ten_mad, '\n')
}

if (! is.null(mad_tree)) {
    if (is_same_root(mad_tree, midpoint_tree)) {
        cat('MAD and midpoint rootings are consistent.\n')
    } else {
        cat('MAD and midpoint rootings are isconsistent.\n')
    }
}
end = proc.time()
cat('elapsed time for species overlap search:', (end - start)[3], 'sec\n')

# check NOTUNG compatibility
files = list.files(dir_notung)
nwk_files = files[grep('[0-9]$', files)]
notung_roots = vector(mode='numeric', length(nwk_files))

do_parallel = FALSE
if (do_parallel) {
    num_parallel = min(nslots, length(nwk_files))
    cluster = makeCluster(num_parallel, 'PSOCK', outfile='')
    registerDoParallel(cluster)
    exp_funs = c("get_phy2_root_in_phy1", "reroot", "is_same_root")
    notung_roots = foreach (i = 1:length(nwk_files), .combine=c, .export=exp_funs) %dopar% {
        library(ape)
        library(phytools)
        tree = read.tree(paste0(dir_notung, nwk_files[i]))
        notung_root = get_phy2_root_in_phy1(phy1=unrooted_tree, phy2=tree, mode="index")
        notung_root
    }
    stopCluster(cluster)    
    mad_root = res[[4]]
    mid_root = get_phy2_root_in_phy1(phy1=unrooted_tree, phy2=midpoint_tree, mode="index")
    cat('NOTUNG root positions:', notung_roots, '\n')
    cat('MAD root position:', mad_root, '\n')
    cat('midpoint root position:', mid_root, '\n')
    is_mid_compatible_with_notung = ifelse(mid_root %in% notung_roots, TRUE, FALSE)
    is_mad_compatible_with_notung = ifelse(mad_root %in% notung_roots, TRUE, FALSE)
} else {
    if (length(input_tree$tip.label)<1000) {
        notung_roots = vector(mode='numeric', length(nwk_files))
        for (i in 1:length(nwk_files)) {
            cat('processing', i, 'th NOTUNG tree.\n')
            tree = read.tree(paste0(dir_notung, nwk_files[i]))
            notung_root = get_phy2_root_in_phy1(phy1=unrooted_tree, phy2=tree, mode="index")
            notung_roots[i] = notung_root
        }
        mad_root = res[[4]]
        mid_root = get_phy2_root_in_phy1(phy1=unrooted_tree, phy2=midpoint_tree, mode="index")
        cat('NOTUNG root positions:', notung_roots, '\n')
        cat('MAD root position:', mad_root, '\n')
        cat('midpoint root position:', mid_root, '\n')
        is_mid_compatible_with_notung = ifelse(mid_root %in% notung_roots, TRUE, FALSE)
        is_mad_compatible_with_notung = ifelse(mad_root %in% notung_roots, TRUE, FALSE)
    } else {
        is_mid_compatible_with_notung = FALSE
        is_mad_compatible_with_notung = FALSE
        for (i in 1:length(nwk_files)) {
            cat('processing', i, 'th NOTUNG tree.\n')
            notung_tree = read.tree(paste0(dir_notung, nwk_files[i]))
            if (is_same_root(midpoint_tree, notung_tree)) {
                is_mid_compatible_with_notung = TRUE
            }
            if (! is.null(mad_tree)) {
                if (is_same_root(mad_tree, notung_tree)) {
                    is_mad_compatible_with_notung = TRUE
                }
            }
        }
    }
}

if (is_mad_compatible_with_notung) {
    cat('the MAD root position is in NOTUNG root positions. Returning the MAD tree.\n')
    out_tree = mad_tree
} else if (is_mid_compatible_with_notung) {
    cat('the MAD root position is not compatible with NOTUNG results.\n')
    cat('the midpoint root position is found to be compatible with NOTUNG resunts. Returning the midpoint tree.\n')
    out_tree = midpoint_tree
} else {
    cat('Neither MAD nor midpoint tree is compatible with NOTUNG results. Returning the first NOTUNG tree.\n')
    out_tree = read.tree(paste0(dir_notung, nwk_files[1]))
}

# Output
out_tree$node.label = NULL
write.tree(out_tree, outfile)
cat('Tree rooting completed.\n')
    

if (mode=='debug') {
    if (check_species_overlap_score) {
        hist(species_overlap_scores)
    }
    par(ces=0.2)
    plot(out_tree)
}
