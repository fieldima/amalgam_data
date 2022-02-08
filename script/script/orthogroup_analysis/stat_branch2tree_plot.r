
mode = ifelse(length(commandArgs(trailingOnly=TRUE))==1, 'debug', 'batch')

if (mode=='debug') {
    #devtools::install_github("kfuku52/rkftools", dep=FALSE)
    #remove.packages('rkftools')
    devtools::install_local(path='/Users/kef74yk/Dropbox_w/repos/rkftools', reload=TRUE, quick=FALSE, local=TRUE, dep=FALSE)
    #options(warn=-1)
}

if (FALSE) {
    devtools::install_git("git://github.com/kfuku52/rkftools.git", branch = "master", dep=FALSE)
}

library(ape)
library(ggplot2)
library(ggimage)
library(ggtree)
library(rkftools)
options(stringsAsFactors=FALSE)

if (mode=="debug") {
    #og = 'OG0002332' # PGK
    og = 'OG0000009'
    #og = 'OG0000020'
    dir_tmp = '/Users/kef74yk/Dropbox_w/db/Ensembl/release-91/orthogroup/tmp/'
    #dir_tmp = '/Users/kef74yk/Dropbox_p/mycoheterotroph/kenji_fukushima/gfe/orthogroup/tmp/'
    dirs = list.files(dir_tmp)
    dir_work = paste0(dir_tmp, dirs[grep(og, dirs)], '/')
    files = list.files(dir_work)
    setwd(dir_work)
    args = c()
    args = c(args, paste0('--stat_branch=', dir_work, og, '.stat.branch.tsv'))
    args = c(args, paste0('--trait_file=', dir_work, og, '.expression.tsv'))
    args = c(args, paste0('--dir_phylopic=', '/Users/kf/Dropbox/kfdata/02_Data/my_db/Ensembl/release-91/phylopic/png_files/'))
    args = c(args, paste0('--max_delta_intron_present=', '-0.5'))
    args = c(args, paste0('--pcm_prefix=', 'l1ou_fpkm_'))
    #args = c(args, paste0('--pcm_prefix=', 'l1ou_cpm_'))
} else if (mode=="batch") {
    args = commandArgs(trailingOnly=TRUE)
}

cat('arguments:\n')
args = rkftools::get_parsed_args(args, print=TRUE)

get_single_branch_tree = function(name, dist) {
    phy = list(
        edge = matrix(c(2,1),1,2),
        tip.label = name,
        edge.length = dist,
        Nnode = 1
    )
    class(phy) = "phylo" 
    return(phy)
}

remove_redundant_root_edge = function(phy) {
    root_num = get_root_num(phy)
    is_root_edge = (phy$edge[,1]!=root_num)
    phy$edge = phy$edge[is_root_edge,]
    phy$edge[phy$edge>root_num] = phy$edge[phy$edge>root_num] - 1
    phy$edge.length = phy$edge.length[is_root_edge]
    phy$Nnode = phy$Nnode - 1
    return(phy)
}

table2phylo = function(df, name_col, dist_col) {
    root_id = max(df$numerical_label)
    df[(df$numerical_label==root_id), 'sister'] = -999
    df[(df$numerical_label==root_id), 'parent'] = -999
    root_name = df[(df$numerical_label==root_id), name_col]
    root_dist = df[(df$numerical_label==root_id), dist_col]
    phy = get_single_branch_tree(root_name, root_dist)
    next_node_ids = sort(df[(df$parent==root_id),'numerical_label'])
    while (length(next_node_ids)>0) {
        for (nni in next_node_ids) {
            nni_name = df[(df$numerical_label==nni), name_col]
            nni_dist = df[(df$numerical_label==nni), dist_col]
            parent_id = df[(df$numerical_label==nni), 'parent']
            parent_name = df[(df$numerical_label==parent_id), name_col]
            parent_num = get_node_num_by_name(phy, parent_name)
            parent_index = (1:nrow(phy$edge))[phy$edge[,2]==parent_num]
            sister_id = df[(df$numerical_label==nni), 'sister']
            sister_name = df[(df$numerical_label==sister_id), name_col]
            sister_dist = df[(df$numerical_label==sister_id), dist_col]
            sister_num = get_node_num_by_name(phy, sister_name)
            sister_index = (1:nrow(phy$edge))[phy$edge[,2]==sister_num]
            branch = get_single_branch_tree(nni_name, nni_dist)
            if (length(parent_num)==0) {
                sister_dist = ifelse(sister_dist<1e-8, 1e-8, sister_dist)
                phy = ape::bind.tree(phy, branch, where=sister_num, position=sister_dist)
            } else {
                phy = ape::bind.tree(phy, branch, where=parent_num, position=0)
            }
        }
        next_node_ids = df[(df$parent %in% next_node_ids),'numerical_label']
    }
    phy = remove_redundant_root_edge(phy)
    phy = ape::ladderize(phy, right=TRUE)
    num_leaf = length(phy$tip.label)
    num_intnode = nrow(phy$edge) - length(phy$tip.label)
    phy$node.label = rep('placeholder', num_intnode)
    next_node_ids = sort(df[(df[[name_col]] %in% phy$tip.label), 'numerical_label'])
    while ((length(next_node_ids)!=1)|(next_node_ids[1]>=0)) {
        tmp_next_node_ids = c()
        for (nni in next_node_ids) {
            if (nni>=0) {
                nni_name = df[(df$numerical_label==nni), name_col]
                nni_num = (1:max(phy$edge))[c(phy$tip.label, phy$node.label)==nni_name]
                parent_num = phy$edge[(phy$edge[,2]==nni_num),1]
                parent_label_index = parent_num - num_leaf
                parent_id = df[(df$numerical_label==nni), 'parent']
                if (parent_id>=0) {
                    parent_name = df[(df$numerical_label==parent_id), name_col]
                    phy$node.label[parent_label_index] = parent_name
                    tmp_next_node_ids = c(tmp_next_node_ids, parent_id)
                }
            }
        }
        next_node_ids = sort(unique(tmp_next_node_ids))
        if (length(next_node_ids)==0) {
            next_node_ids = c(-999)
        }
    }
    if (sum(phy$node.label=='placeholder')>1) {
        warning('Node label "placeholder" appeared more than once.')
    }
    return(phy)
}


attach_neighbor_stats = function(df, neighbor, columns=c()){
    columns = c('numerical_label', columns)
    df_tmp = df[,columns]
    colnames(df_tmp) = paste0(neighbor, '_', colnames(df_tmp))
    df = merge(df, df_tmp, by.x=neighbor, by.y=paste0(neighbor,'_numerical_label'), all.x=TRUE)
    return(df)
}

get_img_info = function(b, dir_phylopic) {
    img_info = data.frame(label=b[b$so_event=='L','node_name'], file=NA)
    phylopic_files = list.files(dir_phylopic)
    for (pf in phylopic_files) {
        sci_name = sub('_', 'PLACEHOLDER', pf)
        sci_name = sub('_.*', '', sci_name)
        sci_name = sub('PLACEHOLDER', '_', sci_name)
        file_path = paste0(dir_phylopic, pf)
        img_info[startsWith(img_info[['label']], sci_name),'file'] = file_path
    }
    img_info = img_info[(!is.na(img_info$file)),]
    return(img_info)
}

geom_divtime = function(g, y=0, size=0.5, font.size=8, unit='') {
    options(scipen=18)
    xmax = max(g$data$x)
    if (xmax>1) {
        ndigit = nchar(round(xmax, digits=0)) - 1
    } else {
        ndigit = -nchar(sub('[1-9].*', '', xmax)) - 1
    }
    xunit = 10^ndigit
    if (xmax%/%(10^ndigit)==1) {
        ndigit = ndigit - 1
        xunit = xunit / 2
    }
    yunit = (max(g$data$y) - min(g$data$y))/100
    out = list()
    out[[1]] = ggplot2::geom_segment(x=0, xend=xmax, y=y, yend=y, size=size)
    count = 2
    for (i in 1:(xmax%/%xunit+1)) {
        x = max(0, xmax-(xunit*i))
        xend = xmax-(xunit*(i-1))
        out[[count+0]] = ggplot2::geom_segment(x=xend, xend=xend, y=y-yunit, yend=y, size=size)
        out[[count+1]] = ggplot2::annotate('text', x=xend, y=y-yunit, label=xunit*(i-1), vjust=1.2, size=font_size)
        count = count+2
    }
    out[[count+1]] = ggplot2::annotate('text', x=xmax/2, y=y-yunit, label=unit, vjust=2.75, size=font_size)
    return(out)
}

suppress_annotation = function(df, suppressed_leaf=c(), target_column='', fill='-') {
    for (node_name in df$node_name) {
        sci_name = sub('_', ' ', node_name)
        sci_name = sub('_.*', '', sci_name)
        if (sci_name %in% suppressed_leaf) {
            df[(df$node_name==node_name),target_column] = fill
        }
    }
    return(df)
}

b = read.table(args[['stat_branch']], header=TRUE, sep='\t', stringsAsFactors=FALSE)
shift_prefix = args[['pcm_prefix']]

# branch and node categories
b[['branch_category']] = ''
b[['is_parent_dup']] = (b[['so_event_parent']]=='D')
b[(!b[['is_parent_dup']]),'branch_category'] = 'S'
if ('delta_intron_present' %in% colnames(b)) {
    cat('Nodes will be categorized into: S, D, R', '\n')
    b = attach_neighbor_stats(df=b, neighbor='sister', columns='delta_intron_present')
    b[['is_retrotransposition']] = (b[['delta_intron_present']] <= args[['max_delta_intron_present']])
    b[['is_retrotransposition']][is.na(b[['is_retrotransposition']])] = FALSE
    b[['is_sister_retrotransposition']] = (b[['sister_delta_intron_present']] <= args[['max_delta_intron_present']])
    b[['is_sister_retrotransposition']][is.na(b[['is_sister_retrotransposition']])] = FALSE
    b[['is_lower_delta_intron_present']] = (b[['delta_intron_present']] <= b[['sister_delta_intron_present']])    
    b[(b[['is_parent_dup']])&(!b[['is_retrotransposition']])&(!b[['is_sister_retrotransposition']]),'branch_category'] = 'D'
    b[(b[['is_parent_dup']])&(b[['is_retrotransposition']])&(!b[['is_sister_retrotransposition']]),'branch_category'] = 'R'
} else {
    cat('Nodes will be categorized into: S, D', '\n')
    b[(b[['is_parent_dup']]),'branch_category'] = 'D'
}

b = attach_neighbor_stats(df=b, neighbor='child1', columns='branch_category')
b = attach_neighbor_stats(df=b, neighbor='child2', columns='branch_category')
b[is.na(b[['child1_branch_category']]),'child1_branch_category'] = ''
b[is.na(b[['child2_branch_category']]),'child2_branch_category'] = ''
b[['node_category']] = NaN
b[(b[['child1_branch_category']]=='S')|(b[['child2_branch_category']]=='S'),'node_category'] = 'S'
b[(b[['child1_branch_category']]=='D')|(b[['child2_branch_category']]=='D'),'node_category'] = 'D'
b[(b[['child1_branch_category']]=='R')|(b[['child2_branch_category']]=='R'),'node_category'] = 'R'
#b[['child1_branch_category']] = NULL
#b[['child2_branch_category']] = NULL
cat('detected node categories:', unique(b[['node_category']][b[['node_category']]!='NaN']), '\n')

# chromosome
is_autosome = (!b[['chromosome']] %in% c('X','Y','W','Z','MT'))
is_autosome[is.na(is_autosome)] = FALSE
b[is_autosome,'chromosome'] = 'A'
nonmammalian_spp = c(
    'Ornithorhynchus anatinus',
    'Anolis carolinensis',
    'Gallus gallus',
    'Xenopus tropicalis',
    'Oreochromis niloticus',
    'Oryzias latipes',
    'Gadus morhua',
    'Astyanax mexicanus',
    'Danio rerio'
)
fragmented_genome_spp = c(
    'Chinchilla lanigera'
)
chromosome_suppressed_spp = c(nonmammalian_spp, fragmented_genome_spp)
b = suppress_annotation(df=b, suppressed_leaf=chromosome_suppressed_spp, target_column='chromosome', fill='-')

# dNdS
omega_method = 'mapdnds'
omega_column = paste0(omega_method, '_omega')
b = attach_neighbor_stats(df=b, neighbor='sister', columns=omega_column)
b[['is_higher_dnds']] = (b[[omega_column]]>=b[[paste0('sister_', omega_column)]])

# Branch width
b[['branch_thickness']] = 0.5
#b[(b[['is_higher_dnds']]&b[['is_parent_dup']]),'branch_thickness'] = 1

# trait table
trait_table = read.table(args[['trait_file']], header=TRUE, row.names=1, sep="\t")
#colnames(trait_table) = toupper(substr(colnames(trait_table),1,1))
#img_info = get_img_info(b, args[['dir_phylopic']])

show_mu_heatmap = FALSE
font_size = 3
branch_thickness = 0.5
pie_size = max(b[['age']])/4
trait_colors = c('#1B9E77','#D95F02','#7570B3','#E7298A','#66A61E','#E6AB02')
node_colors = c('D'='red', 'R'='orange', 'S'='blue')
cat('node_colors:', node_colors, '\n')
cat('names(node_colors):', names(node_colors), '\n')

mu_cols = colnames(b)[startsWith(colnames(b), paste0(shift_prefix, 'mu_'))]
mu_cols = mu_cols[mu_cols!=paste0(shift_prefix, 'mu_complementarity')]

tree = table2phylo(df=b, name_col='node_name', dist_col='bl_dated')
tree = rkftools::force_ultrametric(tree, stop_if_larger_change=0.1)
tree_height = get_node_age(tree, get_root_num(tree))

g = ggtree(tree, size=0, layout='rectangular')
g$data = merge(g$data, b, by.x='label', by.y='node_name', all.x=TRUE)
colnames(g$data) = sub('\\.x$', '', colnames(g$data))
regime_nos = sort(unique(g$data[[paste0(shift_prefix, 'regime')]]))
regime_colors = c('#000000', colorspace::rainbow_hcl(length(regime_nos)-1, c=100))
g$data = merge(g$data, data.frame(regime=regime_nos, regime_color=regime_colors, stringsAsFactors=FALSE), by.x=paste0(shift_prefix, 'regime'), by.y='regime')
g$data = g$data[order(g$data$node),]
g$data[,mu_cols][g$data[,mu_cols]<0] = 0
rownames(g$data) = NULL

shift_node_nums = c(get_root_num(tree), g$data[g$data[[paste0(shift_prefix,'is_shift')]]==1,'node'])
pie_data = g$data[(g$data$node%in%shift_node_nums), c(mu_cols,'node')]
names(trait_colors) = mu_cols
pies = nodepie(data=pie_data, cols=1:length(mu_cols), color=trait_colors, alpha=0.7)
g = inset(tree_view=g, insets=pies, width=pie_size, height=pie_size, x='branch', hjust=0, vjust=0)

g = g + geom_tree(color=g$data$regime_color, size=g$data$branch_thickness)
##g = g + geom_tree(color=g$data$regime_color, size=branch_thickness)
if ('num_intron' %in% colnames(g$data)) {
    g = g + geom_tiplab(label=g$data$num_intron[g$data$isTip], offset=tree_height*0.6, color='black', size=font_size, align=TRUE, linetype='blank', hjust=0.5)
}
if ('chromosome' %in% colnames(g$data)) {
    g = g + geom_tiplab(label=g$data$chromosome[g$data$isTip], offset=tree_height*0.68, color='black', size=font_size, align=TRUE, linetype='blank', hjust=0.5)
}
g = g + geom_tiplab(label=g$data$label[g$data$isTip], offset=tree_height*0.75, color=g$data$regime_color[g$data$isTip], size=font_size, align=TRUE, linetype='blank', linesize=0.5)
g = g + geom_nodelab(mapping=aes(x=branch, y=y, label=support_iqtree), nudge_x=0, nudge_y=0.4, geom='text', hjust=0.5, vjust=0.5, size=font_size)
##g = g + geom_nodelab(mapping=aes(x=branch, y=y, label=node), nudge_x=0, nudge_y=0.4, geom='text', hjust=0.5, vjust=0.5, size=font_size)
g = g + geom_nodepoint(aes(subset=!isTip, color=node_category), position='identity', show.legend=FALSE) + scale_colour_manual(values=node_colors)
## g = annotation_image(g, img_info)
g = gheatmap(g, trait_table, offset=0, width=0.5, low='blue', high='red', colnames_position='top', font.size=font_size)
if (show_mu_heatmap) {
    mu_table = g$data[(g$data$isTip),mu_cols]
    rownames(mu_table) = g$data$label[g$data$isTip]
    colnames(mu_table) = colnames(trait_table)
    g = gheatmap(g, mu_table, offset=max(g$data$x)*0.52, width=0.5, low='blue', high='red', colnames_position='top', font.size=font_size)
}

g = g + ggplot2::xlim(0, tree_height*3.2)
g = g + geom_divtime(g, y=0, size=branch_thickness, font.size=font_size, unit='Million years ago')
if (mode=='debug') {
    plot(g)
}

height = nrow(g$data) / 10
if (mode=='debug') {
    ggsave(g, file='stat_branch2tree_plot.pdf', width=7.2, height=height*0.7, units='in', dpi=300, limitsize=FALSE)
} else {
    ggsave(g, file='stat_branch2tree_plot.pdf', width=7.2, height=height, units='in', dpi=300, limitsize=FALSE)
}



