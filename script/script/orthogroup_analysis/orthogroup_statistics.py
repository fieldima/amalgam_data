#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os, ete3, numpy, pandas, sys, time, colorsys, re, argparse
mode = 'debug' if sys.argv[1] == "-f" else 'batch'
if mode=='debug':
    sys.path.append('/Users/kef74yk/Dropbox_w/repos/kftools')
    sys.path.append('/Users/kef74yk/Dropbox_w/repos/kftools/kftools')

if False:
    get_ipython().system('pip install -U git+git://github.com/kfuku52/kftools@master')

from kftools.kfog import *

pandas.options.display.max_rows=15
pandas.options.display.max_columns=1000


# In[2]:


# Parameters
# mode = "debug" # "debug" or "batch"
print("Setting parameters...\n")
print('mode =', mode, '\n')

if mode=="debug":
    og_id = 'OG0002332'
    #og_id = 'OG0004940'
    dir_ensembl = '/Users/kef74yk/Dropbox_w/db/Ensembl/release-91/'
    dir_og = os.path.join(dir_ensembl, 'orthogroup/tmp/'+og_id+'/')
    os.chdir(dir_og)
    sys.argv = [
        'orthogroup_statistics.py',
        '--species_tree', os.path.join(dir_ensembl, 'timetree/species_timetree.nwk'),
        '--unaligned_aln', dir_og+og_id+'.cds.fasta',
        '--trimal_aln', dir_og+og_id+'.cds.trimal.fasta',
        '--iqtree_tree', dir_og+og_id+'.iqtree.nwk',
        '--iqtree_model', dir_og+og_id+'.model.gz',
        '--root_tree', dir_og+og_id+'.root.nwk',
        '--root_log', dir_og+og_id+'.root.txt',
        '--notung_root_log', dir_og+og_id+'.notung.root/'+og_id+'.iqtree.nwk.rooting.ntglog',
        '--notung_reconcil_stats', dir_og+og_id+'.notung.reconcil/'+og_id+'.root.nwk.reconciled.stats.txt',
        '--dated_tree', dir_og+og_id+'.dated.nwk',
        '--dated_log', dir_og+og_id+'.dated.log.txt',
        '--hyphy_tree_dnds', dir_og+og_id+'.hyphy.dnds.nwk',    

        '--l1ou_prefix', 'tpm', 'fpkm',
        '--l1ou_tree', dir_og+og_id+'.l1ou.tree.tsv', dir_og+og_id+'.fpkm.l1ou.tree.tsv',
        '--l1ou_regime', dir_og+og_id+'.l1ou.regime.tsv', dir_og+og_id+'.fpkm.l1ou.regime.tsv',
        '--l1ou_leaf', dir_og+og_id+'.l1ou.leaf.tsv', dir_og+og_id+'.fpkm.l1ou.leaf.tsv',

        '--phylogeneticem_prefix', 'tpm', 'fpkm',
        '--phylogeneticem_tree', dir_og+og_id+'.PhylogeneticEM.tree.tsv', dir_og+og_id+'.fpkm.PhylogeneticEM.tree.tsv',
        '--phylogeneticem_regime', dir_og+og_id+'.PhylogeneticEM.regime.tsv', dir_og+og_id+'.fpkm.PhylogeneticEM.regime.tsv',
        '--phylogeneticem_leaf', dir_og+og_id+'.PhylogeneticEM.leaf.tsv', dir_og+og_id+'.fpkm.PhylogeneticEM.leaf.tsv',   
        
        '--expression_prefix', 'tpm', 'fpkm',
        '--expression', dir_og+og_id+'.expression.tsv', dir_og+og_id+'.expression.fpkm.tsv',

        '--mapdnds_tree_dn', dir_og+og_id+'.mapdNdS.dN.nwk',
        '--mapdnds_tree_ds', dir_og+og_id+'.mapdNdS.dS.nwk',
        '--scm_intron', dir_og+og_id+'.scm.intron.tsv',
        '--scm_chromosome', dir_og+og_id+'.scm.chromosome.tsv',
    ]


# In[3]:


parser = argparse.ArgumentParser()
parser.add_argument('--species_tree', metavar='PATH', default='', type=str, help='')
parser.add_argument('--unaligned_aln', metavar='PATH', default='', type=str, help='')
parser.add_argument('--trimal_aln', metavar='PATH', default='', type=str, help='')
parser.add_argument('--iqtree_tree', metavar='PATH', default='', type=str, help='')
parser.add_argument('--iqtree_model', metavar='PATH', default='', type=str, help='')
parser.add_argument('--root_tree', metavar='PATH', default='', type=str, help='')
parser.add_argument('--root_log', metavar='PATH', default='', type=str, help='')
parser.add_argument('--notung_root_log', metavar='PATH', default='', type=str, help='')
parser.add_argument('--notung_reconcil_stats', metavar='PATH', default='', type=str, help='')
parser.add_argument('--dated_tree', metavar='PATH', default='', type=str, help='')
parser.add_argument('--dated_log', metavar='PATH', default='', type=str, help='')
parser.add_argument('--hyphy_tree_dnds', metavar='PATH', default='', type=str, help='')

parser.add_argument('--l1ou_prefix', metavar='PATH', default='', type=str, nargs='+', help='')
parser.add_argument('--l1ou_tree', metavar='PATH', default='', type=str, nargs='+', help='')
parser.add_argument('--l1ou_regime', metavar='PATH', default='', type=str, nargs='+', help='')
parser.add_argument('--l1ou_leaf', metavar='PATH', default='', type=str, nargs='+', help='')
parser.add_argument('--phylogeneticem_prefix', metavar='PATH', default='', type=str, nargs='+', help='')
parser.add_argument('--phylogeneticem_tree', metavar='PATH', default='', type=str, nargs='+', help='')
parser.add_argument('--phylogeneticem_regime', metavar='PATH', default='', type=str, nargs='+', help='')
parser.add_argument('--phylogeneticem_leaf', metavar='PATH', default='', type=str, nargs='+', help='')
parser.add_argument('--expression_prefix', metavar='PATH', default='', type=str, nargs='+', help='')
parser.add_argument('--expression', metavar='PATH', default='', type=str, nargs='+', help='')

parser.add_argument('--mapdnds_tree_dn', metavar='PATH', default='', type=str, help='')
parser.add_argument('--mapdnds_tree_ds', metavar='PATH', default='', type=str, help='')
parser.add_argument('--scm_intron', metavar='PATH', default='', type=str, help='')
parser.add_argument('--scm_chromosome', metavar='PATH', default='', type=str, help='')

args = parser.parse_args()  
params = dict()
for attr in [a for a in dir(args) if not a.startswith('_')]:
    params[attr] = getattr(args, attr)

keys = list(params.keys())
for key in keys:
    print(key, '=', params[key])
print('\n')


# In[4]:


df_branch = pandas.DataFrame({'numerical_label':[0,]})

if os.path.exists(params["dated_tree"]):
    df_tmp = nwk2table(tree=params["dated_tree"], mode='node_name')
    df_branch = pandas.merge(df_branch, df_tmp, on='numerical_label', how='outer')
    df_tmp = get_misc_node_statistics(tree_file=params["dated_tree"])
    df_branch = pandas.merge(df_branch, df_tmp, on='numerical_label', how='outer')
if os.path.exists(params["root_tree"]):
    root_tree = ete3.PhyloNode(params['root_tree'], format=1)
    df_tmp = nwk2table(tree=root_tree, mode='branch_length')
    df_tmp.columns = ['numerical_label', 'bl_iqtree']
    df_branch = pandas.merge(df_branch, df_tmp, on='numerical_label', how='outer')
if (os.path.exists(params["iqtree_tree"]))&(os.path.exists(params["root_tree"])):
    iqtree_tree = ete3.PhyloNode(params['iqtree_tree'], format=0)
    iqtree_tree = transfer_root(tree_to=iqtree_tree, tree_from=root_tree)
    df_tmp = nwk2table(tree=iqtree_tree, mode='branch_support')
    df_tmp.columns = ['numerical_label', 'support_iqtree']
    df_branch = pandas.merge(df_branch, df_tmp, on='numerical_label', how='outer')
if os.path.exists(params["dated_tree"]):
    df_tmp = nwk2table(tree=params["dated_tree"], mode='branch_length', age=True)
    df_tmp.columns = ['numerical_label', 'bl_dated', 'age']
    df_branch = pandas.merge(df_branch, df_tmp, on='numerical_label', how='outer')
if (os.path.exists(params["dated_tree"]))&(os.path.exists(params["species_tree"])):
    gene_tree = ete3.PhyloNode(params['dated_tree'], format=1)
    species_tree = ete3.PhyloNode(params['species_tree'], format=1)
    df_tmp = node_gene2species(gene_tree, species_tree)
    df_branch = pandas.merge(df_branch, df_tmp, on='numerical_label', how='outer')
if (all([ os.path.exists(params[key]) for key in ['mapdnds_tree_ds','mapdnds_tree_dn'] ])):
    df_tmp = nwk2table(tree=params["mapdnds_tree_ds"], mode='branch_length')
    df_tmp.columns = ['numerical_label', 'mapdnds_ds']
    df_branch = pandas.merge(df_branch, df_tmp, on='numerical_label', how='outer')
    df_tmp = nwk2table(tree=params["mapdnds_tree_dn"], mode='branch_length')
    df_tmp.columns = ['numerical_label', 'mapdnds_dn']
    df_branch = pandas.merge(df_branch, df_tmp, on='numerical_label', how='outer')
    df_branch['mapdnds_omega'] = df_branch['mapdnds_dn'] / df_branch['mapdnds_ds']
if (all([ os.path.exists(params[key]) for key in ['hyphy_tree_dnds'] ])):
    # 1. E[Syn subs/nucleotide site] tree
    # 2. E[Non-syn subs/nucleotide site] tree
    # 3. dS tree
    # 4. dN tree
    with open(params['hyphy_tree_dnds']) as file:
        hyphy_trees = [ line for line in file.readlines() ]
    hyphy_tree_ds = ete3.PhyloNode(hyphy_trees[2], format=1)
    df_tmp = nwk2table(tree=hyphy_tree_ds, mode='branch_length')
    df_tmp.columns = ['numerical_label', 'hyphy_ds']
    df_branch = pandas.merge(df_branch, df_tmp, on='numerical_label', how='outer')
    hyphy_tree_dn = ete3.PhyloNode(hyphy_trees[3], format=1)
    df_tmp = nwk2table(tree=hyphy_tree_dn, mode='branch_length')
    df_tmp.columns = ['numerical_label', 'hyphy_dn']
    df_branch = pandas.merge(df_branch, df_tmp, on='numerical_label', how='outer')
    df_branch['hyphy_omega'] = df_branch['hyphy_dn'] / df_branch['hyphy_ds']
if (all([ os.path.exists(params[key]) for key in ['scm_intron'] ])):
    df_tmp = pandas.read_csv(params['scm_intron'], sep='\t',  header=0, index_col=None)
    df_tmp.columns = [ c.replace('leaf','node_name') for c in df_tmp.columns ]
    df_branch = pandas.merge(df_branch, df_tmp, on='node_name', how='outer')
    df_branch = df_branch.drop('intron_absent', axis=1)
    df_branch = compute_delta(df_branch, 'intron_present')
if (all([ os.path.exists(params[key]) for key in ['scm_chromosome'] ])):
    df_tmp = pandas.read_csv(params['scm_chromosome'], sep='\t',  header=0, index_col=None)
    df_tmp.columns = [ c.replace('leaf','node_name') for c in df_tmp.columns ]
    df_branch = pandas.merge(df_branch, df_tmp, on='node_name', how='outer')
    df_branch = compute_delta(df_branch, 'A')
    df_branch = compute_delta(df_branch, 'X')
    df_branch = compute_delta(df_branch, 'Y')
    df_branch = compute_delta(df_branch, 'MT')
for i,prefix in enumerate(params['l1ou_prefix']):
    if (os.path.exists(params["l1ou_regime"][i]))&(os.path.exists(params["l1ou_leaf"][i])):
        print('processing l1ou', prefix)
        df_tmp = ou2table(params['l1ou_regime'][i], params['l1ou_leaf'][i], params["dated_tree"])
        df_tmp.columns = [ 'l1ou_'+prefix+'_'+c if c!='numerical_label' else c for c in df_tmp.columns ]
        df_branch = pandas.merge(df_branch, df_tmp, on='numerical_label', how='outer')
for i,prefix in enumerate(params['phylogeneticem_prefix']):    
    if (os.path.exists(params["phylogeneticem_regime"][i]))&(os.path.exists(params["phylogeneticem_leaf"][i])):
        print('processing PhylogeneticEM', prefix)
        df_tmp = ou2table(params['phylogeneticem_regime'][i], params['phylogeneticem_leaf'][i], params["dated_tree"])
        df_tmp.columns = [ 'phylogeneticem_'+prefix+'_'+c if c!='numerical_label' else c for c in df_tmp.columns ]
        df_branch = pandas.merge(df_branch, df_tmp, on='numerical_label', how='outer')
for i,prefix in enumerate(params['expression_prefix']):
    col = params['expression_prefix'][i]+'_clade_min_pearsoncor'
    df_branch.loc[:,col] = numpy.nan
    df_exp = pandas.read_csv(params['expression'][i], sep='\t', index_col=0)
    df_cor = numpy.corrcoef(df_exp)
    gene_tree = ete3.PhyloNode(params['dated_tree'], format=1)
    gene_tree = add_numerical_node_labels(gene_tree)
    for node in gene_tree.traverse(strategy='postorder'):
        if not node.is_leaf():
            leaf_names = node.get_leaf_names()
            leaf_loc = sorted([ df_exp.index.get_loc(ln) for ln in leaf_names ])
            node.min_pearsoncor = df_cor[leaf_loc,:][:,leaf_loc].min()
            df_branch.loc[(df_branch['numerical_label']==node.numerical_label),col] = node.min_pearsoncor
df_branch.to_csv('orthogroup.branch.tsv', sep='\t', index=False)

if (mode=='debug'):
    import IPython.display
    IPython.display.display(df_branch)


# In[5]:


tree_info = dict()
tree_info.update(branch2tree(df_branch))
tree_info['dating_method'] = get_dating_method(params['dated_log'])
if os.path.exists(params["unaligned_aln"]):
    tree_tmp = get_aln_stats(params['unaligned_aln'])
    tree_tmp = add_dict_key_prefix(tree_tmp, 'original')
    tree_info.update(tree_tmp)
if os.path.exists(params["trimal_aln"]):
    tree_tmp = get_aln_stats(params['trimal_aln'])
    tree_tmp = add_dict_key_prefix(tree_tmp, 'cleaned')
    tree_info.update(tree_tmp)
if os.path.exists(params["root_log"]):
    tree_tmp = get_root_stats(params['root_log'])
    tree_info.update(tree_tmp)
if os.path.exists(params["notung_root_log"]):
    tree_tmp = get_notung_root_stats(params['notung_root_log'])
    tree_info.update(tree_tmp)
if os.path.exists(params["notung_reconcil_stats"]):
    tree_tmp = get_notung_reconcil_stats(params['notung_reconcil_stats'])
    tree_info.update(tree_tmp)
if os.path.exists(params["iqtree_model"]):
    tree_tmp = get_iqtree_model_stats(params['iqtree_model'])
    tree_info.update(tree_tmp)
for method in ['l1ou','phylogeneticem']:
    for i,prefix in enumerate(params[method+'_prefix']):
        if os.path.exists(params[method+"_tree"][i]):
            tmp = pandas.read_csv(params[method+'_tree'][i], sep='\t')
            num_shift = tmp['num_shift'].values[0]
            tree_tmp = {method+'_'+prefix+'_num_shift':num_shift,}
            tree_info.update(tree_tmp)
        if os.path.exists(params[method+"_regime"][i]):
            tree_tmp = regime2tree(params[method+'_regime'][i])
            tree_tmp = add_dict_key_prefix(tree_tmp, method+'_'+prefix)
            tree_info.update(tree_tmp)

if mode == 'debug':
    for item in tree_info.items():
        print(item[0], '=', item[1])

df_tree = pandas.DataFrame(tree_info, index=[0, ])
df_tree.to_csv('orthogroup.tree.tsv', sep='\t', index=False)


# In[6]:


print('orthogroup_statistics: done!')


# In[ ]:





# In[ ]:




