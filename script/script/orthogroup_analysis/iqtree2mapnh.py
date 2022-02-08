
# coding: utf-8

# In[1]:


import argparse, pandas, sys, os, re, ete3
mode = 'debug' if sys.argv[1] == "-f" else 'batch'

if mode=='debug':
    sys.path.append('/Users/kf/Dropbox/kfdata/02_Data/my_projects/kftools')
    sys.path.append('/Users/kf/Dropbox/kfdata/02_Data/my_projects/kftools/kftools')

if False:
    get_ipython().system('pip install -U git+git://github.com/kfuku52/kftools@master')

from kftools.kfphylo import *
from kftools.kfseq import *


# In[2]:


if mode=="debug":
    og_id = 'OG0001148'#'OG0005008'
    dir_work = '/Users/kf/Dropbox/kfdata/02_Data/my_db/Ensembl/release-91/orthogroup/tmp/'+og_id+'/'
    os.chdir(dir_work)
    args = dict()
    args['iqtree'] = dir_work+og_id+'.mapdNdS/'+og_id+'.iqtree2mapdNdS.iqtree'
    args['log'] = dir_work+og_id+'.mapdNdS/'+og_id+'.iqtree2mapdNdS.log'
    args['state'] = dir_work+og_id+'.mapdNdS/'+og_id+'.iqtree2mapdNdS.state'
    args['alignment'] = dir_work+og_id+'.cds.trimal.fasta'
    args['treefile'] = dir_work+og_id+'.mapdNdS/'+og_id+'.iqtree2mapdNdS.treefile'
    args['rooted_tree'] = dir_work+og_id+'.root.nwk'
elif mode=="batch":
    parser = argparse.ArgumentParser()
    parser.add_argument('--iqtree', metavar='PATH',type=str, required=True, help='.iqtree output from IQ-TREE')
    parser.add_argument('--log', metavar='PATH',type=str, required=True, help='.log output from IQ-TREE')
    parser.add_argument('--state', metavar='PATH',type=str, required=True, help='.state output from IQ-TREE')
    parser.add_argument('--alignment', metavar='PATH',type=str, required=True, help='input alignment for IQ-TREE')
    parser.add_argument('--treefile', metavar='PATH',type=str, required=True, help='.treefile output from IQ-TREE')
    parser.add_argument('--rooted_tree', metavar='PATH',type=str, required=True, help='A rooted newick tree')
    args = parser.parse_args()
    g = dict()
    for attr in [a for a in dir(args) if not a.startswith('_')]:
        g[attr] = getattr(args, attr)
    args = g
    
print('args:', args)


# In[3]:


with open(args['log'], "r") as file:
    text_log = file.read().split('\n')    
line_omega = [ t for t in text_log if 'Nonsynonymous/synonymous ratio' in t ][0]
line_kappa = [ t for t in text_log if 'Transition/transversion ratio' in t ][0]
line_alpha = [ t for t in text_log if 'Gamma shape alpha' in t ][0]
line_command = [ t for t in text_log if 'Command: ' in t ][0]
value_omega = float(re.sub(r'.*: ', '', line_omega))
value_kappa = float(re.sub(r'.*: ', '', line_kappa))
value_alpha = float(re.sub(r'.*: ', '', line_alpha))
value_model = re.sub(r' .*', '', re.sub(r'.*-m ', '', line_command))
print('omega =', value_omega)
print('kappa =', value_kappa)
print('alpha =', value_alpha)
print('model =', value_model)


# In[4]:


with open(args['iqtree'], "r") as file:
    text_iqtree = file.read().split('\n')
line_pi = [ t for t in text_iqtree if re.search('pi\(...\)', t) ]
value_pi = [ item for line in line_pi for item in line.replace(' ', '').split('pi') if item!='' ]
value_pi = [ item.replace('(','').replace(')','') for item in value_pi ]
codon_freqs = dict()
for pi in value_pi:
    codon_freqs[pi.split('=')[0]] = float(pi.split('=')[1])

nuc_freqs = codon2nuc_freqs(codon_freqs=codon_freqs, model=value_model)
print('equilibrium nucleotide frequency:', nuc_freqs)

thetas = nuc_freq2theta(nuc_freqs=nuc_freqs)
print('equilibrium theta:', thetas)


# In[5]:


rooted_tree = ete3.PhyloNode(args['rooted_tree'])
iqtree_tree = ete3.PhyloNode(args['treefile'], format=1)
iqtree_tree = transfer_root(tree_to=iqtree_tree, tree_from=rooted_tree)

anc_state = pandas.read_csv(args['state'], sep='\t', comment='#')
subroot_node_names = [ n.name for n in iqtree_tree.get_children() ]

subroot_states = dict()
subroot_codon_freqs = dict()
subroot_nuc_freqs = dict()
subroot_thetas = dict()
for snn in subroot_node_names:
    subroot_state = anc_state.loc[(anc_state['Node']==snn),:]
    codon_columns = subroot_state.columns[subroot_state.columns.str.startswith('p_')]
    subroot_codon_freq = subroot_state.loc[:,codon_columns].mean(axis=0)
    subroot_codon_freq.index = codon_columns.str.replace('p_', '')
    subroot_codon_freq = subroot_codon_freq.to_dict()
    subroot_states[snn] = subroot_state
    subroot_codon_freqs[snn] = subroot_codon_freq
    subroot_nuc_freqs[snn] = codon2nuc_freqs(codon_freqs=subroot_codon_freq, model=value_model)
    if (subroot_states[snn].size==0):
        print('a subroot node is a leaf:', snn)
        subroot_nuc_freqs[snn] = alignment2nuc_freqs(leaf_name=snn, alignment_file=args['alignment'], model=value_model)
    subroot_thetas[snn] = nuc_freq2theta(nuc_freqs=subroot_nuc_freqs[snn])

print('subroot nucleotide frequency:', subroot_nuc_freqs)
root_thetas = weighted_mean_root_thetas(subroot_thetas, iqtree_tree, model=value_model)
print('root theta:', root_thetas)


# In[6]:


num_gamma_category = value_model
num_gamma_category = re.sub(r'.*\+G', '', value_model)

out = str()
out = out+'alphabet=Codon(letter=DNA)'+'\n'
out = out+'genetic_code=Standard'+'\n'
out = out+'input.data1=alignment(file=$(SEQ), sites_to_use=all, remove_stop_codons=yes)'+'\n'
out = out+'input.tree1=user(file=$(TREE))'+'\n'
out = out+'model1=YN98(frequencies='+get_mapnh_thetas(model=value_model, thetas=thetas)+',kappa='+str(value_kappa)+',omega='+str(value_omega)+')\n'
out = out+'root_freq1='+get_mapnh_thetas(model=value_model, thetas=root_thetas)+'\n'
out = out+'rate_distribution1=Gamma(n='+num_gamma_category+',alpha='+str(value_alpha)+',Gamma.beta='+str(value_alpha)+')\n'
out = out+'process1=Homogeneous(model=1, rate=1, tree=1, root_freq=1)'+'\n'
out = out+'phylo1=Single(process=1, data=1, useLog=yes)'+'\n'
out = out+'nullProcessParams=YN98.omega*=1'+'\n'
out = out+'map.type=DnDs'+'\n'
out = out+'output.counts=PerBranch(prefix=$(OUT).)'+'\n'    
print(out)

with open('iqtree2mapnh.params', "w") as file:
    file.write(out)
iqtree_tree.write(outfile='iqtree2mapnh.nwk', format=5)
print('iqtree2mapnh, done!')

