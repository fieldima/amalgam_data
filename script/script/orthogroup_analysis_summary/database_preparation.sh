#! /bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -pe def_slot 1
#$ -l s_vmem=32G
#$ -l mem_req=32G
##$ -l epyc
#$ -l d_rt=62:00:00:00
#$ -l s_rt=62:00:00:00
#$ -t 1

source /lustre6/home/lustre1/kfuku/.bashrc
ulimit -s unlimited

echo running on `hostname`
echo "`date`: Starting"

transcriptome_curation='sva'
l1ou_criterion='AICc'
l1ou_alpha_upper='PhylogeneticEM'

dir_script="/lustre6/home/lustre1/kfuku/my_script"
dir_ensembl="/lustre6/home/lustre1/kfuku/my_db/Ensembl/release-91"
dir_og="${dir_ensembl}/orthogroup"
dir_og_stat_branch="${dir_og}/stat.branch.${transcriptome_curation}.${l1ou_criterion}.${l1ou_alpha_upper}"
dir_og_stat_tree="${dir_og}/stat.tree.${transcriptome_curation}.${l1ou_criterion}.${l1ou_alpha_upper}"

python ${dir_script}/generate_orthogroup_database.py \
--overwrite 1 \
--dbpath ${dir_og}/Ensembl.91.orthogroup.db \
--dir_stat_tree ${dir_og_stat_tree} \
--dir_stat_branch ${dir_og_stat_branch}

###################
echo "`date`: Ending"


: <<'#_______________CO_______________'
#_______________CO_______________



