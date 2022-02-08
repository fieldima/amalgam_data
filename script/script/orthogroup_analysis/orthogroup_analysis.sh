#! /bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -pe def_slot 16
#$ -l s_vmem=8G
#$ -l mem_req=8G
##$ -l epyc
#$ -l d_rt=62:00:00:00
#$ -l s_rt=62:00:00:00
#$ -t 1-10

# all: 1-17896
# >3 genes: 1-15610
# 2333 for PGK

source /home/kfuku/.bashrc
ulimit -s unlimited

echo 'Python executable:' `which python`
echo running on `hostname`
echo "`date`: Starting"

MEM_PER_SLOT=`qstat -f -j ${JOB_ID} | grep ",mem_req=" | sed -e "s/.*mem_req=//" -e "s/G,.*//" | sed -e "s/G.*//"`
if [ -z ${MEM_PER_SLOT} ]; then
	MEM_PER_SLOT=3
fi
MEM_PER_HOST=$[${MEM_PER_SLOT}*${NSLOTS}]
echo "MEM_PER_HOST: ${MEM_PER_HOST}"
echo "SGE_TASK_ID: ${SGE_TASK_ID}"
cpu_id=`python -c 'import sys; from numpy import random; a = random.choice(range(64), int(sys.argv[1]), replace=False); print(",".join([str(b) for b in a]))' ${NSLOTS}`
echo CPU IDs for l1ou: ${cpu_id}

date_orthofinder='Mar22'
date_sra='2018_5_1'
transcriptome_curation='sva' # sva, tmmspe_sva, tmmall_sva
max_age=1105
similarity_method='pearson'
similarity_threshold=0.5
l1ou_criterion='AICc' #  "pBIC", "mBIC", "BIC", "AICc"
l1ou_nbootstrap=0
l1ou_use_fit_ind_file=1
l1ou_alpha_upper='PhylogeneticEM'
l1ou_convergence=0
very_large_max_nshift=10
phylogeneticem_use_fit_file=1
mapdnds_codon_freq='F3X4'

intron_gain_rate=0.0001
retrotransposition_rate=0.001 
chr_loc_transition_rate=0.001 # 'empirical' or float
all_params=${l1ou_criterion}_${l1ou_alpha_upper}_${mapdnds_codon_freq}_${intron_gain_rate}_${retrotransposition_rate}_${chr_loc_transition_rate}

exit_if_running=1
if [ ${exit_if_running} -eq 1 ]; then
	delete_preexisting_tmp_dir=1
else
	delete_preexisting_tmp_dir=0
fi
delete_tmp_dir=1

check_fasta_consistency=0
check_dated_tree=1
run_get_fasta=0
run_rps_blast=0
run_mafft=0
run_amas_original=0
run_maxalign=0
run_trimal=0
run_amas_cleaned=0
run_iqtree=0
run_notung_root=0
run_root=0
run_notung_reconcil=0
run_tree_dating=1
run_mapdnds=0
run_hyphy_dnds=0
run_expression=0
run_expression_fpkm=1
run_intron_chromosome_matrix=0
run_scm_intron=1
run_scm_chromosome=1
run_pem_original=0
run_l1ou_original=0
run_pem_fpkm=0
run_l1ou_fpkm=0
run_summary=1
run_tree_plot=0

dir_ensembl="/lustre6/home/lustre1/kfuku/my_db/Ensembl/release-91"
dir_og="${dir_ensembl}/orthogroup"
dir_gff3="${dir_ensembl}/gff3"
dir_transcriptome="${dir_ensembl}/curated_transcriptome/${date_sra}/tpm/tissue_mean"
dir_transcriptome_fpkm="${dir_ensembl}/curated_transcriptome/${date_sra}/tmm_rpkm/tissue_mean"
dir_phylogears="/lustre6/home/lustre1/kfuku/phylogears2-2.0.2016.09.06/bin"
dir_iqtree="/lustre6/home/lustre1/kfuku/iqtree-1.6.5-Linux/bin"
dir_myscript="/lustre6/home/lustre1/kfuku/my_script"
dir_myrepo="/lustre7/home/lustre4/kfuku/my_repository"
dir_amas="/lustre6/home/lustre1/kfuku/AMAS-master/amas"
dir_bpp="/lustre6/home/lustre1/kfuku/bpp/bin"
dir_hyphy="/lustre6/home/lustre1/kfuku/hyphy_2.3.11"
dir_hyphy_out="/lustre6/home/lustre1/kfuku/hyphy_2.3.11/lib/hyphy/TemplateBatchFiles"
file_genecount="${dir_ensembl}/orthofinder/Results_${date_orthofinder}/Orthogroups.GeneCount.csv"
file_og="${dir_ensembl}/orthofinder/Results_${date_orthofinder}/Orthogroups.txt"
species_tree="/lustre6/home/lustre1/kfuku/my_db/Ensembl/release-91/timetree/species_timetree.nwk"

db_cdd="/lustre6/home/lustre1/kfuku/my_db/CDD_NCBI/release_20170327/Cdd_NCBI"

ind=$[${SGE_TASK_ID}-1]
og_id=`python -c "import sys,pandas; df=pandas.read_csv(sys.argv[1],sep='\t',header=0); print(df.loc[int(sys.argv[2]),:].iloc[0])" ${file_genecount} ${ind}`
echo "OrthoGroup ID: ${og_id}"

dir_tmp_main="${dir_og}/tmp"
dir_tmp="${dir_tmp_main}/${SGE_TASK_ID}_${og_id}_${transcriptome_curation}" #_${RANDOM}
dir_cds_mask_longest="${dir_ensembl}/cds.mask.longest"
dir_og_cds_fasta="${dir_og}/cds.fasta"
dir_og_rpsblast="${dir_og}/rpsblast"
dir_og_mafft="${dir_og}/mafft"
dir_og_maxalign="${dir_og}/maxalign"
dir_og_trimal="${dir_og}/trimal"
dir_og_iqtree_tree="${dir_og}/iqtree.tree"
dir_og_iqtree_model="${dir_og}/iqtree.model"
dir_og_notung_root="${dir_og}/notung.root"
dir_og_root="${dir_og}/rooted_tree"
dir_og_root_log="${dir_og}/rooted_tree.log"
dir_og_notung_reconcil="${dir_og}/notung.reconcil"
dir_og_dated_tree="${dir_og}/dated_tree"
dir_og_dated_tree_log="${dir_og}/dated_tree.log"
#dir_og_dated_tree_log2="${dir_og}/dated_tree.log2"
dir_og_mapdnds_dn="${dir_og}/mapdNdS.dN.tree"
dir_og_mapdnds_ds="${dir_og}/mapdNdS.dS.tree"
dir_og_mapdnds_log="${dir_og}/mapdNdS.log"
dir_og_hyphy_dnds="${dir_og}/hyphy.dnds"
dir_og_expression="${dir_og}/character.expression.${transcriptome_curation}"
dir_og_expression_fpkm="${dir_og}/character.expression.fpkm"
dir_og_matrix_intron="${dir_og}/character.intron"
dir_og_matrix_chromosome="${dir_og}/character.chromosome"
dir_og_scm_intron_summary="${dir_og}/scm.intron.summary"
dir_og_scm_intron_plot="${dir_og}/scm.intron.plot"
dir_og_scm_chromosome_summary="${dir_og}/scm.chromosome.summary"
dir_og_scm_chromosome_plot="${dir_og}/scm.chromosome.plot"
# TPM
dir_og_pem_rdata="${dir_og}/PhylogeneticEM.${transcriptome_curation}.rdata"
dir_og_pem_tree="${dir_og}/PhylogeneticEM.${transcriptome_curation}.tree"
dir_og_pem_regime="${dir_og}/PhylogeneticEM.${transcriptome_curation}.regime"
dir_og_pem_leaf="${dir_og}/PhylogeneticEM.${transcriptome_curation}.leaf"
dir_og_pem_plot="${dir_og}/PhylogeneticEM.${transcriptome_curation}.plot"
dir_og_l1ou_fit_rdata="${dir_og}/l1ou.fit.rdata.${transcriptome_curation}.${l1ou_criterion}.${l1ou_alpha_upper}"
dir_og_l1ou_fit_tree="${dir_og}/l1ou.fit.tree.${transcriptome_curation}.${l1ou_criterion}.${l1ou_alpha_upper}"
dir_og_l1ou_fit_regime="${dir_og}/l1ou.fit.regime.${transcriptome_curation}.${l1ou_criterion}.${l1ou_alpha_upper}"
dir_og_l1ou_fit_leaf="${dir_og}/l1ou.fit.leaf.${transcriptome_curation}.${l1ou_criterion}.${l1ou_alpha_upper}"
dir_og_l1ou_fit_plot="${dir_og}/l1ou.fit.plot.${transcriptome_curation}.${l1ou_criterion}.${l1ou_alpha_upper}"
# TMM-FPKM
dir_og_pem_rdata_fpkm="${dir_og}/PhylogeneticEM.${transcriptome_curation}.rdata.fpkm"
dir_og_pem_tree_fpkm="${dir_og}/PhylogeneticEM.${transcriptome_curation}.tree.fpkm"
dir_og_pem_regime_fpkm="${dir_og}/PhylogeneticEM.${transcriptome_curation}.regime.fpkm"
dir_og_pem_leaf_fpkm="${dir_og}/PhylogeneticEM.${transcriptome_curation}.leaf.fpkm"
dir_og_pem_plot_fpkm="${dir_og}/PhylogeneticEM.${transcriptome_curation}.plot.fpkm"
dir_og_l1ou_fit_rdata_fpkm="${dir_og}/l1ou.fit.rdata.${transcriptome_curation}.${l1ou_criterion}.${l1ou_alpha_upper}.fpkm"
dir_og_l1ou_fit_tree_fpkm="${dir_og}/l1ou.fit.tree.${transcriptome_curation}.${l1ou_criterion}.${l1ou_alpha_upper}.fpkm"
dir_og_l1ou_fit_regime_fpkm="${dir_og}/l1ou.fit.regime.${transcriptome_curation}.${l1ou_criterion}.${l1ou_alpha_upper}.fpkm"
dir_og_l1ou_fit_leaf_fpkm="${dir_og}/l1ou.fit.leaf.${transcriptome_curation}.${l1ou_criterion}.${l1ou_alpha_upper}.fpkm"
dir_og_l1ou_fit_plot_fpkm="${dir_og}/l1ou.fit.plot.${transcriptome_curation}.${l1ou_criterion}.${l1ou_alpha_upper}.fpkm"
# Summary
dir_og_amas_original="${dir_og}/amas.original"
dir_og_amas_cleaned="${dir_og}/amas.cleaned"
dir_og_stat_branch="${dir_og}/stat.branch.${transcriptome_curation}.${l1ou_criterion}.${l1ou_alpha_upper}"
dir_og_stat_tree="${dir_og}/stat.tree.${transcriptome_curation}.${l1ou_criterion}.${l1ou_alpha_upper}"
dir_og_tree_plot="${dir_og}/tree_plot.${transcriptome_curation}.${l1ou_criterion}.${l1ou_alpha_upper}"

file_og_cds_fasta="${dir_og_cds_fasta}/${og_id}.cds.fasta"
file_og_rpsblast="${dir_og_rpsblast}/${og_id}.rpsblast.xml"
file_og_mafft="${dir_og_mafft}/${og_id}.cds.aln.fasta"
file_og_maxalign="${dir_og_maxalign}/${og_id}.cds.maxalign.fasta"
file_og_trimal="${dir_og_trimal}/${og_id}.cds.trimal.fasta"
file_og_iqtree_tree="${dir_og_iqtree_tree}/${og_id}.iqtree.nwk"
file_og_iqtree_model="${dir_og_iqtree_model}/${og_id}.model.gz"
file_og_notung_root="${dir_og_notung_root}/${og_id}.notung.root.zip"
file_og_root="${dir_og_root}/${og_id}.root.nwk"
file_og_root_log="${dir_og_root}/${og_id}.root.txt"
file_og_notung_reconcil="${dir_og_notung_reconcil}/${og_id}.notung.reconcil.zip"
file_og_dated_tree="${dir_og_dated_tree}/${og_id}.dated.nwk"
file_og_dated_tree_log="${dir_og_dated_tree_log}/${og_id}.dated.log.txt"
file_og_mapdnds_dn="${dir_og_mapdnds_dn}/${og_id}.mapdNdS.dN.nwk"
file_og_mapdnds_ds="${dir_og_mapdnds_ds}/${og_id}.mapdNdS.dS.nwk"
file_og_mapdnds_log="${dir_og_mapdnds_log}/${og_id}.mapdNdS.log.zip"
file_og_hyphy_dnds="${dir_og_hyphy_dnds}/${og_id}.hyphy.dnds.nwk"
file_og_expression="${dir_og_expression}/${og_id}.expression.tsv"
file_og_expression_fpkm="${dir_og_expression_fpkm}/${og_id}.expression.fpkm.tsv"
file_og_matrix_intron="${dir_og_matrix_intron}/${og_id}.intron.tsv"
file_og_matrix_chromosome="${dir_og_matrix_chromosome}/${og_id}.chromosome.tsv"
file_og_scm_intron_summary="${dir_og_scm_intron_summary}/${og_id}.scm.intron.tsv"
file_og_scm_intron_plot="${dir_og_scm_intron_plot}/${og_id}.scm.intron.pdf"
file_og_scm_chromosome_summary="${dir_og_scm_chromosome_summary}/${og_id}.scm.chromosome.tsv"
file_og_scm_chromosome_plot="${dir_og_scm_chromosome_plot}/${og_id}.scm.chromosome.pdf"
# TPM
file_og_pem_rdata="${dir_og_pem_rdata}/${og_id}.PhylogeneticEM.RData"
file_og_pem_tree="${dir_og_pem_tree}/${og_id}.PhylogeneticEM.tree.tsv"
file_og_pem_regime="${dir_og_pem_regime}/${og_id}.PhylogeneticEM.regime.tsv"
file_og_pem_leaf="${dir_og_pem_leaf}/${og_id}.PhylogeneticEM.leaf.tsv"
file_og_pem_plot="${dir_og_pem_plot}/${og_id}.PhylogeneticEM.pdf"
file_og_l1ou_fit_rdata="${dir_og_l1ou_fit_rdata}/${og_id}.l1ou.RData"
file_og_l1ou_fit_tree="${dir_og_l1ou_fit_tree}/${og_id}.l1ou.tree.tsv"
file_og_l1ou_fit_regime="${dir_og_l1ou_fit_regime}/${og_id}.l1ou.regime.tsv"
file_og_l1ou_fit_leaf="${dir_og_l1ou_fit_leaf}/${og_id}.l1ou.leaf.tsv"
file_og_l1ou_fit_plot="${dir_og_l1ou_fit_plot}/${og_id}.l1ou.pdf"
# FPKM
file_og_pem_rdata_fpkm="${dir_og_pem_rdata_fpkm}/${og_id}.fpkm.PhylogeneticEM.RData"
file_og_pem_tree_fpkm="${dir_og_pem_tree_fpkm}/${og_id}.fpkm.PhylogeneticEM.tree.tsv"
file_og_pem_regime_fpkm="${dir_og_pem_regime_fpkm}/${og_id}.fpkm.PhylogeneticEM.regime.tsv"
file_og_pem_leaf_fpkm="${dir_og_pem_leaf_fpkm}/${og_id}.fpkm.PhylogeneticEM.leaf.tsv"
file_og_pem_plot_fpkm="${dir_og_pem_plot_fpkm}/${og_id}.fpkm.PhylogeneticEM.pdf"
file_og_l1ou_fit_rdata_fpkm="${dir_og_l1ou_fit_rdata_fpkm}/${og_id}.fpkm.l1ou.RData"
file_og_l1ou_fit_tree_fpkm="${dir_og_l1ou_fit_tree_fpkm}/${og_id}.fpkm.l1ou.tree.tsv"
file_og_l1ou_fit_regime_fpkm="${dir_og_l1ou_fit_regime_fpkm}/${og_id}.fpkm.l1ou.regime.tsv"
file_og_l1ou_fit_leaf_fpkm="${dir_og_l1ou_fit_leaf_fpkm}/${og_id}.fpkm.l1ou.leaf.tsv"
file_og_l1ou_fit_plot_fpkm="${dir_og_l1ou_fit_plot_fpkm}/${og_id}.fpkm.l1ou.pdf"
# Summary
file_og_stat_branch="${dir_og_stat_branch}/${og_id}.stat.branch.tsv"
file_og_stat_tree="${dir_og_stat_tree}/${og_id}.stat.tree.tsv"
file_og_amas_original="${dir_og_amas_original}/${og_id}.amas.original.tsv"
file_og_amas_cleaned="${dir_og_amas_cleaned}/${og_id}.amas.cleaned.tsv"
file_og_tree_plot="${dir_og_tree_plot}/${og_id}.tree_plot.pdf"

if [ ${exit_if_running} -eq 1 ]; then
	running_id=( `qstat | grep -v -e " dr " -e "QLOGIN" -e "${JOB_ID} .* ${SGE_TASK_ID}$" | sed -e "s/ qw / qw dummy/" -e "s/\s\+/\t/g" -e "1,2d" | cut -f11 | sort -u` )
	echo "running_id: ${running_id[*]}"
	flag=1
	for r in ${running_id[*]}; do
		if [ `echo ${r} | grep -v ":"` ]; then
			if [ ${r} -eq ${SGE_TASK_ID} ]; then
				flag=0
			fi
		fi
	done
	if [ ${flag} -eq 0 ]; then
		echo "Exit, SGE_TASK_ID=${SGE_TASK_ID}, JOB_ID=${JOB_ID} is running already."
		exit
	fi
fi

directories=( `set | grep "^dir_" | sed -e "s/=.*//"` )
for d in ${directories[@]}; do
	if [ ! -e `eval echo '$'${d}` ] && \
	[ "${d}" != "dir_ensembl" ] && \
	[ "${d}" != "dir_work" ] && \
	[ "${d}" != "dir_hyphy_out" ] && \
	[ "${d}" != "dir_transcriptome" ]; then
		echo creating: `eval echo '$'${d}`
		mkdir `eval echo '$'${d}`
	fi
done

if [ -e ${dir_tmp} ] && [ ${delete_preexisting_tmp_dir} -eq 1 ]; then
	echo "`date`: Deleting preexisting ${dir_tmp}"
	rm -r ${dir_tmp_main}/${SGE_TASK_ID}_*
fi
if [ ! -e ${dir_tmp} ]; then
	echo "Making ${dir_tmp}"
	mkdir ${dir_tmp}
fi
cd ${dir_tmp}
echo "Working in `pwd`"

task="fasta consistency check"
if [ ${check_fasta_consistency} -eq 1 ] && [ -s ${dir_og_cds_fasta}/${og_id}.cds.fasta ]; then
	echo "`date`: Start: ${task}"
	num_gene=`python -c "import sys,pandas; df=pandas.read_csv(sys.argv[1],sep='\t',header=0); print(df.loc[int(sys.argv[2]),'Total'])" ${file_genecount} ${ind}`
	num_seq=`grep -e "^>" ${dir_og_cds_fasta}/${og_id}.cds.fasta | wc -l`
	echo "Number of genes in ${og_id}: ${num_gene}"
	echo "Number of sequences in the fasta: ${num_seq}"
	if [ ${num_gene} -ne ${num_seq} ]; then
		echo "Number of genes and sequences did not match. Deleting existing outputs..."
		directories=( `set | grep "^dir_og_" | sed -e "s/=.*//"` )
		for d in ${directories[@]}; do
			if [ ! -e `eval echo '$'${d}` ] && [ -e ${d}/*${og_id}* ] && [ -n "${og_id}" ] && [ -n "`echo ${og_id} | grep -e OG`" ]; then
				echo 'Deleting the files:'
				echo `ls ${d}/*${og_id}*`
				rm ${d}/*${og_id}*
			fi
		done
	else
		echo "Passed the fasta consistency check: ${dir_og_cds_fasta}/${og_id}.cds.fasta"
	fi
else
	echo "`date`: Skipped: ${task}"
fi

task="check dated tree"
if [ ${check_dated_tree} -eq 1 ] && [ -s ${file_og_dated_tree} ]; then
	echo "`date`: Start: ${task}"
	contain_negative_bl=` grep ":-" ${file_og_dated_tree} | wc -l`
	if [ ${contain_negative_bl} -eq 1 ]; then
		echo "Dated tree has negative branch length. Deleting output files depending on the tree file."
		for key in l1ou pem scm dated stat tree_plot; do
			files=( `set | grep "^file_og_${key}" | sed -e "s/=.*//"` )
			for f in ${files[@]}; do
				if [ -e `eval echo '$'${f}` ]; then
					echo deleting: `eval echo '$'${f}`
					rm `eval echo '$'${f}`
				fi
			done
		done
	else
		echo "Dated tree has no negative branch length. Continue."
	fi
else
	echo "`date`: Skipped: ${task}"
fi

task="fasta generation"
if [ ! -s ${file_og_cds_fasta} ] && [ ${run_get_fasta} -eq 1 ]; then
	echo "`date`: Start: ${task}"
	touch ${og_id}.cds.fasta

	genes=( `cat ${file_og} | grep -e "${og_id}:" | sed -e "s/${og_id}:[[:space:]]//"` )
	spp=( `echo ${genes[@]} | sed -e "s/[[:space:]]/\n/g" | sed -e "s/_/-/" | sed -e "s/_.*//" | sed -e "s/-/_/" | sort -u` )
	echo "species: ${spp[@]}"
	for sp in ${spp[@]}; do
		sp_cds=`find ${dir_cds_mask_longest}/${sp}*`
		sp_genes=( `echo ${genes[@]} | sed -e "s/[[:space:]]/\n/g" | grep -e ${sp}` )
		echo "${sp_genes[@]}"
		echo "${sp_genes[@]}" | sed -e "s/[[:space:]]/\n/g" | fatt extract --stdin ${sp_cds} >> ${og_id}.cds.fasta
	done

	num_gene=`python -c "import sys,pandas; df=pandas.read_csv(sys.argv[1],sep='\t',header=0); print(df.loc[int(sys.argv[2]),'Total'])" ${file_genecount} ${ind}`
	num_seq=`grep -e "^>" ${og_id}.cds.fasta | wc -l`
	echo "Number of genes in ${og_id}: ${num_gene}"
	echo "Number of sequences in the fasta: ${num_seq}"
	if [ ${num_gene} -eq ${num_seq} ]; then
		echo "Number of genes and sequences matched. Fasta generation completed!"
		cp ${og_id}.cds.fasta ${file_og_cds_fasta}
	fi
else
	echo "`date`: Skipped: ${task}"
fi

task="RPS-BLAST"
if [ ! -s ${file_og_rpsblast} ] && [ ${run_rps_blast} -eq 1 ]; then
	echo "`date`: Start: ${task}"
	conda activate py27 > /dev/null

	seqmagick convert \
	--ungap \
	--translate dna2protein \
	--input-format fasta \
	--output-format fasta \
	${file_og_cds_fasta} \
	ungapped_translated_cds.fas

	# makeprofiledb -in ${db_cdd} -dbtype 'rps'

	rpsblast \
	-query ungapped_translated_cds.fas \
	-db ${db_cdd} \
	-out ${og_id}.rpsblast.xml \
	-evalue 0.01 \
	-outfmt 5 \
	-num_threads ${NSLOTS}

	cp ${og_id}.rpsblast.xml ${file_og_rpsblast}

	conda activate base > /dev/null
else
	echo "`date`: Skipped: ${task}"
fi

task="mafft"
if [ ! -s ${file_og_mafft} ] && [ ${run_mafft} -eq 1 ]; then
	echo "`date`: Start: ${task}"

	python ${dir_myscript}/inframe_Nstop2gap.py \
	${file_og_cds_fasta} \
	tmp.cds.fasta \
	"Standard"

	transeq \
	tmp.cds.fasta \
	tmp.pep.fasta

	/lustre6/home/lustre1/kfuku/mafft-7.394/bin/mafft \
	--auto \
	--amino \
	--thread ${NSLOTS} \
	tmp.pep.fasta \
	> tmp.pep.aln.fasta

	tranalign \
	-table 0 \
	-asequence tmp.cds.fasta \
	-bsequence tmp.pep.aln.fasta \
	-outseq tmp.cds.aln.fasta

	python ${dir_myscript}/pad_alignment.py \
	tmp.cds.aln.fasta \
	${og_id}.cds.aln.fasta

	cp ${og_id}.cds.aln.fasta ${file_og_mafft}
else
	echo "`date`: Skipped: ${task}"
fi

task="AMAS for original alignment"
if [ ! -s ${file_og_amas_original} ] && [ ${run_amas_original} -eq 1 ]; then
	echo "`date`: Start: ${task}"
	python ${dir_amas}/AMAS.py summary \
	--in-format fasta \
	--data-type dna \
	--in-files ${file_og_mafft}

	mv summary.txt ${file_og_amas_original}
else
	echo "`date`: Skipped: ${task}"
fi

task="MaxAlign"
if [ ! -s ${file_og_maxalign} ] && [ ${run_maxalign} -eq 1 ]; then
	echo "`date`: Start: ${task}"
	cp ${file_og_mafft} ${og_id}.cds.aln.fasta

	/lustre6/home/lustre1/kfuku/bin/maxalign.pl \
	-i=999 \
	-v=1 \
	-f=${og_id}.maxalign. \
	${og_id}.cds.aln.fasta

	echo Number of sequences before MaxAlign: `grep -e "^>" ${og_id}.cds.aln.fasta | wc -l`
	echo Number of sequences after MaxAlign: `grep -e "^>" ${og_id}.maxalign.heuristic.fsa | wc -l`

	cp ${og_id}.maxalign.heuristic.fsa ${file_og_maxalign}
else
	echo "`date`: Skipped: ${task}"
fi

task="TrimAl"
if [ ! -s ${file_og_trimal} ] && [ ${run_trimal} -eq 1 ]; then
	echo "`date`: Start: ${task}"

	${dir_phylogears}/pgconvseq \
	--output=NEXUS \
	${file_og_maxalign} \
	${og_id}.cds.maxalign.nex \
	> /dev/null 2>&1

	${dir_phylogears}/pgtrimal \
	--frame=1 \
	--method=GAPPYOUT \
	${og_id}.cds.maxalign.nex \
	${og_id}.cds.trimal.nex \
	> /dev/null 2>&1

	${dir_phylogears}/pgconvseq \
	--output=fasta \
	${og_id}.cds.trimal.nex \
	${og_id}.cds.trimal.fasta \
	> /dev/null 2>&1

	cp ${og_id}.cds.trimal.fasta ${file_og_trimal}
else
	echo "`date`: Skipped: ${task}"
fi

task="AMAS for cleaned alignment"
if [ ! -s ${file_og_amas_cleaned} ] && [ ${run_amas_cleaned} -eq 1 ]; then
	echo "`date`: Start: ${task}"
	python ${dir_amas}/AMAS.py summary \
	--in-format fasta \
	--data-type dna \
	--in-files ${file_og_trimal}

	mv summary.txt ${file_og_amas_cleaned}
else
	echo "`date`: Skipped: ${task}"
fi

num_gene_after_maxalign=`grep -e "^>" ${dir_og_maxalign}/${og_id}.cds.maxalign.fasta | wc -l`
if [ ${num_gene_after_maxalign} -lt 3 ]; then
	echo 'Number of genes after MaxAlign is not sufficient (<3) for tree-based analysis. Exitng.'
	exit
else
	echo "Number of genes after MaxAlign: ${num_gene_after_maxalign}"
fi

task="IQ-TREE"
if [[ ( ! -s ${file_og_iqtree_tree} || ! -s ${file_og_iqtree_model} ) && ${run_iqtree} -eq 1 ]]; then
	echo "`date`: Start: ${task}"
	num_seq=`grep -e "^>" ${dir_og_trimal}/${og_id}.cds.trimal.fasta | wc -l`
	if [ ${num_seq} -ge 4 ]; then
		bootstrap_params="-bb 1000 -bnni"
		file_tree="${og_id}.contree"
	else
		bootstrap_params=""
		file_tree="${og_id}.treefile"
	fi

	${dir_iqtree}/iqtree \
	-s ${file_og_trimal} \
	-st DNA \
	-m MFP \
	-pre ${og_id} \
	-nt ${NSLOTS} \
	-mem ${MEM_PER_HOST}G \
	-seed 12345 \
	${bootstrap_params}

	cp ${file_tree} ${file_og_iqtree_tree}
	cp ${og_id}.model.gz ${file_og_iqtree_model}
else
	echo "`date`: Skipped: ${task}"
fi

task="NOTUNG rooting"
if [ ! -s ${file_og_notung_root} ] && [ ${run_notung_root} -eq 1 ]; then
	echo "`date`: Start: ${task}"

	if [ -e ./${og_id}.notung.root ]; then
		rm -r ${og_id}.notung.root
	fi

	memory_notung=$[${MEM_PER_HOST}/2]
	echo "memory_notung: ${memory_notung}"

	java -jar -Xmx${memory_notung}g /lustre6/home/lustre1/kfuku/Notung-2.9/Notung-2.9.jar \
	-s ${species_tree} \
	-g ${file_og_iqtree_tree} \
	--root \
	--infertransfers "false" \
	--treeoutput newick \
	--log \
	--treestats \
	--events \
	--parsable \
	--speciestag prefix \
	--allopt \
	--maxtrees 1000 \
	--nolosses \
	--outputdir ./${og_id}.notung.root

	if [ `ls ${og_id}.notung.root/${og_id}.iqtree.nwk.rooting.0 | wc -l` -gt 0 ]; then
		zip -rq ${og_id}.notung.root.zip ${og_id}.notung.root
		cp ${og_id}.notung.root.zip ${file_og_notung_root}
	fi
else
	echo "`date`: Skipped: ${task}"
fi

task="Tree rooting"
if [[ ( ! -s ${file_og_root} || ! -s ${file_og_root_log} ) && ${run_root} -eq 1 ]]; then
	echo "`date`: Start: ${task}"
	if [ -e ./${og_id}.notung.root ]; then
		rm -r ./${og_id}.notung.root
	fi
	cp ${dir_og_notung_root}/${og_id}.notung.root.zip .
	unzip -q ${og_id}.notung.root.zip

	Rscript ${dir_myscript}/gene_tree_rooting.r \
	./${og_id}.notung.root/ \
	${file_og_iqtree_tree} \
	${og_id}.root.nwk \
	${NSLOTS} \
	2>&1 | tee ${og_id}.root.txt

	if [ -s ${og_id}.root.nwk ]; then
		cp ${og_id}.root.txt ${file_og_root_log}
		cp ${og_id}.root.nwk ${file_og_root}

	fi
else
	echo "`date`: Skipped: ${task}"
fi

task="NOTUNG reconciliation"
if [ ! -s ${file_og_notung_reconcil} ] && [ ${run_notung_reconcil} -eq 1 ]; then
	echo "`date`: Start: ${task}"

	memory_notung=$[${MEM_PER_HOST}/2]
	echo "memory_notung: ${memory_notung}"

	if [ -s ./${og_id}.root.nwk ]; then
		rm ${og_id}.root.nwk
	fi
	if [ -e ./${og_id}.notung.reconcil ]; then
		rm -r ${og_id}.notung.reconcil
	fi
	python ${dir_myscript}/delete_internal_node_name_and_support.py \
	${file_og_root} \
	${og_id}.root.nwk

	java -jar -Xmx${memory_notung}g /lustre6/home/lustre1/kfuku/Notung-2.9/Notung-2.9.jar \
	-s ${species_tree} \
	-g ${og_id}.root.nwk \
	--reconcile \
	--infertransfers "false" \
	--treeoutput newick \
	--log \
	--treestats \
	--events \
	--parsable \
	--speciestag prefix \
	--maxtrees 1 \
	--nolosses \
	--outputdir ./${og_id}.notung.reconcil

	if [ -s ${og_id}.notung.reconcil/${og_id}.root.nwk.reconciled.*parsable.txt ]; then
		zip -rq ${og_id}.notung.reconcil.zip ${og_id}.notung.reconcil
		cp ${og_id}.notung.reconcil.zip ${file_og_notung_reconcil}
	fi
else
	echo "`date`: Skipped: ${task}"
fi

task="Tree dating"
if [[ ( ! -s ${file_og_dated_tree} || ! -s ${file_og_dated_tree_log} ) && ${run_tree_dating} -eq 1 ]]; then
	echo "`date`: Start: ${task}"
	if [ -e ./${og_id}.notung.reconcil ]; then
		rm -r ./${og_id}.notung.reconcil
	fi
	cp ${dir_og_notung_reconcil}/${og_id}.notung.reconcil.zip .
	unzip -q ${og_id}.notung.reconcil.zip
	if [ -s ./${og_id}.notung.reconcil/${og_id}.root.nwk.reconciled.0 ]; then
		cp ./${og_id}.notung.reconcil/${og_id}.root.nwk.reconciled.0 ./${og_id}.notung.reconcil/${og_id}.root.nwk.reconciled
		cp ./${og_id}.notung.reconcil/${og_id}.root.nwk.reconciled.0.parsable.txt ./${og_id}.notung.reconcil/${og_id}.root.nwk.reconciled.parsable.txt
	fi

	${dir_myrepo}/RADTE/radte \
	--species_tree=${species_tree} \
	--gene_tree=./${og_id}.notung.reconcil/${og_id}.root.nwk.reconciled \
	--notung_parsable=./${og_id}.notung.reconcil/${og_id}.root.nwk.reconciled.parsable.txt \
	--max_age=${max_age} \
	--chronos_lambda=1 \
	--chronos_model=discrete \
	--pad_short_edge=0 \
	2>&1 | tee radte.log

	constrained_node=`grep -e "^Calibrated nodes: " radte.log | sed -e "s/Calibrated nodes: //" -e "s/[[:space:]]//g"`
	echo ${constrained_node} > ${og_id}.dated.log.txt
	if [ ! "`echo ${og_id}.dated.log.txt`" = "" ]; then
		cp radte_calibrated_nodes.txt ${file_og_dated_tree_log}
		cp radte_gene_tree_output.nwk ${file_og_dated_tree}
	fi
else
	echo "`date`: Skipped: ${task}"
fi

task="mapdNdS"
if [[ ( ! -s ${file_og_mapdnds_dn} || ! -s ${file_og_mapdnds_ds} || ! -s ${file_og_mapdnds_log} ) && ${run_mapdnds} -eq 1 ]]; then
	echo "`date`: Start: ${task}"

	if [ -e ${dir_tmp}/${og_id}.mapdNdS ]; then
		rm -r ${dir_tmp}/${og_id}.mapdNdS
	fi
	if [ -e ${dir_tmp}/${og_id}.mapdNdS.log.zip ]; then
		rm -r ${dir_tmp}/${og_id}.mapdNdS.log.zip
	fi
	mkdir ${dir_tmp}/${og_id}.mapdNdS
	cd ${dir_tmp}/${og_id}.mapdNdS


	echo "codon_freq: ${mapdnds_codon_freq}"

	${dir_iqtree}/iqtree \
	-s ${file_og_trimal} \
	-st CODON \
	-m GY+${mapdnds_codon_freq}+G4 \
	-pre ${og_id}.iqtree2mapdNdS \
	-nt ${NSLOTS} \
	-mem ${MEM_PER_HOST}G \
	-te ${file_og_root} \
	-asr \
	-seed 12345

	python ${dir_myscript}/iqtree2mapnh.py \
	--iqtree ${og_id}.iqtree2mapdNdS.iqtree \
	--log ${og_id}.iqtree2mapdNdS.log \
	--state ${og_id}.iqtree2mapdNdS.state \
	--alignment ${file_og_trimal} \
	--treefile ${og_id}.iqtree2mapdNdS.treefile \
	--rooted_tree ${dir_og_root}/${og_id}.root.nwk

	LD_LIBRARY_PATH="/lustre6/home/lustre1/kfuku/bpp/lib64/:/lustre6/home/lustre1/kfuku/.pyenv/versions/anaconda3-4.1.0/envs/bpp/lib/:${LD_LIBRARY_PATH}"

	${dir_bpp}/mapnh \
	SEQ=${file_og_trimal} \
	TREE=iqtree2mapnh.nwk \
	OUT=${og_id} \
	param=iqtree2mapnh.params \
	2>&1 | tee mapnh.log.txt

	mv ${og_id}.dN.dnd ${file_og_mapdnds_dn}
	mv ${og_id}.dS.dnd ${file_og_mapdnds_ds}

	cd ${dir_tmp}

	if [ -s ${out1} ] && [ -s ${out2} ]; then
		zip -rq ${og_id}.mapdNdS.log.zip ${og_id}.mapdNdS
		cp ${og_id}.mapdNdS.log.zip ${file_og_mapdnds_log}
	fi
else
	echo "`date`: Skipped: ${task}"
fi

task="Hyphy dN-dS estimation"
if [ ! -s ${file_og_hyphy_dnds} ] && [ ${run_hyphy_dnds} -eq 1 ]; then
	echo "`date`: Start: ${task}"

	LD_LIBRARY_PATH="/lustre6/home/lustre1/kfuku/.pyenv/versions/anaconda3-4.1.0/envs/bpp/lib/:${LD_LIBRARY_PATH}"

	command_hyphy="${dir_hyphy}/bin/HYPHYMP -p ${dir_hyphy}/lib/hyphy/TemplateBatchFiles/AnalyzeCodonData.bf"
	(echo 1; echo "${file_og_trimal}"; echo "MG94W9"; echo 1; echo "${file_og_root}"; echo "y"; echo 9; echo 2; echo "${og_id}.hyphy.dnds.nwk") | ${command_hyphy}

	mv ${dir_hyphy_out}/${og_id}.hyphy.dnds.nwk .
	cp ${og_id}.hyphy.dnds.nwk ${file_og_hyphy_dnds}

	# the output tree order in ${og_id}.hyphy.dnds.nwk
	# 1. E[Syn subs/nucleotide site] tree
	# 2. E[Non-syn subs/nucleotide site] tree
	# 3. dS tree
	# 4. dN tree
else
	echo "`date`: Skipped: ${task}"
fi

task="Expression matrix preparation, TPM"
if [ ! -s ${file_og_expression} ] && [ ${run_expression} -eq 1 ]; then
	echo "`date`: Start: ${task}"

	Rscript ${dir_myscript}/get_expression_matrix.r \
	${file_og_root} \
	${dir_transcriptome}/

	mv expression_matrix.tsv ${og_id}.expression.tsv
	cp ${og_id}.expression.tsv ${file_og_expression}

else
	echo "`date`: Skipped: ${task}"
fi

task="Expression matrix preparation, FPKM"
if [ ! -s ${file_og_expression_fpkm} ] && [ ${run_expression_fpkm} -eq 1 ]; then
	echo "`date`: Start: ${task}"

	Rscript ${dir_myscript}/get_expression_matrix.r \
	${file_og_root} \
	${dir_transcriptome_fpkm}/

	mv expression_matrix.tsv ${og_id}.expression.fpkm.tsv
	cp ${og_id}.expression.fpkm.tsv ${file_og_expression_fpkm}

else
	echo "`date`: Skipped: ${task}"
fi

task="intron and chromosome matrix preparation"
if [[ ( ! -s ${file_og_matrix_intron} || ! -s ${file_og_matrix_chromosome} ) && ${run_intron_chromosome_matrix} -eq 1 ]]; then
	echo "`date`: Start: ${task}"

	# cd ${dir_gff3}
	# wget -r -c -nd -np -A .gff3.gz ftp://ftp.ensembl.org/pub/release-91/gff3/

	Rscript ${dir_myscript}/get_intron_chromosome_matrix.r \
	${file_og_dated_tree} \
	${dir_gff3}

	cp intron_matrix.tsv ${file_og_matrix_intron}
	cp chromosome_matrix.tsv ${file_og_matrix_chromosome}
else
	echo "`date`: Skipped: ${task}"
fi

task="stochastic character mapping of intron evolution"
if [ ! -s ${file_og_scm_intron_summary} ] && [ ${run_scm_intron} -eq 1 ]; then
	echo "`date`: Start: ${task}"

	Rscript ${dir_myscript}/scm_intron_evolution.r \
	--tree_file=${file_og_dated_tree} \
	--trait_file=${file_og_matrix_intron} \
	--intron_gain_rate=${intron_gain_rate} \
	--retrotransposition_rate=${retrotransposition_rate} \
	--nrep=1000 \
	--nslots=${NSLOTS}

	cp intron_evolution_summary.tsv ${file_og_scm_intron_summary}
	if [ -e intron_evolution_plot.pdf ]; then
		cp intron_evolution_plot.pdf ${file_og_scm_intron_plot}
	fi
else
	echo "`date`: Skipped: ${task}"
fi

task="stochastic character mapping of chromosome evolution"
if [ ! -s ${file_og_scm_chromosome_summary} ] && [ ${run_scm_chromosome} -eq 1 ]; then
	echo "`date`: Start: ${task}"

	# Chinchilla lanigera is excluded because of the fragmented genome unanchored to chromosomes
	# Other species are excluded because their sex chromosomes are non-homologous to mammals

	Rscript ${dir_myscript}/scm_chromosome_location_evolution.r \
	--tree_file=${file_og_dated_tree} \
	--trait_file=${file_og_matrix_chromosome} \
	--chr_loc_transition_rate=${chr_loc_transition_rate} \
	--omit_species_prefix='Astyanax_mexicanus|Danio_rerio|Gadus_morhua|Oryzias_latipes|Oreochromis_niloticus|Xenopus_tropicalis|Gallus_gallus|Anolis_carolinensis|Ornithorhynchus_anatinus|Chinchilla_lanigera' \
	--omit_chromosome_prefix='Un' \
	--max_chromosome_nchar=4 \
	--nrep=1000 \
	--nslots=${NSLOTS}

	cp chromosome_location_evolution_summary.tsv ${file_og_scm_chromosome_summary}
	if [ -e chromosome_location_evolution_plot.pdf ]; then
		cp chromosome_location_evolution_plot.pdf ${file_og_scm_chromosome_plot}
	fi
else
	echo "`date`: Skipped: ${task}"
fi

task="PhylogeneticEM, TPM"
if [ ${num_gene_after_maxalign} -le 3 ]; then
	echo "Skipped Phylogenetic EM because num_gene_after_maxalign = ${num_gene_after_maxalign}"
fi
if [[ ( ! -s ${file_og_pem_rdata} || ! -s ${file_og_pem_tree} ||  ! -s ${file_og_pem_regime} || ! -s ${file_og_pem_leaf} ) && ${num_gene_after_maxalign} -gt 3 && ${run_pem_original} -eq 1 ]]; then
	echo "`date`: Start: ${task}"

	pem_fit_file=''
	if [ ${phylogeneticem_use_fit_file} -eq 1 ] && [ -s ${file_og_pem_rdata} ]; then
		pem_fit_file=${file_og_pem_rdata}
	fi

	Rscript ${dir_myscript}/detect_OU_shift_PhylogeneticEM.r \
	--tree_file=${file_og_dated_tree} \
	--trait_file=${file_og_expression} \
	--nslots=${NSLOTS} \
	--fit_file=${pem_fit_file} \
	--similarity_method=${similarity_method} \
	--similarity_threshold=0.99 \
	--ceil_negative=0

	mv PhylogeneticEM.RData ${file_og_pem_rdata}
	mv PhylogeneticEM.tree.tsv ${file_og_pem_tree}
	mv PhylogeneticEM.regime.tsv ${file_og_pem_regime}
	mv PhylogeneticEM.leaf.tsv ${file_og_pem_leaf}
	mv PhylogeneticEM.plot.pdf ${file_og_pem_plot}

else
	echo "`date`: Skipped: ${task}"
fi

task="PhylogeneticEM, FPKM"
if [ ${num_gene_after_maxalign} -le 3 ]; then
	echo "Skipped Phylogenetic EM because num_gene_after_maxalign = ${num_gene_after_maxalign}"
fi
if [[ ( ! -s ${file_og_pem_rdata_fpkm} || ! -s ${file_og_pem_tree_fpkm} ||  ! -s ${file_og_pem_regime_fpkm} || ! -s ${file_og_pem_leaf_fpkm} ) && ${num_gene_after_maxalign} -gt 3 && ${run_pem_fpkm} -eq 1 ]]; then
	echo "`date`: Start: ${task}"

	pem_fit_file=''
	if [ ${phylogeneticem_use_fit_file} -eq 1 ] && [ -s ${file_og_pem_rdata_fpkm} ]; then
		pem_fit_file=${file_og_pem_rdata_fpkm}
	fi

	Rscript ${dir_myscript}/detect_OU_shift_PhylogeneticEM.r \
	--tree_file=${file_og_dated_tree} \
	--trait_file=${file_og_expression_fpkm} \
	--nslots=${NSLOTS} \
	--fit_file=${pem_fit_file} \
	--similarity_method=${similarity_method} \
	--similarity_threshold=0.99 \
	--ceil_negative=0

	mv PhylogeneticEM.RData ${file_og_pem_rdata_fpkm}
	mv PhylogeneticEM.tree.tsv ${file_og_pem_tree_fpkm}
	mv PhylogeneticEM.regime.tsv ${file_og_pem_regime_fpkm}
	mv PhylogeneticEM.leaf.tsv ${file_og_pem_leaf_fpkm}
	mv PhylogeneticEM.plot.pdf ${file_og_pem_plot_fpkm}

else
	echo "`date`: Skipped: ${task}"
fi

num_gene=`grep -e "^>" ${dir_og_cds_fasta}/${og_id}.cds.fasta | wc -l`
if [ ${num_gene} -ge 2000 ]; then
	max_nshift=10
else 
	max_nshift=0
fi

task="l1ou, TPM"
if [[ ( ! -s ${file_og_l1ou_fit_rdata} || ! -s ${file_og_l1ou_fit_tree} || ! -s ${file_og_l1ou_fit_regime} || ! -s ${file_og_l1ou_fit_leaf} ) && ${run_l1ou_original} -eq 1 ]]; then
	echo "`date`: Start: ${task}"

	fit_ind_file=''
	if [ ${l1ou_use_fit_ind_file} -eq 1 ] && [ -s ${file_og_l1ou_fit_rdata} ]; then
		fit_ind_file=${file_og_l1ou_fit_rdata}
	fi

	taskset -c ${cpu_id} Rscript ${dir_myscript}/detect_OU_shift_l1ou.r \
	--max_nshift=${max_nshift} \
	--tree_file=${file_og_dated_tree} \
	--trait_file=${file_og_expression} \
	--nslots=${NSLOTS} \
	--similarity_method=${similarity_method} \
	--similarity_threshold=0.99 \
	--ceil_negative=0 \
	--criterion=${l1ou_criterion} \
	--nbootstrap=${l1ou_nbootstrap} \
	--fit_ind_file=${fit_ind_file} \
	--fit_conv_file='' \
	--alpha_upper=${l1ou_alpha_upper} \
	--detect_convergence=${l1ou_convergence}

	mv fit_ind.RData ${file_og_l1ou_fit_rdata}
	mv l1ou_tree.tsv ${file_og_l1ou_fit_tree}
	mv l1ou_regime.tsv ${file_og_l1ou_fit_regime}
	mv l1ou_leaf.tsv ${file_og_l1ou_fit_leaf}
	mv l1ou_plot.pdf ${file_og_l1ou_fit_plot}
	if [ ${l1ou_nbootstrap} -gt 0 ]; then
		cp l1ou_bootstrap.tsv ${dir_l1ou_bootstrap}/${l1ou_bootstrap}
	fi
else
	echo "`date`: Skipped: ${task}"
fi

task="l1ou, FPKM"
if [[ ( ! -s ${file_og_l1ou_fit_rdata_fpkm} || ! -s ${file_og_l1ou_fit_tree_fpkm} || ! -s ${file_og_l1ou_fit_regime_fpkm} || ! -s ${file_og_l1ou_fit_leaf_fpkm} ) && ${run_l1ou_fpkm} -eq 1 ]]; then
	echo "`date`: Start: ${task}"

	fit_ind_file=''
	if [ ${l1ou_use_fit_ind_file} -eq 1 ] && [ -s ${file_og_l1ou_fit_rdata_fpkm} ]; then
		fit_ind_file=${file_og_l1ou_fit_rdata_fpkm}
	fi

	taskset -c ${cpu_id} Rscript ${dir_myscript}/detect_OU_shift_l1ou.r \
	--max_nshift=${max_nshift} \
	--tree_file=${file_og_dated_tree} \
	--trait_file=${file_og_expression_fpkm} \
	--nslots=${NSLOTS} \
	--similarity_method=${similarity_method} \
	--similarity_threshold=0.99 \
	--ceil_negative=0 \
	--criterion=${l1ou_criterion} \
	--nbootstrap=${l1ou_nbootstrap} \
	--fit_ind_file=${fit_ind_file} \
	--fit_conv_file='' \
	--alpha_upper=${l1ou_alpha_upper} \
	--detect_convergence=${l1ou_convergence}

	mv fit_ind.RData ${file_og_l1ou_fit_rdata_fpkm}
	mv l1ou_tree.tsv ${file_og_l1ou_fit_tree_fpkm}
	mv l1ou_regime.tsv ${file_og_l1ou_fit_regime_fpkm}
	mv l1ou_leaf.tsv ${file_og_l1ou_fit_leaf_fpkm}
	mv l1ou_plot.pdf ${file_og_l1ou_fit_plot_fpkm}
	if [ ${l1ou_nbootstrap} -gt 0 ]; then
		cp l1ou_bootstrap.tsv ${dir_l1ou_bootstrap}/${l1ou_bootstrap_fpkm}
	fi
else
	echo "`date`: Skipped: ${task}"
fi

summary_flag=0
if [ ${run_summary} -eq 1 ]; then
	if [ ! -s ${file_og_stat_branch} ] || [ ! -s ${file_og_stat_tree} ]; then
		echo "Stat file(s) not found. Starting..."
		summary_flag=1
	else
		outfiles=( `set | grep "^file_og_" | grep -v "^file_og_stat_" | grep -v "^file_og_tree_plot" | sed -e "s/=.*//"` )
		for file in ${outfiles[@]}; do
			if [ -e `eval echo '$'${file}` ]; then
				if [ `eval echo '$'${file}` -ot ${file_og_stat_branch} ] && [ `eval echo '$'${file}` -ot ${file_og_stat_tree} ]; then
					echo File older than stat files: `eval echo '$'${file}`
				else
					echo File newer than stat files: `eval echo '$'${file}`
					summary_flag=1					
				fi
			else
				echo File not found: `eval echo '$'${file}`
			fi
		done
	fi
fi

task="summary statistics"
if [ ${summary_flag} -eq 1 ]; then
	echo "`date`: Start: ${task}"

	unzip -q ${file_og_notung_root}
	unzip -q ${file_og_notung_reconcil}

	python ${dir_myscript}/orthogroup_statistics.py \
	--species_tree ${species_tree} \
	--unaligned_aln ${file_og_cds_fasta} \
	--trimal_aln ${file_og_trimal} \
	--iqtree_tree ${file_og_iqtree_tree} \
	--iqtree_model ${file_og_iqtree_model} \
	--root_tree ${file_og_root} \
	--root_log ${file_og_root_log} \
	--notung_root_log ./${og_id}.notung.root/${og_id}.iqtree.nwk.rooting.ntglog \
	--notung_reconcil_stats ./${og_id}.notung.reconcil/${og_id}.root.nwk.reconciled.stats.txt \
	--dated_tree ${file_og_dated_tree} \
	--dated_log ${file_og_dated_tree_log} \
	--hyphy_tree_dnds ${file_og_hyphy_dnds} \
	--l1ou_prefix tpm fpkm \
	--l1ou_tree ${file_og_l1ou_fit_tree} ${file_og_l1ou_fit_tree_fpkm} \
	--l1ou_regime ${file_og_l1ou_fit_regime} ${file_og_l1ou_fit_regime_fpkm} \
	--l1ou_leaf ${file_og_l1ou_fit_leaf} ${file_og_l1ou_fit_leaf_fpkm} \
	--phylogeneticem_prefix tpm fpkm \
	--phylogeneticem_tree ${file_og_pem_tree} ${file_og_pem_tree_fpkm} \
	--phylogeneticem_regime ${file_og_pem_regime} ${file_og_pem_regime_fpkm} \
	--phylogeneticem_leaf ${file_og_pem_leaf} ${file_og_pem_leaf_fpkm} \
	--expression_prefix tpm fpkm \
	--expression ${file_og_expression} ${file_og_expression_fpkm} \
	--mapdnds_tree_dn ${file_og_mapdnds_dn} \
	--mapdnds_tree_ds ${file_og_mapdnds_ds} \
	--scm_intron ${file_og_scm_intron_summary} \
	--scm_chromosome ${file_og_scm_chromosome_summary}

	cp orthogroup.branch.tsv ${file_og_stat_branch}
	cp orthogroup.tree.tsv ${file_og_stat_tree}
else
	echo "`date`: Skipped: ${task}"
fi

task="tree_plot"
if [ ! -s ${file_og_tree_plot} ] && [ ${run_tree_plot} -eq 1 ]; then
	echo "`date`: Start: ${task}"

	Rscript ${dir_myscript}/stat_branch2tree_plot.r \
	--stat_branch=${file_og_stat_branch} \
	--trait_file=${file_og_expression_fpkm} \
	--max_delta_intron_present=-0.5 \
	--pcm_prefix='l1ou_fpkm_'

	mv stat_branch2tree_plot.pdf ${file_og_tree_plot}
else
	echo "`date`: Skipped: ${task}"
fi

if [ ${delete_tmp_dir} -eq 1 ]; then
	if [ -s ${file_og_stat_branch} ] && [ -s ${file_og_stat_tree} ]; then
		echo "Output files detected. Deleting ${dir_tmp}"
		rm -r ${dir_tmp}
	else
		echo "Output files not found. Leaving ${dir_tmp}"
	fi
else
	echo "Leaving ${dir_tmp}"
fi

###################
echo "`date`: Ending"

function dateComp()
{
    ARG1_SECOND=`date -d "$1" '+%s'`
    ARG2_SECOND=`date -d "$2" '+%s'`
    expr $ARG1_SECOND - $ARG2_SECOND
}







: <<'#_______________CO_______________'


#_______________CO_______________



