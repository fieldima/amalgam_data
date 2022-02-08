#! /bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -pe def_slot 1
#$ -l s_vmem=8G
#$ -l mem_req=8G
##$ -l epyc
##$ -l d_rt=1400:00:00
##$ -l s_rt=1400:00:00
#$ -t 12

##1-21

echo `date`: start

dist_method=pearson
mapping_rate_cutoff=0.2

mode='nig'
sra_date=2018_5_1
if [ ${mode} = 'nig' ]; then
	conda activate r-sva
	sra_table="/lustre6/home/lustre1/kfuku/my_db/Ensembl/release-91/sra/${sra_date}/sra_table_mapped_${sra_date}.tsv"
	dir_work="/lustre6/home/lustre1/kfuku/my_project/convergence_duplication/20190429_tmm_normalization"
	dir_out="/lustre6/home/lustre1/kfuku/my_db/Ensembl/release-91/curated_transcriptome/${sra_date}/tmm_rpkm"
	#dir_exp="/lustre6/home/lustre1/kfuku/my_db/Ensembl/release-91/kallisto_summary/${sra_date}/tpm.tmm.kallisto.tsv"
	dir_exp_curated="/lustre6/home/lustre1/kfuku/my_db/Ensembl/release-91/kallisto_summary/${sra_date}/fpkm.tmm.kallisto.gene.log.tsv"
	dir_idm="/lustre6/home/lustre1/kfuku/my_db/Ensembl/release-91/id_mapping"
	dir_script="/lustre6/home/lustre1/kfuku/my_script"
elif [ ${mode} = 'imac' ]; then
	dir_ensembl="/Users/kef74yk/Dropbox_w/db/Ensembl"
	sra_table="${dir_ensembl}/release-91/sra/${sra_date}/sra_table_mapped_${sra_date}.tsv"
	dir_work="/Users/kef74yk/Dropbox_p/data/04_Convergence_Duplication/20190429_tmm_normalization"
	dir_out="${dir_ensembl}/release-91/curated_transcriptome/${sra_date}_tmm"
	#dir_exp="${dir_ensembl}/release-91/kallisto_summary/${sra_date}/tpm.tmm.kallisto.tsv"
	dir_exp_curated="${dir_ensembl}/release-91/kallisto_summary/${sra_date}/fpkm.tmm.kallisto.gene.log.tsv"
	dir_idm="${dir_ensembl}/release-91/id_mapping"
	dir_script="/Users/kef74yk/Dropbox_w/script"
fi

dir_out_tau=${dir_out}/tau
dir_out_tc=${dir_out}/tc
dir_out_r2=${dir_out}/r2
dir_out_sra=${dir_out}/sra
dir_out_uncorrected_tissue_mean=${dir_out}/uncorrected_tissue_mean
dir_out_tissue_mean=${dir_out}/tissue_mean
dir_out_plot=${dir_out}/plot
dir_out_sva=${dir_out}/sva

directories=( `set | grep "^dir_" | sed -e "s/=.*//"` )
for d in ${directories[@]}; do
	if [ ! -e `eval echo '$'${d}` ] && \
	[ "${d}" != "dir_ensembl" ] && \
	[ "${d}" != "dir_work" ] && \
	[ "${d}" != "dir_hyphy_out" ] && \
	[ "${d}" != "dir_transcriptome" ]; then
		echo creating: `eval echo '$'${d}`
		mkdir -p `eval echo '$'${d}`
	fi
done

infiles=( `ls ${dir_exp_curated}` )
echo "all infiles: ${infiles[@]}"
if [ ${mode} = 'nig' ]; then
	ind=$[${SGE_TASK_ID}-1]
	infiles=( ${infiles[${ind}]} )
fi
echo "infiles to be processed: ${infiles[@]}"

for infile in ${infiles[@]}; do
	sci_name=`echo ${infile} | sed -e "s/\..*//"`
	dir_tmp=${dir_work}/${sci_name}
	file_exp=${infile}
	sci_name_idm=`echo ${sci_name:0:1} | awk '{print tolower($0)}'``echo ${sci_name} | sed -e "s/\..*//" -e "s/.*_//"`
	file_idm=`ls ${dir_idm} | grep -e "^${sci_name_idm}"`

	# prepare and set directories
	directories=( `set | grep "^dir_" | sed -e "s/=.*//"` )
	for d in ${directories[@]}; do
		if [ ! -e `eval echo '$'${d}` ]; then
			echo creating: `eval echo '$'${d}`
			mkdir -p `eval echo '$'${d}`
		fi
	done
	cd ${dir_tmp}

#	if [ ! -s ${dir_exp_curated}/${sci_name}.gene.log.tsv ]; then
#		echo `date`: Start: transcriptome_data_preparation.r
#
#		Rscript ${dir_script}/transcriptome_data_preparation.r \
#		${dir_exp}/${file_exp} \
#		${sci_name}.gene.log.tsv \
#		${dir_idm}/${file_idm} \
#		1 \
#		1
#
#		cp ${sci_name}.gene.log.tsv ${dir_exp_curated}/${sci_name}.gene.log.tsv
#	else
#		echo `date`: Skipped: transcriptome_data_preparation.r
#	fi

	if [ ! -s ${dir_out_tc}/${sci_name}.tc.tsv ]; then
		echo `date`: Start: transcriptome_curation2.r
		
		Rscript ${dir_script}/transcriptome_curation2.r \
		${dir_exp_curated}/${sci_name}.gene.log.tsv \
		${sra_table} \
		${dir_tmp} \
		${dist_method} \
		${mapping_rate_cutoff} \
		0 \
		1 \
		'brain|heart|kidney|liver|ovary|testis'

		if [ $? -eq 0 ]; then
			cp *.pdf ${dir_out_plot}
			cp *.sva.*.RData ${dir_out_sva}
			cp *.r2.tsv ${dir_out_r2}
			cp ${sci_name}.tau.tsv ${dir_out_tau}
			cp ${sci_name}.uncorrected.tissue.mean.tsv ${dir_out_uncorrected_tissue_mean}
			cp ${sci_name}.tissue.mean.tsv ${dir_out_tissue_mean}
			cp ${sci_name}.sra.tsv ${dir_out_sra}
			cp ${sci_name}.tc.tsv ${dir_out_tc}
			rm -r ${dir_tmp}
		fi
	else
		echo `date`: Skipped: transcriptome_curation2.r
	fi

done

echo `date`: end





: <<'#_______________CO_______________'



#_______________CO_______________



