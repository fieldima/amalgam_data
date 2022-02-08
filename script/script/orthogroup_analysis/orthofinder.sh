#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -l fat
#$ -l month
#$ -l d_rt=1400:00:00
#$ -l s_rt=1400:00:00
#$ -l s_vmem=4G
#$ -l mem_req=4G
#$ -pe def_slot 256
#$ -t 1-1

#1-66

echo `date`: start

run_of=1

selected_tissues=( brain heart kidney liver ovary testis )
excluded_sci_name=( 'Lepisosteus_oculatus' )

dir_masked_cds='/lustre1/home/kfuku/my_db/Ensembl/release-91/cds.mask'
dir_cds_mask_longest='/lustre1/home/kfuku/my_db/Ensembl/release-91/cds.mask.longest'
dir_pep_mask_longest='/lustre1/home/kfuku/my_db/Ensembl/release-91/pep.mask.longest'
dir_idm="/lustre1/home/kfuku/my_db/Ensembl/release-91/id_mapping"
dir_tissue_mean='/lustre1/home/kfuku/my_db/Ensembl/release-91/curated_transcriptome/tissue_mean'
dir_of='/lustre1/home/kfuku/my_db/Ensembl/release-91/orthofinder'
dir_timetree='/lustre1/home/kfuku/my_db/Ensembl/release-91/timetree'
dir_script="/lustre1/home/kfuku/my_script"
species_tree='/lustre1/home/kfuku/my_db/Ensembl/release-91/timetree/species_timetree.nwk'

if [ ! -e ${dir_pep_mask_longest} ]; then
	mkdir ${dir_pep_mask_longest}
fi
if [ ! -e ${dir_cds_mask_longest} ]; then
	mkdir ${dir_cds_mask_longest}
fi

if [ ! -e ${dir_of} ]; then
	mkdir ${dir_of}
fi
if [ -e ${dir_timetree}/species.txt ]; then
	rm ${dir_timetree}/species.txt
fi
touch ${dir_timetree}/species.txt

echo 'selected sci_name:'
files=( `find ${dir_tissue_mean}/*.tsv` )
for file in ${files[@]}; do
	header=`head -n 1 ${file}`
	flag=1
	for tissue in ${selected_tissues[@]}; do
		if [ ! -n "`echo ${header} | grep -e ${tissue}`" ]; then
			flag=0
		fi
	done
	sci_name=`basename ${file} | sed -e "s/\..*//"`
	sci_name_idm=`echo ${sci_name:0:1} | awk '{print tolower($0)}'``echo ${sci_name} | sed -e "s/\..*//" -e "s/.*_//"`
	if [ ${flag} -eq 1 ]; then
		for es in ${excluded_sci_name[@]}; do
			if [ ${sci_name} = ${es} ]; then
				flag=0
			fi
		done
	fi
	if [ ${flag} -eq 1 ] ; then
		echo ${sci_name}
		echo ${sci_name} | sed -e "s/_/ /" >> ${dir_timetree}/species.txt
		if [ ! -s ${dir_pep_mask_longest}/${sci_name}.fasta ]; then
			
			python ${dir_script}/ensembl_get_longest_transcript.py \
			--file_cds `find ${dir_masked_cds}/${sci_name}*` \
			--file_id  `find ${dir_idm}/${sci_name_idm}*` \
			--file_out ${dir_pep_mask_longest}/tmp.${sci_name}.pep.fasta \
			--translate_table 1 \
			--remove_stop 1

			python ${dir_script}/ensembl_get_longest_transcript.py \
			--file_cds `find ${dir_masked_cds}/${sci_name}*` \
			--file_id  `find ${dir_idm}/${sci_name_idm}*` \
			--file_out ${dir_cds_mask_longest}/tmp.${sci_name}.cds.fasta \
			--translate_table 0 \
			--remove_stop 1

			if [ -s ${dir_pep_mask_longest}/tmp.${sci_name}.pep.fasta ]; then
				sed -e "s/^>/>${sci_name}_/" ${dir_pep_mask_longest}/tmp.${sci_name}.pep.fasta \
				> ${dir_pep_mask_longest}/${sci_name}.fasta
				rm ${dir_pep_mask_longest}/tmp.${sci_name}.pep.fasta
			fi

			if [ -s ${dir_cds_mask_longest}/tmp.${sci_name}.cds.fasta ]; then
				sed -e "s/^>/>${sci_name}_/" ${dir_cds_mask_longest}/tmp.${sci_name}.cds.fasta \
				> ${dir_cds_mask_longest}/${sci_name}.cds.mask.longest.fasta
				rm ${dir_cds_mask_longest}/tmp.${sci_name}.cds.fasta
			fi
		fi
	fi
done

if [ ${run_of} -eq 1 ]; then
	echo 'Start: orthofinder'
	conda activate py27

	orthofinder \
	-t ${NSLOTS} \
	-a $[${NSLOTS}/16] \
	-M dendroblast \
	-S blast \
	-s ${species_tree} \
	-I 1.5 \
	-p ${dir_of} \
	-f ${dir_pep_mask_longest}


	#-fg /lustre1/home/kfuku/my_db/Ensembl/release-91/cds.mask.translate/Results_Mar16_2

	# fg: from groups
fi


echo `date`: end





: <<'#_______________CO_______________'


#_______________CO_______________












