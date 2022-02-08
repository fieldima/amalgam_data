#! /bin/sh

### Start: Modify this block to tailor your analysis ###

run_busco=1

delete_busco_tmp_dir=1
busco_lineage='vertebrata_odb10' # See below for available lineages

### End: Modify this block to tailor your analysis ###

### Modify below if you need to add a new analysis or need to fix some bugs ###

dir_gfe="/gfe_data"
dir_myscript="/script/script"
source /home/.bashrc # Load conda environment
source ${dir_myscript}/gfe_util.sh # Load utility functions
conda activate busco

echo "`date`: Starting Singularity environment"
echo "pwd: "`pwd`
echo "Python executable: `which python`"
echo "dir_gfe: ${dir_gfe}"
echo "NSLOTS: ${NSLOTS}"
echo "JOB_ID: ${JOB_ID}"
echo "MEM_PER_SLOT: ${MEM_PER_SLOT}"
echo "SGE_TASK_ID: ${SGE_TASK_ID}"

if [ -z ${MEM_PER_SLOT} ]; then
	MEM_PER_SLOT=3
fi
MEM_PER_HOST=$[${MEM_PER_SLOT}*${NSLOTS}]
echo "MEM_PER_HOST: ${MEM_PER_HOST}"

dir_og="${dir_gfe}/orthogroup"
dir_sp_cds="${dir_gfe}/species_cds"
dir_busco="${dir_gfe}/species_cds_busco"

infiles=( `ls ${dir_sp_cds}/*.{fa,fasta}` )
infile=${infiles[$[${SGE_TASK_ID}-1]]}
sci_name=`basename ${infile} | sed -e "s/_/|/" -e "s/_.*//" -e "s/|/_/"`
echo "${#infiles[@]} input fasta files were detected in: ${dir_sp_cds}"
echo "Processing ${SGE_TASK_ID}th file: ${infile}"
echo "Scientific name: ${sci_name}"

dir_busco_sp_out="${dir_busco}/${sci_name}"
file_busco_summary="${dir_busco_sp_out}/run_${busco_lineage}/short_summary.txt"
file_full_table="${dir_busco_sp_out}/run_${busco_lineage}/full_table.tsv"
file_full_table_rename="${dir_busco}/${sci_name}.busco.tsv"

directory_generator

task='BUSCO'
if [[ ! -s "${file_busco_summary}" ]] && [[ ! -e "${file_full_table_rename}" ]] && [[ ${run_busco} -eq 1 ]]; then
	echo "`date`: Start : ${task}"

	cd ${dir_busco}

	busco \
	--in ${infile} \
	--out ${sci_name} \
	--cpu ${NSLOTS} \
	--evalue 1e-03 \
	--mode transcriptome \
	--lineage_dataset ${busco_lineage} \
	--limit 3 \
    --force

	echo "`date`: End : ${task}"
else
	echo "`date`: Skipped : ${task}"
fi

if [ -e "${file_busco_summary}" ] && [ -e "${file_full_table}" ]; then
    echo "BUSCO outputs found."
	cat "${file_busco_summary}"
    cp "${file_full_table}" "${file_full_table_rename}"

    if [[ ${delete_busco_tmp_dir} -eq 1 ]]; then
        echo "Deleting tmp directory: ${dir_busco_sp_out}"
        rm -rf ${dir_busco_sp_out}
    fi
else
	echo "BUSCO outputs not found. Exiting."
	exit 1
fi



echo "`date`: Exiting Singularity environment"


: <<'#_______________CO_______________'

Datasets available to be used with BUSCOv4 as of 2019/11/27:

 bacteria_odb10
     - acidobacteria_odb10
     - actinobacteria_phylum_odb10
         - actinobacteria_class_odb10
             - corynebacteriales_odb10
             - micrococcales_odb10
             - propionibacteriales_odb10
             - streptomycetales_odb10
             - streptosporangiales_odb10
         - coriobacteriia_odb10
             - coriobacteriales_odb10
     - aquificae_odb10
     - bacteroidetes-chlorobi_group_odb10
         - bacteroidetes_odb10
             - bacteroidia_odb10
                 - bacteroidales_odb10
             - cytophagia_odb10
                 - cytophagales_odb10
             - flavobacteriia_odb10
                 - flavobacteriales_odb10
             - sphingobacteriia_odb10
         - chlorobi_odb10
     - chlamydiae_odb10
     - chloroflexi_odb10
     - cyanobacteria_odb10
         - chroococcales_odb10
         - nostocales_odb10
         - oscillatoriales_odb10
         - synechococcales_odb10
     - firmicutes_odb10
         - bacilli_odb10
             - bacillales_odb10
             - lactobacillales_odb10
         - clostridia_odb10
             - clostridiales_odb10
             - thermoanaerobacterales_odb10
         - selenomonadales_odb10
         - tissierellia_odb10
             - tissierellales_odb10
     - fusobacteria_odb10
         - fusobacteriales_odb10
     - planctomycetes_odb10
     - proteobacteria_odb10
         - alphaproteobacteria_odb10
             - rhizobiales_odb10
                 - rhizobium-agrobacterium_group_odb10
             - rhodobacterales_odb10
             - rhodospirillales_odb10
             - rickettsiales_odb10
             - sphingomonadales_odb10
         - betaproteobacteria_odb10
             - burkholderiales_odb10
             - neisseriales_odb10
             - nitrosomonadales_odb10
         - delta-epsilon-subdivisions_odb10
             - deltaproteobacteria_odb10
                 - desulfobacterales_odb10
                 - desulfovibrionales_odb10
                 - desulfuromonadales_odb10
             - epsilonproteobacteria_odb10
                 - campylobacterales_odb10
         - gammaproteobacteria_odb10
             - alteromonadales_odb10
             - cellvibrionales_odb10
             - chromatiales_odb10
             - enterobacterales_odb10
             - legionellales_odb10
             - oceanospirillales_odb10
             - pasteurellales_odb10
             - pseudomonadales_odb10
             - thiotrichales_odb10
             - vibrionales_odb10
             - xanthomonadales_odb10
     - spirochaetes_odb10
         - spirochaetia_odb10
             - spirochaetales_odb10
     - synergistetes_odb10
     - tenericutes_odb10
         - mollicutes_odb10
             - entomoplasmatales_odb10
             - mycoplasmatales_odb10
     - thermotogae_odb10
     - verrucomicrobia_odb10
 archaea_odb10
     - thaumarchaeota_odb10
     - thermoprotei_odb10
         - thermoproteales_odb10
         - sulfolobales_odb10
         - desulfurococcales_odb10
     - euryarchaeota_odb10
         - thermoplasmata_odb10
         - methanococcales_odb10
         - methanobacteria_odb10
         - methanomicrobia_odb10
             - methanomicrobiales_odb10
         - halobacteria_odb10
             - halobacteriales_odb10
             - natrialbales_odb10
             - haloferacales_odb10
 eukaryota_odb10
     - alveolata_odb10
         - apicomplexa_odb10
             - aconoidasida_odb10
                 - plasmodium_odb10
             - coccidia_odb10
     - euglenozoa_odb10
     - fungi_odb10
         - ascomycota_odb10
             - dothideomycetes_odb10
                 - capnodiales_odb10
                 - pleosporales_odb10
             - eurotiomycetes_odb10
                 - chaetothyriales_odb10
                 - eurotiales_odb10
                 - onygenales_odb10
             - leotiomycetes_odb10
                 - helotiales_odb10
             - saccharomycetes_odb10
             - sordariomycetes_odb10
                 - glomerellales_odb10
                 - hypocreales_odb10
         - basidiomycota_odb10
             - agaricomycetes_odb10
                 - agaricales_odb10
                 - boletales_odb10
                 - polyporales_odb10
             - tremellomycetes_odb10
         - microsporidia_odb10
         - mucoromycota_odb10
             - mucorales_odb10
     - metazoa_odb10
         - arthropoda_odb10
             - arachnida_odb10
             - insecta_odb10
                 - endopterygota_odb10
                     - diptera_odb10
                     - hymenoptera_odb10
                     - lepidoptera_odb10
                 - hemiptera_odb10
         - mollusca_odb10
         - nematoda_odb10
         - vertebrata_odb10
             - actinopterygii_odb10
                 - cyprinodontiformes_odb10
             - tetrapoda_odb10
                 - mammalia_odb10
                     - eutheria_odb10
                         - euarchontoglires_odb10
                             - glires_odb10
                             - primates_odb10
                         - laurasiatheria_odb10
                             - carnivora_odb10
                             - cetartiodactyla_odb10
                 - sauropsida_odb10
                     - aves_odb10
                         - passeriformes_odb10
     - stramenopiles_odb10
     - viridiplantae_odb10
         - chlorophyta_odb10
         - embryophyta_odb10
             - liliopsida_odb10
                 - poales_odb10
             - eudicots_odb10
                 - brassicales_odb10
                 - fabales_odb10
                 - solanales_odb10

#_______________CO_______________
