#! /bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -pe def_slot 8
#$ -l s_vmem=4G
#$ -l mem_req=4G
##$ -l epyc
#$ -l d_rt=62:00:00:00
#$ -l s_rt=62:00:00:00
#$ -t 1

#SBATCH -J gfe
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -c 1
##SBATCH -p standard
#SBATCH -p test
##SBATCH -t 2:0:0:0
##SBATCH -t 1:0:0
##SBATCH --mail-type=ALL
##SBATCH --mail-user=<aaa@bbb.com>
#SBATCH --mem-per-cpu=4G
#SBATCH --output=gfe_%A_%a.out
#SBATCH --error=gfe_%A_%a.err
##SBATCH --array=1

echo "`date`: Starting"
ulimit -s unlimited

# Change these directories for your custom-made analysis
dir_gfe="${PWD}/../gfe_data" # gfe input and output directory
dir_script="${PWD}" # directory where gfe_util.sh and gfe_*_cmd.sh locate
gfe_image="${PWD}/../gfe.sif" # path to the singularity image

source ${dir_script}/script/gfe_util.sh # loading utility functions
set_singularity_command
variable_SGEnizer
set_singularityenv

cd ${dir_gfe}
${singularity_command} ${gfe_image} < ${dir_script}/gfe_cdsQualityCheck_cmd.sh

echo "`date`: Ending"
