#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=48
#SBATCH --qos=bsc_ls
#SBATCH --output=homfam_SoP_greasy_complete.txt
module load python
echo "filename,original_num_seqs,original_num_cols,clean_num_seqs,clean_num_cols,min_SoP_alignment,"\
    "max_SoP_alignment,median_alignment,std_alignment,alignment_SoP,trimmed_alignment_size_gappyout,"\
    "min_SoP_trimmed_alignment_gappyout,max_SoP_trimmed_alignment_gappyout,median_alignment_gappyout,"\
    "std_alignment_SoP_gappyout,trimmed_alignment_SoP_gappyout,trimmed_alignment_size_strictplus,"\
    "min_SoP_trimmed_alignment_strictplus,max_SoP_trimmed_alignment_strictplus,median_alignment_strictplus,"\
    "std_alignment_SoP_strictplus,trimmed_alignment_SoP_strictplus" > homfam_SoP_stats_complete.txt

FILE=homfam_SoP_calls.txt
/apps/GREASY/2.2/INTEL/IMPI/bin/greasy $FILE >> homfam_SoP_stats_complete.txt
