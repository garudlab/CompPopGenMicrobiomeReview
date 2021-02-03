#!/bin/bash
#$ -N SSU_raxml
#$ -e /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/SSU_raxml_error
#$ -o /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/SSU_raxml_output
#$ -l h_data=8G
#$ -l time=24:00:00
#$ -l highp
#$ -m bea



module load raxml/8.2.4


raxmlHPC -T 4 -f a -m GTRGAMMA -p 12345 -x 12345 -o NC_005042.1.353331-354795 -# autoMRE \
    -s /u/home/w/wrshoema/project-ngarud/negative_selection_microbiome/data/tree/SSU_rRNA_all_species_muscle_clean.frn \
    -n HMP_2 -w /u/home/w/wrshoema/project-ngarud/negative_selection_microbiome/data/tree



# -T = number of threads
# -f = specifies bootstrapping algorithm with ML generating tree at same time
# -m = substitution model, generalized time reversible gamma
# -p = starts tree randomly
# -x = starts tree randomly
# -o = outgroup (name after fasta entry)
# -#  = determines number of bootstrap replicates
# -s = aligned fasta file to be analyzed
# -n = name of output file



raxmlHPC-PTHREADS-SSE3 -T 4 -f a -m GTRGAMMA -p 12345 -x 12345 -o NC_005042.1.353331-354795 -# autoMRE \
    -s /Users/wrshoemaker/GitHub/negative_selection_microbiome/data/tree/SSU_rRNA_all_species_muscle_clean.frn \
    -n HMP_2 -w /Users/wrshoemaker/GitHub/negative_selection_microbiome/data/tree
