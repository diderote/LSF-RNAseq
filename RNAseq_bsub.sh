#!/bin/bash

#BSUB -J <<LSF_RNAseq_JOB>>
#BSUB -R "rusage[mem=3000]"
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -W 120:00
#BSUB -n 1
#BSUB -q general
#BSUB -P <<project>>

### Add a project and job name above before submission ###
### Un-comment one of the following python script commands before submission ###
### Change the path and file names as needed before submission ###

module rm python share-rpms65
source activate RNAseq

# python /path/to/RNAseq.py -f /path/to/RNAseq_experiment_file.yml -t /path/to/RNAseq.ipynb

# python /path/to/RNAseq.py -f /path/to/RNAseq_experiment_file.yml --no-notebook

# python /path/to/RNAseq.py \
#    -f /path/to/RNAseq_experiment_file.yml \
#    -t /path/to/RNAseq.ipynb \
#    -o /path/to/output/newname_RNAseq.ipynb