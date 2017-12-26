# Instructions for using the Nimerlab RNAseq pipeline (version 0.1):

## This pipeline is compatable with UM Pegasus using project nimerlab

1. Download https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh to your nethome folder.
2. Install miniconda: 'bash ~/Miniconda3-latest-Linux-x86_64.sh'
3. After installation type the following commands:
	- module rm python
	- conda env create -f /projects/ctsi/nimerlab/DANIEL/tools/nimerlab-pipelines/RNAseq/installation/RNAseq.yml
	- cp /projects/ctsi/nimerlab/DANIEL/tools/nimerlab-pipelines/RNAseq/installation/fastq_screen.conf ~/miniconda3/envs/RNAseq/share/fastq-screen-0.11.3-0/
4. Copy '/projects/ctsi/nimerlab/DANIEL/tools/nimerlab-pipelines/RNAseq/RNAseq_experiment_file.yml' into your run folder.
5. Copy '/projects/ctsi/nimerlab/DANIEL/tools/nimerlab-pipelines/RNAseq/RNAseq.py' into your run folder.
6. From your run folder, run analysis with this command (replacing 'RNAseq_expiermental_file.yml' with your experimental filename:
	- bsub -q general -n 1 -R 'rusage[mem=1000]' -W 120:00 -o RNAseq.out -e RNAseq.err <<< 'module rm python;source activate RNAseq;python RNAseq.py -f RNAseq_experimental_file.yml' 
