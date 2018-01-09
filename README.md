# Nimerlab RNAseq pipeline

This pipline handles processing and analyses for RNAseq data.  When fully implemented from start to finish it performs the following tasks:
1. Fastq merging
2. Screening for contamination of other genomes
3. Fastq quality measurements
4. Adapter and quality fastq trimming
5. Transcriptome alignment using STAR (with RSEM expected counts) and Kallisto
6. Spike-in quantificaiton (if applicable)
6. Differential expression using DESeq2 (from RSEM and spike-in), and overlap with Sleuth (Kallisto) under most conditions
7. Heatmap of significantly differentially expressed genes
8. GO enrichment of DE genes
9. GSEA from test statistic ranked gene lists from DESeq2
10. PCA analysis
11. Overlap of multiple DE significant gene lists.

The pipeline handles multiple entry/exit points and parse complex experimental designs and compensation types for DE.
In case of error, the pipeline restarts from the last completed step. Progress is tracked in a .log file in the output directory.
All submission scripts, error and output files are saved.

## Instructions for installation and use:

### This pipeline is compatable with UM Pegasus using project nimerlab

1. Download https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh to your nethome folder.
2. Install miniconda: 'bash ~/Miniconda3-latest-Linux-x86_64.sh'
3. After installation type the following commands:
	- module rm python
	- conda env create -f /projects/ctsi/nimerlab/DANIEL/tools/nimerlab-pipelines/RNAseq/environment.yml
	- cp /projects/ctsi/nimerlab/DANIEL/tools/nimerlab-pipelines/RNAseq/fastq_screen.conf ~/miniconda3/envs/RNAseq/share/fastq-screen-0.11.3-0/
4. Copy '/projects/ctsi/nimerlab/DANIEL/tools/nimerlab-pipelines/RNAseq/RNAseq_experiment_file.yml' into your run folder.
5. Rename your experiment file as needed and edit the file to fit your experiment using any text editor.
6. Copy '/projects/ctsi/nimerlab/DANIEL/tools/nimerlab-pipelines/RNAseq/RNAseq' into your run folder.
7. From your run folder, run analysis with this command (replacing 'RNAseq_expiermental_file.yml' with your experimental filename:
	- bsub -q general -n 1 -R 'rusage[mem=1000]' -W 120:00 -o RNAseq.out -e RNAseq.err <<< 'module rm python share-rpms65;source activate RNAseq;./RNAseq -f RNAseq_experimental_file.yml' 
8. In case of error, use the above command to pick up from last completed step.

## To Do:
1. Add splicing analyses.
2. Bam inputs.
4. ICA analysis.
