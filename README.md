# Nimerlab RNAseq pipeline

This pipline handles processing and analyses for RNAseq data.  When fully implemented from start to finish it performs the following tasks:
1. Fastq merging
2. Screening for contamination of other genomes
3. Fastq quality measurements
4. Adapter and quality fastq trimming
5. Spike-in assessment (if applicable)
6. Transcriptome alignment using STAR (with RSEM expected counts) and Kallisto
7. Generation of scaled (reads per million) bigwig files from transcriptome alignment
7. Differential expression using DESeq2
8. Overlap with signifiantly differentially expressed genes (q<0.05, 2 & 1.5 FC) with Sleuth (Kallisto) q<0.05 DE genes 
9. Heatmap of significantly differentially expressed genes using variance stabilized expected counts
10. GO and KEGG enrichment of DE genes
11. GSEA from test statistic ranked gene lists from DESeq2
12. PCA analysis
13. Overlap of multiple gene lists (including from differential expression) and scaled venn diagram output.

The pipeline handles multiple entry/exit points and can parse complex experimental designs and compensation types for DE.
In case of error, the pipeline restarts from the last completed step. Progress is tracked in a .log file in the output directory.
All submission scripts, error and output files are saved.

## Instructions for installation and use:

### This pipeline is compatable with UM Pegasus.

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
	- bsub -q general -n 1 -R 'rusage[mem=1000]' -W 120:00 -o RNAseq.out -e RNAseq.err <<< 'module rm python share-rpms65;source activate RNAseq;./RNAseq.py -f RNAseq_experimental_file.yml' 
8. In case of error, use the above command to pick up from last completed step.

### RNAseq.py can be imported to python as a module with the following attributes:
	New Class:
		- exp = Experiment(scratch='/path/to/scratch/folder/', date='', name='experiment_name', out_dir='/path/to/results/folder', job_folder='/path/to/job/folder',count_matrix=pd.DataFrame(index=genes,columns=samplenames),spike_counts=pd.DataFrame(),stop='',genome='hg38 or mm10',sample_number=int(), samples={1:'sample_name'}, job_id=[],de_groups={},norm='bioinformatic',designs={},overlaps={},tasks_complete=[],de_results={},sig_lists={},overlap_results={},de_sig_overlap={})
	
	Experiment object can be passed to the following functions and returned: exp = function(exp):
	- RNAseq.DESeq2()
	- RNAseq.fastq_cat()       
    - RNAseq.fastq_screen()
    - RNAseq.GO_enrich()
    - RNAseq.fastqc()
    - RNAseq.GSEA()
    - RNAseq.final_qc()
	- RNAseq.PCA()
	- RNAseq.rsem()
	- RNAseq.Sleuth()
	- RNAseq.kallisto()
	- RNAseq.sigs()
	- RNAseq.clustermap()
	- RNAseq.spike()
	- RNAseq.splicing()
	- RNAseq.overlaps()
	- RNAseq.stage()
	- RNAseq.count_matrix()      
	- RNAseq.trim()

	helper functions:
	- RNAseq.parse_yaml() takes required -f yaml file and parses an new experimental object or loads an incomplete one     
	- RNAseq.plot_venn2(Series=pd.Series(), string_name_of_overlap='', folder='') scaled venn of an overlap
	- RNAseq.scaled_bigwig(in_bam='/path/to/bam',out_bw='/path/to/bw',job_log_folder=exp.job_folder,name='',genome=exp.genome) scales bam to bigwig rpm
	- RNAseq.job_wait(id_list=exp.job_id, job_log_folder=exp.job_folder, log_file=exp.log_file) waits for submitted job to finish
	- RNAseq.send_job(command_list=[], job_name='', job_log_folder=exp.job_folder, q='', mem='', log_file=exp.log_file) sends job to LSF resource manager
	- RNAseq.enrichr(gene_list=[], description='', out_dir='')
	- RNAseq.plot_PCA()

## To Do:
1. STAR 2 pass option followed by straberry and DEXseq for splicing.
2. Optional tSNE
3. ICA analysis with chi-square for comparisons for component.
4. Add RUVseq for ERCC for sample visualization. (EDAseq inlane GC normalization?)
5. Change trimming to optional (STAR and Kallisto perform well without trimming) or only trim if flagged as problem.

