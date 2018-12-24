# LSF RNAseq pipeline

This pipline handles processing and analyses for RNAseq data on the University of Miami's Pegasus Computer Cluster using and LSF resource manager.  (Can be readily adpated to other resource managers: SGE, SLURM)  When fully implemented from start to finish it performs the following tasks:

1. Screening for contamination of other genomes
2. Fastq quality measurements
3. Adapter and quality fastq trimming
4. Spike-in alignment and analysis (if applicable)
5. Alignment
	- Transcriptome alignment using STAR with RSEM expected counts and optional Kallisto (not for single-end) with GENCODE annotations
	- Genome alignment using STAR counts with GENCODE annotations
6. Generation of scaled (counts per million per basepair) bigwig files from alignment (stranded or non-stranded)
7. Optional within-sample GC normalization (useful for experiments across multiple sequencing runs).
8. Between sample normalization options: median-ratios (DESeq2 default), removal of unwated variation (RUVSeq) using ERCC spike-ins or empirical negative controls, upper quantile, full or median normalization (EDAseq).
9. Differential expression using DESeq2 (both default and apeglm or ashr shrunken LFC are reported)
10. PCA plots and expression boxplots of all and experimental samples, raw and normalized counts.
11. Output of normalized counts, blinded, with batch compensation, gc normlalized and raw.
12. Optional overlap with signifiantly differentially expressed genes by DESeq2 (q<0.05, 2FC, 1.5FC, no FC filter) with Sleuth (Kallisto) q<0.05 DE genes 
13. Volcano plots generated from DESeq2 results with signifcant genes
14. Heatmap of significantly differentially expressed genes using regularized log2 expresion counts.
15. GO and KEGG enrichment of DE genes
16. GSEA using ranked gene lists from DESeq2
17. Overlap of differentially expressed genes from different tests with scaled venn-diagram and GO/KEGG enrichment of overlaps.  p-value of the overlap is determined using a hypergeometric test with all expressed genes in the the comparison samples as background.

## Option Details:
* Restart: (yes/no) Whether to restart from the beginninng or check for an incomplete pipeline and pickup from last completed step.
* ERCC_spike: align reads to spike index using STAR.
* Stranded: (yes/no) Will generate two stranded signal bigwig files if stranded.  Otherwise, will generate one per sample.
* Sequencer: Nextseq or Hiseq
* Normalization: 
	* Median-Ratios = default DESeq2 method normalization: median of ratios of observed counts. Also reports apeglm lfc shrunken values.
	* ERCC = Account for unwanted variation between samples using spike in counts. (RUVSeq implementation).  Also reports ashr lfc shrunken values.
	* ERCC_Mixed = Same as ERCC except used when ERCC Mix1 and Mix2 are used in the same experiment.  Will reduce the unwanted variance estimation to subgroup B (common concentrations in both mixes).
	* Empirical = Account for unwanted variation between samples using non-differentially expressed genes to estimate unwanted variance (RUVSeq implementation). Also reports ashr lfc shrunken values.
	* Standard upper quartile, full distribution or median normalizaiton based on library size (EDAseq).
* Alignment_Mode: Whether to align to transcriptome with RSEM-STAR (default) or to genome with STAR.
* Signature_Mode: Output signature will either be DESeq2 differentially expressed genes or an overlap of Slueth and DESeq2.
* GC_Normalizaiton: Yes implements within-lane loess normalzation based on gene GC content (EDASeq).  Recommended for samples sequenced in different sequencing runs.  If not 'Nimer' this adds over an hour of computing time for CG content file generatation.
* Sequencing_type: paired or single end sequencing.

Other options:
*  Spike_matrix: if not aligning and using ERCC, path to file containing ERCCs counts.
*  Count_matrix: if not aligning, path to a file containing gene names as rows(first column), sample names as columns(top row), and gene counts in cells.

Lab options if running outside of Nimer Lab access:
* Scratch folder: Folder for staging.
* RSEM_STAR_index: Specify location of previously generated index for alignment.
	* To create the index for transcriptome alignment use rsem-prepare-reference function.  Recommendations: --star, gencocde .gtf, no_alt genome assembly.
* STAR_index: Specify location of previously generated index for alignment.
* Kallisto_index: Specify location of previously generated index for alignment. Create using default Kallisto index function.
* ERCC_STAR_index: Specify location of previously generated STAR index for spike-ins.
	* To create:  donwload ERCC fasta and gtf, use STAR genomeGenerate function
* Project: project for job submission
* GSEA_jar: Specify location of GSEA.jar (tested with GSEA-3.0.jar)
* GSEA_mouse_gmts: Path to folder containing gmts using mouse ensembl gene IDs instead of human gene names. These can be found in the optional file folder.
* Gene_names: optional. provide path to a pickled dictionary to add an attribute such as 'gene_name' to the count_matrix.  Useful for STAR alinging where gene_id is default, or converging mouse GSEA gene ids to gene name.  

The pipeline handles multiple entry/exit points and can parse complex experimental designs and compensation types for DE.  In case of error, the pipeline restarts from the last completed step. Progress is tracked in a .log file in the output directory.

All submission scripts, error and output files are saved.

# Instructions for installation and use:

## This pipeline is compatable with UM Pegasus.

1. Download https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh to your nethome folder.
2. Install miniconda: 'bash ~/Miniconda3-latest-Linux-x86_64.sh'
3. After installation type the following commands:

- If in pegasus:
> module rm python share-rpms65
- To create the environment:
> conda env create -f /path/to/environment.yml

 - If there is a package error, try removing the version of that package the threw the error from the .yml file and try again.

4. ADDITIONAL SETUP:
	- For use of the contamination screen (fastq screen):
		- generate the relevant bowtie2 indices and add the path were indicated in options_files/fastq_screen.conf. (many premade by illumina at https://support.illumina.com/sequencing/sequencing_software/igenome.html)
		- unhash relvant lines before each desired DATABASE
		- cp fastq_screen.conf ~/miniconda3/envs/RNAseq/share/fastq-screen-0.11.3-0/
	- Download GSEA.jar (> 3.0) from http://www.gsea-msigdb.org/gsea/index.jsp.
		- With mm10, the genes have to be converted to mouse IDs.  Provided in the option_files is a folder of preconverted IDs, with a file 'gencode_gene_dict.pkl' to convert the IDs back to mouse gene names.
	- Files for mm10 and hg38 GC Content by gene is found in the options_files folder.  Requires ENSEMBL format annotation if used.
	- ERCC mix file provided in options_files folder.  If using a different control set, mimic this format for use.
5. Copy 'RNAseq_experiment_file.yml' into a newly created folder.
6. Rename your experiment file as needed and edit the file to fit your experiment using a text editor (ie. textEdit).  

#### Run the pipeline:

7. To run the pipeline, use the following code directly in the LSF login-node:

> source activate RNAseq
> python /path/to/RNAseq.py -f /path/to/RNAseq_experimental_file.yml -s -p your_lsf_project

Extra options: 
- Add '-t /path/to/RNAseq.ipynb' if you running as a jupyter notebook
- Add '-o /path/to/output_RNAseq.ipynb' to name your output notebook something other than your experimental file name.
- Add --no-notebook to run as a python script with a log file output.
- For all options: 'python RNAseq.py --help'

8. Alternatively, modify the RNAseq_bsub.sh included in the folder and run:
	
> bsub < RNAseq_bsub.sh

9. In case of error, use the above command to pick up from last completed step if the restart option is 'no' in the experimetnal file.  In notebook mode, restarting from the last completed step will erase all priovous outputs to the notebook and start from the last step.  To avoid this, rename the output notebook.   

Jupyter notebooks can be visualized with nteract (https://nteract.io) if you are not familiar with jupyter notebooks.

## RNAseq.py can be imported to python as a module with the following attributes:
	
Experiment object can be passed to the following functions and returned: exp = function(exp):
- RNAseq.DESeq2(exp)
- RNAseq.fastq_cat(exp)       
- RNAseq.fastq_screen(exp)
- RNAseq.GO_enrich(exp)
- RNAseq.fastqc(exp)
- RNAseq.GSEA(exp)
- RNAseq.final_qc(exp)
- RNAseq.principal_component_analysis(exp)
- RNAseq.star(exp)
- RNAseq.rsem(exp)
- RNAseq.Sleuth(exp)
- RNAseq.kallisto(exp)
- RNAseq.sigs(exp)
- RNAseq.clustermap(exp)
- RNAseq.spike(exp)
- RNAseq.overlaps(exp)
- RNAseq.stage(exp)    
- RNAseq.trim(exp)
- RNAseq.GC_normalization(exp)

- RNAseq.pipeline(experimental_file)

#### Helper functions:
- RNAseq.parse_yaml(experimental_file) takes required yaml file and parses an new experimental object or loads an incomplete one from scratch     
- RNAseq.plot_venn2(Series=pd.Series(), string_name_of_overlap='', folder='') scaled venn of an overlap
- RNAseq.bigwig2bam(in_bam='/path/to/bam',out_bw='/path/to/bw_name',job_log_folder=exp.job_folder,sample='',project=exp.project, stranded=bool) scales bam to bigwig cpm.
- RNAseq.job_wait(id_list=exp.job_id, job_log_folder=exp.job_folder, log_file=exp.log_file) waits for submitted job to finish
- RNAseq.send_job(command_list=[], job_name='', job_log_folder=exp.job_folder, q='', mem='', log_file=exp.log_file) sends job to LSF resource manager
- RNAseq.enrichr(gene_list=[], description='', out_dir='')
- RNAseq.plot_PCA(counts=pd.Dataframe(), colData=pd.Dataframe(), out_dir='', name='')
- RNAseq.volcano(results=pd.Dataframe(DESeq2_results), sig_up=[], sig_down=[], name='', out_dir='')
- RNAseq.gseq_barplot(out_dir,pos_file,neg_file,gmt_name,max_number)
- RNAseq.RUV(data=pd.DataFrame(counts),design='~',colData=pd.DataFrame(),type='ercc'or'empirical',log=exp.log_file, ERCC_counts=pd.DataFrame(), comparison='', plot_dir='')
- RNAseq.plot_exp(data=pd.DataFrame(), plot_dir='', exp_type='', name='')
- RNAseq.rout_write()
- RNAseq.html_header()
- RNAseq.val_folder(folder)
- RNAseq.read_pd(file)
- RNAseq.output(text,log_file)
- RNAseq.image_display(file)
- RNAseq.out_result(image,text)
- RNAseq.validated_run()
- RNAseq.quartile_norm()
