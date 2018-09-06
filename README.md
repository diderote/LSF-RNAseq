# Nimerlab RNAseq pipeline

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
8. Between sample normalization options: median-ratios (DESeq2 default), or removal of unwated variation (RUVSeq) using ERCC spike-ins or empirical negative controls.
9. Differential expression using DESeq2 (both default and apeglm or ashr shrunken LFC are reported)
10. PCA plots and expression boxplots of all and experimental samples, raw and normalized counts
11. Optional overlap with signifiantly differentially expressed genes by DESeq2 (q<0.05, 2FC, 1.5FC, no FC filter) with Sleuth (Kallisto) q<0.05 DE genes 
12. Volcano plots generated from DESeq2 results with signifcant genes
13. Heatmap of significantly differentially expressed genes using regularized log2 expresion counts.
14. GO and KEGG enrichment of DE genes
15. GSEA using ranked gene lists from DESeq2
16. Overlap of multiple gene lists (including from differential expression) and scaled venn diagram output and GO enrichment.

Option Details:
* Restart: (yes/no) Whether to check for an incomplete pipeline and pickup where left off, or restart from scratch.
* ERCC_spike: align reads to spike index using STAR.
* Stranded: (yes/no) Will generate two stranded signal bigwig files if stranded.  Otherwise, will generate one per sample.
* Normalization: 
	* Median-Ratios = default DESeq2 method normalization: median of ratios of observed counts. Also reports apeglm lfc shrunken values.
	* ERCC = Account for unwanted variation between samples using spike in counts. (RUVSeq implementation).  Also reports ashr lfc shrunken values.
	* ERCC_Mixed = Same as ERCC except used when ERCC Mix1 and Mix2 are used in the same experiment.  Will reduce the unwanted variance estimation to subgroup B (common concentrations in both mixes).
	* Empirical = Account for unwanted variation between samples using non-differentially expressed genes to estimate unwanted variance (RUVSeq implementation). Also reports ashr lfc shrunken values.
* Alignment_Mode: Whether to align to transcriptome with RSEM-STAR (default) or to genome with STAR.
* Signature_Mode: Output signature will either be DESeq2 differentially expressed genes or an overlap of Slueth and DESeq2.
* GC_Normalizaiton: Yes implements within-lane loess normalzation based on gene GC content (EDASeq).  Recommended for samples sequenced in different sequencing runs.  If not 'Nimer' this adds over an hour of computing time for CG content file generatation.
* Sequencing_type: paired or single end sequencing.

Other options:
*  Spike_matrix: if not aligning and using ERCC, path to file containing ERCCs counts.
*  Count_matrix: if not aligning, path to a file containing gene names as rows(first column), sample names as columns(top row), and gene counts in cells.
*  Sig_matrix: if solely performing overlaps, path to file with samples names (top row) and gene names to overlap below. Overlap of over two columns will be performed (1v2, 3v4... etc)

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
* Gene_names: optional. provide path to a pickled dictionary to add an attribute such as 'gene_name' to the count_matrix.  Useful for STAR alinging where gene_id is default.  

The pipeline handles multiple entry/exit points and can parse complex experimental designs and compensation types for DE.  In case of error, the pipeline restarts from the last completed step. Progress is tracked in a .log file in the output directory.

All submission scripts, error and output files are saved.

## Instructions for installation and use:

### This pipeline is compatable with UM Pegasus.

1. Download https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh to your nethome folder.
2. Install miniconda: 'bash ~/Miniconda3-latest-Linux-x86_64.sh'
3. After installation type the following commands:

> module rm python
> conda env create -f /projects/ctsi/nimerlab/DANIEL/tools/nimerlab-pipelines/RNAseq/environment.yml

 - If this throws errors, check that miniconda is installed in $HOME/miniconda3/ or change the bottom line of the environment.yml file to match your conda installation.  If there is a package error, try removing the version of that package from the .yml file and trying again.

4. OPTIONAL SETUP:
	- For use of contamination screen (fastq screen):
		- generate the relevant bowtie2 indices and add the path were indicated in options_files/fastq_screen.conf. (many premade by illumina at https://support.illumina.com/sequencing/sequencing_software/igenome.html)
		- unhash relvant lines before each desired DATABASE
		- cp fastq_screen.conf ~/miniconda3/envs/RNAseq/share/fastq-screen-0.11.3-0/
	- Files for mm10 and hg38 GC Content by gene is found in the options_files folder.  Requires ENSEMBL format annotation if used.
	- ERCC mix file provided in options_files folder.  If using a different control set, mimic this format for use.
5. Copy '/projects/ctsi/nimerlab/DANIEL/tools/nimerlab-pipelines/RNAseq/RNAseq_experiment_file.yml' into your run folder.
6. Rename your experiment file as needed and edit the file to fit your experiment using a text editor (ie. textEdit).
7. Copy '/projects/ctsi/nimerlab/DANIEL/tools/nimerlab-pipelines/RNAseq/RNAseq' into your run folder.
8. From your run folder, run analysis with this command (replacing 'RNAseq_expiermental_file.yml' with your experimental filename:
	
> bsub -q general -n 1 -R 'rusage[mem=3000]' -W 120:00 -o RNAseq.out -e RNAseq.err -P <project> <<< 'module rm python share-rpms65;source activate RNAseq;./RNAseq.py -f RNAseq_experimental_file.yml' 

9. In case of error, use the above command to pick up from last completed step.  Until the pipeline is complete, all files are stored and can be accessed in the scratch folder.

### RNAseq.py can be imported to python as a module with the following attributes:
	
Experiment object can be passed to the following functions and returned: exp = function(exp):
- RNAseq.DESeq2()
- RNAseq.fastq_cat()       
- RNAseq.fastq_screen()
- RNAseq.GO_enrich()
- RNAseq.fastqc()
- RNAseq.GSEA()
- RNAseq.final_qc()
- RNAseq.principal_component_analysis()
- RNAseq.star()
- RNAseq.rsem()
- RNAseq.Sleuth()
- RNAseq.kallisto()
- RNAseq.sigs()
- RNAseq.clustermap()
- RNAseq.spike()
- RNAseq.splicing()
- RNAseq.overlaps()
- RNAseq.stage()    
- RNAseq.trim()
- RNAseq.GC_normalization()

#### Helper functions:
- RNAseq.parse_yaml() takes required -f yaml file and parses an new experimental object or loads an incomplete one     
- RNAseq.plot_venn2(Series=pd.Series(), string_name_of_overlap='', folder='') scaled venn of an overlap
- RNAseq.bigwig2bam(in_bam='/path/to/bam',out_bw='/path/to/bw_name',job_log_folder=exp.job_folder,sample='',project=exp.project, stranded=bool) scales bam to bigwig cpm.
- RNAseq.job_wait(id_list=exp.job_id, job_log_folder=exp.job_folder, log_file=exp.log_file) waits for submitted job to finish
- RNAseq.send_job(command_list=[], job_name='', job_log_folder=exp.job_folder, q='', mem='', log_file=exp.log_file) sends job to LSF resource manager
- RNAseq.enrichr(gene_list=[], description='', out_dir='')
- RNAseq.plot_PCA(counts=pd.Dataframe(), colData=pd.Dataframe(), out_dir='', name='')
- RNAseq.volcano(results=pd.Dataframe(DESeq2_results), sig_up=[], sig_down=[], name='', out_dir='')
- RNAseq.RUV(data=pd.DataFrame(counts),design='~',colData=pd.DataFrame(),type='ercc'or'empirical',log=exp.log_file, ERCC_counts=pd.DataFrame(), comparison='', plot_dir='')
- RNAseq.plot_exp(data=pd.DataFrame(), plot_dir='', exp_type='', name='')
- RNAseq.rout_write()
