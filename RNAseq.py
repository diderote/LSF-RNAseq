#!/usr/bin/env python3

'''
University of Miami - Pegasus LSF Cluster RNASeq Pipeline

Reads an experimental design yaml file (Version 0.7).
Requires a conda environment 'RNAseq' made from environment.yml
www.github.com/diderote/Nimerlab-RNAseq/

To do:
    - change designs to allow for non-bindary conditions ('Condition_A: a,a,b,b,c,c')
    - ICA with chi-square with de groups (start as year random)
    - add sleuth for mouse and make compatable with new v.7 strategy
    - optimize STAR
    - add summary statistics at finish
        -fragment length
        -mapping ..etc
    - fine tune DE tests for lfcshrink or GC norm option per test
    - if restarting with papermill... make different ipynb
    - add rlog with design (one blind, other full design)
    - tsne with several perplexities?
    - add single cell (deseq2, min to replace=Inf, other sc options zeroinflation?, LRT)
    - extract log2 instead of rlog with RUVg normcounts or EDAseq normcounts.  ?

'''
import os
import re
import glob
import pickle
import random
import time
from shutil import copy2, copytree, rmtree, move
from datetime import datetime
import yaml
import reprlib

from IPython.display import HTML, display, Image
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib_venn import venn2, venn2_circles
import seaborn as sns
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import rpy2.robjects as ro
import rpy2.rinterface as ri
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from scipy import stats
import gseapy
# import PyPDF2

__author__ = 'Daniel L. Karl'
__license__ = 'MIT'
__version__ = '0.7'

run_main = True if __name__ == '__main__' else False


class Experiment:
    '''
    Experiment object for pipeline
    '''
    def __init__(self):
        self.norm = 'Median-Ratios'
        self.de_sig_overlap = False
        self.gc_norm = False
        self.trim = [0, 0]
        self.tasks_complete = []
        self.job_id = []
        self.designs = {}
        self.overlaps = {}
        self.gene_lists = {}
        self.de_results = {}
        self.sig_lists = {}
        self.overlap_results = {}
        self.genome_indicies = {}

    def __repr__(self):
        exclude = ['overlaps', 'job_id', 'name', 'designs']
        experiment = f'{self.__dict__["name"]}'
        for key, value in self.__dict__.items():
            experiment += f'\n\t{key}: {reprlib.repr(value)}' if key not in exclude else ''
        return f'Experiment({experiment})'


def html_header():
    return ''.join(['<h1>RNAseq Analysis Notebook</h1>',
                    f'<body><b>Experiment Date: {datetime.now():%Y-%m-%d}<br>',
                    f'Pipeline version: {__version__}</b><br>',
                    '<a href="http://www.github.com/diderote/LSF-RNAseq">Pipeline Code</a><br>',
                    f'License: MIT <br> Author: {__author__}'
                    ])


def val_folder(folder):
    folder = folder if folder.endswith('/') else f'{folder}/'
    return '' if folder == '/' else folder


def rout_write(rout):
    '''
    function for setting r_out to print to file
    '''
    print(rout, file=open(f'{os.getcwd()}/R_out_{datetime.now():%Y-%m-%d}.txt', 'a'))


def read_pd(file):
    if (file.split('.')[-1] == 'txt') or (file.split('.')[-1] == 'tab'):
        return pd.read_table(file, header=0, index_col=0)
    elif (file.split('.')[-1] == 'xls') or (file.split('.')[-1] == 'xlsx'):
        return pd.read_excel(file)
    else:
        raise IOError("Cannot parse file.  Make sure it is .txt, .xls, or .xlsx")


def output(text, log_file):
    if run_main:
        print(text, file=open(log_file, 'a'))
    else:
        print(text)


def image_display(file):
    display(Image(file))


def out_result(image, text):
    if not run_main:
        if os.path.isfile(image):
            display(HTML(f'<h2>{text}</h2>'))
            image_display(image)
        else:
            display(HTML(f'<body>No result for {text} found.</body>'))


def parse_yaml(experimental_file):
    '''
    Parse experimental info from yaml file
    '''

    with open(experimental_file, 'r') as file:
        yml = yaml.safe_load(file)

    # Make a new experimental object
    exp = Experiment()

    # Setting Scratch folder
    if yml['Pegasus_Project'] == 'nimerlab':
        exp.scratch = f'/scratch/projects/nimerlab/DANIEL/staging/RNAseq/{yml["Name"]}/'
    else:
        exp.scratch = f'{val_folder(yml["Scratch_folder"])}{yml["Name"]}/'
    os.makedirs(exp.scratch, exist_ok=True)

    # check whether experiment has been attempted
    exp.name = yml['Name']
    filename = f'{exp.scratch}{exp.name}_incomplete.pkl'

    if os.path.isfile(filename):
        if not yml['Restart']:
            with open(filename, 'rb') as experiment:
                exp = pickle.load(experiment)
            os.remove(filename)

            # set new date
            exp.date = f'{datetime.now():%Y-%m-%d}'

            # For output of R logs into job_log_folder
            os.chdir(exp.job_folder)

            output(f'\n{"#" * 30}\nRestarting pipeline on {datetime.now():%Y-%m-%d %H:%M:%S}, from last completed step.', exp.log_file)
            output(f'\nExperimental variables:\n {exp}', exp.log_file)

            return exp
        else:
            os.remove(filename)

    # Passing paramters to new object
    exp.date = f'{datetime.now():%Y-%m-%d}'

    # Make out directory if it doesn't exist
    exp.out_dir = f'{val_folder(yml["Output_directory"])}{exp.name}/'
    os.makedirs(exp.out_dir, exist_ok=True)

    # Log file
    exp.log_file = f'{exp.out_dir}{exp.name}-{exp.date}.log'

    output(f'Pipeline version {str(__version__)} run on {exp.date} \n', exp.log_file)
    output(f'Beginning RNAseq Analysis: {datetime.now():%Y-%m-%d %H:%M:%S}\n', exp.log_file)
    output('Reading experimental file...\n', exp.log_file)
    output(f"Pipeline output folder: {exp.out_dir}\n", exp.log_file)

    # Setting Job Folder
    exp.job_folder = f'{exp.scratch}logs/'
    os.makedirs(exp.job_folder, exist_ok=True)

    # Genome
    exp.genome = yml['Genome'].lower() if yml['Genome'].lower() in ['hg38', 'mm10', 'hg19'] else None
    if exp.genome is None:
        raise ValueError("Genome must be either hg38, hg19, or mm10.")
    output(f'Processing data with: {exp.genome}', exp.log_file)

    # Tasks to complete
    if not yml['Tasks']['Align']:
        exp.tasks_complete = exp.tasks_complete + ['Stage', 'FastQC', 'Fastq_screen', 'Trim', 'STAR', 'RSEM', 'Kallisto', 'Sleuth']
        output('Not performing alignment.', exp.log_file)
        count_matrix_loc = yml['Count_matrix']
        if os.path.exists(count_matrix_loc):
            output(f"Count matrix found at {count_matrix_loc}", exp.log_file)
            output("Performing only DESeq2 on for DE", exp.log_file)
            exp.count_matrix = read_pd(count_matrix_loc)
        else:
            raise IOError("Count Matrix Not Found.")
    elif yml['Tasks']['Align']:
        # Alignment mode. Default is transcript.
        exp.alignment_mode = 'gene' if yml['Alignment_Mode'].lower() == 'gene' else 'transcript'
        if exp.alignment_mode == 'gene':
            exp.tasks_complete.append('RSEM')
        elif exp.alignment_mode == 'transcript':
            exp.tasks_complete.append('STAR')

        # Sequencing type
        if yml['Sequencing_type'].lower() not in ['paired', 'single']:
            raise ValueError("Must specify whether sequence is paired or single end.")
        exp.seq_type = yml['Sequencing_type'].lower()
        exp.tasks_complete = exp.tasks_complete + ['Kallisto', 'Sleuth'] if exp.seq_type == 'single' else exp.tasks_complete

        # Standed
        exp.stranded = True if yml['Stranded'] else False
        output(f'Processing data as {exp.seq_type}-end {("stranded" if exp.stranded else "non-stranded")} sequencing.', exp.log_file)

        # Sequencer for trimming options
        exp.sequencer = yml['Sequencer'].lower() if yml['Sequencer'].lower() in ['nextseq', 'hiseq'] else None
        output(f'Sequencer: {exp.sequencer.capitalize()}', exp.log_file)
    else:
        raise IOError('Please specify whether or not to perform alignment.')

    # Lab specific files
    if yml['Pegasus_Project'] == 'nimerlab':
        exp.genome_indicies['ERCC'] = '/projects/ctsi/nimerlab/DANIEL/tools/genomes/ERCC_spike/STARIndex'
        exp.genome_indicies['GSEA_jar'] = '/projects/ctsi/nimerlab/DANIEL/tools/GSEA/gsea-3.0.jar'
        if exp.genome == 'mm10':
            exp.genome_indicies['RSEM_STAR'] = '/projects/ctsi/nimerlab/DANIEL/tools/genomes/Mus_musculus/mm10/RSEM-STARIndex/mouse'
            exp.genome_indicies['STAR'] = '/projects/ctsi/nimerlab/DANIEL/tools/genomes/Mus_musculus/mm10/STARIndex'
            exp.genome_indicies['Kallisto'] = '/projects/ctsi/nimerlab/DANIEL/tools/genomes/Mus_musculus/mm10/KallistoIndex/GRCm38.transcripts.idx'
            exp.genome_indicies['Gene_names'] = '/projects/ctsi/nimerlab/DANIEL/tools/genomes/Mus_musculus/mm10/gencode_gene_dict.pkl'
            exp.genome_indicies['GMT'] = '/projects/ctsi/nimerlab/DANIEL/tools/GSEA/mouse_gmts/'
        elif exp.genome == 'hg38':
            exp.genome_indicies['RSEM_STAR'] = '/projects/ctsi/nimerlab/DANIEL/tools/genomes/H_sapiens/NCBI/GRCh38/Sequence/RSEM-STARIndex/human'
            exp.genome_indicies['STAR'] = '/projects/ctsi/nimerlab/DANIEL/tools/genomes/H_sapiens/NCBI/GRCh38/Sequence/STARIndex'
            exp.genome_indicies['Kallisto'] = '/projects/ctsi/nimerlab/DANIEL/tools/genomes/H_sapiens/NCBI/GRCh38/Sequence/KallistoIndex/GRCh38.transcripts.idx'
            exp.genome_indicies['Gene_names'] = '/projects/ctsi/nimerlab/DANIEL/tools/genomes/H_sapiens/NCBI/GRCh38/Annotation/Archives/archive-2015-08-11-09-31-31/Genes.gencode/ENSG_NAME_dict.pkl'
        elif exp.genome == 'hg19':
            exp.genome_indicies['RSEM_STAR'] = '/projects/ctsi/nimerlab/DANIEL/tools/genomes/H_sapiens/Hg19/NCBI-RNAseq/RSEM-STAR/human'
            exp.genome_indicies['STAR'] = '/projects/ctsi/nimerlab/DANIEL/tools/genomes/H_sapiens/Hg19/NCBI-RNAseq/STAR'
            exp.genome_indicies['Kallisto'] = '/projects/ctsi/nimerlab/DANIEL/tools/genomes/H_sapiens/Hg19/NCBI-RNAseq/Kallisto/GRCh37.transcripts.idx'
            exp.genome_indicies['Gene_names'] = '/projects/ctsi/nimerlab/DANIEL/tools/genomes/H_sapiens/Hg19/gencode_gene_dict.pkl'
    else:
        exp.genome_indicies['RSEM_STAR'] = yml['RSEM_STAR_index']
        exp.genome_indicies['STAR'] = yml['STAR_index']
        exp.genome_indicies['Kallisto'] = yml['Kallisto_index']
        exp.genome_indicies['ERCC'] = yml['ERCC_STAR_index']
        exp.genome_indicies['GSEA_jar'] = yml['GSEA_jar']
        exp.genome_indicies['Gene_names'] = yml['Gene_names']
        if exp.genome == 'mm10':
            exp.genome_indicies['GMT'] = val_folder(yml['GSEA_mouse_gmx_folder'])

    # GC_normalizaton
    if yml['GC_Normalization']:
        exp.gc_norm = True
    else:
        exp.tasks_complete.append('GC')

    # Support Files:
    if yml['Pegasus_Project'] == 'nimerlab':
        exp.genome_indicies['ERCC_Mix'] = '/projects/ctsi/nimerlab/DANIEL/tools/genomes/ERCC_spike/cms_095046.txt'
        if exp.genome == 'mm10':
            exp.genome_indicies['GC_Content'] = '/projects/ctsi/nimerlab/DANIEL/tools/genomes/Mus_musculus/mm10/mm10_GC_Content.txt'
        elif exp.genome == 'hg38':
            exp.genome_indicies['GC_Content'] = '/projects/ctsi/nimerlab/DANIEL/tools/genomes/H_sapiens/hg38_GC_Content.txt'
        elif exp.genome == 'hg19':
            exp.genome_indicies['GC_Content'] = '/projects/ctsi/nimerlab/DANIEL/tools/genomes/H_sapiens/hg38_GC_Content.txt'
    else:
        exp.genome_indicies['ERCC_Mix'] = yml['ERCC_Mix_file']
        exp.genome_indicies['GC_Content'] = yml['GC_Content_file']

    # No DE option
    if not yml['Tasks']['Differential_Expression']:
        exp.tasks_complete = exp.tasks_complete + ['GC', 'DESeq2', 'Sleuth', 'Sigs', 'Heatmaps', 'GO_enrich', 'GSEA', 'PCA']
        output('Not performing differential expression analyses.', exp.log_file)

    # Spike
    if not yml['Tasks']['ERCC_align']:
        exp.tasks_complete.append('Spike')

    # Fastq Folder
    exp.fastq_folder = val_folder(yml['Fastq_directory'])
    if (os.path.isdir(exp.fastq_folder) is False) and (yml['Tasks']['Align']):
        raise IOError("Can't Find Fastq Folder.")

    # Hard clip
    # exp.trim = yml['trim'] if yml['trim'] is not None else exp.trim

    # Project
    exp.project = yml['Pegasus_Project']

    # Sample Names
    exp.samples = {key: name for key, name in yml['Samples'].items() if (name != 'Name') and (name is not None)}
    exp.sample_number = len(exp.samples)
    output("Samples: ", exp.log_file)
    for number, sample in exp.samples.items():
        output(f'{number}: {sample}', exp.log_file)
    output(f'\nProcessing {exp.sample_number} samples.\n', exp.log_file)

    # Differential Expression Groups
    if yml['Tasks']['Differential_Expression']:
        output("Parsing experimental design for differential expression...\n", exp.log_file)

        exp.conditions = {key: condition.split(',') for key, condition in yml['Conditions'].items() if condition is not None}
        exp.de_tests = {key: test for key, test in yml['Designs'].items() if test['Test_condition'] is not None}
        for key, test in exp.de_tests.items():
            all_conditions = test['All_conditions'].split(',')
            all_samples = list({exp.samples[int(sample)] for sample in test['All_samples'].split(',')})
            exp.designs[key] = {'all_samples': all_samples}
            for condition in all_conditions:
                exp.designs[key][condition] = list({exp.samples[int(sample)] for sample in exp.conditions[condition]})
            exp.designs[key]['Test_condition'] = test['Test_condition'].split(',')

            if len(exp.designs[key]['Test_condition']) == 1:
                exp.designs[key]['reduced'] = '~' + ' + {}'.join([condition for condition in all_conditions if condition not in exp.designs[key]['Test_condition']])
                exp.designs[key]['design'] = f'{exp.designs[key]["reduced"]} + {exp.designs[key]["Test_condition"][0]}' if len(all_conditions) > 1 else f'~{exp.designs[key]["Test_condition"][0]}'
            elif len(exp.designs[key]['Test_condition']) == 2:
                exp.designs[key]['reduced'] = '~' + ' + '.join([f'{condition}' for condition in all_conditions])
                intersection = f'{exp.designs[key]["Test_condition"][0]}:{exp.designs[key]["Test_condition"][1]}'
                exp.designs[key]['design'] = f'{exp.designs[key]["reduced"]} + {intersection}'
            else:
                raise ValueError('Cannot handle this experimental design.')
            exp.designs[key]['reduced'] = '~1' if exp.designs[key]['reduced'] == '~' else exp.designs[key]['reduced']

            exp.designs[key]['colData'] = pd.DataFrame({f'{condition}': ['yes' if sample in exp.designs[key][f'{condition}'] else 'no' for sample in all_samples] for condition in all_conditions}, index=all_samples)
            exp.designs[key]['Test_type'] = exp.de_tests[key]['Test_type'].lower()

        all_colData = pd.DataFrame(index=[sample for sample in exp.samples.values()])
        for condition in exp.conditions.keys():
            all_colData[condition] = ['yes' if str(sample) in exp.conditions[condition] else 'no' for sample in exp.samples.keys()]

        exp.designs['complete'] = {'design': '~' + ' + '.join([condition for condition in exp.conditions.keys()]),
                                   'colData': all_colData
                                   }

        for name, items in exp.designs.items():
            output(f'\n{name}:', exp.log_file)
            output(str(items['colData']), exp.log_file)

        # Normalization method
        if yml['Normalization'].lower() == 'ercc' or yml['Normalization'].lower() == 'ercc_mixed':
            exp.norm = yml['Normalization'].lower()
            if not yml['Tasks']['ERCC_align']:
                spike_matrix_loc = yml['Spike_matrix']
                if os.path.exists(spike_matrix_loc):
                    output(f"Spike Count matrix found at {spike_matrix_loc}", exp.log_file)
                    exp.spike_counts = read_pd(spike_matrix_loc)
                else:
                    output("Cannot find spike matrix.", exp.log_file)
            output('\nNormalizing samples for differential expression analysis using ERCC spike-in variance.\n', exp.log_file)
        elif yml['Normalization'].lower() == 'empirical':
            output('\nNormalizing samples for differential expression analysis using empirical negative controls for variance.\n', exp.log_file)
            exp.norm = 'empirical'
        elif yml['Normalization'].lower() == 'median-ratios':
            output('\nNormalizing samples for differential expression analysis using deseq2 size factors determined using default median of ratios method.\n', exp.log_file)
        else:
            output(f"\nI don't know the {yml['Normalization']} normalization method.  Using default median-ratios.\n", exp.log_file)

        exp.lfcshrink = True if yml['LFCshrink'] else False
        if exp.lfcshrink:
            output("\nReporting additional file with log fold change shrinkage for differential expression tests.", exp.log_file)

    # Initialize DE sig overlaps
    exp.de_sig_overlap = True if yml['Signature_Mode'].lower() == 'combined' else False
    if exp.de_sig_overlap is False or exp.alignment_mode.lower() == 'transcript':
        exp.tasks_complete = exp.tasks_complete + ['Kallisto', 'Sleuth']

    # DE Overlaps
    for key, item in yml['Overlaps'].items():
        if bool(item):
            exp.overlaps[key] = item.split(':')
    if str(len(list(exp.overlaps.keys()))) != 0:
        output(f'\nOverlapping {len(list(exp.overlaps.keys()))} differential analysis comparison(s).', exp.log_file)
        output(f'{exp.overlaps}\n', exp.log_file)
    else:
        exp.tasks_complete.append('Overlaps')
        output('Not performing signature overlaps', exp.log_file)

    # Initialized Process Complete List
    exp.tasks_complete.append('Parsed')

    # For output of R logs into job_log_folder
    os.chdir(exp.job_folder)

    output(f'Experiment file parsed: {datetime.now():%Y-%m-%d %H:%M:%S}\n', exp.log_file)

    return exp


def send_job(command_list, job_name, job_log_folder, q, mem, log_file, project, cores=1):
    '''
    Sends job to LSF pegasus.ccs.miami.edu
    '''

    os.makedirs(job_log_folder, exist_ok=True)

    rand_id = str(random.randint(0, 100000))
    str_comd_list = '\n'.join(command_list)
    cmd = f'''
#!/bin/bash

#BSUB -J JOB_{job_name}_ID_{rand_id}
#BSUB -R "rusage[mem={mem}]"
#BSUB -R "span[ptile={cores}]"
#BSUB -o {job_log_folder}{job_name}_logs_{rand_id}.stdout.%J
#BSUB -e {job_log_folder}{job_name}_logs_{rand_id}.stderr.%J
#BSUB -W 120:00
#BSUB -n {cores}
#BSUB -q {q}
#BSUB -P {project}

{str_comd_list}'''

    job_path_name = f'{job_log_folder}{job_name}.sh'
    write_job = open(job_path_name, 'w')
    write_job.write(cmd)
    write_job.close()
    os.system(f'bsub < {job_path_name}')
    output(f'sending {job_name} as ID_{rand_id}...', log_file)
    time.sleep(2)  # too many conda activations at once sometimes leads to inability to activate during a job.

    return rand_id


def job_wait(id_list, log_file):
    '''
    Waits for jobs sent by send job to finish.
    '''
    waiting = True
    while waiting:
        with os.popen('bhist -w') as stream:
            jobs_list = stream.read()
        current = []
        for rand_id in id_list:
            if len([j for j in re.findall(r'ID_(\d+)', jobs_list) if j == rand_id]) != 0:
                current.append(rand_id)
        if len(current) == 0:
            waiting = False
        else:
            output(f'Waiting for jobs to finish... {datetime.now():%Y-%m-%d %H:%M:%S}', log_file)
            time.sleep(60)


def stage(exp):
    '''
    Stages files in Pegasus Scratch
    '''

    # Stage Experiment Folder in Scratch
    output(f'Staging in {exp.scratch}\n', exp.log_file)

    # Copy Fastq to scratch fastq folder
    if os.path.exists(f'{exp.scratch}Fastq'):
        rmtree(f'{exp.scratch}Fastq')

    os.makedirs(f'{exp.scratch}Fastq', exist_ok=True)

    files = [file for file in glob.glob(f'{exp.fastq_folder}*') for sample in exp.samples.values() if file.endswith('.fastq.gz') & (sample in file)]
    for file in files:
        copy2(file, f'{exp.scratch}Fastq')

    exp.fastq_folder = f'{exp.scratch}Fastq/'
    exp.tasks_complete.append('Stage')

    output(f'Staging complete: {datetime.now():%Y-%m-%d %H:%M:%S}\n', exp.log_file)

    return exp


def fastqc(exp):
    '''
    Performs fastq spec analysis with FastQC
    '''
    output('Assessing fastq quality. \n', exp.log_file)

    # Make QC folder
    exp.qc_folder = f'{exp.scratch}QC/'
    os.makedirs(exp.qc_folder, exist_ok=True)

    for number, sample in exp.samples.items():
        command_list = ['module rm python',
                        'module rm perl',
                        'source activate RNAseq',
                        f'fastqc {exp.fastq_folder}{exp.fastq_folder} *'
                        ]

        exp.job_id.append(send_job(command_list=command_list,
                                   job_name=f'{sample}_fastqc',
                                   job_log_folder=exp.job_folder,
                                   q='general',
                                   mem=5000,
                                   log_file=exp.log_file,
                                   project=exp.project
                                   ))

    # Wait for jobs to finish
    job_wait(exp.job_id, exp.log_file)

    # move to qc folder
    fastqc_files = glob.glob(f'{exp.fastq_folder}*.zip')
    fastqc_files = fastqc_files + glob.glob(f'{exp.fastq_folder}*.html')
    for f in fastqc_files:
        copy2(f, exp.qc_folder)
        os.remove(f)

    exp.tasks_complete.append('FastQC')

    output(f'FastQC complete: {datetime.now():%Y-%m-%d %H:%M:%S}\n', exp.log_file)

    return exp


def fastq_screen(exp):
    '''
    Checks fastq files for contamination with alternative genomes using Bowtie2
    '''

    output(f'Screening for contamination during sequencing: {datetime.now():%Y-%m-%d %H:%M:%S}\n', exp.log_file)

    # Make QC folder
    exp.qc_folder = f'{exp.scratch}QC/'
    os.makedirs(exp.qc_folder, exist_ok=True)
    cwd = val_folder(os.getcwd())
    os.chdir(exp.fastq_folder)

    if exp.seq_type == 'paired':
        fastq_end = '_R1.fastq.gz'
    elif exp.seq_type == 'single':
        fastq_end = '.fastq.gz'

    # Submit fastqc and fastq_screen jobs for each sample
    for number, sample in exp.samples.items():
        command_list = ['module rm python',
                        'module rm perl',
                        'source activate RNAseq',
                        f'fastq_screen --threads 4 --aligner bowtie2 {exp.fastq_folder}{sample}{fastq_end}'
                        ]

        exp.job_id.append(send_job(command_list=command_list,
                                   job_name=f'{sample}_fastq_screen',
                                   job_log_folder=exp.job_folder,
                                   q='general',
                                   mem=3000,
                                   log_file=exp.log_file,
                                   project=exp.project,
                                   cores=2
                                   ))
        time.sleep(1)

    # Wait for jobs to finish
    job_wait(exp.job_id, exp.log_file)

    # move to qc folder
    fastqs_files = glob.glob(f'{exp.fastq_folder}*screen*')
    for f in fastqs_files:
        copy2(f, exp.qc_folder)
        os.remove(f)

    # change to experimental directory in scratch
    os.chdir(cwd)

    exp.tasks_complete.append('Fastq_screen')
    output(f'Screening complete: {datetime.now():%Y-%m-%d %H:%M:%S}\n', exp.log_file)

    return exp


def trim(exp):
    '''
    Trimming based on standard UM SCCC Core Nextseq 500 technical errors.  Cudadapt can hard clip both ends, but may ignore 3' in future.
    '''

    output(f'Beginning fastq trimming: {datetime.now():%Y-%m-%d %H:%M:%S}\n', exp.log_file)

    scan = 0
    while scan < 2:

        quality = '--nextseq-trim=20' if exp.sequencer == 'nextseq' else '-q 20'

        # Submit trimming files for each sample
        for number, sample in exp.samples.items():

            move_on = True
            if exp.seq_type == 'paired':
                if f'{exp.fastq_folder}{sample}_trim_R2.fastq.gz' in glob.glob(f'{exp.fastq_folder}*.gz'):
                    move_on = False
            elif exp.seq_type == 'single':
                if f'{exp.fastq_folder}{sample}_trim.fastq.gz' in glob.glob(f'{exp.fastq_folder}*.gz'):
                    move_on = False

            if move_on:
                output(f'Trimming {sample}: ', exp.log_file)

                if exp.seq_type == 'paired':
                    trim_u = str(exp.trim[0])
                    trim_U = str(exp.trim[1])
                    cutadapt = f'cutadapt -j 4 -a AGATCGGAAGAGC -A AGATCGGAAGAGC --cores=10 {quality} -u {trim_u} -u -{trim_u} -U {trim_U} -U -{trim_U} -m 18 -o {exp.fastq_folder}{sample}_trim_R1.fastq.gz -p {exp.fastq_folder}{sample}_trim_R2.fastq.gz {exp.fastq_folder}{sample}_R1.fastq.gz {exp.fastq_folder}{sample}_R2.fastq.gz'
                elif exp.seq_type == 'single':
                    cutadapt = f'cutadapt -j 4 -a AGATCGGAAGAGC --cores=10 {quality} -m 18 -o {exp.fastq_folder}{sample}_trim.fastq.gz {exp.fastq_folder}{sample}.fastq.gz'

                command_list = ['module rm python',
                                'module rm perl',
                                'source activate RNAseq',
                                cutadapt
                                ]

                exp.job_id.append(send_job(command_list=command_list,
                                           job_name=f"{sample}_trim",
                                           job_log_folder=exp.job_folder,
                                           q='general',
                                           mem=1000,
                                           log_file=exp.log_file,
                                           project=exp.project,
                                           cores=2
                                           ))

        # Wait for jobs to finish
        job_wait(exp.job_id, exp.log_file)

        scan += 1

    # move logs to qc folder
    output('\nTrimming logs are found in stdout files from bsub.  Cutadapt does not handle log files in multi-core mode.', exp.log_file)

    for number, sample in exp.samples.items():
        if exp.seq_type == 'paired':
            if f'{exp.fastq_folder}{sample}_trim_R2.fastq.gz' not in glob.glob(f'{exp.fastq_folder}*.gz'):
                raise ValueError('Not all samples were trimmed.')
        elif exp.seq_type == 'single':
            if f'{exp.fastq_folder}{sample}_trim.fastq.gz' not in glob.glob(f'{exp.fastq_folder}*.gz'):
                raise ValueError('Not all samples were trimmed.')

    exp.tasks_complete.append('Trim')
    output(f'Trimming complete: {datetime.now():%Y-%m-%d %H:%M:%S}\n', exp.log_file)

    return exp


def spike(exp):
    '''
    If calling from jupyter.

    Align sequencing files to ERCC index using STAR aligner.
    '''

    plt.clf()
    output(f"Processing with ERCC spike-in: {datetime.now():%Y-%m-%d %H:%M:%S}\n", exp.log_file)

    ERCC_folder = f'{exp.scratch}ERCC/'
    os.makedirs(ERCC_folder, exist_ok=True)

    if not os.path.exists(f'{exp.scratch}Fastq'):
        copytree(exp.fastq_folder, f'{exp.scratch}Fastq')
        exp.fastq_folder = f'{exp.scratch}Fastq/'

    scan = 0
    while scan < 2:
        for number, sample in exp.samples.items():
            # Scan if succesful during second loop.
            if f'{ERCC_folder}{sample}_ERCCReadsPerGene.out.tab' not in glob.glob(f'{ERCC_folder}*.tab'):
                # Submit STAR alingment for spike-ins for each sample
                output(f'Aligning {sample} to spike-in.\n', exp.log_file)

                if os.path.isfile(f'{exp.fastq_folder}{sample}_trim_R1.fastq.gz') or os.path.isfile(f'{exp.fastq_folder}{sample}_trim.fastq.gz'):
                    fname = f'{exp.fastq_folder}{sample}_trim'
                elif os.path.isfile(f'{exp.fastq_folder}{sample}_R1.fastq.gz') or os.path.isfile(f'{exp.fastq_folder}{sample}.fastq.gz'):
                    fname = f'{exp.fastq_folder}{sample}'
                else:
                    output('Cannot find fastq files for spike-in alignment. \n', exp.log_file)
                    raise IOError('Cannot find fastq files for spike-in alignment.')

                if exp.seq_type == 'paired':
                    spike = f'STAR --runThreadN 4 --genomeDir {exp.genome_indicies["ERCC"]} --readFilesIn {fname}_R1.fastq.gz {fname}_R2.fastq.gz --readFilesCommand zcat --outFileNamePrefix {ERCC_folder}{sample}_ERCC --quantMode GeneCounts'
                elif exp.seq_type == 'single':
                    spike = f'STAR --runThreadN 4 --genomeDir {exp.genome_indicies["ERCC"]} --readFilesIn {fname}.fastq.gz --readFilesCommand zcat --outFileNamePrefix {ERCC_folder}{sample}_ERCC --quantMode GeneCounts'

                command_list = ['module rm python',
                                'module rm perl',
                                'source activate RNAseq',
                                spike
                                ]

                exp.job_id.append(send_job(command_list=command_list,
                                           job_name=f'{sample}_ERCC',
                                           job_log_folder=exp.job_folder,
                                           q='general',
                                           mem=5000,
                                           log_file=exp.log_file,
                                           project=exp.project,
                                           cores=2
                                           ))

        # Wait for jobs to finish
        job_wait(exp.job_id, exp.log_file)

        scan += 1

    for number, sample in exp.samples.items():
        sam_file = f'{ERCC_folder}{sample}_ERCCAligned.out.sam'
        if os.path.isfile(sam_file):
            os.remove(sam_file)

    output('Spike-in alignment jobs finished.', exp.log_file)

    # Generate one matrix for all spike_counts
    try:
        ERCC_counts = glob.glob(f'{ERCC_folder}*_ERCCReadsPerGene.out.tab')
        if len(ERCC_counts) != exp.sample_number:
            output('At least one ERCC alignment failed.', exp.log_file)
            raise RuntimeError('At least one ERCC alignment failed. Check scripts and resubmit.')
        else:
            exp.spike_counts = pd.DataFrame(index=pd.read_csv(ERCC_counts[1], header=None, index_col=0, sep="\t").index)

            for number, sample in exp.samples.items():
                exp.spike_counts[sample] = pd.read_csv(f'{ERCC_folder}{sample}_ERCCReadsPerGene.out.tab', header=None, index_col=0, sep="\t")[[3]]
            exp.spike_counts = exp.spike_counts.iloc[4:, :]
            exp.spike_counts.to_csv(f'{ERCC_folder}ERCC.count.matrix.txt', header=True, index=True, sep="\t")

    except:
        output('Error generating spike_count matrix.', exp.log_file)
        raise RuntimeError('Error generating spike_count matrix. Make sure the file is not empty.')

    # check to see if there were any spike in reads, if not, change
    if exp.spike_counts.loc['ERCC-00002', :].sum(axis=0) < 50:
        output('ERCC has low or no counts, skipping further spike-in analysis.', exp.log_file)
        return exp

    if exp.genome_indicies['ERCC_Mix'] is not None:
        # Filtering for counts with more than 5 counts in two samples
        spike_counts = exp.spike_counts.copy()
        spike_counts = spike_counts[spike_counts[spike_counts > 5].apply(lambda x: len(x.dropna()) > 1, axis=1)]
        mix = pd.read_csv(exp.genome_indicies['ERCC_Mix'], header=0, index_col=1, sep="\t")
        mix = mix.rename(columns={'concentration in Mix 1 (attomoles/ul)': 'Mix_1',
                                  'concentration in Mix 2 (attomoles/ul)': 'Mix_2'})
        names = list(spike_counts.columns)
        spike_counts = spike_counts.join(mix)

        merged_spike = pd.DataFrame(columns=['value', 'Mix_1', 'Mix_2', 'subgroup'])
        name = []
        length = len(spike_counts)
        for sample in names:
            merged_spike = pd.concat([merged_spike,
                                     spike_counts[[sample, 'Mix_1', 'Mix_2', 'subgroup']].rename(columns={sample: 'value'})],
                                     ignore_index=True)
            name = name + [sample] * length
        merged_spike['Sample'] = name
        merged_spike['log'] = merged_spike.value.apply(lambda x: np.log2(x))
        merged_spike['log2_Mix_1'] = np.log2(merged_spike.Mix_1)
        merged_spike['log2_Mix_2'] = np.log2(merged_spike.Mix_2)

        # Plot ERCC spike.
        plt.clf()
        sns.set(context='paper', font_scale=2, style='white', rc={'figure.dpi': 300, 'figure.figsize': (6, 6)})
        M1 = sns.lmplot(x='log2_Mix_1', y='log', hue='Sample', data=merged_spike, size=10, aspect=1, ci=None)
        M1.set_ylabels(label='spike-in counts (log2)')
        M1.set_xlabels(label='ERCC Mix (log2(attamoles/ul))')
        plt.title("ERCC Mix 1 Counts per Sample")
        sns.despine()
        M1.savefig(f'{ERCC_folder}ERCC_Mix_1_plot.png')
        out_result(f'{ERCC_folder}ERCC_Mix_1_plot.png', 'ERCC Mix1 Plot')
        if run_main:
            plt.close()

        plt.clf()
        sns.set(context='paper', font_scale=2, style='white', rc={'figure.dpi': 300, 'figure.figsize': (6, 6)})
        M2 = sns.lmplot(x='log2_Mix_2', y='log', hue='Sample', data=merged_spike, size=10, aspect=1, ci=None)
        M2.set_ylabels(label='spike-in counts (log2)')
        M2.set_xlabels(label='ERCC Mix (log2(attamoles/ul))')
        plt.title("ERCC Mix 2 Counts per Sample")
        sns.despine()
        M2.savefig(f'{ERCC_folder}ERCC_Mix_2_plot.png')
        out_result(f'{ERCC_folder}ERCC_Mix_2_plot.png', 'ERCC Mix2 Plot')
        if run_main:
            plt.close()

        plt.clf()
        sns.set(context='paper', font_scale=2, style='white', rc={'figure.dpi': 300, 'figure.figsize': (6, 6)})
        setB = sns.lmplot(x='log2_Mix_1', y='log', hue='Sample', data=merged_spike[merged_spike.subgroup == 'B'], size=10, aspect=1, ci=None)
        setB.set_ylabels(label='spike-in counts (log2)')
        setB.set_xlabels(label='ERCC Mix (log2(attamoles/ul))')
        plt.title("ERCC Subgroup B Counts per Sample")
        sns.despine()
        setB.savefig(f'{ERCC_folder}ERCC_Subgroup_B_plot.png')
        out_result(f'{ERCC_folder}ERCC_Subgroup_B_plot.png', 'ERCC Mix1-Mix2 Common (Group B) Plot')
        if run_main:
            plt.close()

    output(f"ERCC spike-in processing complete: {datetime.now():%Y-%m-%d %H:%M:%S}\n", exp.log_file)

    exp.tasks_complete.append('Spike')
    return exp


def bam2bw(in_bam, out_bw, job_log_folder, sample, project, stranded, log_file):

    if stranded:
        command_list = ['module rm python share-rpms65',
                        'source activate RNAseq',
                        f'bamCoverage -p 4 --filterRNAstrand forward -b {in_bam} --normalizeUsing CPM -bs 1 -o {out_bw}.cpm.fwd.bw',
                        f'bamCoverage -p 4 --filterRNAstrand reverse --scaleFactor -1 -b {in_bam} --normalizeUsing CPM -bs 1 -o {out_bw}.cpm.rev.bw'
                        ]
    else:
        command_list = ['module rm python share-rpms65',
                        'source activate RNAseq',
                        f'bamCoverage -p 4 -b {in_bam} --normalizeUsing CPM -bs 1 -o {out_bw}.cpm.bw'
                        ]

    return (send_job(command_list=command_list,
                     job_name=f'{sample}_bw',
                     job_log_folder=job_log_folder,
                     q='general',
                     mem='10000',
                     cores=2,
                     log_file=log_file,
                     project=project
                     ))


def star(exp):
    '''
    Alignment to genome using STAR
    '''
    output(f'\n Beginning genome alignments with STAR: {datetime.now():%Y-%m-%d %H:%M:%S}', exp.log_file)
    out_dir = f'{exp.scratch}STAR_results/'
    os.makedirs(out_dir, exist_ok=True)

    scan = 0
    while scan < 2:  # Loop twice to make sure source activate didn't fail the first time
        for number, sample in exp.samples.items():
            if f'{out_dir}{sample}_Aligned.sortedByCoord.out.bam' in glob.glob(f'{out_dir}*.bam'):
                pass
            else:
                output(f'Aligning using STAR to genome for {sample}.\n', exp.log_file)
                fname = f'{exp.fastq_folder}{sample}_trim'

                if exp.seq_type == 'paired':
                    align = f'STAR --runThreadN 4 --outSAMtype BAM SortedByCoordinate --genomeDir {exp.genome_indicies["STAR"]} --readFilesIn {fname}_R1.fastq.gz {fname}_R2.fastq.gz --readFilesCommand zcat --outFileNamePrefix {out_dir}{sample}_ --quantMode GeneCounts'
                elif exp.seq_type == 'single':
                    align = f'STAR --runThreadN 4 --outSAMtype BAM SortedByCoordinate --genomeDir {exp.genome_indicies["STAR"]} --readFilesIn {fname}.fastq.gz --readFilesCommand zcat --outFileNamePrefix {out_dir}{sample}_ --quantMode GeneCounts'

                command_list = ['module rm python share-rpms65',
                                'source activate RNAseq',
                                align,
                                f'samtools index {sample}_Aligned.sortedByCoord.out.bam'
                                ]

                exp.job_id.append(send_job(command_list=command_list,
                                           job_name=f'{sample}_STAR',
                                           job_log_folder=exp.job_folder,
                                           q='bigmem',
                                           mem=50000,
                                           log_file=exp.log_file,
                                           project=exp.project,
                                           cores=2))

                time.sleep(5)

        # Wait for jobs to finish
        job_wait(exp.job_id, exp.log_file)

        scan += 1

    output('STAR alignment to genome finished.', exp.log_file)

    # Generate signal files

    for sample in exp.samples.values():

        output(f'Generating bigwig singal file for {sample}.\n', exp.log_file)

        exp.job_id.append(bam2bw(in_bam=f'{out_dir}{sample}_Aligned.sortedByCoord.out.bam',
                                 out_bw=f'{out_dir}{sample}.star.',
                                 job_log_folder=exp.job_folder,
                                 sample=sample,
                                 project=exp.project,
                                 log_file=exp.log_file,
                                 stranded=exp.stranded
                                 ))

    # Generate one matrix for all counts
    try:
        counts_glob = glob.glob(f'{out_dir}*_ReadsPerGene.out.tab')
        if len(counts_glob) != exp.sample_number:
            output('At least one STAR alignment failed.', exp.log_file)
            raise RuntimeError('At least one STAR alignment failed. Check scripts and resubmit.')
        else:
            exp.count_matrix = pd.DataFrame(index=read_pd(counts_glob[1]).index)

            for number, sample in exp.samples.items():
                exp.count_matrix[sample] = read_pd(f'{out_dir}{sample}_ReadsPerGene.out.tab')[[3]]
            exp.count_matrix = exp.count_matrix.iloc[4:, :]
            exp.count_matrix.to_csv(f'{out_dir}ALL_STAR.count.matrix.txt', header=True, index=True, sep="\t")
            if os.path.isfile(exp.genome_indicies['Gene_names']):
                with open(exp.genome_indicies['Gene_names'], 'rb') as file:
                    gene_dict = pickle.load(file)
                exp.count_matrix['name'] = exp.count_matrix.index
                exp.count_matrix = exp.count_matrix[exp.count_matrix.name.isin(gene_dict.keys())]
                exp.count_matrix['name'] = exp.count_matrix.name.apply(lambda x: f'{x}_{gene_dict[x]}')
                exp.count_matrix.index = exp.count_matrix.name
                exp.count_matrix = exp.count_matrix.drop(columns=['name'])
                exp.count_matrix.to_csv(f'{out_dir}Filtered_STAR.count.matrix.txt', header=True, index=True, sep="\t")

    except:
        output('Error generating count matrix.', exp.log_file)
        raise RuntimeError('Error generating STARcount matrix. Make sure the file is not empty.')

    exp.tasks_complete.append('STAR')
    output(f'STAR alignemnt and count generation complete: {datetime.now():%Y-%m-%d %H:%M:%S}\n', exp.log_file)

    return exp


def rsem(exp):
    '''
    Alignment to transcriptome using STAR and estimating expected counts using EM with RSEM

    '''
    output(f'\n Beginning RSEM-STAR transcriptome alignments: {datetime.now():%Y-%m-%d %H:%M:%S}', exp.log_file)

    out_dir = f'{exp.scratch}RSEM_results/'
    os.makedirs(out_dir, exist_ok=True)

    cur_dir = val_folder(os.getcwd())
    os.chdir(out_dir)

    scan = 0
    while scan < 2:  # Loop twice to make sure source activate didn't fail the first time
        for number, sample in exp.samples.items():
            if f'{out_dir}{sample}.genome.sorted.bam' in glob.glob(f'{out_dir}*.bam'):
                pass
            else:
                output(f'Aligning using STAR and counting transcripts using RSEM for {sample}.\n', exp.log_file)

                if exp.seq_type == 'paired':
                    align = f'rsem-calculate-expression --star --star-gzipped-read-file --paired-end --append-names --output-genome-bam --sort-bam-by-coordinate -p 4 {exp.fastq_folder}{sample}_trim_R1.fastq.gz {exp.fastq_folder}{sample}_trim_R2.fastq.gz {exp.genome_indicies["RSEM_STAR"]} {sample}'
                elif exp.seq_type == 'single':
                    align = f'rsem-calculate-expression --star --star-gzipped-read-file --append-names --output-genome-bam --sort-bam-by-coordinate -p 4 {exp.fastq_folder}{sample}_trim.fastq.gz {exp.genome_indicies["RSEM_STAR"]} {sample}'

                plot_model = f'rsem-plot-model {sample} {sample}.models.pdf'

                command_list = ['module rm python share-rpms65',
                                'source activate RNAseq',
                                align,
                                plot_model
                                ]

                exp.job_id.append(send_job(command_list=command_list,
                                           job_name=f'{sample}_RSEM',
                                           job_log_folder=exp.job_folder,
                                           q='bigmem',
                                           mem=50000,
                                           log_file=exp.log_file,
                                           project=exp.project,
                                           cores=2
                                           ))

                time.sleep(5)

        # Wait for jobs to finish
        job_wait(exp.job_id, exp.log_file)

        scan += 1

    # Generate signal files

    for sample in exp.samples.values():

        output(f'Generating bigwig singal file for {sample}.\n', exp.log_file)

        exp.job_id.append(bam2bw(in_bam=f'{out_dir}{sample}.genome.sorted.bam',
                                 out_bw=f'{out_dir}{sample}.rsem.bw',
                                 job_log_folder=exp.job_folder,
                                 sample=sample,
                                 project=exp.project,
                                 log_file=exp.log_file,
                                 stranded=exp.stranded))

    remove_files = ['genome.bam', 'transcript.bam', 'transcript.sorted.bam', 'transcrpt.sorted.bam.bai', 'wig']
    for number, sample in exp.samples.items():
        for file in remove_files:
            del_file = f'{out_dir}{sample}.{file}'
            if os.path.isfile(del_file):
                os.remove(del_file)
            pdf = f'{out_dir}{sample}.models.pdf'
            if os.path.isdir(exp.qc_folder) and os.path.isfile(pdf):
                move(pdf, f'{exp.qc_folder}{sample}.models.pdf')

    # Generate one matrix for all expected_counts
    output(f'Generating Sample Matrix from RSEM.gene.results: {datetime.now():%Y-%m-%d %H:%M:%S}\n', exp.log_file)
    matrix = 'rsem-generate-data-matrix '
    columns = []
    for number, sample in exp.samples.items():
        matrix = f'{matrix}{out_dir}{sample}.genes.results '
        columns.append(sample)

    matrix = matrix + f'> {out_dir}RSEM.count.matrix.txt'

    command_list = ['module rm python',
                    'source activate RNAseq',
                    matrix
                    ]

    exp.job_id.append(send_job(command_list=command_list,
                               job_name='Generate_Count_Matrix',
                               job_log_folder=exp.job_folder,
                               q='general',
                               mem=1000,
                               log_file=exp.log_file,
                               project=exp.project
                               ))

    # Wait for jobs to finish
    job_wait(exp.job_id, exp.log_file)

    counts = read_pd(f'{out_dir}RSEM.count.matrix.txt')
    counts.columns = columns
    counts.to_csv(f'{out_dir}RSEM.count.matrix.txt', header=True, index=True, sep="\t")

    exp.count_matrix = counts
    exp.tasks_complete.append('RSEM')
    output(f'STAR alignemnt and RSEM counts complete: {datetime.now():%Y-%m-%d %H:%M:%S}\n', exp.log_file)

    os.chdir(cur_dir)

    return exp


def kallisto(exp):
    '''
    Second/alternate alignment to transcriptome using kallisto
    '''
    # make Kallisto_results folder

    os.makedirs(f'{exp.scratch}Kallisto_results/', exist_ok=True)

    scan = 0
    while scan < 2:

        output(f'Beginning Kallisto alignments: {datetime.now():%Y-%m-%d %H:%M:%S}\n', exp.log_file)

        # Submit kallisto for each sample
        for number, sample in exp.samples.items():

            kal_out = f'{exp.scratch}Kallisto_results/{sample}/'
            os.makedirs(kal_out, exist_ok=True)

            if f'{kal_out}abundance.tsv' not in glob.glob(f'{kal_out}*.tsv'):
                if exp.seq_type == 'paired':
                    align = f'kallisto quant --index={exp.genome_indicies["Kallisto"]} --output-dir={kal_out} --threads=4 --bootstrap-samples=100 {exp.fastq_folder}{sample}_trim_R1.fastq.gz {exp.fastq_folder}{sample}_trim_R2.fastq.gz'

                output(f'Aligning {sample} using Kallisto.\n', exp.log_file)

                command_list = ['module rm python',
                                'module rm perl',
                                'source activate RNAseq',
                                align
                                ]

                exp.job_id.append(send_job(command_list=command_list,
                                           job_name=f'{sample}_Kallisto',
                                           job_log_folder=exp.job_folder,
                                           q='general',
                                           mem=10000,
                                           log_file=exp.log_file,
                                           project=exp.project,
                                           cores=2))

        # Wait for jobs to finish
        job_wait(exp.job_id, exp.log_file)

        scan += 1

    exp.tasks_complete.append('Kallisto')

    return exp


def plot_PCA(counts, colData, out_dir, name, test_condition):
    '''
    Inputs
    ------
    counts: dataframe of counts
    colData: dataframe of colData (DESeq2 format). If '', then skip test_condition.
    out_dir: string of output directory
    name: name of PCA plot
    test_condition: name of colData column corresponding to test condition

    Outputs
    -------
    None

    prints PCA plot to out_dir

    '''

    try:
        to_remove = ['gene_name', 'id', 'name']
        for x in to_remove:
            if x in list(counts.columns):
                counts = counts.drop(x, axis=1)

        pca = PCA(n_components=2)
        bpca = pca.fit_transform(counts.T)
        pca_score = pca.explained_variance_ratio_
        bpca_df = pd.DataFrame(bpca)
        bpca_df.index = counts.T.index
        bpca_df['name'] = bpca_df.index

        plt.clf()
        fig = plt.figure(figsize=(8, 8), dpi=300)
        ax = fig.add_subplot(111)
        if len(colData) == 0:
            ax.scatter(bpca_df[0], bpca_df[1], marker='o', color='black')
        else:
            bpca_df['group'] = colData[test_condition].tolist()
            ax.scatter(bpca_df[bpca_df.group == 'yes'][0], bpca_df[bpca_df.group == 'yes'][1], marker='o', color='blue')
            ax.scatter(bpca_df[bpca_df.group == 'no'][0], bpca_df[bpca_df.group == 'no'][1], marker='o', color='red')
            red_patch = mpatches.Patch(color='red', alpha=.4, label=f'Not {test_condition}')
            blue_patch = mpatches.Patch(color='blue', alpha=.4, label=test_condition)

        ax.set_xlabel(f'PCA Component 1: {pca_score[0]:0.1%} variance')
        ax.set_ylabel(f'PCA Component 2: {pca_score[1]:0.1%} varinace')

        for i, sample in enumerate(bpca_df['name'].tolist()):
            xy = (bpca_df.iloc[i, 0], bpca_df.iloc[i, 1])
            xytext = tuple([sum(x) for x in zip(xy, ((sum(abs(ax.xaxis.get_data_interval())) * .01), (sum(abs(ax.yaxis.get_data_interval())) * .01)))])
            ax.annotate(sample, xy=xy, xytext=xytext)

        if len(colData) != 0:
            ax.legend(handles=[blue_patch, red_patch], loc=1)

        sns.despine()
        plt.tight_layout()
        plt.subplots_adjust(right=0.8, top=.8)

        os.makedirs(val_folder(out_dir), exist_ok=True)
        ax.figure.savefig(f'{out_dir}{name}_PCA.png')
        ax.figure.savefig(f'{out_dir}{name}_PCA.svg')
        plt.close()
        out_result(f'{out_dir}{name}_PCA.png', f'PCA Plot: {name}')

    except:
        raise RuntimeError('Error during plot_PCA. Fix problem then resubmit with same command to continue from last completed step.')


def GC_normalization(exp):
    '''
    Within lane loess GC normalization using EDAseq
    '''
    pandas2ri.activate()

    ri.set_writeconsole_regular(rout_write)
    ri.set_writeconsole_warnerror(rout_write)

    edaseq = importr('EDASeq')
    normCounts = ro.r('normCounts')
    as_df = ro.r("as.data.frame")

    output(f'Beginning within-lane GC length/content loess normalization for all samples: {datetime.now():%Y-%m-%d %H:%M:%S}\n', exp.log_file)

    GC_content = read_pd(exp.genome_indicies['GC_Content'])
    raw_counts = exp.count_matrix
    raw_counts['id'] = raw_counts.index
    raw_counts['id'] = raw_counts.id.apply(lambda x: x.split("_")[0].split(".")[0])
    GC_genes = set(GC_content.split.tolist())

    # Keep only counts with GC data (based on latest ensembl biomart).  see EDAseq package and use 'biomart' after dropping ensembl name version.
    GC_counts = round(raw_counts[raw_counts.id.apply(lambda x: x in GC_genes)])
    GC_gene_set = set(GC_counts.id.tolist())
    GC_content = GC_content[GC_content.split.apply(lambda x: x in GC_gene_set)]
    EDA_set = edaseq.newSeqExpressionSet(counts=GC_counts.drop(columns='id').values, featureData=GC_content)
    gcNorm = edaseq.withinLaneNormalization(EDA_set, 'gc', 'loess')
    data_norm = ro.pandas2ri.ri2py(as_df(normCounts(gcNorm)))
    data_norm.index = GC_counts.index
    data_norm.columns = GC_counts.drop(columns='id').columns
    exp.gc_count_matrix = data_norm

    output(f'Finished GC normalization: {datetime.now():%Y-%m-%d %H:%M:%S}\n', exp.log_file)
    exp.tasks_complete.append('GC')

    return exp


def RUV(RUV_data, test_type, design, reduced, colData, norm_type, log, ERCC_counts, comparison, test_condition, plot_dir, de, lfcshrink):

    '''
    Performs normalization method for removing unwanted variance across samples using controls sequences.

    Inputs
    ------
    RUV_data = pandas dataframe of counts
    type_type = 'lrt' or 'wald'
    design = string of design (ie '~Condition_A')
    reduced = string of reduced design (for use if using lrt)
    colData = pandas dataframe of DESeq2 format colData or empty list if compensating with all.
    norm_type = string 'ercc' or 'empirical'
    log = log file for output printing
    ERCC_counts = unsused if 'empirical', else a dataframe of ERCC_counts
    comparison = string of comparison name
    out_dir = directory for pca plots
    de = whether or not to perform differential expression
    test_condition = name of colData column associated with the test_condition
    '''
    try:
        pandas2ri.activate()

        ri.set_writeconsole_regular(rout_write)
        ri.set_writeconsole_warnerror(rout_write)

        deseq = importr('DESeq2')
        ruvseq = importr('RUVSeq')
        edaseq = importr('EDASeq')
        as_df = ro.r("as.data.frame")
        as_cv = ro.r('as.character')
        counts = ro.r("counts")
        normCounts = ro.r('normCounts')
        pdata = ro.r('pData')
        results = ro.r('results')

        os.makedirs(plot_dir, exist_ok=True)

        if comparison == 'ALL':
            plot_PCA(counts=RUV_data, colData=[], out_dir=plot_dir, test_condition=[], name=f'{comparison}_pre-{norm_type}_RUV_ALL_counts')
        else:
            plot_PCA(counts=RUV_data, colData=colData, out_dir=plot_dir, test_condition=test_condition, name=f'{comparison}_pre-{norm_type}_RUV_counts')

        # retain gene name
        RUV_data['name'] = RUV_data.index

        # RUVseq
        if norm_type.lower() == 'empirical':
            output(f'Performing Normalization by removing unwatned variance of empirical negative control genes for {comparison}: {datetime.now():%Y-%m-%d %H:%M:%S}\n', log)

            # determining non differentially expressed genes to use as empirical negative controls
            dds_emp = deseq.DESeqDataSetFromMatrix(countData=RUV_data.drop(columns='name').values,
                                                   colData=colData,
                                                   design=ro.Formula(design))
            dds_emp = deseq.DESeq(dds_emp)
            results_emp = pandas2ri.ri2py(as_df(deseq.results(dds_emp)))
            results_emp.index = RUV_data.index
            results_emp.sort_values(by='padj', inplace=True)
            top_de = list(results_emp.head(10000).index)

            # rename indices to reflect rpy2 conversion to R dataframe
            RUV_data.index = range(1, (len(RUV_data) + 1))

            # empirical negative controls
            empirical = list(RUV_data[RUV_data.name.apply(lambda x: x not in top_de)].drop(columns='name').index)

            # generate normalization scaling based on unwanted variance from empirical negative controls
            data_set = edaseq.newSeqExpressionSet(RUV_data.drop(columns='name').values, phenoData=colData)
            RUVg_set = ruvseq.RUVg(x=data_set, cIdx=as_cv(empirical), k=1)

            output(f'\nEmpirical negative control normalization complete for {comparison}: {datetime.now():%Y-%m-%d %H:%M:%S}\n', log)

        elif norm_type.lower() == 'ercc':
            output(f'Performing Normalization by removing unwanted variance using ERCC spike-ins for {comparison}: {datetime.now():%Y-%m-%d %H:%M:%S}\n', log)

            # rename ERCC join ERCC counts to gene counts and reindex for rpy2 R dataframe
            ERCC_counts['name'] = ERCC_counts.index
            ERCC_counts['name'] = ERCC_counts.name.apply(lambda x: f'{x}_{x}')

            RUV_data = RUV_data.append(ERCC_counts).dropna()  # added dropna
            RUV_data.index = range(1, (len(RUV_data) + 1))

            # generate index locations of ERCC spikes
            spike_list = list(RUV_data[RUV_data.name.apply(lambda x: x in list(ERCC_counts.name))].index)

            # normalize samples based on unwanted variance between ERCC spike in controls
            data_set = edaseq.newSeqExpressionSet(RUV_data.drop(columns='name').values, phenoData=colData)
            RUVg_set = ruvseq.RUVg(x=data_set, cIdx=as_cv(spike_list), k=1)
            output(f'\nERCC normalization complete for {comparison}: {datetime.now():%Y-%m-%d %H:%M:%S}\n', log)

        else:
            raise ValueError('RUV() takes only "ercc" or "empirical" as options.')

        # Differential expression (wald or lrt) to account for scaled variances between samples
        if comparison == 'ALL':
            RUV_dds = deseq.DESeqDataSetFromMatrix(countData=counts(RUVg_set), colData=pdata(RUVg_set), design=ro.Formula('~W_1'))
            RUV_dds = deseq.estimateSizeFactors_DESeqDataSet(RUV_dds)
        else:
            RUV_dds = deseq.DESeqDataSetFromMatrix(countData=counts(RUVg_set), colData=pdata(RUVg_set), design=ro.Formula(f'~W_1 + {design.split("~")[-1]}'))
            RUV_dds = deseq.DESeq(RUV_dds)
            if test_type == 'lrt':
                RUV_dds = deseq.DESeq(RUV_dds, test='LRT', reduced=ro.Formula('~W1' if reduced == '~1' else f'~W1 + {reduced.split("~")[-1]}'))

        # generate normalized counts for pca
        counts_df = pandas2ri.ri2py(as_df(normCounts(RUVg_set)))
        counts_df.columns = RUV_data.drop(columns='name').columns

        if comparison == 'ALL':
            plot_PCA(counts=counts_df.dropna(), colData=[], out_dir=plot_dir, test_condition=[], name=f'{comparison}_post-{norm_type}_RUV_raw_counts')
        else:
            plot_PCA(counts=counts_df.dropna(), colData=colData, out_dir=plot_dir, test_condition=test_condition, name=f'{comparison}_post-{norm_type}_RUV_raw_counts')

        if de:
            # extract results and relabel samples and genes
            output(f'{comparison} results type: ', log)

            results = pandas2ri.ri2py(as_df(deseq.results(RUV_dds, contrast=as_cv([f'{design.split(" ")[-1].split("~")[-1]}', 'yes', 'no']))))
            results.index = RUV_data.name

            if lfcshrink:
                output(f'Perfomring log fold change shrinking for {comparison}.  Switched to ashr method for lfc shrinkage with RUV normalizaiton.', log)
                lfc = pandas2ri.ri2py(as_df(deseq.lfcShrink(RUV_dds, contrast=as_cv([f'{design.split(" ")[-1].split("~")[-1]}', 'yes', 'no']), type='ashr')))
                lfc.index = RUV_data.name
            else:
                lfc = None

        # Change normlaized counts
        RUV_normcounts = pandas2ri.ri2py(as_df(counts(RUV_dds, normalized=True)))
        RUV_normcounts.columns = RUV_data.drop(columns='name').columns
        RUV_normcounts.index = RUV_data.name
        if comparison == 'ALL':
            plot_PCA(counts=RUV_normcounts, colData=[], out_dir=plot_dir, test_condition=[], name=f'{comparison}_post-{norm_type}_RUV_normalized_counts')
        else:
            plot_PCA(counts=RUV_normcounts, colData=colData, out_dir=plot_dir, test_condition=test_condition, name=f'{comparison}_post-{norm_type}_RUV_normalized_counts')

        output(f'Unwanted variance normalization complete for {comparison} using RUVSeq: {datetime.now():%Y-%m-%d %H:%M:%S}', log)

        if de:
            return results, lfc, RUV_normcounts
        else:
            return RUV_normcounts

    except:
        raise RuntimeError('Error during RUVseq.')


def quartile_norm(norm_type, count_matrix, colData, design, reduced, log, comparison, plot_dir, lfcshrink, test_type, test_condition):
    '''
    type = 'median','upper','full'

    '''
    output(f'Beginning differential expresion with DESeq2 and {norm_type} normlization for {comparison}', log)

    pandas2ri.activate()

    ri.set_writeconsole_regular(rout_write)
    ri.set_writeconsole_warnerror(rout_write)

    edaseq = importr('EDASeq')
    deseq = importr('DESeq2')
    normCounts = ro.r('normCounts')
    as_df = ro.r("as.data.frame")
    as_cv = ro.r('as.character')

    count_matrix['name'] = count_matrix.index
    colData_size = colData
    colData_size['sizeFactor'] = 1.0

    plot_PCA(counts=count_matrix.drop(['name']), colData=colData, out_dir=plot_dir, test_condition=test_condition, name=f'{comparison}_pre-{norm_type}_counts')

    EDAset = edaseq.newSeqExpressionSet(counts=count_matrix.values, phenoData=colData)
    EDAset = edaseq.betweenLaneNormalization(EDAset, which=f'{norm_type}')
    norm_counts = pandas2ri.ri2py(as_df(normCounts(EDAset)))

    norm_dds = deseq.DESeqDataSetFromMatrix(countData=norm_counts.values, colData=colData_size, design=ro.Fomula(design))
    norm_dds = deseq.estimateDispersions_DESeqDataSet(norm_dds)

    if test_type == 'lrt':
        norm_dds = deseq.nbinomLRT(norm_dds, full=ro.Formula(design), reduced=ro.Formula(reduced))
    else:
        norm_dds = deseq.nbinomWaldTest(norm_dds)
        results = pandas2ri.ri2py(as_df(deseq.results(norm_dds), contrast=as_cv([f'{design.split(" ")[-1].split("~")[-1]}', 'yes', 'no'])))
        results.index = count_matrix.name

    normcounts = pandas2ri.ri2py(as_df(norm_counts))
    normcounts.columns = count_matrix.drop(columns='name').columns
    normcounts.index = count_matrix.name
    log2_normcounts = np.log2(normcounts + 1)

    if lfcshrink:
        output(f'Perfomring log fold change shrinking for {comparison}.', log)
        lfc = pandas2ri.ri2py(as_df(deseq.lfcShrink(norm_dds, contrast=as_cv([f'{design.split(" ")[-1].split("~")[-1]}', 'yes', 'no']), type='apeglm')))
        lfc.index = count_matrix.name
    else:
        lfc = None

    plot_PCA(counts=log2_normcounts, colData=colData, out_dir=plot_dir, test_condition=test_condition, name=f'{comparison}_post-{norm_type}_counts')

    return results, lfc, log2_normcounts


def DESeq2(exp):
    '''
    Differential Expression using DESeq2

    Inputs
    ------
    exp.job_folder: '/path/to/job/log/folder'
    exp.log_file: 'log_file.txt'
    exp.scratch: out_dir folder
    exp.alignment_mode: 'gene' or 'transcript'
    exp.designs['colData']: pd.DataFrame()
    exp.designs['all_samples']: ['list','of','samples']
    exp.designs['design']: ex. '~ConditionA'
    exp.norm: 'ERCC','Median-Ratios','ERCC_Mixed','Empirical'
    exp.tasks_complete: []
    exp.count_matrix: pd.DataFrame()

    Optional
    --------
    exp.gc_norm: bool
    exp.gc_count_matrix
    exp.spike_counts
    exp.genome_indicies['ERCC_Mix']

    Outputs
    -------
    exp.de_results: {}
    saves file to exp.scratch/DESeq2_Results/...

    '''

    output(f'Beginning DESeq2 differential expression analysis: {datetime.now():%Y-%m-%d %H:%M:%S}\n', exp.log_file)

    pandas2ri.activate()

    ri.set_writeconsole_regular(rout_write)
    ri.set_writeconsole_warnerror(rout_write)

    deseq = importr('DESeq2')
    as_df = ro.r("as.data.frame")
    as_cv = ro.r('as.character')
    assay = ro.r("assay")
    session = ro.r("sessionInfo")

    out_dir = f'{exp.scratch}DESeq2_results/'
    os.makedirs(out_dir, exist_ok=True)

    if exp.gc_norm:
        output('Using GC normalized counts for differential expression.\n', exp.log_file)
        count_matrix = exp.gc_count_matrix[list(exp.samples.values())]
    else:
        output('Using STAR or STAR-RSEM aligned counts for differential expression.\n', exp.log_file)
        count_matrix = round(exp.count_matrix[list(exp.samples.values())])

    dds = {}

    for comparison, designs in exp.designs.items():
        output(f'Beginning {comparison}: {datetime.now():%Y-%m-%d %H:%M:%S}\n', exp.log_file)
        colData = designs['colData']
        design = ro.Formula(designs['design'])
        data = count_matrix[designs['all_samples']]

        # filtering for genes with more than 1 count in two samples
        data = round(data[data[data > 1].apply(lambda x: len(x.dropna()) > 1, axis=1)])

        if exp.norm.lower() in ['ercc', 'ercc_mixed']:
            output('Determining ERCC scaling vs Sample scaling using median of ratios of counts for rough comparison.  This may point out potentially problematic samples.\n', exp.log_file)

            dds[comparison] = deseq.DESeqDataSetFromMatrix(countData=data.values,
                                                           colData=colData,
                                                           design=design
                                                           )

            ERCC_data = round(exp.spike_counts[designs['all_samples']])
            ERCC_dds = deseq.DESeqDataSetFromMatrix(countData=ERCC_data.values, colData=colData, design=design)
            ERCC_size = deseq.estimateSizeFactors_DESeqDataSet(ERCC_dds)
            deseq2_size = deseq.estimateSizeFactors_DESeqDataSet(dds[comparison])
            sizeFactors = ro.r("sizeFactors")

            # compare size factors from DESeq2 and ERCC for inconsistencies
            ERCC_vector = pandas2ri.ri2py_vector(sizeFactors(ERCC_size))
            deseq2_vector = pandas2ri.ri2py_vector(sizeFactors(deseq2_size))
            if len(ERCC_vector) == len(deseq2_vector):
                for x in range(len(ERCC_vector)):
                    if abs((ERCC_vector[x] - deseq2_vector[x]) / (ERCC_vector[x] + deseq2_vector[x])) > 0.1:
                        output(f'ERCC spike ({x + 1} in list) is greater than 10 percent different than deseq2 size factor for {comparison}. \n', exp.log_file)
                output(f'Samples: {designs["all_samples"]}\n', exp.log_file)
                output(f'ERCC size factors: {ERCC_vector}', exp.log_file)
                output(f'DESeq2 size factors: {deseq2_vector}\n', exp.log_file)
            else:
                output(f'\nERCC and deseq2 column lengths are different for {comparison}', exp.log_file)

        # Differential Expression
        if exp.norm.lower() == 'median-ratios':
            output('Using DESeq2 standard normalization of scaling by median of the ratios of observed counts.', exp.log_file)
            output(f'Performing {designs["Test_type"]} test for differential expression for {comparison}\n', exp.log_file)

            dds[comparison] = deseq.DESeqDataSetFromMatrix(countData=data.values,
                                                           colData=colData,
                                                           design=design
                                                           )

            dds[comparison] = deseq.DESeq(dds[comparison])

            if designs['Test_type'] == 'lrt':
                reduced = ro.Formula(designs['reduced'])
                dds[comparison] = deseq.DESeq(dds[comparison], test='LRT', reduced=reduced)

            output(f'{comparison} results type: ', exp.log_file)

            # get results
            exp.de_results[f'DE2_{comparison}'] = pandas2ri.ri2py(as_df(deseq.results(dds[comparison], contrast=as_cv([f"{designs['design'].split(' ')[-1].split('~')[-1]}", 'yes', 'no']))))
            exp.de_results[f'DE2_{comparison}'].index = data.index

            # get shrunken lfc (apeglm) method)
            if exp.lfcshrink:
                output(f'Perfomring log fold change shrinkage for {comparison} using the "apeglm" method.', exp.log_file)
                coef = as_cv(f'{designs["design"].split(" ")[-1].split("~")[-1]}_yes_vs_no')
                exp.de_results[f'shrunkenLFC_{comparison}'] = pandas2ri.ri2py(as_df(deseq.lfcShrink(dds[comparison], coef=coef, type='apeglm')))
                exp.de_results[f'shrunkenLFC_{comparison}'].index = data.index

            # regularized log transformed
            exp.de_results[f'{comparison}_rlog_counts'] = pandas2ri.ri2py_dataframe(assay(deseq.rlog(dds[comparison], blind=False)))
            exp.de_results[f'{comparison}_rlog_counts'].columns = data.columns
            exp.de_results[f'{comparison}_rlog_counts'].index = data.index

        elif exp.norm.lower() == 'ercc':
            exp.de_results[f'DE2_{comparison}'], \
                exp.de_results[f'shrunkenLFC_{comparison}'], \
                exp.de_results[f'{comparison}_log2_normCounts'] = RUV(RUV_data=data,
                                                                      test_type=designs['Test_type'],
                                                                      design=designs['design'],
                                                                      reduced=designs['reduced'],
                                                                      test_condition=f'{designs["Test_condition"][-1]}',
                                                                      colData=colData,
                                                                      norm_type='ERCC',
                                                                      ERCC_counts=exp.spike_counts[designs['all_samples']],
                                                                      log=exp.log_file,
                                                                      comparison=comparison,
                                                                      plot_dir=f'{exp.scratch}PCA/{comparison}/',
                                                                      de=True,
                                                                      lfcshrink=exp.lfsshink
                                                                      )

        elif exp.norm.lower() == 'ercc_mixed':
            full_counts = exp.spike_counts[designs['all_samples']]
            mix = pd.read_csv(exp.genome_indicies['ERCC_Mix'], header=0, index_col=1, sep="\t")
            subgroupB = mix[mix.subgroup == 'B'].index.tolist()
            exp.de_results[f'DE2_{comparison}'], \
                exp.de_results[f'shrunkenLFC_{comparison}'], \
                exp.de_results[f'{comparison}_log2_normCounts'] = RUV(RUV_data=data,
                                                                      test_type=designs['Test_type'],
                                                                      design=designs['design'],
                                                                      reduced=designs['reduced'],
                                                                      test_condition=f'{designs["Test_condition"][-1]}',
                                                                      colData=colData,
                                                                      norm_type='ERCC',
                                                                      ERCC_counts=full_counts.loc[subgroupB],
                                                                      log=exp.log_file,
                                                                      comparison=comparison,
                                                                      plot_dir=f'{exp.scratch}PCA/{comparison}/',
                                                                      de=True,
                                                                      lfcshrink=exp.lfcshrink
                                                                      )

        elif exp.norm.lower() == 'empirical':
            exp.de_results[f'DE2_{comparison}'], \
                exp.de_results[f'shrunkenLFC_{comparison}'], \
                exp.de_results[f'{comparison}_log2_normCounts'] = RUV(RUV_data=data,
                                                                      test_type=designs['Test_type'],
                                                                      design=designs['design'],
                                                                      reduced=designs['reduced'],
                                                                      test_condition=f'{designs["Test_condition"][-1]}',
                                                                      colData=colData,
                                                                      norm_type='empirical',
                                                                      ERCC_counts=None,
                                                                      log=exp.log_file,
                                                                      comparison=comparison,
                                                                      plot_dir=f'{exp.scratch}PCA/{comparison}/',
                                                                      de=True,
                                                                      lfcshrink=exp.lfcshrink
                                                                      )

        elif exp.norm.lower() in ['upper', 'median', 'full']:
            exp.de_results[f'DE2_{comparison}'], \
                exp.de_results[f'shrunkenLFC_{comparison}'], \
                exp.de_results[f'{comparison}_log2_normCounts'] = quartile_norm(count_matrix=data,
                                                                                colData=colData,
                                                                                norm_type=exp.norm.lower(),
                                                                                design=designs['design'],
                                                                                reduced=designs['reduced'],
                                                                                log=exp.log_file,
                                                                                comparison=comparison,
                                                                                plot_dir=f'{exp.scratch}PCA/{comparison}/',
                                                                                lfcshrink=exp.lfchrink,
                                                                                test_type=designs['Test_type'],
                                                                                test_condition=f'{designs["Test_condition"][-1]}'
                                                                                )

        else:
            raise ValueError('Can only use "median-ratios", "ercc", "ercc_mixed", "empirical", "upper", "median", or "full" for normalization of DESeq2.')

        # DESeq2 results
        exp.de_results[f'DE2_{comparison}'].sort_values(by='padj', ascending=True, inplace=True)
        exp.de_results[f'DE2_{comparison}']['gene_name'] = exp.de_results[f'DE2_{comparison}'].index
        exp.de_results[f'DE2_{comparison}']['gene_name'] = exp.de_results[f'DE2_{comparison}'].gene_name.apply(lambda x: x.split("_")[1])
        exp.de_results[f'DE2_{comparison}'].to_csv(f'{out_dir}{comparison}-DESeq2-results.txt', header=True, index=True, sep="\t")
        # Shrunken LFC using apeglm or ashr method
        if exp.lfcshrink:
            exp.de_results[f'shrunkenLFC_{comparison}'].sort_values(by='log2FoldChange', ascending=False, inplace=True)
            exp.de_results[f'shrunkenLFC_{comparison}']['gene_name'] = exp.de_results[f'shrunkenLFC_{comparison}'].index
            exp.de_results[f'shrunkenLFC_{comparison}']['gene_name'] = exp.de_results[f'shrunkenLFC_{comparison}'].gene_name.apply(lambda x: x.split("_")[1])
            exp.de_results[f'shrunkenLFC_{comparison}'].to_csv(f'{out_dir}{comparison}-DESeq2-shrunken-LFC.txt', header=True, index=True, sep="\t")

        if exp.norm == 'median-ratios':
            count_type = f'{comparison}_rlog_counts'
        else:
            count_type = f'{comparison}_log2_normCounts'

        # Normlized counts
        exp.de_results[count_type].to_csv(f'{out_dir}{comparison}-rlog-counts.txt', header=True, index=True, sep="\t")

    # blind rlog count matrix for all samples.
    colData = pd.DataFrame(index=count_matrix.columns, data={'condition': ['A'] * exp.sample_number})
    design = ro.Formula("~1")

    # count_matrix = round(count_matrix[count_matrix[count_matrix > 5].apply(lambda x: len(x.dropna()) > 1 , axis=1)])
    dds_blind = deseq.DESeqDataSetFromMatrix(countData=count_matrix.values, colData=colData, design=design)
    exp.de_results['blind_rlog'] = pandas2ri.ri2py_dataframe(assay(deseq.rlog(dds_blind)))
    exp.de_results['blind_rlog'].index = count_matrix.index
    exp.de_results['blind_rlog'].columns = count_matrix.columns
    exp.de_results['blind_rlog']['gene_name'] = exp.de_results['blind_rlog'].index
    exp.de_results['blind_rlog']['gene_name'] = exp.de_results['blind_rlog'].gene_name.apply(lambda x: x.split("_")[1])
    exp.de_results['blind_rlog'].to_csv(f'{out_dir}ALL-samples-blind-rlog-counts.txt', header=True, index=True, sep="\t")

    '''
    dds_complete = deseq.DESeqDataSetFromMatrix(countData=count_matrix.values,
                                                colData=exp.designs['complete']['colData'],
                                                design=ro.Formula(exp.designs['complete']['design'])
                                                )

    exp.de_results['complete_rlog'] = pandas2ri.ri2py_dataframe(assay(deseq.rlog(dds_complete)))
    exp.de_results['complete_rlog'].index = count_matrix.index
    exp.de_results['complete_rlog'].columns = count_matrix.columns
    exp.de_results['complete_rlog']['gene_name'] = exp.de_results['complete_rlog'].index
    exp.de_results['complete_rlog']['gene_name'] = exp.de_results['complete_rlog'].gene_name.apply(lambda x: x.split("_")[1])
    exp.de_results['complete_rlog'].to_csv('{}ALL-samples-complete-design-rlog-counts.txt'.format(out_dir),
                                           header=True,
                                           index=True,
                                           sep="\t"
                                           )
    '''

    if exp.norm.lower() == 'ercc':
        exp.de_results['all_ERCC_normCounts'] = RUV(RUV_data=count_matrix,
                                                    test_type='',
                                                    design=[],
                                                    reduced=[],
                                                    test_condition=[],
                                                    colData=colData,
                                                    norm_type='ERCC',
                                                    ERCC_counts=round(exp.spike_counts),
                                                    log=exp.log_file,
                                                    comparison='ALL',
                                                    plot_dir=f'{exp.scratch}PCA/ALL/',
                                                    de=False,
                                                    lfcshrink=exp.lfcshrink
                                                    )

        exp.de_results['all_ERCC_normCounts'].to_csv(f'{out_dir}ALL-samples-ERCC_DE2_normCounts.txt',
                                                     header=True,
                                                     index=True,
                                                     sep="\t"
                                                     )

    output(session(), exp.log_file)
    exp.tasks_complete.append('DESeq2')
    output(f'DESeq2 differential expression complete: {datetime.now():%Y-%m-%d %H:%M:%S}\n', exp.log_file)

    return exp


def plot_exp(data, plot_dir, ylabel, name, log_file=''):
    '''
    Inputs
    ------
    data: dataframe.  samples in columns, counts in rows
    plot_dir: output directory string
    ylabel: string of data type for ylabel (ex 'Normalized Log$_2$')
    name : title for plot and for file
    log_file: ouput log file

    Ouputs
    ------
    None
    Prints boxplot to plot_dir

    '''
    output('Starting global sample expression comparisons.', log_file)

    plt.clf()
    sns.set(context='paper', font='Arial', style='white', rc={'figure.dpi': 300, 'figure.figsize': (4, 4)})
    pl = sns.boxplot(data=data, color='darkgrey', medianprops={'color': 'red'})
    pl.set_ylabel(f'Expression Counts\n({ylabel})')
    pl.set_title(name)
    for tick in pl.xaxis.get_ticklabels():
            tick.set_rotation(90)
    plt.tight_layout()
    sns.despine()
    plt.savefig(f'{plot_dir}{name}_expression_barplot.png', dpi=300)
    plt.close()
    out_result(f'{plot_dir}{name}_expression_barplot.png', f'Expression Barplot: {name}')


def Principal_Component_Analysis(exp):

    out_dir = f'{exp.scratch}PCA/'
    all_out = f'{out_dir}ALL/'
    os.makedirs(out_dir, exist_ok=True)

    # PCA on raw data
    output('Starting PCA analysis for all raw counts.', exp.log_file)
    plot_PCA(counts=exp.count_matrix[list(exp.samples.values())],
             colData=[],
             out_dir=all_out,
             test_condition=[],
             name='all_raw_counts'
             )

    plot_exp(data=exp.count_matrix[list(exp.samples.values())],
             plot_dir=all_out,
             exp_type='raw counts',
             name='all_raw_counts',
             log_file=exp.log_file
             )

    output('Starting PCA analysis for DESeq2 regularized log counts of all samples blinded.', exp.log_file)
    plot_PCA(counts=exp.de_results['blind_rlog'],
             colData=[],
             out_dir=all_out,
             test_condition=[],
             name='all_samples_blind_rlog'
             )

    plot_exp(data=exp.de_results['blind_rlog'],
             plot_dir=all_out,
             exp_type='Normalized log$_2$',
             name='all_samples_blind_rlog',
             log_file=exp.log_file
             )

    '''
    output('Starting PCA analysis for DESeq2 regularized log counts of all samples with full design.', exp.log_file)
    plot_PCA(counts=exp.de_results['complete_rlog'],
             colData=[],
             out_dir=all_out,
             test_condition=[],
             name='all_samples_complete_rlog'
             )

    plot_exp(data=exp.de_results['complete_rlog'],
             plot_dir=all_out,
             exp_type='Normalized log$_2$',
             name='all_samples_comoplete_rlog',
             log_file=exp.log_file
             )
    '''

    if exp.norm.lower() == 'ercc':
        output('starting PCA analysis for ALL ERCC-normalized counts.', exp.log_file)
        plot_PCA(counts=exp.de_results['all_ERCC_normCounts'],
                 colData=[],
                 out_dir=all_out,
                 test_condition=[],
                 name='all_ercc_normalized_normCounts'
                 )

        plot_exp(data=exp.de_results['all_ERCC_normCounts'],
                 plot_dir=all_out,
                 exp_type='Normalized',
                 name='all_ercc_normalized_normCounts',
                 log_file=exp.log_file
                 )

    for comparison, design in exp.designs.items():
        if exp.norm == 'median-ratios':
            count_type = f'{comparison}_rlog_counts'
        else:
            count_type = f'{comparison}_log2_normCounts'

        output(f'Starting DESeq2 PCA analysis for {comparison}: {datetime.now():%Y-%m-%d %H:%M:%S}\n', exp.log_file)
        plot_PCA(counts=exp.de_results[count_type],
                 colData=design['colData'],
                 out_dir=f'{out_dir}{comparison}/',
                 test_condition=f'{design["Test_condition"][-1]}',
                 name=comparison
                 )

        plot_exp(data=exp.de_results[count_type],
                 plot_dir=f'{out_dir}{comparison}/',
                 exp_type='Normalized log$_2$',
                 name=comparison,
                 log_file=exp.log_file
                 )

    if exp.gc_norm:
        output('starting PCA analysis for gc normalized raw counts.', exp.log_file)
        plot_PCA(counts=exp.gc_count_matrix,
                 colData=[],
                 out_dir=all_out,
                 test_condition=[],
                 name='gc_normalized_raw_counts'
                 )

        plot_exp(data=exp.gc_count_matrix,
                 plot_dir=all_out,
                 exp_type='GC Normalized Raw',
                 name='gc_normalized_raw_counts',
                 log_file=exp.log_file
                 )

    exp.tasks_complete.append('PCA')
    output(f'PCA for DESeq2 groups complete: {datetime.now():%Y-%m-%d %H:%M:%S}\n', exp.log_file)

    return exp


def Sleuth(exp):
    '''
    NEEDS TO BE UPDATED WITH 0.v7 Strategy

    Differential expression using sleuth from the Pachter lab: https://pachterlab.github.io/sleuth/
    '''

    '''
    pandas2ri.activate()

    ri.set_writeconsole_regular(rout_write)
    ri.set_writeconsole_warnerror(rout_write)

    sleuth = importr('sleuth')
    biomart = importr('biomaRt')
    dplyr = importr('dplyr', on_conflict="warn")
    session=r("sessionInfo")
    out_dir= '{}Sleuth_results/'.format(exp.scratch)
    kal_dir= '{}Kallisto_results/'.format(exp.scratch)
    os.makedirs(out_dir, exist_ok=True)

    for comparison,design in exp.designs.items():
        output('Beginning Sleuth differential expression analysis for {}: {:%Y-%m-%d %H:%M:%S}\n'.format(comparison, datetime.now()), exp.log_file)

        path = []
        for name in design['colData'].index.tolist():
            path.append(kal_dir + name)

        if 'compensation' in design['colData'].columns.tolist():
            s2c = pd.DataFrame({'sample': design['colData'].index.tolist(),
                                'compensation': design['colData'].compensation.tolist(),
                                'condition': design['colData'].main_comparison.tolist(),
                                'path': path
                               },
                               index=range(1, len(path)+1)
                              )
            s2c = s2c[['sample','compensation','condition','path']]
            condition=Formula('~ compensation + condition')
            reduced = Formula('~compensation')
        else:
            s2c = pd.DataFrame({'sample': design['colData'].index.tolist(),
                                'condition': design['colData'].main_comparison.tolist(),
                                'path': path
                               },
                               index=range(1, len(path)+1)
                              )
            s2c = s2c[['sample','condition','path']]
            condition=Formula('~ condition')
            reduced = Formula('~1')

        globalenv["s2c"] = s2c
        r('s2c$path = as.character(s2c$path)')
        s2c = globalenv["s2c"]

        if exp.genome == 'mm10':
            mart = biomart.useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset = "mmusculus_gene_ensembl",host = "useast.ensembl.org")
        elif exp.genome == 'hg38':
            mart = biomart.useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset = "hsapiens_gene_ensembl", host = 'useast.ensembl.org')
        elif exp.genome == 'hg19':
            mart = biomart.useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset = "hsapiens_gene_ensembl", host = 'useast.ensembl.org')
        else:
            raise ValueError('Error in sleuth, pipeline only handles hg38,hg19, and mm10')

        t2g = biomart.getBM(attributes = ro.StrVector(("ensembl_transcript_id_version", "ensembl_gene_id","external_gene_name")), mart=mart)
        t2g = dplyr.rename(t2g, target_id = 'ensembl_transcript_id_version', ens_gene = 'ensembl_gene_id', ext_gene = 'external_gene_name')

        so = sleuth.sleuth_prep(s2c, target_mapping = t2g, num_cores=1, aggregation_column = 'ens_gene')
        so = sleuth.sleuth_fit(so, condition, 'full')
        so = sleuth.sleuth_fit(so, reduced, 'reduced')
        so = sleuth.sleuth_lrt(so, 'reduced', 'full')
        output(sleuth.models(so), exp.log_file)
        sleuth_table=sleuth.sleuth_results(so, 'reduced:full','lrt',show_all=True)
        exp.de_results['SL_{}'.format(comparison)] = pandas2ri.ri2py(sleuth_table)
        exp.de_results['SL_{}'.format(comparison)].to_csv('{}{}_slueth_results.txt'.format(out_dir,comparison), header=True, index=True, sep="\t")

        output(session(), exp.log_file)
        output('Sleuth differential expression complete: {:%Y-%m-%d %H:%M:%S}\n'.format(datetime.now()), exp.log_file)
    '''
    exp.tasks_complete.append('Sleuth')
    return exp


def volcano(results, sig_up, sig_down, name, out_dir):
    '''
    Generate volcano plot from deseq2 results dataframe and significant genes

    Inputs
    ------
    results: deseq2 results as a dataframe
    sig_up: set or list of genes up to be highlithed
    sig_down: set or list of genes down to be highlighted
    name: string of name of plot
    out_dir: string of output directory

    Outputs
    -------
    None

    Saves plot to file in out_dir

    '''

    plt.clf()
    sns.set(context='paper', style='white', font_scale=1)
    fig = plt.figure(figsize=(6, 6), dpi=300)
    ax = fig.add_subplot(111)

    results['logp'] = results.pvalue.apply(lambda x: -np.log10(x))

    scatter = ax.scatter(results.log2FoldChange, results.logp, marker='o', color='gray', alpha=0.1, s=10, label='_nolegend_')

    sig_results = results[results.padj < 0.05]
    scatter = ax.scatter(sig_results[sig_results.gene_name.apply(lambda x: x in sig_up)].log2FoldChange,
                         sig_results[sig_results.gene_name.apply(lambda x: x in sig_up)].logp,
                         marker='o', alpha=0.3, color='firebrick', s=10, label='Genes UP'
                         )

    scatter = ax.scatter(sig_results[sig_results.gene_name.apply(lambda x: x in sig_down)].log2FoldChange,
                         sig_results[sig_results.gene_name.apply(lambda x: x in sig_down)].logp,
                         marker='o', alpha=0.3, color='steelblue', s=10, label='Genes DOWN'
                         )

    ax.axes.set_xlabel('Fold Change (log$_2$)')
    ax.axes.set_ylabel('p-value (-log$_10$)')

    ax.legend(loc='upper left', markerscale=3)
    fig.suptitle(name)

    sns.despine()
    plt.tight_layout()
    plt.savefig(f'{out_dir}/{name}-Volcano-Plot.png', dpi=200)
    plt.savefig(f'{out_dir}/{name}-Volcano-Plot.svg', dpi=200)
    out_result(f'{out_dir}/{name}-Volcano-Plot.png', f'Volcano Plot: {name}')
    if run_main:
        plt.close()

    return


def sigs(exp):
    '''
    Identifies significantly differentially expressed genes at 2 fold and 1.5 fold cutoffs with q<0.05. Generates Volcano Plots of results.
    '''
    out_dir = f'{exp.scratch}Sigs_and_volcano_plots/'
    os.makedirs(out_dir, exist_ok=True)

    for comparison, design in exp.designs.items():

        if exp.de_sig_overlap:
            output(f'Performing overlaps of signifcant genes from Kallisto/Sleuth and STAR/RSEM/DESeq2 for {comparison}.', exp.log_file)

            exp.sig_lists[comparison] = {}
            DE_results = exp.de_results[f'DE2_{comparison}']
            SL_results = exp.de_results[f'SL_{comparison}']
            SL_sig = set(SL_results[SL_results.qval < 0.05].ext_gene.tolist())

            DE2_2UP = set(DE_results[(DE_results.padj < 0.05) & (DE_results.log2FoldChange > 1)].gene_name.tolist())
            DE2_2DN = set(DE_results[(DE_results.padj < 0.05) & (DE_results.log2FoldChange < -1)].gene_name.tolist())
            DE2_15UP = set(DE_results[(DE_results.padj < 0.05) & (DE_results.log2FoldChange > .585)].gene_name.tolist())
            DE2_15DN = set(DE_results[(DE_results.padj < 0.05) & (DE_results.log2FoldChange < -.585)].gene_name.tolist())
            DE2_UP = set(DE_results[(DE_results.padj < 0.05) & (DE_results.log2FoldChange > 0)].gene_name.tolist())
            DE2_DN = set(DE_results[(DE_results.padj < 0.05) & (DE_results.log2FoldChange < 0)].gene_name.tolist())

            exp.sig_lists[comparison]['2FC_UP'] = DE2_2UP & SL_sig
            exp.sig_lists[comparison]['2FC_DN'] = DE2_2DN & SL_sig
            exp.sig_lists[comparison]['15FC_UP'] = DE2_15UP & SL_sig
            exp.sig_lists[comparison]['15FC_DN'] = DE2_15DN & SL_sig
            exp.sig_lists[comparison]['All_UP'] = DE2_UP & SL_sig
            exp.sig_lists[comparison]['All_DN'] = DE2_DN & SL_sig

        else:
            output(f'Only using significant genes called from DESeq2 for {comparison} analyses.', exp.log_file)

            DE_results = exp.de_results[f'DE2_{comparison}']

            exp.sig_lists[comparison] = {}

            DE2_2UP = set(DE_results[(DE_results.padj < 0.05) & (DE_results.log2FoldChange > 1)].gene_name.tolist())
            DE2_2DN = set(DE_results[(DE_results.padj < 0.05) & (DE_results.log2FoldChange < -1)].gene_name.tolist())
            DE2_15UP = set(DE_results[(DE_results.padj < 0.05) & (DE_results.log2FoldChange > .585)].gene_name.tolist())
            DE2_15DN = set(DE_results[(DE_results.padj < 0.05) & (DE_results.log2FoldChange < -.585)].gene_name.tolist())
            DE2_UP = set(DE_results[(DE_results.padj < 0.05) & (DE_results.log2FoldChange > 0)].gene_name.tolist())
            DE2_DN = set(DE_results[(DE_results.padj < 0.05) & (DE_results.log2FoldChange < 0)].gene_name.tolist())

            exp.sig_lists[comparison]['2FC_UP'] = DE2_2UP
            exp.sig_lists[comparison]['2FC_DN'] = DE2_2DN
            exp.sig_lists[comparison]['15FC_UP'] = DE2_15UP
            exp.sig_lists[comparison]['15FC_DN'] = DE2_15DN
            exp.sig_lists[comparison]['All_UP'] = DE2_UP
            exp.sig_lists[comparison]['All_DN'] = DE2_DN

        # volcano_plot
        volcano_out = f'{out_dir}{comparison}/'
        os.makedirs(volcano_out, exist_ok=True)

        output('Generating Volcano Plots using DESeq2 results for significance', exp.log_file)
        volcano(results=DE_results, sig_up=DE2_2UP, sig_down=DE2_2DN, name=f'{comparison}_2_FC', out_dir=volcano_out)
        volcano(results=DE_results, sig_up=DE2_15UP, sig_down=DE2_15DN, name=f'{comparison}_1.5_FC', out_dir=volcano_out)
        volcano(results=DE_results, sig_up=DE2_UP, sig_down=DE2_DN, name=f'{comparison}_noFC_filter', out_dir=volcano_out)

    for comparison, sigs in exp.sig_lists.items():
        sig_out = f'{out_dir}{comparison}/'
        os.makedirs(sig_out, exist_ok=True)
        for sig, genes in sigs.items():
            with open(f'{sig_out}{sig}.txt', 'w') as file:
                for gene in genes:
                    file.write(f'{gene}\n')

    output(f'Signature and Volcano Plot generation complete: {datetime.now():%Y-%m-%d %H:%M:%S}\n', exp.log_file)
    exp.tasks_complete.append('Sigs')
    return exp


def clustermap(exp):
    '''
    Generate heatmap of differentially expressed genes using regularized log2counts.
    '''

    heat_dir = f'{exp.scratch}Heatmaps/'

    for comparison, design in exp.designs.items():
        out_dir = f'{heat_dir}{comparison}/'
        os.makedirs(out_dir, exist_ok=True)

        rlog = exp.de_results[f'{comparison}_rlog_counts']
        rlog['gene_name'] = [name.split("_")[1] for name in rlog.index.tolist()]

        sig = set(exp.sig_lists[comparison]['2FC_UP'] | exp.sig_lists[comparison]['2FC_DN'])
        if len(sig) < 2:
            output(f'There are not enough significantly differentially expressed genes with 2 fold chagnes in {comparison}.  Ignoring heatmap for this group. \n', exp.log_file)
        else:
            plt.clf()
            CM = sns.clustermap(rlog[rlog.gene_name.apply(lambda x: x in sig)].drop('gene_name', axis=1), z_score=0, method='complete', cmap='RdBu_r', yticklabels=False)
            CM.savefig(f'{out_dir}{comparison}_2FC_Heatmap.png', dpi=300)
            CM.savefig(f'{out_dir}{comparison}_2FC_Heatmap.svg', dpi=300)
            out_result(f'{out_dir}{comparison}_2FC_Heatmap.png', f'Heatmap (2FC): {comparison}')
            if run_main:
                plt.close()

        sig15 = set(exp.sig_lists[comparison]['15FC_UP'] | exp.sig_lists[comparison]['15FC_DN'])
        if len(sig15) < 2:
            output(f'There are not enough significantly differentially expressed genes with 1.5 fold chagnes in {comparison}.  Ignoring heatmap for this group. \n', exp.log_file)
        else:
            plt.clf()
            CM15 = sns.clustermap(rlog[rlog.gene_name.apply(lambda x: x in sig15)].drop('gene_name', axis=1), z_score=0, method='complete', cmap='RdBu_r', yticklabels=False)
            CM15.savefig(f'{out_dir}{comparison}_1.5FC_Heatmap.png', dpi=300)
            CM15.savefig(f'{out_dir}{comparison}_1.5FC_Heatmap.svg', dpi=300)
            out_result(f'{out_dir}{comparison}_1.5FC_Heatmap.png', f'Heatmap (1.5 FC): {comparison}')
            if run_main:
                plt.close()

        sigAll = set(exp.sig_lists[comparison]['All_UP'] | exp.sig_lists[comparison]['All_DN'])
        if len(sigAll) < 2:
            output(f'There are not enough significantly differentially expressed genes without a fold change in {comparison}.  Ignoring heatmap for this group. \n', exp.log_file)
        else:
            plt.clf()
            CM15 = sns.clustermap(rlog[rlog.gene_name.apply(lambda x: x in sigAll)].drop('gene_name', axis=1), z_score=0, method='complete', cmap='RdBu_r', yticklabels=False)
            CM15.savefig(f'{out_dir}{comparison}_noFCfilter_Heatmap.png', dpi=300)
            CM15.savefig(f'{out_dir}{comparison}_noFCfilter_Heatmap.svg', dpi=300)
            out_result(f'{out_dir}{comparison}_noFCfilter_Heatmap.png', f'Heatmap (no filter): {comparison}')
            if run_main:
                plt.close()

    exp.tasks_complete.append('Heatmaps')
    output('Heatmaps for DESeq2 differentially expressed genes complete: {:%Y-%m-%d %H:%M:%S}\n'.format(datetime.now()), exp.log_file)

    return exp


def enrichr(gene_list, description, out_dir):
    '''
    Perform GO enrichment and KEGG enrichment Analysis using Enrichr: http://amp.pharm.mssm.edu/Enrichr/
    '''
    gene_sets = 'KEGG_2016'
    gseapy.enrichr(gene_list=gene_list, description=description, gene_sets='KEGG_2016', outdir=out_dir)

    out_result(f'{out_dir}{gene_sets}.{description}.enrichr.reports.png', f'Enrichr: {gene_sets} for {description}')

    gene_sets = 'GO_Biological_Process_2017b'
    gseapy.enrichr(gene_list=gene_list, description=description, gene_sets=gene_sets, outdir=out_dir)

    out_result(f'{out_dir}{gene_sets}.{description}.enrichr.reports.png', f'Enrichr: {gene_sets} for {description}')

    gene_sets = 'GO_Molecular_Function_2017b'
    gseapy.enrichr(gene_list=gene_list, description=description, gene_sets='GO_Molecular_Function_2017b', outdir=out_dir)

    out_result(f'{out_dir}{gene_sets}.{description}.enrichr.reports.png', f'Enrichr: {gene_sets} for {description}')

    return


def GO_enrich(exp):
    '''
    Perform GO enrichment analysis on significanttly differentially expressed genes.
    '''
    GO_dir = f'{exp.scratch}GO_enrichment/'
    os.makedirs(GO_dir, exist_ok=True)

    for comparison, design in exp.designs.items():
        output(f'Beginning GO enrichment for {comparison}: {datetime.now():%Y-%m-%d %H:%M:%S}\n', exp.log_file)

        for name, sig in exp.sig_lists[comparison].items():
            if len(sig) == 0:
                output(f'There are no significantly differentially expressed genes in {name} {comparison}.  Ignoring gene enrichment. \n', exp.log_file)
            else:
                GO_out = f'{GO_dir}{comparison}/'
                os.makedirs(GO_out, exist_ok=True)
                enrichr(gene_list=list(sig), description=f'{comparison}_{name}', out_dir=GO_out)

    exp.tasks_complete.append('GO_enrich')
    output(f'GO Enrichment analysis for DESeq2 differentially expressed genes complete: {datetime.now():%Y-%m-%d %H:%M:%S}\n', exp.log_file)

    return exp


def gsea_barplot(out_dir, pos_file, neg_file, gmt_name, max_number=20):
    '''
    Inputs
    ------
    out_dir: directory output or '' for current directory
    pos_file: GSEA positive enrichment .xls file
    neg_file: GSEA negative enrichment .xls file
    gmt_name: name of enrichment (ex: Hallmarks)
    max_number: max number of significant sets to report (default 20)

    Returns
    -------
    string of save file
    top pos term
    top neg term

    '''

    out_dir = val_folder(out_dir)
    pos = pd.read_table(pos_file).head(max_number) if os.path.isfile(pos_file) else pd.DataFrame(columns=['FDR q-val'])
    top_pos = pos.NAME.tolist()[0]
    pos[gmt_name] = [' '.join(name.split('_')[1:]) for name in pos.NAME.tolist()]
    neg = pd.read_table(neg_file).head(max_number) if os.path.isfile(neg_file) else pd.DataFrame(columns=['FDR q-val'])
    top_neg = neg.NAME.tolist()[0]
    neg[gmt_name] = [' '.join(name.split('_')[1:]) for name in neg.NAME.tolist()]

    plt.clf()
    sns.set(context='paper', font='Arial', font_scale=.9, style='white', rc={'figure.dpi': 300, 'figure.figsize': (8, 6)})
    fig, (ax1, ax2) = plt.subplots(ncols=1, nrows=2)
    fig.suptitle(f'{gmt_name} GSEA enrichment\n(q<0.05, max {max_number})')

    if len(pos[pos['FDR q-val'] < 0.05]) > 0:
        UP = sns.barplot(data=pos[pos['FDR q-val'] < 0.05], x='NES', y=gmt_name, color='firebrick', ax=ax1)
        UP.set_title('Positive Enrichment')
        sns.despine()
    else:
        ax1.text(0.5, 0.5, 'No Significant Positive Enrichments.',
                 horizontalalignment='center',
                 verticalalignment='center',
                 transform=ax1.transAxes
                 )

    if len(neg[neg['FDR q-val'] < 0.05]) > 0:
        DN = sns.barplot(data=neg[neg['FDR q-val'] < 0.05], x='NES', y=gmt_name, color='steelblue', ax=ax2)
        DN.set_title('Negative Enrichment')
        sns.despine()
    else:
        ax2.text(0.5, 0.5, 'No Significant Netagive Enrichments.',
                 horizontalalignment='center',
                 verticalalignment='center',
                 transform=ax2.transAxes
                 )

    try:
        plt.tight_layout(h_pad=1, w_pad=1)
    except ValueError:
        pass

    plt.subplots_adjust(top=0.88)
    file = f'{out_dir}{gmt_name}_GSEA_NES_plot.png'
    fig.savefig(file, dpi=300)
    out_result(file, f'GSEA Normalized Enrichment Plot: {gmt_name}')
    if run_main:
        plt.close()

    pos_png = glob.glob(f'{out_dir}*/enplot*{top_pos}*.png')
    out_result(pos_png[0], f'Top positive {gmt_name} GSEA')

    neg_png = glob.glob(f'{out_dir}*/enplot*{top_neg}*.png')
    out_result(neg_png[0], f'Top negative {gmt_name} GSEA')

    return file


def GSEA(exp):
    '''
    Perform Gene Set Enrichment Analysis using gsea 3.0 from the Broad Institute.
    '''

    output('Starting GSEA enrichment.', exp.log_file)

    out_dir = f'{exp.scratch}DESeq2_GSEA'
    os.makedirs(out_dir, exist_ok=True)
    cwd = val_folder(os.getcwd())

    if exp.genome == 'mm10':
        gmt_list = glob.glob(f'{exp.genome_indicies["GMT"]}*.gmt')
        gmts = {'Hallmarks': [gmt for gmt in gmt_list if 'h.all' in gmt][0],
                'KEGG': [gmt for gmt in gmt_list if 'c2.cp.kegg' in gmt][0],
                'GO_Biological_Process': [gmt for gmt in gmt_list if 'c5.bp' in gmt][0],
                'GO_Molecular_Function': [gmt for gmt in gmt_list if 'c5.mf' in gmt][0],
                'Curated_Gene_Sets': [gmt for gmt in gmt_list if 'c2.cgp' in gmt][0]
                }
    else:
        gmts = {'Hallmarks': 'gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/h.all.v6.2.symbols.gmt',
                'KEGG': 'gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c2.cp.kegg.v6.2.symbols.gmt',
                'GO_Biological_Process': 'gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c5.bp.v6.2.symbols.gmt',
                'GO_Molecular_Function': 'gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c5.mf.v6.2.symbols.gmt',
                'Curated_Gene_Sets': 'gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c2.cgp.v6.2.symbols.gmt'
                }

    for comparison, design in exp.designs.items():

        indexpath = glob.glob(f'{out_dir}/{comparison}/Hallmarks/*/index.html')
        if len(indexpath) > 0:
            output(f'GSEA Hallmarks already done for {comparison}. Skipping GSEA for {comparison}\n', exp.log_file)
        else:
            output(f'Running GSEA for {comparison}. Results found in {exp.out_dir}DESeq2_GSEA/{comparison}. \n', exp.log_file)
            out_compare = f'{out_dir}/{comparison}'
            os.makedirs(out_compare, exist_ok=True)
            os.chdir(out_compare)

            results = exp.de_results[f'DE2_{comparison}'].dropna()

            # generate ranked list based on shrunken log2foldchange
            if exp.lfcshrink:
                lfc = exp.de_results[f'shrunkenLFC_{comparison}'].dropna()
                lfc.sort_values(by='log2FoldChange', ascending=False, inplace=True)
                lfc.index = lfc.gene_name
                lfc = lfc.log2FoldChange.dropna()
                lfc.to_csv(f'{out_compare}/{comparison}_shrunkenLFC.rnk', header=False, index=True, sep="\t")
                rnk2 = f'{out_compare}/{comparison}_shrunkenLFC.rnk'

            output('Using Wald statistic for gene preranking.', exp.log_file)
            rnk = f'{comparison}_stat.rnk'

            if exp.genome == 'mm10':
                results['Ens_ID'] = [ID.split('.')[0] for ID in results.index.tolist()]
                results[['Ens_ID', 'stat']].dropna(subset=['stat']).to_csv(f'{out_compare}/{comparison}_stat.rnk', header=None, index=None, sep="\t")
            else:
                results[['gene_name', 'stat']].dropna(subset=['stat']).to_csv(f'{out_compare}/{comparison}_stat.rnk', header=None, index=None, sep="\t")

            output(f'Beginning GSEA enrichment for {comparison} using preranked genes: {datetime.now():%Y-%m-%d %H:%M:%S}', exp.log_file)
            output('Genes with positive LFC (to the left left in GSEA output graph) are upregulated in experimental vs control conditions. Genes with negative LFC (on right) are downregulated genes in experimental samples vs controls.\n', exp.log_file)

            for name, gset in gmts.items():
                set_dir = f'{out_compare}/{name}'
                os.makedirs(set_dir, exist_ok=True)

                command_list = ['module rm python java perl share-rpms65',
                                'source activate RNAseq',
                                f'java -cp {exp.genome_indicies["GSEA_jar"]} -Xmx2048m xtools.gsea.GseaPreranked -gmx {gset} -norm meandiv -nperm 1000 -rnk "{rnk}" -scoring_scheme weighted -rpt_label {comparison}_{name}_wald -create_svgs false -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 1000 -set_min 10 -zip_report false -out {name} -gui false'
                                ]

                if exp.lfcshrink:
                    command_list.append(f'java -cp {exp.genome_indicies["GSEA_jar"]} -Xmx2048m xtools.gsea.GseaPreranked -gmx gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/{gset}.v6.1.symbols.gmt -norm meandiv -nperm 1000 -rnk {rnk2} -scoring_scheme weighted -rpt_label {comparison}_{gset}_shrunkenLFC -create_svgs false -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 1000 -set_min 10 -zip_report false -out {name} -gui false')

                exp.job_id.append(send_job(command_list=command_list,
                                           job_name=f'{comparison}_{name}_GSEA',
                                           job_log_folder=exp.job_folder,
                                           q='general',
                                           mem=3000,
                                           log_file=exp.log_file,
                                           project=exp.project))
                time.sleep(1)

    # Wait for jobs to finish
    job_wait(exp.job_id, exp.log_file)

    for comparison, design in exp.designs.items():
        for name, gset in gmts.items():
            path = glob.glob(f'{out_dir}/{comparison}/{name}/*')[0]
            if 'index.html' == f'{path}/index.html'.split('/')[-1]:
                new_dir = f'{out_dir}/{comparison}/{name}/'

                pos_file = glob.glob(f'{path}/gsea_report_for_na_pos*.xls')[0]
                neg_file = glob.glob(f'{path}/gsea_report_for_na_neg*.xls')[0]

                display(HTML(f'<h1>{comparison} GSEA Summary</h1>'))
                barplot = gsea_barplot(out_dir=new_dir, pos_file=pos_file, neg_file=neg_file, gmt_name=name)

                msg = f'Open "index.html" in subfolder for all results.\nOpen {barplot} for barplot summary of top enrichments.'

                with open(f'{new_dir}README.txt', 'w') as fp:
                    fp.write(msg)

                output(f'{comparison} GSEA enrichment barplot for {name} can be found here:\n {barplot}', exp.log_file)

            else:
                output(f'GSEA did not complete {name} for {comparison}.', exp.log_file)

    os.chdir(cwd)
    exp.tasks_complete.append('GSEA')
    output(f'GSEA analysis complete: {datetime.now():%Y-%m-%d %H:%M:%S}\n', exp.log_file)

    return exp


def plot_venn2(Series, string_name_of_overlap, folder):
    '''
    Series with with overlaps 10,01,11
    Plots a 2 way venn.
    Saves to file.
    '''

    os.makedirs(folder, exist_ok=True)
    plt.clf()
    plt.figure(figsize=(7, 7))

    font = {'family': 'sans-serif',
            'weight': 'normal',
            'size': 16,
            }

    plt.rc('font', **font)

    # make venn
    venn_plot = venn2(subsets=(Series.iloc[0], Series.iloc[1], Series.iloc[2]), set_labels=Series.index.tolist())
    patch = ['10', '01', '11']
    colors = ['green', 'blue', 'teal']
    for patch, color in zip(patch, colors):
        venn_plot.get_patch_by_id(patch).set_color('none')
        venn_plot.get_patch_by_id(patch).set_alpha(0.4)
        venn_plot.get_patch_by_id(patch).set_edgecolor('none')

    c = venn2_circles(subsets=(Series.iloc[0], Series.iloc[1], Series.iloc[2]))
    colors_circle = ['green', 'blue']
    for circle, color in zip(c, colors_circle):
        circle.set_edgecolor(color)
        circle.set_alpha(0.8)
        circle.set_linewidth(3)

    plt.title(f'{string_name_of_overlap.replace("_", " ")} Overlaps')
    plt.tight_layout()
    plt.savefig(f'{folder}{string_name_of_overlap}-overlap-{datetime.now():%Y-%m-%d}.svg')
    plt.savefig(f'{folder}{string_name_of_overlap}-overlap-{datetime.now():%Y-%m-%d}.png', dpi=300)
    out_result(f'{folder}{string_name_of_overlap}-overlap-{datetime.now():%Y-%m-%d}.png', f'Overlap Venn: {string_name_of_overlap.replace("_", " ")}')
    if run_main:
        plt.close()


def overlaps(exp):
    '''
    Performs overlaps of two or more de_sig lists.
    '''
    out_dir = f'{exp.scratch}Overlaps/'
    os.makedirs(out_dir, exist_ok=True)

    if len(exp.overlaps) != 0:
        names = ['2FC_UP', '2FC_DN', '15FC_UP', '15FC_DN', 'All_UP', 'All_DN']
        output(f'Beginning overlap of significant genes: {datetime.now():%Y-%m-%d %H:%M:%S}\n', exp.log_file)

        for overlap, comparison_list in exp.overlaps.items():
            if len(comparison_list) != 0:
                for name in names:
                    key = f'{overlap}_{name}'
                    exp.overlap_results[f'{key}_overlap'] = exp.sig_lists[comparison_list[0]][name] & exp.sig_lists[comparison_list[1]][name]
                    exp.overlap_results[f'{key}_uniqueA'] = exp.sig_lists[comparison_list[0]][name] - exp.sig_lists[comparison_list[1]][name]
                    exp.overlap_results[f'{key}_uniqueB'] = exp.sig_lists[comparison_list[1]][name] - exp.sig_lists[comparison_list[0]][name]

                    if len(exp.overlap_results[f'{key}_overlap']) == 0:
                        output(f'{overlap}_{name} have no overlapping genes', exp.log_file)
                    else:
                        venn = pd.Series([len(exp.sig_lists[comparison_list[0]][name]) - len(exp.overlap_results[f'{key}_overlap']),
                                          len(exp.sig_lists[comparison_list[1]][name]) - len(exp.overlap_results[f'{key}_overlap']),
                                          len(exp.overlap_results[f'{key}_overlap'])
                                          ],
                                         index=comparison_list + ['Overlap']
                                         )
                        plot_venn2(venn, key, f'{out_dir}{overlap}_{name}/')

    for name, sig in exp.overlap_results.items():
        sig_out = f'{out_dir}{name}_GEA/'
        os.makedirs(sig_out, exist_ok=True)

        if len(sig) == 0:
            output(f'Not performing GO enrichment for {name} overlaps since there are no overlapping genes.\n', exp.log_file)
        else:
            output(f'Performing GO enrichment for {name} overlaps: {datetime.now()} \n', exp.log_file)
            enrichr(gene_list=list(sig), description=f'{name}_overlap', out_dir=sig_out)

            with open(f'{sig_out}{name}.txt', 'w') as file:
                for gene in list(sig):
                    file.write(f'{gene}\n')

    exp.tasks_complete.append('Overlaps')
    output(f'Overlap analysis complete: {datetime.now():%Y-%m-%d %H:%M:%S}\n', exp.log_file)

    return exp


def plot_col(df, title, ylabel, out='', xy=(None, None), xticks=[''], plot_type=['violin', 'swarm'], pvalue=False, compare_tags=None, log_file=None):
    '''
    One or two column boxplot from dataframe.  Titles x axis based on column names.

    Inputs
    ------
    df: dataframe (uses first two columns)
    title: string of title
    ylabel: string of y label
    xy: If specified, will x is the label column and y is the data column. (default: (None,None): Data separated into two columns).
    xticks: list of xtick names (default is none)
    pvalue: bool to perform ttest (default is False).  Will only work if xy=(None,None) or ther are only two labels in x.
    plot_type: list of one or more: violin, box, swarm (default=violin)
    compare_tags:  if xy and pvalue is specified and there are more than two tags in x, specify the tags to compare. eg. ['a','b']
    out: out parent directory.  if none returns into colplot/
    log_file: log_file

    Returns
    ------
    None
    '''

    out = f'{val_folder(out)}/colplot/' if len(out) != 0 else 'colplot/'
    os.makedirs(out, exist_ok=True)

    plt.clf()
    sns.set(context='paper', font='Arial', font_scale=2, style='white', rc={'figure.dpi': 300, 'figure.figsize': (5, 6)})

    if type(plot_type) != list:
        plot_type = plot_type.split()
    lower_plot_type = [x.lower() for x in plot_type]

    if len(lower_plot_type) == 0:
        raise IOError('Input a plot type.')
    elif True not in {x in lower_plot_type for x in ['violin', 'box', 'swarm']}:
        raise IOError('Did not recognize plot type.')

    if 'swarm' in lower_plot_type:
        if xy == (None, None):
            fig = sns.swarmplot(data=df, color='black', s=4)
        else:
            fig = sns.swarmplot(data=df, x=xy[0], y=xy[1], color='black', s=4)
    if 'violin' in lower_plot_type:
        if xy == (None, None):
            fig = sns.violinplot(data=df)
        else:
            fig = sns.violinplot(data=df, x=xy[0], y=xy[1])
    if 'box' in lower_plot_type:
        if xy == (None, None):
            fig = sns.boxplot(data=df)
        else:
            fig = sns.boxplot(data=df, x=xy[0], y=xy[1])

    fig.yaxis.set_label_text(ylabel)
    fig.set_title(title)
    if xticks:
        fig.xaxis.set_ticklabels(xticks)
        fig.xaxis.set_label_text('')
        for tick in fig.xaxis.get_ticklabels():
            tick.set_fontsize(12)

    if pvalue:
        if xy == (None, None):
            _, pvalue = stats.ttest_ind(a=df.iloc[:, 0], b=df.iloc[:, 1])
            compare_tags = df.columns
        else:
            _, pvalue = stats.ttest_ind(a=df[df[xy[0]] == compare_tags[0]][xy[1]], b=df[df[xy[0]] == compare_tags[1]][xy[1]])
        fig.text(s=f'p-value = {pvalue:.03g}, {compare_tags[0]} v {compare_tags[1]}', x=0, y=-.12, transform=fig.axes.transAxes, fontsize=12)

    sns.despine()
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.17, top=0.9)
    plt.savefig(f'{out}{title.replace(" ", "_")}.png', dpi=300)
    if run_main:
        plt.close()

    out_result(f'{out}{title.replace(" ", "_")}.png', f'{title} Plot')
    output(f'{title.replace(" ", "_")}.png found in {out}', log_file)


def final_qc(exp):
    try:
        output(f'Beginning final qc: {datetime.now():%Y-%m-%d %H:%M:%S}\n', exp.log_file)

        command_list = ['module rm python', 'source activate RNAseq', f'cd {exp.scratch}', 'multiqc *']

        exp.job_id.append(send_job(command_list=command_list,
                                   job_name='MultiQC',
                                   job_log_folder=exp.job_folder,
                                   q='general',
                                   mem=1000,
                                   log_file=exp.log_file,
                                   project=exp.project))

        # Wait for jobs to finish
        job_wait(exp.job_id, exp.log_file)

        if os.path.isdir(f'{exp.scratch}multiqc_data'):
            copytree(f'{exp.scratch}multiqc_data', f'{exp.scratch}/QC/multiqc_data')
            rmtree(f'{exp.scratch}multiqc_data')

        log_file = None if run_main else exp.log_file
        samples = [sample for sample in exp.samples.values()]

        # Summary plots for RSEM alignment
        rsem_file = f'{exp.scratch}/QC/multiqc_data/multiqc_rsem.txt'
        if os.path.isfile(rsem_file):
            rsem_stats = read_pd(rsem_file)
            plot_col(df=rsem_stats.Alignable / 1e6,
                     title='Aligned Reads per Sample',
                     ylabel='Reads (Millions)',
                     log_file=log_file
                     )
            plot_col(df=rsem_stats.alignable_percent,
                     title='Percent Aligned per Sample',
                     ylabel='Percentage Aligned',
                     log_file=log_file
                     )

            '''
            fragment_series = pd.Series()
            for sample in samples:
                modelfile = '{}/QC/{}.models.pdf'.format(exp.scratch,sample)
                if os.path.isfile(modelfile):
                    with open(modelfile, 'rb') as fp:
                        text = PyPDF2.PdfFileReader(fp).getPage(0).extractText()
                        mean = float(text.split('\n')[1].split[','][1].split(' ')[-1])
                fragment_series[sample] = mean
            plot_col(df=fragment_series,
                     title='Mean Fragment Lengths per Sample',
                     ylable='Fragment Length',
                     log_file=log_file
                    )
            '''

        # Summary plots for FastQC data
        fastqc_file = f'{exp.scratch}/QC/multiqc_data/multiqc_fastqc.txt'
        if os.path.isfile(fastqc_file):
            gen_stats = read_pd(f'{exp.scratch}/QC/multiqc_data/multiqc_general_stats.txt')
            if exp.seq_type == 'paired':
                samples = [f'{sample}_R2' for sample in samples]
            plot_col(df=gen_stats.loc[samples, 'FastQC_mqc-generalstats-fastqc-total_sequences'] / 1e6,
                     title='Total Sequencer Reads per Sample',
                     ylabel='Reads (Millions)',
                     log_file=log_file
                     )

            plot_col(df=gen_stats.loc[samples, 'FastQC_mqc-generalstats-fastqc-percent_gc'],
                     title='Percent GC Content per Sample',
                     ylabel='Percentage of Reads with GC Content',
                     log_file=log_file
                     )

        display(HTML('<h1>Final QC Summary</h1>'))
        display(HTML(f'{exp.scratch}/multiqc_report.html'))

        exp.tasks_complete.append('MultiQC')

        return exp

    except:
        output('Error during MultiQC.', exp.log_file)
        filename = f'{exp.scratch}{exp.name}_incomplete.pkl'
        with open(filename, 'wb') as experiment:
            pickle.dump(exp, experiment)
        raise RuntimeError('Error during MultiQC. Fix problem then resubmit with same command to continue from last completed step.')


def finish(exp):

    job_wait(exp.job_id, exp.log_file)

    try:

        if os.path.isdir(f'{exp.scratch}/Fastq'):
            rmtree(f'{exp.scratch}/Fastq')

        output(f'\nConda environment file: {exp.job_folder}{exp.name}_environmnet.yml\nPackage versions: ', exp.log_file)

        os.system(f'conda env export > {exp.job_folder}{exp.name}_environmnet.yml')
        with open(f'{exp.job_folder}{exp.name}_environmnet.yml', 'r') as fp:
            versions = yaml.load(fp)
        for package in versions['dependencies']:
            output(package, exp.log_file)

        output(f'\n{exp.name} analysis complete! \n', exp.log_file)
        output(f'Copying all results into {exp.out_dir}: {datetime.now():%Y-%m-%d %H:%M:%S}\n', exp.log_file)

        scratch_log = f'{exp.scratch}{exp.log_file.split("/")[-1]}'
        if run_main:
            copy2(exp.log_file, scratch_log)

        rmtree(exp.out_dir)
        copytree(exp.scratch, exp.out_dir)

        exp.tasks_complete.append('Finished')

        filename = f'{exp.out_dir}{exp.name}_{exp.date}.pkl'
        with open(filename, 'wb') as experiment:
            pickle.dump(exp, experiment)

        output(f'Python Experiment: \n{exp}', exp.log_file)
        output(f'Moved all files into {exp.out_dir}: {datetime.now():%Y-%m-%d %H:%M:%S}\n', exp.log_file)
        output("\n Finger's Crossed!!!", exp.log_file)

        return exp

    except:
        output('Error while finishing pipeline.', exp.log_file)
        filename = f'{exp.scratch}{exp.name}_incomplete.pkl'
        with open(filename, 'wb') as experiment:
            pickle.dump(exp, experiment)
        raise RuntimeError('Error finishing pipeline. Fix problem then resubmit with same command to continue from last completed step.')


def validated_run(task, func, exp):
    try:
        if task in exp.tasks_complete:
            output(f'Skipping {task}...', exp.log_file)
            return exp
        else:
            return func(exp)
    except:
        output(f'Error in {task}.', exp.log_file)
        filename = f'{exp.scratch}{exp.name}_incomplete.pkl'
        with open(filename, 'wb') as experiment:
            pickle.dump(exp, experiment)
        raise RuntimeError(f'Error in {task}. Fix problem then resubmit with same command to continue from last completed step.')


def pipeline(experimental_file):
        exp = parse_yaml(experimental_file)
        exp = validated_run('Stage', stage, exp)
        exp = validated_run('Fastq_screen', fastq_screen, exp)
        exp = validated_run('Trim', trim, exp)
        exp = validated_run('FastQC', fastqc, exp)
        exp = validated_run('Spike', spike, exp)
        exp = validated_run('STAR', star, exp)
        exp = validated_run('RSEM', rsem, exp)
        exp = validated_run('Kallisto', kallisto, exp)
        exp = validated_run('GC', GC_normalization, exp)
        exp = validated_run('DESeq2', DESeq2, exp)
        exp = validated_run('Sleuth', Sleuth, exp)
        exp = validated_run('PCA', Principal_Component_Analysis, exp)
        exp = validated_run('Sigs', sigs, exp)
        exp = validated_run('Heatmaps', clustermap, exp)
        exp = validated_run('GO_enrich', GO_enrich, exp)
        exp = validated_run('GSEA', GSEA, exp)
        exp = validated_run('Overlaps', overlaps, exp)
        # exp =validated_run('decomp',decomposition,exp)
        exp = validated_run('MultiQC', final_qc, exp)
        exp = validated_run('Finished', finish, exp)


if run_main:
    import argparse
    plt.switch_backend('Agg')

    parser = argparse.ArgumentParser()
    parser.add_argument('--experimental_file', '-f', required=True, help='experimental yaml file', type=str)
    args = parser.parse_args()

    pipeline(args.experimental_file)
