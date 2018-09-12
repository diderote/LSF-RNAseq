#!/usr/bin/env python3

'''
University of Miami - Pegasus Cluster RNASeq Pipeline

Reads an experimental design yaml file (Version 0.7).
Requires a conda environment 'RNAseq' made from environment.yml
www.github.com/diderote/Nimerlab-RNAseq/

To do:
    - change designs to allow for non-bindary conditions ('Condition_A: a,a,b,b,c,c')
    - ICA with chi-square with de groups
    - t-SNE (add as option)
    - add sleuth for mouse
    - optimize STAR and updated rout_write
    - hiseq option
    - add papermill notebook implemenatation

'''
__author__ = 'Daniel L. Karl'
__license__ = 'MIT'
__version__ = '0.7'

import os
import re
import glob
import pickle
import math
import random
import time
from shutil import copy2,copytree,rmtree,move
from datetime import datetime
import subprocess as sub
import yaml
import reprlib

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib_venn import venn2, venn2_circles
import seaborn as sns
from sklearn.decomposition import PCA
import rpy2.robjects as ro
import rpy2.rinterface as ri
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri, r, globalenv, Formula
import gseapy


#make decorator validator for file endings for parser.
#def for pd read to read pandas

class Experiment:
    '''
    Experiment object for pipeline
    '''

    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        self.norm = 'Median-Ratios'
        self.spike = False
        self.de_sig_overlap = False
        self.gc_norm=False
        self.trim=[0,0]
        self.tasks_complete=[]
        self.job_id=[]
        self.designs={}
        self.overlaps = {}
        self.gene_lists={}
        self.de_results = {}
        self.sig_lists={}
        self.overlap_results={}
        self.genome_indicies={}

    def __repr__(self):
        exclude = ['overlaps','job_id','name','designs']
        experiment = '{}'.format(self.__dict__['name'])
        for key,value in self.__dict__.items():
            experiment += '\n\t{}: {}'.format(key,reprlib.repr(value)) if key not in exclude else ''
        return 'Experiment({})'.format(experiment)   
    
class RaiseError(Exception):
    pass
    
def read_pd(file):
    if (file.split('.')[-1] == 'txt') or (file.split('.')[-1] == 'tab'):
        return pd.read_table(file, header= 0, index_col=0)
    elif (file.split('.')[-1] == 'xls') or (file.split('.')[-1] == 'xlsx'):
        return pd.read_excel(file)
    else:
        raise IOError("Cannot parse count matrix.  Make sure it is .txt, .xls, or .xlsx")

def parse_yaml(experimental_file):
    '''
    Parse experimental info from yaml file
    '''    
    
    with open(experimental_file,'r') as file:
        yml = yaml.safe_load(file)

    #Make a new experimental object
    exp = Experiment()
    
    #Setting Scratch folder
    if yml['Lab'].lower() == 'nimer':
        exp.scratch = '/scratch/projects/nimerlab/DANIEL/staging/RNAseq/{}/'.format(yml['Name'])
    else:
        exp.scratch ='{}{}/'.format(yml['Scratch_folder'],yml['Name']) if yml['Scratch_folder'].endswith('/') else '{}/{}/'.format(yml['Scratch_folder'],yml['Name'])
    os.makedirs(exp.scratch, exist_ok=True)
    
    #check whether experiment has been attempted
    exp.name = yml['Name']   
    filename= '{}{}_incomplete.pkl'.format(exp.scratch,exp.name)
    
    if os.path.isfile(filename):
        if yml['Restart'] == False:
            with open(filename, 'rb') as experiment:
                exp = pickle.load(experiment)
            os.remove(filename)

            #set new date
            exp.date = format(datetime.now(), '%Y-%m-%d') 

            print('\n#############\nRestarting pipeline on {:%Y-%m-%d %H:%M:%S}, from last completed step.'.format(datetime.now()), file=open(exp.log_file,'a'))

            return exp 
        else:
            os.remove(filename)

    #Passing paramters to new object
    exp.date = format(datetime.now(), '%Y-%m-%d') 
    
    #Make out directory if it doesn't exist
    exp.out_dir = '{}{}/'.format(yml['Output_directory'],exp.name) if yml['Output_directory'].endswith('/') else '{}/{}/'.format(yml['Output_directory'], exp.name)
    os.makedirs(exp.out_dir, exist_ok=True)

    #Log file
    exp.log_file = '{}{}-{}.log'.format(exp.out_dir,exp.name,exp.date)
    
    print('Pipeline version {} run on {} \n'.format(str(__version__),exp.date), file=open(exp.log_file, 'w'))
    print('Beginning RNAseq Analysis: {:%Y-%m-%d %H:%M:%S}\n'.format(datetime.now()), file=open(exp.log_file, 'a'))
    print('Reading experimental file...\n', file=open(exp.log_file, 'a'))
    print("Pipeline output folder: {}\n".format(exp.out_dir), file=open(exp.log_file, 'a'))

    #Setting Job Folder
    exp.job_folder = '{}logs/'.format(exp.scratch)
    os.makedirs(exp.job_folder, exist_ok=True)

    #Genome
    exp.genome = yml['Genome'].lower() if yml['Genome'].lower() in ['hg38', 'mm10', 'hg19'] else None
    if exp.genome is None:
        raise ValueError("Genome must be either hg38, hg19, or mm10.")
    print('Processing data with: ' + str(exp.genome), file=open(exp.log_file, 'a'))

    #Sequencing type
    if yml['Sequencing_type'].lower() not in ['paired','single']:
        raise ValueError("Must specify whether sequence is paired or single end.")
    exp.seq_type = yml['Sequencing_type'].lower()
    exp.tasks_complete = exp.tasks_complete + ['Kallisto','Sleuth'] if exp.seq_type == 'single' else exp.tasks_complete

    #Standed
    exp.stranded = True if yml['Stranded'] else False
    print('Processing data as {}-end {} sequencing.'.format(exp.seq_type, ('stranded' if exp.stranded else 'non-stranded')), file =open(exp.log_file,'a'))

    #Tasks to complete
    if yml['Tasks']['Align'] == False:
        exp.tasks_complete = exp.tasks_complete + ['Stage','FastQC','Fastq_screen','Trim','STAR','Kallisto', 'Sleuth']
        print('Not performing alignment.', file=open(exp.log_file,'a'))
        count_matrix_loc=yml['Count_matrix']
        if os.path.exists(count_matrix_loc):
            print("Count matrix found at {}".format(count_matrix_loc), file=open(exp.log_file, 'a'))
            print("Performing only DESeq2 on for DE", file=open(exp.log_file, 'a'))
            exp.count_matrix = read_pd(count_matrix_loc)
        else:
            raise IOError("Count Matrix Not Found.")
    elif yml['Tasks']['Align'] == True:
        #Alignment mode. Default is transcript.
        exp.alignment_mode = 'gene' if yml['Tasks']['Alignment_Mode'].lower() == 'gene' else 'transcript'
    else:
        raise IOError('Please specify whether or not to perform alignment.')   
    
    #Lab specific files
    if yml['Lab'].lower() == 'other':
        exp.genome_indicies['RSEM_STAR'] = yml['RSEM_STAR_index']
        exp.genome_indicies['STAR'] = yml['STAR_index']
        exp.genome_indicies['Kallisto'] = yml['Kallisto_index']
        exp.genome_indicies['ERCC'] = yml['ERCC_STAR_index']
        exp.genome_indicies['GSEA_jar'] = yml['GSEA_jar']
        exp.genome_indicies['Gene_names'] = yml['Gene_names']
        if exp.genome == 'mm10':
            exp.genome_indicies['GMT'] = yml['GSEA_mouse_gmx_folder'] if yml['GSEA_mouse_gmx_folder'].endswith('/') else '{}/'.format(yml['GSEA_mouse_gmx_folder'])
    elif yml['Lab'].lower() == 'nimer':
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
            
    #GC_normalizaton
    if yml['Tasks']['GC_Normalization']:
        exp.gc_norm = True
    else:
        exp.tasks_complete.append('GC')

    #Support Files:
    if yml['Lab'].lower() == 'nimer':
        exp.genome_indicies['ERCC_Mix'] = '/projects/ctsi/nimerlab/DANIEL/tools/genomes/ERCC_spike/cms_095046.txt'
        if exp.genome == 'mm10':
            exp.genome_indicies['GC_Content'] = '/projects/ctsi/nimerlab/DANIEL/tools/genomes/Mus_musculus/mm10/mm10_GC_Content.txt'
        elif exp.genome == 'hg38':
            exp.genome_indicies['GC_Content'] = '/projects/ctsi/nimerlab/DANIEL/tools/genomes/H_sapiens/hg38_GC_Content.txt'
        elif exp.genome == 'hg19':
            exp.genome_indicies['GC_Content'] = '/projects/ctsi/nimerlab/DANIEL/tools/genomes/H_sapiens/hg38_GC_Content.txt'
    elif yml['Lab'].lower() == 'other':
        exp.genome_indicies['ERCC_Mix'] = yml['ERCC_Mix_file']
        exp.genome_indicies['GC_Content'] = yml['GC_Content_file']

    #No DE option
    if yml['Tasks']['Differential_Expression'] == False:
        exp.tasks_complete = exp.tasks_complete + ['GC','DESeq2','Sleuth','Sigs','Heatmaps','GO_enrich','GSEA','PCA']
        print('Not performing differential expression analyses.', file=open(exp.log_file,'a'))

    #Spike
    if yml['Tasks']['ERCC_align'] == False:
        exp.tasks_complete.append('Spike')   

    #Fastq Folder
    exp.fastq_folder= yml['Fastq_directory'] if yml['Fastq_directory'].endswith('/') else '{}/'.format(yml['Fastq_directory'])
    if (os.path.isdir(exp.fastq_folder) == False) and (yml['Tasks']['Align']):
        raise IOError("Can't Find Fastq Folder.")
    
    #Hard clip
    #exp.trim = yml['trim'] if yml['trim'] is not None else exp.trim

    #Project
    exp.project = 'nimerlab' if yml['Lab'].lower() == 'nimer' else yml['Pegasus_Project']
    
    #Sample Names
    exp.samples={key:name for key,name in yml['Samples'].items() if (name != 'Name') and (name is not None)}
    exp.sample_number = len(exp.samples)
    print("Samples: ", file=open(exp.log_file, 'a'))
    for number,sample in exp.samples.items():
        print('{number}: {sample}'.format(number=number,sample=sample), file=open(exp.log_file, 'a'))
    print('\nProcessing {} samples.\n'.format(exp.sample_number), file=open(exp.log_file, 'a'))

    #Differential Expression Groups
    if yml['Tasks']['Differential_Expression']:
        print("Parsing experimental design for differential expression...\n", file=open(exp.log_file, 'a'))
        
        exp.conditions = {key:condition.split(',') for key, condition in yml['Conditions'].items() if condition is not None}
        exp.de_tests = {key:test for key,test in yml['Designs'].items() if test['Test_condition'] is not None}
        for key,test in exp.de_tests.items():
            all_conditions = test['All_conditions'].split(',')
            all_samples = [exp.samples[int(sample)] for sample in test['All_samples'].split(',')]
            exp.designs[key] = {'all_samples': all_samples}
            for condition in all_conditions:
                exp.designs[key]['Condition_{}'.format(condition)] = [exp.samples[int(sample)] for group in condition for sample in exp.conditions[group]]
            exp.designs[key]['Test_condition'] = test['Test_condition'].split(',') 

            if len(exp.designs[key]['Test_condition']) == 1:
                exp.designs[key]['reduced'] = '~' + ' + '.join(['Condition_{}'.format(condition) for condition in all_conditions if condition not in exp.designs[key]['Test_condition']])
                exp.designs[key]['design'] = '{} + Condition_{}'.format(exp.designs[key]['reduced'],exp.designs[key]['Test_condition'][0]) if len(all_conditions) > 1 else '~Condition_{}'.format(exp.designs[key]['Test_condition'][0])  
            elif len(exp.designs[key]['Test_condition']) == 2:
                exp.designs[key]['reduced'] = '~' + ' + '.join(['Condition_{}'.format(condition) for condition in all_conditions])
                intersection='Condition_{}:Condition_{}'.format(exp.designs[key]['Test_condition'][0],exp.designs[key]['Test_condition'][1])
                exp.designs[key]['design'] = '{} + {}'.format(exp.designs[key]['reduced'],intersection) if len(all_conditions) > 2 else '~{}'.format(intersection)
            else:
                raise ValueError('Cannot handle this experimental design.')
            exp.designs[key]['reduced'] = '~1' if exp.designs[key]['reduced'] == '~' else exp.designs[key]['reduced']

            exp.designs[key]['colData'] = pd.DataFrame({'Condition_{}'.format(condition): ['condition' if sample in exp.designs[key]['Condition_{}'.format(condition)] else 'not_condition' for sample in all_samples] for condition in all_conditions},
                                                        index=all_samples) 
            exp.designs[key]['Test_type'] = exp.de_tests[key]['Test_type'].lower()
        
        for name,items in exp.designs.items():
            print('\n{}:'.format(name), file=open(exp.log_file,'a'))
            print(str(items['colData']), file=open(exp.log_file,'a'))


        #Normalization method
        if yml['Normalization'].lower() == 'ercc' or yml['Normalization'].lower() == 'ercc_mixed':
            exp.norm = yml['Normalization'].lower() 
            if yml['Tasks']['ERCC_align'] == False:
                spike_matrix_loc = yml['Spike_matrix']
                if os.path.exists(spike_matrix_loc):
                    print("Spike Count matrix found at {}".format(spike_matrix_loc), file=open(exp.log_file, 'a'))
                    exp.spike_counts = read_pd(spike_matrix_loc)
                else:
                    print("Cannot find spike matrix.", file=open(exp.log_file, 'a'))
            print('\nNormalizing samples for differential expression analysis using ERCC spike-in variance'+ '\n', file=open(exp.log_file, 'a'))
        elif yml['Normalization'].lower() == 'empirical':
            print('\nNormalizing samples for differential expression analysis using empirical negative controls for variance'+ '\n', file=open(exp.log_file, 'a'))
            exp.norm = 'empirical'
        elif yml['Normalization'].lower() == 'median-ratios':
            print('\nNormalizing samples for differential expression analysis using deseq2 size factors determined using default median of ratios method.'+ '\n', file=open(exp.log_file, 'a'))
        else:
            print("\nI don't know the {} normalization method.  Using default median-ratios.\n".format(yml['Normalization']), file=open(exp.log_file, 'a'))

    #Initialize DE sig overlaps
    exp.de_sig_overlap = True if yml['Tasks']['Signature_Mode'].lower() == 'combined' else False
    if exp.de_sig_overlap == False or exp.alignment_mode.lower() == 'transcript':
        exp.tasks_complete = exp.tasks_complete + ['Kallisto','Sleuth']
        
    #DE Overlaps
    for key, item in yml['Overlaps'].items():
        if bool(item):   
            exp.overlaps[key] = item.split('v')
    if str(len(list(exp.overlaps.keys()))) != 0:
        print('\nOverlapping ' + str(len(list(exp.overlaps.keys()))) + ' differential analysis comparison(s).', file=open(exp.log_file, 'a'))
        print(str(exp.overlaps)+ '\n', file=open(exp.log_file, 'a'))
    else:
        exp.tasks_complete.append('Overlaps')
        print('Not performing signature overlaps', file=open(exp.log_file,'a'))
    
    
    #Initialized Process Complete List
    exp.tasks_complete.append('Parsed')

    print('Experiment file parsed: {:%Y-%m-%d %H:%M:%S}\n'.format(datetime.now()), file=open(exp.log_file, 'a'))
    
    return exp

def send_job(command_list, job_name, job_log_folder, q, mem, log_file, project, threads=1):
    '''
    Sends job to LSF pegasus.ccs.miami.edu
    '''
    
    os.makedirs(job_log_folder, exist_ok=True)

    rand_id = str(random.randint(0, 100000))
    str_comd_list =  '\n'.join(command_list)
    cmd = '''
#!/bin/bash

#BSUB -J JOB_{job_name}_ID_{random_number}
#BSUB -R "rusage[mem={mem}]"
#BSUB -R "span[ptile={threads}]"
#BSUB -o {job_log_folder}{job_name_o}_logs_{rand_id}.stdout.%J
#BSUB -e {job_log_folder}{job_name_e}_logs_{rand_id}.stderr.%J
#BSUB -W 120:00
#BSUB -n {threads}
#BSUB -q {q}
#BSUB -P {project}

{commands_string_list}'''.format(job_name = job_name,
                                 job_log_folder=job_log_folder,
                                 job_name_o=job_name,
                                 job_name_e=job_name,
                                 commands_string_list=str_comd_list,
                                 random_number=rand_id,
                                 rand_id=rand_id,
                                 q=q,
                                 mem=mem,
                                 project=project,
                                 threads=threads
                                )
    
    job_path_name = job_log_folder + job_name+'.sh'
    write_job = open(job_path_name, 'w')
    write_job.write(cmd)
    write_job.close()
    os.system('bsub < {}'.format(job_path_name))
    print('sending job ID_{rand_id}...'.format(rand_id=str(rand_id)), file=open(log_file, 'a'))
    time.sleep(1) #too many conda activations at once sometimes leads to inability to activate during a job.
    
    return rand_id

def job_wait(id_list, job_log_folder, log_file):
    '''
    Waits for jobs sent by send job to finish.
    '''
    running = True
    while running:
        jobs_list = os.popen('sleep 60|bhist -w').read()
        current=[]
        for rand_id in id_list:
            if len([j for j in re.findall('ID_(\d+)', jobs_list) if j == rand_id]) != 0:
                current.append(rand_id)
        if len(current) == 0:
            running = False
        else:
            print('Waiting for jobs to finish... {:%Y-%m-%d %H:%M:%S}'.format(datetime.now()), file=open(log_file, 'a'))

#def rout_write(x, folder):
#    '''
#    function for setting r_out to print to file instead of jupyter
#    rpy2.rinterface.set_writeconsole_regular(rout_write)
#    rpy2.rinterface.set_writeconsole_warnerror(rout_write)
#    '''
#    print(x, file=open('{}/R_out_{:%Y-%m-%d}.txt'.format(folder, datetime.now()), 'a'))

def stage(exp):
    '''
    Stages files in Pegasus Scratch
    '''
    
    #Stage Experiment Folder in Scratch
    print('Staging in {}\n'.format(exp.scratch), file=open(exp.log_file, 'a'))
    
    #Copy Fastq to scratch fastq folder
    if os.path.exists('{}Fastq'.format(exp.scratch)):
        rmtree('{}Fastq'.format(exp.scratch))
    copytree(exp.fastq_folder, '{}Fastq'.format(exp.scratch))

    #change to experimental directory in scratch
    os.chdir(exp.scratch)
    
    exp.fastq_folder= '{}Fastq/'.format(exp.scratch)
    
    exp.tasks_complete.append('Stage')
    
    print('Staging complete: {:%Y-%m-%d %H:%M:%S}\n'.format(datetime.now()), file=open(exp.log_file, 'a'))

    return exp

def fastqc(exp):
    '''
    Performs fastq spec analysis with FastQC
    '''
    print('Assessing fastq quality. \n', file=open(exp.log_file, 'a'))

    #Make QC folder
    exp.qc_folder = exp.scratch + 'QC/'
    os.makedirs(exp.qc_folder, exist_ok=True)
    
        
    for number,sample in exp.samples.items():
        command_list = ['module rm python',
                        'module rm perl',
                        'source activate RNAseq',
                        'fastqc ' + exp.fastq_folder + sample + '*',
                       ]

        exp.job_id.append(send_job(command_list=command_list, 
                                   job_name= sample + '_fastqc',
                                   job_log_folder=exp.job_folder,
                                   q= 'general',
                                   mem=5000,
                                   log_file=exp.log_file,
                                   project=exp.project
                                  )
                         )

    #Wait for jobs to finish
    job_wait(id_list=exp.job_id, job_log_folder=exp.job_folder, log_file=exp.log_file)
    
    #move to qc folder
    fastqc_files = glob.glob(exp.fastq_folder + '*.zip')
    fastqc_files = fastqc_files + glob.glob(exp.fastq_folder + '*.html')
    for f in fastqc_files:
        copy2(f,exp.qc_folder)
        os.remove(f)
        
    exp.tasks_complete.append('FastQC')
    
    print('FastQC complete: {:%Y-%m-%d %H:%M:%S}\n'.format(datetime.now()), file=open(exp.log_file, 'a'))
    
    return exp

def fastq_screen(exp):
    '''
    Checks fastq files for contamination with alternative genomes using Bowtie2
    '''

    print('Screening for contamination during sequencing: {:%Y-%m-%d %H:%M:%S}\n'.format(datetime.now()), file=open(exp.log_file, 'a'))
    
    #Make QC folder
    exp.qc_folder = exp.scratch + 'QC/'
    os.makedirs(exp.qc_folder, exist_ok=True)

    #change to experimental directory in scratch
    os.chdir(exp.fastq_folder)
    
    if exp.seq_type == 'paired':
        fastq_end = '_R1.fastq.gz'
    elif exp.seq_type == 'single':
        fastq_end = '.fastq.gz'

    #Submit fastqc and fastq_screen jobs for each sample
    for number,sample in exp.samples.items():
        command_list = ['module rm python',
                        'module rm perl',
                        'source activate RNAseq',
                        'fastq_screen --threads 4 --aligner bowtie2 ' + exp.fastq_folder + sample + fastq_end
                       ]

        exp.job_id.append(send_job(command_list=command_list, 
                                   job_name= sample + '_fastq_screen',
                                   job_log_folder=exp.job_folder,
                                   q= 'general',
                                   mem=3000,
                                   log_file=exp.log_file,
                                   project=exp.project,
                                   threads=2
                                  )
                         )
        time.sleep(1)
    
    #Wait for jobs to finish
    job_wait(id_list=exp.job_id, job_log_folder=exp.job_folder, log_file=exp.log_file)
    
    #move to qc folder        
    fastqs_files = glob.glob(exp.fastq_folder + '*screen*')
    for f in fastqs_files:
        copy2(f,exp.qc_folder)
        os.remove(f)

    #change to experimental directory in scratch
    os.chdir(exp.scratch)
    exp.tasks_complete.append('Fastq_screen')
    print('Screening complete: {:%Y-%m-%d %H:%M:%S}\n'.format(datetime.now()), file=open(exp.log_file, 'a'))
    
    return exp

def trim(exp):
    '''
    Trimming based on standard UM SCCC Core Nextseq 500 technical errors.  Cudadapt can hard clip both ends, but may ignore 3' in future.
    '''

    print('Beginning fastq trimming: {:%Y-%m-%d %H:%M:%S}\n'.format(datetime.now()), file=open(exp.log_file, 'a'))
        
    #change to experimental directory in scratch
    os.chdir(exp.fastq_folder)
    
    scan=0
    while scan < 2:

        #Submit trimming files for each sample
        for number,sample in exp.samples.items():

            move_on = True
            if exp.seq_type == 'paired':
                if '{loc}{sample}_trim_R2.fastq.gz'.format(loc=exp.fastq_folder,sample=sample) in glob.glob(exp.fastq_folder + '*.gz'):
                    move_on = False
            elif exp.seq_type == 'single':
                if '{loc}{sample}_trim.fastq.gz'.format(loc=exp.fastq_folder,sample=sample) in glob.glob(exp.fastq_folder + '*.gz'):
                    move_on = False
            
            if move_on:
                print('Trimming {sample}: '.format(sample=sample), file=open(exp.log_file, 'a'))

                if exp.seq_type == 'paired':
                    trim_u=str(exp.trim[0])
                    trim_U=str(exp.trim[1])
                    cutadapt = 'cutadapt -j 4 -a AGATCGGAAGAGC -A AGATCGGAAGAGC --cores=10 --nextseq-trim=20 -u {trim_u} -u -{trim_u} -U {trim_U} -U -{trim_U} -m 18 -o {loc}{sample}_trim_R1.fastq.gz -p {loc}{sample}_trim_R2.fastq.gz {loc}{sample}_R1.fastq.gz {loc}{sample}_R2.fastq.gz'.format(qc=exp.qc_folder,loc=exp.fastq_folder,sample=sample,trim_u=trim_u,trim_U=trim_U)
                elif exp.seq_type == 'single':
                    cutadapt = 'cutadapt -j 4 -a AGATCGGAAGAGC --cores=10 --nextseq-trim=20 -m 18 -o {loc}{sample}_trim.fastq.gz {loc}{sample}.fastq.gz'.format(qc=exp.qc_folder,loc=exp.fastq_folder,sample=sample)
                
                command_list = ['module rm python',
                                'module rm perl',
                                'source activate RNAseq',
                                cutadapt
                               ]

                exp.job_id.append(send_job(command_list=command_list, 
                                           job_name= sample + "_trim",
                                           job_log_folder=exp.job_folder,
                                           q= 'general',
                                           mem=1000,
                                           log_file=exp.log_file,
                                           project=exp.project,
                                           threads=2
                                          )
                                 )
            
        #Wait for jobs to finish
        job_wait(id_list=exp.job_id, job_log_folder=exp.job_folder, log_file=exp.log_file)

        scan += 1
    
    #move logs to qc folder        
    print('\nTrimming logs are found in stdout files from bsub.  Cutadapt does not handle log files in multi-core mode.', file=open(exp.log_file, 'a'))

    for number,sample in exp.samples.items():
        if exp.seq_type == 'paired':
            if '{loc}{sample}_trim_R2.fastq.gz'.format(loc=exp.fastq_folder,sample=sample) not in glob.glob(exp.fastq_folder + '*.gz'):
                raise RaiseError('Not all samples were trimmed.')
        elif exp.seq_type == 'single':
            if '{loc}{sample}_trim.fastq.gz'.format(loc=exp.fastq_folder,sample=sample) not in glob.glob(exp.fastq_folder + '*.gz'):
                raise RaiseError('Not all samples were trimmed.')

    #change to experimental directory in scratch
    os.chdir(exp.scratch)

    exp.tasks_complete.append('Trim')
    print('Trimming complete: {:%Y-%m-%d %H:%M:%S}\n'.format(datetime.now()), file=open(exp.log_file, 'a'))

    return exp

def spike(exp, backend='Agg'):
    '''
    If calling from jupyter.  Change backend as needed.
    
    Align sequencing files to ERCC index using STAR aligner.
    '''
    plt.switch_backend(backend)

    print("Processing with ERCC spike-in: {:%Y-%m-%d %H:%M:%S}\n".format(datetime.now()), file=open(exp.log_file, 'a'))
        
    ERCC_folder=exp.scratch + 'ERCC/'
    os.makedirs(ERCC_folder, exist_ok=True)

    if not os.path.exists('{}Fastq'.format(exp.scratch)):
        copytree(exp.fastq_folder, '{}Fastq'.format(exp.scratch))
        exp.fastq_folder= '{}Fastq/'.format(exp.scratch)

    scan = 0
    while scan < 2:
        for number,sample in exp.samples.items():
            #Scan if succesful during second loop.
            if '{loc}{sample}_ERCCReadsPerGene.out.tab'.format(loc=ERCC_folder,sample=sample) not in glob.glob(ERCC_folder + '*.tab'):
                #Submit STAR alingment for spike-ins for each sample
                print('Aligning {sample} to spike-in.'.format(sample=sample)+ '\n', file=open(exp.log_file, 'a'))

                if os.path.isfile('{floc}{sample}_trim_R1.fastq.gz'.format(floc=exp.fastq_folder, sample=sample)) or os.path.isfile('{floc}{sample}_trim.fastq.gz'.format(floc=exp.fastq_folder, sample=sample)):
                    fname = '{floc}{sample}_trim'.format(floc=exp.fastq_folder, sample=sample)
                elif os.path.isfile('{floc}{sample}_R1.fastq.gz'.format(floc=exp.fastq_folder, sample=sample)) or os.path.isfile('{floc}{sample}.fastq.gz'.format(floc=exp.fastq_folder, sample=sample)):
                    fname = '{floc}{sample}'.format(floc=exp.fastq_folder, sample=sample)
                else:
                    print('Cannot find fastq files for spike-in alignment. \n', file=open(exp.log_file, 'a'))
                    raise IOError('Cannot find fastq files for spike-in alignment.')

                if exp.seq_type == 'paired':
                    spike='STAR --runThreadN 4 --genomeDir {index} --readFilesIn {fname}_R1.fastq.gz {fname}_R2.fastq.gz --readFilesCommand zcat --outFileNamePrefix {loc}{sample}_ERCC --quantMode GeneCounts'.format(index=exp.genome_indicies['ERCC'],fname=fname,loc=ERCC_folder,sample=sample)
                elif exp.seq_type == 'single':
                    spike='STAR --runThreadN 4 --genomeDir {index} --readFilesIn {fname}.fastq.gz --readFilesCommand zcat --outFileNamePrefix {loc}{sample}_ERCC --quantMode GeneCounts'.format(index=exp.genome_indicies['ERCC'],fname=fname,loc=ERCC_folder,sample=sample)

                command_list = ['module rm python',
                                'module rm perl',
                                'source activate RNAseq',
                                spike
                               ]

                exp.job_id.append(send_job(command_list=command_list, 
                                           job_name= sample + '_ERCC',
                                           job_log_folder=exp.job_folder,
                                           q= 'general',
                                           mem=5000,
                                           log_file=exp.log_file,
                                           project=exp.project,
                                           threads=2
                                          )
                                 )

        #Wait for jobs to finish
        job_wait(id_list=exp.job_id, job_log_folder=exp.job_folder, log_file=exp.log_file)

        scan += 1

    for number,sample in exp.samples.items():
        sam_file='{ERCC_folder}{sample}_ERCCAligned.out.sam'.format(ERCC_folder=ERCC_folder,sample=sample)
        if os.path.isfile(sam_file):
            os.remove(sam_file)

    print('Spike-in alignment jobs finished.', file=open(exp.log_file, 'a'))
    
    ### Generate one matrix for all spike_counts
    try:
        ERCC_counts = glob.glob(ERCC_folder + '*_ERCCReadsPerGene.out.tab')
        if len(ERCC_counts) != exp.sample_number:
            print('At least one ERCC alignment failed.', file=open(exp.log_file,'a'))
            raise RaiseError('At least one ERCC alignment failed. Check scripts and resubmit.')
        else:
            exp.spike_counts = pd.DataFrame(index=pd.read_csv(ERCC_counts[1], header=None, index_col=0, sep="\t").index)
        
            for number,sample in exp.samples.items():
                exp.spike_counts[sample] = pd.read_csv('{loc}{sample}_ERCCReadsPerGene.out.tab'.format(loc=ERCC_folder, sample=sample),header=None, index_col=0, sep="\t")[[3]]
            exp.spike_counts = exp.spike_counts.iloc[4:,:]
            exp.spike_counts.to_csv('{loc}ERCC.count.matrix.txt'.format(loc=ERCC_folder), header=True, index=True, sep="\t")

    except:
        print('Error generating spike_count matrix.', file=open(exp.log_file,'a'))
        raise RaiseError('Error generating spike_count matrix. Make sure the file is not empty.')
    
    #check to see if there were any spike in reads, if not, change
    if exp.spike_counts.loc['ERCC-00002',:].sum(axis=0) < 50:
        print('ERCC has low or no counts, skipping further spike-in analysis.', file=open(exp.log_file,'a'))
        return exp 

    if exp.genome_indicies['ERCC_Mix'] != None:
        # Filtering for counts with more than 5 counts in two samples
        spike_counts = exp.spike_counts.copy()
        spike_counts = spike_counts[spike_counts[spike_counts > 5].apply(lambda x: len(x.dropna()) > 1 , axis=1)]
        mix = read_pd(exp.genome_indicies['ERCC_Mix'])
        mix = mix.rename(columns={'concentration in Mix 1 (attomoles/ul)': 'Mix_1',
                                  'concentration in Mix 2 (attomoles/ul)': 'Mix_2'})
        names = list(spike_counts.columns)
        spike_counts = spike_counts.join(mix)
        
        merged_spike = pd.DataFrame(columns=['value','Mix_1','Mix_2','subgroup'])
        name = []
        length = len(spike_counts)
        for sample in names:
            merged_spike = pd.concat([merged_spike,
                                     spike_counts[[sample,'Mix_1','Mix_2','subgroup']].rename(columns={sample:'value'})],
                                    ignore_index=True)
            name=name + [sample]*length
        merged_spike['Sample']=name
        merged_spike['log'] = merged_spike.value.apply(lambda x: np.log2(x))
        merged_spike['log2_Mix_1']=np.log2(merged_spike.Mix_1)
        merged_spike['log2_Mix_2']=np.log2(merged_spike.Mix_2)

        # Plot ERCC spike.
        sns.set(context='paper', font_scale=2, style='white',rc={'figure.dpi': 300, 'figure.figsize':(6,6)})
        M1 = sns.lmplot(x='log2_Mix_1', y='log', hue='Sample', data=merged_spike, size=10, aspect=1)
        M1.set_ylabels(label='spike-in counts (log2)')
        M1.set_xlabels(label='ERCC Mix (log2(attamoles/ul))')
        plt.title("ERCC Mix 1 Counts per Sample")
        sns.despine()
        M1.savefig(ERCC_folder + 'ERCC_Mix_1_plot.png')
        if __name__ == "__main__":
            plt.close()

        sns.set(context='paper', font_scale=2, style='white',rc={'figure.dpi': 300, 'figure.figsize':(6,6)})
        M2 = sns.lmplot(x='log2_Mix_2', y='log', hue='Sample', data=merged_spike, size=10, aspect=1)
        M2.set_ylabels(label='spike-in counts (log2)')
        M2.set_xlabels(label='ERCC Mix (log2(attamoles/ul))')
        plt.title("ERCC Mix 2 Counts per Sample")
        sns.despine()
        M2.savefig(ERCC_folder + 'ERCC_Mix_2_plot.png')
        if __name__ == "__main__":
            plt.close()

        sns.set(context='paper', font_scale=2, style='white',rc={'figure.dpi': 300, 'figure.figsize':(6,6)})
        setB = sns.lmplot(x='log2_Mix_1', y='log', hue='Sample', data=merged_spike[merged_spike.subgroup == 'B'], size=10, aspect=1)
        setB.set_ylabels(label='spike-in counts (log2)')
        setB.set_xlabels(label='ERCC Mix (log2(attamoles/ul))')
        plt.title("ERCC Subgroup B Counts per Sample")
        sns.despine()
        setB.savefig(ERCC_folder + 'ERCC_Subgroup_B_plot.png')
        if __name__ == "__main__":
            plt.close()


    print("ERCC spike-in processing complete: {:%Y-%m-%d %H:%M:%S}\n".format(datetime.now()), file=open(exp.log_file, 'a'))
      
    exp.tasks_complete.append('Spike')
    return exp 

def bam2bw(in_bam,out_bw,job_log_folder,sample,project,stranded,log_file):
    
    if stranded:
        command_list = ['module rm python share-rpms65',
                        'source activate RNAseq',
                        'bamCoverage --filterRNAstrand forward -b {in_bam} --normalizeUsing CPM -bs 1 -o {out_bw}.cpm.fwd.bw'.format(in_bam,out_bw),
                        'bamCoverage --filterRNAstrand reverse --scaleFactor -1 -b {in_bam} --normalizeUsing CPM -bs 1 -o {out_bw}.cpm.rev.bw'.format(in_bam,out_bw)
                        ]
    else:
        command_list = ['module rm python share-rpms65',
                        'source activate RNAseq',
                        'bamCoverage -b {in_bam} --normalizeUsing CPM -bs 1 -o {out_bw}.cpm.bw'.format(in_bam,out_bw)
                        ]
    send_job(command_list=command_list,
             job_name= '{}_stranded_bw'.format(sample),
             job_log_folder=job_log_folder,
             q='general',
             mem='10000',
             threads=2,
             log_file=log_file,
             project=project
            )
                 

def star(exp):
    '''
    Alignment to genome using STAR
    '''
    print('\n Beginning genome alignments with STAR: {:%Y-%m-%d %H:%M:%S}'.format(datetime.now()), file=open(exp.log_file, 'a'))
    out_dir = '{}STAR_results/'.format(exp.scratch)
    os.makedirs(out_dir, exist_ok=True)
    os.chdir(out_dir) 

    scan=0
    while scan < 2: #Loop twice to make sure source activate didn't fail the first time
        for number,sample in exp.samples.items():      
            if '{}{}_Aligned.sortedByCoord.out.bam'.format(out_dir,sample) in glob.glob(out_dir + '*.bam'):
                pass
            else:
                print('Aligning using STAR to genome for {}.\n'.format(sample), file=open(exp.log_file, 'a'))
                fname = '{floc}{sample}_trim'.format(floc=exp.fastq_folder, sample=sample)

                if exp.seq_type == 'paired':
                    align='STAR --runThreadN 4 --outSAMtype BAM SortedByCoordinate --genomeDir {index} --readFilesIn {fname}_R1.fastq.gz {fname}_R2.fastq.gz --readFilesCommand zcat --outFileNamePrefix {loc}{sample}_ --quantMode GeneCounts'.format(index=exp.genome_indicies['STAR'],fname=fname,loc=out_dir,sample=sample)
                elif exp.seq_type == 'single':
                    align='STAR --runThreadN 4 --outSAMtype BAM SortedByCoordinate --genomeDir {index} --readFilesIn {fname}.fastq.gz --readFilesCommand zcat --outFileNamePrefix {loc}{sample}_ --quantMode GeneCounts'.format(index=exp.genome_indicies['STAR'],fname=fname,loc=out_dir,sample=sample)

                command_list = ['module rm python share-rpms65',
                                'source activate RNAseq',
                                align,
                                'samtools index {}_Aligned.sortedByCoord.out.bam'.format(sample)
                                ]

                exp.job_id.append(send_job(command_list=command_list, 
                                            job_name= sample + '_STAR',
                                            job_log_folder=exp.job_folder,
                                            q= 'bigmem',
                                            mem=50000,
                                            log_file=exp.log_file,
                                            project=exp.project,
                                            threads=2
                                            )
                                  )

                time.sleep(5)

        #Wait for jobs to finish
        job_wait(id_list=exp.job_id, job_log_folder=exp.job_folder, log_file=exp.log_file)

        scan += 1

    print('STAR alignment to genome finished.', file=open(exp.log_file, 'a'))

     #Generate signal files

    for sample in exp.samples.values():

        print('Generating bigwig singal file for {}.\n'.format(sample), file=open(exp.log_file,'a'))

        bam2bw(in_bam='{}{}_Aligned.sortedByCoord.out.bam'.format(out_dir,sample),
               out_bw='{}{}.star.'.format(out_dir,sample),
               job_log_folder=exp.job_folder,
               sample=sample,
               project=exp.project,
               log_file=exp.log_file,
               stranded=exp.stranded
              )
    
    ### Generate one matrix for all counts
    try:
        counts_glob = glob.glob(out_dir + '*_ReadsPerGene.out.tab')
        if len(counts_glob) != exp.sample_number:
            print('At least one STAR alignment failed.', file=open(exp.log_file,'a'))
            raise RaiseError('At least one STAR alignment failed. Check scripts and resubmit.')
        else:
            exp.count_matrix = pd.DataFrame(index=read_pd(counts_glob[1]).index)
        
            for number,sample in exp.samples.items():
                exp.count_matrix[sample] = read_pd('{}{}_ReadsPerGene.out.tab'.format(out_dir,sample))[[3]]
            exp.count_matrix = exp.count_matrix.iloc[4:,:]
            exp.count_matrix.to_csv('{}ALL_STAR.count.matrix.txt'.format(out_dir), header=True, index=True, sep="\t")
            if os.path.isfile(exp.genome_indicies['Gene_names']):
                with open(exp.genome_indicies['Gene_names'], 'rb') as file:
                    gene_dict = pickle.load(file)
                exp.count_matrix['name'] = exp.count_matrix.index
                exp.count_matrix = exp.count_matrix[exp.count_matrix.name.isin(gene_dict.keys())]
                exp.count_matrix['name'] = exp.count_matrix.name.apply(lambda x: '{}_{}'.format(x, gene_dict[x]))
                exp.count_matrix.index=exp.count_matrix.name
                exp.count_matrix = exp.count_matrix.drop(columns=['name'])
                exp.count_matrix.to_csv('{}Filtered_STAR.count.matrix.txt'.format(out_dir), header=True, index=True, sep="\t")

    except:
        print('Error generating count matrix.', file=open(exp.log_file,'a'))
        raise RaiseError('Error generating STARcount matrix. Make sure the file is not empty.')

    exp.tasks_complete.append('STAR')
    print('STAR alignemnt and count generation complete: {:%Y-%m-%d %H:%M:%S}\n'.format(datetime.now()), file=open(exp.log_file, 'a'))

    return exp

def rsem(exp):
    '''
    Alignment to transcriptome using STAR and estimating expected counts using EM with RSEM

    '''  
    print('\n Beginning RSEM-STAR transcriptome alignments: {:%Y-%m-%d %H:%M:%S}'.format(datetime.now()), file=open(exp.log_file, 'a'))
    
    out_dir = '{}RSEM_results/'.format(exp.scratch)
    
    os.makedirs(out_dir, exist_ok=True)
    os.chdir(out_dir)        

    scan=0
    while scan < 2: #Loop twice to make sure source activate didn't fail the first time
        for number,sample in exp.samples.items():      
            if '{}{}.genome.sorted.bam'.format(out_dir,sample) in glob.glob(out_dir + '*.bam'):
                pass
            else:
                print('Aligning using STAR and counting transcripts using RSEM for {}.\n'.format(sample), file=open(exp.log_file, 'a'))

                if exp.seq_type == 'paired':
                    align='rsem-calculate-expression --star --star-gzipped-read-file --paired-end --append-names --output-genome-bam --sort-bam-by-coordinate -p 4 {loc}{sample}_trim_R1.fastq.gz {loc}{sample}_trim_R2.fastq.gz {index} {sample}'.format(loc=exp.fastq_folder,index=exp.genome_indicies['RSEM_STAR'],sample=sample)
                elif exp.seq_type == 'single':
                    align= 'rsem-calculate-expression --star --star-gzipped-read-file --append-names --output-genome-bam --sort-bam-by-coordinate -p 4 {loc}{sample}_trim.fastq.gz {index} {sample}'.format(loc=exp.fastq_folder,index=exp.genome_indicies['RSEM_STAR'],sample=sample)

                plot_model='rsem-plot-model {sample} {sample}.models.pdf' .format(sample=sample)  

                command_list = ['module rm python share-rpms65',
                                'source activate RNAseq',
                                align,
                                plot_model
                                ]

                exp.job_id.append(send_job(command_list=command_list, 
                                            job_name= sample + '_RSEM',
                                            job_log_folder=exp.job_folder,
                                            q= 'bigmem',
                                            mem=50000,
                                            log_file=exp.log_file,
                                            project=exp.project,
                                            threads=2
                                            )
                                  )

                time.sleep(5)

        #Wait for jobs to finish
        job_wait(id_list=exp.job_id, job_log_folder=exp.job_folder, log_file=exp.log_file)
    
        scan += 1

    #Generate signal files

    for sample in exp.samples.values():

        print('Generating bigwig singal file for {}.\n'.format(sample), file=open(exp.log_file,'a'))

        bam2bw(in_bam='{}{}.genome.sorted.bam'.format(out_dir,sample),
               out_bw='{}{}.rsem.'.format(out_dir,sample),
               job_log_folder=exp.job_folder,
               sample=sample,
               project=exp.project,
               log_file=exp.log_file,
               stranded=exp.stranded
              )

    remove_files = ['genome.bam','transcript.bam','transcript.sorted.bam','transcrpt.sorted.bam.bai','wig']
    for number,sample in exp.samples.items():
        for file in remove_files:
            del_file='{}{}.{}'.format(out_dir, sample,file)
            if os.path.isfile(del_file):
                os.remove(del_file)
            pdf = '{}{}.models.pdf'.format(out_dir,sample)
            if os.path.isdir(exp.qc_folder) and os.path.isfile(pdf):
                move(pdf, '{}{}.models.pdf'.format(exp.qc_folder,sample))
    
    ### Generate one matrix for all expected_counts
    print('Generating Sample Matrix from RSEM.gene.results: {:%Y-%m-%d %H:%M:%S}\n'.format(datetime.now()), file=open(exp.log_file, 'a'))
    matrix='rsem-generate-data-matrix '
    columns=[]
    for number,sample in exp.samples.items():
        matrix = '{}{}{}.genes.results '.format(matrix,out_dir,sample)
        columns.append(sample)
        
    matrix = matrix + '> {}RSEM.count.matrix.txt'.format(out_dir)
        
    command_list = ['module rm python',
                    'source activate RNAseq',
                    matrix
                   ]

    exp.job_id.append(send_job(command_list=command_list, 
                               job_name= 'Generate_Count_Matrix',
                               job_log_folder=exp.job_folder,
                               q= 'general',
                               mem=1000,
                               log_file=exp.log_file,
                               project=exp.project
                              )
                     )
    
    #Wait for jobs to finish
    job_wait(id_list=exp.job_id, job_log_folder=exp.job_folder, log_file=exp.log_file)
    
    counts = read_pd('{}RSEM.count.matrix.txt'.format(out_dir))
    counts.columns = columns
    counts.to_csv('{}RSEM.count.matrix.txt'.format(out_dir), header=True, index=True, sep="\t")

    exp.count_matrix = counts
    exp.tasks_complete.append('STAR')
    print('STAR alignemnt and RSEM counts complete: {:%Y-%m-%d %H:%M:%S}\n'.format(datetime.now()), file=open(exp.log_file, 'a'))
    
    return exp
    
def kallisto(exp):
    '''
    Second/alternate alignment to transcriptome using kallisto
    '''
    #make Kallisto_results folder

    os.makedirs(exp.scratch + 'Kallisto_results/', exist_ok=True)

    scan = 0
    while scan < 2:

        print('Beginning Kallisto alignments: {:%Y-%m-%d %H:%M:%S}\n'.format(datetime.now()), file=open(exp.log_file, 'a'))

        #Submit kallisto for each sample
        for number,sample in exp.samples.items():

            kal_out = '{}Kallisto_results/{}/'.format(exp.scratch, sample)
            os.makedirs(kal_out, exist_ok=True)

            if '{loc}abundance.tsv'.format(loc=kal_out) not in glob.glob(kal_out + '*.tsv'):
                if exp.seq_type == 'paired':
                    align = 'kallisto quant --index={index} --output-dir={out} --threads=4 --bootstrap-samples=100 {loc}{sample}_trim_R1.fastq.gz {loc}{sample}_trim_R2.fastq.gz'.format(index=exp.genome_indicies['Kallisto'],out=kal_out,loc=exp.fastq_folder,sample=sample)

                print('Aligning {sample} using Kallisto.'.format(sample=sample)+ '\n', file=open(exp.log_file, 'a'))

                command_list = ['module rm python',
                                'module rm perl',
                                'source activate RNAseq',
                                align
                               ]

                exp.job_id.append(send_job(command_list=command_list, 
                                           job_name= sample + '_Kallisto',
                                           job_log_folder=exp.job_folder,
                                           q= 'general',
                                           mem=10000,
                                           log_file=exp.log_file,
                                           project=exp.project,
                                           threads=2
                                          )
                                 )

        #Wait for jobs to finish
        job_wait(id_list=exp.job_id, job_log_folder=exp.job_folder, log_file=exp.log_file)
            
        scan += 1

    exp.tasks_complete.append('Kallisto')
    
    return exp

def plot_PCA(counts, colData, out_dir, name, backend='Agg'):
    '''
    Inputs
    ------
    counts: dataframe of counts
    colData: dataframe of colData (DESeq2 format)
    out_dir: string of output directory
    name: name of PCA plot
    backend: Agg if not interactive (default), else use 'Qt5Agg'
    
    Outputs
    -------
    None

    prints PCA plot to out_dir

    '''
    plt.switch_backend(backend)

    try:
        to_remove=['gene_name','id', 'name']
        for x in to_remove:
            if x in list(counts.columns):
                counts = counts.drop(x, axis=1)

        pca = PCA(n_components=2)
        bpca = pca.fit_transform(counts.T)
        pca_score = pca.explained_variance_ratio_
        bpca_df = pd.DataFrame(bpca)
        bpca_df.index = counts.T.index
        bpca_df['name']= bpca_df.index

        fig = plt.figure(figsize=(8,8), dpi=300)
        ax = fig.add_subplot(111)
        if len(colData) == 0:
            ax.scatter(bpca_df[0], bpca_df[1], marker='o', color='black')
        else:
            bpca_df['group']= colData['main_comparison'].tolist()
            ax.scatter(bpca_df[bpca_df.group == 'Experimental'][0],bpca_df[bpca_df.group == 'Experimental'][1], marker='o', color='blue')
            ax.scatter(bpca_df[bpca_df.group == 'Control'][0],bpca_df[bpca_df.group == 'Control'][1], marker='o', color='red')
            red_patch = mpatches.Patch(color='red', alpha=.4, label='Control')
            blue_patch = mpatches.Patch(color='blue', alpha=.4, label='Experimental')

        ax.set_xlabel('PCA Component 1: {:0.1%} variance'.format(pca_score[0]))
        ax.set_ylabel('PCA Component 2: {:0.1%} varinace'.format(pca_score[1]))


        for i,sample in enumerate(bpca_df['name'].tolist()):
            xy=(bpca_df.iloc[i,0], bpca_df.iloc[i,1])
            xytext=tuple([sum(x) for x in zip(xy, ((sum(abs(ax.xaxis.get_data_interval()))*.01),(sum(abs(ax.yaxis.get_data_interval()))*.01)))])
            ax.annotate(sample, xy= xy, xytext=xytext)             
        
        if len(colData) != 0:
            ax.legend(handles=[blue_patch, red_patch], loc=1)
        
        sns.despine()
        plt.tight_layout()
        plt.subplots_adjust(right=0.8, top=.8)

        os.makedirs(out_dir, exist_ok=True)
        ax.figure.savefig('{}{}_PCA.png'.format(out_dir,name))
        ax.figure.savefig('{}{}_PCA.svg'.format(out_dir,name))
        if __name__ == "__main__":
            plt.close()
    except:
        raise RaiseError('Error during plot_PCA. Fix problem then resubmit with same command to continue from last completed step.')

def GC_normalization(exp):
    '''
    Within lane loess GC normalization using EDAseq
    '''
    pandas2ri.activate()

    #ri.set_writeconsole_regular(rout_write(folder=exp.job_folder))
    #ri.set_writeconsole_warnerror(rout_write(folder=exp.job_folder))

    edaseq = importr('EDASeq')
    as_df=ro.r("as.data.frame")
    assay=ro.r("assay")
    as_cv = ro.r('as.character')
    counts = ro.r("counts")
    fdata=ro.r('fData')
    normCounts=ro.r('normCounts')

    print('Beginning within-lane GC length/content loess normalization for all samples: {:%Y-%m-%d %H:%M:%S}\n'.format(datetime.now()), file=open(exp.log_file,'a'))

    GC_content = read_pd(exp.genome_indicies['GC_Content'])
    raw_counts = exp.count_matrix
    raw_counts['id']=raw_counts.index
    raw_counts['id']=raw_counts.id.apply(lambda x: x.split("_")[0].split(".")[0])
    GC_genes = set(GC_content.split.tolist())

    #Keep only counts with GC data (based on latest ensembl biomart).  see EDAseq package and use 'biomart' after dropping ensembl name version.
    GC_counts = round(raw_counts[raw_counts.id.apply(lambda x: x in GC_genes)])
    GC_gene_set = set(GC_counts.id.tolist())
    GC_content = GC_content[GC_content.split.apply(lambda x: x in GC_gene_set)]
    EDA_set = edaseq.newSeqExpressionSet(counts=GC_counts.drop(columns='id').values,featureData=GC_content)
    gcNorm = edaseq.withinLaneNormalization(EDA_set, 'gc','loess')
    data_norm = ro.pandas2ri.ri2py_dataframe(normCounts(gcNorm))
    data_norm.index = GC_counts.index
    data_norm.columns = GC_counts.drop(columns='id').columns
    exp.gc_count_matrix = data_norm

    print('Finished GC normalization: {:%Y-%m-%d %H:%M:%S}\n'.format(datetime.now()), file=open(exp.log_file,'a'))
    exp.tasks_complete.append('GC')

    return exp 

def RUV(RUV_data,test_type,design,reduced,colData,norm_type,log, ERCC_counts, comparison, plot_dir,de):

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
    '''
    try:
        pandas2ri.activate()
        
        #ri.set_writeconsole_regular(rout_write(folder=exp.job_folder))
        #ri.set_writeconsole_warnerror(rout_write(folder=exp.job_folder))

        deseq = importr('DESeq2')
        ruvseq = importr('RUVSeq')
        edaseq = importr('EDASeq')
        as_df=ro.r("as.data.frame")
        assay=ro.r("assay")
        as_cv = ro.r('as.character')
        counts = ro.r("counts")
        normCounts=ro.r('normCounts')
        pdata=ro.r('pData')
        results = ro.r('results')
        head = ro.r('head')

        os.makedirs(plot_dir, exist_ok=True)

        if comparison == 'ALL':
            plot_PCA(counts = RUV_data, colData=[], out_dir= plot_dir, name= '{}_pre-{}_RUV_ALL_counts'.format(comparison, norm_type))
        else:
            plot_PCA(counts = RUV_data, colData=colData, out_dir= plot_dir, name= '{}_pre-{}_RUV_counts'.format(comparison, norm_type))

        #retain gene name
        RUV_data['name'] = RUV_data.index

        #RUVseq
        if norm_type.lower() == 'empirical':
            print('Performing Normalization by removing unwatned variance of empirical negative control genes for {}: {:%Y-%m-%d %H:%M:%S}\n'.format(comparison,datetime.now()) , file=open(log,'a'))
            
            #determining non differentially expressed genes to use as empirical negative controls
            dds_emp = deseq.DESeqDataSetFromMatrix(countData = RUV_data.drop(columns='name').values,
                                                   colData=colData,
                                                   design=ro.Formula(design)
                                                  )
            dds_emp = deseq.DESeq(dds_emp)
            results_emp = pandas2ri.ri2py(as_df(deseq.results(dds_emp)))
            results_emp.index=RUV_data.index
            results_emp.sort_values(by='padj', inplace=True)
            top_de = list(results_emp.head(10000).index)
            
            #rename indices to reflect rpy2 conversion to R dataframe
            RUV_data.index= range(1,(len(RUV_data)+1))
            
            #empirical negative controls
            empirical = list(RUV_data[RUV_data.name.apply(lambda x: x not in top_de)].drop(columns='name').index)
            
            #generate normalization scaling based on unwanted variance from empirical negative controls
            data_set = edaseq.newSeqExpressionSet(RUV_data.drop(columns='name').values, phenoData=colData)
            RUVg_set = ruvseq.RUVg(x=data_set, cIdx=as_cv(empirical), k=1)

            print('\nEmpirical negative control normalization complete for {}: {:%Y-%m-%d %H:%M:%S}\n'.format(comparison,datetime.now()), file=open(log, 'a'))

        elif norm_type.lower() == 'ercc':
            print('Performing Normalization by removing unwanted variance using ERCC spike-ins for {}: {:%Y-%m-%d %H:%M:%S}\n'.format(comparison,datetime.now()), file=open(log,'a'))
            
            #rename ERCC join ERCC counts to gene counts and reindex for rpy2 R dataframe
            ERCC_counts['name'] = ERCC_counts.index
            ERCC_counts['name'] = ERCC_counts.name.apply(lambda x: '{}_{}'.format(x,x))

            RUV_data = RUV_data.append(ERCC_counts).dropna() ####added dropna
            RUV_data.index= range(1,(len(RUV_data)+1))
            
            #generate index locations of ERCC spikes
            spike_list = list(RUV_data[RUV_data.name.apply(lambda x: x in list(ERCC_counts.name))].index)
            
            #normalize samples based on unwanted variance between ERCC spike in controls
            data_set = edaseq.newSeqExpressionSet(RUV_data.drop(columns='name').values, phenoData=colData)
            RUVg_set = ruvseq.RUVg(x=data_set, cIdx=as_cv(spike_list), k=1)
            print('\nERCC normalization complete for {}: {:%Y-%m-%d %H:%M:%S}\n'.format(comparison, datetime.now()), file=open(log, 'a'))

        else:
            raise RaiseError('RUV() takes only "ercc" or "empirical" as options.')

        #Differential expression (wald or lrt) to account for scaled variances between samples
        if comparison == 'ALL':
            RUV_dds = deseq.DESeqDataSetFromMatrix(countData=counts(RUVg_set), colData=pdata(RUVg_set), design=ro.Formula('~W_1'))
            RUV_dds = deseq.estimateSizeFactors_DESeqDataSet(RUV_dds)
        else:
            RUV_dds = deseq.DESeqDataSetFromMatrix(countData=counts(RUVg_set), colData=pdata(RUVg_set), design=ro.Formula('~W_1 + {}'.format(design.split('~')[-1])))
            RUV_dds = deseq.DESeq(RUV_dds)
            if test_type == 'lrt':
                RUV_dds = deseq.DESeq(RUV_dds,test='LRT',reduced=ro.Formula('~W1' if reduced == '~1' else '~W1 + {}'.format(reduced.split('~')[-1])))

        #generate normalized counts for pca
        counts_df = pandas2ri.ri2py(as_df(normCounts(RUVg_set)))
        counts_df.columns = RUV_data.drop(columns='name').columns

        if comparison == 'ALL':
            plot_PCA(counts = counts_df.dropna(), colData=[], out_dir=plot_dir, name='{}_post-{}_RUV_raw_counts'.format(comparison, norm_type))
        else:
            plot_PCA(counts = counts_df.dropna(), colData= colData, out_dir=plot_dir, name='{}_post-{}_RUV_raw_counts'.format(comparison, norm_type))

        if de:
            #extract results and relabel samples and genes
            print('{} results type: '.format(comparison), file = open(log, 'a'))

            results = pandas2ri.ri2py(as_df(deseq.results(RUV_dds, contrast=as_cv(['{}'.format(design.split(' ')[-1]),'condition','not_condition']))))
            results.index = RUV_data.name

            lfc = pandas2ri.ri2py(as_df(deseq.lfcShrink(RUV_dds, contrast=as_cv(['{}'.format(design.split(' ')[-1]),'condition','not_condition']), type='ashr')))
            lfc.index = RUV_data.name
            print('Switched to ashr method for lfcShrinkage for {} with RUV normalizaiton.'.format(comparison), file=open(log,'a'))

        RUV_normcounts = pandas2ri.ri2py(as_df(counts(RUV_dds, normalized=True)))
        RUV_normcounts.columns = RUV_data.drop(columns='name').columns
        RUV_normcounts.index = RUV_data.name
        if comparison == 'ALL':
            plot_PCA(counts = RUV_normcounts, colData=[], out_dir= plot_dir, name= '{}_post-{}_RUV_normalized_counts'.format(comparison, norm_type))
        else:
            plot_PCA(counts = RUV_normcounts, colData= colData, out_dir=plot_dir, name='{}_post-{}_RUV_normalized_counts'.format(comparison, norm_type))


        rld = pandas2ri.ri2py_dataframe(assay(deseq.rlog(RUV_dds, blind=False)))
        rld.columns = RUV_data.drop(columns='name').columns
        rld.index = RUV_data.name

        print('Unwanted variance normalization complete for {} using RUVSeq: {:%Y-%m-%d %H:%M:%S}'.format(comparison, datetime.now()), file=open(log,'a'))

        if de:
            return results, rld, lfc, RUV_normcounts
        else:
            return rld, RUV_normcounts

    except: 
        raise RaiseError('Error during RUVseq.')

def DESeq2(exp):
        
    '''
    Differential Expression using DESeq2
    
    Inputs
    ------


    '''

    print('Beginning DESeq2 differential expression analysis: {:%Y-%m-%d %H:%M:%S}\n'.format(datetime.now()), file=open(exp.log_file, 'a'))
    
    pandas2ri.activate()
    
    #ri.set_writeconsole_regular(rout_write(folder=exp.job_folder))
    #ri.set_writeconsole_warnerror(rout_write(folder=exp.job_folder))

    deseq = importr('DESeq2')
    as_df=ro.r("as.data.frame")
    as_cv=ro.r('as.character')
    assay=ro.r("assay")
    session=ro.r("sessionInfo")
    results = ro.r('results')
    counts = ro.r("counts")
    head = ro.r('head')
    
    out_dir= exp.scratch + 'DESeq2_results/'
    os.makedirs(out_dir, exist_ok=True)
    
    if exp.gc_norm:
        if exp.alignment_mode == 'gene':
            print('Using GC normalized counts for differential expression.\n', file=open(exp.log_file,'a'))
        else:
            print('Using GC normalized RSEM expected counts for differential expression.\n', file=open(exp.log_file,'a'))
        count_matrix = exp.gc_count_matrix
    else:
        if exp.alignment_mode == 'gene':
            print('Using STAR counts for differential expression.\n', file=open(exp.log_file,'a'))
        else:
            print('Using rounded RSEM expected counts for differential expression.\n', file=open(exp.log_file,'a'))
        count_matrix = round(exp.count_matrix)
    
    dds={}
    
    for comparison,designs in exp.designs.items():
        print('Beginning {}: {:%Y-%m-%d %H:%M:%S}\n'.format(comparison, datetime.now()), file=open(exp.log_file, 'a'))
        colData=designs['colData']
        design=ro.Formula(designs['design'])
        data=count_matrix[designs['all_samples']]

        # filtering for genes with more than 1 count in two samples
        data = round(data[data[data > 1].apply(lambda x: len(x.dropna()) > 1 , axis=1)]) 

        dds[comparison] = deseq.DESeqDataSetFromMatrix(countData = data.values,
                                                       colData=colData,
                                                       design=design
                                                      )

        if exp.norm.lower() == 'ercc':
            print('Determining ERCC scaling vs Sample scaling using median of ratios of counts for rough comparison.  This may point out potentially problematic samples.\n', file=open(exp.log_file, 'a'))
            ERCC_data = round(exp.spike_counts[designs['all_samples']])
            ERCC_dds = deseq.DESeqDataSetFromMatrix(countData = ERCC_data.values, colData=colData, design=design)
            ERCC_size = deseq.estimateSizeFactors_DESeqDataSet(ERCC_dds)
            deseq2_size = deseq.estimateSizeFactors_DESeqDataSet(dds[comparison])
            sizeFactors=ro.r("sizeFactors")
            
            #Legacy:  Do not scale by ERCC size factors using DESeq2.
            #dds[comparison].do_slot('colData').do_slot('listData')[1] = sizeFactors(ERCC_size)
            #dds[comparison] = deseq.DESeq(dds[comparison])

            #compare size factors from DESeq2 and ERCC for inconsistencies
            ERCC_vector=pandas2ri.ri2py_vector(sizeFactors(ERCC_size))
            deseq2_vector=pandas2ri.ri2py_vector(sizeFactors(deseq2_size))
            if len(ERCC_vector) == len(deseq2_vector):
                for x in range(len(ERCC_vector)):
                    if abs((ERCC_vector[x]-deseq2_vector[x])/(ERCC_vector[x]+deseq2_vector[x])) > 0.1:
                        print('ERCC spike ({x} in list) is greater than 10 percent different than deseq2 size factor for {comparison}. \n'.format(x=x+1,comparison=comparison), file=open(exp.log_file,'a'))
                print('Samples: {}\n'.format(str(designs['all_samples'])), file=open(exp.log_file,'a'))
                print('ERCC size factors: {}'.format(str(ERCC_vector)), file=open(exp.log_file,'a'))
                print('DESeq2 size factors: {}\n'.format(str(deseq2_vector)), file=open(exp.log_file,'a'))
            else:
                print('\nERCC and deseq2 column lengths are different for {}'.format(comparison), file=open(exp.log_file,'a'))

        #Differential Expression
        if exp.norm.lower() == 'median-ratios':
            print('Using DESeq2 standard normalization of scaling by median of the ratios of observed counts.', file=open(exp.log_file, 'a'))
            
            print('Performing {} test for differential expression for {}\n'.format(design['Test_type'],comparison), file=open(exp.log_file, 'a'))
            dds[comparison] = deseq.DESeq(dds[comparison])
            
            if designs['Test_type'] == 'lrt':
                reduced = ro.Formula(designs['reduced'])
                dds[comparison] = deseq.DESeq(dds[comparison],test = 'LRT',reduced=reduced)
            
            print('{} results type: '.format(comparison), file = open(exp.log_file, 'a'))

            #get results
            exp.de_results['DE2_{}'.format(comparison)] = pandas2ri.ri2py(as_df(deseq.results(dds[comparison], contrast=as_cv(['{}'.format(design.split(' ')[-1]),'condition','not_condition']))))
            exp.de_results['DE2_{}'.format(comparison)].index = data.index
            
            #get shrunken lfc (apeglm) method)
            exp.de_results['shrunkenLFC_{}'.format(comparison)] = pandas2ri.ri2py(as_df(deseq.lfcShrink(dds[comparison], coef=as_cv('{}_{}_vs_{}'.format('{}'.format(design.split(' ')[-1]),'condition','not_condition')), type='apeglm')))
            exp.de_results['shrunkenLFC_{}'.format(comparison)].index = data.index
            print('Using apeglm method for lfc shrinkage for {}.'.format(comparison), file=open(exp.log_file, 'a'))

            #regularized log transformed
            exp.de_results['{}_rlog_counts'.format(comparison)] = pandas2ri.ri2py_dataframe(assay(deseq.rlog(dds[comparison], blind=False)))
            exp.de_results['{}_rlog_counts'.format(comparison)].columns = data.columns
            exp.de_results['{}_rlog_counts'.format(comparison)].index = data.index

            #normalized counts
            exp.de_results['{}_DE2_normCounts'.format(comparison)] = pandas2ri.ri2py_dataframe(counts(dds[comparison], normalized=True))
            exp.de_results['{}_DE2_normCounts'.format(comparison)].columns = data.columns
            exp.de_results['{}_DE2_normCounts'.format(comparison)].index = data.index

        elif exp.norm.lower() == 'ercc':
            exp.de_results['DE2_{}'.format(comparison)],exp.de_results['{}_rlog_counts'.format(comparison)], exp.de_results['shrunkenLFC_{}'.format(comparison)], exp.de_results['{}_DE2_normCounts'.format(comparison)]  = RUV(RUV_data = data, 
                                                                          test_type=designs['Test_type'],
                                                                          design=designs['design'],
                                                                          reduced=designs['reduced'],
                                                                          colData=colData, 
                                                                          norm_type='ERCC', 
                                                                          ERCC_counts = exp.spike_counts[designs['all_samples']], 
                                                                          log=exp.log_file,
                                                                          comparison=comparison,
                                                                          plot_dir = '{}PCA/'.format(exp.scratch),
                                                                          de= True
                                                                         )

        elif exp.norm.lower() == 'ercc_mixed':
            full_counts = exp.spike_counts[designs['all_samples']]
            mix = read_pd(exp.genome_indicies['ERCC_Mix'])
            subgroupB=mix[mix.subgroup == 'B'].index.tolist()
            exp.de_results['DE2_{}'.format(comparison)],exp.de_results['{}_rlog_counts'.format(comparison)], exp.de_results['shrunkenLFC_{}'.format(comparison)], exp.de_results['{}_DE2_normCounts'.format(comparison)]  = RUV(RUV_data = data,
                                                                          test_type=designs['Test_type'],
                                                                          design=designs['design'],
                                                                          reduced=designs['reduced'],
                                                                          colData=colData, 
                                                                          norm_type='ERCC', 
                                                                          ERCC_counts = full_counts.loc[subgroupB], 
                                                                          log=exp.log_file,
                                                                          comparison=comparison,
                                                                          plot_dir = '{}PCA/'.format(exp.scratch),
                                                                          de= True
                                                                         )
    
        elif exp.norm.lower() == 'empirical':
            exp.de_results['DE2_{}'.format(comparison)],exp.de_results['{}_rlog_counts'.format(comparison)], exp.de_results['shrunkenLFC_{}'.format(comparison)], exp.de_results['{}_DE2_normCounts'.format(comparison)]  = RUV(RUV_data = data,
                                                                          test_type=designs['Test_type'],
                                                                          design=designs['design'], 
                                                                          reduced=designs['reduced'],
                                                                          colData=colData, 
                                                                          norm_type='empirical', 
                                                                          ERCC_counts = None, 
                                                                          log=exp.log_file,
                                                                          comparison=comparison,
                                                                          plot_dir ='{}PCA/'.format(exp.scratch),
                                                                          de=True
                                                                         )
        else:
            raise RaiseError('Can only use "median-ratios", "ercc", "ercc_mixed" or "empirical" for normalization of DESeq2.')

        #DESeq2 results
        exp.de_results['DE2_{}'.format(comparison)].sort_values(by='padj', ascending=True, inplace=True)
        exp.de_results['DE2_{}'.format(comparison)]['gene_name']=exp.de_results['DE2_{}'.format(comparison)].index
        exp.de_results['DE2_{}'.format(comparison)]['gene_name']=exp.de_results['DE2_{}'.format(comparison)].gene_name.apply(lambda x: x.split("_")[1])
        exp.de_results['DE2_{}'.format(comparison)].to_csv('{}{}-DESeq2-results.txt'.format(out_dir,comparison), 
                                                           header=True, 
                                                           index=True, 
                                                           sep="\t"
                                                          )
        ##Shrunken LFC using apeglm or ashr method
        exp.de_results['shrunkenLFC_{}'.format(comparison)].sort_values(by='log2FoldChange', ascending=False, inplace=True)
        exp.de_results['shrunkenLFC_{}'.format(comparison)]['gene_name']=exp.de_results['shrunkenLFC_{}'.format(comparison)].index
        exp.de_results['shrunkenLFC_{}'.format(comparison)]['gene_name']=exp.de_results['shrunkenLFC_{}'.format(comparison)].gene_name.apply(lambda x: x.split("_")[1])
        exp.de_results['shrunkenLFC_{}'.format(comparison)].to_csv('{}{}-DESeq2-shrunken-LFC.txt'.format(out_dir,comparison), 
                                                                   header=True, 
                                                                   index=True, 
                                                                   sep="\t"
                                                                  )

        #Regularized log2 expected counts.
        exp.de_results['{}_rlog_counts'.format(comparison)].to_csv('{}{}-rlog-counts.txt'.format(out_dir,comparison), 
                                                                   header=True, 
                                                                   index=True, 
                                                                   sep="\t"
                                                                  )

        #Normalized counts (not-log transformed)
        exp.de_results['{}_DE2_normCounts'.format(comparison)].to_csv('{}{}-normalized-counts.txt'.format(out_dir,comparison),
                                                                      header=True,
                                                                      index=True,
                                                                      sep="\t"
                                                                     )

    #rlog count matrix for all samples.
    colData = pd.DataFrame(index=count_matrix.columns, data={'condition': ['A']*exp.sample_number})
    design=ro.Formula("~1")
    #count_matrix = round(count_matrix[count_matrix[count_matrix > 5].apply(lambda x: len(x.dropna()) > 1 , axis=1)]) 
    dds_all = deseq.DESeqDataSetFromMatrix(countData = count_matrix.values,
                                           colData=colData,
                                           design=design
                                          )
    exp.de_results['all_rlog'] = pandas2ri.ri2py_dataframe(assay(deseq.rlog(dds_all)))
    exp.de_results['all_rlog'].index=count_matrix.index
    exp.de_results['all_rlog'].columns=count_matrix.columns
    exp.de_results['all_rlog']['gene_name']=exp.de_results['all_rlog'].index
    exp.de_results['all_rlog']['gene_name']=exp.de_results['all_rlog'].gene_name.apply(lambda x: x.split("_")[1])
    exp.de_results['all_rlog'].to_csv('{}ALL-samples-blind-rlog-counts.txt'.format(out_dir), 
                                      header=True, 
                                      index=True, 
                                      sep="\t"
                                     )

    if exp.norm.lower() == 'ercc':
        exp.de_results['all_ERCC_rlog'],exp.de_results['all_ERCC_normCounts']  = RUV(RUV_data = count_matrix, 
                                                     test_type='',
                                                     design=design,
                                                     reduced=design,
                                                     colData=colData, 
                                                     norm_type='ERCC', 
                                                     ERCC_counts = round(exp.spike_counts), 
                                                     log=exp.log_file,
                                                     comparison='ALL',
                                                     plot_dir = exp.scratch + 'PCA/',
                                                     de=False
                                                    )
        exp.de_results['all_ERCC_rlog'].to_csv('{}ALL-samples-blind-ERCC-rlog-counts.txt'.format(out_dir), 
                                               header=True, 
                                               index=True, 
                                               sep="\t"
                                               )
        exp.de_results['all_ERCC_normCounts'].to_csv('{}ALL-samples-ERCC_DE2_normCounts.txt'.format(out_dir), 
                                                     header=True, 
                                                     index=True, 
                                                     sep="\t"
                                                    )

    print(session(), file=open(exp.log_file, 'a'))    
    exp.tasks_complete.append('DESeq2')
    print('DESeq2 differential expression complete: {:%Y-%m-%d %H:%M:%S}\n'.format(datetime.now()), file=open(exp.log_file, 'a'))
    
    return exp

def plot_exp(data, plot_dir,exp_type, name, log_file, backend='Agg'):
    '''
    Inputs
    ------
    data: dataframe.  samples in columns, counts in rows
    plot_dir: output directory string
    exp_type: string of data type for ylabel (ex 'Normalized Log$_2$')
    name : title for plot and for file
    log_file: ouput log file
    backend: matplotlib backedn to use.  Agg default if not interactive, else use 'Qt5Agg'

    Ouputs
    ------
    None
    Prints boxplot to plot_dir

    '''
    print('Starting global sample expression comparisons.', file=open(log_file, 'a'))

    if __name__ == "__main__":
        plt.clf()
    sns.set(context='paper', font='Arial', style='white',rc={'figure.dpi':300,'figure.figsize':(4,4)})
    pl=sns.boxplot(data=data, color='darkgrey', medianprops={'color':'red'})
    pl.set_ylabel('Expression Counts\n({})'.format(exp_type))
    pl.set_title(name)
    for tick in pl.xaxis.get_ticklabels():
            tick.set_rotation(90)
    plt.tight_layout()
    sns.despine()
    plt.savefig('{}{}_expression_barplot.png'.format(plot_dir, name), dpi=300)
    if __name__ == "__main__":
        plt.close()

def Principal_Component_Analysis(exp):

    out_dir = exp.scratch + 'PCA/'
    os.makedirs(out_dir, exist_ok=True)
    
    #PCA on raw data
    print('Starting PCA analysis for all raw counts.', file=open(exp.log_file, 'a'))
    plot_PCA(counts=exp.count_matrix,
             colData=[],
             out_dir=out_dir,
             name='all_raw_counts'
             )
    plot_exp(data=exp.count_matrix,
             plot_dir=out_dir,
             exp_type= 'raw counts',
             name='all_raw_counts',
             log_file=exp.log_file
             )

    print('Starting PCA analysis for DESeq2 regularized log counts of all samples.', file=open(exp.log_file, 'a'))
    plot_PCA(counts=exp.de_results['all_rlog'],
             colData=[],
             out_dir=out_dir,
             name='all_samples_rlog'
             )

    plot_exp(data=exp.de_results['all_rlog'],
             plot_dir=out_dir,
             exp_type='Normalized log$_2$',
             name='all_samples_rlog',
             log_file=exp.log_file
             )

    if exp.norm.lower() == 'ercc':
        print('starting PCA analysis for ALL ERCC-normalized counts.', file=open(exp.log_file, 'a'))
        plot_PCA(counts = exp.de_results['all_ERCC_normCounts'],
                 colData=[],
                 out_dir=out_dir,
                 name='all_ercc_normalized_normCounts'
                )
        plot_exp(data=exp.de_results['all_ERCC_normCounts'],
                 plot_dir=out_dir,
                 exp_type='Normalized',
                 name='all_ercc_normalized_normCounts',
                 log_file=exp.log_file
                 )


        print('starting PCA analysis for all ERCC regularlized log transformed counts.', file=open(exp.log_file,'a'))
        plot_PCA(counts = exp.de_results['all_ERCC_rlog'],
                 colData=[],
                 out_dir=out_dir,
                 name='all_ercc_rlog_+counts'
                )
        plot_exp(data=exp.de_results['all_ERCC_rlog'],
                 plot_dir=out_dir,
                 exp_type='Normalized log$_2$',
                 name='all_ercc_rlog_+counts',
                 log_file=exp.log_file
                 )

    for comparison,design in exp.designs.items():
        print('Starting DESeq2 PCA analysis for {}: {:%Y-%m-%d %H:%M:%S}\n'.format(comparison, datetime.now()), file=open(exp.log_file, 'a'))
        plot_PCA(counts=exp.de_results['{}_rlog_counts'.format(comparison)],
                 colData= design['colData'],
                 out_dir=out_dir,
                 name=comparison
                )
        plot_exp(data=exp.de_results['{}_rlog_counts'.format(comparison)],
                 plot_dir=out_dir,
                 exp_type='Normalized log$_2$',
                 name=comparison,
                 log_file=exp.log_file
                 )

    if exp.gc_norm:
        print('starting PCA analysis for gc normalized raw counts.', file=open(exp.log_file, 'a'))
        plot_PCA(counts = exp.gc_count_matrix,
                 colData=[],
                 out_dir=out_dir,
                 name='gc_nromalized_raw_counts'
                )
        plot_exp(data=exp.gc_count_matrix,
                 plot_dir=out_dir,
                 exp_type='GC Normalized Raw',
                 name='gc_nromalized_raw_counts',
                 log_file=exp.log_file
                 )

    exp.tasks_complete.append('PCA')
    print('PCA for DESeq2 groups complete: {:%Y-%m-%d %H:%M:%S}\n'.format(datetime.now()), file=open(exp.log_file, 'a'))

    return exp

def Sleuth(exp):
    '''
    Differential expression using sleuth from the Pachter lab: https://pachterlab.github.io/sleuth/
    '''

    pandas2ri.activate()

    #ri.set_writeconsole_regular(rout_write(folder=exp.job_folder))
    #ri.set_writeconsole_warnerror(rout_write(folder=exp.job_folder))
    
    sleuth = importr('sleuth') 
    biomart = importr('biomaRt')
    dplyr = importr('dplyr', on_conflict="warn")
    session=r("sessionInfo")
    out_dir= exp.scratch + 'Sleuth_results/'
    kal_dir= exp.scratch + 'Kallisto_results/'
    os.makedirs(out_dir, exist_ok=True)

    for comparison,design in exp.designs.items():
        print('Beginning Sleuth differential expression analysis for {}: {:%Y-%m-%d %H:%M:%S}\n'.format(comparison, datetime.now()), file=open(exp.log_file, 'a'))

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
            raise RaiseError('Error in sleuth, pipeline only handles hg38,hg19, and mm10')

        t2g = biomart.getBM(attributes = ro.StrVector(("ensembl_transcript_id_version", "ensembl_gene_id","external_gene_name")), mart=mart)
        t2g = dplyr.rename(t2g, target_id = 'ensembl_transcript_id_version', ens_gene = 'ensembl_gene_id', ext_gene = 'external_gene_name')

        so = sleuth.sleuth_prep(s2c, target_mapping = t2g, num_cores=1, aggregation_column = 'ens_gene')
        so = sleuth.sleuth_fit(so, condition, 'full')
        so = sleuth.sleuth_fit(so, reduced, 'reduced')
        so = sleuth.sleuth_lrt(so, 'reduced', 'full')
        print(sleuth.models(so), file=open(exp.log_file,'a'))
        sleuth_table=sleuth.sleuth_results(so, 'reduced:full','lrt',show_all=True)
        exp.de_results['SL_{}'.format(comparison)] = pandas2ri.ri2py(sleuth_table)
        exp.de_results['SL_{}'.format(comparison)].to_csv('{}{}_slueth_results.txt'.format(out_dir,comparison), header=True, index=True, sep="\t")

        print(session(), file=open(exp.log_file, 'a'))    
        print('Sleuth differential expression complete: {:%Y-%m-%d %H:%M:%S}\n'.format(datetime.now()), file=open(exp.log_file, 'a'))
    

    exp.tasks_complete.append('Sleuth')
    return exp

def volcano(results, sig_up, sig_down, name, out_dir, backend='Agg'):
    '''
    Generate volcano plot from deseq2 results dataframe and significant genes

    Inputs
    ------
    results: deseq2 results as a dataframe
    sig_up: set or list of genes up to be highlithed
    sig_down: set or list of genes down to be highlighted
    name: string of name of plot
    out_dir: string of output directory
    backend: matplotlib backend to use.  Default 'Agg' if not interactive, else switch to 'Qt5Agg'

    Outputs
    -------
    None

    Saves plot to file in out_dir

    '''
    plt.switch_backend(backend)


    sns.set(context='paper', style='white', font_scale=1)
    fig = plt.figure(figsize=(6,6), dpi=300)
    ax = fig.add_subplot(111)

    results['logp'] = results.pvalue.apply(lambda x: -np.log10(x))

    scatter = ax.scatter(results.log2FoldChange, results.logp, marker='o', color='gray', alpha=0.1, s=10, label='_nolegend_')

    scatter = ax.scatter(results[results.gene_name.apply(lambda x: x in sig_up)].log2FoldChange, 
                         results[results.gene_name.apply(lambda x: x in sig_up)].logp,
                         marker = 'o', alpha = 0.3, color='firebrick', s=10, label= 'Genes UP'
                         )

    scatter = ax.scatter(results[results.gene_name.apply(lambda x: x in sig_down)].log2FoldChange,
                         results[results.gene_name.apply(lambda x: x in sig_down)].logp,
                         marker='o', alpha = 0.3, color='steelblue', s=10, label = 'Genes DOWN'
                        )

    ax.axes.set_xlabel('Fold Change (log$_2$)')
    ax.axes.set_ylabel('p-value (-log$_10$)')

    ax.legend(loc = 'upper left', markerscale=3)
    fig.suptitle(name)

    sns.despine()
    plt.tight_layout()
    plt.savefig('{}/{}-Volcano-Plot.png'.format(out_dir,name), dpi=200)
    plt.savefig('{}/{}-Volcano-Plot.svg'.format(out_dir,name), dpi=200)
    if __name__ == "__main__":
        plt.close()

    return

def sigs(exp):
    '''
    Identifies significantly differentially expressed genes at 2 fold and 1.5 fold cutoffs with q<0.05. Generates Volcano Plots of results.
    '''
    out_dir = exp.scratch + 'Sigs_and_volcano_plots/'
    os.makedirs(out_dir, exist_ok=True)

    for comparison,design in exp.designs.items():

        if exp.de_sig_overlap:
            print('Performing overlaps of signifcant genes from Kallisto/Sleuth and STAR/RSEM/DESeq2 for {comparison}.'.format(comparison=comparison), file=open(exp.log_file,'a'))
            
            exp.sig_lists[comparison] = {}
            DE_results=exp.de_results['DE2_{}'.format(comparison)]
            SL_results=exp.de_results['SL_{}'.format(comparison)]
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
            print('Only using significant genes called from DESeq2 for {} analyses.'.format(comparison), file=open(exp.log_file, 'a'))
        
            DE_results=exp.de_results['DE2_{}'.format(comparison)]

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

        #volcano_plot    
        volcano_out = '{}{}/'.format(out_dir,comparison)
        os.makedirs(volcano_out, exist_ok=True)

        print('Generating Volcano Plots using DESeq2 results for significance', file=open(exp.log_file, 'a'))
        volcano(results = DE_results, sig_up=DE2_2UP, sig_down=DE2_2DN, name='{}_2_FC'.format(comparison), out_dir=volcano_out)
        volcano(results = DE_results, sig_up=DE2_15UP, sig_down=DE2_15DN, name='{}_1.5_FC'.format(comparison), out_dir=volcano_out)
        volcano(results = DE_results, sig_up=DE2_UP, sig_down=DE2_DN, name='{}_noFC_filter'.format(comparison), out_dir=volcano_out)

    for comparison, sigs in exp.sig_lists.items():
        sig_out='{}{}/'.format(out_dir,comparison)
        os.makedirs(sig_out, exist_ok=True)
        for sig, genes in sigs.items():
            with open('{}{}.txt'.format(sig_out,sig), 'w') as file:
                for gene in genes:
                    file.write('{}\n'.format(gene))

    print('Signature and Volcano Plot generation complete: {:%Y-%m-%d %H:%M:%S}\n'.format(datetime.now()), file=open(exp.log_file, 'a'))
    exp.tasks_complete.append('Sigs')
    return exp

def clustermap(exp):
    '''
    Generate heatmap of differentially expressed genes using regularized log2counts.
    '''
    
    out_dir=exp.scratch + 'Heatmaps/'
    os.makedirs(out_dir, exist_ok=True)
    
    for comparison,design in exp.designs.items():
        rlog = exp.de_results['{}_rlog_counts'.format(comparison)]
        rlog['gene_name']=[name.split("_")[1] for name in rlog.index.tolist()]

        sig = set(exp.sig_lists[comparison]['2FC_UP'] | exp.sig_lists[comparison]['2FC_DN'])
        if len(sig) < 2:
            print('There are not enough significantly differentially expressed genes with 2 fold chagnes in {comparison}.  Ignoring heatmap for this group. \n'.format(comparison=comparison), file=open(exp.log_file,'a'))
        else:
            CM = sns.clustermap(rlog[rlog.gene_name.apply(lambda x: x in sig)].drop('gene_name',axis=1), z_score=0, method='complete', cmap='RdBu_r', yticklabels=False)
            CM.savefig('{}{}_2FC_Heatmap.png'.format(out_dir,comparison), dpi=300)
            CM.savefig('{}{}_2FC_Heatmap.svg'.format(out_dir,comparison), dpi=300)

        sig15 = set(exp.sig_lists[comparison]['15FC_UP'] | exp.sig_lists[comparison]['15FC_DN'])
        if len(sig15) < 2:
            print('There are not enough significantly differentially expressed genes with 1.5 fold chagnes in {comparison}.  Ignoring heatmap for this group. \n'.format(comparison=comparison), file=open(exp.log_file,'a'))
        else:
            CM15 = sns.clustermap(rlog[rlog.gene_name.apply(lambda x: x in sig15)].drop('gene_name',axis=1), z_score=0, method='complete', cmap='RdBu_r', yticklabels=False)
            CM15.savefig('{}{}_1.5FC_Heatmap.png'.format(out_dir,comparison), dpi=300)
            CM15.savefig('{}{}_1.5FC_Heatmap.svg'.format(out_dir,comparison), dpi=300)

        sigAll = set(exp.sig_lists[comparison]['All_UP'] | exp.sig_lists[comparison]['All_DN'])
        if len(sigAll) < 2:
            print('There are not enough significantly differentially expressed genes without a fold change in {comparison}.  Ignoring heatmap for this group. \n'.format(comparison=comparison), file=open(exp.log_file,'a'))
        else:
            CM15 = sns.clustermap(rlog[rlog.gene_name.apply(lambda x: x in sigAll)].drop('gene_name',axis=1), z_score=0, method='complete', cmap='RdBu_r', yticklabels=False)
            CM15.savefig('{}{}_noFCfilter_Heatmap.png'.format(out_dir,comparison), dpi=300)
            CM15.savefig('{}{}_noFCfilter_Heatmap.svg'.format(out_dir,comparison), dpi=300)
    
    exp.tasks_complete.append('Heatmaps')
    print('Heatmaps for DESeq2 differentially expressed genes complete: {:%Y-%m-%d %H:%M:%S}\n'.format(datetime.now()), file=open(exp.log_file, 'a'))
    
    return exp

def enrichr(gene_list, description, out_dir):
    '''
    Perform GO enrichment and KEGG enrichment Analysis using Enrichr: http://amp.pharm.mssm.edu/Enrichr/
    '''
    
    gseapy.enrichr(gene_list=gene_list,
                   description='{description}_KEGG'.format(description=description),
                   gene_sets='KEGG_2016', 
                   outdir=out_dir
                   )
    gseapy.enrichr(gene_list=gene_list,
                   description='{description}_GO_biological_process'.format(description=description),
                   gene_sets='GO_Biological_Process_2017b', 
                   outdir=out_dir
                  )
    gseapy.enrichr(gene_list=gene_list, 
                   description='{description}_GO_molecular_function'.format(description=description),
                   gene_sets='GO_Molecular_Function_2017b', 
                   outdir=out_dir
                  )
    return

def GO_enrich(exp):
    '''
    Perform GO enrichment analysis on significanttly differentially expressed genes.
    '''
    GO_dir='{}GO_enrichment/'.format(exp.scratch)
    os.makedirs(GO_dir, exist_ok=True)
    
    for comparison,design in exp.designs.items():
        print('Beginning GO enrichment for {}: {:%Y-%m-%d %H:%M:%S}\n'.format(comparison, datetime.now()), file=open(exp.log_file, 'a'))
        
        for name,sig in exp.sig_lists[comparison].items():
            if len(sig) == 0:
                print('There are no significantly differentially expressed genes in {} {}.  Ignoring GO enrichment. \n'.format(name,comparison), file=open(exp.log_file,'a'))
            else:
                GO_out = '{}{}/'.format(GO_dir, comparison)
                os.makedirs(GO_out,exist_ok=True)
                enrichr(gene_list=list(sig), description='{}_{}'.format(comparison,name),out_dir=GO_out)

    exp.tasks_complete.append('GO_enrich')
    print('GO Enrichment analysis for DESeq2 differentially expressed genes complete: {:%Y-%m-%d %H:%M:%S}\n'.format(datetime.now()), file=open(exp.log_file, 'a'))
    
    return exp

def gsea_barplot(out_dir,pos_file,neg_file,gmt_name,max_number=20):
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

    '''

    out_dir = out_dir if out_dir.endswith('/') else '{}/'.format(out_dir)
    out_dir = '' if out_dir == '/' else out_dir 
    pos = pd.read_table(pos_file).head(max_number) if os.path.isfile(pos_file) else pd.DataFrame(columns=['FDR q-val'])
    pos[gmt_name] = [' '.join(name.split('_')[1:]) for name in pos.NAME.tolist()]
    neg = pd.read_table(neg_file).head(max_number) if os.path.isfile(neg_file) else pd.DataFrame(columns=['FDR q-val'])
    neg[gmt_name] = [' '.join(name.split('_')[1:]) for name in neg.NAME.tolist()]
    
    sns.set(context='paper', font='Arial',font_scale=.9, style='white', rc={'figure.dpi': 300, 'figure.figsize':(8,6)})
    fig,(ax1,ax2) = plt.subplots(ncols=1, nrows=2)
    fig.suptitle('{} GSEA enrichment\n(q<0.05, max {})'.format(gmt_name, max_number))
    
    if len(pos[pos['FDR q-val'] < 0.05]) > 0:
        UP = sns.barplot(data=pos[pos['FDR q-val'] < 0.05], x = 'NES', y=gmt_name, color='firebrick', ax=ax1)
        UP.set_title('Positive Enrichment')
        sns.despine()
       
    if len(neg[neg['FDR q-val'] < 0.05]) > 0:
        DN = sns.barplot(data=neg[neg['FDR q-val'] < 0.05], x = 'NES', y=gmt_name, color='steelblue', ax=ax2)
        DN.set_title('Negative Enrichment')
        sns.despine()
    
    try:
        plt.tight_layout(h_pad=1,w_pad=1)
    except ValueError:
        pass

    plt.subplots_adjust(top=0.88)
    file='{}{}_GSEA_NES_plot.png'.format(out_dir,gmt_name)
    fig.savefig(file, dpi=300)
    plt.close()

    return file

def GSEA(exp):
    '''
    Perform Gene Set Enrichment Analysis using gsea 3.0 from the Broad Institute.
    '''

    print('Starting GSEA enrichment.', file=open(exp.log_file, 'a'))

    out_dir = '{}DESeq2_GSEA'.format(exp.scratch)
    os.makedirs(out_dir, exist_ok=True)

    if exp.genome == 'mm10':
        gmt_list = glob.glob('{}*.gmt'.format(exp.genome_indicies['GMT']))
        gmts={'Hallmarks': [gmt for gmt in gmt_list if 'h.all' in gmt][0],
              'KEGG':[gmt for gmt in gmt_list if 'c2.cp.kegg' in gmt][0],
              'GO_Biological_Process':[gmt for gmt in gmt_list if 'c5.bp' in gmt][0],
              'GO_Molecular_Function':[gmt for gmt in gmt_list if 'c5.mf' in gmt][0],
              'Curated_Gene_Sets':[gmt for gmt in gmt_list if 'c2.cgp' in gmt][0]
              }
    else:
        gmts={'Hallmarks': 'gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/h.all.v6.2.symbols.gmt',
              'KEGG':'gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c2.cp.kegg.v6.2.symbols.gmt',
              'GO_Biological_Process':'gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c5.bp.v6.2.symbols.gmt',
              'GO_Molecular_Function':'gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c5.mf.v6.2.symbols.gmt',
              'Curated_Gene_Sets':'gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c2.cgp.v6.2.symbols.gmt'
              }

    for comparison,design in exp.designs.items():
        
        indexpath=glob.glob('{}/{}/Hallmarks/*/index.html'.format(out_dir,comparison))
        if len(indexpath) > 0:
            print('GSEA Hallmarks already done for {name} for {comparison}. Skipping GSEA for {comparison}\n'.format(name=name,comparison=comparison), file=open(exp.log_file,'a'))
        else:
            print('Running GSEA for {comparison}. Results found in {out}DESeq2_GSEA/{comparison}. \n'.format(comparison=comparison, out=exp.out_dir), file=open(exp.log_file, 'a'))
            out_compare = '{}/{}'.format(out_dir,comparison)
            os.makedirs(out_compare, exist_ok=True)

            results=exp.de_results['DE2_{}'.format(comparison)].dropna()

            #generate ranked list based on shrunken log2foldchange
            lfc = exp.de_results['shrunkenLFC_' + comparison].dropna()
            lfc.sort_values(by='log2FoldChange', ascending=False, inplace=True)
            lfc.index = lfc.gene_name
            lfc = lfc.log2FoldChange.dropna()
            lfc.to_csv('{}/{}_shrunkenLFC.rnk'.format(out_compare,comparison), header=False, index=True, sep="\t")

            print('Using Wald statistic for gene preranking.', file = open(exp.log_file,'a'))
            rnk = '{}_stat.rnk'.format(comparison)

            if exp.genome == 'mm10':
                results['Ens_ID'] = [ID.split('.')[0] for ID in results.index.tolist()]
                results[['Ens_ID','stat']].dropna(subset=['stat']).to_csv('{}/{}_stat.rnk'.format(out_compare,comparison), header=None, index=None, sep="\t")
            else:
                results[['gene_name','stat']].dropna(subset=['stat']).to_csv('{}/{}_stat.rnk'.format(out_compare,comparison), header=None, index=None, sep="\t")

            os.chdir(out_compare)

            print('Beginning GSEA enrichment for {} using preranked genes: {:%Y-%m-%d %H:%M:%S}'.format(comparison, datetime.now()), file=open(exp.log_file, 'a'))
            print('Genes with positive LFC (to the left left in GSEA output graph) are upregulated in experimental vs control conditions. Genes with negative LFC (on right) are downregulated genes in experimental samples vs controls.\n', file=open(exp.log_file,'a'))
                  
            for name,gset in gmts.items():
                set_dir='{}/{}'.format(out_compare,name) 
                os.makedirs(set_dir, exist_ok=True)

                command_list = ['module rm python java perl share-rpms65',
                                'source activate RNAseq',
                                'java -cp {jar} -Xmx2048m xtools.gsea.GseaPreranked -gmx {gset} -norm meandiv -nperm 1000 -rnk {rnk} -scoring_scheme weighted -rpt_label {comparison}_{name}_wald -create_svgs false -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 1000 -set_min 10 -zip_report false -out {name} -gui false'.format(jar=exp.genome_indicies['GSEA_jar'],gset=gset,comparison=comparison,name=name,rnk=rnk)
                               ] #for shrunken lfc: 'java -cp {jar} -Xmx2048m xtools.gsea.GseaPreranked -gmx gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/{gset}.v6.1.symbols.gmt -norm meandiv -nperm 1000 -rnk {rnk2} -scoring_scheme weighted -rpt_label {comparison}_{gset}_shrunkenLFC -create_svgs false -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 1000 -set_min 10 -zip_report false -out {name} -gui false'.format(jar=exp.genome_indicies['GSEA_jar'],gset=gset,comparison=comparison,name=name,rnk2=rnk2)

                exp.job_id.append(send_job(command_list=command_list, 
                                           job_name='{}_{}_GSEA'.format(comparison,name),
                                           job_log_folder=exp.job_folder,
                                           q= 'general',
                                           mem=3000,
                                           log_file=exp.log_file,
                                           project=exp.project
                                          )
                                 )
                time.sleep(1)

    #Wait for jobs to finish
    job_wait(id_list=exp.job_id, job_log_folder=exp.job_folder, log_file=exp.log_file)

    for comparison,design in exp.designs.items():
        for name,gset in gmts.items():
            path=glob.glob('{}/{}/{}/*'.format(out_dir,comparison,name))[0]
            if 'index.html' == '{}/index.html'.format(path).split('/')[-1]:
                new_dir='{}/{}/{}'.format(out_dir,comparison,name)
                os.chdir(new_dir)
                
                pos_file = glob.glob('{}/gsea_report_for_na_pos*.xls'.format(path))[0]
                neg_file = glob.glob('{}/gsea_report_for_na_neg*.xls'.format(path))[0]

                barplot = gsea_barplot(out_dir = new_dir, pos_file=pos_file, neg_file=neg_file, gmt_name=name)

                msg = 'Open "index.html" in subfolder for all results.\nOpen {} for barplot summary of top enrichments.'.format(barplot)

                with open('README.txt','w') as fp:
                    fp.write(msg)

                print('{} GSEA enrichment barplot for {} can be found here:\n {}'.format(comparison,name,barplot), open(exp.log_file,'a'))

            else:
                print('GSEA did not complete {} for {}.'.format(name,comparison), file=open(exp.log_file,'a'))            

    os.chdir(exp.scratch)
    exp.tasks_complete.append('GSEA')
    print('GSEA analysis complete: {:%Y-%m-%d %H:%M:%S}\n'.format(datetime.now()), file=open(exp.log_file, 'a'))
    
    return exp

def plot_venn2(Series, string_name_of_overlap, folder, backend='Agg'):
    '''
    Series with with overlaps 10,01,11
    Plots a 2 way venn.
    Saves to file.
    '''
    plt.switch_backend(backend)

    os.makedirs(folder, exist_ok=True)

    plt.figure(figsize=(7,7))
    
    font = {'family': 'sans-serif',
            'weight': 'normal',
            'size': 16,
           }
    
    plt.rc('font', **font)
    
    #make venn
    venn_plot = venn2(subsets=(Series.iloc[0], Series.iloc[1], Series.iloc[2]), set_labels = Series.index.tolist())
    patch=['10','01','11']
    colors=['green','blue','teal']
    for patch,color in zip(patch,colors):
        venn_plot.get_patch_by_id(patch).set_color('none')
        venn_plot.get_patch_by_id(patch).set_alpha(0.4)
        venn_plot.get_patch_by_id(patch).set_edgecolor('none')   

    c= venn2_circles(subsets=(Series.iloc[0], Series.iloc[1], Series.iloc[2]))
    colors_circle=['green','blue']
    for circle,color in zip(c,colors_circle): 
        circle.set_edgecolor(color)
        circle.set_alpha(0.8)
        circle.set_linewidth(3)
        
    plt.title(string_name_of_overlap + " Overlaps")
    plt.tight_layout()
    plt.savefig('{}{}-overlap-{:%Y-%m-%d}.svg'.format(folder,string_name_of_overlap, datetime.now()))
    plt.savefig('{}{}-overlap-{:%Y-%m-%d}.png'.format(folder,string_name_of_overlap, datetime.now()), dpi=300)
    if __name__ == "__main__":
        plt.close()

def overlaps(exp):
    '''
    Performs overlaps of two or more de_sig lists.
    '''
    out_dir = exp.scratch + 'Overlaps/'
    os.makedirs(out_dir, exist_ok=True)
    
    if len(exp.overlaps) != 0:
        names=['2FC_UP', '2FC_DN', '15FC_UP','15FC_DN', 'All_UP', 'All_DN']
        print('Beginning overlap of significant genes: {:%Y-%m-%d %H:%M:%S}\n'.format(datetime.now()), file=open(exp.log_file, 'a'))

        for overlap,comparison_list in exp.overlaps.items():
            if len(comparison_list) != 0:
                for name in names:
                    key= '{}_{}'.format(overlap,name)
                    exp.overlap_results[key] = exp.sig_lists[comparison_list[0]][name] & exp.sig_lists[comparison_list[1]][name] 
                    
                    if len(exp.overlap_results[key]) == 0:
                        print('{}_{} have no overlapping genes'.format(overlap,name), file=open(exp.log_file,'a'))
                    else:
                        venn = pd.Series([len(exp.sig_lists[comparison_list[0]][name])-len(exp.overlap_results[key]),
                                          len(exp.sig_lists[comparison_list[1]][name])-len(exp.overlap_results[key]),
                                          len(exp.overlap_results[key])
                                         ],
                                         index= comparison_list + ['Overlap']
                                        )
                        plot_venn2(venn, key, '{}{}/'.format(out_dir,name))
            
    elif len(exp.gene_lists) != 0:
        for name, gene_list in exp.gene_lists.items():
            exp.overlap_results[name]= gene_list[0] & gene_list[1]
            if len(exp.overlap_results[name]) == 0:
                    print('{} has no overlapping genes'.format(name), file=open(exp.log_file,'a'))
            else:
                list_names = gene_list.keys()
                venn = pd.Series([len(gene_list[list_names[0]])-len(exp.overlap_results[name]),
                                  len(gene_list[list_names[1]])-len(exp.overlap_results[name]),
                                  len(overlap_results[name])
                                 ],
                                 index= list_names + ['Overlap']
                                )
                plot_venn2(venn,name,'{}{}/'.format(out_dir,name))

    for name,sig in exp.overlap_results.items():
        sig_out='{}{}/'.format(out_dir,name)
        os.makedirs(sig_out, exist_ok=True)

        if len(sig) == 0:
            print('Not performing GO enrichment for {} overlaps since there are no overlapping genes.\n'.format(name), file=open(exp.log_file, 'a'))
        else:
            print('Performing GO enrichment for {} overlaps: {} \n'.format(name,datetime.now()), file=open(exp.log_file, 'a'))                    
            enrichr(gene_list=list(sig), description='{}_overlap'.format(name),out_dir=sig_out)

            with open('{}{}.txt'.format(sig_out,name), 'w') as file:
                for gene in list(sig):
                    file.write('{}\n'.format(gene))

    exp.tasks_complete.append('Overlaps')
    print('Overlap analysis complete: {:%Y-%m-%d %H:%M:%S}\n'.format(datetime.now()), file=open(exp.log_file, 'a'))
                   
    return exp

def final_qc(exp):
    try:
        print('Beginning final qc: {:%Y-%m-%d %H:%M:%S}\n'.format(datetime.now()), file=open(exp.log_file, 'a'))
        
        os.chdir(exp.scratch)

        command_list = ['module rm python',
                        'source activate RNAseq',
                        'multiqc *'
                       ]
        
        exp.job_id.append(send_job(command_list=command_list, 
                                   job_name= 'MultiQC',
                                   job_log_folder=exp.job_folder,
                                   q= 'general',
                                   mem=1000,
                                   log_file=exp.log_file,
                                   project=exp.project
                                  )
                         )
        
        #Wait for jobs to finish
        job_wait(id_list=exp.job_id, job_log_folder=exp.job_folder, log_file=exp.log_file)
        
        if os.path.isdir(exp.scratch + '/multiqc_data'):
            rmtree(exp.scratch + '/multiqc_data')

        exp.tasks_complete.append('MultiQC')
        
        return exp

    except:
        print('Error during MultiQC.', file=open(exp.log_file,'a'))
        filename= '{}{}_incomplete.pkl'.format(exp.scratch,exp.name)
        with open(filename, 'wb') as experiment:
            pickle.dump(exp, experiment)
        raise RaiseError('Error during MultiQC. Fix problem then resubmit with same command to continue from last completed step.')

def finish(exp):
    try:

        os.chdir(exp.scratch)
        
        if exp.seq_type == 'paired':
            for number,sample in exp.samples.items():
                R_list = ['{}{}_R1.fastq.gz'.format(exp.fastq_folder,sample),
                          '{}{}_R2.fastq.gz'.format(exp.fastq_folder,sample),
                          '{}{}_trim_R1.fastq.gz'.format(exp.fastq_folder,sample),
                          '{}{}_trim_R2.fastq.gz'.format(exp.fastq_folder,sample)
                         ]
                for R in R_list:
                    if os.path.isfile(R):
                        os.remove(R)
        elif exp.seq_type == 'single':
            for number,sample in exp.samples.items():
                R_list = ['{}{}.fastq.gz'.format(exp.fastq_folder,sample),
                          '{}{}_trim.fastq.gz'.format(exp.fastq_folder,sample)                          
                         ]
                for R in R_list:
                    if os.path.isfile(R):
                        os.remove(R)
        
        print('\nPackage versions: ', file=open(exp.log_file, 'a'))
        with open('/projects/ctsi/nimerlab/DANIEL/tools/nimerlab-pipelines/RNAseq/environment.yml','r') as file:
            versions = yaml.load(file)
        for package in versions['dependencies']:
            print(package, file=open(exp.log_file, 'a'))

        print('\n{} analysis complete!  Performed the following tasks: '.format(exp.name)+ '\n', file=open(exp.log_file, 'a'))
        print(str(exp.tasks_complete) + '\n', file=open(exp.log_file, 'a'))
        
        scratch_log= exp.scratch + exp.log_file.split("/")[-1]
        copy2(exp.log_file, scratch_log)
        rmtree(exp.out_dir)
        copytree(exp.scratch, exp.out_dir)

        exp.tasks_complete.append('Finished')

        filename= '{}{}_{}.pkl'.format(exp.out_dir,exp.name,exp.date)
        with open(filename, 'wb') as experiment:
            pickle.dump(exp, experiment) 

        print('Moved all files into {}: {:%Y-%m-%d %H:%M:%S}\n'.format(exp.out_dir, datetime.now()), file=open(exp.log_file, 'a'))
        print("\n Finger's Crossed!!!", file=open(exp.log_file, 'a'))

    except:
        print('Error while finishing pipeline.', file=open(exp.log_file,'a'))
        filename= '{}{}_incomplete.pkl'.format(exp.scratch, exp.name)
        with open(filename, 'wb') as experiment:
            pickle.dump(exp, experiment)
        raise RaiseError('Error finishing pipeline. Fix problem then resubmit with same command to continue from last completed step.')

def preprocess(exp):
    try:
        pipe_stage='preprocessing'
        #exp=fastq_cat(exp)
        if 'Stage' not in exp.tasks_complete:
            pipe_stage = 'staging'
            exp=stage(exp) 
        if 'Fastq_screen' not in exp.tasks_complete:
            pipe_stage = 'contamination screening'
            exp=fastq_screen(exp)
        if 'Trim' not in exp.tasks_complete:
            pipe_stage = 'fastq trimming'
            exp=trim(exp)
        if 'FastQC' not in exp.tasks_complete:
            pipe_stage = 'FastQC'
            exp=fastqc(exp)  
        return exp
    except:
        print('Error in {}.'.format(pipe_stage), file=open(exp.log_file,'a'))
        filename= '{}{}_incomplete.pkl'.format(exp.scratch,exp.name)
        with open(filename, 'wb') as experiment:
            pickle.dump(exp, experiment)
        raise RaiseError('Error in {}. Fix problem then resubmit with same command to continue from last completed step.'.format(pipe_stage))

def align(exp):
    try:
        pipe_stage='alignment'
        if 'Spike' not in exp.tasks_complete:
            pipe_stage = 'spike in processing'
            exp=spike(exp)
        if 'STAR' not in exp.tasks_complete:
            pipe_stage = 'STAR alignment'
            if exp.alignment_mode == 'gene':
                exp = star(exp)
            else:
                exp=rsem(exp)
        if 'Kallisto' not in exp.tasks_complete:
            pipe_stage = 'Kallisto alignment'
            exp=kallisto(exp)
        return exp
    except:
        print('Error in {}.'.format(pipe_stage), file=open(exp.log_file,'a'))
        filename= '{}{}_incomplete.pkl'.format(exp.scratch,exp.name)
        with open(filename, 'wb') as experiment:
            pickle.dump(exp, experiment)
        raise RaiseError('Error in {}. Fix problem then resubmit with same command to continue from last completed step.'.format(pipe_stage))

def diff_exp(exp):
    try:
        pipe_stage='differential expression'
        if 'GC' not in exp.tasks_complete:
            pipe_stage = 'GC Normalization'
            exp = GC_normalization(exp)
        if 'DESeq2' not in exp.tasks_complete:
            pipe_stage = 'DESeq2'
            exp = DESeq2(exp)
        if 'PCA' not in exp.tasks_complete:
            pipe_stage = 'PCA'
            exp = Principal_Component_Analysis(exp)
        if 'Sleuth' not in exp.tasks_complete:
            pipe_stage = 'Sleuth'
            exp = Sleuth(exp)
        if 'Sigs' not in exp.tasks_complete:
            pipe_stage = 'signature generation'
            exp = sigs(exp)
        if 'Heatmaps' not in exp.tasks_complete:
            pipe_stage = 'heatmap generation'
            sep = clustermap(exp)
        if 'GO_enrich' not in exp.tasks_complete:
            pipe_stage = 'GO enrichment'
            exp = GO_enrich(exp)
        if 'GSEA' not in exp.tasks_complete:
            pipe_stage = 'GSEA'
            exp = GSEA(exp)
        if 'Overlaps' not in exp.tasks_complete:
            pipe_stage = 'signature overlaps'
            exp = overlaps(exp)
        #exp = decomposition(exp)  
        return exp
    except:
        print('Error in {}.'.format(pipe_stage), file=open(exp.log_file,'a'))
        filename= '{}{}_incomplete.pkl'.format(exp.scratch,exp.name)
        with open(filename, 'wb') as experiment:
            pickle.dump(exp, experiment)
        raise RaiseError('Error in {}. Fix problem then resubmit with same command to continue from last completed step.'.format(pipe_stage))

def pipeline(experimental_file):
    exp = parse_yaml(experimental_file)
    exp = preprocess(exp)
    exp = align(exp)
    exp = diff_exp(exp)
    if 'MultiQC' not in exp.tasks_complete:
        exp = final_qc(exp)
    finish(exp)

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--experimental_file', '-f', required=True, help='experimental yaml file', type=str)
    args = parser.parse_args()
    
    pipeline(args.experimental_file)

