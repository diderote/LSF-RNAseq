#!/usr/bin/env python3

'''
Nimerlab RNASeq Pipeline v0.5
python3/utf-8

Reads an experimetnal design yaml file (Version 0.5).
Requires a conda environment 'RNAseq' made from environment.yml

To do:
    - ICA with chi-square with de groups
    - t-SNE (add as option)
    - if no chrname - skip bigwig generation

Author: Daniel Karl

'''

import os
import re
import glob
import pickle
import math
import random
import time
import yaml
from shutil import copy2,copytree,rmtree,move
from datetime import datetime
import subprocess as sub
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib_venn import venn2, venn2_circles
import seaborn as sns
from sklearn.decomposition import PCA
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri, r, globalenv, Formula
import gseapy

version=0.5

class Experiment(object):
    '''
    Experiment object for pipeline
    '''
    def __init__(self, scratch='', date='', name='', out_dir='', job_folder='', qc_folder='', 
                  log_file='',fastq_folder='',spike=False, count_matrix=pd.DataFrame(), trim=[0,0],
                  spike_counts=pd.DataFrame(),genome='',sample_number=int(), samples={}, 
                  job_id=[],de_groups={},norm='Median-Ratios',designs={}, overlaps={}, gene_lists={},
                  tasks_complete=[],de_results={},sig_lists={},overlap_results={},de_sig_overlap={},
                  genome_indicies={},project='', gc_norm=False, gc_count_matrix=pd.DataFrame(),
                  seq_type='', alignment_mode='transcript'
                 ):
        self.scratch = scratch
        self.date = date
        self.name = name
        self.out_dir =out_dir
        self.job_folder=job_folder
        self.qc_folder=qc_folder
        self.log_file=log_file
        self.fastq_folder=fastq_folder
        self.spike = spike
        self.count_matrix = count_matrix
        self.trim=trim
        self.spike_counts = spike_counts
        self.genome = genome
        self.sample_number =sample_number
        self.samples = samples
        self.job_id=job_id
        self.de_groups = de_groups
        self.norm = norm
        self.designs=designs
        self.overlaps = overlaps
        self.gene_lists=gene_lists
        self.tasks_complete=tasks_complete
        self.de_results = de_results
        self.sig_lists=sig_lists
        self.overlap_results=overlap_results
        self.de_sig_overlap = de_sig_overlap
        self.genome_indicies=genome_indicies
        self.project=project
        self.gc_norm=gc_norm
        self.gc_count_matrix = gc_count_matrix
        self.seq_type=seq_type
        self.alignment_mode=alignment_mode

class RaiseError(Exception):
    pass

def parse_yaml():
    '''
    Parse experimental info from yaml file
    '''    
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--experimental_file', '-f', required=True, help='experimental yaml file', type=str)
    args = parser.parse_args()
    
    with open(args.experimental_file,'r') as file:
        yml = yaml.safe_load(file)

    #Make a new experimental object
    exp = Experiment()
    
    #Setting Scratch folder
    if yml['Lab'].lower() == 'nimer':
        exp.scratch = '/scratch/projects/nimerlab/DANIEL/staging/RNAseq/{}/'.format(yml['Name'])
        os.makedirs(exp.scratch, exist_ok=True)
    else:
        try:
            exp.scratch = '{}{}/'.format(yml['Scratch_folder'],yml['Name'])
            os.makedirs(exp.scratch, exist_ok=True)
        except:
            raise RaiseError('Error making scratch/staging directory', file=open(exp.log_file,'a'))
    
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
    
    if yml['Output_directory'].endswith('/'):
        exp.out_dir = '{}{}/'.format(yml['Output_directory'],exp.name)
    else:
        exp.out_dir = '{}/{}/'.format(yml['Output_directory'], exp.name)
    
    #Make out directory if it doesn't exist
    os.makedirs(exp.out_dir, exist_ok=True)

    #Log file
    exp.log_file = '{}{}-{}.log'.format(exp.out_dir,exp.name,exp.date)
    
    print('Pipeline version {} run on {} \n'.format(str(version),exp.date), file=open(exp.log_file, 'w'))
    print('Beginning RNAseq Analysis: {:%Y-%m-%d %H:%M:%S}\n'.format(datetime.now()), file=open(exp.log_file, 'a'))
    print('Reading experimental file...\n', file=open(exp.log_file, 'a'))
    print("Pipeline output folder: {}\n".format(exp.out_dir), file=open(exp.log_file, 'a'))

    #Setting Job Folder
    exp.job_folder = '{}logs/'.format(exp.scratch)
    os.makedirs(exp.job_folder, exist_ok=True)

    #Genome
    if yml['Genome'].lower() not in ['hg38', 'mm10', 'hg19']:
        raise ValueError("Genome must be either hg38, hg19, or mm10.")
    else:
        exp.genome = yml['Genome'].lower()
        print('Processing data with: ' + str(exp.genome), file=open(exp.log_file, 'a'))

    #Sequencing type
    if yml['Sequencing_type'].lower() not in ['paired','single']:
        raise ValueError("Must specify whether sequence is paired or single end.")
    else:
        exp.seq_type = yml['Sequencing_type'].lower()
        if exp.seq_type == 'single':
            exp.tasks_complete = exp.tasks_complete + ['Kallisto','Sleuth']
        print('Processing data as {}-end sequencing.'.format(exp.seq_type), file =open(exp.log_file,'a'))

    #Tasks to complete
    if yml['Tasks']['Align'] == False:
        exp.tasks_complete = exp.tasks_complete + ['Stage','FastQC','Fastq_screen','Trim','STAR','Kallisto', 'Sleuth']
        print('Not performing alignment.', file=open(exp.log_file,'a'))
        count_matrix_loc=yml['Count_matrix']
        if os.path.exists(count_matrix_loc):
            print("Count matrix found at {}".format(count_matrix_loc), file=open(exp.log_file, 'a'))
            print("Performing only DESeq2 on for DE", file=open(exp.log_file, 'a'))
            if count_matrix_loc.split('.')[-1] == 'txt':
                exp.count_matrix = pd.read_csv(count_matrix_loc, header= 0, index_col=0, sep="\t")
            elif (count_matrix_loc.split('.')[-1] == 'xls') or (count_matrix_loc.split('.')[-1] == 'xlsx'):
                exp.count_matrix = pd.read_excel(count_matrix_loc)
            else:
                raise IOError("Cannot parse count matrix.  Make sure it is .txt, .xls, or .xlsx")
        else:
            raise IOError("Count Matrix Not Found.")
    elif yml['Tasks']['Align'] == True:
        #Alignment mode. Default is transcript.
        if yml['Tasks']['Alignment_Mode'].lower() == 'gene':
            exp.alignment_mode = 'gene'
    else:
        raise IOError('Please specify whether or not to perform alignment.')   
    
    #Lab specific files
    if yml['Lab'].lower() == 'other':
        exp.genome_indicies['RSEM_STAR'] = yml['RSEM_STAR_index']
        exp.genome_indicies['STAR'] = yml['STAR_index']
        exp.genome_indicies['Kallisto'] = yml['Kallisto_index']
        exp.genome_indicies['ERCC'] = yml['ERCC_STAR_index']
        exp.genome_indicies['chrLen'] = yml['ChrNameLength_file']
    elif yml['Lab'].lower() == 'nimer':
        exp.genome_indicies['ERCC'] = '/projects/ctsi/nimerlab/DANIEL/tools/genomes/ERCC_spike/STARIndex'
        if exp.genome == 'mm10':
            exp.genome_indicies['RSEM_STAR'] = '/projects/ctsi/nimerlab/DANIEL/tools/genomes/Mus_musculus/mm10/RSEM-STARIndex/mouse'
            exp.genome_indicies['STAR'] = '/projects/ctsi/nimerlab/DANIEL/tools/genomes/Mus_musculus/mm10/STARIndex'
            exp.genome_indicies['chrLen'] = '/projects/ctsi/nimerlab/DANIEL/tools/genomes/Mus_musculus/mm10/RSEM-STARIndex/chrNameLength.txt'
            exp.genome_indicies['Kallisto'] = '/projects/ctsi/nimerlab/DANIEL/tools/genomes/Mus_musculus/mm10/KallistoIndex/GRCm38.transcripts.idx'
        elif exp.genome == 'hg38':
            exp.genome_indicies['RSEM_STAR'] = '/projects/ctsi/nimerlab/DANIEL/tools/genomes/H_sapiens/NCBI/GRCh38/Sequence/RSEM-STARIndex/human'
            exp.genome_indicies['STAR'] = '/projects/ctsi/nimerlab/DANIEL/tools/genomes/H_sapiens/NCBI/GRCh38/Sequence/STARIndex'
            exp.genome_indicies['chrLen'] = '/projects/ctsi/nimerlab/DANIEL/tools/genomes/H_sapiens/NCBI/GRCh38/Sequence/RSEM-STARIndex/chrNameLength.txt'
            exp.genome_indicies['Kallisto'] = '/projects/ctsi/nimerlab/DANIEL/tools/genomes/H_sapiens/NCBI/GRCh38/Sequence/KallistoIndex/GRCh38.transcripts.idx'
        elif exp.genome == 'hg19':
            exp.genome_indicies['RSEM_STAR'] = '/projects/ctsi/nimerlab/DANIEL/tools/genomes/H_sapiens/Hg19/NCBI-RNAseq/RSEM-STAR/human'
            exp.genome_indicies['STAR'] = '/projects/ctsi/nimerlab/DANIEL/tools/genomes/H_sapiens/Hg19/NCBI-RNAseq/STAR'
            exp.genome_indicies['chrLen'] = '/projects/ctsi/nimerlab/DANIEL/tools/genomes/H_sapiens/Hg19/NCBI-RNAseq/RSEM-STAR/chrNameLength.txt'
            exp.genome_indicies['Kallisto'] = '/projects/ctsi/nimerlab/DANIEL/tools/genomes/H_sapiens/Hg19/NCBI-RNAseq/Kallisto/GRCh37.transcripts.idx'
 
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
    if yml['ERCC_spike']:
        exp.spike = True
        if 'Stage' in exp.tasks_complete:
            if yml['Fastq_directory'][-1] != '/':
                yml['Fastq_directory'] = yml['Fastq_directory'] + '/'
            if os.path.isdir(yml['Fastq_directory']):
                exp.fastq_folder=yml['Fastq_directory']
            else:
                raise IOError("Can't Find Fastq Folder.")

    if (yml['Normalization'].lower() == 'ercc') and (yml['Tasks']['Align'] == False):
        spike_matrix_loc = yml['Spike_matrix']
        if os.path.exists(spike_matrix_loc):
            print("Spike Count matrix found at {}".format(spike_matrix_loc), file=open(exp.log_file, 'a'))
            if spike_matrix_loc.split('.')[-1] == 'txt':
                exp.spike_counts = pd.read_csv(spike_matrix_loc, header= 0, index_col=0, sep="\t")
            elif (spike_matrix_loc.split('.')[-1] == 'xls') or (spike_matrix_loc.split('.')[-1] == 'xlsx'):
                exp.spike_counts = pd.read_excel(spike_matrix_loc)
            else:
                raise IOError("Cannot parse spike count matrix.  Make sure it is .txt, .xls, or .xlsx")
        else:
            raise IOError("Spike Count Matrix Not Found. ")

    #Fastq Folder
    if 'Stage' not in exp.tasks_complete:
        if yml['Fastq_directory'][-1] != '/':
            yml['Fastq_directory'] = yml['Fastq_directory'] +'/'
        if os.path.isdir(yml['Fastq_directory']):
            exp.fastq_folder=yml['Fastq_directory']
        else:
            raise IOError("Can't Find Fastq Folder.")
    
    #Hard clip
    if exp.trim != [0,0]:
        exp.trim = yml['trim']

    #Project
    if yml['Lab'].lower()=='nimer':
        exp.project = 'nimerlab'
    elif yml['Lab'].lower() == 'other':
        if len(yml['Pegasus_Project']) == 0:
            exp.project = ''
        else:
            exp.project = yml['Pegasus_Project']
    
    #Counts
    if not 0 < yml['Total_sample_number'] < 33:
        raise ValueError("This pipeline is only set up to handle up to 32 samples.")
    else:
        exp.sample_number = yml['Total_sample_number']
        print('Processing ' + str(exp.sample_number) + ' samples.'+ '\n', file=open(exp.log_file, 'a'))
    
    #Sample Names
    count = 1
    for key,name in yml['Samples'].items():
        if count <= exp.sample_number:
            exp.samples[key]=name
            count += 1
        else:
            break
    print("Samples: ", file=open(exp.log_file, 'a'))
    for number,sample in exp.samples.items():
        print('{number}: {sample}'.format(number=number,sample=sample), file=open(exp.log_file, 'a'))
    
    #Differential Expression Groups
    if yml['Tasks']['Differential_Expression']:
        for key, item in yml['Groups'].items():
            if bool(item):
                exp.de_groups[key] = [exp.samples[int(x)] for x in item.split(',')]
            
        print("Parsing experimental design for differential expression...\n", file=open(exp.log_file, 'a'))
        
        #Normalization method
        if yml['Normalization'].lower() == 'ercc':
            exp.norm = 'ERCC' 
            print('Normalizing samples for differential expression analysis using ERCC spike-in variance'+ '\n', file=open(exp.log_file, 'a'))
        elif yml['Normalization'].lower() == 'empirical':
            print('Normalizing samples for differential expression analysis using empirical negative controls for variance'+ '\n', file=open(exp.log_file, 'a'))
            exp.norm = 'empirical'
        elif yml['Normalization'].lower() == 'median-ratios':
            print('Normalizing samples for differential expression analysis using deseq2 size factors determined using default median of ratios method.'+ '\n', file=open(exp.log_file, 'a'))
        else:
            print("I don't know the " + yml['Normalization'] + ' normalization method.  Using default median-ratios.'+ '\n', file=open(exp.log_file, 'a'))
    
        for key, comparison in yml['Comparisons'].items():
            if bool(comparison):

                #determine comparison strategy for DESeq2
                exp_groups = comparison.split('v')[0]
                if len(exp_groups.split('-')) == 1:
                    exp_group_type = 'simple'
                elif len(exp_groups.split('-')) == 2:
                    Exp1, Exp2 = exp_groups.split('-')[0], exp_groups.split('-')[1]
                    exp_group_type = 'compensation'
                else:
                    raise ValueError("Cannot process {}. Check formating for E1 or E1-E2. Or too many groups for current pipline.".format(key))
                
                cont_groups = comparison.split('v')[1]
                if len(cont_groups.split('-')) == 1:
                    cont_group_type = 'simple'
                elif len(cont_groups.split('-'))== 2:
                    Cont1,Cont2 = cont_groups.split('-')[0], cont_groups.split('-')[1]
                    cont_group_type = 'compensation'
                else:
                    raise ValueError("Cannot process {}. Check formating for C1 or C1-C2. Or too man comparisons.".format(key))
            
                #Check comparison for group consistency.
                error = "Can't make a comparison with an unspecified group. Make sure your Comparisons match your Groups for DE"
                groups = list(exp.de_groups.keys())
                
                exp.designs[key]={}
                exp.designs[key]['all_samples']=[]
                exp.designs[key]['main_comparison']=[]
                exp.designs[key]['compensation']=[]
                
                if exp_group_type == 'simple':
                    if exp_groups not in groups:
                        raise ValueError(error)
                    else:
                        exp.designs[key]['all_samples'].extend(exp.de_groups[exp_groups])
                        exp.designs[key]['main_comparison'].extend(['Experimental']*len(exp.de_groups[exp_groups]))
                        if cont_group_type == 'simple':
                            if cont_groups not in groups:
                                raise ValueError(error)
                            else:
                                exp.designs[key]['all_samples'].extend(exp.de_groups[cont_groups])
                                exp.designs[key]['main_comparison'].extend(['Control']*len(exp.de_groups[cont_groups]))
                        elif cont_group_type == 'compensation':
                            raise ValueError("Cannot batch compensate 1 Experimental group with 2 Control groups")
                        else:
                            raise ValueError(error)
                        
                        exp.designs[key]['design'] = "~main_comparison"
                        exp.designs[key]['colData'] = pd.DataFrame({"sample_names": exp.designs[key]['all_samples'],
                                                                    "main_comparison": exp.designs[key]['main_comparison']})
                     
                elif exp_group_type == 'compensation':
                    if (Exp1 not in groups) or (Exp2 not in groups):
                        raise ValueError(error)
                    else:
                        exp.designs[key]['all_samples'].extend(exp.de_groups[Exp1])
                        exp.designs[key]['all_samples'].extend(exp.de_groups[Exp2])
                        exp.designs[key]['main_comparison'].extend(['Experimental']*len(exp.de_groups[Exp1] + exp.de_groups[Exp2]))
                        if cont_group_type == 'simple':
                            raise ValueError("Cannot batch compensate 2 Experimental groups with 1 Control groups.")
                        elif cont_group_type == 'compensation':
                            if (Cont1 not in groups) or (Cont2 not in groups):
                                raise ValueError(error)
                            else:
                                exp.designs[key]['all_samples'].extend(exp.de_groups[Cont1])
                                exp.designs[key]['all_samples'].extend(exp.de_groups[Cont2])
                                exp.designs[key]['main_comparison'].extend(['Control']*len(exp.de_groups[Cont1] + exp.de_groups[Cont2]))
                        else:
                            raise ValueError(error)                                     
                        
                        exp.designs[key]['compensation'].extend((['Group_1']*len(exp.de_groups[Exp1]) +
                                                                 ['Group_2']*len(exp.de_groups[Cont1]) +
                                                                 ['Group_1']*len(exp.de_groups[Exp2]) +
                                                                 ['Group_2']*len(exp.de_groups[Cont2])))
                        exp.designs[key]['design'] = "~compensation + main_comparison"
                        exp.designs[key]['colData']= pd.DataFrame({"sample_names": exp.designs[key]['all_samples'],
                                                                   "main_comparison": exp.designs[key]['main_comparison'],
                                                                   "compensation": exp.designs[key]['compensation']})
                else:
                    raise ValueError(error)  

        for name,items in exp.designs.items():
            print('\n{}:'.format(name), file=open(exp.log_file,'a'))
            print(str(items['colData']), file=open(exp.log_file,'a'))

    #Initialize DE sig overlaps
    for comparison, design in exp.designs.items():
        exp.de_sig_overlap[comparison] = (True if yml['Tasks']['Signature_Mode'].lower() == 'combined' else False)
        
    #Overlaps
    if yml['Tasks']['Overlap_of_genes'] == False:
        exp.tasks_complete.append('Overlaps')
        print('\nNot performing signature overlaps', file=open(exp.log_file,'a'))
    elif (yml['Tasks']['Differential_Expression'] == False) and yml['Tasks']['Overlap_of_genes']:
        gene_file=yml['Sig_matrix']
        if os.path.exists(gene_file):
            print("Gene lists found at {}".format(gene_file), file=open(exp.log_file,'a'))
            if gene_file.split('.')[-1] == 'txt':
                gene_fileDF = pd.read_csv(gene_file, header= 0, index_col=None, sep="\t")
            elif (gene_file.split('.')[-1] == 'xls') or (gene_file.split('.')[-1] == 'xlsx'):
                gene_fileDF = pd.read_excel(gene_file)
            else:
                raise IOError("Cannot parse gene lists file.  Make sure it is .txt, .xls, or .xlsx")
            
            overlap_number = 1
            count = 0
            if len(gene_fileDF.columns)%2 == 0:
                for x in range(len(gene_fileDF.columns)/2):
                    overlap_name='Gene_Overlap_{}'.format(str(overlap_number))
                    list_1 = gene_fileDF.columns[count]
                    list_2 = gene_fileDF.columns[count + 1]
                    exp.gene_lists[overlap_name] = {}
                    exp.gene_lists[overlap_name][list_1] = set(gene_fileDF[list_1].tolist())
                    exp.gene_lists[overlap_name][list_2] = set(gene_fileDF[list_2].tolist())
                    count += 2
                    overlap_number = count/2
                print('\nPerforming {} overlaps.'.format(str(count/2)), file=open(exp.log_file,'a'))
            else:
                raise IOError("Cannot parse gene lists file. Requires an even number of gene lists.")

        else:
            raise IOError("Gene List not found. If not doing differential expression, you need to provide an list of genes for overlaps.", file=opne(exp.log_file, 'a'))
    #DE Overlaps
    elif yml['Tasks']['Overlap_of_genes']:
        for key, item in yml['Overlaps'].items():
            if bool(item):   
                exp.overlaps[key] = item.split('v')
        print('\nOverlapping ' + str(len(list(exp.overlaps.keys()))) + ' differential analysis comparison(s).', file=open(exp.log_file, 'a'))
        if str(len(list(exp.overlaps.keys()))) != 0:
            print(str(exp.overlaps)+ '\n', file=open(exp.log_file, 'a'))
    else:
        print("Can't process design for overlaps.  Continuing without overlap analyses.", file=open(exp.log_file, 'a'))
        exp.tasks_complete.append('Overlaps')

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
#BSUB -n 1
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

def stage(exp):
    '''
    Stages files in Pegasus Scratch
    '''
    
    #Stage Experiment Folder in Scratch
    print('Staging in ' + exp.scratch+ '\n', file=open(exp.log_file, 'a'))
    
    #Copy Fastq to scratch fastq folder
    if os.path.exists(exp.scratch + 'Fastq'):
        rmtree(exp.scratch + 'Fastq')
    copytree(exp.fastq_folder, exp.scratch + 'Fastq')

    #change to experimental directory in scratch
    os.chdir(exp.scratch)
    
    exp.fastq_folder= exp.scratch + 'Fastq/'
    
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
                                   mem=1000,
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
                                   threads=4
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
                print('\nTrimming {sample}: '.format(sample=sample), file=open(exp.log_file, 'a'))

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
                                           threads=4
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

def spike(exp):
    '''
    Align sequencing files to ERCC index using STAR aligner.
    '''
    if exp.spike:
        print("Processing with ERCC spike-in: {:%Y-%m-%d %H:%M:%S}\n".format(datetime.now()), file=open(exp.log_file, 'a'))
            
        ERCC_folder=exp.scratch + 'ERCC/'
        os.makedirs(ERCC_folder, exist_ok=True)

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
                                               threads=4
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
            exp.spike = False

        if exp.spike:
            # Prep spike counts for plot (only if Nimer)
            if exp.genome_indicies['ERCC_Mix'] != None:
                # Filtering for counts with more than 5 counts in two samples
                spike_counts = exp.spike_counts.copy()
                spike_counts = spike_counts[spike_counts[spike_counts > 5].apply(lambda x: len(x.dropna()) > 1 , axis=1)]
                mix = pd.read_csv(exp.genome_indicies['ERCC_Mix'], header=0, index_col=1, sep="\t")
                mix = mix.rename(columns={'concentration in Mix 1 (attomoles/ul)': 'Mix_1',
                                          'concentration in Mix 2 (attomoles/ul)': 'Mix_2'})
                names = list(spike_counts.columns)
                spike_counts = spike_counts.join(mix)
                
                merged_spike = pd.DataFrame(columns=['value','Mix_1','Mix_2'])
                name = []
                length = len(spike_counts)
                for sample in names:
                    merged_spike = pd.concat([merged_spike,
                                             spike_counts[[sample,'Mix_1','Mix_2']].rename(columns={sample:'value'})],
                                            ignore_index=True)
                    name=name + [sample]*length
                merged_spike['Sample']=name
                merged_spike['log'] = merged_spike.value.apply(lambda x: np.log2(x))
                merged_spike['log2_Mix_1']=np.log2(merged_spike.Mix_1)
                merged_spike['log2_Mix_2']=np.log2(merged_spike.Mix_2)

                # Plot ERCC spike.
                sns.set(context='paper', font_scale=2, style='white')
                M1 = sns.lmplot(x='log2_Mix_1', y='log', hue='Sample', data=merged_spike, size=10, aspect=1)
                M1.set_ylabels(label='spike-in counts (log2)')
                M1.set_xlabels(label='ERCC Mix (log2(attamoles/ul))')
                plt.title("ERCC Mix 1 Counts per Sample")
                sns.despine()
                M1.savefig(ERCC_folder + 'ERCC_Mix_1_plot.png')
                if __name__ == "__main__":
                    plt.close()

                sns.set(context='paper', font_scale=2, style='white')
                M2 = sns.lmplot(x='log2_Mix_2', y='log', hue='Sample', data=merged_spike, size=10, aspect=1)
                M2.set_ylabels(label='spike-in counts (log2)')
                M2.set_xlabels(label='ERCC Mix (log2(attamoles/ul))')
                plt.title("ERCC Mix 2 Counts per Sample")
                sns.despine()
                M2.savefig(ERCC_folder + 'ERCC_Mix_2_plot.png')
                if __name__ == "__main__":
                    plt.close()

            else:
                print('Not plotting ERCC counts for other labs.', file=open(exp.log_file,'a'))

        print("ERCC spike-in processing complete: {:%Y-%m-%d %H:%M:%S}\n".format(datetime.now()), file=open(exp.log_file, 'a'))
    
    else:
        print("No ERCC spike-in processing.\n", file=open(exp.log_file, 'a'))
    
    exp.tasks_complete.append('Spike')
    return exp 

def bam2bw(in_bam,out_bw,job_log_folder,name,genome):

    script='{job_log_folder}{name}.py'.format(job_log_folder=job_log_folder,name=name)
    print('#!/usr/bin/env python\nimport pybedtools\nimport subprocess', file=open(script,'w'))
    print('kwargs=dict(bg=True,split=True,g="{genome}")'.format(genome=genome), file=open(script,'a'))
    print('readcount=pybedtools.contrib.bigwig.mapped_read_count("{in_bam}")'.format(in_bam=in_bam), file=open(script,'a'))
    print('_scale = 1 / (readcount / 1e6)\nkwargs["scale"] = _scale', file=open(script,'a'))
    print('x = pybedtools.BedTool("{in_bam}").genome_coverage(**kwargs)'.format(in_bam=in_bam), file=open(script,'a'))
    print('cmds = ["bedGraphToBigWig", x.fn, "{genome}", "{out_bw}"]'.format(genome=genome,out_bw=out_bw), file=open(script,'a'))
    print('p = subprocess.Popen(cmds, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)', file=open(script,'a'))
    print('stdout, stderr = p.communicate()', file=open(script,'a'))

    return script

def star(exp):
    '''
    Alignment to genome using STAR and htseq
    '''
    print('\n Beginning RSEM-STAR transcriptome alignments: {:%Y-%m-%d %H:%M:%S}'.format(datetime.now()), file=open(exp.log_file, 'a'))
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
                    spike='STAR --runThreadN 4 --outSAMtype BAM SortedByCoordinate --genomeDir {index} --readFilesIn {fname}_R1.fastq.gz {fname}_R2.fastq.gz --readFilesCommand zcat --outFileNamePrefix {loc}{sample}_ --quantMode GeneCounts'.format(index=exp.genome_indicies['STAR'],fname=fname,loc=out_dir,sample=sample)
                elif exp.seq_type == 'single':
                    spike='STAR --runThreadN 4 --outSAMtype BAM SortedByCoordinate --genomeDir {index} --readFilesIn {fname}.fastq.gz --readFilesCommand zcat --outFileNamePrefix {loc}{sample}_ --quantMode GeneCounts'.format(index=exp.genome_indicies['STAR'],fname=fname,loc=out_dir,sample=sample)

                genome=exp.genome_indicies['chrLen']

                scaled=bam2bw(in_bam='{}{}_Aligned.sortedByCoord.out.bam'.format(out_dir,sample),
                              out_bw='{}{}.star.rpm.bw'.format(out_dir,sample),
                              job_log_folder=exp.job_folder,
                              name='{}_to_bigwig'.format(sample),
                              genome=genome
                             )

                command_list = ['module rm python share-rpms65',
                                'source activate RNAseq',
                                align,
                                'python {}'.format(scaled),
                                ]

                exp.job_id.append(send_job(command_list=command_list, 
                                            job_name= sample + '_STAR',
                                            job_log_folder=exp.job_folder,
                                            q= 'bigmem',
                                            mem=60000,
                                            log_file=exp.log_file,
                                            project=exp.project,
                                            threads=4
                                            )
                                  )

                time.sleep(5)

        #Wait for jobs to finish
        job_wait(id_list=exp.job_id, job_log_folder=exp.job_folder, log_file=exp.log_file)

        scan += 1

    for number,sample in exp.samples.items():
        del_file='{}{}.wig'.format(out_dir, sample,file)
        if os.path.isfile(del_file):
            os.remove(del_file)

    print('STAR alignment to genome finished.', file=open(exp.log_file, 'a'))
    
    ### Generate one matrix for all counts
    try:
        counts_glob = glob.glob(out_dir + '*_ReadsPerGene.out.tab')
        if len(counts_glob) != exp.sample_number:
            print('At least one STAR alignment failed.', file=open(exp.log_file,'a'))
            raise RaiseError('At least one STAR alignment failed. Check scripts and resubmit.')
        else:
            counts = pd.DataFrame(index=pd.read_csv(counts_glob[1], header=None, index_col=0, sep="\t").index)
        
            for number,sample in exp.samples.items():
                exp.counts[sample] = pd.read_csv('{}{}_ReadsPerGene.out.tab'.format(out_dir,sample),header=None, index_col=0, sep="\t")[[3]]
            exp.count_matrix = counts.iloc[4:,:]
            exp.count_matrix.to_csv('{}STAR.count.matrix.txt'.format(out_dir), header=True, index=True, sep="\t")
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
                genome=exp.genome_indicies['chrLen']

                scaled=bam2bw(in_bam='{}{}.genome.sorted.bam'.format(out_dir,sample),
                              out_bw='{}{}.rsem.rpm.bw'.format(out_dir,sample),
                              job_log_folder=exp.job_folder,
                              name='{}_to_bigwig'.format(sample),
                              genome=genome
                             )

                command_list = ['module rm python share-rpms65',
                                'source activate RNAseq',
                                align,
                                'python {}'.format(scaled),
                                plot_model
                                ]

                exp.job_id.append(send_job(command_list=command_list, 
                                            job_name= sample + '_RSEM',
                                            job_log_folder=exp.job_folder,
                                            q= 'bigmem',
                                            mem=60000,
                                            log_file=exp.log_file,
                                            project=exp.project,
                                            threads=4
                                            )
                                  )


                time.sleep(5)



        #Wait for jobs to finish
        job_wait(id_list=exp.job_id, job_log_folder=exp.job_folder, log_file=exp.log_file)
    
        scan += 1

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
        matrix = '{}{}{}.gene.results '.format(matrix,out_dir,sample)
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
    
    counts = pd.read_csv('{}RSEM.count.matrix.txt'.format(out_dir), header=0, index_col=0, sep="\t")
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
                    align = 'kallisto quant --index={index} --output-dir={out} --threads=15 --bootstrap-samples=100 {loc}{sample}_trim_R1.fastq.gz {loc}{sample}_trim_R2.fastq.gz'.format(index=exp.genome_indicies['Kallisto'],out=kal_out,loc=exp.fastq_folder,sample=sample)

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
                                           project=exp.project
                                          )
                                 )

        #Wait for jobs to finish
        job_wait(id_list=exp.job_id, job_log_folder=exp.job_folder, log_file=exp.log_file)
            
        scan += 1

    exp.tasks_complete.append('Kallisto')
    
    return exp

def plot_PCA(counts, colData, out_dir, name):
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

        fig = plt.figure(figsize=(8,8), dpi=100)
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
        plt.subplots_adjust(right=0.9, top=.9)

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

    edaseq = importr('EDASeq')
    as_df=ro.r("as.data.frame")
    assay=ro.r("assay")
    as_cv = ro.r('as.character')
    counts = ro.r("counts")
    fdata=ro.r('fData')
    normCounts=ro.r('normCounts')

    print('Beginning within-lane GC length/content loess normalization for all samples: {:%Y-%m-%d %H:%M:%S}\n'.format(datetime.now()), file=open(exp.log_file,'a'))

    GC_content = pd.read_csv(exp.genome_indicies['GC_Content'], header=0, index_col=0, sep="\t")
    raw_counts = exp.count_matrix
    raw_counts['id']=raw_counts.index
    raw_counts['id']=raw_counts.id.apply(lambda x: x.split("_")[0].split(".")[0])
    GC_genes = GC_content.split.tolist()

    #Keep only counts with GC data (based on latest ensembl biomart).  see EDAseq package and use 'biomart' after dropping ensembl name version.
    GC_counts = round(raw_counts[raw_counts.id.apply(lambda x: x in GC_genes)].drop(columns='id'))
    EDA_set = edaseq.newSeqExpressionSet(counts=GC_counts.values,featureData=GC_content)
    gcNorm = edaseq.withinLaneNormalization(EDA_set, 'gc','loess')
    data_norm = ro.pandas2ri.ri2py_dataframe(normCounts(gcNorm))
    data_norm.index = GC_counts.index
    data_norm.columns = GC_counts.columns
    exp.gc_count_matrix = data_norm

    print('Finished GC normalization: {:%Y-%m-%d %H:%M:%S}\n'.format(datetime.now()), file=open(exp.log_file,'a'))
    exp.tasks_complete.append('GC')

    return exp 

def RUV(RUV_data,design,colData,norm_type,log, ERCC_counts, comparison, plot_dir,de):

    '''
    Performs normalization method for removing unwanted variance across samples using controls sequences.

    Inputs
    ------
    RUV_data = pandas dataframe of counts
    design = string of design (ie '~main_comparison')
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
            plot_PCA(counts = RUV_data, colData=[], out_dir= plot_dir, name= '{}_pre-{}_RUV_ALL_raw_counts_PCA'.format(comparison, norm_type))
        else:
            plot_PCA(counts = RUV_data, colData=colData, out_dir= plot_dir, name= '{}_pre-{}_RUV_raw_counts_PCA'.format(comparison, norm_type))

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

        #Differential expression (Wald DESeq2) to account for scaled variances between samples
        if design == '~main_comparison':
            RUV_dds = deseq.DESeqDataSetFromMatrix(countData=counts(RUVg_set), colData=pdata(RUVg_set), design=ro.Formula('~W_1 + main_comparison'))
            RUV_dds = deseq.DESeq(RUV_dds)
        #Differential expression (Wald DESeq2) to account for scaled variances between samples
        elif design == '~compensation + main_comparison':
            RUV_dds = deseq.DESeqDataSetFromMatrix(countData=counts(RUVg_set), colData=pdata(RUVg_set), design=ro.Formula('~W_1 + compensation + main_comparison'))
            RUV_dds = deseq.DESeq(RUV_dds)
        elif comparison == 'ALL':
            RUV_dds = deseq.DESeqDataSetFromMatrix(countData=counts(RUVg_set), colData=pdata(RUVg_set), design=ro.Formula('~W_1'))


        #generate normalized counts for pca
        counts_df = pandas2ri.ri2py(as_df(normCounts(RUVg_set)))
        counts_df.columns = RUV_data.drop(columns='name').columns

        if comparison == 'ALL':
            plot_PCA(counts = counts_df.dropna(), colData= [], out_dir=plot_dir, name='{}_post-()_RUV_raw_counts_PCA'.format(comparison, norm_type))
        else:
            plot_PCA(counts = counts_df.dropna(), colData= colData, out_dir=plot_dir, name='{}_post-{}_RUV_raw_counts_PCA'.format(comparison, norm_type))

        if de:
            #extract results and relabel samples and genes
            print('{} results type: '.format(comparison), file = open(log, 'a'))
            print(head(results(RUV_dds), 0), file = open(log,'a'))

            results = pandas2ri.ri2py(as_df(deseq.results(RUV_dds, contrast=as_cv(['main_comparison','Experimental','Control']))))
            results.index = RUV_data.name

            lfc = pandas2ri.ri2py(as_df(deseq.lfcShrink(RUV_dds, contrast=as_cv(['main_comparison','Experimental','Control']), type='ashr')))
            lfc.index = RUV_data.name
            print('Switched to ashr method for lfcShrinkage for {}.'.format(comparison), file=open(log,'a'))

        vst = pandas2ri.ri2py_dataframe(assay(deseq.varianceStabilizingTransformation(RUV_dds)))
        vst.columns = RUV_data.drop(columns='name').columns
        vst.index = RUV_data.name

        print('Unwanted variance normalization complete for {} using RUVSeq: {:%Y-%m-%d %H:%M:%S}'.format(comparison, datetime.now()), file=open(log,'a'))

        if de:
            return results, vst, lfc
        else:
            return vst

    except: 
        raise RaiseError('Error during RUVseq.')

def DESeq2(exp):
        
    '''
    Differential Expression using DESeq2
    '''
    print('Beginning DESeq2 differential expression analysis: {:%Y-%m-%d %H:%M:%S}\n'.format(datetime.now()), file=open(exp.log_file, 'a'))
    
    pandas2ri.activate()
    
    deseq = importr('DESeq2')
    as_df=ro.r("as.data.frame")
    as_cv=ro.r('as.character')
    assay=ro.r("assay")
    session=ro.r("sessionInfo")
    results = ro.r('results')
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

        if exp.spike:
            print('Determining ERCC scaling vs Sample scaling using median of ratios of counts for rough comparison.  This may point out potentially problematic samples.'+ '\n', file=open(exp.log_file, 'a'))
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
        else:
            pass

        #Differential Expression
        if exp.norm.lower() == 'median-ratios':
            print('Using DESeq2 standard normalization of scaling by median of the ratios of observed counts.', file=open(exp.log_file, 'a'))
            
            print('Performing Wald Test for differential expression for {}\n'.format(comparison), file=open(exp.log_file, 'a'))
            dds[comparison] = deseq.DESeq(dds[comparison])
            
            print('{} results type: '.format(comparison), file = open(exp.log_file, 'a'))
            print(head(results(dds[comparison]), 0), file = open(exp.log_file,'a'))

            #get results
            exp.de_results['DE2_' + comparison] = pandas2ri.ri2py(as_df(deseq.results(dds[comparison], contrast=as_cv(['main_comparison','Experimental','Control']))))
            exp.de_results['DE2_' + comparison].index = data.index
            
            #get shrunken lfc (apeglm) method)
            exp.de_results['shrunkenLFC_' + comparison] = pandas2ri.ri2py(as_df(deseq.lfcShrink(dds[comparison], coef=as_cv('main_comparison_Experimental_vs_Control'), type='apeglm')))
            exp.de_results['shrunkenLFC_' + comparison].index = data.index
            print('Using apeglm method for lfc shrinkage for {}.'.format(comparison), file=open(exp.log_file, 'a'))

            #variance stabalized transcript counts (also size scaled)
            exp.de_results[comparison + '_vst'] = pandas2ri.ri2py_dataframe(assay(deseq.varianceStabilizingTransformation(dds[comparison])))
            exp.de_results[comparison + '_vst'].columns = data.columns
            exp.de_results[comparison + '_vst'].index = data.index



        elif exp.norm.lower() == 'ercc':
            exp.de_results['DE2_' + comparison], exp.de_results[comparison + '_vst'], exp.de_results['shrunkenLFC_' + comparison]  = RUV(RUV_data = data, 
                                                                                                                                         design=designs['design'], 
                                                                                                                                         colData=colData, 
                                                                                                                                         norm_type='ERCC', 
                                                                                                                                         ERCC_counts = round(exp.spike_counts[designs['all_samples']]), 
                                                                                                                                         log=exp.log_file,
                                                                                                                                         comparison=comparison,
                                                                                                                                         plot_dir = exp.scratch + 'PCA/',
                                                                                                                                         de= True
                                                                                                                                        )
    
        elif exp.norm.lower() == 'empirical':
            exp.de_results['DE2_' + comparison], exp.de_results[comparison + '_vst'], exp.de_results['shrunkenLFC_' + comparison] = RUV(RUV_data = data, 
                                                                                                                                        design=designs['design'], 
                                                                                                                                        colData=colData, 
                                                                                                                                        norm_type='empirical', 
                                                                                                                                        ERCC_counts = None, 
                                                                                                                                        log=exp.log_file,
                                                                                                                                        comparison=comparison,
                                                                                                                                        plot_dir = exp.scratch + 'PCA/',
                                                                                                                                        de=True
                                                                                                                                       )
        else:
            raise RaiseError('Can only use "median-ratios", "ercc", or "empirical" for normalization of DESeq2.')

        #DESeq2 results
        exp.de_results['DE2_' + comparison].sort_values(by='padj', ascending=True, inplace=True)
        exp.de_results['DE2_' + comparison]['gene_name']=exp.de_results['DE2_'+comparison].index
        exp.de_results['DE2_' + comparison]['gene_name']=exp.de_results['DE2_' + comparison].gene_name.apply(lambda x: x.split("_")[1])
        exp.de_results['DE2_' + comparison].to_csv(out_dir + comparison + '-DESeq2-results.txt', 
                                                   header=True, 
                                                   index=True, 
                                                   sep="\t"
                                                  )
        ##Shrunken LFC using apeglm or ashr method
        exp.de_results['shrunkenLFC_' + comparison].sort_values(by='log2FoldChange', ascending=False, inplace=True)
        exp.de_results['shrunkenLFC_' + comparison]['gene_name']=exp.de_results['shrunkenLFC_' + comparison].index
        exp.de_results['shrunkenLFC_' + comparison]['gene_name']=exp.de_results['shrunkenLFC_' + comparison].gene_name.apply(lambda x: x.split("_")[1])
        exp.de_results['shrunkenLFC_' + comparison].to_csv(out_dir + comparison + '-DESeq2-shrunken-LFC.txt', 
                                                           header=True, 
                                                           index=True, 
                                                           sep="\t"
                                                          )

        #Variance Stabilized log2 expected counts.
        exp.de_results[comparison + '_vst'].to_csv(out_dir + comparison + '-VST-counts.txt', 
                                                   header=True, 
                                                   index=True, 
                                                   sep="\t"
                                                  )

    #Variance Stabalized count matrix for all samples.
    colData = pd.DataFrame(index=count_matrix.columns, data={'condition': ['A']*exp.sample_number})
    design=ro.Formula("~1")
    #count_matrix = round(count_matrix[count_matrix[count_matrix > 5].apply(lambda x: len(x.dropna()) > 1 , axis=1)]) 
    dds_all = deseq.DESeqDataSetFromMatrix(countData = count_matrix.values,
                                           colData=colData,
                                           design=design
                                          )
    exp.de_results['all_vst'] = pandas2ri.ri2py_dataframe(assay(deseq.varianceStabilizingTransformation(dds_all)))
    exp.de_results['all_vst'].index=count_matrix.index
    exp.de_results['all_vst'].columns=count_matrix.columns
    exp.de_results['all_vst']['gene_name']=exp.de_results['all_vst'].index
    exp.de_results['all_vst']['gene_name']=exp.de_results['all_vst'].gene_name.apply(lambda x: x.split("_")[1])
    exp.de_results['all_vst'].to_csv(out_dir +'ALL-samples-VST-counts.txt', 
                                     header=True, 
                                     index=True, 
                                     sep="\t"
                                    )

    if exp.norm.lower() == 'ercc':
        exp.de_results['all_ERCC_vst']  = RUV(RUV_data = count_matrix, 
                                              design=design, 
                                              colData=colData, 
                                              norm_type='ERCC', 
                                              ERCC_counts = round(exp.spike_counts[designs['all_samples']]), 
                                              log=exp.log_file,
                                              comparison='ALL',
                                              plot_dir = exp.scratch + 'PCA/',
                                              de=False
                                            )

    print(session(), file=open(exp.log_file, 'a'))    
    exp.tasks_complete.append('DESeq2')
    print('DESeq2 differential expression complete: {:%Y-%m-%d %H:%M:%S}\n'.format(datetime.now()), file=open(exp.log_file, 'a'))
    
    return exp

def Principal_Component_Analysis(exp):

    out_dir = exp.scratch + 'PCA/'
    os.makedirs(out_dir, exist_ok=True)
    
    for comparison,design in exp.designs.items():
        print('Starting DESeq2 VST PCA analysis for {}: {:%Y-%m-%d %H:%M:%S}\n'.format(comparison, datetime.now()), file=open(exp.log_file, 'a'))
        plot_PCA(counts=exp.de_results[comparison + '_vst'],
                 colData= design['colData'],
                 out_dir=out_dir,
                 name=comparison
                )

    print('Starting DESeq2 VST PCA analysis for all samples.', file=open(exp.log_file, 'a'))
    plot_PCA(counts=exp.de_results['all_vst'],
             colData=[],
             out_dir=out_dir,
             name='all_samples'
             )

    print('Starting PCA analysis for all raw counts.', file=open(exp.log_file, 'a'))
    plot_PCA(counts=exp.count_matrix,
             colData=[],
             out_dir=out_dir,
             name='all_raw_counts'
             )

    if exp.gc_norm:
        print('starting PCA analysis for gc normalized raw counts.', file=open(exp.log_file, 'a'))
        plot_PCA(counts = exp.gc_count_matrix,
                 colData=[],
                 out_dir=out_dir,
                 name='gc_nromalized_raw_counts'
                )
    if exp.norm.lower() == 'ercc':
        print('starting PCA analysis for ercc normalized vst counts.', file=open(exp.log_file, 'a'))
        plot_PCA(counts = exp.de_results['all_ERCC_vst'],
                 colData=[],
                 out_dir=out_dir,
                 name='ercc_normalized_vst_counts'
                )

    exp.tasks_complete.append('PCA')
    print('PCA for DESeq2 groups complete: {:%Y-%m-%d %H:%M:%S}\n'.format(datetime.now()), file=open(exp.log_file, 'a'))

    return exp

def Sleuth(exp):
    '''
    Differential expression using sleuth from the Pachter lab: https://pachterlab.github.io/sleuth/
    '''

    pandas2ri.activate()
 
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
        for name in design['colData'].sample_names.tolist():
            path.append(kal_dir + name)

        if 'compensation' in design['colData'].columns.tolist():
            s2c = pd.DataFrame({'sample': design['colData'].sample_names.tolist(),
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
            s2c = pd.DataFrame({'sample': design['colData'].sample_names.tolist(),
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

def volcano(results, sig_up, sig_down, name, out_dir):
    '''
    Generate volcano plot from deseq2 results dataframe and significant genes
    '''

    sns.set(context='paper', style='white', font_scale=1)
    fig = plt.figure(figsize=(6,6), dpi=200)
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

        if exp.de_sig_overlap[comparison]:
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
    Generate heatmap of differentially expressed genes using variance stablized transfrmed log2counts.
    '''
    
    out_dir=exp.scratch + 'Heatmaps/'
    os.makedirs(out_dir, exist_ok=True)
    
    for comparison,design in exp.designs.items():
        vst = exp.de_results['{}_vst'.format(comparison)]
        vst['gene_name']=vst.index
        vst['gene_name']=vst.gene_name.apply(lambda x: x.split("_")[1])

        sig = list(exp.sig_lists[comparison]['2FC_UP'] | exp.sig_lists[comparison]['2FC_DN'])
        if len(sig) == 0:
            print('There are no significantly differentially expressed genes with 2 fold chagnes in {comparison}.  Ignoring heatmap for this group. \n'.format(comparison=comparison), file=open(exp.log_file,'a'))
        else:
            CM = sns.clustermap(vst[vst.gene_name.apply(lambda x: x in sig)].drop('gene_name',axis=1), z_score=0, method='complete', cmap='RdBu_r', yticklabels=False)
            CM.savefig('{}{}_2FC_Heatmap.png'.format(out_dir,comparison), dpi=300)
            CM.savefig('{}{}_2FC_Heatmap.svg'.format(out_dir,comparison), dpi=300)

        sig15 = list(exp.sig_lists[comparison]['15FC_UP'] | exp.sig_lists[comparison]['15FC_DN'])
        if len(sig15) == 0:
            print('There are no significantly differentially expressed genes with 1.5 fold chagnes in {comparison}.  Ignoring heatmap for this group. \n'.format(comparison=comparison), file=open(exp.log_file,'a'))
        else:
            CM15 = sns.clustermap(vst[vst.gene_name.apply(lambda x: x in sig15)].drop('gene_name',axis=1), z_score=0, method='complete', cmap='RdBu_r', yticklabels=False)
            CM15.savefig('{}{}_1.5FC_Heatmap.png'.format(out_dir,comparison), dpi=300)
            CM15.savefig('{}{}_1.5FC_Heatmap.svg'.format(out_dir,comparison), dpi=300)

        sigAll = list(exp.sig_lists[comparison]['All_UP'] | exp.sig_lists[comparison]['All_DN'])
        if len(sigAll) == 0:
            print('There are no significantly differentially expressed genes without a fold change in {comparison}.  Ignoring heatmap for this group. \n'.format(comparison=comparison), file=open(exp.log_file,'a'))
        else:
            CM15 = sns.clustermap(vst[vst.gene_name.apply(lambda x: x in sigAll)].drop('gene_name',axis=1), z_score=0, method='complete', cmap='RdBu_r', yticklabels=False)
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

def GSEA(exp):
    '''
    Perform Gene Set Enrichment Analysis using gsea 3.0 from the Broad Institute.
    '''

    out_dir = '{}DESeq2_GSEA'.format(exp.scratch)
    os.makedirs(out_dir, exist_ok=True)

    print('Starting GSEA enrichment.', file=open(exp.log_file, 'a'))

    for comparison,design in exp.designs.items():
        gsea_file='{out}DESeq2/{comparison}/Hallmarks'

        #check if comparison already done.

        print('GSEA for {comparison} found in {out}DESeq2_GSEA/{comparison}. \n'.format(comparison=comparison, out=exp.out_dir), file=open(exp.log_file, 'a'))
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
        results.sort_values(by='stat', ascending=False, inplace=True)
        results.index = results.gene_name
        results = results.stat.dropna()
        results.to_csv('{}/{}_stat.rnk'.format(out_compare,comparison), header=False, index=True, sep="\t")
        
        rnk = '{}_stat.rnk'.format(comparison)
        rnk2 = '{}_shrunkenLFC.rnk'.format(comparison)

        os.chdir(out_compare)

        print('Beginning GSEA enrichment for {} using preranked genes: {:%Y-%m-%d %H:%M:%S}'.format(comparison, datetime.now()), file=open(exp.log_file, 'a'))
        print('Genes with positive LFC (to the left left in GSEA output graph) are upregulated in experimental vs control conditions. Genes with negative LFC (on right) are downregulated genes in experimental samples vs controls.\n', file=open(exp.log_file,'a'))

        gmts={'h.all': 'Hallmarks',
              'c2.cp.kegg': 'KEGG',
              'c5.bp': 'GO_Biological_Process',
              'c5.mf': 'GO_Molecular_Function',
              'c2.cgp': 'Curated_Gene_Sets'
              }
              
        for gset,name in gmts.items():
            set_dir='{}/{}'.format(out_compare,name) 
            os.makedirs(set_dir, exist_ok=True)

            command_list = ['module rm python java perl share-rpms65',
                            'source activate RNAseq',
                            'java -cp /projects/ctsi/nimerlab/DANIEL/tools/GSEA/gsea-3.0.jar -Xmx2048m xtools.gsea.GseaPreranked -gmx gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/{gset}.v6.1.symbols.gmt -norm meandiv -nperm 1000 -rnk {rnk} -scoring_scheme weighted -rpt_label {comparison}_{gset}_wald -create_svgs false -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 1000 -set_min 10 -zip_report false -out {name} -gui false'.format(gset=gset,comparison=comparison,name=name,rnk=rnk)
                           ] #for shrunken lfc: 'java -cp /projects/ctsi/nimerlab/DANIEL/tools/GSEA/gsea-3.0.jar -Xmx2048m xtools.gsea.GseaPreranked -gmx gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/{gset}.v6.1.symbols.gmt -norm meandiv -nperm 1000 -rnk {rnk2} -scoring_scheme weighted -rpt_label {comparison}_{gset}_shrunkenLFC -create_svgs false -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 1000 -set_min 10 -zip_report false -out {name} -gui false'.format(gset=gset,comparison=comparison,name=name,rnk2=rnk2)

            exp.job_id.append(send_job(command_list=command_list, 
                                       job_name='{}_{}_GSEA'.format(comparison,gset),
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
        for gset,name in gmts.items():
            path=glob.glob('{}/{}/{}/*'.format(out_dir,comparison,name))[0]
            if 'index.html' == '{}/index.html'.format(path).split('/')[-1]:
                os.chdir('{}/{}/{}'.format(out_dir,comparison,name))
                open('Within each folder click "index.html" for results','w')
            else:
                print('GSEA did not complete {} for {}.'.format(name,comparison), file=open(exp.log_file,'a'))            

    os.chdir(exp.scratch)
    exp.tasks_complete.append('GSEA')
    print('GSEA analysis complete: {:%Y-%m-%d %H:%M:%S}\n'.format(datetime.now()), file=open(exp.log_file, 'a'))
    
    return exp

def plot_venn2(Series, string_name_of_overlap, folder):
    '''
    Series with with overlaps 10,01,11
    Plots a 2 way venn.
    Saves to file.
    '''
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

def pipeline():
    exp = parse_yaml()
    exp = preprocess(exp)
    exp = align(exp)
    exp = diff_exp(exp)
    if 'MultiQC' not in exp.tasks_complete:
        exp = final_qc(exp)
    finish(exp)

if __name__ == "__main__":
    plt.switch_backend('agg')
    pipeline()

