
#!/usr/bin/env python
# coding: utf-8

'''

Nimerlab-RNASeq-pipeline-v0.1

Copyright © 2017-2018 Daniel L. Karl

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation 
files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, 
modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the 
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE 
WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR 
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


Reads an experimetnal design yaml file (Version 0.1).
Requires a conda environment 'RNAseq' made from RNAseq.yml
 
To do:
    - Set up for mm10 alignment
    - fastq cat (cat.py preliminary)
    - start and stop points


'''

import os,re,datetime
version=0.1

#### Parse Experimental File
def parse_yaml():
    
    import argparse,yaml
    import pandas as pd
    
    
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--experimental_file', '-f', required=True, help='experimental yaml file', type=str)
    args = parser.parse_args()
    exp_input = open(args.experimental_file,'r')


    yml=yaml.safe_load(exp_input)
    exp_input.close()
    
    #Create an new object to store experimental variables
    exp = type('', (), {})()
    
    #Setting Scratch folder
    exp.scratch = '/scratch/projects/nimerlab/DANIEL/staging/RNAseq/' + yml['Name'] + '/'

    #Setting Job Folder
    exp.job_folder = exp.scratch + 'logs/'
    os.makedirs(exp.job_folder, exist_ok=True)

    #Passing paramters to new object
    exp.date = yml['Rundate']   
    exp.name = yml['Name']
    exp.out_dir = yml['Output_directory']
    
    #Log file
    exp.log_file = exp.job_folder + exp.name + "-" + exp.date + '.out'
    
    print('Pipeline version ' + str(version) + ' run on ' + datetime.datetime.today().strftime('%Y-%m-%d') + '\n', file=open(exp.log_file, 'w'))
    print('Beginning RNAseq Analysis: ' + str(datetime.datetime.now()) + '\n', file=open(exp.log_file, 'a'))
    print('Reading experimental file...' + '\n', file=open(exp.log_file, 'a'))

    #Start Point
    start=[]
    if yml['Startpoint']['Fastq']['Start'] : start.append('Fastq') 
    if yml['Startpoint']['Bam']['Start'] : start.append('Bam')
    if yml['Startpoint']['Gene_Counts']['Start'] : start.append('Counts')
    if len(start) != 1:
        raise ValueError("There are more than one startpoints in this experimental file.  Please fix the file and resubmit.")
    else:
        exp.start = start[0]
        print('Pipeline starting with: ' + str(exp.start)+ '\n', file=open(exp.log_file, 'a'))
   
    #Start Fastq
    exp.fastq = {'Run': False}
    if yml['Startpoint']['Fastq']['Start']:
        exp.fastq['Run'] = True
        exp.fastq['catenated']=yml['Startpoint']['Fastq']['Pre-combined']
        if os.path.isdir(yml['Startpoint']['Fastq']['Fastq_directory']):
            exp.fastq['Folder']=yml['Startpoint']['Fastq']['Fastq_directory']
        else:
            raise IOError("Can't Find Fastq Folder.")

    #Spike
    exp.spike = False
    if yml['ERCC_spike']:
        exp.spike = True
        
    #Initialize trimmed
    exp.trimmed = False
            
    #Start Bam
    exp.bam = {'Run': False}
    if yml['Startpoint']['Bam']['Start']:
        exp.bam['Run'] = True
        if os.path.isdir(yml['Startpoint']['Bam']['Bam_folder']):
            exp.bam['Folder'] = yml['Startpoint']['Bam']['Bam_folder']
        else: 
            raise IOError("Can't Find Bam Folder.")    
        print("Count matrix found at " + yml['Startpoint']['Gene_Counts']['Count_matrix_location']+ '\n', file=open(exp.log_file, 'a'))

    #Start Gene Counts
    exp.counts = {'Run': False}
    if yml['Startpoint']['Gene_Counts']['Start']:
        exp.counts['Run'] = True
        if os.path.exists(yml['Startpoint']['Gene_Counts']['Count_matrix_location']):
            exp.count_matrix = pd.read_csv(yml['Startpoint']['Gene_Counts']['Count_matrix_location'], 
                                           header= 0, index_col=0, sep="\t")
        else:
            raise IOError("Count Matrix Not Found.")    
        print("Count matrix found at " + yml['Startpoint']['Gene_Counts']['Count_matrix_location']+ '\n', file=open(exp.log_file, 'a'))
    
    #End Point
    if yml['Stop']['Alignment']:
        exp.stop = 'Alignment'
        print('Pipeline stopping after alignment.'+ '\n', file=open(exp.log_file, 'a'))
    elif yml['Stop']['Differential_Expression']:
        exp.stop = 'DE'
        print('Pipeline stopping after differential expression analysis.'+ '\n', file=open(exp.log_file, 'a'))
    else:
        exp.stop = 'END'
        print('Pipeline stopping after full analysis.'+ '\n', file=open(exp.log_file, 'a'))
    
    #Genome
    if yml['Genome'] not in ['hg38', 'mm10']:
        raise ValueError("Genome must be either hg38 or mm10.")
    else:
        exp.genome = yml['Genome']
        print('Processing data with: ' + str(exp.genome)+ '\n', file=open(exp.log_file, 'a'))
    
    #Counts
    if not 0 < yml['Total_sample_number'] < 19:
        raise ValueError("This pipeline is only set up to handle up to 18 samples.")
    else:
        exp.sample_number = yml['Total_sample_number']
        print('Processing ' + str(exp.sample_number) + ' samples.'+ '\n', file=open(exp.log_file, 'a'))
    
    #Sample Names
    exp.samples = {}
    count = 1
    for key,name in yml['Samples'].items():
        if count <= exp.sample_number:
            exp.samples[key]=name
            count += 1
        else:
            break
    print("Samples: "+ '\n', file=open(exp.log_file, 'a'))
    print(str(exp.samples) + '\n', file=open(exp.log_file, 'a'))
    
    #Initialize Job_ID list
    exp.job_id=[]
    
    #Out Folder
    os.makedirs(exp.out_dir, exist_ok=True)
    print("Pipeline output folder: " + str(exp.out_dir)+ '\n', file=open(exp.log_file, 'a'))
    
    #Differential Expression Groups
    if exp.stop == 'Alignment':
        pass
    else: 
        exp.de_groups = {}
        for key, item in yml['Differential_Expression_Groups'].items():
            if item == None:
                pass
            else:
                temp=item.split(',')
                exp.de_groups[key] = []
                for x in temp:
                    exp.de_groups[key].append(exp.samples[int(x)])
            
    #Differential Expression Design
    if exp.stop == 'Alignment':
        pass
    else:
        print("Parsing experimental design for differential expression..."+ '\n', file=open(exp.log_file, 'a'))
        
        #Normalization method
        exp.norm = 'bioinformatic'
        if yml['Differential_Expression_Normalizaiton'] == 'ERCC':
            exp.norm = 'ERCC'
            print('Normalizing samples for differential expression analysis using ERCC spike-ins'+ '\n', file=open(exp.log_file, 'a'))
        elif yml['Differential_Expression_Normalizaiton'] == 'bioinformatic':
            print('Normalizing samples for differential expression analysis using conventional size factors'+ '\n', file=open(exp.log_file, 'a'))
        else:
            print("I don't know the " + yml['Differential_Expression_Normalizaiton'] + ' normalization method.  Using size factors.'+ '\n', file=open(exp.log_file, 'a'))
    
        exp.designs={}
        for key, comparison in yml['Differential_Expression_Comparisons'].items():
            if comparison == None:
                pass
            else:
                exp.designs[key]={}
                E = comparison.split('v')[0]
                if len(E.split('-')) == 1:
                    E_type = 1
                elif len(E.split('-')) == 2:
                    E1, E2 = E.split('-')[0], E.split('-')[1]
                    E_type = 2
                else:
                    raise ValueError("Cannot process " + str(key) + ".  Check format E1 or E1-E2. Or too many Groups for pipline.")
                C = comparison.split('v')[1]
                if len(C.split('-')) == 1:
                    C_type = 1
                elif len(C.split('-'))== 2:
                    C1,C2 = C.split('-')[0], C.split('-')[1]
                    C_type = 2
                else:
                    raise ValueError("Cannot process " + str(key) + ".  Check format C1 or C1-C2.  Or too man comparisons.")
            
                #Check comparison for group consistency.
                error = "Can't make a comparison with an unspecified group. Make sure your Comparisons match your Groups for DE"
                groups = list(exp.de_groups.keys())
                
                exp.designs[key]['all_samples']=[]
                exp.designs[key]['main_comparison']=[]
                exp.designs[key]['compensation']=[]
                
                if E_type == 1:
                    if E not in groups:
                        raise ValueError(error)
                    else:
                        exp.designs[key]['all_samples'].extend(exp.de_groups[E])
                        exp.designs[key]['main_comparison'].extend(['Experimental']*len(exp.de_groups[E]))
                        if C_type == 1:
                            if C not in groups:
                                raise ValueError(error)
                            else:
                                exp.designs[key]['all_samples'].extend(exp.de_groups[C])
                                exp.designs[key]['main_comparison'].extend(['Control']*len(exp.de_groups[C]))
                        elif C_type == 2:
                            raise ValueError("Cannot batch compensate 1 Experimental group with 2 Control groups")
                        else:
                            raise ValueError(error)
                        
                        exp.designs[key]['design'] = "~main_comparison"
                        exp.designs[key]['colData'] = pd.DataFrame({"sample_names": exp.designs[key]['all_samples'],
                                                                    "main_comparison": exp.designs[key]['main_comparison']})
                    print('DE design: '+ '\n', file=open(exp.log_file, 'a'))
                    print(str(exp.designs) + '\n', file=open(exp.log_file, 'a')) 
                elif E_type == 2:
                    if E1 not in groups or E2 not in groups:
                        raise ValueError(error)
                    else:
                        exp.designs[key]['all_samples'].extend(exp.de_groups[E1])
                        exp.designs[key]['all_samples'].extend(exp.de_groups[E2])
                        exp.designs[key]['main_comparison'].extend(['Experimental']*len(exp.de_groups[E1] + exp.de_groups[E2]))
                        if C_type == 1:
                            raise ValueError("Cannot batch compensate 2 Experimental groups with 1 Control groups.")
                        elif C_type == 2:
                            if C1 not in groups or C2 not in groups:
                                raise ValueError(error)
                            else:
                                exp.designs[key]['all_samples'].extend(exp.de_groups[C1])
                                exp.designs[key]['all_samples'].extend(exp.de_groups[C2])
                                exp.designs[key]['main_comparison'].extend(['Control']*len(exp.de_groups[C1] + exp.de_groups[C2]))
                        else:
                            raise ValueError(error)                                     
                        
                        exp.designs[key]['compensation'].extend((['Group_1']*len(exp.de_groups[E1]) +
                                                                 ['Group_2']*len(exp.de_groups[C1]) +
                                                                 ['Group_1']*len(exp.de_groups[E2]) +
                                                                 ['Group_2']*len(exp.de_groups[C2])))
                        exp.designs[key]['design'] = "~compensation + main_comparison"
                        exp.designs[key]['colData']= pd.DataFrame({"sample_names": exp.designs[key]['all_samples'],
                                                                   "main_comparison": exp.designs[key]['main_comparison'],
                                                                   "compensation": exp.designs[key]['compensation']})
                    print('DE design: '+ '\n', file=open(exp.log_file, 'a'))
                    print(str(exp.designs)+ '\n', file=open(exp.log_file, 'a')) 
                else:
                    raise ValueError(error)
    
    #DE overlaps
    exp.overlaps = {}
    if yml['Overlaps'] == None:
        print('There are no overlaps to process for muliptle differential expression analyses.'+ '\n', file=open(exp.log_file, 'a'))
        exp.run_overlap = False
    else:
        exp.run_overlap = True
        for key, item in yml['Overlaps'].items():
            if item == None:
                pass    
            else:
                exp.overlaps[key] = item.split('v')
        print('Overlapping ' + str(len(list(exp.overlaps.keys()))) + ' differential analysis comparison(s).'+ '\n', file=open(exp.log_file, 'a'))
        print(str(exp.overlaps)+ '\n', file=open(exp.log_file, 'a'))
        
    #Initialized Process Complete List
    exp.tasks_completed=['Parsed ' + str(datetime.datetime.now())]
    
    print('Experiment file parsed: ' + str(datetime.datetime.now())+ '\n', file=open(exp.log_file, 'a'))
    
    return exp


# Sends job to LSF resource manager queues
def send_job(command_list, job_name, job_log_folder, q, mem):
    
    import os,re,random
    
    os.makedirs(job_log_folder, exist_ok=True)

    '''
    Send job to LSF pegasus.ccs.miami.edu
    Example:
    print(send_job(command_list=['module rm python,
                                  'source activate RNAseq' ,
                                  'fastqc...  ',
                                  'fastq_screen.... ' 
                                 ], 
                   job_name='fastqc',
                   job_log_folder=exp.job_folder,
                   q='bigmem',
                   mem=60000
                  )
         )
    '''

    rand_id = str(random.randint(0, 10000))
    str_comd_list =  '\n'.join(command_list)
    cmd = '''

    #!/bin/bash

    #BSUB -J JOB_{job_name}_ID_{random_number}
    #BSUB -P nimerlab
    #BSUB -o {job_log_folder}/{job_name_o}_logs_{rand_id}.stdout.%J
    #BSUB -e {job_log_folder}/{job_name_e}_logs_{rand_id}.stderr.%J
    #BSUB -W 120:00
    #BSUB -n 1
    #BSUB -q {q}
    #BSUB -R "rusage[mem={mem}]"

    {commands_string_list}'''.format(job_name = job_name,
                                     job_log_folder=job_log_folder,
                                     job_name_o=job_name,
                                     job_name_e=job_name,
                                     commands_string_list=str_comd_list,
                                     random_number=rand_id,
                                     rand_id=rand_id,
                                     q=q,
                                     mem=mem
                                    )
    
    job_path_name = job_log_folder + job_name+'.sh'
    write_job = open(job_path_name, 'a')
    write_job.write(cmd)
    write_job.close()
    print(cmd+ '\n', file=open(exp.log_file, 'a'))
    os.system('bsub < {}'.format(job_path_name))
    print('sending job ID_' + str(rand_id) + '...'+ '\n', file=open(exp.log_file, 'a'))
   
    return rand_id


# waits for LSF jobs to finish
def job_wait(rand_id, job_log_folder):
    
    running = True
    time = 0
    while running:
        jobs_list = os.popen('sleep 10|bhist -w').read()
        print('job  {}'.format(rand_id)+ '\n', file=open(exp.log_file, 'a'))

        if len([j for j in re.findall('ID_(\d+)', jobs_list) if j == rand_id]) == 0:
            running = False
        else:
            time += 10
    

def fastq_cat(exp):
    
    '''
    Not scripted.  Using previously catenated fastq files.
    '''
    
    return exp

# Stages experiment in a scratch folder
def stage(exp):
    
    import subprocess as sub
    from shutil import copytree,rmtree
    
    set_temp='/scratch/projects/nimerlab/tmp'
    sub.run('export TMPDIR=' + set_temp, shell=True)
    print('TMP directory set to ' + set_temp+ '\n', file=open(exp.log_file, 'a'))
    
    #Stage Experiment Folder in Scratch
    os.makedirs(exp.scratch, exist_ok=True)
    print('Staging in ' + exp.scratch+ '\n', file=open(exp.log_file, 'a'))
    
    #Copy Fastq to scratch fastq folder
    if os.path.exists(exp.scratch + 'Fastq'):
        rmtree(exp.scratch + 'Fastq')
    copytree(exp.fastq['Folder'], exp.scratch + 'Fastq')

    #change to experimental directory in scratch
    os.chdir(exp.scratch)
    
    exp.fastq['Folder']= exp.scratch + 'Fastq/'
    
    exp.tasks_completed.append('Stage ' + str(datetime.datetime.now()))
    
    print('Staging complete: ' + str(datetime.datetime.now())+ '\n', file=open(exp.log_file, 'a'))
    
    return exp


def fastqc(exp):
    
    from shutil import move
    import glob
    
    print('Assessing fastq quality.'+ '\n', file=open(exp.log_file, 'a'))
    
    #Make QC folder
    exp.qc_folder = exp.fastq['Folder'] + 'qc/'
    os.makedirs(exp.qc_folder, exist_ok=True)
    
    #Submit fastqc and fastq_screen jobs for each sample
    if exp.trimmed == False:
        
        for number,sample in exp.samples.items():
            command_list = ['module rm python',
                            'module rm perl',
                            'source activate RNAseq',
                            'fastqc ' + exp.fastq['Folder'] + sample + '*',
                           ]

            exp.job_id.append(send_job(command_list=command_list, 
                                       job_name= sample + '_fastqc',
                                       job_log_folder=exp.job_folder,
                                       q= 'general',
                                       mem=1000
                                      )
                             )
    elif exp.trimmed:

        for number,sample in exp.samples.items():
            command_list = ['module rm python',
                            'module rm perl',
                            'source activate RNAseq',
                            'fastqc ' + exp.fastq['Folder'] + sample + '_trim_*',
                           ]

            exp.job_id.append(send_job(command_list=command_list, 
                                       job_name= sample + '_fastqc_trim',
                                       job_log_folder=exp.job_folder,
                                       q= 'general',
                                       mem=1000
                                      )
                             )
    else:
        raise ValueError("Error processing trimming status.")
    
    #Wait for jobs to finish
    for rand_id in exp.job_id:
        job_wait(rand_id=rand_id, job_log_folder=exp.job_folder)
    
    #move to qc folder
    fastqc_files = glob.glob(exp.fastq['Folder'] + '*.zip')
    for f in fastqc_files:
        move(f,exp.qc_folder)
     
    exp.tasks_completed.append('FastQC ' + str(datetime.datetime.now()))
    
    print('FastQC complete: ' + str(datetime.datetime.now())+ '\n', file=open(exp.log_file, 'a'))
    
    return exp


def fastq_screen(exp):
    
    from shutil import move
    import glob
    
    print('Screening for contamination during sequencing: '  + str(datetime.datetime.now())+ '\n', file=open(exp.log_file, 'a'))
    
    #Make QC folder
    exp.qc_folder = exp.fastq['Folder'] + 'qc/'
    os.makedirs(exp.qc_folder, exist_ok=True)
    
    #Submit fastqc and fastq_screen jobs for each sample
    for number,sample in exp.samples.items():
        command_list = ['module rm python',
                        'module rm perl',
                        'source activate RNAseq',
                        'fastq_screen --aligner bowtie2 ' + exp.fastq['Folder'] + sample + '_R1.fastq.gz'
                       ]

        exp.job_id.append(send_job(command_list=command_list, 
                                   job_name= sample + '_fastq_screen',
                                   job_log_folder=exp.job_folder,
                                   q= 'general',
                                   mem=1000
                                  )
                         )
    
    #Wait for jobs to finish
    for rand_id in exp.job_id:
        job_wait(rand_id=rand_id, job_log_folder=exp.job_folder)
    
    #move to qc folder        
    fastqs_files = glob.glob(exp.fastq['Folder'] + '*screen.txt')
    for f in fastqs_files:
        move(f,exp.qc_folder)
    
    exp.tasks_completed.append('Fastq_screen ' + str(datetime.datetime.now()))
    print('Screening complete: ' + str(datetime.datetime.now())+ '\n', file=open(exp.log_file, 'a'))
    
    return exp


# Trimming based on standard UM SCCC Core Nextseq 500 technical errors
def trim(exp):
    
    print('Beginning fastq trimming: '  + str(datetime.datetime.now())+ '\n', file=open(exp.log_file, 'a'))
    
    #Submit fastqc and fastq_screen jobs for each sample
    for number,sample in exp.samples.items():
        print('Trimming {sample}: '.format(sample=sample)+ '\n', file=open(exp.log_file, 'a'))
        
        trim_galore= 'trim_galore --clip_R1 2 --clip_R2 2 --paired --three_prime_clip_R1 4 --three_prime_clip_R2 4 {loc}{sample}_R1.fastq.gz {loc}{sample}_R2.fastq.gz'.format(loc=exp.fastq['Folder'],sample=sample) 
        skewer='skewer --mode pe --end-quality 20 --compress --min 18 --threads 15 -n {loc}{sample}_R1_val_1.fq.gz {loc}{sample}_R2_val_2.fq.gz'.format(loc=exp.fastq['Folder'],sample=sample)
        
        command_list = ['module rm python',
                        'module rm perl',
                        'source activate RNAseq',
                        trim_galore,
                        skewer
                       ]

        exp.job_id.append(send_job(command_list=command_list, 
                                   job_name= sample + '_trim',
                                   job_log_folder=exp.job_folder,
                                   q= 'general',
                                   mem=1000
                                  )
                         )
    
    #Wait for jobs to finish
    for rand_id in exp.job_id:
        job_wait(rand_id=rand_id, job_log_folder=exp.job_folder)
    
    
    for number,sample in exp.samples.items():
        os.rename('{loc}{sample}_R1_val_1.fq-trimmed-pair1.fastq.gz'.format(loc=exp.fastq['Folder'],sample=sample),
                  '{loc}{sample}_trim_R1.fastq.gz'.format(loc=exp.fastq['Folder'],sample=sample)
                 )
        
        os.rename('{loc}{sample}_R1_val_1.fq-trimmed-pair2.fastq.gz'.format(loc=exp.fastq['Folder'],sample=sample),
                  '{loc}{sample}_trim_R2.fastq.gz'.format(loc=exp.fastq['Folder'],sample=sample)
                 )
    
    #move logs to qc folder        
    logs = glob.glob(exp.fastq['Folder'] + '*.txt')
    logs.append(glob.glob(exp.fastq['Folder'] + '*.log'))
    for l in logs:
        move(l,exp.qc_folder) 
    
    exp.tasks_completed.append('Trim ' + str(datetime.datetime.now()))
    print('Trimming complete: ' + str(datetime.datetime.now())+ '\n', file=open(exp.log_file, 'a'))
    
    return exp

def preprocess(exp):
    
    #exp=fastq_cat(exp)
    exp=stage(exp)
    exp=fastqc(exp)
    exp=fastq_screen(exp)
    exp=trim(exp)
    exp=fastqc(exp)
    
    return exp


def spike(exp):
    
    if exp.spike:
        print("Processing with ERCC spike-in: " + str(datetime.datetime.now())+ '\n', file=open(exp.log_file, 'a'))
        
        #Submit STAR alingment for spike-ins for each sample
        for number,sample in exp.samples.items():
            print('Aligning {sample} to spike-in.'.format(sample=sample)+ '\n', file=open(exp.log_file, 'a'))

            spike='STAR --runThreadN 10 --genomeDIR /projects/ctsi/nimerlab/DANIEL/tools/genomes/ERCC_spike/STARIndex --readFilesIn {loc}{sample}_trim_R1.fastq.gz {loc}{sample}_trim_R2.fastq.gz --readFilesCommand zcat --outFileNamePrefix {loc}{sample}_ERCC --quantMode GeneCounts'.format(loc=exp.fastq['Folder'],sample=sample)

            command_list = ['module rm python',
                            'module rm perl',
                            'source activate RNAseq',
                            spike
                           ]

            exp.job_id.append(send_job(command_list=command_list, 
                                       job_name= sample + '_ERCC',
                                       job_log_folder=exp.job_folder,
                                       q= 'general',
                                       mem=5000
                                      )
                             )

        #Wait for jobs to finish
        for rand_id in exp.job_id:
            job_wait(rand_id=rand_id, job_log_folder=exp.job_folder)
        
        #move results to ERCC folder
        os.makedirs(exp.scratch + 'ERCC/', exist_ok=True)
        files = glob.glob(exp.scratch + '*ERCC*.tab')
        for file in files:
            move(file,exp.scratch + 'ERCC/')
        
        ### Generate one matrix for all spike_counts
        matrix='rsem-generate-data-matrix '
        columns=[]
        for number,sample in exp.samples.items():
            matrix = matrix + exp.scratch + 'ERCC/' + '{loc}{sample}_ERCCReadsPerGene.out.tab '.format(loc=exp.scratch + 'ERCC/', sample=sample)
            columns.append(sample)
        
        matrix = matrix + '> {loc}ERCC.count.matrix'.format(loc=exp.scratch + 'ERCC/')
        
        command_list = ['module rm python',
                        'source activate RNAseq',
                        matrix
                       ]

        exp.job_id.append(send_job(command_list=command_list, 
                                   job_name= 'ERCC_Count_Matrix',
                                   job_log_folder=exp.job_folder,
                                   q= 'general',
                                   mem=1000
                                  )
                         )
        
        #Wait for jobs to finish
        for rand_id in exp.job_id:
            job_wait(rand_id=rand_id, job_log_folder=exp.job_folder)
        
        exp.spike_counts = pd.read_csv('{loc}ERCC.count.matrix'.format(loc=exp.scratch + 'ERCC/'),
                                       header=0,
                                       index_col=0,
                                       sep="\t"
                                      )
        
        print("ERCC spike-in processing complete: " + str(datetime.datetime.now())+ '\n', file=open(exp.log_file, 'a'))
        exp.tasks_completed.append('Spike ' + str(datetime.datetime.now()))
        
    else:
        print("No ERCC spike-in processing."+ '\n', file=open(exp.log_file, 'a'))
    
    return(exp)



def rsem(exp):
    
    print('Beginning RSEM-STAR alignments: '  + str(datetime.datetime.now())+ '\n', file=open(exp.log_file, 'a'))
    
    if exp.genome == 'hg38':
        #Submit RSEM-STAR for each sample
        for number,sample in exp.samples.items():
            print('Aligning using STAR and counting transcripts using RSEM for {sample}.'.format(sample=sample)+ '\n', file=open(exp.log_file, 'a'))

            align='rsem-calculate-expression --star --star-gzipped-read-file --paired-end --append-names --output-genome-bam --sort-bam-by-coordinate -p 15 {loc}{sample}_trim_R1.fastq.gz {loc}{sample}_trim_R2.fastq.gz /projects/ctsi/nimerlab/DANIEL/tools/genomes/H_sapiens/Ensembl/GRCh38/Sequence/RSEM_STARIndex/human {sample}'.format(loc=exp.fastq['Folder'],sample=sample)
            bam2wig='rsem-bam2wig {sample}.genome.sorted.bam {sample}.wig {sample}'.format(sample=sample)
            wig2bw='wigToBigWig {sample}.wig /projects/ctsi/nimerlab/DANIEL/tools/genomes/H_sapiens/Ensembl/GRCh38/Sequence/RSEM_STARIndex/chrNameLength.txt {sample}.rsem.bw'.format(sample=sample)
            plot_model='rsem-plot-model {sample} {sample}.models.pdf'

            command_list = ['module rm python',
                            'module rm perl',
                            'source activate RNAseq',
                            align,
                            bam2wig,
                            wig2bw,
                            plot_model
                            ]

            exp.job_id.append(send_job(command_list=command_list, 
                                        job_name= sample + '_RSEM',
                                        job_log_folder=exp.job_folder,
                                        q= 'bigmem',
                                        mem=60000
                                        )
                              )
    else:
        ## MM10 here.
        raise ValueError("This pipeline is not set up for mm10.")
    
    
    #Wait for jobs to finish
    for rand_id in exp.job_id:
        job_wait(rand_id=rand_id, job_log_folder=exp.job_folder)
     
    #make RSEM_results folder
    os.makedirs(exp.scratch + 'RSEM_results/', exist_ok=True)
    
    #move results to folder        
    results = glob.glob(exp.scratch + '*.models.pdf')
    results.append(glob.glob(exp.scratch + '*.genes.results'))
    results.append(glob.glob(exp.scratch + '*.isoforms.results'))
    results.append(glob.glob(exp.scratch + '*.genome.sorted.bam'))
    results.append(glob.glob(exp.scratch + '*.genome.sorted.bam.bai'))
    results.append(glob.glob(exp.scratch + '*.rsem.bw RSEM_results'))
    for file in results:
        move(file,exp.scratch + 'RSEM_results/')

    exp.tasks_completed.append('RSEM ' + str(datetime.datetime.now()))
    print('STAR alignemnt and RSEM counts complete: ' + str(datetime.datetime.now())+ '\n', file=open(exp.log_file, 'a'))
    
    return exp

def kallisto(exp):
    
    print('Beginning Kallisto alignments: '  + str(datetime.datetime.now())+ '\n', file=open(exp.log_file, 'a'))
    
    #make Kallisto_results folder
    os.makedirs(exp.scratch + 'Kallisto_results/', exist_ok=True)
    
    if exp.genome == 'hg38':
        #Submit kallisto for each sample
        for number,sample in exp.samples.items():
            print('Aligning {sample} using Kallisto.'.format(sample=sample)+ '\n', file=open(exp.log_file, 'a'))

            align='kallisto quant --index=/projects/ctsi/nimerlab/DANIEL/tools/genomes/H_sapiens/Ensembl/GRCh38/Sequence/KallistoIndex/GRCh38.transcripts.idx --output-dir={out}Kallisto_results --threads=15 --bootstrap-samples=100 {loc}{sample}_trim_R1.fastq.gz {loc}{sample}_trim_R2.fastq.gz'.format(out=exp.scratch,loc=exp.fastq['Folder'],sample=sample)

            command_list = ['module rm python',
                            'module rm perl',
                            'source activate RNAseq',
                            align
                            ]

            exp.job_id.append(send_job(command_list=command_list, 
                                        job_name= sample + '_Kallisto',
                                        job_log_folder=exp.job_folder,
                                        q= 'bigmem',
                                        mem=60000
                                        )
                              )
    else:
        ## MM10 here.
        raise ValueError("This pipeline is not set up for mm10.")

    #Wait for jobs to finish
    for rand_id in exp.job_id:
        job_wait(rand_id=rand_id, job_log_folder=exp.job_folder)
    
    exp.tasks_completed.append('Kallisto ' + str(datetime.datetime.now()))
    
    return exp


def align(exp):
    
    exp=spike(exp)
    exp=rsem(exp)
    exp=kalliso(exp)
    
    return exp


def count_matrix(exp):
    
    import pandas as pd
    
    print('Generating Sample Matrix from RSEM.gene.results: ' + str(datetime.datetime.now())+ '\n', file=open(exp.log_file, 'a'))

    ### Generate one matrix for all expected_counts
    matrix='rsem-generate-data-matrix '
    columns=[]
    for number,sample in exp.samples.items():
        matrix = matrix + exp.scratch + 'RSEM_results/' + sample + '.genes.results '
        columns.append(sample)
        
    matrix = matrix + '> {loc}RSEM.count.matrix'.format(loc=exp.scratch + 'RSEM_results/')
        
    command_list = ['module rm python',
                    'source activate RNAseq',
                    matrix
                   ]

    exp.job_id.append(send_job(command_list=command_list, 
                               job_name= 'Generate_Count_Matrix',
                               job_log_folder=exp.job_folder,
                               q= 'general',
                               mem=1000
                              )
                     )
    
    #Wait for jobs to finish
    for rand_id in exp.job_id:
        job_wait(rand_id=rand_id, job_log_folder=exp.job_folder)
    
    counts = pd.read_csv('{loc}RSEM.count.matrix'.format(loc=(exp.scratch + 'RSEM_results/')), header=0, index_col=0, sep="\t")
    counts.columns = columns
    
    exp.count_matrix = counts
    exp.tasks_completed.append('Count_Matrix ' + str(datetime.datetime.now()))
    print('Sample count matrix complete: ' + str(datetime.datetime.now())+ '\n', file=open(exp.log_file, 'a'))
    
    return exp
    


def DESeq2(exp):
        
    '''
    Differential Expression using DESeq2
    '''
    
    print('Beginning DESeq2 differential expression analysis: ' + str(datetime.datetime.now())+ '\n', file=open(exp.log_file, 'a'))
    
    import pandas as pd
    import numpy as np
    import rpy2.robjects as ro
    ro.pandas2ri.activate()
    
    deseq = ro.packages.importr('DESeq2')
    as_df=ro.r("as.data.frame")
    assay=ro.r("assay")
    
    out_dir= exp.scratch + 'DESeq2_results/'
    os.makedirs(out_dir, exist_ok=True)
    
    count_matrix = exp.count_matrix
    exp.de_results = {}
    dds={}
    
    for comparison,design in exp.design.items():
        print('Beginning ' + comparison + ': ' + str(datetime.datetime.now())+ '\n', file=open(exp.log_file, 'a'))
        colData=design['colData']
        design=design['design']
        data=count_matrix[design['all_samples']]
        dds[comparison] = deseq.DESeqDataSetFromMatrix(countData = data.values,
                                                       colData=colData,
                                                       design=design
                                                      )
        
        print('Performing differential expression with DESeq2: ' + str(datetime.datetime.now())+ '\n', file=open(exp.log_file, 'a'))
        dds[comparison] = deseq.DESeq(dds[comparison])
        
        if exp.norm == 'ERCC':
            print('Determining ERCC size factors...'+ '\n', file=open(exp.log_file, 'a'))
            ERCC_data = exp.spike_counts[design['all_samples']]
            ERCC_dds = deseq.DESeqDataSetFromMatrix(countData = ERCC_data.values, colData=colData, design=design)
            ERCC_size = deseq.estimateSizeFactors_DESeqDataSet(ERCC_dds)
            sizeFactors=robjects.r("sizeFactors")
            dds[comparison].do_slot('colData').do_slot('listData')[1] = sizeFactors(exp.ERCC_size)
            dds[comparison] = deseq.DESeq(dds[comparison])
        
        #DESeq2 results
        exp.de_results[comparison] = ro.pandas2ri.ri2py(as_df(deseq.results(dds[comparison])))
        exp.de_results[comparison].index = data.index
        exp.de_results[comparison].sort_values(by='padj', ascending=True, inplace=True)
        exp.de_results[comparison]['gene_name']=exp.de_results[comparison].index
        exp.de_results[comparison].to_csv(out_dir + comparison + '-DESeq2-results.txt', 
                                          header=True, 
                                          index=True, 
                                          sep="\t"
                                         )
        #Variance Stabilized log2 expected counts.
        exp.de_results[comparison + '_vst'] = ro.pandas2ri.ri2py_dataframe(assay(deseq.varianceStabilizingTransformation(dds[comparison])))
        exp.de_results[comparison + '_vst'].columns = data.columns
        exp.de_results[comparison + '_vst'].index = data.index
        exp.de_results[comparison + '_vst'].to_csv(out_dir + comparison + '-VST-counts.txt', 
                                                   header=True, 
                                                   index=True, 
                                                   sep="\t"
                                                  )
        
    exp.tasks_completed.append('DESeq2 ' + str(datetime.datetime.now()))
    print('DESeq2 differential expression complete: ' + str(datetime.datetime.now())+ '\n', file=open(exp.log_file, 'a'))
    
    return exp


def clustermap(exp):
    
    import pandas as pd
    import seaborn as sns
    
    os.makedirs(exp.scratch + 'DESeq2_results/Heatmaps/', exist_ok=True)
    
    for comparison,design in exp.design.items():
        results=exp.de_results[comparison]
        vst = exp.de_results[comparison + '_vst']
        sig = results[(results.padj < 0.05) & ((results.log2FoldChange > 1) | (results.log2FoldChange < -1))].gene_name.tolist()
        CM = sns.clustermap(vst[vst.gene_name.apply(lambda x: x in sig)], z_score=0, method='complete', cmap='RdBu_r')
        CM.savefig(exp.scratch + 'DESeq2_results/Heatmaps/{comparison}_2FC_Heatmap.png'.format(comparison=comparison), dpi=200)
        CM.savefig(exp.scratch + 'DESeq2_results/Heatmaps/{comparison}_2FC_Heatmap.png'.format(comparison=comparison), dpi=200)
    
        sig15 = results[(results.padj < 0.05) & ((results.log2FoldChange > 0.585) | (results.log2FoldChange < -0.585))].gene_name.tolist()
        CM15 = sns.clustermap(vst[vst.gene_name.apply(lambda x: x in sig15)], z_score=0, method='complete', cmap='RdBu_r')
        CM15.savefig(exp.scratch + 'DESeq2_results/Heatmaps/{comparison}_1.5FC_Heatmap.png'.format(comparison=comparison), dpi=200)
        CM15.savefig(exp.scratch + 'DESeq2_results/Heatmaps/{comparison}_1.5FC_Heatmap.svg'.format(comparison=comparison), dpi=200)
    
    exp.tasks_completed.append('DESeq2_Heatmaps ' + str(datetime.datetime.now()))
    print('Heatmaps for DESeq2 differentially expressed genes complete: ' + str(datetime.datetime.now())+ '\n', file=open(exp.log_file, 'a'))
    
    return exp


def enrichr_de(exp):
    import pandas as pd
    import gseapy
    
    out_dir = exp.scratch + 'DESeq2_results/enrichr'
    os.makedirs(out_dir, exist_ok=True)
    
    for comparison,design in exp.design.items():
        print('Beginning GO enrichment for {comparison}: '.format(comparison=comparison) + str(datetime.datetime.now())+ '\n', file=open(exp.log_file, 'a'))
        
        results=exp.de_results[comparison]
        results['gene_name']=results.gene_name.apply(lambda x: x.split("_")[1])
        
        exp.sig_lists={}
        exp.sig_lists[comparison] = {}
        exp.sig_lists[comparison]['2FC_UP'] = set(results[(results.padj < 0.05) & (results.log2FoldChange > 1)].gene_name.tolist())
        exp.sig_lists[comparison]['2FC_DN'] = set(results[(results.padj < 0.05) & (results.log2FoldChange < -1)].gene_name.tolist())
        exp.sig_lists[comparison]['15FC_UP'] = set(results[(results.padj < 0.05) & (results.log2FoldChange > .585)].gene_name.tolist())
        exp.sig_lists[comparison]['15FC_DN'] = set(results[(results.padj < 0.05) & (results.log2FoldChange < -.585)].gene_name.tolist())

        for name,sig in exp.sig_lists[comparison].items():
            gseapy.enrichr(gene_list=list(sig), 
                           description='{comparison}_{name}_KEGG'.format(comparison=comparison,name=name),
                           gene_sets='KEGG_2016', 
                           outdir=out_dir
                          )
            gseapy.enrichr(gene_list=list(sig), 
                           description='{comparison}_{name}_GO_biological_process'.format(comparison=comparison,name=name), 
                           gene_sets='GO_Biological_Process_2017b', 
                           outdir=out_dir
                          )
            gseapy.enrichr(gene_list=list(sig), 
                           description='{comparison}_{name}_GO_molecular_function'.format(comparison=comparison,name=name), 
                           gene_sets='GO_Molecular_Function_2017b', 
                           outdir=out_dir
                          )
                
    
    exp.tasks_completed.append('Enrichr_DE ' + str(datetime.datetime.now()))
    print('Enrichment analysis for DESeq2 differentially expressed genes complete: ' + str(datetime.datetime.now())+ '\n', file=open(exp.log_file, 'a'))
    
    return exp


def GSEA(exp):
    
    import pandas as pd
    import geseapy
    
    out_dir = exp.scratch + 'DESeq2_results/GSEA'
    os.makedirs(out_dir, exist_ok=True)
    
    for comparison,design in exp.design.items():
        print('Beginning GSEA enrichment for {comparison}: '.format(comparison=comparison) + str(datetime.datetime.now())+ '\n', file=open(exp.log_file, 'a'))
        results=exp.de_results[comparison]
        results['gene_name']=results.gene_name.apply(lambda x: x.split("_")[1])
        results.sort_values(by='stat', ascending=False, inplace=True)
        
        print('Beginning GSEA:Hallmark enrichment for {comparison}: '.format(comparison=comparison) + str(datetime.datetime.now())+ '\n', file=open(exp.log_file, 'a'))
        gseapy.prerank(rnk= results.stat, gene_sets= '/projects/ctsi/nimerlab/DANIEL/tools/gene_sets/h.all.v6.1.symbols.gmt', outdir=out_dir)
        print('Beginning GSEA:KEGG enrichment for {comparison}: '.format(comparison=comparison) + str(datetime.datetime.now())+ '\n', file=open(exp.log_file, 'a'))
        gseapy.prerank(rnk= results.stat, gene_sets= '/projects/ctsi/nimerlab/DANIEL/tools/gene_sets/c2.cp.kegg.v6.1.symbols.gmt', outdir=out_dir)
        print('Beginning GSEA:GO biological process enrichment for {comparison}: '.format(comparison=comparison) + str(datetime.datetime.now())+ '\n', file=open(exp.log_file, 'a'))
        gseapy.prerank(rnk= results.stat, gene_sets= '/projects/ctsi/nimerlab/DANIEL/tools/gene_sets/c5.bp.v6.1.symbols.gmt', outdir=out_dir)
        print('Beginning GSEA:GO molecular function enrichment for {comparison}: '.format(comparison=comparison) + str(datetime.datetime.now())+ '\n', file=open(exp.log_file, 'a'))
        gseapy.prerank(rnk= results.stat, gene_sets= '/projects/ctsi/nimerlab/DANIEL/tools/gene_sets/c5.mf.v6.1.symbols.gmt', outdir=out_dir)
        print('Beginning GSEA:Perturbation enrichment for {comparison}: '.format(comparison=comparison) + str(datetime.datetime.now())+ '\n', file=open(exp.log_file, 'a'))
        gseapy.prerank(rnk= results.stat, gene_sets= '/projects/ctsi/nimerlab/DANIEL/tools/gene_sets/c2.cgp.v6.1.symbols.gmt', outdir=out_dir)

    exp.tasks_completed.append('GSEA_DESeq2 ' + str(datetime.datetime.now()))
    print('GSEA for DESeq2 stat preranked genes complete: ' + str(datetime.datetime.now())+ '\n', file=open(exp.log_file, 'a'))
    
    return exp


def PCA(exp):
    
    import pandas as pd
    from sklearn.decomposition import PCA
    import matplotlib.pyplot as plt 
    import matplotlib.patches as mpatches
    
    out_dir = exp.scratch + 'DESeq2_results/PCA/'
    os.makedirs(out_dir, exist_ok=True)
    
    for comparison,design in exp.design.items():
        print('Starting PCA analysis for {comparison}: '.format(comparison=comparison) + str(datetime.datetime.now())+ '\n', file=open(exp.log_file, 'a'))
        pca = PCA(n_components=2)
        bpca = bpca=pca.fit_transform(exp.de_results[comparison + '_vst'].T)
        pca_score = pca.explained_variance_ratio_
        bpca_df = pd.DataFrame(bpca)
        bpca_df.index = exp.de_results[comparison + '_vst'].T.index
        bpca_df['group']= design['colData']['main_comparison'].tolist()
        bpca_df['name']= design['colData']['sample_names'].tolist()
            
        plt.clf()
        fig = plt.figure(figsize=(8,8), dpi=100)
        ax = fig.add_subplot(111)
        ax.scatter(bpca_df[bpca_df.group == 'Experimental'][0],bpca_df[bpca_df.group == 'Experimental'][1], marker='o', color='blue')
        ax.scatter(bpca_df[bpca_df.group == 'Control'][0],bpca_df[bpca_df.group == 'Control'][1], marker='o', color='red')
        ax.set_xlabel('PCA Component 1: {var}% variance'.format(var=int(pca_score[0]*100))) 
        ax.set_ylabel('PCA Component 2: {var}% varinace'.format(var=int(pca_score[1]*100)))
        red_patch = mpatches.Patch(color='red', alpha=.4, label='Control')
        blue_patch = mpatches.Patch(color='blue', alpha=.4, label='Experimental')

        for i,sample in enumerate(bpca_df['name'].tolist()):
            ax.annotate(sample, (bpca_df.iloc[i,0], bpca_df.iloc[i,1]), textcoords='offset points')             
        ax.legend(handles=[blue_patch, red_patch], loc=1)
        ax.figure.savefig(out_dir, '{comparison}_PCA.png'.format(comparison=comparison))
        ax.figure.savefig(out_dir, '{comparison}_PCA.svg'.format(comparison=comparison))

    exp.tasks_completed.append('PCA ' + str(datetime.datetime.now()))
    print('PCA for DESeq2 groups complete: ' + str(datetime.datetime.now())+ '\n', file=open(exp.log_file, 'a'))

    return exp


def diff_exp(exp):
    
    exp = count_matrix(exp)
    exp = spike_norm(exp)
    exp = DESeq2(exp)
    exp = clustermap(exp)
    exp = enrichr_de(exp)
    exp = GSEA(exp)
    exp = PCA(exp)

    #Sleuth
    #ICA  

    return exp


def plot_venn2(Series, string_name_of_overlap, folder):
    '''
    Series with with overlaps 10,01,11
    Plots a 2 way venn.
    Saves to file.
    '''
    
    from matplotlib_venn import venn2
    import matplotlib.pyplot as plt
    
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
        venn_plot.get_patch_by_id(patch).set_color(color)
        venn_plot.get_patch_by_id(patch).set_alpha(0.4)
        venn_plot.get_patch_by_id(patch).set_edgecolor('none')    
     
    plt.title(string_name_of_overlap + " Overlaps")
    plt.tight_layout()
    plt.savefig(folder + string_name_of_overlap + "-overlap-" + datetime.datetime.today().strftime('%Y-%m-%d') + ".svg", dpi=200)
    plt.savefig(folder + string_name_of_overlap + "-overlap-" + datetime.datetime.today().strftime('%Y-%m-%d') + ".png", dpi=200)


def overlaps(exp):
    
    import pandas as pd
    import gseapy
    
    out_dir = exp.scratch + 'Overlaps/'
    os.makedirs(out_dir, exist_ok=True)
    
    exp.overlap_results={}
    names=['2FC_UP', '2FC_DN', '15FC_UP','15FC_DN']
      
    print('Beginning overlap of significant genes: ' + str(datetime.datetime.now())+ '\n', file=open(exp.log_file, 'a'))

    for overlap,comparison_list in exp.overlaps.items():
        for name in names:
            key= '{overlap}_{name}'.format(overlap=overlap,name=name)
            exp.overlap_results[key] = exp.sig_lists[comparison_list[0]][name] & exp.sig_lists[comparison_list[1]][name] 
            venn = pd.Series([len([comparison_list[0]][name])-len(exp.overlap_results[key]),
                              len([comparison_list[1]][name])-len(exp.overlap_results[key]),
                              len(exp.overlap_results[key])
                             ],
                             index= [comparison_list] + ['Overlap']
                            )
            plot_venn2(venn, key, out_dir)
                   
    for name,sig in exp.overlap_results.items():
        print('Perfomring GO enrichment for ' + name + ' overlaps: ' + str(datetime.datetime.now())+ '\n', file=open(exp.log_file, 'a'))
        gseapy.enrichr(gene_list=list(sig),
                       description='{name}_overlap_KEGG'.format(name=name),
                       gene_sets='KEGG_2016', 
                       outdir=out_dir
                      )
        gseapy.enrichr(gene_list=list(sig),
                       description='{name}_overlap_GO_biological_process'.format(name=name),
                       gene_sets='GO_Biological_Process_2017b', 
                       outdir=out_dir
                      )
        gseapy.enrichr(gene_list=list(sig), 
                       description='{name}_overlap_GO_molecular_function'.format(name=name),
                       gene_sets='GO_Molecular_Function_2017b', 
                       outdir=out_dir
                      )
    exp.tasks_completed.append('Overlaps ' + str(datetime.datetime.now()))
    print('Overlap analysis complete: ' + str(datetime.datetime.now())+ '\n', file=open(exp.log_file, 'a'))
                   
    return exp


def final_qc(exp):
    
    print('Beginning final qc: ' + str(datetime.datetime.now())+ '\n', file=open(exp.log_file, 'a'))
    
    command_list = ['module rm python',
                    'source activate RNAseq',
                    'multiqc {folders}'.format(folders=exp.scratch)
                   ]
    
    exp.job_id.append(send_job(command_list=command_list, 
                               job_name= 'MultiQC',
                               job_log_folder=exp.job_folder,
                               q= 'general',
                               mem=1000
                              )
                     )
    
    #Wait for jobs to finish
    for rand_id in exp.job_id:
        job_wait(rand_id=rand_id, job_log_folder=exp.job_folder)
    
    exp.tasks_completed.append('MultiQC ' + str(datetime.datetime.now()))
    
    return exp


def finish(exp):
    
    import watermark
    import pickle
    import shutil
    ## add wattermark, python and r versions to exp.
    ## print them here.
    
    filename= '{out}{name}_{date}.pkl'.format(out=exp.scratch, name=exp.name, date=exp.date)
    experiment = open(filename, 'a') 
    pickle.dump(exp, experiment) 
    
    #move all files from scratch to out directory
    
    shutil.move(exp.scratch, exp.outdir)
    
    print('{name} analysis complete!  Performed the following tasks: '.format(name=exp.name)+ '\n', file=open(exp.log_file, 'a'))
    print(str(exp.tasks_completed) + '\n', file=open(exp.log_file, 'a'))
    print('Moved all files into {out}: '.format(out=exp.outdir) + str(datetime.datetime.now())+ '\n', file=open(exp.log_file, 'a'))
    print("Finger's Crossed!!!"+ '\n', file=open(exp.log_file, 'a'))
    

### Pipeline:

exp = parse_yaml()
exp = preprocess(exp)
exp = align(exp)
exp = diff_exp(exp) #finetune for cluster map
if exp.run_overlap:
    exp=Overlaps(exp)
exp = final_qc(exp)
finish(exp)
