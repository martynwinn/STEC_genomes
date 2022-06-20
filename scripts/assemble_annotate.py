# Use velvet for assembly
# Use prokka for annotation
# Prepare for roary

import os, sys, glob
import numpy as np
import pandas as pd

# isolates
# options: deer, mince, cheese, clinical
isolate_reads = 'clinical'

# best to create pandas dataframe in one go after loop
if isolate_reads == 'deer':
    reads_dir = 'deer_moredun/reads_batch_1/'
elif isolate_reads == 'mince':
    reads_dir = 'mince/'
elif isolate_reads == 'cheese': 
    reads_dir = 'errington/'
elif isolate_reads == 'clinical':
    reads_dir = 'clinical_claire/batch_1/'

for isolate_file in glob.glob(reads_dir + 'SRR2036124*_1.fastq'):
    #DEER file_stem = os.path.basename(isolate_file)[:-16]
    #isolate = file_stem[:file_stem.index('_S')]
    #MINCE file_stem = os.path.basename(isolate_file)[:-18]
    #file_stem = os.path.basename(isolate_file)[:-18]
    #CHEESE file_stem = os.path.basename(isolate_file)[:-8]
    #file_stem = os.path.basename(isolate_file)[:-8]
    #CLINICAL file_stem = os.path.basename(isolate_file)[:-8]
    file_stem = os.path.basename(isolate_file)[:-8]
    isolate = file_stem
    print('Assembling ',isolate)

    pair_1 = reads_dir + file_stem + '_1.fastq'
    pair_2 = reads_dir + file_stem + '_2.fastq'
    velvet_outdir = reads_dir + file_stem + '_velvet'

    # run QC, how is this used?
    os.system("fastqc %s %s" % (pair_1, pair_2))

    # run fastq interlacer? is this necessary if use -separate option for velveth?

    # run velveth
    # Nic suggests hash length 29. Threads?
    os.system("velveth %s 29 -fastq -shortPaired -separate %s %s" % (velvet_outdir, pair_1, pair_2))
    # run velvetg
    os.system("velvetg %s" % (velvet_outdir))

    # run quast

    # run prokka

    # run roary

