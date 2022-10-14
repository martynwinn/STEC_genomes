# Use velvet or shovill for assembly
# Use prokka for annotation
# Prepare for roary

import os, sys, glob
import numpy as np
import pandas as pd

# isolates
# options: deer, mince, cheese, clinical
isolate_reads = 'deer'

# best to create pandas dataframe in one go after loop
if isolate_reads == 'deer':
    reads_dir = 'deer_moredun/reads_batch_2/'
elif isolate_reads == 'mince':
    reads_dir = 'mince/'
elif isolate_reads == 'cheese': 
    reads_dir = 'errington/'
elif isolate_reads == 'clinical':
    reads_dir = 'clinical_claire/batch_1/'

do_fastqc = False
do_shovill = True

for isolate_file in glob.glob(reads_dir + '24131wH3*_1_trimmed.fastq'):
    #DEER 
    file_stem = os.path.basename(isolate_file)[:-16]
    #MINCE file_stem = os.path.basename(isolate_file)[:-18]
    #file_stem = os.path.basename(isolate_file)[:-18]
    #CHEESE file_stem = os.path.basename(isolate_file)[:-8]
    #file_stem = os.path.basename(isolate_file)[:-8]
    #CLINICAL file_stem = os.path.basename(isolate_file)[:-8]
    #file_stem = os.path.basename(isolate_file)[:-8]
    isolate = file_stem

    pair_1 = reads_dir + file_stem + '_1_trimmed.fastq'
    pair_2 = reads_dir + file_stem + '_2_trimmed.fastq'
    trim_tag = '_shov'
    pair_trimmed = reads_dir + file_stem + trim_tag + '.fastq'
    trimmomatic_log = reads_dir + file_stem + trim_tag + '.log'
    pair_1_trimmed = reads_dir + file_stem + trim_tag + '_1P.fastq'
    pair_2_trimmed = reads_dir + file_stem + trim_tag + '_2P.fastq'
    assembly_outdir = reads_dir + file_stem + '_assembly' + trim_tag
    quast_outdir = reads_dir + file_stem + '_quast' + trim_tag
    prokka_outdir = reads_dir + file_stem + '_prokka' + trim_tag

    # run FastQC, this produces a report to be manually inspected
    if do_fastqc:
        print('Running FastQC on ',isolate)
        os.system("fastqc %s %s" % (pair_1, pair_2))

    if do_shovill:
        print('Running shovill on ',isolate)
        # assembly with Seemann's shovill
        ###os.system("shovill --outdir %s --tmpdir /work4/scd/scarf228/genome_STEC/tmp --force --R1 %s --R2 %s" % (assembly_outdir,pair_1,pair_2))

    else:
        # try Trimmomatic instead
        print('Running Trimmomatic on ',isolate)
        # ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 
        os.system("java -jar ~/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar PE -basein %s -baseout %s -trimlog %s SLIDINGWINDOW:4:20" % (pair_1,pair_trimmed,trimmomatic_log))

        # run fastq interlacer? is this necessary if use -separate option for velveth?

        # run velveth
        # Nic suggests hash length 29. Threads?
        print('Running Velvet on ',isolate)
        os.system("velveth %s 29 -fastq -shortPaired -separate %s %s" % (assembly_outdir, pair_1_trimmed, pair_2_trimmed))
        # run velvetg
        os.system("velvetg %s" % (assembly_outdir))

    # run quast for assembly quality
    print('Running Quast on ',isolate)
    ###os.system("quast -o %s %s" % (quast_outdir,assembly_outdir+'/contigs.fa'))

    # run prokka
    print('Running prokka on ',isolate)
    # --compliant necessary because contig names from velvet too long
    os.system("prokka --usegenus --genus Escherichia --species coli --outdir %s --force --compliant %s" % (prokka_outdir,assembly_outdir+'/contigs.fa'))
    os.system("mv %s %s" % (prokka_outdir+'/PROKKA_08122022.gff',prokka_outdir+'/'+isolate+'.gff'))
    os.system("ln -s %s %s" % ('../../'+prokka_outdir+'/'+isolate+'.gff','roary/gff_files/deer_b2_'+isolate+'.gff')) 

    # run roary
    print('Running Roary on ',isolate)
    ###os.system("roary -p 8 %s" % (gff_file))

