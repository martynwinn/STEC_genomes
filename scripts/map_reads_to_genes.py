# Use bowtie2 to map ENA reads to stx and virulence genes
# Should quality control first
# samtools for analysis

import os, sys, glob
import numpy as np
import pandas as pd

# genes to search for
#ref_indices = 'stx_genes'
ref_indices = 'virulence_all'
os.environ['BOWTIE2_INDEXES'] = '/work4/scd/scarf228/genome_STEC/bowtie_indices'

# best to create pandas dataframe in one go after loop
#isolate_cols = ['name','serotype','stx_Anne','stx2a','stx2c','stx1a','stx1c','stx1d','stx2d','stx2b','stx2e','stx2f','stx2g']
isolate_cols = ['name','serotype','simpPatho','tot_reads']
#isolate_cols = ['name','seqtype','stx_Tim','tot_reads']
#isolate_cols = ['name','serotype','stx_Anne','tot_reads']
#isolate_cols = ['name','Mairi_row','serotype','tot_reads']
for igene in range(957):
    label = 'v'+str(igene+1)
    isolate_cols.append(label)
ncols = len(isolate_cols)

# read in strain metadata
metadata = {}
with open('strain_metadata_clinical2.tab','r') as f:
    line = f.readline()
    while line:
        metadata[line.split()[0]] = [line.split()[1],line.split()[2]]
        line = f.readline()

# STEC samples to search through
#reads_dir = 'deer_moredun/reads_batch_3/'
#reads_dir = 'mince/'
#reads_dir = 'errington/'
reads_dir = 'clinical_claire/batch_2/'

data = []
for isolate_file in glob.glob(reads_dir + '*_1.fastq'):
    #MINCE file_stem = os.path.basename(isolate_file)[:-18]
    file_stem = os.path.basename(isolate_file)[:-8]
    #DEER file_stem = os.path.basename(isolate_file)[:-16]
    #isolate = file_stem[:file_stem.index('_S')]
    isolate = file_stem
    print('Mapping ',isolate)

    pair_1 = reads_dir + file_stem + '_1.fastq'
    pair_2 = reads_dir + file_stem + '_2.fastq'
    out_sam = reads_dir + file_stem + '_bowtie_vall.sam'
    out_bam_sorted = reads_dir + file_stem + '_bowtie_vall_sorted.bam'
    out_log = reads_dir + file_stem + '_bowtie_vall.log'
    samtools_log = reads_dir + file_stem + '_samtools_vall.log'

    isolate_data = [isolate,metadata[isolate][0],metadata[isolate][1]]

    # options
    # --no-unal
    # -a output all, bigger .sam file
    os.system("bowtie2 -x %s -1 %s -2 %s -S %s --no-unal >& %s" % (ref_indices, pair_1, pair_2, out_sam, out_log))

# export LD_LIBRARY_PATH=/work4/scd/scarf228/genome_STEC/
## for libcrypto.so

    # "I had a hunch if I did samtools depth and added up the coverages it would equal the 
    #  output [of samtools bedcov]"
    # also consider samtools coverage, but need newer version (currently v1.7)
    # want to filter by quality
    # samtools idxstats gives similar numbers (slightly less) to grepping the .sam file
    # it will count reads multiple times if use bowtie -a
    os.system("samtools sort -o %s %s" % (out_bam_sorted, out_sam))
    os.system("samtools index %s" % (out_bam_sorted))
    os.system("samtools idxstats %s > %s" % (out_bam_sorted, samtools_log))

    # extract total number of read pairs for use in normalisation
    with open(out_log,'r') as f:
        line = f.readline()
        isolate_data.append(line.split()[0])

    # extract reads aligned to each gene
    with open(samtools_log,'r') as f:
        line = f.readline()
        while line:
            isolate_data.append(line.split()[2])
            line = f.readline()
    print(isolate_data)
    data.append(isolate_data[:ncols])

print(data)
df = pd.DataFrame(data,columns=isolate_cols)

df.to_csv('mapped_reads_virulence_all_clinical2.csv')
