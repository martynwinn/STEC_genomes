# Run roary

import os, sys, glob
import numpy as np
import pandas as pd

os.system("roary -v -e --mafft -p 8 -f roary/deer_b1_b2 roary/gff_files/deer_b*.gff")

# pangenome plots 

### create_pan_genome_plots.R 


# build core phylogeny

### fasttree -nt core_gene_alignment.aln > core.tre

# more plots

###python roary_plots.py core.tre gene_presence_absence.csv

# phandango

### drag n drop core.tre gene_presence_absence.csv
