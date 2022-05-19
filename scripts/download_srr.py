# Download reads from ...
# see https://erilu.github.io/python-fastq-downloader/

# SCARF /work4 "1-250TB per dept", not backed up.
# df says /work4/scd is 17% of 227 TB. 

# deer batch 1 zip: 12GB unpacked: 44GB  current: so far, still 44GB
# deer batch 2 zip: 9GB unpacked: 32GB
# deer batch 3 zip: 15GB unpacked: 69GB
# mince zip:  unpacked: 50GB
# errington zip: 
# clinical (489)  guess ca. 800 GB

import os, sys, glob
import numpy as np
import pandas as pd

srr_codes = ['SRR3574026','SRR3242001','SRR3242000','SRR3581423','SRR3581420','SRR3581418','SRR3742209','SRR3581400','SRR3581399','SRR3581397','SRR3581395','SRR3581394','SRR3581393','SRR3581389','SRR3581387','SRR3581384','SRR3581380','SRR3574302','SRR3581377','SRR3581376','SRR3581375','SRR3581373']

#next []

for code in srr_codes:
   print("prefetching ",code)
   os.system("~/SRA_toolkit/sratoolkit.3.0.0-centos_linux64/bin/prefetch %s" % (code))
   print("extracting ",code)
   os.system("~/SRA_toolkit/sratoolkit.3.0.0-centos_linux64/bin/fastq-dump -I --split-files %s/%s.sra" % (code,code))
   print("finished ",code)
