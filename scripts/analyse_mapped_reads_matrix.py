# This follows on from map_reads_to_genes.py

import os, sys, glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")   # don't need DISPLAY
import matplotlib.pyplot as plt
import matplotlib.collections as collections

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# read virulence families
vfam = pd.read_csv('UniqueVFs_List.csv')
vir_gene_families = {}
for index,row in vfam.iterrows():
   if index > 0:
      vir_gene_families[name] = (number,row["MW_VF_Numbering"]-1)
   name = row["Name"]
   number = row["MW_VF_Numbering"]
vir_gene_families[name] = (number,957)
print(vir_gene_families)

#vir_gene_families = {
#   'cba':(18,32),
#   'intimin':(96,140),
#   'efa1':(144,154),
#   'espI':(209,210),
#   'gad':(249,318),
#   'ehxA':(319,330),
#   'pic':(530,535),
#   'sat':(537,542),
#   'senB':(543,545),
#   'sepA':(546,552),
#   'stx1A':(561,578),
#   'stx1B':(579,592),
#   'stx2A':(593,706),
#   'stx2B':(707,749),
#   'tccP':(756,789),
#   'tir':(790,825)
#}

frames = []
# want input files in particular order
infile_list = ['mapped_reads_virulence_all_reads1.csv','mapped_reads_virulence_all_reads2.csv','mapped_reads_virulence_all_reads3.csv','mapped_reads_virulence_all_mince.csv','mapped_reads_virulence_all_cheese.csv','mapped_reads_virulence_all_clinical1.csv','mapped_reads_virulence_all_clinical2.csv']
keys=['deer1','deer2','deer3','mince','cheese','clinical1','clinical2']
for infile in infile_list:
   print('Reading mapped reads file %s' % infile)
   frames.append(pd.read_csv(infile))
df_mapped = pd.concat(frames, keys=keys)
#print(df_mapped['v96'].to_numpy())

# normalise each row by the total number of reads and drop metadata
# multiply by a million to make nice numbers
### also normalise by length of ref gene? 
df_norm = df_mapped.loc[:,'v1':'v957'].div(df_mapped.loc[:,'tot_reads']*0.000001,axis='index')
print(df_norm)
print(df_norm['v799'].sum(axis='index'))

# new columns for families
for family in vir_gene_families.keys():
   df_norm[family] = df_norm.iloc[:,vir_gene_families[family][0]-1:vir_gene_families[family][1]].sum(axis=1)
print(df_norm)

#selected correlations
#TODO need to swap with drop column bit
cor_matrix = df_norm.loc[:,'astA':'air'].corr()
pd.set_option('display.max_columns', 100)
print(cor_matrix)
pd.reset_option('display.max_columns')

# drop columns where no reads map for any sample
for col in df_norm.columns:
   if col[0] == 'v' and df_norm[col].sum(axis='index') < 0.001:
      df_norm.drop([col], axis = 1, inplace=True)

print(df_norm)

# plot some columns
fam1 = 'eae'
fam2 = 'tir'
fig, (ax1) = plt.subplots(1,1)
g = np.arange(0,229)
ax1.plot(g,df_norm[fam1].to_numpy(),label=fam1)
ax1.plot(g,df_norm[fam2].to_numpy(),label=fam2)
ax1.set_xlabel('Samples')
ax1.set_ylabel('Mapped reads')
ax1.legend()

ax1.text(30,1800,'Deer',color='red')
ax1.text(95,1900,'Mince',color='red')
ax1.text(110,2150,'Cheese',color='red')
ax1.text(150,2350,'Clinical',color='red')
collection = collections.BrokenBarHCollection(
    [(90,31),(125,104)], (0,3000), facecolor='green', alpha=0.1)
ax1.add_collection(collection)

#for i,label in enumerate(df_mapped['name']):
#   ax1.annotate(label,(x_pca[i,0],x_pca[i,1]))

plt.savefig("vall_reads.jpg")
plt.close()

exit(1)

# scale mapped reads to unit variance
#scaler = StandardScaler()
#scaler.fit(df_mapped.iloc[:,5:])
#scaled_data = scaler.transform(df_mapped.iloc[:,5:])

# do PCA
pca = PCA(n_components=10)
pca.fit(df_norm.loc[:,'v1':'v957'])
print(pca.explained_variance_ratio_)
x_pca = pca.transform(df_norm.loc[:,'v1':'v957'])
#x_pca_deer = pca.transform(df_norm.loc[:90,'v1':'v957'])
#x_pca_mince = pca.transform(df_norm.loc[90:121,'v1':'v957'])
#x_pca_cheese = pca.transform(df_norm.loc[121:125,'v1':'v957'])
#x_pca_clinical = pca.transform(df_norm.loc[126:,'v1':'v957'])

fig, (ax1) = plt.subplots(1,1)
ax1.scatter(x_pca[:90,0],x_pca[:90,1],c='red')
ax1.scatter(x_pca[90:121,0],x_pca[90:121,1],c='blue')
ax1.scatter(x_pca[121:125,0],x_pca[121:125,1],c='green')
ax1.scatter(x_pca[126:,0],x_pca[126:,1],c='yellow')
ax1.set_xlabel('First principal component')
ax1.set_ylabel('Second Principal Component')
#for i,label in enumerate(df_mapped['name']):
#   ax1.annotate(label,(x_pca[i,0],x_pca[i,1]))

plt.savefig("pca_vall.jpg")
plt.close()
