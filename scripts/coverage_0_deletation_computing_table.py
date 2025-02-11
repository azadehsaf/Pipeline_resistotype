```{} 
import pandas as pd 
import numpy as np
import csv,re
import os,io,sys
```


coverage_base=sys.argv[1]
coverage_quality_control=sys.argv[2]

df1=pd.read_csv(coverage_base, delimiter="\t", header=None)

#computing length for each gene
df1["len"]=df1[2]-df1[1]
df1_genes=df1[[3,5,"len"]]
df1_genes=df1_genes.rename(columns={ 5 : 'profondueur', 3 : 'gène', 4:'base_position'})
pos_deb=0
pos_fin=0


df1_genes['zero_profondeur']=df1_genes['profondueur'].apply(lambda x :  1  if x == 0 else 0  )
pos_deb

for_coverage=df1_genes.groupby(['gène']).min()



#for_coverage["coverage"]=for_coverage['zero_profondeur']/for_coverage['len']

for_coverage_median=df1_genes.groupby('gène').agg(lambda x: np.median(x))
for_coverage_median=for_coverage_median.rename(columns={'profondueur' : 'mediane'})
for_coverage["mediane"]=for_coverage_median['mediane']
for_coverage.to_csv(coverage_quality_control, sep='\t', index=True)


