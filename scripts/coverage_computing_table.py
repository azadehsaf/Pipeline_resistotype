import pandas as pd 
import numpy as np
import csv,re
import os,io,sys

coverage_base=sys.argv[1]
coverage_quality_control=sys.argv[2]

df1=pd.read_csv(coverage_base, delimiter="\t", header=None)

#computing length for each gene
df1["len"]=df1[2]-df1[1]
df1_genes=df1[[3,5,"len"]]
df1_genes=df1_genes.rename(columns={ 5 : 'profondueur' , 3 : 'gène'})
df1_genes['30X_len']=df1_genes['profondueur'].apply(lambda x :  1  if int(x) >= 30 else 0  )
for_coverage=df1_genes.groupby('gène').agg({'len':'count' , 'profondueur' : 'median' , '30X_len': 'sum'})
for_coverage["coverage"]=for_coverage['30X_len']/for_coverage['len']
for_coverage=for_coverage.rename(columns={'profondueur' : 'mediane'})
for_coverage=for_coverage[['len', '30X_len', 'coverage', 'mediane']]

for_coverage.to_csv(coverage_quality_control, sep='\t', index=True)