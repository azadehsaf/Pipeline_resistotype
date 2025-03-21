import pandas as pd 
import numpy as np
import csv,re
import os,io,sys

coverage_base=sys.argv[1]
coverage_quality_control=sys.argv[2]
coverage_deleted_area=sys.argv[3]

path_to_write, filename = os.path.split(sys.argv[2])
all_table=filename.split('_')[0]+"_"+filename.split('_')[1]

path_table_toWriteExcel=path_to_write+"/"+all_table+"_mapping_QC.xlsx"
path_table1_toWrite=sys.argv[2]
path_table2_toWrite=sys.argv[3]

PCmapped=sys.argv[4]
medCoveGenome=sys.argv[5]


#print(coverage_quality_control)
#print()


#echo "Pourcentage de reads mappées: $pc_mapped_read \n" >> {output.mapping_quality}
#echo "Profondeur médiane sur le génome: $pc_mediane_coverage_chrom \n" >> {output.mapping_quality}

dfmapping = pd.DataFrame.from_dict({'Pourcentage de reads mappées:': [PCmapped] , 'Profondeur médiane sur le génome:' :[medCoveGenome] }, orient='index')

with pd.ExcelWriter (path_table_toWriteExcel) as writer:
    dfmapping.to_excel(writer,sheet_name='mapping_quality', index=True)  


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
with pd.ExcelWriter (path_table_toWriteExcel, mode='a') as writer:
    for_coverage.to_excel(writer, sheet_name='couverture_gènes')

#deleted area of genes, result in genomic positions
df1=pd.read_csv(coverage_base,delimiter="\t", header=None)
df1["longueur_gène"]=df1[2]-df1[1]+1
df1_genes=df1[[1,2,3,4,5,"longueur_gène"]]
df1_genes=df1_genes.rename(columns={ 1: 'début' , 2 : 'fin', 5 : 'profondueur', 3 : 'gène', 4:'base_position'})

df1_genes=df1_genes[df1_genes.profondueur==0]
df1_genes["position_génomique"]=df1_genes.apply(lambda x : x[0]+x[3], axis=1)
table=df1_genes.groupby(['gène', 'début', 'fin' , 'longueur_gène' ])

second_table=table.agg(position_début_délétion=( 'base_position', min ),   
                        longueur=('base_position',np.ptp) , 
                        position_génomique_délétion=('position_génomique' , np.min))


second_table.to_csv(coverage_deleted_area, sep='\t', index=True)
with pd.ExcelWriter (path_table_toWriteExcel, mode='a') as writer:
    second_table.to_excel(writer, sheet_name='Deletion_gènes')