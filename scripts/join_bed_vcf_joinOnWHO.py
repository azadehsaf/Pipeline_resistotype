import pandas as pd
import csv
import io,sys

# read the datasets

#python scripts/join_bed_vcf.py mutation_annote_nom_gene gene.bed fiche_jont.tab
#annotations=pd.read_csv(r"data/interval2023_with_Header.bed",delimiter="\t")
df1_extracted_mutations = pd.read_csv(sys.argv[1], delimiter="\t")
annotations=pd.read_csv(sys.argv[2],delimiter="\t")

# print the datasets

concat_data = pd.merge(df1_extracted_mutations,annotations, how='left')
concat_data.to_csv(sys.argv[3], sep='\t', index=False)


