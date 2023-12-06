import pandas as pd
import csv
import io,sys

# read the datasets



df1_extracted_mutations = pd.read_csv(sys.argv[1], delimiter="\t")
annotations=pd.read_csv(sys.argv[2],delimiter="\t")

# print the datasets

concat_data = pd.merge(df1_extracted_mutations,annotations, how='left')
concat_data.to_csv(sys.argv[3], sep='\t', index=False)


