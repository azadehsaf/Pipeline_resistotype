import pandas as pd
import csv
import io,sys

#python remove_redundancy_varscan.py tab_to_cleaning final_table
df1= pd.read_csv(sys.argv[1], delimiter="\t")
df_final=df1.drop_duplicates(subset=['POS'], keep='first', low_memory=False)
df_final.to_csv(sys.argv[2], sep='\t', index=False)
