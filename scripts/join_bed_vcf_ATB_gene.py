import pandas as pd
import csv
import io,sys

# read the datasets


df1_extracted_mutations = pd.read_csv(sys.argv[1], delimiter="\t")
antibio=pd.read_csv(sys.argv[2],delimiter="\t")

# print the datasets

d_drug={'Amikacin':'AMK' , 'Bedaquiline':'BDQ', 'Capreomycin' :'CAP', 'Clofazimine' :'CFZ', 'Delamanid':'DLM', 'Ethambutol':'EMB', 'Ethionamide':'ETO' , 'Isoniazid':'INH' , 'Kanamycin':'KAN', 'Levofloxacin' : 'LFX' , 'Linezolid':'LNZ', 'Moxifloxacin':'MFX', 'Pyrazinamide':'PZA', 'Rifampicin':'RIF' , 'Streptomycin':'STR' , 'Fluoroquinolones' : 'FQ' , 'Ciprofloxacin' : 'CIP', 'Ofloxacin' :'OFX' , 'Para-aminosalicylic': 'PAS' , 'Aminoglycosides' : 'AMG',  'Cycloserine' : 'CYS', 'Para-aminosalicylic_acid':'PAS' }
d_drug_complet={'AMK':'Amikacin', 'BDQ':'Bedaquiline',  'CAP':'Capreomycin',  'CFZ':'Clofazimine',  'DLM':'Delamanid',  'EMB':'Ethambutol', 'ETO':'Ethionamide',  'INH':'Isoniazid',  'KAN':'Kanamycin',  'LFX':'Levofloxacin',  'LNZ':'Linezolid',  'MFX':'Moxifloxacin',  'PZA':'Pyrazinamide',  'RIF':'Rifampicin', 'STR':'Streptomycin',  'FQ':'Fluoroquinolones' , 'CIP':'Ciprofloxacin', 'OFX':'Ofloxacin',  'PAS':'Para-aminosalicylic' , 'AMG':'Aminoglycosides',  'CYS':'Cycloserine', 'PAS':'Para-aminosalicylic_acid' , 'CYC':'Cycloserine' ,'SQ109' : 'SQ109'  }


antibio['completName']= antibio['Antibio'].apply(lambda x:d_drug_complet[x])
concat_data = pd.merge(df1_extracted_mutations,antibio)
#x[18] is completName
#x[10] is drug

print ("hereeeee")
concat_data['keep']=concat_data.apply(lambda x:  ((x[18] in (str(x[10]).split(','))) or pd.isna(x[10])) , axis=1)


filter=concat_data[concat_data['keep']].drop(columns="keep")
filter.to_csv(sys.argv[3], sep='\t', index=False)


