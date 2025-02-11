import pandas as pd
import numpy as np
import csv
import io,sys

# read the datasets


AA_abreviation= { 'C':'Cys', 'D':'Asp', 'S':'Ser',  'Q' : 'Gln',  'K' : 'Lys',  'I' : 'Ile', 'P' : 'Pro', 'T' : 'Thr', 'F' : 'Phe',  'N' : 'Asn', 'G' : 'Gly', 'H' : 'His', 'L' : 'Leu', 'R' : 'Arg',  'W' : 'Trp',  'A' : 'Ala',  'V' : 'Val', 'E' : 'Glu',  'Y' : 'Tyr', 'M' : 'Met', 'X' : 'Undetermined'}

vartype = pd.read_csv(sys.argv[1], delimiter="\t")
all_annotations=pd.read_csv(sys.argv[2],delimiter="\t")


def using_d_Abreviation_to_AA_translate(x):
    print("now")
    print("x")
    print(x)
   
    if x=="Nan" or x=="nan" :
        return "NaN"
    if "ext" in x :
        return(x)
    if  "fs" in x :
        return(AA_abreviation[x[0]]+x[1:len(x)-2]+"fs")
    if  "del" in x  or "ins" in x:
        return(x)
    if x[len(x)-1]=="s" or x[len(x)-1]=="X" or x[len(x)-1]=="*" or x[len(x)-1]=="_" or  x[len(x)-1]=="?"  :
        
        return(AA_abreviation[x[0]]+x[1:len(x)-1])
    elif not any(c.isalpha() for c in x):
        return (x)
    elif not x[0].isalpha():
        return(x[0:len(x)-1]+AA_abreviation[x[len(x)-1]])
    elif "dup" in x :
        return(x)
    
    return(AA_abreviation[x[0]]+x[1:len(x)-1]+AA_abreviation[x[len(x)-1]])
    
dataf=all_annotations
if not dataf.empty:
    

#dataf['EVENT']=np.where(dataf['EVENT'].isna(),  str(dataf['ANN[*].GENE']) + "_p." + using_d_Abreviation_to_AA_translate(str(dataf['ANN[*].HGVS_P']).replace('p.' ,"")), dataf['EVENT'].apply(str).str.split(','))

    print(pd.isnull(dataf['EVENT']))
    dataf['EVENT']=dataf.apply(lambda x: x[9].split(',') if pd.isnull(x[9])==False  else x[6] + "_p." + using_d_Abreviation_to_AA_translate(str(x[7]).replace("p.","")) , axis=1)
    print(dataf)
#dataf['EVENT'] = dataf['EVENT'].apply(str).str.split(',')
    
    
df_exploded = dataf.explode('EVENT', ignore_index=True)

concat_data = pd.merge(vartype,df_exploded)
concat_data.to_csv(sys.argv[3], sep='\t', index=False)


