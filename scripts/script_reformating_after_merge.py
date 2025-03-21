#!/usr/bin/env python
# coding: utf-8

# python scripts/script_reformating_after_merge.py  {input.all_merge_notfiltered} {output.table1} {output.table2} {output.table3} {dossier} {espece} {lignee}



import pandas as pd 
import numpy as np
import csv,re
import os,io,sys


d_drug_french={'AMK':'Amikacin' ,  'BDQ':'Bédaquiline','CAP':'Capréomycin','CFZ':'Clofazimine', 'DLM':'Délamanid','Pa':'Prétomanide', 'EMB':'Ethambutol', 'ETO':'Ethionamide', 'INH':'Isoniazid','KAN':'Kanamycin','LFX':'Levofloxacin','LNZ':'Linézolid','MFX':'Moxifloxacin','PZA':'Pyrazinamide','RIF':'Rifampicin' , 'STR':'Streptomycin', 'FQ':'Fluoroquinolones', 'CIP':'Ciprofloxacin', 'OFX':'Ofloxacin' , 'Para-aminosalicylic':'PAS', 'AMG':'Aminoglycosides',  'CYS':'Cyclosérine', 'Para-aminosalicylic_acid' :'PAS', 'PAS':'PAS' }

d_drug_english={'AMK':'Amikacin' ,  'BDQ':'Bedaquiline','CAP':'Capreomycin','CFZ':'Clofazimine', 'DLM':'Delamanid','Pa':'Pretomanide', 'EMB':'Ethambutol', 'ETO':'Ethionamide', 'INH':'Isoniazid','KAN':'Kanamycin','LFX':'Levofloxacin','LNZ':'Linezolid','MFX':'Moxifloxacin','PZA':'Pyrazinamide','RIF':'Rifampicin' , 'STR':'Streptomycin', 'FQ':'Fluoroquinolones', 'CIP':'Ciprofloxacin', 'OFX':'Ofloxacin', 'AMG':'Aminoglycosides',  'CYS':'Cycloserine', 'PAS':'PAS', 'CYC':'Cycloserine' }
d_drug_french['BDQ']



#df1= pd.read_csv("/home/azadeh/workspace/snakemake_tutorial/snakemake-snakemake-tutorial-data-f17b467/output_Run_57/2412067826KIA_lib/2412067826KIA_lib_join_table_dailyPipeline_Tbprofiler_MTBseq_NOTFiltered.csv", delimiter="\t") #pipeline
#/home/azadeh/workspace/snakemake_tutorial/snakemake-snakemake-tutorial-data-f17b467/script_reformating_graphicInterface_after_merge/2408056426BEL_S48_join_table_dailyPipeline_Tbprofiler_MTBseq_NOTFiltered.csv
#place_to_write=os.path.dirname(sys.argv[1])+"/"
path_to_write, filename = os.path.split(sys.argv[1])
all_table=filename.split('_')[0]+"_"+filename.split('_')[1]
path_table1_toWrite=sys.argv[2]
path_table_toWriteExcel=path_to_write+"/"+ all_table +".xlsx"
path_table2_toWrite=sys.argv[3]
path_table3_toWrite=sys.argv[4]
dossier=sys.argv[5]
espece=sys.argv[6]
lignee=sys.argv[7]



dfDossier = pd.DataFrame.from_dict({'N° Dossier:': [dossier] , 'Espèce:' :[espece.replace('_',' ')] , 'lignée Coll':[lignee.replace('_',' ')] } , orient='index')

with pd.ExcelWriter (path_table_toWriteExcel) as writer:
    dfDossier.to_excel(writer,sheet_name='dossier_lignée', index=True)  



df1=pd.read_csv(sys.argv[1], delimiter="\t")

df1['HGVS_C']=df1.apply(lambda x:  x[2]+str(x[8])+x[3]  if (x[7]=="rrs" or x[7]=="rrl" )  else  x[10], axis=1)
df1=df1.drop(columns=['N_WHOALL_R', 'N_WHOALL_S','GENEID', 'REF', 'ALT'])
def mut_synonym(x):
    return(x[0]==x[-1] and x!="NaN" and x[0] in ['C','D','S','Q','K','I', 'P', 'T','F', 'N','G', 'H', 'L', 'R', 'W', 'A','V','E','Y','M'])

df1['ANTIBIO']=df1['ANTIBIO'].apply(lambda x: d_drug_english[x] if  pd.notnull(x) else "Na")
dfreformat=df1.rename(columns={'ANTIBIO':'Antibiotique' , 'GENE':'Gène' , 'POS':'Position génomique' , 'HGVS_C':'Mutation (nucléotide)', 'HGVS_P':'Mutation (AA 1 lettre)', 'HGVS_P_3L':'Mutation (AA 3 lettres)','Freq_Mutated':'Fréquence','DP': 'Profondeur', 'P_INTER':'PhyResSE', 'W_INTER':'Walker', 'Res_TBprof':'Tbprofiler','Res_MTBseq':'MTBseq','WHOALL_GRADING':'OMS'})



#removing ddla T365A
dfreformat=dfreformat.loc[(dfreformat["Gène"]!="ddlA" ) | (dfreformat["Mutation (AA 1 lettre)"]!="T365A") ]
#removing synonym mutation #keeping only L203L fabG1
dfreformat['Mutation (AA 1 lettre)']=dfreformat.apply( lambda x: "True" if (mut_synonym(str(x[5])) and not (str(x[4])=="fabG1" and str(x[5])=="L203L")) else x[5] , axis=1)
dfreformat=dfreformat.loc[dfreformat['Mutation (AA 1 lettre)']!="True"]

#table_3 renaming all columns 
dfreformat_table=dfreformat[['Antibiotique','Gène', 'Position génomique','Mutation (nucléotide)','Mutation (AA 1 lettre)', 'Mutation (AA 3 lettres)', 'Fréquence', 'Profondeur', 'PhyResSE', 'Walker' , 'Tbprofiler','MTBseq', 'OMS' ]]
dfreformat_table_3=dfreformat_table
dfreformat_table_3['MTBseq']=dfreformat_table_3['MTBseq'].replace("ND" , "Incertain")
dfreformat_table_3['MTBseq']=dfreformat_table_3['MTBseq'].replace("Res" , "Résistant")
dfreformat_table_3['MTBseq']=dfreformat_table_3['MTBseq'].apply(lambda x : "Phylogénie" if re.match(r"phylo \(.*\)", str(x)) or x=="Polymorphism" else x)

dfreformat_table_3['Tbprofiler']=dfreformat_table_3['Tbprofiler'].replace("ND" , "Incertain")
dfreformat_table_3['Tbprofiler']=dfreformat_table_3['Tbprofiler'].replace("Res" , "Résistant")

dfreformat_table_3['Walker']=dfreformat_table_3['Walker'].replace("Resistant" , "Résistant")
dfreformat_table_3['Walker']=dfreformat_table_3['Walker'].replace("Uncharacterised" , "Incertain")
dfreformat_table_3['Walker']=dfreformat_table_3['Walker'].replace("phylo" , "Phylogénie")
#dfreformat_table_3['Walker']=dfreformat_table_3['Walker'].apply(lambda x : "Phylogénie" if x in ["Benign", "phylo" ] else x)



dfreformat_table_3['PhyResSE']=dfreformat_table_3['PhyResSE'].replace("Resistant" , "Résistant")
dfreformat_table_3['PhyResSE']=dfreformat_table_3['PhyResSE'].replace("CONF","Résistant")
dfreformat_table_3['PhyResSE']=dfreformat_table_3['PhyResSE'].apply(lambda x : "Phylogénie" if re.match(r"phylo\_*\(.*\)", str(x)) or (x in ["Polymorphism" , "phylo"]) else x)
dfreformat_table_3['Mutation (nucléotide)']=dfreformat_table_3['Mutation (nucléotide)'].apply(lambda x : str(x).lower())



dfreformat_table_3.to_csv(path_table3_toWrite, sep='\t', index=False)
#path_table_toWriteExcel="test.xlsx"


################

    



with pd.ExcelWriter (path_table_toWriteExcel, mode='a') as writer:
    dfreformat_table_3.to_excel(writer,sheet_name='table_3', index=False)
#dfreformat_table_2
#dfreformat=dfreformat1




def tbprof_res(x) :
    tbprof=x[10]
       
    if tbprof=="Res":
        return ("Résistant" , "tbprofiler")
    else:
        return ("INCONNU", "tbprofiler")

def mtbseq_res(x) :
    mtbseq=x[11]    
    if mtbseq=="Res":
        return ("Résistant" , "mtbseq")
    if re.match(r"phylo \(.*\)", str(mtbseq)) or mtbseq=="Polymorphism" :
        return ("Sensible" , "mtbseq")
    else: #ND or NAN 
        return ("INCONNU" , "mtbseq")
        

def walker_res(x) :
    walker=x[9]    
    if walker=="Resistant":
        return ("Résistant" , "walker")
    if walker=="phylo":
        return ("Sensible" , "walker")
    else:
        return ("INCONNU" , "walker")

def phy_res(x) :
    phyr=str(x[8])
    if phyr in ["Resistant" , "CONF"]:
        return ("Résistant" , "PhyResSE")
    if re.match(r"phylo\_*\(.*\)", str(phyr)) or (phyr in ["Polymorphism" , "phylo"]) :
        return ("Sensible" , "PhyResSE")
    else:
        return ("INCONNU" , "PhyResSE")


#a=test4
#tbprof_res(a)
#mtbseq_res(a)
#walker_res(a)
#phy_res("Resistant")



# In[ ]:


def interpretation_source(x):
    
    oms=x[12]
    #print(oms)
    oms
    match oms: 
        case "1)_Assoc_w_R" | "2)_Assoc_w_R_-_Interim" : 
            return ("Résistant" , "OMS")
            
        case "4)_Not_assoc_w_R_-_Interim" | "5)_Not_assoc_w_R" : 
            return ("Sensible" , "OMS") 
          
        case _:
            tbprofres=tbprof_res(x)
            mtbseqres= mtbseq_res(x)
            phyres=phy_res(x)
            walkerres=walker_res(x)

            if  "Résistant" in tbprofres:
                return tbprofres
            if  "Résistant" in mtbseqres:
                return mtbseqres
            if  "Résistant" in phyres:
                return phyres 
            if  "Résistant" in walkerres:
                return walkerres
            else:
                if "Sensible" in mtbseqres:
                    return mtbseqres
                if "Sensible" in phyres:
                    return phyres
                if "Sensible" in walkerres:
                    return walkerres
               
                else:
                    if oms =="3)_Uncertain_significance" :
                        return ("Indeterminé" , "all") 
                    else:
                        return ("INCONNU" , "all")
                                   
                

#list of gene for printing

lis_of_gene_for_print=['eis', 'rrs', 'atpE' , 'mmpL5' , 'mmpS5' , 'pepQ', 'Rv0678', 'tlyA' , 'ddn', 'fbiA' , 'fbiB' ,'fbiC', 'folC', 'fgd1' ,'fbiD', 'Rv2983','embB' , 'ethA', 'inhA','ahpC', 'dnaA', 'glpK', 'hadA', 'katG', 'mshA', 'ndh', 'Rv0010c', 'Rv1129c', 'Rv1258c', 'Rv2752c', 'gyrA','gyrB', 'rplC' , 'rrl', 'clpC1' , 'panD', 'pncA', 'PPE35', 'rpsA', 'Rv3236c', 'sigE', 'rpoB', 'gid', 'rpsL']


dfreformat=df1.rename(columns={'ANTIBIO':'Antibiotique' , 'GENE':'Gène' , 'POS':'Position génomique' , 'HGVS_C':'Mutation (nucléotide)', 'HGVS_P':'Mutation (AA 1 lettre)', 'HGVS_P_3L':'Mutation (AA 3 lettres)','Freq_Mutated':'Fréquence','DP': 'Profondeur', 'P_INTER':'PhyResSE', 'W_INTER':'Walker', 'Res_TBprof':'Tbprofiler','Res_MTBseq':'MTBseq','WHOALL_GRADING':'OMS'})
dfreformat=dfreformat.loc[(dfreformat["Gène"]!="ddlA" ) | (dfreformat["Mutation (AA 1 lettre)"]!="T365A") ]
dfreformat['Mutation (AA 1 lettre)']=dfreformat.apply( lambda x: "True" if (mut_synonym(str(x[5])) and not (str(x[4])=="fabG1" and str(x[5])=="L203L")) else x[5] , axis=1)
dfreformat=dfreformat.loc[dfreformat['Mutation (AA 1 lettre)']!="True"]


dfreformat['Mutation (nucléotide)']=dfreformat['Mutation (nucléotide)'].apply(lambda x : str(x).lower())
dfreformat['Interprétation']=dfreformat.apply(lambda x: interpretation_source(x)[0],axis=1)
dfreformat['Source']=dfreformat.apply(lambda x: interpretation_source(x)[1],axis=1)
dfreformat_table2=dfreformat.drop(columns=['Mutation (AA 3 lettres)','Profondeur', 'PhyResSE',	'Walker',	'Tbprofiler',	'MTBseq', 'OMS' ])
dfreformat_table2




dfreformat_table2.to_csv(path_table2_toWrite, sep='\t', index=False) #tableau 2
with pd.ExcelWriter (path_table_toWriteExcel, mode='a') as writer :
    dfreformat_table2.to_excel(writer,sheet_name='table_2',  index=False)




table1_step1=dfreformat_table2.drop(columns=['Position génomique', 'Fréquence',	'Source'])
table1_step1_sansNa=table1_step1[~ table1_step1.Antibiotique.isin(['Fluoroquinolones' , 'Na']) ] # removing Fluoroquinolones
table1_step1_sansNa=table1_step1 [ table1_step1.Gène.isin(lis_of_gene_for_print)  ] # removing all genes wich are not conserning resistance



def present_mutation(x):
    if pd.notna(x[2]) :
        return str(x[1])+" "+str(x[2])
    else:
        return str(x[1])+" "+str(x[3])

table1_step1_sansNa['Mutation']=table1_step1_sansNa.apply(lambda x: present_mutation(x) if str(x[4])!="Sensible" else "" , axis=1)
table1_step1_sansNa
#table1_NoSensible=table1_afterNA[table1_afterNA['Interprétation']!="SENSIBLE"]
#table1_NoSensible




table_for_interpretation_temp=table1_step1_sansNa.drop(columns=['Gène','Mutation (AA 1 lettre)', 'Mutation (nucléotide)'])
table_for_interpretation_temp.reset_index(drop=True, inplace=True) #reset all indices 
table_for_interpretation_temp




antibio_diag_sorted={'Rifampicin':'Sensible','Isoniazid':'Sensible','Pyrazinamide':'Sensible','Ethambutol':'Sensible','Streptomycin':'Sensible','Moxifloxacin':'Sensible','Levofloxacin':'Sensible','Bedaquiline':'Sensible','Clofazimine':'Sensible','Linezolid':'Sensible','Cycloserine' :'Sensible','Delamanid':'Sensible','Pretomanide':'Sensible', 'Kanamycin': 'Sensible', 'Amikacin':'Sensible', 'Capreomycin':'Sensible', 'Ethionamide':'Sensible','PAS':'Sensible'}
order=list(antibio_diag_sorted.keys())

def adding_sensible_antibiotiques(x):
    #print(x.Antibiotique)
    for element in order:
        
        if element not in x.Antibiotique.to_list():
            ind=order.index(element)
            print(element)
            print(ind)
            new_row=[element, "Sensible",'']
            print(new_row)
            x.loc[len(x)]=new_row
            

adding_sensible_antibiotiques(table_for_interpretation_temp)

table_for_interpretation_temp


table_for_interpretation_to_order=table_for_interpretation_temp.drop_duplicates()
table_for_interpretation_to_order


table_for_interpretation_to_order['Antibiotique'] = pd.Categorical(table_for_interpretation_to_order['Antibiotique'], categories = order)
table_for_interpretation=table_for_interpretation_to_order.sort_values(by = 'Antibiotique')
table_for_interpretation


def formating_for_print(x):
    #print(order)
    for cle in order:
        if cle in x.Antibiotique.unique() :
            print(cle)
            if "Résistant" in x.loc[x['Antibiotique']==cle].Interprétation.unique():
                x.drop(x[(x['Antibiotique']==cle) & (x['Interprétation']!="Résistant")].index, inplace=True )
            elif "INCONNU" in x.loc[x['Antibiotique']==cle].Interprétation.unique():
                x.drop(x[(x['Antibiotique']==cle) & (x['Interprétation']=="Sensible")].index, inplace=True )
            elif "Indeterminé" in x.loc[x['Antibiotique']==cle].Interprétation.unique(): 
                x.drop(x[(x['Antibiotique']==cle) & (x['Interprétation']=="Sensible")].index, inplace=True )
    print(x.Antibiotique.unique())

formating_for_print(table_for_interpretation)
#table_for_interpretation




#table_for_interpretation.loc[(table_for_interpretation['Antibiotique']=="Cycloserine") & (table_for_interpretation['Interprétation']=="INCONNU")].Mutation
#table1_step1_sansNa_sorted.loc[table1_step1_sansNa_sorted[table1_step1_sansNa_sorted.Antibiotique=="Streptomycin"].Interprétation=="RESISTANT"]




table1_resistance_Mutaggregate=table_for_interpretation.groupby(['Antibiotique']).agg(lambda x: '; '.join(x))
table1_resistance_Mutaggregate




#table1_resistance_Mutaggregate.columns
table1_to_print =table1_resistance_Mutaggregate.rename(columns={'Interprétation':'Résistance génotypique' , 'Mutation' : 'Mutations' })
table1_to_print.to_csv(path_table1_toWrite, sep='\t', index=True)
with pd.ExcelWriter (path_table_toWriteExcel, mode='a') as writer :
    table1_to_print.to_excel(writer,sheet_name='table_1')


