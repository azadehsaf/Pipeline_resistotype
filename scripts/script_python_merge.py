
import pandas as pd
import numpy as np
import csv,re
import os,io,sys

#python script_python_merge.py 2310061520DEL_S13_finale_table_AntiBiotics_noRodundant.txt  table_tbProfiler_a_merger.txt mtbseq_mut_table.txt current_folder/all_merge.csv


place_to_write=os.path.dirname(sys.argv[4])+"/"
#excel_table_to_write=sys.argv[4].replace(".csv",".xlsx")

print(place_to_write)
#pipeline
df1= pd.read_csv(sys.argv[1], delimiter="\t") #pipeline
df1_tb= pd.read_csv(sys.argv[2], delimiter="\t") # TBprofiler
df1_mtbseqBFilt= pd.read_csv(sys.argv[3], delimiter="\t") #MTBseq



# renaming columns TB


df1_tb_renamed=df1_tb.rename(columns={'Genome_Position': 'POS', 'Locus_Tag' : 'GENEID', 'Gene':'GENE' , 'Change':'EVENT_TB', 'Estimated_Fraction':'Freq_Mutated_TB', 'Drug':'Res_TBprof'})


#renaming columns MTBSEQ
#removing consequensive deleted bases

if len (df1_mtbseqBFilt.loc[df1_mtbseqBFilt.Allel=="GAP"]) > 1 :
    bases_del="".join(df1_mtbseqBFilt.loc[np.where( (df1_mtbseqBFilt.Pos.diff(1) == 1) & (df1_mtbseqBFilt.Allel=="GAP") )]['Ref'])
    if len(bases_del) > 1 :
        start_del=(df1_mtbseqBFilt.loc[np.where((df1_mtbseqBFilt.Pos.diff(1) == 1) & (df1_mtbseqBFilt.Allel=="GAP"))]['Pos']).iloc[0]-1
        df1_mtbseq=df1_mtbseqBFilt.loc[np.where((df1_mtbseqBFilt.Pos.diff(1) != 1) & (df1_mtbseqBFilt.Allel=="GAP") | (df1_mtbseqBFilt.Allel!="GAP"))]
        bases_to_replace=" ".join(df1_mtbseq[df1_mtbseq['Pos']==start_del]['Ref'])+bases_del
        df1_mtbseq.loc[df1_mtbseq['Pos']==start_del,'Ref']=bases_to_replace
    else:
        df1_mtbseq=df1_mtbseqBFilt
else:
    df1_mtbseq=df1_mtbseqBFilt
    

df1_mtbseq.to_csv(place_to_write+"mtbseq_filtered.csv", sep='\t', index=False)


df1_mtbseq_renamed=df1_mtbseq.rename(columns={'Pos':'POS', 'Ref':'REF', 'Allel':'ALT', 'Freq':'Freq_Mutated_MTBSEQ', 'Cov':'DP_MTBSEQ', 'Subst':'EVENT_MTBSEQ', 'Gene':'GENEID', 'GeneName':'GENE', 'Product': 'Product_MTBSEQ','ResistanceSNP':'Res_MTBseq',  'PhyloSNP':'PhyloSNP_MTBSEQ',  'InterestingRegion' : 'ANTIBIO_mtbseq'})

df1_mtbseq.to_csv(place_to_write+"mtbseq_filtered_after_rename.csv", sep='\t', index=False)

d_AA = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K', 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

AA_abreviation= { 'C':'CYS', 'D':'ASP', 'S':'SER',  'Q' : 'GLN',  'K' : 'LYS',  'I' : 'ILE', 'P' : 'PRO', 'T' : 'THR', 'F' : 'PHE',  'N' : 'ASN', 'G' : 'GLY', 'H' : 'HIS', 'L' : 'LEU', 'R' : 'ARG',  'W' : 'TRP',  'A' : 'ALA',  'V' : 'VAL', 'E' : 'GLU',  'Y' : 'TYR', 'M' : 'MET', 'X' : 'Undetermined'}




d_drug={'Amikacin':'AMK' , 'Bedaquiline':'BDQ', 'Capreomycin' :'CAP', 'Clofazimine' :'CFZ', 'Delamanid':'DLM', 'Ethambutol':'EMB', 'Ethionamide':'ETO' , 'Isoniazid':'INH' , 'Kanamycin':'KAN', 'Levofloxacin' : 'LFX' , 'Linezolid':'LNZ', 'Moxifloxacin':'MFX', 'Pyrazinamide':'PZA', 'Rifampicin':'RIF' , 'Streptomycin':'STR' , 'Fluoroquinolones' : 'FQ' , 'Ciprofloxacin' : 'CIP', 'Ofloxacin' :'OFX' , 'Para-aminosalicylic': 'PAS' , 'Aminoglycosides' : 'AMG',  'Cycloserine' : 'CYS', 'Para-aminosalicylic_acid' :'PAS' }
 


def shorten(x):
    if len(str(x)) % 3 != 0: 
        raise ValueError('Input length should be a multiple of three')

    y = ''
    for i in range(len(x) // 3):
        y += d[x[3 * i : 3 * i + 3]]
    return y



#translate three letters AA to one 
def using_d_AA_to_translate(x):
    #print(x)
    y=''
    X=str(x).replace('p.','').upper()

    if X[0]=="_" and X[len(X)-1]=='_':
        return X

    if X[0]=="_" :
        y=X[1:len(X)-3]+d_AA[X[-3:len(X)]]
        return y

    first=d_AA[X[0:3]]   
    if X[len(X)-1]=="*" or  X[len(X)-1]=="_":
        
       y+=first+X[3:len(X)]
    else:
        second=d_AA[X[-3:len(X)]]
        y+=first+X[3:len(X)-3]+second
    return y


def using_d_Abreviation_to_AA_translate(x):
    #print("now")
    #print(x)
    #print (x[len(x)-5:len(x)])
    if x[0]=="_" and x[len(x)-1]=="_":
        return (x)
    if "ext" in x :
        return (x)
    if "fs" in x:
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
    


#split to take protein codon in MTBseq Event column
def col_split(x):
    return x.split(' ')[0]

#converting this c.309C>T to this  C309T
#Chromosome:g.2714130_2714450del

def change_event_nucleo_tbProfiler(x):
    if 'c.' in str(x) :
        x_c_removed=str(x).replace('c.','')
    if 'Chromosome:g.' in str(x) :
        x_c_removed=str(x).replace('Chromosome:g.','')
    if "del" in x_c_removed :
        print("herreeeeeeeeeeeeeeeeeeeeeeEVENT TP ")
        print(x)
        print(x_c_removed)
        return x
    
    index_car=x_c_removed.find('>')
    res=x_c_removed[index_car-1]+x_c_removed[0:index_car-1]+x_c_removed[index_car+1]
    return res



def  change_drug_to_abbreviation(x):
    #print(str(x))
    if x in ["Nan", "Na", '']:
       return ("NA",)
   
    if len(x.split(','))  > 1:
        table=x
        
        a=()
        for t in table.split(','):
           
            a=a+(d_drug[t.capitalize()],)
          
        return a
    #print(str(x))
    return (d_drug[str(x).capitalize()],)


def fill_freq_mutated(x):
    if (str(x[5]) in  ["NA" ,"" ,"nan"]) and (not(str(x[19]) in ["", "nan", "NA"])):
        return (str(float(x[19])*100)+"%").replace(',', '.')
    elif (str(x[5]) in ["NA" , " " , "nan" ]) and (not (str(x[25]) in ["NA" , " " , "nan"])):
        return ((str(x[25])+"%").replace(',', '.')) 
    return (str(x[5]).replace(',', '.') )                                                                 

def fill_depth_mutated(x):
    if ( str(x[4]) in[ "NA" , " " , "nan"]) and (not (str(x[26])in ["NA" ," " ,"nan"])):
        return (x[26])
    return (x[4])

def fill_HGVS_P(x) :
    if (str(x[8])in [ "NA" , " " , "nan"]) and (not (str(x[22]) in ["NA" ," " ,"nan"])):
         return (str(x[22]))
    elif  (str(x[8])in [ "NA" , " " , "nan"])  and  (not (str(x[27]) in ["NA" ," " ,"nan"])) :
        return (str(x[27]))
    return (x[8])


def fill_HGVS_C(x) :
    if (str(x[9])in [ "NA" , " " , "nan"]) and (not (str(x[21]) in ["NA" ," " ,"nan"])):
         return (str(x[21]))
    return (x[9])



def fill_antibio(x):
    if (str(x[0]) in ["NA" ,"" ,"nan"])  and  (not(str(x[23]) in ["", "nan", "NA"])) :
        return str(x[23])
    elif (str(x[0]) in ["NA" ,"" ,"nan"])  and  (not(str(x[31]) in ["", "nan", "NA"])):
        return str(x[31])
    return "NA"




def res_phylo(x) :
    #resistance or phylo, if resistance is not identified then return phylogeny 
    if  str(x[11])!="NA" and str(x[11]) !=" " and str(x[11])!= "nan":
        return(x[11])
    return(x[10])

    
#TB_profiler table reconfiguration
df1_tb_renamed['GENEID']=df1_tb_renamed['GENEID'].apply(lambda x: "EBG00000313339" if (x=="rrl") else x)
df1_tb_renamed['GENEID']=df1_tb_renamed['GENEID'].apply(lambda x: "EBG00000313325" if (x=="rrs") else x)
df1_tb_renamed['HGVS_C_TB']=df1_tb_renamed['EVENT_TB'].apply(lambda x: change_event_nucleo_tbProfiler(x) if ("c." in str(x))  or ("Chromosome:g." in str(x))   else "NA")
df1_tb_renamed['HGVS_P_TB']=df1_tb_renamed['EVENT_TB'].apply(lambda x: x if ("p." in str(x) ) else "NA")
df1_tb_renamed['HGVS_P_TB']=df1_tb_renamed['HGVS_P_TB'].apply(lambda x: using_d_AA_to_translate(x) if str(x) !="NA" else x)
df1_tb_renamed['ANTIBIO_TB']=df1_tb_renamed['Res_TBprof'].apply(lambda x: change_drug_to_abbreviation(str(x).capitalize())[0] if not  x == "NaN" else "NA")


#translating AA three letters to one letter and removing p.
df1_mtbseq_renamed['GENE']=df1_mtbseq_renamed.apply(lambda x:  x[7] if x[8]=="-" else x[8], axis=1 )
#print(df1_mtbseq_renamed['GENE'])
df1_mtbseq_renamed['EVENT_MTBSEQ']=df1_mtbseq_renamed['EVENT_MTBSEQ'].transform(col_split)
df1_mtbseq_renamed['EVENT_MTBSEQ']=df1_mtbseq_renamed['EVENT_MTBSEQ'].apply(lambda x: using_d_AA_to_translate(x) if x !="" else "NA")
df1_mtbseq_renamed['ANTIBIO_mtbseq']=df1_mtbseq_renamed['ANTIBIO_mtbseq'].apply(lambda x: change_drug_to_abbreviation(col_split(x).capitalize())[0]  if x !=" " else "NA")
df1_mtbseq_renamed['Res_MTBseq']=df1_mtbseq_renamed.apply(lambda x: res_phylo(x) if (str(x)!="NA") else "NA" , axis=1)
 


#df1
df1['HGVS_C']=df1['HGVS_C'].apply(lambda x: str(x).replace('c.','').replace('n.','') if ("c." in x ) else "NA")
df1['HGVS_P']=df1['HGVS_P'].apply(lambda x: str(x).replace('p.','') if (x!=""or x!="nan") else "NA")



#dfmerged_df1_tbprof=df1.merge(df1_tb_renamed, how="outer", on=["POS","GENEID","GENE","ANTIBIO"])
dfmerged_df1_tbprof=df1.merge(df1_tb_renamed, how="outer", on=["POS","GENE", "GENEID"])
dfmerged_df1_tbprof.to_csv(place_to_write+"dfmerged_df1_tbprof_mtbseq_firsTest.csv", sep='\t', index=False)
dfmerged_df1_tbprof_mtbseq=dfmerged_df1_tbprof.merge(df1_mtbseq_renamed, how="outer", on=["POS","GENE","GENEID","REF","ALT"]) #,"ANTIBIO" 
dfmerged_df1_tbprof_mtbseq.to_csv(place_to_write+"dfmerged_df1_tbprof_mtbseq_all_test.csv", sep='\t', index=False)

###############
dfmerged_df1_tbprof_mtbseq['Freq_Mutated']=dfmerged_df1_tbprof_mtbseq.apply(lambda x: fill_freq_mutated(x) if(str(x)!="NA" and str(x)!=" " and str(x)!="nan") else "NA", axis=1)
dfmerged_df1_tbprof_mtbseq['DP']=dfmerged_df1_tbprof_mtbseq.apply(lambda x:  fill_depth_mutated(x) if(str(x)!="NA" and str(x)!=" " and str(x)!="nan") else "NA", axis=1)

dfmerged_df1_tbprof_mtbseq['HGVS_P']=dfmerged_df1_tbprof_mtbseq.apply(lambda x: fill_HGVS_P(x) if(str(x)!="NA" and str(x)!=" " and str(x)!="nan") else "NA", axis=1)
dfmerged_df1_tbprof_mtbseq['HGVS_C']=dfmerged_df1_tbprof_mtbseq.apply(lambda x: fill_HGVS_C(x) if(str(x)!="NA" and str(x)!=" " and str(x)!="nan") else "NA", axis=1)
dfmerged_df1_tbprof_mtbseq.to_csv(place_to_write+"dfmerged_df1_tbprof_mtbseq_all_test_motre.csv", sep='\t', index=False)

dfmerged_df1_tbprof_mtbseq['ANTIBIO']=dfmerged_df1_tbprof_mtbseq.apply(lambda x: fill_antibio(x) if((not (str(x) in ["NA" , " ","nan"])) and (str(x[0]) in ["NA" ,"" ,"nan", "Na"])) else str(x[0]) , axis=1)
#merge_after_drop=dfmerged_df1_tbprof_mtbseq.drop(columns=['EVENT', 'PROG', 'Type', 'Freq_Mutated_MTBSEQ', 'DP_MTBSEQ', 'EVENT_MTBSEQ', 'Product_MTBSEQ','PhyloSNP_MTBSEQ', 'Freq_Mutated_TB'])
merge_after_drop=dfmerged_df1_tbprof_mtbseq.drop(columns=['EVENT', 'PROG', 'Type', 'Freq_Mutated_MTBSEQ', 'DP_MTBSEQ', 'Product_MTBSEQ','PhyloSNP_MTBSEQ', 'Freq_Mutated_TB','ANTIBIO_mtbseq'])
 #'EVENT_MTBSEQ'
#print("here")

merge_after_drop['DRUG']=merge_after_drop['DRUG'].apply(lambda x: change_drug_to_abbreviation(x) if (str(x)!="NA" and str(x) != "nan" ) else "Na")

merge_after_drop.to_csv(place_to_write+"dfmerged_df1_tbprof_mtbseq_all_afterDroping.csv", sep='\t', index=False)


def validate_antibio_grading(x):
    a=x[10]
  
    if x[0]=='FQ' or a=="Na":
        return "NA"
    if len(a)>1:
        ind,val=0,""
        for ind,val in enumerate(a, start=1):
            if val==x[0]:
                return ind
        return "NA"
    else:
        return 0

def validate_antibio_number_sensible(x):
    ind=x[23]
    a=x[12]

    if ind in ['NA','NAN','nan','None','NaN',"", '', ' '] :
        return "NA"
    if ind > len(str(a).split(',')):
        return "NA"
    if len(str(a).split(','))>1:
        return str(a).split(',')[int(ind)-1].replace('(','').replace(')','')
    else:
        return str(x[12])
  
    

def validate_antibio_number_resistance(x):
    ind=x[23]
    r=x[11]
    if ind in ['NA','NAN','nan','NaN','None',"", '']:
        return "NA"
    if ind > len(str(r).split(',')):
        return "NA"
    if len(str(r).split(','))>1:
        return str(r).split(',')[int(ind)-1].replace('(','').replace(')','')
    else:   
        return str(x[11])



def validate_antibio_all_grading(x):
    ind=x[23]
    g=x[13]
    if str(ind) in  ['NA','NAN','nan','None',"", '']:
        return "NA"
    if ind > len(str(g).split(',')):
        return "NA"
    if len(str(g).split(','))>1:
        return str(g).split(',')[int(ind)-1]
    else:   
        return str(x[13])
    

    
def validate_p_inter(x):
    ind=x[23]
    g=x[14]
    if str(ind) in ['NA','NAN','nan','None',"", '']:
       return str(g)
    if g=="NA":
        return "NA"
    if re.match(r"phylo_\(.*\)", str(g)):
        return g
    if ind > 0 and len(str(g).split(','))>1 and ind <= len(str(g).split(',')):
        return str(g).split(',')[int(ind)-1]
    if ind > len(str(g).split(',')):
        return str(g).split(',')[0] 
    else:   
        return str(g)


def validate_w_inter(x):
   
    ind=x[23]
    g=x[15]
   
    if str(ind) in ["NA","NAN","nan","None","", " "]:
        return str(g)
    if g=="NA":
        return "NA"
    if ind >0 and len(str(g).split(','))>1 and ind <= len(str(g).split(',')):
        return str(g).split(',')[int(ind)-1]
    if ind > len(str(g).split(',')):
        return str(g).split(',')[0]   
    else:   
        return str(g)
    

def validate_tbProf_resist(x):
    if str(x[16])=="nan":
        return " "
    elif (str(x[17])=="nan" and str(x[19])==str(x[8])) or ( str(x[17])=="nan" and str(x[18])==str(x[9])):
        return "ND"
    elif (x[0] in change_drug_to_abbreviation(str(x[17]).capitalize())) or (str(x[17])!="nan" and str(x[0])=="nan") :
        print(change_drug_to_abbreviation(str(x[17]).capitalize()))
        return "Res"
    elif str(x[16])!="nan":
        return "ND"
    else:
        return " "


def validate_mtbseq_resist(x):
    if str(x[21])==" " or str(x[21]).capitalize()=="NAN":
        return " " 
    if  (str(x[21])==str(x[8]) and str(x[22])==" "):
        return "ND"
    if "phylo" in str(x[22]):
        return str(x[22])
    elif (change_drug_to_abbreviation(str(x[22]).split(' ')[0].capitalize())[0]==x[0]) or (str(x[22])!="nan" and str(x[22])!=" " and str(x[0])=="nan"):
        return "Res"
    elif str(x[21])=="NA" :
        return "ND"
    else:
        return " "
    




merge_after_drop['INDEXING']=merge_after_drop.apply(lambda x: validate_antibio_grading(x) if (str(x)!="NA") else "NA" , axis=1)
merge_after_drop['N_WHOALL_S']=merge_after_drop.apply(lambda x: validate_antibio_number_sensible(x) if (str(x)!="NaN") else "NA" , axis=1)
merge_after_drop['N_WHOALL_R']=merge_after_drop.apply(lambda x: validate_antibio_number_resistance(x) if (str(x)!="NA") else "NA" , axis=1)

merge_after_drop['WHOALL_GRADING']=merge_after_drop.apply(lambda x: validate_antibio_all_grading(x) if (str(x)!="NA") else "NA" , axis=1)
merge_after_drop['P_INTER']=merge_after_drop.apply(lambda x:  validate_p_inter(x) if (str(x)!="NA") else "NA" , axis=1)
merge_after_drop['W_INTER']=merge_after_drop.apply(lambda x:  validate_w_inter(x) if (str(x)!="NA") else "NA" , axis=1)
merge_after_drop['Res_TBprof']=merge_after_drop.apply(lambda x: validate_tbProf_resist(x) if (str(x)!="NA")  else "NA" , axis=1)
merge_after_drop['Res_MTBseq']=merge_after_drop.apply(lambda x: validate_mtbseq_resist(x) if (str(x)!="NA") else "NA" , axis=1)
merge_after_drop['HGVS_P_3L']=merge_after_drop['HGVS_P'].apply(lambda x: using_d_Abreviation_to_AA_translate(x) if (str(x)!=" "and  str(x)!="nan") else "NA")

merge_after_drop=merge_after_drop.drop(columns=['INDEXING', 'DRUG' , 'HGVS_C_TB' , 'HGVS_P_TB' , 'EVENT_TB', 'ANTIBIO_TB','EVENT_MTBSEQ' ])

#merge_after_drop.to_csv("before_antibio.csv", sep='\t', index=False)


merge_after_drop=merge_after_drop.loc[(merge_after_drop['ANTIBIO']!='NA')]
merge_after_dropAntibioNa=merge_after_drop[merge_after_drop['ANTIBIO'].notnull()]
#merge_after_dropAntibioNa=merge_after_drop


merge_after_dropAntibioNa = merge_after_dropAntibioNa[['ANTIBIO', 'POS', 'REF', 'ALT', 'DP', 'Freq_Mutated', 'GENEID', 'GENE' ,'HGVS_P','HGVS_P_3L','HGVS_C', 'P_INTER', 'W_INTER','Res_TBprof','Res_MTBseq', 'N_WHOALL_R','N_WHOALL_S','WHOALL_GRADING']].replace("nan", '', regex=True).replace("NA", '', regex=True).sort_values(by=['ANTIBIO','POS'])
#merge_after_dropAntibioNa = merge_after_dropAntibioNa[['ANTIBIO', 'POS', 'REF', 'ALT', 'DP', 'Freq_Mutated', 'GENEID', 'GENE' ,'HGVS_P','HGVS_P_3L','HGVS_C', 'P_INTER', 'W_INTER','Res_TBprof','Res_MTBseq', 'N_WHOALL_R','N_WHOALL_S','WHOALL_GRADING']].replace("nan", '', regex=True).replace("NA", '', regex=True).sort_values(by=['ANTIBIO','POS'])

merge_after_dropAntibioNa.to_csv(sys.argv[4], sep='\t', index=False)
#merge_after_dropAntibioNa.to_excel(excel_table_to_write)



#kraken2 --db /home/admin-rac/database/kraken2db  -threads 8 --report res_kraken2_seq.report --gzip-compressed mtbseq_results_files/2407065649LAVD_S48_R1.fastq.gz mtbseq_results_files/2407065649LAVD_S48_R2.fastq.gz 

