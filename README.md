# Pipeline_resistotype
This repository contains a simple pipeline, a snakemake file with related scripts and data. The pipeline is a tool to identify all SNV SNP/InDel within a group of antibiotics resistance genes belonging to Mycobacterium tuberculosis reference H37rv .
This pipeline also allows to identify the taxonomy and lineage of a given MTB strain.

# Usage 
This pipline is useable for short paired_end read of Mycobacterium Tuberculosis strain in ".fastq.gz" or "fastq.dsrc" file.

The sequence have to be called as follow:
([A-Za-z]+|[0-9]+)+_([A-Za-z]+|[0-9]+)+_R[12].fastq.gz

# The pipline bases on:

-A homemade pipeline wich call for each SNPs and InDels through resistance genes. 

-MTBseq https://github.com/ngs-fzb/MTBseq_source

-TB-profiler 2.8 https://github.com/jodyphelan/TBProfiler

-Kraken2 https://github.com/DerrickWood/kraken2

# INSTALLATION:

To run the pipeline you need:

 -Python and pandas (2.0.3) 
 
 -install snakemake.

The list of packages is given by "spec-file-resistotype-pipeline.txt"
The easiest way to install all needed tools (including tbProfiler and kraken2) is to clone the file "spec-file-resistotype-pipeline.txt" and then create a new conda env and install all mandatory software in the same conda envs.  

conda create --name myenv --file spec-file-resistotype-pipeline.txt
conda install --name myenv --file spec-file.txt

Next you need to install MTBseq, by clonning and using "mtbseq-env.yaml"

conda env create -f mtbseq-env.yaml

If you use DSRC compression, then before using pipeline please install DSRC in your home or in conda env that you created for running pipeline.  

Once conda environment is activated then you can run the pipeline 


# Running the Pipelie: 

snakemake -s  snakefile_WHO2023_Withcomment_removing_pindel_modifEtoile_graphicTables --use-conda --core 15

The pipeline activates MTBseq conda environment while running MTBseq rule.

# Calling SNVs SNPs InDel:

Raw data were analyzed using in house pipeline (BWA-MEM 0.7.17-r1188, Samtools 1.15.1 mpileup 1.15.1, Freebayes V1.0.2 and Varscan V2.4.6). 
All SNPs and SNVs were annotated by OMS 2023 Catalogue of mutations in MTB.  
For lineage identification, MTBC_BC barcoding method ("10.1128/msphere.00169-23") was used.

Furthermore two more pipelines: TB-profiler 2.8.12 and MTBseq 1.1.0  were used  for cross validation of our findings. Coll lineage reported by both pipelines was used to confirm our lineage report. 

[SNV_calling_pipeline.pdf](https://github.com/user-attachments/files/19707089/SNV_calling_pipeline.pdf)
![alt_text](https://github.com/user-attachments/assets/780383d6-5fc8-42f0-a582-dc119619fcee)


# Resistance Reporting:

For identifying all variants involving in resistance the pipeline uses Walker and PhyresSE database of mutation in addition to OMS 2023 Catalogue of mutations in MTB.
Furthermore Tbprofiler and MTBseq report their own findins. To sum it all up pipeline uses the following steps persuming a variant "Resistance"/"Sensitive"/"Undetermined" or "Inconnu". 

![alt_text](https://github.com/user-attachments/assets/270f36a2-d440-4f1a-b349-6be3cf5921b5)



# Config file .yaml

datas: "/home/admin-rac/Documents/run_87" 

output_folder: "/home/admin-rac/workspace/Pipeline_resistotype-master/output_Run_87"

reference: "/home/admin-rac/workspace/Pipeline_resistotype-master/data"


pair:

 [tab]R1 : R1
 
 [tab]R2 : R2
    
samples:

 [tab]2503005535OUK_lib : 2503005535OUK_lib
 
 [tab]2412067017SHA_lib : 2412067017SHA_lib,
 
 [tab]2412072999BEG_lib : 2412072999BEG_lib,
 
 [tab]2502034815OLL_lib : 2502034815OLL_lib






# Output Files and Format
All result is availabe in "output_folder" given in config.yaml file.

For each sample the is a result folder placed in "output_folder".

The output of SNV calling with taxonomy and lineage is generated in a  excel table "sampleName.xlsx"


The statistics of gene coverage and mapping along with report of deleted fragments is genated  "sampleName_mapping_QC.xlsx   




# Scripts and Data Folder

# Runnin on Computing Cloud: 
Need a config File 

