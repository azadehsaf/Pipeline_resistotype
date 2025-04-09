# Pipeline_resistotype
This repository contains a simple pipeline, snakemake file with related scripts and data. The pipeline is a tool to identify all SNV SNP/InDel within a group of antibiotics resistance genes belonging to Mycobacterium tuberculosis (reference H37rv).
This pipeline also allows to identify the taxonomy and lineage of a given MTB strain.
This pipline is useable for short read Mycobacterium Tuberculosis strain in paired_end "Fastq" file given in [A-Z][0-9]_[A-Z][0-9]_R[12].fastq 

#The pipline uses:

-A homemade pipeline wich call for each SNPs and InDels through resistance genes. 
-MTBseq https://github.com/ngs-fzb/MTBseq_source
-TB-profiler 2.8
-kraken 2



# INSTALLATION:

To run the pipeline you need:

 -Python and pandas (2.0.3) 
 
 -install snakemake.

The list of package is given by "spec-file-resistotype-pipeline.txt"
The easiest way to install all needed tools is to clone the file "spec-file-resistotype-pipeline.txt" and then create the conda envs and use conda to install all mandatory software in conda envs.  

conda create --name myenv --file spec-file-resistotype-pipeline.txt
conda install --name myenv --file spec-file.txt

Next you need to install MTBseq, by clonning and using "mtbseq-env.yaml"

conda env create -f mtbseq-env.yaml


You may need to install pdflatex.

Once conda environment is activated then you can run the pipeline 


# Running the Pipelie: 

snakemake -s  snakefile_WHO2023_Withcomment_removing_pindel_modifEtoile_graphicTables --use-conda --core 15
The pipeline activates MTBseq conda environment while running MTBseq rule.

