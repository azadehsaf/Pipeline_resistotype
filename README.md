# Pipeline_resistotype
This repository contains a simple pipeline, snakemake file with related scripts and data. The pipeline is a tool to identify all SNV SNP/InDel within a group of antibiotics resistance genes belonging to Mycobacterium tuberculosis (reference H37rv).


INSTALLATION:
To run teh pipeline you need to install snakemake.

To run the pipeline you need to create a conda environment. The list of package is given by conda_env_routine_resistotype_pipeline_spec_file.txt
Once conda environment is activated then you can run the pipeline 


Running the Pipelie: 

snakemake -s snakefile_WHO2023_Withcomment --use-conda --core 15
The pipeline activates MTBseq conda environment while running MTBseq rule.

