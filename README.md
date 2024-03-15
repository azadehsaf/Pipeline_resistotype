Pipeline_resistotype
This repository contains a pipeline (snakemake file) aiming at identifing
all SNP/SNV within a group of antibiotic-resistant genes in Mycobacterium tuberculosis(reference H37rv).

To run the pipeline you need to create a conda environment. The list of package is given by conda_env_routine_resistotype_pipeline_spec_file.txt

Once conda environment is activated then you can run the pipeline :

snakemake -s snakefile_WHO2023_Withcomment --use-conda --core 15

The pipeline activates MTBseq conda environment while running MTBseq rule. 


