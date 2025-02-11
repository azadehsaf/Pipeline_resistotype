


#tbProfilerout_results=outfolder+"/{sample}/tb_profiler_results_files/results/{sample}_tbprofiler_report.results.txt", 
#        big_table_join_sorted_antiBiotics=outfolder+"/{sample}/final_results_assembling_files/{sample}_finale_table_AntiBiotics_noRodundant.txt",
#        final_lineage=outfolder+"/{sample}/final_results_assembling_files/{sample}_Final_lineage.txt",
#        mtbseq_lin=outfolder+"/{sample}/{sample}_MTBseq_final_report.tab",
#        mtbseq_call=outfolder+"/{sample}/{sample}_called_mutation_name.tab",
#        spoligo=outfolder+"/{sample}/spolpred_spoligoTyping/{sample}_spoligoType_DRall",
#        non_coveredGenes_out=outfolder+"/{sample}/coverage/{sample}_coverageResistGenes_selectedNonCovered.txt"

run=$1
sample=$2

coveragefile="$run/$sample/coverage/$sample"_"coverageResistGenes_selectedNonCovered.txt"

spoligofile="$run/$sample/spolpred_spoligoTyping/$sample"_"spoligoType_DRall"

respipe="$run/$sample/final_results_assembling_files/$sample"_"finale_table_AntiBiotics_noRodundant.txt"
restbprof="$run/$sample/$sample"_"tbprofiler_to_merge.txt"
resmtbseq="$run/$sample/$sample"_"mtbseq_to_merge.txt"
alljoinNofilter="$run/$sample/$sample"_"join_table_dailyPipeline_Tbprofiler_MTBseq_NOTFiltered.csv"
allmerge="$run/$sample/$sample"_"join_table_dailyPipeline_Tbprofiler_MTBseq.csv"
reslintbprof="$run/$sample/tb_profiler_results_files/results/$sample"_"tbprofiler_report.results.txt"
MTBClineage="$run/$sample/final_results_assembling_files/$sample"_"Final_lineage.txt"
mtbseqlineage="$run/$sample/$sample"_"MTBseq_final_report.tab"
krakenreport="$run/$sample/$sample"_"kraken2_result_report.report"

python scripts/script_python_merge.py $respipe $restbprof $resmtbseq $alljoinNofilter




rm  $allmerge

cat $alljoinNofilter  | awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$18}' >> $allmerge

        echo -e "\n\nNon covered genes : \n" >> $allmerge
        cat $coveragefile >> $allmerge
     
	echo -e "\n\nKraken2 taxonomy : \n" >>$allmerge
        cat $krakenreport | cut -f1,4,6 | awk -F'\t' '($2=="S"|| $2=="S1"){print $1"\t"$3}'| sort -n -r -k1 | sed -n 1,2p >>$allmerge
	
        echo -e "\n\nSpoligotyping :\n" >> $allmerge
        echo -e "(Spolpred): \n" >> $allmerge

	cat $spoligofile| cut -d'/' -f6 |cut -f1,2 | sed 's/.fastq//g' >>  $allmerge
        echo -e "\n\nTB_profiler lineage" >> $allmerge

	cat $reslintbprof | sed -n  '/^Lineage report/,/Resistance report/p' | grep -v 'report'  | grep -v '^-' | grep -v '^$' | sort -n -k 1,1 | sed 1d | tr ' ' '_' >> $allmerge
        echo -e "\n\n" >> $allmerge


	echo -e "\n\nMTBC_BC_lineage\n10.1128/msphere.00169-23 \n" >> $allmerge
        cat  $MTBClineage >> $allmerge

	echo -e "\n\nMTBseq lineage:\n" >> $allmerge
        cat  $mtbseqlineage >> $allmerge
