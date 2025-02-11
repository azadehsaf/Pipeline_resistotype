  extractedFiledDB=$1 #$1=#{output.extractedField_snpEffFreebayes_WHO2023}
  extractedFieldVartype=$2  #$2=#{output.extractedField_varType_WHO2023}

  ficheDBallReformated=$3 #{output.all_who_reformated}
  ficheVartypeallReformated=$4 #{output.vartype_who_reformated}
  mergeBothFiche=$5 #({output.merge_who})

        cat $extractedFiledDB | awk -F'\t' 'BEGIN{{ print "POS\tvarType\tGEN[*].DP\tGEN[*].RO\tGEN[*].AO\tANN[*].GENEID\tANN[*].GENE\tANN[*].HGVS_P\tANN[*].HGVS_C\tEVENT" }}; {{ if ($10 ~ /del/) {{ print $1"\tDEL\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11 }} else {{ if (($10 ~ /ins/) || ($10 ~ /dup/)) {{print $1"\tINS\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11 }}}}}}'  > $ficheDBallReformated
        cat $extractedFieldVartype | awk -F'\t' 'BEGIN{{ print "POS\tREF\tALT\tvarType\tEVENT\tDRUG\tN_WHOALL_R\tN_WHOALL_S\tWHOALL_GRADING" }}; {{if (($4 ~ /DEL/)|| ($4 ~ /INS/)) {{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}}}}'   > $ficheVartypeallReformated
         
        python scripts/join_bed_vcf_vartype_all_WHO.py   $ficheVartypeallReformated $ficheDBallReformated $mergeBothFiche
        len_merge=`cat mergeBothFiche |wc -l`
        if [ "$len_merge" -gt "1" ];
        then
            sed  -i '/del/d' $extractedFiledDB 
            sed  -i '/ins/d' $extractedFiledDB
            cat $mergeBothFiche | sed 1d |  awk -F'\t' '{{print $1"\t"$2"\t"$3"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}}' >> $extractedFiledDB
        fi
