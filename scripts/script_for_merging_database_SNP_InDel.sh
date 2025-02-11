  extractedFiledDB=$1 #$1=#{output.extractedField_snpEffFreebayes_DB}
  extractedFieldVartype=$2  #$2=#{output.extractedField_varType_DB}

  ficheDBallReformated=$3 #({output.all_DB_reformated})
  ficheVartypeallReformated=$4 #({output.vartype_DB_reformated})
  mergeBothFiche=$5 #({output.merge_db})

        cat $extractedFiledDB | awk -F'\t' 'BEGIN{{ print "POS\tvarType\tGEN[*].DP\tGEN[*].RO\tGEN[*].AO\tANN[*].GENEID\tANN[*].GENE\tANN[*].HGVS_P\tANN[*].HGVS_C" }}; {{if ($10 ~ /del/){{ print $1"\tDEL\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10 }} else {{if (($10 ~ /ins/) || ($10 ~ /dup/)) {{print $1"\tINS\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10 }}}}}}' > $ficheDBallReformated
        cat $extractedFieldVartype | awk -F'\t' 'BEGIN{{ print "POS\tREF\tALT\tvarType\tP_INTER\tW_INTER"}};{{if (($4 ~ /DEL/)|| ($4 ~ /INS/))  {{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$7 }}}}'   > $ficheVartypeallReformated
         
        python scripts/join_bed_vcf.py  $ficheVartypeallReformated $ficheDBallReformated $mergeBothFiche
        len_merge=`cat mergeBothFiche |wc -l`
        if [ "$len_merge" -gt "1" ];
        then
            sed  -i '/del/d' $extractedFiledDB 
            sed  -i '/ins/d' $extractedFiledDB
            cat $mergeBothFiche | sed 1d | awk -F'\t'  '{{print $1"\t"$2"\t"$3"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$5"\t"$6}}' >> $extractedFiledDB
        fi
