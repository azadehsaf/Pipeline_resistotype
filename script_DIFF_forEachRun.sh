
new_list=$1

while read -r LINE
do

    echo " HERRRRRRRRRRRREEEEEEEEEEEE $LINE"
    A=`diff output_Run_106/$LINE/$LINE"_"join_table_dailyPipeline_Tbprofiler_MTBseq.csv output_Run_106/ALL_CSV_old/$LINE"_"join_table_dailyPipeline_Tbprofiler_MTBseq.csv`

    echo "$A"
    
done<$1

