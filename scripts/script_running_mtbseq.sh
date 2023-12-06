gzip -c $1 > $3
gzip -c $2 > $4
cd $5
echo "running MTBseq"
MTBseq --step TBfull --threads $6 --lowfreq_vars 1 --minfreq 2 --mincovf 2 --mincovr 2 --minbqual 13 2> ../mtbseq_all_error
echo "MTBseq finished"
cd 
