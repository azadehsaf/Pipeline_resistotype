
new_list=$1
run=$2

while read -r LINE
do

    echo "erreur $LINE"
    ./script_for_regenrating_pipelineResultFile.sh $run $LINE
   
done<$1

