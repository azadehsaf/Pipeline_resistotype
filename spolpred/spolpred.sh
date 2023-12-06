#/bin/bash

name=$1
shift

echo $@

./SpolPred $@

sed -i "s/^[^\t]*/$name/" /tmp/output.txt

exit 0
