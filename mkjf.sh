#!/bin/sh

# compute k-mer database for givin k using jellyfish

if [ $# -ne 2 ]; then
    echo "Usage: $0 <target dir> <k>" 1>&2
    exit 1
fi

dir=${1}
k=${2}

for file in `find ${dir} -name *.fasta`; do
    if [ ! -e ${file%.fasta}.${k}.jf ]; then
	jellyfish count -m ${k} -s 100M -t 6 ${file} -o ${file%.fasta}.${k}.jf &
	#echo "jellyfish count -m ${k} -s 100M -t 6 ${file} -o ${file%.fasta}.${k}.jf"
    fi
done

wait

