#!/bin/sh

# compute k-mer database for givin k using jellyfish

if [ $# -ne 3 ]; then
    echo "Usage: $0 <target dir> <k> <parallel>" 1>&2
    exit 1
fi

dir=${1}
k=${2}
P=${3}
python="${HOME}/bin/python"
jf_to_all="/work2/yt/chromcode/jf_to_all.py"

find ${dir} -type f -name "*.${k}.jf" | xargs -n1 -P${P} -I@ ${python} ${jf_to_all} @ ${k}

