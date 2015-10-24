#!/bin/sh

# compute k-mer database for givin k using jellyfish

if [ $# -ne 3 ]; then
    echo "Usage: $0 <target dir> <k> <parallel>" 1>&2
    exit 1
fi

dir=${1}
k=${2}
P=${3}
jellyfish="${HOME}/bin/jellyfish"

find ${dir} -type f -name '*.fasta' | xargs -n1 -P${P} -I@ ${jellyfish} count -m ${k} -s 100M -t 6 @ -o @.${k}.jf
