#!/bin/bash
SDIR="$( cd "$( dirname "$0" )" && pwd )"
GENOME=$1
MAPPING=$2

if [ "$#" -ne 2 ]; then
    echo "usage runPEMapperMultiDirectories.sh GENOME MAPPING_FILE"
    exit
fi

for sample in $(cat $MAPPING | cut -f2 | sort | uniq); do
    echo $sample;
    cat $MAPPING | awk -v S=$sample '$2==S{print $4}' \
        | xargs $SDIR/pipe.sh -s $sample $GENOME
done
