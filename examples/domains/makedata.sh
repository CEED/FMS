#!/bin/bash

if [[ ! -e FMS_test_data ]] ; then
    mkdir FMS_test_data
fi

for PROTOCOL in ascii yaml json hdf5; do
    ORDER=1
    for NAME in one two three ; do
        FILENAME="FMS_test_data/domains_${PROTOCOL}_order_${ORDER}.fms"
        FILENAME2="FMS_test_data/domains_${PROTOCOL}_order_${ORDER}.3D"
        ./domains $PROTOCOL $ORDER > $FILENAME2
        mv domains.fms $FILENAME
        ORDER=$((ORDER+1))
    done
done
