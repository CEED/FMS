#!/bin/bash

if [[ ! -e FMS_test_data ]] ; then
    mkdir FMS_test_data
fi

for PROTOCOL in ascii yaml json hdf5; do
    ORDER=1
    for NAME in one two three ; do
        # Use number names in the filename so VisIt will not group them.
        if [[ "$ORDER" == "1" ]] ; then
            ORDERSTR="one"
        elif [[ "$ORDER" == "2" ]] ; then
            ORDERSTR="two"
        elif [[ "$ORDER" == "3" ]] ; then
            ORDERSTR="three"
        fi
        FILENAME="FMS_test_data/domains_${PROTOCOL}_order_${ORDERSTR}.fms"
        FILENAME2="FMS_test_data/domains_${PROTOCOL}_order_${ORDERSTR}.3D"
        echo ./domains $PROTOCOL $ORDER \> $FILENAME2
        ./domains $PROTOCOL $ORDER > $FILENAME2
        mv domains.fms $FILENAME
        ORDER=$((ORDER+1))
    done
done
