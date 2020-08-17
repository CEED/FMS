#!/bin/bash

if [[ ! -e FMS_test_data ]] ; then
    mkdir FMS_test_data
fi

for PROTOCOL in ascii yaml json hdf5; do
    ORDER=1
    for ORDERSTR in one two three ; do
        FILENAME="FMS_test_data/domains_${PROTOCOL}_order_${ORDERSTR}.fms"
        FILENAME2="FMS_test_data/domains_${PROTOCOL}_order_${ORDERSTR}.3D"
        echo ./domains $PROTOCOL $ORDER \> $FILENAME2
        ./domains $PROTOCOL $ORDER > $FILENAME2
        mv domains.fms $FILENAME
        ORDER=$((ORDER+1))
    done
done

for PROTOCOL in ascii yaml hdf5; do
    ORDER=1
    for ORDERSTR in one two three four; do
        FILENAME="FMS_test_data/quads_${PROTOCOL}_order_${ORDERSTR}.fms"
        FILENAME2="FMS_test_data/quads_${PROTOCOL}_order_${ORDERSTR}.3D"
        echo ./quads $PROTOCOL $ORDER \> $FILENAME2
        ./quads $PROTOCOL $ORDER > $FILENAME2
        mv quads.fms $FILENAME
        ORDER=$((ORDER+1))
    done
done

for PROTOCOL in ascii hdf5; do
    ORDER=1
    for ORDERSTR in one two three four five; do
        FILENAME="FMS_test_data/hex_${PROTOCOL}_order_${ORDERSTR}.fms"
        FILENAME2="FMS_test_data/hex_${PROTOCOL}_order_${ORDERSTR}.3D"
        FILENAME3="FMS_test_data/hex_${PROTOCOL}_order_${ORDERSTR}.lines"
        echo ./hexes $PROTOCOL $ORDER 6 6 6
        ./hexes $PROTOCOL $ORDER 6 6 6
        mv hex.fms $FILENAME
        mv hex.3D $FILENAME2
        mv hex.lines $FILENAME3
        ORDER=$((ORDER+1))
    done
done
