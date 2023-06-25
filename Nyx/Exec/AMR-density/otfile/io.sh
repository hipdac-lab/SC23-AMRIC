#!/bin/bash

echo "---------- Writing Time for No Compression ----------"
for i in {0..2}
do
    echo "***** run $i *****"
    python3 readnocomp.py nocomp-$i.txt
done
echo "------------------------ END ------------------------"
echo " "

echo "---------- Writing Time for original AMReX Compression ----------"
for i in {0..2}
do
    echo "***** run $i *****"
    python3 readAmrex.py comp-$i.txt
done
echo "------------------------ END ------------------------"
echo " "

echo "---------- Writing Time for AMRIC-SZ_L/R ----------"
for i in {0..2}
do
    echo "***** run $i *****"
    python3 readSZLR.py sz2-$i.txt
done
echo "------------------------ END ------------------------"
echo " "

echo "---------- Writing Time for AMRIC-SZInterp ----------"
for i in {0..2}
do
    echo "***** run $i *****"
    python3 readSZInterp.py sz3-$i.txt
done
echo "------------------------ END ------------------------"

