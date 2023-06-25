


echo "---------- Data Quality for original AMReX Compression ----------"

./compareData -d gold.raw comp-005.raw

echo "------------------------ END ------------------------"
echo " "

echo "---------- Data Quality for AMRIC-SZ_L/R ----------"

./compareData -d gold.raw sz2-001.raw

echo "------------------------ END ------------------------"
echo " "

echo "---------- Data Quality for AMRIC-SZInterp ----------"

./compareData -d gold.raw sz3-001.raw

echo "------------------------ END ------------------------"
echo " "

echo "---------- CR for original AMReX Compression ----------"

python3 cr.py nocomp.h5 comp-005.h5

#echo "------------------------ END ------------------------"
echo " "

echo "---------- CR for AMRIC-SZ_L/R ----------"

python3 cr.py nocomp.h5 sz2-001.h5

#echo "------------------------ END ------------------------"
echo " "

echo "---------- CR for AMRIC-SZInterp ----------"

python3 cr.py nocomp.h5 sz3-001.h5
echo " "

#echo "------------------------ END ------------------------"


