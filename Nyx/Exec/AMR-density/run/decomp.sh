. ../base.sh
./ss comp-005.h5
mv ss.raw comp-005.raw

./ss nocomp.h5
mv ss.raw gold.raw

. ../sz2.sh
./des sz2-001.h5
mv test.raw sz2-001.raw

. ../sz3.sh
./stack sz3-001.h5
mv test.raw sz3-001.raw
