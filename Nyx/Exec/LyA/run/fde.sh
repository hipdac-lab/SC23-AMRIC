# . ../64-base.sh
# ./ss comp-005.h5
# mv ss.raw comp-005.raw

./ss-adp nocomp-241.h5
mv ss.raw gold-241.raw

# # . ../64-sz2.sh
# # ./des sz2-001.h5
# # mv test.raw sz2-001.raw

# . ../64-sz3.sh
# ./stack-f sz3-001.h5
# mv test.raw sz3-001.raw

. ../sz3-pad.sh
./pad-adp pad-281-001.h5
mv test.raw pad-281-001.raw

. ../sz3-adp.sh
./pad-adp-c test.h5
mv test.raw pad-offset-241.raw

./pad-adp pad-241.h5
mv test.raw pad-241.raw