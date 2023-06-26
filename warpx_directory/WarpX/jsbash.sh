source ~/.bashrc
for i in {0..2}
do
    echo "Running iteration $i..."
    echo "0.005" > eb.txt
    echo "sz2"
    source ./64-sz2.sh
    jsrun -n 8 ./sz2 inputs > otfile/sz2-$i.txt
    mv diags/plt001001.h5 diags/sz2-001.h5

    echo "0.005" > eb.txt
    echo "sz3"
    source ./64-sz3.sh
    jsrun -n 8 ./sz3 inputs > otfile/sz3-$i.txt
    mv diags/plt001001.h5 diags/sz3-001.h5

    echo "0.01" > eb.txt
    source ./64-base.sh
    echo "nocomp"
    jsrun -n 8 ./nocomp inputs > otfile/nocomp-$i.txt
    mv diags/plt001001.h5 diags/nocomp.h5

    echo "comp"
    jsrun -n 8 ./comp inputs > otfile/comp-$i.txt
    mv diags/plt001001.h5 diags/comp-005.h5
done
rm -rf diags/*old*
