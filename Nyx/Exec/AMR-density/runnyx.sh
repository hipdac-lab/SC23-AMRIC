source ~/.bashrc
for i in {0..2}
do
    echo "Running iteration $i..."
    echo "0.001" > eb.txt
    echo "sz2"
    source ./sz2.sh
    mpirun -np 8 ./sz2 inputs > otfile/sz2-$i.txt
    mv run/plt00211_multi.h5 run/sz2-001.h5

    echo "0.001" > eb.txt
    echo "sz3"
    source ./sz3.sh
    mpirun -np 8 ./sz3 inputs > otfile/sz3-$i.txt
    mv run/plt00211_multi.h5  run/sz3-001.h5

    echo "0.005" > eb.txt
    source ./base.sh
    echo "nocomp"
    mpirun -np 8 ./nocomp inputs > otfile/nocomp-$i.txt
    mv run/plt00211_multi.h5  run/nocomp.h5

    echo "comp"
    mpirun -np 8 ./comp inputs > otfile/comp-$i.txt
    mv run/plt00211_multi.h5  run/comp-005.h5
done
rm -rf run/*old*
