cd $AMRIC_HOME
git clone https://github.com/FabioGrosso/AMRIC-AMReX.git amrex

cd $AMRIC_HOME
git clone https://github.com/FabioGrosso/AMRIC-Nyx.git Nyx
cd Nyx/subprojects
. build-sun.sh

cd $AMRIC_HOME/Nyx/Exec/AMR-density
make clean
cp $AMRIC_HOME/amrex/Src/Extern/HDF5/AMReX_PlotFileUtilHDF5.nast $AMRIC_HOME/amrex/Src/Extern/HDF5/AMReX_PlotFileUtilHDF5.cpp
make -j 16 USE_HDF5_SZ3=TRUE SZ3_HOME=$SZ_SLE_HOME
mv Nyx3d* sz2

make clean
cp $AMRIC_HOME/amrex/Src/Extern/HDF5/AMReX_PlotFileUtilHDF5.stack $AMRIC_HOME/amrex/Src/Extern/HDF5/AMReX_PlotFileUtilHDF5.cpp
make -j 16 USE_HDF5_SZ3=TRUE SZ3_HOME=$SZ3_HOME
mv Nyx3d* sz3

make clean
cp $AMRIC_HOME/amrex/Src/Extern/HDF5/AMReX_PlotFileUtilHDF5.ori $AMRIC_HOME/amrex/Src/Extern/HDF5/AMReX_PlotFileUtilHDF5.cpp
make -j 16 USE_HDF5_SZ3=FALSE 
mv Nyx3d* nocomp

make clean
cp $AMRIC_HOME/amrex/Src/Extern/HDF5/AMReX_PlotFileUtilHDF5.ori $AMRIC_HOME/amrex/Src/Extern/HDF5/AMReX_PlotFileUtilHDF5.cpp
make -j 16 USE_HDF5_SZ3=TRUE SZ3_HOME=$ORI_SZ3_HOME
mv Nyx3d* comp

cd $AMRIC_HOME

