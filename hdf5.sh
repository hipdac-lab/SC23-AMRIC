wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.12/hdf5-1.12.2/src/hdf5-1.12.2.tar.gz
tar -xvf hdf5-1.12.2.tar.gz
cd hdf5-1.12.2
mkdir install
./configure --enable-parallel --prefix=$AMRIC_HOME/hdf5-1.12.2/install MPI_CC=mpicc MPI_CXX=mpicxx
make -j 16
make install
cd $AMRIC_HOME
export AMRIC_H5_HOME=$AMRIC_HOME/hdf5-1.12.2/install
echo export AMRIC_H5_HOME=$AMRIC_HOME/hdf5-1.12.2/install >> ~/.bashrc
