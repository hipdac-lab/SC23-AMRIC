# SC23-AMRIC Artifacts

AMRIC is a novel in-situ lossy compression framework that leverages the HDF5 filter to enhance both I/O efficiency and compression quality for Adaptive Mesh Refinement (AMR) applications. AMRIC was implemented into the AMReX framework and evaluated on two real-world AMR applications, Nyx and WarpX.

While preparing the artifacts, we executed them on a single node from the Chameleon Cloud, equipped with two Intel Xeon Gold 6242 CPUs and 192 GB of memory (specifically, ```compute_skylake``` configuration). We recommend that reviewers also use the Chameleon Cloud for artifact evaluation.

## Method 1: Use Singularity Image (Recommended)
The entire workflow takes approximately 10 minutes to execute, including downloading container image and preparing environment (3 mins), running WarpX simulation (3 mins), running Nyx simulation (3 mins), and evaluating compression performance (1 min).

### Minimum system requirements
OS: Ubuntu (20.04 is recommended)

Memory: >= 16 GB RAM

Processor: >= 8 cores

### Step 1: Install Singularity
Install [Singularity](https://singularity-tutorial.github.io/01-installation/)

### Step 2: Download the pre-built Singularity image file via gdown
Press Enter after finishing.
```
sudo pip3 install gdown
gdown https://drive.google.com/uc?id=14v_xUmET-HvCFO3LqmD4sNJL65jBcd0L&export=download
```
or via GitHub
```
git clone https://github.com/hipdac-lab/AMRIC-image.git
cat AMRIC-image/img/amric.sif-* > amric.sif
```
### Step 3: Build and run the image file (need root privilege)
```
sudo singularity build --sandbox artiAmr amric.sif
sudo singularity shell --writable artiAmr
```

### Step 4: Set up environmental variables
```
export OMPI_DIR=/opt/ompi 
export OMPI_VERSION=4.1.1
export PATH=$OMPI_DIR/bin:$PATH
export LD_LIBRARY_PATH=$OMPI_DIR/lib:$LD_LIBRARY_PATH
export MANPATH=$OMPI_DIR/share/man:$MANPATH
export C_INCLUDE_PATH=/opt/ompi/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=/opt/ompi/include:$CPLUS_INCLUDE_PATH
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
```

### Step 5: Run WarpX simulation with no compression, AMReX's original compression, and AMRIC
```
cd /home/wpx256/
. bash.sh
```

### Step 6: Run NYX simulation with no compression, AMReX's original compression, and AMRIC
```
cd /home/nyx128/
. bash.sh
```

### Step 7: Evaluate WarpX's data quality and compression ratio for original AMReX compression and our AMRIC
```
cd /home/wpx256/diags/
. decomp.sh > temp.txt
. qualityCR.sh
```

### Step 8: Evaluate NYX's data quality and compression ratio for original AMReX compression and our AMRIC
```
cd /home/nyx128/run/
. decomp.sh > temp.txt
. qualityCR.sh
```

### Step 9: Compare I/O perf for baselines (i.e., no compression and ori AMReX compression) and AMRIC in WarpX

```
cd /home/wpx256/otfile/
. io.sh
```

### Step 10: Compare I/O performance between baselines and AMRIC in NYX

```
cd /home/nyx128/otfile/
. io.sh
```

## Method 2: Build From Source
### Minimum system & software libraries requirements
OS: Linux (Ubuntu is recommended)

Memory: >= 16 GB RAM

Processor: >= 8 cores

gcc/9.4.0 (or 9.3.0)

cmake (>= 3.23)

OpenMPI/4.1.1 (install scripts provided, or spectrum-mpi)

python/3.8

hdf5/1.12.2 (install scripts provided)


### Step 1: Download AMRIC, checkpoint files, and set up environmental variables (2 mins)
```
git clone https://github.com/hipdac-lab/SC23-AMRIC.git
cd SC23-AMRIC
export AMRIC_HOME=$(pwd)
echo "# start of AMRIC env" >> ~/.bashrc
echo export AMRIC_HOME=$(pwd) >> ~/.bashrc
```
### Step 2: Load or install CMake and numpy. For example, in Ubuntu
```
pip3 install numpy
sudo snap install cmake
sudo snap install cmake --classic
```
### Step 3: Load or install OpenMPI. For example, in Ubuntu (7 mins)
```
sudo bash openmpi.sh 
. mpi_env.sh
```
### Step 4: Download and install the HDF5 library (4 mins)
```
. hdf5.sh
```
### Step 5: Install optimized SZ3 (de)compressor and H5Z-SZ3 (de)compression filter (5 mins)
```
. compressor.sh
```
### Step 6: Install AMReX and Nyx with AMRIC (8 mins)
```
. nyx.sh
```
### Step 7: Install WarpX with AMRIC (9 mins)
```
. warpx.sh
```
### Step 8: Download qcat (compression analysis tool, 1 min)
```
. qcat.sh
```
### Step 9: Run WarpX with no compression, AMReX’s original compression, and AMRIC (3 mins)
```
cd $AMRIC_HOME/warpx_directory/WarpX
. bash.sh
```
### Step 10: Run NYX with no compression, AMReX’s original compression, and AMRIC (3 mins).
```
cd $AMRIC_HOME/Nyx/Exec/AMR-density
. bash.sh
```
### Step 11: Evaluate WarpX’s data quality and compression ratio for original AMReX compression and our AMRIC.
```
cd $AMRIC_HOME/warpx_directory/WarpX/diags
cp $AMRIC_HOME/SZ_SLE/build/tools/H5Z-SZ3/test/des-w .
cp $AMRIC_HOME/orisz3/build/tools/H5Z-SZ3/test/ss-w .
cp $AMRIC_HOME/orisz3/build/tools/H5Z-SZ3/test/stack-w .
cp $AMRIC_HOME/qcat/install/bin/compareData .
. decomp.sh > out.txt
. qualityCR.sh
```
### Step 12: Evaluate NYX’s data quality and compression ratio for original AMReX compression and our AMRIC.
```
cd $AMRIC_HOME/Nyx/Exec/AMR-density/run
cp $AMRIC_HOME/qcat/install/bin/compareData .
cp $AMRIC_HOME/SZ_SLE/build/tools/H5Z-SZ3/test/des .
cp $AMRIC_HOME/orisz3/build/tools/H5Z-SZ3/test/ss .
cp $AMRIC_HOME/orisz3/build/tools/H5Z-SZ3/test/stack .
. decomp.sh > out.txt
. qualityCR.sh
```
### Step 13: Compare I/O performance between baselines (i.e., no compression and ori AMReX compression) and AMRIC in WarpX.
```
cd $AMRIC_HOME/warpx_directory/WarpX/otfile
. io.sh
```
### Step 14: Compare I/O performance between baselines (i.e., no compression and ori AMReX compression) and AMRIC in Nyx.
```
cd $AMRIC_HOME/Nyx/Exec/AMR-density/otfile
. io.sh
```
## Method 2: Build from source
### Minimum system & software libraries requirements


