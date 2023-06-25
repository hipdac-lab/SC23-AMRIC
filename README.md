# SC23-AMRIC Artifacts

## Method 1: Use Singularity Image (Recommended)
The entire workflow takes approximately 10 minutes to execute, including downloading container image and preparing environment (3 mins), running WarpX simulation (3 mins), running Nyx simulation (3 mins), and evaluating compression performance (1 min).
### Minimum system requirements
OS: Ubuntu (20.04 is recommended)

Memory: >= 16 GB RAM

Processor: >= 8 cores

### Step 1: Install Singularity
Install [Singularity](https://singularity-tutorial.github.io/01-installation/)

### Step 2: Download the pre-built Singularity image file
Press Enter after finishing.
```
sudo pip3 install gdown
gdown https://drive.google.com/uc?id=14v_xUmET-HvCFO3LqmD4sNJL65jBcd0L&export=download
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

## Method 2: Build from source
### Minimum system & software libraries requirements
OS: Linux

Memory: >= 16 GB RAM

Processor: >= 8 cores

gcc/9.3.0

cmake (>= 3.23)

spectrum-mpi/10.4.0.3

python/3.8

hdf5/1.12.2 (install scripts provided)


