git clone https://github.com/szcompressor/qcat.git
cd qcat
mkdir install 
./configure --prefix=$AMRIC_HOME/qcat/install
make -j 8
make install
cd $AMRIC_HOME
