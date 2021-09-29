#!/bin/sh
# Change this to point to the right MPI C,C++ compiler
export CC=gcc
export CXX=g++
export FC=gfortran
NUM_PROCS=2

##############################################
############  DO NOT MODIFY #################
##############################################
INSTALL_PATH=$PWD/DIST-ompss-2
ROOT_DIR=$PWD
EXTRA=$PWD/extra

# Install metis
cd $ROOT_DIR
tar -zxf metis-5.1.0.tar.gz
cd metis-5.1.0
make config cc=$CC prefix=$INSTALL_PATH
make -j$NUM_PROCS install && cd $ROOT_DIR && rm -rf metis-5.1.0

#HDF5
cd $ROOT_DIR
tar -jxvf hdf5-1.8.21.tar.bz2
cd hdf5-1.8.21
./configure --prefix=$INSTALL_PATH --enable-fortran --disable-shared
make -j$NUM_PROCS install && cd $ROOT_DIR && rm -rf hdf5-1.8.21

