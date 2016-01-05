#!/bin/bash

ffile=$1
export ffile
echo $ffile

EXECUTABLE=$2
export EXECUTABLE


gfortran -fopenmp -g -O2 -save-temps -I/usr/local/netcdf-4.1.2/f90 -I/usr/local/hdf-4.2.5/include -o $EXECUTABLE -lnetcdff -lnetcdf -L/usr/local/hdf5-1.8.5/lib -lhdf5 -lhdf5_hl -L/usr/local/lib -L/usr/lib -lcurl -lz $ffile
# -L/usr/local/netcdf-4.1.2/lib -L/usr/local/hdf-4.2.5/lib -lmfhdf -ldf -ljpeg