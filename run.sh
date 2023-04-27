#!/bin/bash

make
cd $PWD
./LBM > result.txt
make clean