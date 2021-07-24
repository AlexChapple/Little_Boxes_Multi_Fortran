#!/bin/sh

gfortran-9 Little_Boxes_Multi_RK.f08 
./a.out
python3 plot_results.py


