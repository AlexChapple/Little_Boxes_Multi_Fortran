#!/bin/sh

gfortran-9 testing.f08 
./a.out
python3 plot_results.py


