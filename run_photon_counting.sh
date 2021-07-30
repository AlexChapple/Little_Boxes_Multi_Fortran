#!/bin/sh

gfortran-9  photon_counting.f08 -o photon_counting.out
./photon_counting.out
python3 photon_counting.py 


