#!/bin/sh

gfortran-9  waiting_time_dist.f08 
./a.out
python3 plot_waiting_distribution.py


