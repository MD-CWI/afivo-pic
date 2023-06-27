#! /bin/bash

export OMP_NUM_THREADS=16

./apic CO2_3D.cfg | tee CO2_50mbar_E0125_out.txt