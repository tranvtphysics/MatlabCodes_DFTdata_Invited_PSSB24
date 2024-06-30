#! /bin/bash

QE_path='/usr/QE/qe-6.7/bin' # Path on the Danymede machine

mpirun -np 30 $QE_path/pw.x < RUN_vc-relax_pw_QE.in > 5AGNR.vc-relax.out
