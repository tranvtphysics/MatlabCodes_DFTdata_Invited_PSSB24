#! /bin/bash

QE_path='/usr/QE/qe-6.7/bin' # Path on the Danymede machine

mpirun -np 10 $QE_path/pw.x < RUN_relax_pw_QE.in > Chain2atoms.relax.out
