#! /bin/bash

QE_path='/usr/QE/qe-6.7/bin' # Path on the Danymede machine

# Run 1: Converge: force: 1E-3, energy: 1E-8
# mpirun -np 10 $QE_path/pw.x < RUN1_vc-relax_pw_QE.in > Silicene.vc-relax.out1

# Run 2: Converge: force: 1E-3, energy: 1E-10. With the input atomic positions from the final atomic positions in Run 1
mpirun -np 20 $QE_path/pw.x < RUN_vc-relax_pw_QE.in > vc-relax.out


# Run 3: Converge: force: 1E-4, energy: 1E-10