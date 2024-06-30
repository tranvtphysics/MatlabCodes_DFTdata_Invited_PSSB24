#! /bin/bash
QE_path='/usr/QE/qe-6.7/bin' # Path on the Danymede machine

# mpirun -np 10 $QE_path/pw.x < RUN_scf_pw_QE.in > MoS2.scf.out
mpirun -np 10 $QE_path/pw.x < RUN_nscf_pw_QE.in > MoS2.nscf.out
mpirun -np 10 $QE_path/bands.x < MoS2.bands.in > MoS2.bands.out
$QE_path/plotband.x < MoS2.plotband
# mpirun -np 10 $QE_path/dos.x < RUN_dos_QE.in> MoS2.dos.out
