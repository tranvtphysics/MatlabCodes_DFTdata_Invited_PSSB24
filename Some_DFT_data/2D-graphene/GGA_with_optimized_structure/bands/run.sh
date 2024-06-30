#! /bin/bash
QE_path='/usr/QE/qe-6.7/bin' # Path on the Danymede machine

mpirun -np 20 $QE_path/pw.x < RUN_scf_pw_QE.in > 2Dgraphene.scf.out
mpirun -np 20 $QE_path/pw.x < RUN_nscf_pw_QE.in > 2Dgraphene.nscf.out
mpirun -np 20 $QE_path/bands.x < 2Dgraphene.bands.in > 2Dgraphene.bands.out
$QE_path/plotband.x < 2Dgraphene.plotband
mpirun -np 20 $QE_path/dos.x < RUN_dos_QE.in> 2Dgraphene.dos.out
