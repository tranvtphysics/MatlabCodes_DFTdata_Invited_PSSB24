#! /bin/bash
QE_path='/usr/QE/qe-6.7/bin' # Path on the Danymede machine

mpirun -np 30 $QE_path/pw.x < RUN_scf_pw_QE.in > 5AGNR.scf.out
mpirun -np 30 $QE_path/pw.x < RUN_nscf_pw_QE.in > 5AGNR.nscf.out
mpirun -np 20 $QE_path/bands.x < 5AGNR.bands.in > 5AGNR.bands.out
$QE_path/plotband.x < 5AGNR.plotband
mpirun -np 1 $QE_path/dos.x < RUN_dos_QE.in> 5AGNR.dos.out
