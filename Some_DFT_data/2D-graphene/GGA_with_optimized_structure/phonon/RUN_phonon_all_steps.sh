#! /bin/bash
QE_path='/usr/QE/qe-6.7/bin' # Path on the Danymede machine

# step 0: run scf to have the converged charge density
mpirun -np 20 $QE_path/pw.x < RUN_scf_pw_QE.in > 2Dgraphene.scf.out

# step1_Phonon_QE to have .dyn matrices (the most expensive step: should use many cores with this step)
mpirun -np 30 $QE_path/ph.x < RUN_step1_Dyn_on_q_ph_QE.in > print_step1_phonon.out

# step2_Dyn_realspace_Phonon_QE
$QE_path/q2r.x < RUN_step2_Dyn_on_real_space_q2r_QE.in > print_step2_q2r_phonon.out

# step3_PhoBands_QE.sh
$QE_path/matdyn.x < RUN_step3_Phobands_matdyn_QE.in > print_step3_Phonon_bands.out
$QE_path/plotband.x < RUN_plotband_Phonon.in

# step4_PhoDOS_QE
$QE_path/matdyn.x < RUN_step4_PhoDOS_matdyn_QE.in > print_step4_phonon_DOS.out
