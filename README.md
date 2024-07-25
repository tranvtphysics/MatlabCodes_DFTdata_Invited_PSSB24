# DFT data and Matlab codes for the invited paper pssb 2024

## DFT data
### 1. 1D structures:
	- Linear_chain_2atoms_B_N
	- 5_AGNR_H_passive
	- 6_AGNR_H_passive
### 2. 2D structures:
	- 2D-graphene
	- MoS2-monolayer

### 3. 3D structures:
	- 3D_Silicon

## Matlab codes
### 1. Read_DynMat_QE
	Codes to read data from ph.x calculations and construct dynamical matrices of cells
	- Read_DynMat_QE.m
	Matlab function helps to read data of force constant in a single file *.dyn1 or *.dyn2, ... and construct dynamical matrix at a specific q-point.
	- Read_DynMat_QE_in_individual_files_from_ph_calculation.m
	Matlab function helps to read data from many individual files *.dyn and construct dynamical matrices at different q-points.
	- Construct_DynMat_by_reading_fc_file_from_q2r_QE.m 
	Matlab funcion help to read interatomic force constant from the .fc file generated by q2r.x calculations and constructure dynamical matrices in real space.

### 2. Code_to_Plot_Phonon_Bands
	Codes to digonalize the dynamical matrices and compute eigen frequencies to plot phonon bands
	- Compute_Phonon_bands_QE_from_dynmat_file.m
	This Matlab script helps to diagonalize the Dynamical matrix (q) 3Nax3Na constructed from data of force constants in the file *.dynmat (by using the Matlab function: Read_DynMat_QE.m)
	- Read_freq_file_by_QE_and_plot_Phonon_band.m
	This Matlab script helps to read data of phonon frequency modes in the file *.freq generated by phonon QE calculations and then plot the phonon bands.	
### 3. Fourier_transform_for_dynamical_matrix
	Codes of Fourier transform to construct dynamical matrices in real space of cells from dynamical matrices in q-space. Codes are specific for individual systems
	- Acoustic_Sum_check_diaginal_blocks_of_DynMat_real_space.m
	Matlab function to correct diagonal block lines following the acoustic sum rule
	- Fourier_Transform_read_DynMat_from_individual_files_QE_chain.m
	- Fourier_Transform_read_DynMat_from_individual_files_QE_Graphen.m
	- Fourier_Transform_read_DynMat_from_individual_files_QE_MoS2.m
	- Fourier_Transform_read_DynMat_from_individual_files_QE_Silicon.m
### 4. Construct_DynMat_conventional_from_DynMat_primary_cell
	Codes to construct dynamical matrices of conventional cells directly from that of the primary cells
	- DynMat_conventional_from_DynMat_primary_2DGraphene.m
	- DynMat_conventional_from_DynMat_primary_2DMoS2.m