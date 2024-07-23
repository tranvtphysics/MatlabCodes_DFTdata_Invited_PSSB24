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
### 2. Code_to_Plot_Phonon_Bands
	Codes to digonalize the dynamical matrices and compute eigen frequencies to plot phonon bands
### 3. Fourier_transform_for_dynamical_matrix
	Codes of Fourier transform to construct dynamical matrices in real space of cells from dynamical matrices in q-space. Codes are specific for individual systems
### 4. Construct_DynMat_conventional_from_DynMat_primary_cell
	Codes to construct dynamical matrices of conventional cells directly from that of the primary cells