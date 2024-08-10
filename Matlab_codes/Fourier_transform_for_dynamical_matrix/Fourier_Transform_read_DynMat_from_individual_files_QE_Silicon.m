%% The copyright by Van-Truong TRAN
% 1st version: 11/2021

%% THIS CODE TRANSFORS DYNAMYCAL MATRIX FROM RECIPROCAL SPACE (q) to REAL SPACE (r):

% We know the total dynamical matrix D in reciprocal space (D as a function of q vector)
% We need to construct D00, D01, D0-1, D02, D0-2, ... in real space for the interaction
% between cell 0 and cells 1, -1, 2, -2, ...:

% METHOD: FOURIER TRANSFORM
% D(q)= sum_over_R(D_0_beta*exp(i*q*R_0_beta)) with D_0_beta: is the
% dynammycal matrix interactions between cell 0 and cell beta-th.
% R_0_beta=R_beta-R0: vector from cell 0 to cell beta-th.

% D_0_beta = 1/Nq*sum_over_q[D(q)*exp(-i*q*R_0_beta)], Here beta=0,1,-1. R_0_beta=0 if beta=0,

% REQUIRE: WITH nat atoms per unit cell, D should be constructed in the size
% 3*natx3*nat
%% Data input: given the set of data {D(q)} and {q} corresponding to Nq q-points


clear all

% Add the path of the tool:
addpath(genpath('/media/tranvt/Importance/DATA/RESEARCH/My papers/Drafts/2024.Paper2_DFT_phonon_transport/Paper_PSSA2024_Matlab_codes_and_DFT_data/Matlab_codes/Read_DynMat_QE'));
addpath(genpath('/media/tranvt/Importance/DATA/RESEARCH/My papers/Drafts/2024.Paper2_DFT_phonon_transport/Paper_PSSA2024_Matlab_codes_and_DFT_data/Matlab_codes/Fourier_transform_for_dynamical_matrix'));

%%

cd('/media/tranvt/Importance/DATA/RESEARCH/My papers/Drafts/2024.Paper2_DFT_phonon_transport/Paper_PSSA2024_Matlab_codes_and_DFT_data/Some_DFT_data/3D_Silicon/phonon_Nq_odd/Results/') 


nat=2;% number of atoms per cell
amass=[28.0855,28.0855];%[207.2,127.60]; % amass(Pb)=207.2;amss(Te)=127.60; % Mass for MoS2: [95.94,32.065,32.065 ]

number_of_individual_files=10; % number of file .dyn1, .dyn2, ...

filename_dynmat_with_index=sprintf('Si.dyn');

[qpoints,DynaMat_full]=Read_DynMat_QE_in_individual_files_from_ph_calculation(number_of_individual_files,filename_dynmat_with_index,nat,amass);

% qpoints is in unit of 2pi/alat
% DynaMat_full is in SI unit: N/(kg.m)=1/s^2


Nq=size(qpoints,1);
D_total_q=DynaMat_full;

%% Data input: given the set of {R_0_beta} for considered neighbor cells around the cell 0
% if N1_shells, N2_shells, N3_shells are the cut-off numbers of neighbor cells in interaction along a1, a2, a3 
%(with a=([a1;a2;a3]) is the lattice vector

% R_0_beta=i1*a1 + i2*a2 + i3*a3
       

% beta corresponding to the cell index {n1,n2,n3}

au_of_length_to_Angstrom=0.529177249;

alat_au=  10.2623; % a.u unit of length
alat=alat_au*au_of_length_to_Angstrom
Lattice_vector=alat*[-0.500000   0.000000   0.500000
                     0.000000   0.500000   0.500000
                     -0.500000   0.500000   0.000000];
                 
                 
a1=Lattice_vector(1,:);
a2=Lattice_vector(2,:);
a3=Lattice_vector(3,:);                 
                 
q_points=qpoints*2*pi/alat;

figure(10)
XX=0;
for n=1:Nq
XX=XX+1/Nq*exp(1i*dot(q_points(n,:),Lattice_vector(1,:)));
hold on
plot(q_points(n,1),q_points(n,2),'or')
end
XX % this number must go to zeros
xlabel('q_x')
ylabel('q_y')
%% Check if q along a1, a2,a3 satisfy: -pi<=q1*a1<=pi. Or q1=n1*b1, so -1/2<= n1 <=1/2
% Lattice vectors in reciprical space
[b1,b2,b3] = Reciprocal_Lattice_Vectors(a1,a2,a3);

q1_set=zeros(Nq,3);
n1_set=zeros(Nq,1);

q2_set=zeros(Nq,3);
n2_set=zeros(Nq,1);

q3_set=zeros(Nq,3);
n3_set=zeros(Nq,1);
for n=1:Nq

q1_n=q_points(n,1);
q2_n=q_points(n,2);    
q3_n=q_points(n,3);

cos_q_b1=dot(q_points(n,:),b1)/norm(q_points(n,:))/norm(b1);
q1_set(n,:)=norm(q_points(n,:))*cos_q_b1*b1/norm(b1);   
n1_set(n)=dot(q1_set(n,:),b1)/dot(b1,b1);

cos_q_b2=dot(q_points(n,:),b2)/norm(q_points(n,:))/norm(b2);
q2_set(n,:)=norm(q_points(n,:))*cos_q_b2*b2/norm(b2);   
n2_set(n)=dot(q2_set(n,:),b2)/dot(b2,b2);   

cos_q_b3=dot(q_points(n,:),b3)/norm(q_points(n,:))/norm(b3);
q3_set(n,:)=norm(q_points(n,:))*cos_q_b3*b3/norm(b3);   
n3_set(n)=dot(q3_set(n,:),b3)/dot(b3,b3);   
end
[min(n1_set), max(n1_set)]  % this number must go to zeros

[min(n2_set), max(n2_set)] 

[min(n3_set), max(n3_set)] 
%% COMPUTING D_0_beta  BY FOURIER TRANSFORM:

% Check cut-off neighboring shells (N1, N2, N3)

% [N1_shells,N2_shells,N3_shells] = Check_CutOff_neighbors_in_Fourier_transform(Lattice_vector,q_points,D_total_q);


N1_shells=2; % equal Nq1/2 with Nq1 in ph.x calculation
N2_shells=2; % equal Nq2/2 with Nq2 in ph.x calculation
N3_shells=2; % equal Nq3/2 with Nq3 in ph.x calculation
% Compute the full set of D_0_beta
    n0_1=N1_shells+1;
    n0_2=N2_shells+1;
    n0_3=N3_shells+1;

D_0_beta_ALL=zeros(3*nat,3*nat,2*N1_shells+1,2*N2_shells+1,2*N3_shells+1);

for i1=1:2*N1_shells+1
    for i2=1:2*N2_shells+1
        for i3=1:2*N3_shells+1
            
            R_0_beta=(i1-n0_1)*a1 + (i2-n0_2)*a2 + (i3-n0_3)*a3;
            
            D_0_beta=zeros(3*nat,3*nat);
            for n=1:Nq
                q=q_points(n,:);
                D_q=D_total_q{n};
                D_0_beta=D_0_beta+1/Nq*D_q*exp(-1i*dot(q,R_0_beta)); % ALSO RE-COMPUATE D_0_0
            end
            D_0_beta_ALL(:,:,i1,i2,i3)=D_0_beta;
        end
    end
end    

%% Applying Acoustic Sum rule:
  [D_0_beta_ALL]=Acoustic_Sum_check_diaginal_blocks_of_DynMat_real_space(D_0_beta_ALL,amass,N1_shells,N2_shells,N3_shells);

save D_0_beta_ALL D_0_beta_ALL
%% TEST
% D_0_0
% 
% D_0_beta_ALL(:,:,n0_1,n0_2,n0_3)
% 
% D_0_100=D_0_beta_ALL(:,:,n0_1+1,n0_2,n0_3)
% 
% D_0_200=D_0_beta_ALL(:,:,n0_1+2,n0_2,n0_3)

%% Plot phonon bands from D_0_beta in real space obtained from Furier transform

% q_points to plot bands: Plot only along x direction

% Nq_1=100;
% 
% Nq=Nq_1;
% q_points=zeros(Nq,3);
% n=0;
% for n1=1:Nq_1
%     qx=(n1-1)*(pi/a1(1))/(Nq_1-1);
%     n=n+1;
%     q_points(n,:)=[qx,0,0];
% 
% end

%  cd ..
%  cd ..
prefix='Si.dynmat';%'PbTe.dynmat';
% read dynamic matrix file: prefix.dynmat
[q_points_path_0,DynaMat_full_0]=Read_DynMat_QE(prefix,nat,amass);

q_points_path=q_points_path_0*2*pi/alat;
Nq_path=size(q_points_path,1);

Omega=zeros(Nq_path,3*nat);
for n=1:Nq_path
    q=q_points_path(n,:);
    
    
    D_q=zeros(3*nat,3*nat);
    
    for i1=1:2*N1_shells+1
        for i2=1:2*N2_shells+1
            for i3=1:2*N3_shells+1

                R_0_beta=(i1-n0_1)*a1 +(i2-n0_2)*a2 + (i3-n0_3)*a3;

                D_0_beta=D_0_beta_ALL(:,:,i1,i2,i3);
                
                D_q=D_q+D_0_beta*exp(1i*dot(q,R_0_beta)) ;
                
                %D_q;

            end
        end
    end
    
    
    Omega2=sort(real(eig(D_q)));
    
    for i=1:3*nat
        if (real(Omega2(i))>=0)
        Omega(n,i)=sqrt(real(Omega2(i))); % unit 1/s=Hz
        else
        Omega(n,i)=-sqrt(-real(Omega2(i))); % unit 1/s=Hz 
        end
    end    
    
    
end

%% PLOT RESULTS:
cmm1_to_Hz=2.99793*1e10; % equivalent cmm1_to_secondm1

%Freq=Omega/(2*pi)*1e-12;% convert to Thz

Freq=Omega/(2*pi)/cmm1_to_Hz;% convert to cm-1

qx_points=q_points_path(:,1);

hold on
figure(1)
set(gca,'Fontsize',24)
plot(([1:Nq_path]-1)/(Nq_path-1),Freq,'k:','linewidth',2);
legend(sprintf('Using Dyn(r) obtained from FourierTrans, N_1xN_2xN_3=%dx%dx%d',N1_shells,N2_shells,N3_shells))
title('Si: with correcting diagonal blocks (or Acoustic Sum rule)')
box on
xlim([0,Nq_path-1]/(Nq_path-1))
xlabel('q-path')
%ylabel('frequency (Thz)')
ylabel(' Frequency (cm^{-1})','Fontsize',24);

%% save figures
string1=sprintf('Phonon_bands_%s_from_FourTrans.png',filename_dynmat_with_index);
saveas(gcf,string1);% save figure
string2=sprintf('Phonon_bands_%s_from_FourTrans.fig',filename_dynmat_with_index);
saveas(gcf,string2);% save figure
