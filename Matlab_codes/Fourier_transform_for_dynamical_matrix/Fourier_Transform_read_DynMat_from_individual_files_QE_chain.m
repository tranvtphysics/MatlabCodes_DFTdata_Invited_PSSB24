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

cd('/media/tranvt/Importance/DATA/RESEARCH/My papers/Drafts/2024.Paper2_DFT_phonon_transport/Paper_PSSA2024_Matlab_codes_and_DFT_data/Some_DFT_data/Linear_chain_2atoms_B_N/phonon/Result/')

nat=2;% number of atoms per cell
amass=[10.811,14.0067];%[207.2,127.60]; % amass(Pb)=207.2;amss(Te)=127.60; % Mass for MoS2: [95.94,32.065,32.065 ]


number_of_individual_files=5; % number of file .dyn1, .dyn2, ...

filename_dynmat_with_index=sprintf('Chain2atoms.dyn');

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

alat=  4.9525;% for chain C12-C48: 6.61404144; % a.u unit of length
Lattice_vector=alat*[1.000000   0.000000   0.000000
                     0.000000   7.631422   0.000000 % for chain C12-C48: 0.000000   5.714286   0.000000
                     0.000000   0.000000   7.631422]; % for chain C12-C48: 0.000000   0.000000   5.714286];
                 
Lattice_vector*au_of_length_to_Angstrom
                 
a1=Lattice_vector(1,:);
a2=Lattice_vector(2,:);
a3=Lattice_vector(3,:);                 
                 
q_points=qpoints*2*pi/alat;

XX=0;
for n=1:Nq
XX=XX+1/Nq*exp(1i*dot(q_points(n,:),Lattice_vector(1,:)));
end
XX % this number must go to zeros
%% COMPUTING D_0_beta  BY FOURIER TRANSFORM:

% Check cut-off neighboring shells (N1, N2, N3)

% [N1_shells,N2_shells,N3_shells] = Check_CutOff_neighbors_in_Fourier_transform(Lattice_vector,q_points,D_total_q);


N1_shells=4;
N2_shells=0;
N3_shells=0;

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

    
%save D_0_beta_ALL            
%% TEST
% D_0_0
% 
% D_0_beta_ALL(:,:,n0_1,n0_2,n0_3)
% 
% D_0_100=D_0_beta_ALL(:,:,n0_1+1,n0_2,n0_3)
% 
% D_0_200=D_0_beta_ALL(:,:,n0_1+2,n0_2,n0_3)


%% Applying Acoustic Sum rule:
[D_0_beta_ALL]=Acoustic_Sum_check_diaginal_blocks_of_DynMat_real_space(D_0_beta_ALL,amass,N1_shells,N2_shells,N3_shells);

%% Plot phonon bands from D_0_beta in real space obtained from Furier transform

% q_points to plot bands: Plot only along x direction

Nq_1=100;

Nq=Nq_1;
q_points=zeros(Nq,3);
n=0;
for n1=1:Nq_1
    qx=(n1-1)*(pi/a1(1))/(Nq_1-1);
    n=n+1;
    q_points(n,:)=[qx,0,0];

end

%===================================================================================
%  %Test using shorter N_shell cut-off
% 
% 
% N1_shells_cut=N1_shells-2; % Cut-off considered range N1_shells_cut<=N1_shells
% N2_shells_cut=N2_shells; % Cut-off considered range N1_shells_cut<=N1_shells
% N3_shells_cut=N3_shells; % Cut-off considered range N1_shells_cut<=N1_shells
% 
% 
% 
% D_0_beta_ALL_cut=zeros(3*nat,3*nat,2*N1_shells_cut+1,2*N2_shells_cut+1,2*N3_shells_cut+1);
%     for i1=1:2*N1_shells_cut+1
%         i1_old=(n0_1-N1_shells_cut)-1+i1;
%         for i2=1:2*N2_shells_cut+1
%             i2_old=(n0_2-N2_shells_cut)-1+i2;
%             for i3=1:2*N3_shells_cut+1
%                 i3_old=(n0_3-N3_shells_cut)-1+i3;
%                     D_0_beta_ALL_cut(:,:,i1,i2,i3)=D_0_beta_ALL(:,:,i1_old,i2_old,i3_old);
% 
%             end
%         end
%     end
% 
% [D_0_beta_ALL_cut]=Acoustic_Sum_check_diaginal_blocks_of_DynMat_real_space(D_0_beta_ALL_cut,amass,N1_shells_cut,N2_shells_cut,N3_shells_cut);
%     
%     
% 
%     n0_1_new=N1_shells_cut+1;
%     n0_2_new=N2_shells_cut+1;
%     n0_3_new=N3_shells_cut+1;
% 
% 
% Omega=zeros(Nq,3*nat);
% for n=1:Nq
%     q=q_points(n,:);
%     
%     
%     D_q=zeros(3*nat,3*nat);
%     
%     for i1=1:2*N1_shells_cut+1
%         for i2=1:2*N2_shells_cut+1
%             for i3=1:2*N3_shells_cut+1
% 
%                 R_0_beta=(i1-n0_1_new)*a1 + (i2-n0_2_new)*a2 + (i3-n0_3_new)*a3;
% 
%                 D_0_beta=D_0_beta_ALL_cut(:,:,i1,i2,i3);
%                 
%                 D_q=D_q+D_0_beta*exp(1i*dot(q,R_0_beta)); 
% 
%             end
%         end
%     end
%     
%     
%     Omega2=sort(real(eig(D_q)));
%     
%     for i=1:3*nat
%         if (real(Omega2(i))>=0)
%         Omega(n,i)=sqrt(real(Omega2(i))); % unit 1/s=Hz
%         else
%         Omega(n,i)=-sqrt(-real(Omega2(i))); % unit 1/s=Hz 
%         end
%     end    
%     
%     
% end
%  % End Test shorter N_shell cut-off
%===================================================================================


Omega=zeros(Nq,3*nat);
for n=1:Nq
    q=q_points(n,:);
    
    
    D_q=zeros(3*nat,3*nat);
    
    for i1=1:2*N1_shells+1
        for i2=1:2*N2_shells+1
            for i3=1:2*N3_shells+1

                R_0_beta=(i1-n0_1)*a1 + (i2-n0_2)*a2 + (i3-n0_3)*a3;

                D_0_beta=D_0_beta_ALL(:,:,i1,i2,i3);
                
                D_q=D_q+D_0_beta*exp(1i*dot(q,R_0_beta)); 

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
%D_q
 
 
%% PLOT RESULTS:
cmm1_to_Hz=2.99793*1e10; % equivalent cmm1_to_secondm1

% Freq=Omega/(2*pi)*1e-12;% convert to Thz

Freq=Omega/(2*pi)/cmm1_to_Hz;% convert to cm-1

qx_points=q_points(:,1);

%=============================================================================
npath=1;
segment_qpoints=[50];%Graphene: [40,20,40,20], {'G','K', 'M', 'G'};
Labels={'G','X'}; %{'G','X'};%{'L','G', 'X', 'W', 'G'};%{'G','X'};%
%Labels={'', '', '', '', ''};

Labels_position=zeros(1,npath+1);
Labels_position(1)=1;
for i=1:npath
Labels_position(i+1)=sum(segment_qpoints(1:i))+1;
end


Nq_path=Nq; % all q-point in the q-paths considered
q_distance=1:Nq_path;%
%=============================================================================
figure(1)
h1=plot(([1:Nq_path]-1)/(Nq_path-1),Freq,':k','linewidth',3); % with the band-crossing treatment

for i=2:npath
ln(i-1)=line([q_distance(Labels_position(i))-1 q_distance(Labels_position(i))-1]/(Nq_path-1), [-10 1.1*max(max(Freq))],'LineStyle','--');
end
ylabel(' Frequency (cm^{-1})','Fontsize',24);
xlabel('q-paths ','Fontsize',24);
set(gca,'Xtick',(q_distance(Labels_position)-1)/(Nq_path-1),'XTickLabel',Labels,'Fontsize',24,'Color','none');

set(gca,'Fontsize',24)
box on
%xlim([0,max(q_distance)])
ylim([-50,2200])

legend(sprintf('Using Dyn(r) obtained from FourierTrans, N_1xN_2xN_3=%dx%dx%d',N1_shells,N2_shells,N3_shells))
title('with correcting diagonal blocks (or Acoustic Sum rule)')
%% save figures
% string1=sprintf('Phonon_bands_with_D_realspace_from_FourierTransf_%dx%dx%d.png',N1_shells,N2_shells,N3_shells);
% saveas(gcf,string1);% save figure
% string2=sprintf('Phonon_bands_with_D_realspace_from_FourierTransf_%dx%dx%d.png.fig',N1_shells,N2_shells,N3_shells);
% saveas(gcf,string2);% save figure



