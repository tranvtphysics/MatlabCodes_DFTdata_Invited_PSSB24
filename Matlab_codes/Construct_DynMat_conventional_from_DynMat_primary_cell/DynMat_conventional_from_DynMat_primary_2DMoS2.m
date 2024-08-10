% THE CODE ALLOWS US TO CONSTRUCT DYNAMICAL MATRICES OF THE
% CONVENTIONA-CELL LATTICE FROM THE INFORMATION OF THE PRIMARY-CELL LATTICE

clear all

clc

%% =========================================================================
% INPUT: 
au_of_length_to_Angstrom=0.529177249;
%--------------------------------------------------------------------------
% THE INFORMATION OF THE PRIMARY-CELLS LATTICE

% Force constants of atoms in the primary cells
cd('/media/tranvt/Importance/DATA/RESEARCH/My papers/Drafts/2024.Paper2_DFT_phonon_transport/Paper_PSSA2024_Matlab_codes_and_DFT_data/Some_DFT_data/MoS2-monolayer/NEW_with_vc_relax_forces1Em4/phonon_Nq_odd/Result/')
load D_0_beta_ALL

 D_0_neighbors=D_0_beta_ALL;
 
 
 N1_shells=(size(D_0_neighbors,3)-1)/2;
 N2_shells=(size(D_0_neighbors,4)-1)/2;
 N3_shells=(size(D_0_neighbors,5)-1)/2;
 
 n0_1=N1_shells+1;
 n0_2=N2_shells+1;
 n0_3=N3_shells+1;
 

% lattice vectors of the primary-cell lattice
alat_primary=5.9830;% a.u unit
alat_primary=alat_primary*au_of_length_to_Angstrom;
a1=alat_primary*[1.000000   0.000000   0.000000 ];
a2=alat_primary*[-0.504034   0.872904   0.000000];
a3=alat_primary*[0.000000   0.000000   3.565287];

A=[a1' a2' a3'];% matrix 3x3 containing a1, a2, a3

% Coordinates of atoms in the primary-cell lattice: 1: S1, 2: S2, 3: Mo
X0=[+1.583032996   +4.584306830   -1.566417180 1   0 0 0% this atom in cell [m1,m2,m3]=[0,0,0]
    +1.583032996   +4.584306830   +1.566417180 2   0 0 0
    +1.583032996   +2.741893578   +0.000000000 2   0 0 0]; % this atom in cell [m1,m2,m3]=[0,0,0]


amass=[32.065,32.065,95.94];
mass=amass*1.660540199E-27;% convert masses to kg
%--------------------------------------------------------------------------
% THE INFORMATION OF THE CONVENTIONAL-CELLS LATTICE

% Cordinates of the cell {m1'=0,m2'=0,m3'=0}  
X0_c=[
+1.583032996   +4.584306830   -1.566417180 1
+1.583032996   +4.584306830   +1.566417180 2
+1.583032996   +2.741893578   +0.000000000 3
+3.166065991   +1.842413252   -1.566417180 1
+3.166065991   +1.842413252   +1.566417180 2
+4.749098987   +2.741893578   +0.000000000 3
+4.749098987   +4.584306830   -1.566417180 1
+4.749098987   +4.584306830   +1.566417180 2
+3.166065992   +5.483787156   +0.000000000 3
+6.332131983   +5.483787156   +0.000000000 3
+6.332131982   +1.842413252   -1.566417180 1
+6.332131982   +1.842413252   +1.566417180 2]


% lattice vectors of the conventional-cell lattice
a1_c=[6.3321   0.000000   0.000000];
a2_c=[0.000000   5.4838   0.000000];
a3_c=[0.000000   0.000000   20.00000];
%% =========================================================================
% STEP 0: DETERMINE [m1,m2,m3] INDECES of atoms in cell X0_c
Nat_c=size(X0_c,1);

for i=1:Nat_c
    if X0_c(i,4)==1 % atoms 1
        R0=X0_c(i,1:3)-X0(1,1:3);
        delta_m0=(A\eye(3))*R0';
        
    else
        if X0_c(i,4)==2 % atoms 2
            R0=X0_c(i,1:3)-X0(2,1:3);
            delta_m0=(A\eye(3))*R0';
        else % atoms 3
               R0=X0_c(i,1:3)-X0(3,1:3);
            delta_m0=(A\eye(3))*R0';  
        end
    end
    X0_c(i,5:7)=round(delta_m0);
end

% now X0_c contains the information of [m1,m2,m3] for all atoms 



%% Analyzing K between atoms in the primary-cells lattice:
D_0_0=D_0_beta_ALL(:,:,n0_1,n0_2,n0_3);

nat=size(X0,1);

K33_all=zeros(3,3,2*N1_shells+1,2*N2_shells+1,2*N3_shells+1,nat,nat);
for na=1:nat
    
    D33_n1n1=zeros(3,3);
    
    for nb=1:nat
    mass_n1=mass(na);mass_n2=mass(nb);   
        
        for n1=1:2*N1_shells+1
            for n2=1:2*N2_shells+1
                for n3=1:2*N3_shells+1            

                    D_0_beta=D_0_beta_ALL(:,:,n1,n2,n3);


                            D33_n1_n2=D_0_beta(1+3*(na-1):3*na,1+3*(nb-1):3*nb);                 
                      
                            
                            K33_n1_n2=D33_n1_n2*sqrt(mass_n1*mass_n2);

                            K33_all(:,:,n1,n2,n3,na,nb)=real(K33_n1_n2);
                end
            end
        end
    end
end

%% =========================================================================
% PERFORMANCE: STEP 2: CONSTRUCTING K(m1',m2',m3',na',nb') based on K(m1,m2,m3,na,nb)
N1_shells_c=1;
N2_shells_c=1;
N3_shells_c=0;

K33_c_all=zeros(3,3,2*N1_shells_c+1,2*N2_shells_c+1,2*N3_shells_c+1,Nat_c,Nat_c);

D_0_beta_ALL_c=zeros(3*Nat_c,3*Nat_c,2*N1_shells_c+1,2*N2_shells_c+1,2*N3_shells_c+1);

    
for m1_c=-N1_shells_c:N1_shells_c
    n1_c=m1_c+(N1_shells_c+1);
    for m2_c=-N2_shells_c:N2_shells_c
        n2_c=m2_c+(N2_shells_c+1);
        for m3_c=-N3_shells_c:N3_shells_c
            n3_c=m3_c+(N3_shells_c+1);

           %============================================================================================        
           % PERFORMANCE: STEP 1: DETERMINE X'OF THE CELL {m1',m2',m3'} FROM X'0 OF THE CELL {0,0,0}

            %m1_c=1; m2_c=0; m3_c=0;

            R_c=m1_c*a1_c+m2_c*a2_c+m3_c*a3_c;

            xyz_m1m2m3_c=zeros(Nat_c,3);
            order_m1m2m3_c=zeros(Nat_c,1);
            for i=1:Nat_c
            xyz_m1m2m3_c(i,:)=X0_c(i,1:3)+R_c;

            order_m1m2m3_c(i)=X0_c(i,4);
            end





            %
            delta_m=(A\eye(3))*R_c';

            m1m2m3_m1m2m3_c=zeros(Nat_c,3);
            for i=1:Nat_c
            m1m2m3_m1m2m3_c(i,:)=X0_c(i,5:7)+delta_m';
            end

            %m1m2m3_m1m2m3_c

            X_m1m2m3_c=[xyz_m1m2m3_c order_m1m2m3_c m1m2m3_m1m2m3_c];                 
            %============================================================================================                       
                    
            %============================================================================================                      
            % CONSTRUCTING K(m1',m2',m3',na',nb') or K(n1',n2',n3',na',nb')
            % based on K(m1,m2,m3,na,nb)or K(n1,n2,n3,na,nb)
            
            D_0_beta_c=zeros(3*Nat_c,3*Nat_c);

            for na_c=1:Nat_c  % cell {m1_c=0,m2_c=0,m3_c0=0}
                order_na_in_prim=X0_c(na_c,4);
                m1m2m3_na_in_prim=X0_c(na_c,5:7);

                for nb_c=1:Nat_c % cell {m1_c,m2_c,m3_c0}
                    order_nb_in_prim=X_m1m2m3_c(nb_c,4);
                    m1m2m3_nb_in_prim=X_m1m2m3_c(nb_c,5:7);

                    m1m2m3_nb_compared_na_in_prim=m1m2m3_nb_in_prim-m1m2m3_na_in_prim;
                    m1_nb_compared_na_in_prim=m1m2m3_nb_compared_na_in_prim(1);
                    m2_nb_compared_na_in_prim=m1m2m3_nb_compared_na_in_prim(2);
                    m3_nb_compared_na_in_prim=m1m2m3_nb_compared_na_in_prim(3);

                    n1_nb_compared_na_in_prim=round(m1_nb_compared_na_in_prim+(N1_shells+1));
                    n2_nb_compared_na_in_prim=round(m2_nb_compared_na_in_prim+(N2_shells+1));
                    n3_nb_compared_na_in_prim=round(m3_nb_compared_na_in_prim+(N3_shells+1));

                    if  (abs(m1_nb_compared_na_in_prim)> N1_shells) || (abs(m2_nb_compared_na_in_prim)> N2_shells) || (abs(m3_nb_compared_na_in_prim)> N3_shells)
                        K33_c_all(:,:,n1_c,n2_c,n3_c,na_c,nb_c)=zeros(3);
                    else    
                    K33_c_all(:,:,n1_c,n2_c,n3_c,na_c,nb_c)=K33_all(:,:,n1_nb_compared_na_in_prim,n2_nb_compared_na_in_prim,n3_nb_compared_na_in_prim,order_na_in_prim,order_nb_in_prim);
                    end
                    
                    D_0_beta_c(1+3*(na_c-1):3*na_c,1+3*(nb_c-1):3*nb_c)=K33_c_all(:,:,n1_c,n2_c,n3_c,na_c,nb_c)/sqrt(mass(order_na_in_prim)*mass(order_nb_in_prim));


                end
            end

                D_0_beta_ALL_c(:,:,n1_c,n2_c,n3_c)=D_0_beta_c;

        end
    end
end
%% save data D_0_beta_ALL_c
save D_0_beta_ALL_c D_0_beta_ALL_c