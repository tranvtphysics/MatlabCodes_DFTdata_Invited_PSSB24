
function [D_0_beta_ALL]=Acoustic_Sum_check_diaginal_blocks_of_DynMat_real_space(D_0_beta_ALL,amass,N1_shells,N2_shells,N3_shells)


    n0_1=N1_shells+1;
    n0_2=N2_shells+1;
    n0_3=N3_shells+1;
    
    
    nat=length(amass); % number of atoms per cell
%% Applied Acoustic sum:

mass=amass*1.660540199E-27;% convert masses to kg

D_0_0=D_0_beta_ALL(:,:,n0_1,n0_2,n0_3);

D33_nn=zeros(3,3,nat);
for n1=1:nat
    
    D33_n1n1=zeros(3,3);
    
    for n2=1:nat
    mass_n1=mass(n1);mass_n2=mass(n2);   
        
        for i1=1:2*N1_shells+1
            for i2=1:2*N2_shells+1
                for i3=1:2*N3_shells+1            

                    D_0_beta=D_0_beta_ALL(:,:,i1,i2,i3);



                            D33_n1_n2=D_0_beta(1+3*(n1-1):3*n1,1+3*(n2-1):3*n2);
                            
                            
                            if (n1==n2) && (i1==n0_1) && (i2==n0_2) && (i3==n0_3) 
                                K33_n1_n2=zeros(3,3);
                            else
                            

                            
                            K33_n1_n2=D33_n1_n2*sqrt(mass_n1*mass_n2);
                            end
                            
                            D33_n1n1=D33_n1n1-K33_n1_n2/mass_n1; % acoustic sum

                end
            end
        end
    end
    
D33_nn(:,:,n1)=D33_n1n1;
 end
    
    
    
for n=1:nat
    n1=n;n2=n;
    D_0_0(1+3*(n1-1):3*n1,1+3*(n2-1):3*n2)=D33_nn(:,:,n);
end
 
D_0_beta_ALL(:,:,n0_1,n0_2,n0_3)=D_0_0;