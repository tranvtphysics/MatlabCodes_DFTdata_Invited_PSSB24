function [D_0_beta_ALL,R_0_beta_ALL,Position_cells]=Construct_DynMat_by_reading_fc_file_from_q2r_QE(filename_fc,amass,Lattice_vector)
% The copyright by Van-Truong TRAN.
% Code version: V1.0 (02/2022)

% READING interatomic force constant from the .fc file generated by q2r.x (QE) calculculation
% filename_fc: file data is "prefix.fc". Format: string
% nat: number of atoms in the unit cell in QE calculation. Format: integer
% amass: atomic mass of each kind of atoms. Format: real

%%
% Dir='/media/tranvt/Importance/DATA/BACKUP/Onedrive/RESEARCH/C2N_2021/SIMULATIONS/DFT_QE/2D-graphene/GGA_new/RESULTS/phonon';
% cd(Dir);
% 
% filename_fc='gr.fc';
% nat=2;
% amass=[12,12];%[10.811,14.0067];

mass=amass*1.660540199E-27;% convert masses to kg


Ryd_to_Nm=2.1798741E-18;% 1 Ryd=2.1798741E-18 N*m;

au_of_length_to_meter=5.29177249E-11;
%% Reading the file "prefix.fc"
% the content of the file looks like this:
% 2    2  8  4.9524876  7.6314219  7.6314219  0.0000000  0.0000000  0.0000000
%            1  'B  '    9853.62371224769     
%            2  'N  '    12766.3260799500     
%     1    1     -0.1411180033      0.0000000000      0.0000000000
%     2    2      0.3587817939      0.0000000000      0.0000000000
%  T
%           1.187375343203          0.000000000000          0.000000000000
%           0.000000000000          1.020407615077          0.000000000000
%           0.000000000000          0.000000000000          1.020407615077
%     1
%       7.3740412      0.0000000      0.0000000
%       0.0000000      0.4561708      0.0000000
%       0.0000000      0.0000000      0.4561708
%     2
%      -7.3740412      0.0000000      0.0000000
%       0.0000000     -0.4561708      0.0000000
%       0.0000000      0.0000000     -0.4561708
%    8   1   1
%    1   1   1   1
%    1   1   1   2.28533168892E+00
%    2   1   1  -1.63435062033E-01
%    3   1   1  -7.04399109463E-02
%    4   1   1  -4.45511288173E-02
%    5   1   1  -3.72102409688E-02
%    6   1   1  -4.45511288173E-02
%    7   1   1  -7.04399109463E-02
%    8   1   1  -1.63435062033E-01
%    1   1   1   2
%    1   1   1  -1.04115111417E+00
%    2   1   1  -1.04551837301E+00
%    3   1   1   1.06689616276E-01
%    4   1   1   5.62198889597E-02
%    5   1   1   4.05771663398E-02
%    6   1   1   4.05705591369E-02
%    7   1   1   5.62674923719E-02
%    8   1   1   1.06758539743E-01


file_data=sprintf('%s',filename_fc);

fileID = fopen(file_data,'r');

[status1,cmdout1]=system(sprintf('wc -l < %s',file_data)); % Counting number of data lines in the data file
Nlines=str2num(cmdout1);

number_lines_scaned=Nlines; % Look at txt file to see how many lines need to be scaned
Intro = textscan(fileID,'%s',number_lines_scaned,'Delimiter','\n');
IF=Intro{1};
fclose(fileID);


%% read information of the first line:
Line_1=strtrim(IF{1});
Line_1_values=str2num(Line_1);

Nb_species=Line_1_values(1);
Nb_atoms=Line_1_values(2);

Type_of_lattice=Line_1_values(3);

Lattice_parameters=Line_1_values(4:end);

%% From line 2nd to 1 + Nb_species: the information of species and their atomistic mass in a.u



%% From line 1 + Nb_species + 1 to 1 + Nb_species + Nb_atoms: the information of atoms: order, species-index, coordinates in a0 units
n=0;
species_xyz_atoms=zeros(Nb_atoms,4);
for i=1 + Nb_species + 1:1 + Nb_species + Nb_atoms
Line_i=strtrim(IF{i});
Line_i_values=str2num(Line_i);
n=n+1;
species_xyz_atoms(n,:)=Line_i_values(2:end);

end
%species_xyz_atoms ;

%% Line 1 + Nb_species + Nb_atoms + 1
% T if the file contains epsilon and Z*, F otherwise

% If the letter is "T":
% from line (1 + Nb_species + Nb_atoms + 1) + 1 to (1 + Nb_species +
% Nb_atoms + 1) + 3 : The information of epsilon
% From line (1 + Nb_species + Nb_atoms + 1 + 3) + 1 : (1 + Nb_species + Nb_atoms + 1 + 3) + 4*Nb_atoms : the information of Z*

%% Information of grid q: matrix 3x1 : nq1xnq2xnq3
% If line Line 1 + Nb_species + Nb_atoms + 1 contains the letter "T":
% line (1 + Nb_species + Nb_atoms + 1 + 3 + 4*Nb_atoms) +1 shows the gird-q

% If If line Line 1 + Nb_species + Nb_atoms + 1 contains the letter "F":
% line (1 + Nb_species + Nb_atoms + 1) +1 shows the gird-q

Line_letter_F_or_T=strtrim(IF{1 + Nb_species + Nb_atoms + 1});
if strcmp(Line_letter_F_or_T,'T')==1
    index_line_grid_q=(1 + Nb_species + Nb_atoms + 1 + 3 + 4*Nb_atoms) +1;
    grid_q=str2num(IF{index_line_grid_q});
end

if strcmp(Line_letter_F_or_T,'F')==1
    index_line_grid_q=(1 + Nb_species + Nb_atoms + 1) +1;
    grid_q=str2num(IF{index_line_grid_q});
end

%grid_q;

Nq1=grid_q(1);
Nq2=grid_q(2);
Nq3=grid_q(3);
%% After line of the grid q: force constant between atoms in different directions
%     1   1   1   1  --> ??
% indices: i,j,na,nb (polarization1 , polarization 2, atom 1 , atom 2)

% m1 m2 m3  force_constant(m1,m2,m3,i,j,na,nb) 
% where m1=1,...,n1, m2=1,...,n2, m3=1,...,n3, define a lattice vector R:
%    R = (m1-1)*tau1 + (m2-1)*tau2 + (m3-1)*tau3
% (if m1-1 > n1/2 and so on, refold m1-1 => m1-1-n1 and so on)
%---------------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------------
% My format:
IFC_all=zeros(Nq1,Nq2,Nq3,Nb_atoms,Nb_atoms,3,3);
% note: Nq1,Nq2,Nq3: indicate the cell neighbors: m1=1:Nq1,m2=1:Nq2,m3=1:Nq3
% Nb_atoms,Nb_atoms: indicate atom index: na,nb: atom na in cell 0, nb in
% cell (m1,m2,m3)
% 3,3: indicate direction x, y, z (1:x, 2:y, 3:z)

% IFC_all(m1,m2,m3,na,nb,i,j)

% Now read information of force constant between atoms in different
% directions xx,xy,xz,... and put into this format matrix

% Line index_line_grid_q + 1: i,j,na,nb (polarization1 , polarization 2, atom 1 , atom 1)
% From Line (index_line_grid_q + 1)+1 to (index_line_grid_q + 1)+Nq1xNq2xNq3: m1 m2 m3  force_constant(m1,m2,m3,i,j,na,nb) 

% Thus: there are 3x3xNb_atomsxNb_atoms blocks, each block has (Nq1xNq2xNq3+1) rows

Nb_blocks=3*3*Nb_atoms*Nb_atoms;

Nb_lines_each_block=Nq1*Nq2*Nq3+1;

for i_bl=1:Nb_blocks
    index_line_ij_nanb=index_line_grid_q+1+(i_bl-1)*Nb_lines_each_block;
    
    ij_nanb=str2num(IF{index_line_ij_nanb});
    i=ij_nanb(1);
    j=ij_nanb(2);
    na=ij_nanb(3);
    nb=ij_nanb(4);
    
    for m=1:Nq1*Nq2*Nq3
      index_line_m1_m2_m3_FC=index_line_ij_nanb+m; 
    
      m1_m2_m3_FC=str2num(IF{index_line_m1_m2_m3_FC});
      
      m1=m1_m2_m3_FC(1);
      m2=m1_m2_m3_FC(2);
      m3=m1_m2_m3_FC(3);
      FC=m1_m2_m3_FC(4); % units of Ry/au^2
      
      IFC_all(m1,m2,m3,na,nb,i,j)=FC;% units of Ry/au^2
      
    end
    
end


%% Constructing Matrix force constant between atoms:
IFC_all_3x3=cell(Nq1,Nq2,Nq3,Nb_atoms,Nb_atoms);
for m1=1:Nq1
    for m2=1:Nq2
        for m3=1:Nq3
            for na=1:Nb_atoms
                for nb=1:Nb_atoms
                    
                    IFC_na_nb_3x3=zeros(3);
                    for i=1:3                  
                       
                        for j=1:3                        
                        IFC_na_nb_3x3(i,j)=IFC_all(m1,m2,m3,na,nb,i,j); % Note: IFC_all(m1,m2,m3,na,nb,i,j) in units of Ry/au^2
                        end     
                                       
                    end
                    
                    IFC_all_3x3{m1,m2,m3,na,nb}=IFC_na_nb_3x3*Ryd_to_Nm/au_of_length_to_meter^2;  % Convert to SI unit: N/m;
                end
            end
        end
    end
end
            

%% Constructing the Dynamical matrices: interaction between cell 0 and the neighbor cells

% au_of_length_to_Angstrom=0.529177249;
% 
% alat_au=  4.6531; % a.u unit of length
% alat=alat_au*au_of_length_to_Angstrom
% Lattice_vector=alat*[1.000000   0.000000   0.000000
%                      -0.500000   0.866025   0.000000
%                      0.000000   0.000000   4.061244];

% NOTE: Lattice vector in the form 3x3, can be checked in output of scf (pw.x) calculations

                 
a1=Lattice_vector(1,:);
a2=Lattice_vector(2,:);
a3=Lattice_vector(3,:);    



% 
D_0_beta_ALL=zeros(3*Nb_atoms,3*Nb_atoms,Nq1,Nq2,Nq3);
R_0_beta_ALL=zeros(3,Nq1,Nq2,Nq3);
Position_cells=cell(Nq1,Nq2,Nq3);

for m1=1:Nq1
    for m2=1:Nq2
        for m3=1:Nq3
   
            % NOTE: important to define R(m1,m2,m3)
            % m1 m2 m3  force_constant(m1,m2,m3,i,j,na,nb) 
            % where m1=1,...,n1, m2=1,...,n2, m3=1,...,n3, define a lattice vector R:
            %    R = (m1-1)*tau1 + (m2-1)*tau2 + (m3-1)*tau3
            % (if m1-1 > n1/2 and so on, refold m1-1 => m1-1-n1 and so on)

               if (m1-1)>Nq1/2
                   position_along_a1=m1-1-Nq1;
               else
                  position_along_a1=m1-1;
               end
               if (m2-1)>Nq2/2
                   position_along_a2=m2-1-Nq2;
               else
                  position_along_a2=m2-1;
               end
               if (m3-1)>Nq3/2
                   position_along_a3=m3-1-Nq3;
               else
                  position_along_a3=m3-1;
               end               
         Position_cells{m1,m2,m3}=[position_along_a1,position_along_a2,position_along_a3];
         
         R_0_beta_ALL(:,m1,m2,m3)=position_along_a1*a1 + position_along_a2*a2 + position_along_a3*a3;
         
            D_0_beta=zeros(3*Nb_atoms,3*Nb_atoms);            
              for n1=1:Nb_atoms
                for n2=1:Nb_atoms
                    mass_n1=mass(n1);mass_n2=mass(n2);
                    
                    D_0_beta(1+3*(n1-1):3*n1,1+3*(n2-1):3*n2)= IFC_all_3x3{m1,m2,m3,n1,n2}/sqrt(mass_n1*mass_n2); % so, in unit of N/(kg.m)=1/s^2
                end
              end
   
            D_0_beta_ALL(:,:,m1,m2,m3)=D_0_beta; %unit of N/(kg.m)=1/s^2
        end
    end
end    


