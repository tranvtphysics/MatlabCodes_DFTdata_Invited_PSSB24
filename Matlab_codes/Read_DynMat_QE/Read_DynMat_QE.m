function [qpoints,DynaMat_full]=Read_DynMat_QE(filename_dynmat,nat,amass)
% The copyright by Van-Truong TRAN.
% Code version: V2.0 (01/2022)
% change compared V1.0: Scanning from the line 'Dynamical  Matrix in
% cartesian axes', so the code can read also dynamical matrices from
% individual files .dyn1, .dyn2, ...


% READING Dynamical matrix file generated by QE calculculation
% filename_dynmat: name of the system: file data is "prefix.dynmat". Format: string
% nat: number of atoms in the unit cell in QE calculation. Format: integer
% amass: atomic mass of each kind of atoms. Format: real

%%
% prefix='PbTe';
% nat=2;
% amass=[1,1];

mass=amass*1.660540199E-27;% convert masses to kg


Ryd_to_Nm=2.1798741E-18;% 1 Ryd=2.1798741E-18 N*m;
Bohr_to_m=5.29177e-11; %1bohr=5.29177e-11 m
%% Reading the file "prefix.dynmat"
% the content of the file looks like this:
%     Dynamical  Matrix in cartesian axes
% 
%      q = (    0.000000000   0.300000000   0.000000000 ) 
% 
%     1    1
%   0.02451659  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000
%   0.00000000  0.00000000    0.09561959  0.00000000    0.00000000  0.00000000
%   0.00000000  0.00000000    0.00000000  0.00000000    0.02451659  0.00000000
%     1    2
%  -0.01467754  0.02020189    0.00000000  0.00000000    0.00000000  0.00000000
%   0.00000000  0.00000000   -0.03361601  0.04626847    0.00000000  0.00000000
%   0.00000000  0.00000000    0.00000000  0.00000000   -0.01467754  0.02020189
%     2    1
%  -0.01467754 -0.02020189    0.00000000  0.00000000    0.00000000  0.00000000
%   0.00000000  0.00000000   -0.03361601 -0.04626847    0.00000000  0.00000000
%   0.00000000  0.00000000    0.00000000  0.00000000   -0.01467754 -0.02020189
%     2    2
%   0.03236203  0.00000000    0.00000000  0.00000000    0.00000000  0.00000000
%   0.00000000  0.00000000    0.10257058  0.00000000    0.00000000  0.00000000
%   0.00000000  0.00000000    0.00000000  0.00000000    0.03236203  0.00000000


file_data=sprintf('%s',filename_dynmat);

fileID = fopen(file_data,'r');

[status1,cmdout1]=system(sprintf('wc -l < %s',file_data)); % Counting number of data lines in the data file
Nlines=str2num(cmdout1);

number_lines_scaned=Nlines; % Look at txt file to see how many lines need to be scaned
Intro = textscan(fileID,'%s',number_lines_scaned,'Delimiter','\n');
IF=Intro{1};

N_q=0;
for i=1:Nlines
    
    Line_i=IF{i};
    Line_i=strtrim(Line_i);
    
    if length(Line_i)>5 && strcmp(Line_i,'Dynamical  Matrix in cartesian axes')==1
        
        i_q=i+2; % index the of line start with q =
        
        Line_q=strtrim(IF{i_q});
        %Line_i
        N_q=N_q+1; % counting number of q-points
        qpoints(N_q,1:3)=str2num(Line_q(6:end-1)); % q-points in the .dynmat file are in the unit of 2pi/alat 
        
        % data for this q-point start from the line iq+2 to i + 2 + (4 *nat*nat)
          start_line_data = i_q + 2;
          end_line_data = i_q + 2 + (4*nat*nat)-1;
          
          data_q=cell(4 * nat * nat);
          for j=1:(4 * nat * nat)
              data_q{j}=IF{i_q+2+j-1};
          end
          
          
          for j=1:4*nat*nat
              if mod(j-1,4)==0
                  n1_n2=str2num(data_q{j});
                  n1=n1_n2(1); n2=n1_n2(2);
                                    
                  Row1=str2num(data_q{j+1});
                  Row2=str2num(data_q{j+2});
                  Row3=str2num(data_q{j+3});
                  
                    Dn1n2(1,1)=Row1(1)+1i*Row1(2);Dn1n2(1,2)=Row1(3)+1i*Row1(4); Dn1n2(1,3)=Row1(5)+1i*Row1(6) ; 

                    Dn1n2(2,1)=Row2(1)+1i*Row2(2);Dn1n2(2,2)=Row2(3)+1i*Row2(4); Dn1n2(2,3)=Row2(5)+1i*Row2(6) ; 

                    Dn1n2(3,1)=Row3(1)+1i*Row3(2);Dn1n2(3,2)=Row3(3)+1i*Row3(4); Dn1n2(3,3)=Row3(5)+1i*Row3(6) ;
                  
                    D{N_q,n1,n2}=Dn1n2; % NOTE: Unit of Dynamical matrix in the .dynmat file is in the unit of Ryd/bohr^2. This is acctually Force-constant matrix rather than dynamical matrix 

                        %n1_n2
                        %D{N_q,n1,n2}
              end
          end
                            
          
          
    end
    
    
    %Line_i=str2num(IF{ii+2});
   
     
end


fclose(fileID) ;
%% Constructing the full Dynamical matrix for the whole system (at each q-point)
DynaMat_full=cell(N_q); % at each q-point, the size of the DynMat_ful is 3natX3nat. Convert to SI unit


for i_q=1:N_q
    for n1=1:nat
        for n2=1:nat
            mass_n1=mass(n1);mass_n2=mass(n2);
            
            D_q_n1n2=D{i_q,n1,n2}; % matrix 3x3. unit Ryd/bohr^2 
            D_q_n1n2=D_q_n1n2*Ryd_to_Nm/Bohr_to_m^2; % so, unit N/m
            D_q(1+3*(n1-1):3*n1,1+3*(n2-1):3*n2)= D_q_n1n2/sqrt(mass_n1*mass_n2); % so, in unit of N/(kg.m)=1/s^2
        end
    end
DynaMat_full{i_q}=D_q;
end

