clear all

%% Add the path of the tool:
addpath(genpath('/media/tranvt/Importance/DATA/RESEARCH/Works in time/2012-2015.In France.PhD/Dropbox/My papers/Drafts/2024.Paper2_DFT_phonon_transport/Matlab_codes_and_DFT_data/Matlab_codes/Read_DynMat_QE'));

%%

%cd('/media/tranvt/Importance/DATA/BACKUP/Onedrive/RESEARCH/IMPMC_Paris6_Postdoc_2019_2021/SIMULATIONS/Others/grupy-master/Example_Truong_practice/1.00')
cd('/media/tranvt/Importance/DATA/RESEARCH/Works in time/2012-2015.In France.PhD/Dropbox/My papers/Drafts/2024.Paper2_DFT_phonon_transport/Matlab_codes_and_DFT_data/Some_DFT_data/MoS2-monolayer/NEW_with_vc_relax_forces1Em4/phonon/Result')

nat=3;
amass=[95.94,32.065,32.065];%[207.2,127.60]; % amass(Pb)=207.2;amss(Te)=127.60; % Mass for MoS2: [95.94,32.065,32.065 ]

prefix='MoS2.dynmat';%'PbTe.dynmat';
% read dynamic matrix file: prefix.dynmat
[qpoints,DynaMat_full]=Read_DynMat_QE(prefix,nat,amass);

% qpoints is in unit of 2pi/alat
% DynaMat_full is in SI unit: N/(kg.m)=1/s^2

N_q=size(qpoints,1);
Omega=zeros(N_q,3*nat);
for i_q=1:N_q
     D_q = DynaMat_full{i_q};

    Omega2=eig(D_q);
    
    for n=1:3*nat
        if (real(Omega2(n))>=0)
        Omega(i_q,n)=sqrt(real(Omega2(n))); % unit 1/s=Hz
        else
        Omega(i_q,n)=-sqrt(-real(Omega2(n))); % unit 1/s=Hz 
        end        
    end
   Omega(i_q,:)=sort(Omega(i_q,:));     
    
end

%% 

cmm1_to_Hz=2.99793*1e10; % equivalent cmm1_to_secondm1

Freq=Omega/(2*pi)/cmm1_to_Hz;% convert to cm^-1


%% PLOT
%=============================================================================
npath=3;
segment_qpoints=[40,20,40,20];%Graphene: [40,20,40,20], {'G','K', 'M', 'G'};
Labels={'G','K', 'M', 'G'}; %{'G','X'};%{'L','G', 'X', 'W', 'G'};%{'G','X'};%
%Labels={'', '', '', '', ''};

Labels_position=zeros(1,npath+1);
Labels_position(1)=1;
for i=1:npath
Labels_position(i+1)=sum(segment_qpoints(1:i))+1;
end


Nq_path=N_q; % all q-point in the q-paths considered
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
legend(sprintf('%s',prefix))
set(gca,'Fontsize',24)
box on
%xlim([0,max(q_distance)])
ylim([-10,500])
title(sprintf('read from %s QE',prefix))



%% save figures
string1=sprintf('Phonon_bands_used_data_dynmat_file.png');
saveas(gcf,string1);% save figure
string2=sprintf('Phonon_bands_used_data_dynmat_file.fig');
saveas(gcf,string2);% save figure