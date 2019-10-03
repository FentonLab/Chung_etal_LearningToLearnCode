%%%%%%%%%%%%%%%%%%%%%%%%%%CSD_analysis_AUC_32ch_HM4_load.m%%%%%%%%%%%%%%%%%
%
% Title : CSD_analysis_AUC_32ch_HM4.m
% Detail : CSD data analysis of 32channel experiments
% copied from EEPP_Analysis_32channel.m
% ANIMALS : HM4
% Detail : Only one animal loading. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%CSD_analysis_AUC_32ch_HM4_load.m%%%%%%%%%%%%%%%%%%
%
%   Author  : Patrick H
%   Date    : 03/16/2017
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%CSD_analysis_AUC_32ch_HM4_load.m%%%%%%%%%%%%%%%%%

%% Clear the variable and Load the data
clear all; 

 



%% Initialization and Data aquisition
channel32=1; % ACTIVATE CHANNEL32
Total_number_session_truncated=10; % First 10 session
unit_scale = 1e-3; % A/m^3 -> muA/mm^3


Total_num_animal=1;



    filename={ 'PPx_STIM_PRETRAIN_BEFORE_stim',...
                'PPx_STIM_PRETRAIN_AFTER_stim',...
                'PPx_STIM_17_TR1_BEFORE_stim',...
                'PPx_STIM_21_TR1_AFTER_stim',...
                'PPx_STIM_23_TR2_BEFORE_stim',...
                'PPx_STIM_26_TR2_AFTER_stim',...
                'PPx_STIM_28_TR3_BEFORE_stim',...
                'PPx_STIM_31_TR3_AFTER_stim',...
                'PPx_STIM_56_RET_BEFORE_stim',...
                'PPx_STIM_59_RET_AFTER_stim'
              };


          correction_flag=zeros(1,length(filename));
           %correction_flag(9)=1;

          Intensity=250;

Total_number_session =length(filename);

if channel32==1
    for i=1:Total_number_session
        DataS(i,1)=getData32(filename{i},correction_flag(i));
    end
else
    for i=1:Total_number_session
        DataS(i,1)=getData(filename{i},correction_flag(i));
    end
    
end


%% CSD using CSDplotter based on Dino code
Intensity = 300;
%Intensity_step=DataS(1,1).stimulus_step;
Intensity_step=25;
Intensity_index=Intensity/Intensity_step+1;
Effective_intensity_max32=DataS(1,1).num_stimulus_different;
for j=1:Effective_intensity_max32
    for k=1:Total_num_animal
        figure;
 
        CSDToT(:,:,:,k,j)=plotCSDPlotter(DataS(:,k),Total_number_session_truncated,(j-1)*Intensity_step,channel32);
        % 
        % --Inveseback to the usual ordering
        
    end
end

CSDToT=CSDToT*unit_scale;
CSD=CSDToT(:,:,:,:,Intensity_index);


%%
save('CSDToTAUC_031319_PPx_training','CSDToT','Total_number_session_truncated', 'DataS')
disp('CSD data for 32ch saved')

%% Validity check CSD time plot for one channel
figure;
Intensity = 500;
ch_test=11;
session_test=1;
animal_test =1; % 1: only one animal
sampling_freq_test= DataS(session_test,animal_test).sampling_freq;   
time_test=(1:size(CSD(ch_test,:,session_test,animal_test),2))/sampling_freq_test/unit_scale;
plot(time_test,CSD(ch_test,:,session_test,animal_test))
