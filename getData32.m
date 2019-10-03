%%%%%%%%%%%%%%%%%%%%%%%%%%getData32.m%%%%%%%%%%%%%%%%%
%
% Title : getData32.m
% Detail : Construct the consistent structure data from the 32 channel mat
% file to incoorporate with 16 channel data (getData.m) 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%getData32.m%%%%%%%%%%%%%%%%%%
%
%   Author  : Ain Chung
%   Date    : 09/18/2016
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%getData32.m%%%%%%%%%%%%%%%%%

function DataS=getData32(filename,correction_flag, ppi_flag)

if nargin<3
   ppi_flag=0; 
end

 Raw_data1=load(filename);
 if correction_flag==1
 Raw_data.data=Raw_data1.stimData(:,:,:);
 Raw_data.data(:,:,2:end+1)=Raw_data.data(:,:,1:end);
 
 elseif correction_flag==11
 Raw_data.data=Raw_data1.stimData(:,:,:);
  Raw_data.data(:,:,21:end+1)=Raw_data.data(:,:,20:end);

 
 elseif correction_flag==2
 Raw_data.data=Raw_data1.stimData(:,:,:);
 Raw_data.data(:,:,3:end+2)=Raw_data.data(:,:,1:end);

   elseif correction_flag==22 
 Raw_data.data=Raw_data1.stimData(:,:,:);
  Raw_data.data(:,:,2:end+1)=Raw_data.data(:,:,1:end); 
 Raw_data.data(:,:,35:end+1)=Raw_data.data(:,:,34:end);
 elseif correction_flag==-3
 Raw_data.data=Raw_data1.stimData(:,:,2:end-3);  
 elseif correction_flag==3
    
 Raw_data.data=Raw_data1.stimData(:,:,:);
 Raw_data.data (:,:,4:end+3)=Raw_data1.stimData(:,:,1:end);  
  
    elseif correction_flag==33 
 Raw_data.data=Raw_data1.stimData(:,:,:);
  Raw_data.data(:,:,2:end+1)=Raw_data.data(:,:,1:end); 
 Raw_data.data(:,:,21:end+1)=Raw_data.data(:,:,20:end);
 Raw_data.data(:,:,35:end+1)=Raw_data.data(:,:,34:end);
 elseif correction_flag==44
 Raw_data.data=Raw_data1.stimData(:,:,:);
 Raw_data.data (:,:,5:end+4)=Raw_data1.stimData(:,:,1:end);
 
     elseif correction_flag==444
 Raw_data.data=Raw_data1.stimData(:,:,:);
 Raw_data.data(:,:,2:end+1)=Raw_data.data(:,:,1:end); 
 Raw_data.data(:,:,21:end+1)=Raw_data.data(:,:,20:end);
 Raw_data.data(:,:,35:end+1)=Raw_data.data(:,:,34:end);
 Raw_data.data (:,:,end+1)=Raw_data1.stimData(:,:,end);

 
 elseif correction_flag==4
 Raw_data.data=Raw_data1.stimData(:,:,:);
 Raw_data.data (:,:,2:end+1)=Raw_data1.stimData(:,:,1:end);
 Raw_data.data (:,:,end+1)=Raw_data1.stimData(:,:,end);
 Raw_data.data (:,:,end+1)=Raw_data1.stimData(:,:,end);
 Raw_data.data (:,:,end+1)=Raw_data1.stimData(:,:,end);
 
 elseif correction_flag==55
 Raw_data.data=Raw_data1.stimData(:,:,:);
 Raw_data.data (:,:,6:end+5)=Raw_data.data(:,:,1:end);
 
 elseif correction_flag==16
 Raw_data.data=Raw_data1.stimData(:,:,2:end-16);
 
  elseif correction_flag==20
 Raw_data.data=Raw_data1.stimData(:,:,:);
 Raw_data.data (:,:,41:45)=Raw_data.data(:,:,36:end);
  Raw_data.data (:,:,46:55)=Raw_data.data(:,:,36:45);

 elseif correction_flag==5
 Raw_data.data=Raw_data1.stimData(:,:,2:end-2);
 elseif correction_flag==6
 Raw_data.data=Raw_data1.stimData(:,:,9:end-9); 
 elseif correction_flag==8
 Raw_data.data=Raw_data1.stimData(:,:,2:end-8);
 elseif correction_flag==-1
 Raw_data.data=Raw_data1.stimData(:,:,:);
 Raw_data.data(:,:,end+1)=Raw_data.data(:,:,end);
  elseif correction_flag==-2 
 Raw_data.data=Raw_data1.stimData(:,:,:);
 Raw_data.data(:,:,end+1)=Raw_data.data(:,:,end);
  Raw_data.data(:,:,end+1)=Raw_data.data(:,:,end);
  elseif correction_flag==10 % in the middle 175 missing
 Raw_data.data=Raw_data1.stimData(:,:,:);
 Raw_data.data(:,:,22:end+1)=Raw_data.data(:,:,21:end);
  elseif correction_flag==-5
 Raw_data.data=Raw_data1.stimData(:,:,:);
 Raw_data.data(:,:,end+1)=Raw_data.data(:,:,end);
  Raw_data.data(:,:,end+1)=Raw_data.data(:,:,end);
  Raw_data.data(:,:,end+1)=Raw_data.data(:,:,end);
  Raw_data.data(:,:,end+1)=Raw_data.data(:,:,end);
  Raw_data.data(:,:,end+1)=Raw_data.data(:,:,end);
 
   elseif correction_flag==-6
 Raw_data.data=Raw_data1.stimData(:,:,:);
 Raw_data.data(:,:,end+1)=Raw_data.data(:,:,end);
  Raw_data.data(:,:,end+1)=Raw_data.data(:,:,end);
  Raw_data.data(:,:,end+1)=Raw_data.data(:,:,end);
  Raw_data.data(:,:,end+1)=Raw_data.data(:,:,end);
  Raw_data.data(:,:,end+1)=Raw_data.data(:,:,end);
  Raw_data.data(:,:,end+1)=Raw_data.data(:,:,end);
 

 elseif correction_flag==11 % 03/13/2017: Dino: HM4x2_CNO_10mg_before_SAL_STIM_stim.mat
 Raw_data.data=Raw_data1.stimData(:,:,:);
 Raw_data.data(:,:,2:end+1)=Raw_data.data(:,:,1:end); 
 Raw_data.data(:,:,end+1:end+5)=Raw_data.data(:,:,end-4:end); % 275 intensity
 Raw_data.data(:,:,end+1:end+5)=Raw_data.data(:,:,end-4:end); % 300 intensity 
 elseif correction_flag==12 % 03/13/2017: Dino: HM4x2_CNO_10mg_before_SAL_STIM_stim.mat
 Raw_data.data=Raw_data1.stimData(:,:,:);
 Raw_data.data(:,:,end+1:end+5)=Raw_data.data(:,:,end-4:end); % 275 intensity
 Raw_data.data(:,:,end+1:end+5)=Raw_data.data(:,:,end-4:end); % 300 intensity 
  elseif correction_flag==13 % 03/13/2017: Dino: pp12 LTP experiment 225 to 300 X
 Raw_data.data=Raw_data1.stimData(:,:,:);
 Raw_data.data(:,:,end+1:end+3)=Raw_data.data(:,:,end-2:end); % 225 intensity
 Raw_data.data(:,:,end+1:end+3)=Raw_data.data(:,:,end-2:end); % 250 intensity 
 Raw_data.data(:,:,end+1:end+3)=Raw_data.data(:,:,end-2:end); % 275 intensity 
 Raw_data.data(:,:,end+1:end+3)=Raw_data.data(:,:,end-2:end); % 300 intensity 
 else
 Raw_data.data=Raw_data1.stimData(:,:,1:end);    % ERROR_stimulus_0_exp1 =1; correction term
 end
 
%% Initialization and Data aquisition
num_channel = size(Raw_data.data,1);            % number of channel : 16 channel
num_stimulus = size(Raw_data.data,3);           % number of total experiments.
num_time_sampled = size(Raw_data.data,2); % number of sampled data for each experiment
num_stimulus_different=size(Raw_data1.stimLevels,2); % number of different stimulus

sampling_freq= Raw_data1.stimFS;                                 % Sampling frequency 48kHz
milisecond = 1000;
measurement_time = num_time_sampled/sampling_freq;      % Experiment time on each experiment


stimulus_offset_time =0.005;                            % Inacitve period to find the offset level (4ms)
stimulus_offset_sample=sampling_freq*stimulus_offset_time; % number of samples in the offest time
time_sample =((1/sampling_freq:1/sampling_freq: measurement_time)-stimulus_offset_time)*milisecond; 

stimulus_valid_time_start = 0.0065;                            % Valid time after the stimulus (7ms)
stimulus_valid_time_start_sample=sampling_freq*stimulus_valid_time_start; 

stimulus_start=Raw_data1.stimLevels(1);                                       % Starting level of the stimulus
stimulus_final=Raw_data1.stimLevels(end);                                     % Final lelvel of the stimulus

if ppi_flag==0
    stimulus_step=Raw_data1.stimLevels(end)/(num_stimulus_different-1);           % Stimulus step increments
elseif ppi_flag==1
    stimulus_step=Raw_data1.stimLevels(end)/(num_stimulus_different);           % PPI Stimulus step increments 100/10 not 100/9
end
stimulus_experiment_number = num_stimulus/num_stimulus_different;  % number of experiment in each stimulus
stimulus_experiment_void_number =0;                                 % Keep every experiment
stimulus_experiment_effective_number= stimulus_experiment_number-stimulus_experiment_void_number; % number of VALID experiment in each stimulus

valid_data = zeros(num_channel,num_time_sampled, stimulus_experiment_effective_number*num_stimulus_different);

% Figure mode
% valid_data_offset_flag=1;

% Analysis mode
 valid_data_offset_flag=1;


for j=1:num_stimulus_different
    for k=(stimulus_experiment_effective_number)*(j-1)+1:(stimulus_experiment_effective_number)*j
         valid_data(:,:,k)=Raw_data.data(:,:,stimulus_experiment_void_number*j+k);
         % OFFSET compensation
         valid_data(:,:,k)=valid_data(:,:,k)-mean(valid_data(:,1:stimulus_offset_sample,k),2)*ones(1,num_time_sampled);
    end
end

if valid_data_offset_flag==0
    truncated_valid_data=valid_data(:,stimulus_valid_time_start_sample+1:end,:);
    
else
    truncated_valid_data=valid_data;
    
end
DataS = struct('data',truncated_valid_data,'stimulus_step',stimulus_step,'stimulus_start',stimulus_start,'stimulus_final',stimulus_final,...
                     'num_stimulus_different',num_stimulus_different,'stimulus_experiment_effective_number',stimulus_experiment_effective_number,...
                     'sampling_freq',sampling_freq,'stimulus_experiment_number',stimulus_experiment_number,'stimulus_offset_sample',stimulus_offset_sample,...
                     'stimulus_valid_time_start_sample',stimulus_valid_time_start_sample,'time_sample',time_sample,'num_channel',num_channel);
                 


end