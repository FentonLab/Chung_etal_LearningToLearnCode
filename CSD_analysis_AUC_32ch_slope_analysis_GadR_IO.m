f%%%%%%%%%%%%%%%%%%%%%%%%%%CSD_analysis_AUC_32ch_slope_analysis_PPI.m%%%%%%%%%%%%%%%%%
%
% Title : CSD_analysis_AUC_32ch_slope_analysis.m
% Detail : CSD data Temporal (Slope/Depth measure) analysis of 32channel experiments 
% copied from EEPP_Analysis_32channel.m
% ANIMALS : Single animal 
% Detail : read from CSD_analysis_AUC_32ch_load.m file and analyze.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%CSD_analysis_AUC_32ch_HM4_slope_analysis.m%%%%%%%%%%%%%%%%%%
%
%   Author  : Ain Chung
%   Date    : Nov 2018
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%CSD_analysis_AUC_32ch_HM4_slope_analysis.m%%%%%%%%%%%%%%%%%

%% Clear the variable and Load the data
clear all; 


%% Initialization and Data aquisition

channel32=1; % ACTIVATE CHANNEL32


%------------------- GADR --------------
%strCSD=load('CSDToTAUC_20181003_GadR1_IO.mat');
%strCSD=load('CSDToTAUC20180508_GadR1_IO_optic.mat');


%strCSD=load('CSDToTAUC_20181003_GadR3_wo_optic.mat');
%strCSD=load('CSDToTAUC_20181003_GadR3_IO_w_feedbackoptic.mat');

%strCSD=load('CSDToTAUC_20181003_GadR4_wo_optic.mat');
%strCSD=load('CSDToTAUC_20181003_GadR4_IO_w_feedbackoptic.mat');

%strCSD=load('CSDToTAUC_20181003_GadR5_wo_optic.mat');
%strCSD=load('CSDToTAUC_20181003_GadR5_IO_w_feedbackoptic.mat');

%strCSD=load('CSDToTAUC_20181003_GadR9_wo_optic.mat');
%strCSD=load('CSDToTAUC_20181003_GadR9_IO_w_feedbackoptic.mat');

 
%strCSD=load('CSDToTAUC_20181003_GadR10_wo_optic.mat');
%strCSD=load('CSDToTAUC_20181003_GadR10_IO_w_feedbackoptic.mat');

%strCSD=load('CSDToTAUC_20181003_GadR11_wo_optic.mat');
%strCSD=load('CSDToTAUC_20181003_GadR11_IO_w_feedbackoptic.mat');
%strCSD=load('CSDToTAUC_20181003_GadR11_IOcurve_w_optic.mat');

%strCSD=load('CSDToTAUC_20181003_GadR12_wo_optic.mat');
%strCSD=load('CSDToTAUC_20181003_GadR12_IO_w_feedbackoptic.mat');
%strCSD=load('CSDToTAUC_20181003_GadR12_IOcurve_w_optic.mat');


%strCSD=load('CSDToTAUC_20181003_GadR13_wo_optic.mat');
%strCSD=load('CSDToTAUC_20181003_GadR13_IO_w_feedbackoptic.mat');
%strCSD=load('CSDToTAUC_20181003_GadR13_IOcurve_w_optic.mat');


%strCSD=load('CSDToTAUC_20181003_GadR14_wo_optic.mat');
%strCSD=load('CSDToTAUC_20181003_GadR14_IO_w_feedbackoptic.mat');
%strCSD=load('CSDToTAUC_20181003_GadR14_IOcurve_w_optic.mat');

%strCSD=load('CSDToTAUC_20181003_GadR16_wo_optic.mat');
%strCSD=load('CSDToTAUC_20181003_GadR16_IO_w_feedbackoptic.mat');
%strCSD=load('CSDToTAUC_20181003_GadR16_IOcurve_w_optic.mat');


%strCSD=load('CSDToTAUC_20181003_GadR20_wo_optic.mat');
%strCSD=load('CSDToTAUC_20181003_GadR20_IO_w_feedbackoptic.mat');
%strCSD=load('CSDToTAUC_20181003_GadR20_IOcurve_w_optic.mat');


%strCSD=load('CSDToTAUC_20181003_GadR22_wo_optic.mat');
%strCSD=load('CSDToTAUC_20181003_GadR22_IO_w_feedbackoptic.mat');%
%strCSD=load('CSDToTAUC_20181003_GadR22_IOcurve_w_optic.mat');



%strCSD=load('CSDToTAUC_20181003_GadR23_wo_optic.mat');
%strCSD=load('CSDToTAUC_20181003_GadR23_IO_w_feedbackoptic.mat');
%strCSD=load('CSDToTAUC_20181003_GadR23_IOcurve_w_optic.mat');


%strCSD=load('CSDToTAUC_20181003_GadR24_wo_optic.mat');
%strCSD=load('CSDToTAUC_20181003_GadR24_IO_w_feedbackoptic.mat');
%strCSD=load('CSDToTAUC_20181003_GadR24_IOcurve_w_optic.mat');

%strCSD=load('CSDToTAUC_20181003_GadR25_wo_optic.mat');
%strCSD=load('CSDToTAUC_20181003_GadR25_IO_w_feedbackoptic.mat');
%strCSD=load('CSDToTAUC_20181003_GadR25_IOcurve_w_optic.mat');

strCSD=load('CSDToTAUC_20181003_GadR26_wo_optic.mat');
%strCSD=load('CSDToTAUC_20181003_GadR26_IO_w_feedbackoptic.mat');
%strCSD=load('CSDToTAUC_20181003_GadR26_IOcurve_w_optic.mat');

%strCSD=load('CSDToTAUC_20181003_GadR27_wo_optic.mat');
%strCSD=load('CSDToTAUC_20181003_GadR27_IO_w_feedbackoptic.mat');
%strCSD=load('CSDToTAUC_20181003_GadR27_IOcurve_w_optic.mat');

%strCSD=load('CSDToTAUC_20181003_GadR28_wo_optic.mat');
%strCSD=load('CSDToTAUC_20181003_GadR28_IO_w_feedbackoptic.mat');
%strCSD=load('CSDToTAUC_20181003_GadR28_IOcurve_w_optic.mat');


Total_number_session_truncated=strCSD.Total_number_session_truncated; % First 10 session

Total_num_section = 6; % Anatomical section (CA1/ SLM/ DGU_dend/ DGU_Soma/ DGL_Soma/ DGL_Dend)

CSDToT=strCSD.CSDToT;
DataS=strCSD.DataS;
num_channel=DataS(1,1).num_channel;
unit_scale = 1e-3; % A/m^3 -> muA/mm^3

Total_num_animal=1;
Effective_intensity_max32=DataS(1,1).num_stimulus_different;

if channel32==0
    Effective_intensity_max=100;
    Effective_intensity_max32=floor(Effective_intensity_max/DataS(1,1).stimulus_step)+1;
end
Total_num_channel_analysis32=DataS(1,1).num_channel;
Intensity_step=DataS(1,1).stimulus_step;

% Sampling freq and time_test for AUC integration in this code.
% REMEMBER that it is integrated with the unit milisecond. 
sampling_freq_auc= DataS(1,1).sampling_freq;   
time_test_auc=(1:size(CSDToT(1,:,1),2))/sampling_freq_auc/unit_scale;

%Validity check CSD time plot for one channel
figure;
Intensity = 400;
Intensity_step=DataS(1,1).stimulus_step;
Intensity_index=Intensity/Intensity_step;
ch_test=5;
session_test=1;
animal_test =1; % 1: only one animal
CSD=strCSD.CSDToT(:,:,:,:,Intensity_index);
sampling_freq_test= DataS(session_test,animal_test).sampling_freq;   
time_test=(1:size(CSD(:,:,session_test,animal_test),2))/sampling_freq_test/unit_scale;

plot(time_test,CSD(ch_test,:,session_test,animal_test));
hold on
plot(time_test,CSD(ch_test,:,session_test,animal_test));

Intensity = 500;

Intensity_index=Intensity/Intensity_step;
figure;
for x=1:Total_number_session_truncated
    subplot (1,Total_number_session_truncated,x);
    
    session_test=x;

    for j =1:DataS(session_test,1).num_channel
            plot(time_test,CSDToT(j,:,session_test,animal_test,Intensity_index)*2+(num_channel-j)*300)
%           plotNoStimulusErrorbar(DataS(before),getAverage(DataS(before),j,Intensity,channel32)*1+(DataS(before).num_channel-j)*0.0090,j);

        hold on; 
        xlim([0,50]);
        ylim([-1000,10000]);
    end
    
        hold off 
end

%%
figure;

for w = 1:10
    
    for j =1:DataS(1,1).num_channel
            plot(time_test,CSDToT(j,:,1,animal_test,w)*4+(num_channel-j)*300);
%           plotNoStimulusErrorbar(DataS(before),getAverage(DataS(before),j,Intensity,channel32)*1+(DataS(before).num_channel-j)*0.0090,j);

        hold on; 
        xlim([0,100]);
      % ylim([-0.002,0.27]);
    end
    
    hold on

   
    
end



%% 2. Temporal analysis 0) (INITIALIZATION) OR (LOAD DATA) for slope analysis

choice = questdlg('Is it your first time for this animal?', ...
	'Initialization', ...
	'Yes. First time. I need to initialize','No. I want to load the existing file','I do not know','I do not know');
% Handle response
switch choice
    case 'Yes. First time. I need to initialize'
        disp(' Initializing')
        initialization_flag = 1;
    case 'No. I want to load the existing file'
        disp('Loading the data')
        initialization_flag = 2;
    case 'I do not know'
        error('Make sure of what you are doing');
    case ''
        error('You need to select one')
end

if initialization_flag==1
    Slope_data=zeros(Total_number_session_truncated,Effective_intensity_max32,num_channel,Total_num_animal);
    Slope_data_max=zeros(Total_number_session_truncated,Effective_intensity_max32,num_channel,Total_num_animal);
    Slope_data_curv_max=zeros(Total_number_session_truncated,Effective_intensity_max32,num_channel,Total_num_animal);
    
    Depth_data=zeros(Total_number_session_truncated,Effective_intensity_max32,num_channel,Total_num_animal);
    Time_data=zeros(Total_number_session_truncated,Effective_intensity_max32,num_channel,Total_num_animal,4);
    Section_Index_data = zeros(Total_number_session_truncated,Effective_intensity_max32,num_channel,Total_num_animal);
    PopSpike_Channel_data = zeros(Total_number_session_truncated,Effective_intensity_max32,Total_num_animal, Total_num_section);
    
   CSD_AUC_DGU_data=zeros(Total_number_session_truncated, Effective_intensity_max32, num_channel, Total_num_animal);
   CSD_Integral_data=zeros(Total_number_session_truncated, Effective_intensity_max32, num_channel, Total_num_animal);
   
else
    Slope_load=load('CSDSlope_112118_GadR26_IO.mat');
    Depth_load=load('CSDDepth_112118_GadR26_IO.mat');
    Time_load=load('CSDTIME_112118_GadR26_IO.mat');
    AUC_load=load('CSDAUC_112118_GadR26_IO.mat');
    
    %-- 11/6/17 new data structure
    Section_Index_load=load('CSDSecIndex_112118_GadR26_IO.mat');
    PopSpike_Channel_load=load('CSDPopChannel_112118_GadR26_IO.mat');
    AUC_DGU_load=load('CSDAUC_DGU_112118_GadR26_IO.mat');
% %     
    Slope_data = Slope_load.Slope_data;
    Slope_data_max = Slope_load.Slope_data_max;
    Slope_data_curv_max = Slope_load.Slope_data_curv_max;
    
    Depth_data = Depth_load.Depth_data;
    Time_data = Time_load.Time_data;
    CSD_Integral_data = AUC_load.CSD_Integral_data;
    
    %-- 11/6/17 new data strucutre
    Section_Index_data = Section_Index_load.Section_Index_data;
    PopSpike_Channel_data = PopSpike_Channel_load.PopSpike_Channel_data;
    CSD_AUC_DGU_data=AUC_DGU_load.CSD_AUC_DGU_data;

end


%% 2. Temporal analysis 1) Get critical time manually. (GET)
% JUST FOR THE TEST WARNING!!!!!!!!!%1) Get critical time manually. (INITIALIZATION) OR (LOAD DATA) for slope analysis
% figure(120120);
global ResetTest t1 t2 t3 t4;

ResetTest= 0;
% close 120120;
milisecond=1/unit_scale;    %%% ATTENTION: DataS.time_sample is in the milisecond unit

test_animal=1;
Intensity_step=DataS(1,1).stimulus_step;
% Start number fix : R emove the first experiment. In result: 7 out of 8
start_experiment= DataS(1,1).stimulus_experiment_number-DataS(1,1).stimulus_experiment_effective_number+1;



%%%%%%%%%%%%%%%%%******* MODE SELECTION ***********************************

flag_manual=1;    % 2:Semi-Manual 1: Manual 0: Automatic

flag_second_peak = 0;    %0: inactive, 1:active (second peak)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-- Change here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---Section division---
ch_SLM_start = 2;
ch_DGU_neg_start = 4;          % DG upper level starting channel
ch_DGU_pos_start = 7;       % DG upper level - Dendritic ending channel
ch_DGL_pos_start = 10;            %  DG Lower level soma starting channel
ch_DGL_neg_start = 12;            %  DG lower level dendritic starting channel

%---Population Spike reference channel
PopSpike_reference_CA1_ch=1;
PopSpike_reference_DGU_ch= 7;
PopSpike_reference_DGL_ch= 11;
%----------------------

Intensity=0;    
test_session=1;           % Session number (before after...etc)
test_ch=11;                % Channel number %pp5: 20 24 29 %pp7: 18 22 26 %pp12: 19 24 30 %pp13: 
min_max_flag =1;        %1:in semi-mode, detact popspike automatically
 

t1=14;
t2=t1+0.5; 
t3=10.3;
t4=11.1;

%%%%%%%%%%%%%%% FIRST PEAK TIMING %%%%%%%%%%%%%%%
start_window_time = (t1+0.1)/milisecond;
end_window_time  =start_window_time+5/milisecond;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


start_session =1;
end_session = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Intensity_index=Intensity/Intensity_step+1;
CSD=CSDToT(:,:,:,:,Intensity_index);


sampling_freq_test= DataS(test_session,test_animal).sampling_freq;   
time_test=(1:size(CSD(test_ch,:,test_session),2))/sampling_freq_test/unit_scale;



%%%%%%%%%%%%%%%%%******* Plot Property ***********************************

xlim_min=0;
xlim_max=150;              % 2:Semi-Manual 1: Manual 0: Automatic

ylim_min=0;
ylim_max=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%---Test property display

disp('%%-----Test Property-----%%');

fprintf('Time : %s\n',datestr(datetime('now')))
if flag_manual==1
   (fprintf('Mode : Manual\n'));
   (fprintf('Current Session : %d\n', test_session));
elseif flag_manual==2
    (fprintf('Mode : Semi-manual\n'));
end
(fprintf('Current Intensity : %d\n', Intensity));
(fprintf('Current channel : %d\n', test_ch));
if (test_ch>=ch_DGU_neg_start) && (test_ch<ch_DGU_pos_start)
    (fprintf('Current Section : DG Upper Dendritic (%d, %d)\n',ch_DGU_neg_start,ch_DGU_pos_start-1));
    test_section= 3;
    
elseif (test_ch>=ch_DGU_pos_start) && (test_ch<ch_DGL_pos_start)
    (fprintf('Current Section : DG Upper Soma (%d, %d)\n',ch_DGU_pos_start,ch_DGL_pos_start-1));
    test_section= 4;

elseif (test_ch>=ch_SLM_start)&&(test_ch<ch_DGU_neg_start)
    fprintf('Current Section : SLM (%d, %d)\n', ch_SLM_start, ch_DGU_neg_start-1)
    test_section= 2;

elseif test_ch<ch_SLM_start
    fprintf('Current Section : CA1 (%d, %d)\n', 1, ch_SLM_start-1)
    test_section= 1;

elseif (test_ch>=ch_DGL_pos_start)&&(test_ch<ch_DGL_neg_start)
    fprintf('Current Section : DG Lower Soma (%d, %d)\n', ch_DGL_pos_start, ch_DGL_neg_start-1)
    test_section= 5;
    
elseif (test_ch>=ch_DGL_neg_start)
    fprintf('Current Section : DG Lower Dendritic (%d, %d)\n', ch_DGL_neg_start, num_channel)
    test_section= 6;
end

%----------------------------

  if flag_manual==1
 
 
      Data=CSD(test_ch,:,test_session,test_animal);
            
        if min_max_flag==1
            start_window_pts=floor(sampling_freq_test*start_window_time);
            end_window_pts=floor(sampling_freq_test*end_window_time);
            [~,Imax]=max(Data(start_window_pts:end_window_pts));
            [~,Imin]=min(Data(start_window_pts:end_window_pts));
            Imax=Imax+start_window_pts-1;
            Imin=Imin+start_window_pts-1;

        end
        figure(5112015);
        set(gcf, 'Position', [100, 100, 1000, 1000])
        hold on;
        plot(time_test,Data,'bo-');
        
        plot(time_test,zeros(1,length(time_test)),'k-');
        
        xlim([xlim_min, xlim_max])
        
        if min_max_flag==1
            
            t3 = time_test(Imin);
            tempcritical_pts =round(sampling_freq_test*t3/milisecond); 
            temp_data=CSD(test_ch,tempcritical_pts,test_session,test_animal);
            plot(tempcritical_pts*milisecond/sampling_freq_test,temp_data,'ro')
        
            t4 = time_test(Imax);
            tempcritical_pts =round(sampling_freq_test*t4/milisecond); 
            temp_data=CSD(test_ch,tempcritical_pts,test_session,test_animal);
            plot(tempcritical_pts*milisecond/sampling_freq_test,temp_data,'ro')
        
            h = zoom;
            set(h,'Motion','horizontal','Enable','on');

            pause;

            [t1, ~] = ginput(1);
            t1=round(t1,1);
            tempcritical_pts =round(sampling_freq_test*t1/milisecond);
            temp_data=CSD(test_ch,tempcritical_pts,test_session,test_animal);
            plot(tempcritical_pts*milisecond/sampling_freq_test,temp_data,'ro')

            [t2, ~] = ginput(1);
            t2=round(t2,1);
            tempcritical_pts =round(sampling_freq_test*t2/milisecond); 
            temp_data=CSD(test_ch,tempcritical_pts,test_session,test_animal);
            plot(tempcritical_pts*milisecond/sampling_freq_test,temp_data,'ro')

            critical_time_pts=[t1;t2;t3;t4]/milisecond;
            critical_pts=round(sampling_freq_test*critical_time_pts);
            temp_data=CSD(test_ch,critical_pts,test_session,test_animal);
        else
        

        h = zoom;
        set(h,'Motion','horizontal','Enable','on');

        pause;
     
        [t1, ~] = ginput(1);
        t1=round(t1,1);
        tempcritical_pts =round(sampling_freq_test*t1/milisecond);
        temp_data=CSD(test_ch,tempcritical_pts,test_session,test_animal);
        plot(tempcritical_pts*milisecond/sampling_freq_test,temp_data,'ro')

        [t2, ~] = ginput(1);
        t2=round(t2,1);
        tempcritical_pts =round(sampling_freq_test*t2/milisecond); 
        temp_data=CSD(test_ch,tempcritical_pts,test_session,test_animal);
        plot(tempcritical_pts*milisecond/sampling_freq_test,temp_data,'ro')
        
        [t3, ~] = ginput(1);
        t3=round(t3,1);
        tempcritical_pts =round(sampling_freq_test*t3/milisecond); 
        temp_data=CSD(test_ch,tempcritical_pts,test_session,test_animal);
        plot(tempcritical_pts*milisecond/sampling_freq_test,temp_data,'ro')
        
        [t4, ~] = ginput(1);
        t4=round(t4,1);
        tempcritical_pts =round(sampling_freq_test*t4/milisecond); 
        temp_data=CSD(test_ch,tempcritical_pts,test_session,test_animal);
        plot(tempcritical_pts*milisecond/sampling_freq_test,temp_data,'ro')
        
        critical_time_pts=[t1;t2;t3;t4]/milisecond;
        critical_pts=round(sampling_freq_test*critical_time_pts);
        temp_data=CSD(test_ch,critical_pts,test_session,test_animal);

        end
          plot(critical_time_pts*milisecond,temp_data,'r*')
  
        btn = uicontrol('Style', 'pushbutton', 'String', 'Close',...
        'Position', [20 20 650 20],...
        'CallBack','close(gcf)');
        %-- Define the reference channel for each cell (CA1, DGU, DGL)
        if test_section==1 || test_section==2
            PopSpike_reference_ch=PopSpike_reference_CA1_ch;
        elseif test_section==3 || test_section==4
            PopSpike_reference_ch=PopSpike_reference_CA1_ch;
        elseif test_section==5 || test_section==6
            PopSpike_reference_ch=PopSpike_reference_CA1_ch;
        end
        
        %-- Mark the PopSpike reference time t3
        if PopSpike_reference_ch==test_ch
             plot(critical_time_pts(3)*milisecond,temp_data(3),'ro','MarkerFaceColor','r');
        end 
      
        hold off;
        
       fprintf('t1, t2, t3, t4 : (%d, %d, %d, %d)\n', t1,t2,t3,t4)
       fprintf('PopSpike reference Channel : %d\n', PopSpike_reference_ch)
     
        % Max slope and Max curvature 
            critical_pts_index = critical_pts(1):critical_pts(2);
            critical_data = CSD(test_ch,critical_pts_index,test_session,test_animal);
            critical_slope_index = diff(critical_data);
            [der,max_slope_candi,min_slope_candi]=slope_max_min(critical_data,sampling_freq_test);
            [der_curv,max_slope_curv_candi,min_slope_curv_candi]=slope_max_min(der,sampling_freq_test);

            if abs(max_slope_candi)>abs(min_slope_candi)
            max_slope=max_slope_candi;
            else
            max_slope=min_slope_candi;
            end
            if abs(max_slope_curv_candi)>abs(min_slope_curv_candi)
            max_slope_curv=max_slope_curv_candi;
            else
            max_slope_curv=min_slope_curv_candi;
            end

        % AUC (integral)

        Effective_time_test=time_test(critical_pts(1):critical_pts(2));
        Effective_temp_data=CSD(test_ch,critical_pts(1):critical_pts(2),test_session,test_animal);

        CSD_integral=@(t) interp1(Effective_time_test,Effective_temp_data,t,'pchip');
        auc=integral(CSD_integral,Effective_time_test(1),Effective_time_test(end));


        %--------------------------------------------------------------------------------
        % Save MANUALLY
        [slope, ~]= getSlopeDepthThreepoint(critical_time_pts(1:3,1),temp_data(1:3));
        Slope_data(test_session,Intensity_index,test_ch,test_animal)=slope(1);
        Slope_data_max(test_session,Intensity_index,test_ch,test_animal)=max_slope;
        Slope_data_curv_max(test_session,Intensity_index,test_ch,test_animal)=max_slope_curv;
        CSD_Integral_data(test_session,Intensity_index,test_ch,test_animal)=auc;

        [~, depth]= getSlopeDepthThreepoint(critical_time_pts(2:4,1),temp_data(2:4));
        Depth_data(test_session,Intensity_index,test_ch,test_animal)=depth;
        Time_data(test_session,Intensity_index,test_ch,test_animal,:)=critical_time_pts;
        
            %-- CA1: 1 SLM :2 DGU-Dend :3 DGU-Soma :4 DGL-Soma :5 DGL-Dend :6
        Section_Index_data(test_session,Intensity_index,1:ch_SLM_start-1,test_animal)= 1;
        Section_Index_data(test_session,Intensity_index,ch_SLM_start:ch_DGU_neg_start-1,test_animal)= 2;
        Section_Index_data(test_session,Intensity_index,ch_DGU_neg_start:ch_DGU_pos_start-1,test_animal)= 3;
        Section_Index_data(test_session,Intensity_index,ch_DGU_pos_start:ch_DGL_pos_start-1,test_animal)= 4;
        Section_Index_data(test_session,Intensity_index,ch_DGL_pos_start:ch_DGL_neg_start-1,test_animal)= 5;
        Section_Index_data(test_session,Intensity_index,ch_DGL_neg_start:end,test_animal)= 6;    
            %-- Popsulation Spike reference Channel
        PopSpike_Channel_data(test_session,Intensity_index,test_animal,1)=PopSpike_reference_CA1_ch;
        PopSpike_Channel_data(test_session,Intensity_index,test_animal,2)=PopSpike_reference_CA1_ch;
        PopSpike_Channel_data(test_session,Intensity_index,test_animal,3)=PopSpike_reference_DGU_ch;
        PopSpike_Channel_data(test_session,Intensity_index,test_animal,4)=PopSpike_reference_DGU_ch;
        PopSpike_Channel_data(test_session,Intensity_index,test_animal,5)=PopSpike_reference_DGL_ch;
        PopSpike_Channel_data(test_session,Intensity_index,test_animal,6)=PopSpike_reference_DGL_ch;
    
  elseif flag_manual==2

 %%%%%*************SEMI MANULA MODE 4 TIME POINTS ****************************

    figure;
     xlim([xlim_min, xlim_max])
    figure_gap =60;
    
    for kk=start_session:end_session
         
        test_session=kk;
        
     
          Data=CSD(test_ch,:,test_session,test_animal);
             
        if min_max_flag==1
            
            start_window_pts=floor(sampling_freq_test*start_window_time);
            end_window_pts=floor(sampling_freq_test*end_window_time);
            [~,Imax]=max(Data(start_window_pts:end_window_pts));
            [~,Imin]=min(Data(start_window_pts:end_window_pts));
            Imax=Imax+start_window_pts-1;
            Imin=Imin+start_window_pts-1;

        end
        
        plot(time_test,Data-(kk-1)*figure_gap);
        
        hold on;    

        if min_max_flag==1
            t3 = time_test(Imin);
            t4 = time_test(Imax);
            critical_time_pts=[t1;t2;t3;t4]/milisecond;
            critical_pts=round(sampling_freq_test*critical_time_pts);
            temp_data=CSD(test_ch,critical_pts,test_session,test_animal);
        else

            critical_time_pts=[t1;t2;t3;t4]/milisecond;
            critical_pts=round(sampling_freq_test*critical_time_pts);
            temp_data=CSD(test_ch,critical_pts,test_session,test_animal);

        end
        plot(critical_time_pts*milisecond,temp_data-(kk-1)*figure_gap,'ro')

        %-- Mark the PopSpike reference time t3
        if PopSpike_reference_DGU_ch==test_ch
            plot(critical_time_pts(3)*milisecond,temp_data(3)-(kk-1)*figure_gap,'ro','MarkerFaceColor','r');
        end

        % Max slope and Max curvature 
            critical_pts_index = critical_pts(1):critical_pts(2);
            critical_data = CSD(test_ch,critical_pts_index,test_session,test_animal);
            critical_slope_index = diff(critical_data);
            [der,max_slope_candi,min_slope_candi]=slope_max_min(critical_data,sampling_freq_test);
            [der_curv,max_slope_curv_candi,min_slope_curv_candi]=slope_max_min(der,sampling_freq_test);

            if abs(max_slope_candi)>abs(min_slope_candi)
            max_slope=max_slope_candi;
            else
            max_slope=min_slope_candi;
            end
            if abs(max_slope_curv_candi)>abs(min_slope_curv_candi)
            max_slope_curv=max_slope_curv_candi;
            else
            max_slope_curv=min_slope_curv_candi;
            end

        % AUC (integral)

        Effective_time_test=time_test(critical_pts(1):critical_pts(2));
        Effective_temp_data=CSD(test_ch,critical_pts(1):critical_pts(2),test_session,test_animal);

        CSD_integral=@(t) interp1(Effective_time_test,Effective_temp_data,t,'pchip');
        auc=integral(CSD_integral,Effective_time_test(1),Effective_time_test(end));


        %--------------------------------------------------------------------------------
        % Save MANUALLY
        [slope, ~]= getSlopeDepthThreepoint(critical_time_pts(1:3,1),temp_data(1:3));
        Slope_data(test_session,Intensity_index,test_ch,test_animal)=slope(1);
        Slope_data_max(test_session,Intensity_index,test_ch,test_animal)=max_slope;
        Slope_data_curv_max(test_session,Intensity_index,test_ch,test_animal)=max_slope_curv;
        CSD_Integral_data(test_session,Intensity_index,test_ch,test_animal)=auc;

        [~, depth]= getSlopeDepthThreepoint(critical_time_pts(2:4,1),temp_data(2:4));
        Depth_data(test_session,Intensity_index,test_ch,test_animal)=depth;
        Time_data(test_session,Intensity_index,test_ch,test_animal,:)=critical_time_pts;

            %-- CA1: 1 SLM :2 DGU-Dend :3 DGU-Soma :4 DGL-Soma :5 DGL-Dend :6
        Section_Index_data(test_session,Intensity_index,1:ch_SLM_start-1,test_animal)= 1;
        Section_Index_data(test_session,Intensity_index,ch_SLM_start:ch_DGU_neg_start-1,test_animal)= 2;
        Section_Index_data(test_session,Intensity_index,ch_DGU_neg_start:ch_DGU_pos_start-1,test_animal)= 3;
        Section_Index_data(test_session,Intensity_index,ch_DGU_pos_start:ch_DGL_pos_start-1,test_animal)= 4;
        Section_Index_data(test_session,Intensity_index,ch_DGL_pos_start:ch_DGL_neg_start-1,test_animal)= 5;
        Section_Index_data(test_session,Intensity_index,ch_DGL_neg_start:end,test_animal)= 6;    
            %-- Popsulation Spike reference Channel
        PopSpike_Channel_data(test_session,Intensity_index,test_animal,1)=PopSpike_reference_CA1_ch;
        PopSpike_Channel_data(test_session,Intensity_index,test_animal,2)=PopSpike_reference_CA1_ch;
        PopSpike_Channel_data(test_session,Intensity_index,test_animal,3)=PopSpike_reference_DGU_ch;
        PopSpike_Channel_data(test_session,Intensity_index,test_animal,4)=PopSpike_reference_DGU_ch;
        PopSpike_Channel_data(test_session,Intensity_index,test_animal,5)=PopSpike_reference_DGL_ch;
        PopSpike_Channel_data(test_session,Intensity_index,test_animal,6)=PopSpike_reference_DGL_ch;
    
    end
    
  end

  disp('%%%%%%%%%%%%%%%%%%%%%%');
  
%% 3. Temporal analysis 1) save slope and depth
choice = questdlg('Are you sure you want to save?', ...
	'Saving', ...
	'Yes. I do want to save the current work','No. I do not want to save', 'No. I do not want to save');
% Handle response
saving_flag=0;
switch choice
    case 'Yes. I do want to save the current work'
        disp(' Saving start')
        saving_flag=1;
    case 'No. I do not want to save'
        disp('Terminate the process')
        error('It has been stop since the user did not want to save it');
end

if saving_flag==1
    save('CSDSlope_112118_GadR26_IO','Slope_data', 'Slope_data_max', 'Slope_data_curv_max')
    save('CSDDepth_112118_GadR26_IO','Depth_data')
    save('CSDTIME_112118_GadR26_IO','Time_data')
    save('CSDAUC_112118_GadR26_IO','CSD_Integral_data')
    save('CSDAUC_DGU_112118_GadR26_IO','CSD_AUC_DGU_data'); 
    save('CSDSecIndex_112118_GadR26_IO','Section_Index_data')
    save('CSDPopChannel_112118_GadR26_IO','PopSpike_Channel_data');
% % %  %
%     save('CSDSlope_112118_GadR28_IO_w_0ms_optic','Slope_data', 'Slope_data_max', 'Slope_data_curv_max')
%     save('CSDDepth_112118_GadR28_IO_w_0ms_optic','Depth_data')
%     save('CSDTIME_112118_GadR28_IO_w_0ms_optic','Time_data')
%     save('CSDAUC_112118_GadR28_IO_w_0ms_optic','CSD_Integral_data')
%     save('CSDAUC_DGU_112118_GadR28_IO_w_0ms_optic.mat','CSD_AUC_DGU_data'); %why there is only 1?
%     save('CSDSecIndex_112118_GadR28_IO_w_0ms_optic','Section_Index_data')
%     save('CSDPopChannel_112118_GadR28_IO_w_0ms_optic','PopSpike_Channel_data');
% %     
    
    
    
    disp('It was saved successfully.')
else
    disp('It was not saved.')
end

    
   %% 2. Temporal analysis 2) Get Area Under the Curve. (GET)
choice = questdlg('Measure the AUC?', ...
	'AUC for each channel', ...
	'Yes','No', 'No');
% Handle response

test_animal =1;

%%%%%%%%%%% Change this %%%%%%%%%%%%

AUCsession_start=1;
AUCsession_end=10;

AUCintensity_start= 150;
AUCintensity_end =250;
AUCintensity_step =50;
 
DG_start_reference =1; % 1 : select the start time from reference
                       % 0 : select the start time from t1 on each channel

% AUCchannel_list = 1:num_channel; % Compute all AUC
AUCchannel_list = [12];  % List of channel to compute the AUC

                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AUCintensity_step_size = AUCintensity_step/Intensity_step;
AUCintensity_start_index= AUCintensity_start/Intensity_step+1;
AUCintensity_end_index =AUCintensity_end/Intensity_step+1;

AUCintensity_list = AUCintensity_start_index:AUCintensity_step_size:AUCintensity_end_index;


switch choice
    case 'Yes'
        fprintf('\n----AUC Calculation----\n')
       fprintf('Stimulus Intensity:');
        fprintf('%d ',(AUCintensity_list-1)*Intensity_step);
        fprintf('\nSessions :');
        fprintf('%d ', AUCsession_start:AUCsession_end);
       
        
    case 'No'
        disp('Terminate the process')
        error('It has been stop since the user did not want to measure the AUC');
end


for AUCintensity_id=AUCintensity_list
  
    CSD=CSDToT(:,:,:,:,AUCintensity_id);

    for AUCsession_id=AUCsession_start:AUCsession_end
        
        sampling_freq_test= DataS(AUCsession_id,test_animal).sampling_freq;   

        for AUCchannel_id=AUCchannel_list
            AUCsection_id = Section_Index_data(AUCsession_id,AUCintensity_id,AUCchannel_id,test_animal);
            
            %--Get Reference PopSpike time
            PopSpike_reference_ch = PopSpike_Channel_data(AUCsession_id,AUCintensity_id,test_animal,AUCsection_id);
            PopSpike_time_start = Time_data(AUCsession_id,AUCintensity_id,PopSpike_reference_ch,test_animal,1);
            PopSpike_time_end = Time_data(AUCsession_id,AUCintensity_id,PopSpike_reference_ch,test_animal,3);
            
            %--Get critical points
            time_test=(1:size(CSD(AUCchannel_id,:,AUCsession_id),2))/sampling_freq_test/unit_scale;
            critical_pts=round(sampling_freq_test*squeeze(Time_data(AUCsession_id,AUCintensity_id,AUCchannel_id,test_animal,:)));

            if AUCsection_id<=2
                PopSpike_pts=round(sampling_freq_test*PopSpike_time_end);
                Integral_start_pts = round(sampling_freq_test*PopSpike_time_start);
                Integral_end_pts = PopSpike_pts;
            elseif AUCsection_id>=3
                PopSpike_pts=round(sampling_freq_test*PopSpike_time_end);
                if DG_start_reference==1
                    Integral_start_pts = round(sampling_freq_test*PopSpike_time_start);
                elseif DG_start_reference==0
                    Integral_start_pts = critical_pts(1);
                end
                Integral_end_pts = PopSpike_pts;
            end
            fprintf('\n Int%d| S%d| Ch%d | PopSpk ch:', (AUCintensity_id-1)*Intensity_step, AUCsession_id,AUCchannel_id);
            fprintf('%d |', PopSpike_reference_ch);
            
            

            Effective_time_test=time_test(Integral_start_pts:Integral_end_pts);
            Effective_temp_data=CSD(AUCchannel_id,Integral_start_pts:Integral_end_pts,AUCsession_id,test_animal);

            CSD_integral=@(t) interp1(Effective_time_test,Effective_temp_data,t,'pchip');
            auc=integral(CSD_integral,Effective_time_test(1),Effective_time_test(end));
        
            CSD_AUC_DGU_data(AUCsession_id,AUCintensity_id,AUCchannel_id,test_animal)=auc;
            fprintf('AUC: %d |', auc);
        end
    end
end
fprintf('\n-----------------------\n')

%3) Temporal analysis 3) save AUC
choice = questdlg('Are you sure you want to save AUC?', ...
	'Saving', ...
	'Yes. I do want to save the current AUC','No. I do not want to save', 'No. I do not want to save');
% Handle response
saving_flag=0;
switch choice
    case 'Yes. I do want to save the current AUC'
        disp(' Saving start for AUC')
        saving_flag=1;
    case 'No. I do not want to save'
        disp('Terminate the process')
        error('It has been stop since the user did not want to save it');
end

if saving_flag==1
    save('CSDAUC_DGL_100918_PPw5.mat','CSD_AUC_DGL_data'); %why there is only 1?
    
    disp('It was saved successfully.')
else
    disp('It was not saved.')
end


%% Slope Depth verification
%--Read the saved data of slope and depth.
Slope_load=load('CSDSlope_100918_PP10_training.mat');
Depth_load= load('CSDDepth_100918_PP10_training.mat');
Time_load=load('CSDTIME_100918_PP10_training.mat');
%AUC_load=load('CSDAUC_100918_PP10_training.mat');


chid=21;
intensity_list = Intensity_step*(0:(Effective_intensity_max32-1));
%--Slope plot 
figure;
cnt=0;
for i=1:Total_number_session_truncated
    cnt=cnt+1;
    subplot(1,ceil(Total_number_session_truncated),cnt)
hold on;

test_error=plot(intensity_list,Slope_load.Slope_data(i,1:Effective_intensity_max32,chid));
xlim([0, 250]);
ylim([-1*10^5, 0]);
%ylim([0, 2*10^5 ]);
axis square;

end
hold off;
figure;
cnt=0;
for i=1:Total_number_session_truncated 
    cnt=cnt+1;
    subplot(1,ceil(Total_number_session_truncated),cnt)
hold on;

test_error=plot(intensity_list,Depth_load.Depth_data(i,1:Effective_intensity_max32,chid));
% xlim([0, 300]);
% ylim([-0.1,3]);
axis square;
  
end
hold off;



%% Verify

figure;
test_session=5;
Intensity_index=6;
test_ch=1;
test_animal=3;
critical_time_pts=Time_load.Time_data(test_session,Intensity_index,test_ch,test_animal,:);
critical_pts=floor(sampling_freq_test*critical_time_pts);
temp_data=CSDToT(chid32(test_animal,test_ch,test_session),:,test_session,test_animal,Intensity_index);
plot(time_test, temp_data);
hold on;
temp_data=CSDToT(chid32(test_animal,test_ch,test_session),critical_pts,test_session,test_animal,Intensity_index);
plot(reshape(critical_time_pts*milisecond,1,4),temp_data,'ro');



CSD_Integral_data(test_session,Intensity_index,test_ch,test_animal)

