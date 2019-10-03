%%%%%%%%%%%%%%%%%%%%%%%%%%CSD_SlopeStatistics.m%%%%%%%%%%%%%%%%%
%
% Title : CSD_SlopeStatistics.m
% Detail : Originated from EEPP_SlopeStatistics.m 
%
%           
%
%   Author  : Ain Chung
%   Date    : 04/23/2017
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%CSD_SlopeStatistics.m%%%%%%%%%%%%%%%%%

%% clear and close data

clear all;

%% Load slope and Depth for all animals

Animal_name = {'PP2','PP3','PP4','PP1w','PP3w','PP5w'};  %Initial training

Total_num_animal =6;

%- [0-0] Initialization 13 intensity / 25 step size %%%%%%%%%%%

Effective_intensity_max = 6;
IntensityStep =50;
Intensity_list=IntensityStep*(0:Effective_intensity_max-1);
Intensity_list_index50=1:2:11; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%-[0] Load the data
 
raw_slope_PP2=load('CSDSlope_073019_PP2_2m.mat');  
raw_slope_PP3=load('CSDSlope_073019_PP3_2m.mat'); 
raw_slope_PP4=load('CSDSlope_073019_PP4_2m.mat');
raw_slope_PP1w=load('CSDSlope_073019_PP1w_2m.mat');
raw_slope_PP3w=load('CSDSlope_073019_PP3w_2m.mat');
raw_slope_PP5w=load('CSDSlope_073019_PP5w_2m.mat');

%- Truncate PP10 PP14 to 10 sessions as in PP17--------------------------------
Total_number_session =4; % Total number to truncate
% Num_session_16_start=1;
Num_session_start =1; % SESSION START NUMBER 9 to 18

Session_selection = Num_session_start:(Num_session_start+Total_number_session-1);

slope_PP2=raw_slope_PP2.Slope_data(Session_selection,:,:,:);
slope_PP3=raw_slope_PP3.Slope_data(Session_selection,:,:,:);
slope_PP4=raw_slope_PP4.Slope_data(Session_selection,:,:,:);
slope_PP1w=raw_slope_PP1w.Slope_data(Session_selection,:,:,:);
slope_PP3w=raw_slope_PP3w.Slope_data(Session_selection,:,:,:);
slope_PP5w=raw_slope_PP5w.Slope_data(Session_selection,:,:,:);







% size_equal =(size(slope_PP5)==size(slope_PP7)).*(size(slope_PP7)==size(slope_PP12)).*(size(slope_PP12)==size(slope_PP15));
% flag_equal=1;
% for i=1:length(size_equal)
%     flag_equal=size_equal(i)*flag_equal;
% % end
flag_equal=1;
if flag_equal
[Total_number_session, num_Stimulus_different, num_channel, stimulus_experiment_number]=size(slope_PP4);
num_Stimulus_different16=size(slope_PP2,3);   
else 
    error('Size are different in the data. Either number of session, number of stimulus, number of channel');
end

%%
%- Initialization of the data set

chid_PP2 = [2,6,8,10];%pp2:2/6/8/10 %pp3:1/3/8/10 %pp4:6/9/12/14
chid_PP3 = [1,3,8,10];
chid_PP4 = [7,9,12,14]; %PP1w:9/11/12/13 %ppw3:9/10/12/14 %ppw5:10/11/12/14
chid_PP1w = [9,11,12,14]; 
chid_PP3w = [9,10,12,14]; 
chid_PP5w = [10,11,12,14]; 

raw_chid=[chid_PP2; chid_PP3; chid_PP4; chid_PP1w; chid_PP3w; chid_PP5w];
Total_num_channel_analysis = size(raw_chid,2);

Slope_data = zeros(Total_number_session,6,32,Total_num_animal);


%%
Slope_data(:,:,1:16,1)=slope_PP2(:,1:6,:);
Slope_data(:,:,1:16,2)=slope_PP3(:,1:6,:);
Slope_data(:,:,1:16,3)=slope_PP4(:,1:6,:);
Slope_data(:,:,:,4)=slope_PP1w(:,1:2:11,:);
Slope_data(:,:,:,5)=slope_PP3w(:,1:2:11,:);
Slope_data(:,:,:,6)=slope_PP5w(:,1:2:11,:);
% Slope_data(:,:,:,6)=slope_PP15(:,1:2:11,:);



total_avg_slope=zeros(Total_number_session,num_Stimulus_different,Total_num_channel_analysis);
total_std_slope=zeros(Total_number_session,num_Stimulus_different,Total_num_channel_analysis);
total_avg_depth=zeros(Total_number_session,num_Stimulus_different,Total_num_channel_analysis);
total_std_depth=zeros(Total_number_session,num_Stimulus_different,Total_num_channel_analysis);

% Validity check


%% Initalization
chid=zeros(Total_num_animal,Total_num_channel_analysis,Total_number_session);
%Chid size : (Total_num_animal,Effective_stimulus_experiment_number,Total_number_session)
for k=1:Total_number_session

    chid(:,:,k)=raw_chid;  
end





%% [3] Inter-animal Statistics
slope_animal= zeros(Total_number_session,Effective_intensity_max,Total_num_channel_analysis,Total_num_animal);
depth_animal= zeros(Total_number_session,Effective_intensity_max,Total_num_channel_analysis,Total_num_animal);
avg_slope_animal = zeros(Total_number_session,Effective_intensity_max,Total_num_channel_analysis);
std_slope_animal = zeros(Total_number_session,Effective_intensity_max,Total_num_channel_analysis);
avg_depth_animal = zeros(Total_number_session,Effective_intensity_max,Total_num_channel_analysis);
std_depth_animal = zeros(Total_number_session,Effective_intensity_max,Total_num_channel_analysis);

% Normalized data 2018/11/12
% slope_animal_normalized=zeros(Total_number_session,Effective_intensity_max,Total_num_channel_analysis,Total_num_animal);
% avg_slope_animal_normalized = zeros(Total_number_session,Effective_intensity_max,Total_num_channel_analysis);
% std_slope_animal_normalized = zeros(Total_number_session,Effective_intensity_max,Total_num_channel_analysis);
% avg_depth_animal_normalized = zeros(Total_number_session,Effective_intensity_max,Total_num_channel_analysis);
% std_depth_animal_normalized = zeros(Total_number_session,Effective_intensity_max,Total_num_channel_analysis);

normal_session_num = 1;


% STEP 1 Plot

for i=1:Total_number_session
    for j=1:Total_num_channel_analysis
        for l= 6:-1:1

             animal_index=[];  
             
             cnt=1;

             for k=1:Total_num_animal
%                 %%%%%%%%%%%%%%%%%%%%%%%%%% Context C PP2 non-data removal%%%%%%%%%%%%%%%%
%                 if (k==2)&&context_index(i)==3  % If PP2 and Context C(Conflict is selected)
% 
%                 else
%                 %%%%%%%%%%%%%%%%%%%%%%%%%% Context C PP2 non-data removal%%%%%%%%%%%%%%%%

                animal_index=[animal_index,k]; % Exclude PP2 with Context C : usually 5 or 4 in this case. 
                slope_animal(i,l,j,k)=Slope_data(i,l,chid(k,j,i),k);
              
                 slope_animal_normalized(i,l,j,k)=slope_animal(i,l,j,k)/slope_animal(normal_session_num,6,j,k);

                
                cnt=cnt+1;
                
%                 %%%%%%%%%%%%%%%%%%%%%%%%%% Context C PP2 non-data removal%%%%%%%%%%%%%%%%
%                 end
%                 %%%%%%%%%%%%%%%%%%%%%%%%%% Context C PP2 non-data removal%%%%%%%%%%%%%%%%
             end
                
                avg_slope_animal(i,l,j)=mean(slope_animal(i,l,j,1:cnt-1));
                std_slope_animal(i,l,j)=std(slope_animal(i,l,j,1:cnt-1))/sqrt(length(animal_index));
                avg_slope_animal_normalized(i,l,j)=mean(slope_animal_normalized(i,l,j,1:cnt-1));
                std_slope_animal_normalized(i,l,j)=std(slope_animal_normalized(i,l,j,1:cnt-1))/sqrt(length(animal_index));
     
        end

    end
end

%% [3-1] Plot the inter-animal statistics

figure;
    for j=1:Total_num_channel_analysis
        
        for i=1:2:Total_number_session

        subplot(1,Total_num_channel_analysis,j)

        hold on
%         test_error=errorbar(Intensity_list,avg_slope_animal(i,:,j),std_slope_animal(i,:,j));
        test_error=errorbar(Intensity_list,avg_slope_animal_normalized(i,:,j),std_slope_animal_normalized(i,:,j));
        xlim([0, (Effective_intensity_max-1)*IntensityStep]);
        %ylim([-5 5]);
        axis square;
        hold off;    
        end
    end
    


%%
figure;
    for j=1:Total_num_channel_analysis
        for i=1:Total_number_session

        subplot(Total_num_channel_analysis,Total_number_session/2,Total_number_session/2*(j-1)+ceil(i/2))


        hold on;
        if rem(i,2)
           colr='b-'; 
        else
            colr='r-';
        end
        test_error=errorbar(Intensity_list,avg_slope_animal(i,:,j),std_slope_animal(i,:,j),colr);
        xlim([0, (Effective_intensity_max-1)*IntensityStep]);
       % ylim([-5 5]);
        axis square;
        hold off;    
        end
    end
    name_t='Sum/ w/o normalization';
    text(15,22,name_t);

%%

max_slope           = zeros(Total_number_session,Total_num_channel_analysis,Total_num_animal); % 1 case
max_slope_avg       = zeros(Total_number_session,Total_num_channel_analysis); % 1 case avg
max_slope_std       = zeros(Total_number_session,Total_num_channel_analysis); % 1 case std

max_slope_half      = zeros(Total_number_session,Total_num_channel_analysis,Total_num_animal); % 2 case
max_slope_half_avg       = zeros(Total_number_session,Total_num_channel_analysis); % 2 case avg
max_slope_half_std       = zeros(Total_number_session,Total_num_channel_analysis); % 2 case std

max_slope_variation = zeros(Total_number_session,Total_num_channel_analysis,Total_num_animal);  % 3 case : SAVE Maximum and its intensity index V(I+1)-V(I)
max_slope_variation_avg       = zeros(Total_number_session,Total_num_channel_analysis); % 3 case avg
max_slope_variation_std       = zeros(Total_number_session,Total_num_channel_analysis); % 3 case std

max_slope_variation30 = zeros(Total_number_session,Total_num_channel_analysis,Total_num_animal,2);  % 3 case : SAVE Maximum and its intensity index V(I+1)-V(I)
max_slope_variation30_avg       = zeros(Total_number_session,Total_num_channel_analysis); % 3 case avg
max_slope_variation30_std       = zeros(Total_number_session,Total_num_channel_analysis); % 3 case std

max_slope_auc           = zeros(Total_number_session,Total_num_channel_analysis,Total_num_animal); % 1 case
max_slope_auc_avg       = zeros(Total_number_session,Total_num_channel_analysis); % 1 case avg
max_slope_auc_std       = zeros(Total_number_session,Total_num_channel_analysis); % 1 case std


for i=1:Total_number_session
    for j=1:Total_num_channel_analysis
        
        animal_index=[];
        for k=1:Total_num_animal %%%%%%%%%%%%%-WARNING----Total_num_animal16 was chosen JUST TO SEE THE CONTEXT B CASE 11/09/2016 
            %%%%--Should be always Total_num_animal
            
%             %%%%%%%%%%%%%%%%%%%%%%%%%% Context C PP2 non-data removal%%%%%%%%%%%%%%%%
%             if (k==2)&&context_index(i)==3  % If PP2 and Context C(Conflict is selected)
% 
%             else
%             %%%%%%%%%%%%%%%%%%%%%%%%%% Context C PP2 non-data
%             %%%%%%%%%%%%%%%%%%%%%%%%%% removal%%%%%%%%%%%%%%%%
            
            animal_index=[animal_index,k]; % Exclude PP2 with Context C : usually 5 or 4 in this case. 
            
            
            max_slope(i,j,k)=slope_animal(i,Effective_intensity_max,j,k); % maximum slope
            slope_variation = slope_animal(i,2:Effective_intensity_max,j,k)...
                               -slope_animal(i,1:Effective_intensity_max-1,j,k);
            [~, I_max]=max(abs(slope_variation));
            max_slope_variation30(i,j,k,:)=[slope_variation(I_max); I_max];
            
            avg_func =@(x) 1*(interp1(Intensity_list,slope_animal(i,:,j,k),x,'linear')-max_slope(i,j,k)/2)^2;
            slope_50_init=Intensity_list(floor(Effective_intensity_max/2));
            slope_50=lsqnonlin(avg_func,slope_50_init,Intensity_list(1),Intensity_list(end));
            cnt=0;
            while avg_func(slope_50)>0.01&& slope_50_init<Intensity_list(Effective_intensity_max)
                cnt=cnt+1;
            slope_50_init=Intensity_list(floor(Effective_intensity_max/2)+cnt);    
                slope_50=lsqnonlin(avg_func,slope_50_init,Intensity_list(1),Intensity_list(end));
                
            end
            cnt=0;
            while avg_func(slope_50)>0.01&& slope_50_init>Intensity_list(1)
                cnt=cnt-1;
            slope_50_init=Intensity_list(floor(Effective_intensity_max/2)+cnt);    
                slope_50=lsqnonlin(avg_func,slope_50_init,Intensity_list(1),Intensity_list(end));
                
            end
            if avg_func(slope_50)>0.01
                error('Exact 50% cannot be found');
               
            end
            slope_50_int=(slope_50-mod(slope_50,IntensityStep))/IntensityStep+1;
            max_slope_half(i,j,k)=slope_50;
            max_slope_variation(i,j,k)=slope_animal(i,slope_50_int+1,j,k)-slope_animal(i,slope_50_int,j,k);
            
            avg_func_auc =@(x) interp1(Intensity_list,slope_animal(i,:,j,k),x,'linear');
            max_slope_auc(i,j,k)=integral(avg_func_auc,Intensity_list(1),Intensity_list(end));
%             %%%%%%%%%%%%%%%%%%%%%%%%%% Context C PP2 non-data removal%%%%%%%%%%%%%%%%
%             end
%             %%%%%%%%%%%%%%%%%%%%%%%%%% Context C PP2 non-data removal%%%%%%%%%%%%%%%%
        end
        
            num_animal_index=length(animal_index);
            max_slope_avg(i,j)=mean(max_slope(i,j,animal_index));
            max_slope_std(i,j)=std(max_slope(i,j,animal_index))/sqrt(num_animal_index);

            max_slope_half_avg(i,j)=mean(max_slope_half(i,j,animal_index));
            max_slope_half_std(i,j)=std(max_slope_half(i,j,animal_index))/sqrt(num_animal_index);

            max_slope_variation_avg(i,j)=mean(max_slope_variation(i,j,animal_index,1));
            max_slope_variation_std(i,j)=std(max_slope_variation(i,j,animal_index,1))/sqrt(num_animal_index);
            
            max_slope_variation30_avg(i,j)=mean(max_slope_variation30(i,j,animal_index,1));
            max_slope_variation30_std(i,j)=std(max_slope_variation30(i,j,animal_index,1))/sqrt(num_animal_index);
        
            max_slope_auc_avg(i,j)=mean(max_slope_auc(i,j,animal_index));
            max_slope_auc_std(i,j)=std(max_slope_auc(i,j,animal_index))/sqrt(num_animal_index);

        
    end
end


%% [4-1] Plot 1. Max 2. Max 50%(X) 3. Max 50%(Y)
figure;
Total_num_max_analysis=5; % case 1 2 3
Session_list= 1:2:10;
show_session_list = 1:length(Session_list);
for i=1:Total_num_max_analysis
    for j=1:Total_num_channel_analysis

        subplot(Total_num_channel_analysis,Total_num_max_analysis,Total_num_max_analysis*(j-1)+i)

        hold on;
        
        cnt=1;
        for k=Session_list
            switch i 
            case 1
                bar(cnt,max_slope_avg(k,j));
            case 2
                bar(cnt,max_slope_half_avg(k,j));
            case 3
                bar(cnt,max_slope_variation_avg(k,j));
            case 4
                bar(cnt,max_slope_variation30_avg(k,j));
            case 5
                bar(cnt,max_slope_auc_avg(k,j));
            end 
            cnt=cnt+1;
        end
        
        
        colr='ro';
        switch i 
        case 1
           % bar(Session_list,max_slope_avg(Session_list,j));
            test_error=errorbar(show_session_list,max_slope_avg(Session_list,j),max_slope_std(Session_list,j),colr);
        case 2
           % bar(Session_list,max_slope_half_avg(Session_list,j));
            test_error=errorbar(show_session_list,max_slope_half_avg(Session_list,j),max_slope_half_std(Session_list,j),colr);
        case 3
           % bar(Session_list,max_slope_variation_avg(Session_list,j));
            test_error=errorbar(show_session_list,max_slope_variation_avg(Session_list,j),max_slope_variation_std(Session_list,j),colr);
        case 4
           % bar(Session_list,max_slope_variation_avg(Session_list,j));
            test_error=errorbar(show_session_list,max_slope_variation30_avg(Session_list,j),max_slope_variation30_std(Session_list,j),colr);
        case 5
           % bar(Session_list,max_slope_variation_avg(Session_list,j));
            test_error=errorbar(show_session_list,max_slope_auc_avg(Session_list,j),max_slope_auc_std(Session_list,j),colr);
        
        end
        
        
        hold off;
        
    end
end


%%
upper_max_slope_variation30=reshape(max_slope_variation30(:,1,:,1),Total_number_session,Total_num_animal);
lower_max_slope_variation30=reshape(max_slope_variation30(:,4,:,1),Total_number_session,Total_num_animal);
upper_max_slope=reshape(max_slope(:,1,:),Total_number_session,Total_num_animal);
lower_max_slope=reshape(max_slope(:,4,:),Total_number_session,Total_num_animal);
upper_max_slope_half=reshape(max_slope_half(:,1,:),Total_number_session,Total_num_animal);
lower_max_slope_half=reshape(max_slope_half(:,4,:),Total_number_session,Total_num_animal);
upper_max_slope_variation=reshape(max_slope_variation(:,1,:,1),Total_number_session,Total_num_animal);
lower_max_slope_variation=reshape(max_slope_variation(:,4,:,1),Total_number_session,Total_num_animal);
upper_max_slope_auc=reshape(max_slope_auc(:,1,:),Total_number_session,Total_num_animal);
lower_max_slope_auc=reshape(max_slope_auc(:,4,:),Total_number_session,Total_num_animal);
