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

Animal_name = {'PPw1','PPw2','PPw3','PPw4','PPw5'};
Total_num_animal =size(Animal_name,2);

%- [0-0] Initialization 13 intensity / 25 step size %%%%%%%%%%%

Effective_intensity_max = 11;
IntensityStep =25;
Intensity_list=IntensityStep*(0:Effective_intensity_max-1);
Intensity_list_index50=1:11;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%-[0] Load the data
%------------------------
raw_slope_PPw1=load('CSDSlope_100918_PPw1_water.mat');   
raw_slope_PPw2=load('CSDSlope_100918_PPw2_water.mat');   
raw_slope_PPw3=load('CSDSlope_100918_PPw3_water.mat');   
raw_slope_PPw4=load('CSDSlope_100918_PPw4_water.mat');    
raw_slope_PPw5=load('CSDSlope_100918_PPw5_water.mat');    

% Load max slope

%- Truncate PP10 PP14 to 10 sessions as in PP17--------------------------------
Total_number_session =10; % Total number to truncate
% Num_session_16_start=1;
Num_session_start =1; % SESSION START NUMBER 9 to 18

Session_selection = Num_session_start:(Num_session_start+Total_number_session-1);

slope_PPw1=raw_slope_PPw1.Slope_data(Session_selection,:,:,:);
slope_PPw2=raw_slope_PPw2.Slope_data(Session_selection,:,:,:);
slope_PPw3=raw_slope_PPw3.Slope_data(Session_selection,:,:,:);
slope_PPw4=raw_slope_PPw4.Slope_data(Session_selection,:,:,:);
slope_PPw5=raw_slope_PPw5.Slope_data(Session_selection,:,:,:);



[Total_number_session, num_Stimulus_different, num_channel, stimulus_experiment_number]=size(slope_PPw1);


%%
%- Initialization of the data set

chid_PPw1 = [4, 8, 9, 10, 13]; 
chid_PPw2 = [4, 8, 10, 11, 16];
chid_PPw3 = [3, 7, 8, 9, 14]; 
chid_PPw4 = [3, 7, 8, 9, 12];
chid_PPw5 = [3, 7, 9, 10, 14];

raw_chid=[chid_PPw1; chid_PPw2; chid_PPw3;chid_PPw4;chid_PPw5];
Total_num_channel_analysis = size(raw_chid,2);

Slope_data = zeros(Total_number_session,num_Stimulus_different,num_channel,Total_num_animal);




Slope_data(:,:,:,1)=slope_PPw1;
Slope_data(:,:,:,2)=slope_PPw2;
Slope_data(:,:,:,3)=slope_PPw3;
Slope_data(:,:,:,4)=slope_PPw4;
Slope_data(:,:,:,5)=slope_PPw5;


total_avg_slope=zeros(Total_number_session,num_Stimulus_different,Total_num_channel_analysis);
total_std_slope=zeros(Total_number_session,num_Stimulus_different,Total_num_channel_analysis);
total_avg_depth=zeros(Total_number_session,num_Stimulus_different,Total_num_channel_analysis);
total_std_depth=zeros(Total_number_session,num_Stimulus_different,Total_num_channel_analysis);

% Validity check


% Initalization
chid=zeros(Total_num_animal,Total_num_channel_analysis,Total_number_session);

% Chid size : (Total_num_animal,Effective_stimulus_experiment_number,Total_number_session)
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
slope_animal_normalized=zeros(Total_number_session,Effective_intensity_max,Total_num_channel_analysis,Total_num_animal);
avg_slope_animal_normalized = zeros(Total_number_session,Effective_intensity_max,Total_num_channel_analysis);
std_slope_animal_normalized = zeros(Total_number_session,Effective_intensity_max,Total_num_channel_analysis);
avg_depth_animal_normalized = zeros(Total_number_session,Effective_intensity_max,Total_num_channel_analysis);
std_depth_animal_normalized = zeros(Total_number_session,Effective_intensity_max,Total_num_channel_analysis);

normal_session_num = 1;


% STEP 1 Plot

for i=1:Total_number_session
    for j=1:Total_num_channel_analysis
        for l= 11:-1:1

             animal_index=[];  
             
             cnt=1;

             for k=1:Total_num_animal


                animal_index=[animal_index,k]; % 
                slope_animal(i,l,j,k)=Slope_data(i,l,chid(k,j,i),k);
              
                 slope_animal_normalized(i,l,j,k)=slope_animal(i,l,j,k)/slope_animal(normal_session_num,11,j,k);

                
                cnt=cnt+1;
                

             end
                
                avg_slope_animal(i,l,j)=mean(slope_animal(i,l,j,1:cnt-1));
                std_slope_animal(i,l,j)=std(slope_animal(i,l,j,1:cnt-1))/sqrt(length(animal_index));
                avg_slope_animal_normalized(i,l,j)=mean(slope_animal_normalized(i,l,j,1:cnt-1));
                std_slope_animal_normalized(i,l,j)=std(slope_animal_normalized(i,l,j,1:cnt-1))/sqrt(length(animal_index));
     
        end

    end
end
% STEP 1 Plot



%% [3-1] Plot the inter-animal statistics

figure;
    for j=1:Total_num_channel_analysis
        
        for i=1:2:Total_number_session

        subplot(1,Total_num_channel_analysis,j)

        hold on
        test_error=errorbar(Intensity_list,avg_slope_animal(i,1:11,j),std_slope_animal(i,1:11,j));
        xlim([0, (Effective_intensity_max-1)*IntensityStep]);
        %ylim([-5 5]);
        axis square;
        hold off;    
        end
    end
    name_t='Sum/ w/o normalization';
    text(15,22,name_t);

figure;
    for j=1:Total_num_channel_analysis
        
        for i=1:2:Total_number_session

        subplot(1,Total_num_channel_analysis,j)

        hold on
%         test_error=errorbar(Intensity_list,avg_slope_animal(i,:,j),std_slope_animal(i,:,j));
        test_error=errorbar(Intensity_list,avg_slope_animal_normalized(i,1:11,j),std_slope_animal_normalized(i,1:11,j));
        xlim([0, (Effective_intensity_max-1)*IntensityStep]);
        %ylim([-5 5]);
        axis square;
        hold off;    
        end
    end
    name_t='Sum/ w/o normalization';
    text(15,22,name_t);

    
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
        test_error=errorbar(Intensity_list,avg_slope_animal(i,1:2:11,j),std_slope_animal(i,1:2:11,j),colr);
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
        for k=1:Total_num_animal %

            
            animal_index=[animal_index,k]; % 
            
            
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
upper_max_slope_variation30=reshape(max_slope_variation30(:,3,:,1),Total_number_session,Total_num_animal);
lower_max_slope_variation30=reshape(max_slope_variation30(:,5,:,1),Total_number_session,Total_num_animal);
upper_max_slope=reshape(max_slope(:,3,:),Total_number_session,Total_num_animal);
lower_max_slope=reshape(max_slope(:,5,:),Total_number_session,Total_num_animal);
upper_max_slope_half=reshape(max_slope_half(:,3,:),Total_number_session,Total_num_animal);
lower_max_slope_half=reshape(max_slope_half(:,5,:),Total_number_session,Total_num_animal);
upper_max_slope_variation=reshape(max_slope_variation(:,3,:,1),Total_number_session,Total_num_animal);
lower_max_slope_variation=reshape(max_slope_variation(:,5,:,1),Total_number_session,Total_num_animal);
upper_max_slope_auc=reshape(max_slope_auc(:,3,:),Total_number_session,Total_num_animal);
lower_max_slope_auc=reshape(max_slope_auc(:,5,:),Total_number_session,Total_num_animal);
