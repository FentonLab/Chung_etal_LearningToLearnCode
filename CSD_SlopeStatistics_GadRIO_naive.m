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

Animal_name = {'GadR20','Gad22','Gad23','GadR25', 'GadR26'};  %Initial training

Total_num_animal =size(Animal_name,2);
Total_num_PPI = 2;

%- [0-0] Initialization 13 intensity / 25 step size %%%%%%%%%%%

Effective_intensity_max = 11;
IntensityStep =50;
Intensity_list=IntensityStep*(0:Effective_intensity_max-1);
% Intensity_list_index50=1:2:13; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%-[0] Load the data


raw_slope_GadR20_R1=load('CSDSlope_112118_GadR20_IO.mat'); % 
raw_slope_GadR22_R1=load('CSDSlope_112118_GadR22_IO.mat'); %  
raw_slope_GadR23_R1=load('CSDSlope_112118_GadR23_IO.mat'); %  
raw_slope_GadR25_R1=load('CSDSlope_112118_GadR25_IO.mat'); %  
raw_slope_GadR26_R1=load('CSDSlope_112118_GadR26_IO.mat'); %  
raw_slope_GadR20_R2=load('CSDSlope_112118_GadR20_IO_w_optic.mat'); % 
raw_slope_GadR22_R2=load('CSDSlope_112118_GadR22_IO_w_optic.mat'); %  
raw_slope_GadR23_R2=load('CSDSlope_112118_GadR23_IO_w_optic.mat'); %  
raw_slope_GadR25_R2=load('CSDSlope_112118_GadR25_IO_w_optic.mat'); %
raw_slope_GadR26_R2=load('CSDSlope_112118_GadR26_IO_w_optic.mat'); %

raw_depth_GadR20_R1=load('CSDDepth_112118_GadR20_IO.mat'); % 
raw_depth_GadR22_R1=load('CSDDepth_112118_GadR22_IO.mat'); %  
raw_depth_GadR23_R1=load('CSDDepth_112118_GadR23_IO.mat'); %  
raw_depth_GadR25_R1=load('CSDDepth_112118_GadR25_IO.mat'); % 
raw_depth_GadR26_R1=load('CSDDepth_112118_GadR26_IO.mat'); % 
raw_depth_GadR20_R2=load('CSDDepth_112118_GadR20_IO_w_optic.mat'); % 
raw_depth_GadR22_R2=load('CSDDepth_112118_GadR22_IO_w_optic.mat'); %  
raw_depth_GadR23_R2=load('CSDDepth_112118_GadR23_IO_w_optic.mat'); %  
raw_depth_GadR25_R2=load('CSDDepth_112118_GadR25_IO_w_optic.mat'); %  
raw_depth_GadR26_R2=load('CSDDepth_112118_GadR26_IO_w_optic.mat'); %  


Total_number_session =1; % Total number to truncate

Num_session_start =1; 

Session_selection = Num_session_start:(Num_session_start+Total_number_session-1);


slope_GadR20_R1=raw_slope_GadR20_R1.Slope_data(Session_selection,:,:);
slope_GadR20_R2=raw_slope_GadR20_R2.Slope_data(Session_selection,:,:);
slope_GadR22_R1=raw_slope_GadR22_R1.Slope_data(Session_selection,:,:);
slope_GadR22_R2=raw_slope_GadR22_R2.Slope_data(Session_selection,:,:);
slope_GadR23_R1=raw_slope_GadR23_R1.Slope_data(Session_selection,:,:);
slope_GadR23_R2=raw_slope_GadR23_R2.Slope_data(Session_selection,:,:);
slope_GadR25_R1=raw_slope_GadR25_R1.Slope_data(Session_selection,:,:);
slope_GadR25_R2=raw_slope_GadR25_R2.Slope_data(Session_selection,:,:);
slope_GadR26_R1=raw_slope_GadR26_R1.Slope_data(Session_selection,:,:);
slope_GadR26_R2=raw_slope_GadR26_R2.Slope_data(Session_selection,:,:);

depth_GadR20_R1=raw_depth_GadR20_R1.Depth_data(Session_selection,:,:);
depth_GadR20_R2=raw_depth_GadR20_R2.Depth_data(Session_selection,:,:);
depth_GadR22_R1=raw_depth_GadR22_R1.Depth_data(Session_selection,:,:);
depth_GadR22_R2=raw_depth_GadR22_R2.Depth_data(Session_selection,:,:);
depth_GadR23_R1=raw_depth_GadR23_R1.Depth_data(Session_selection,:,:);
depth_GadR23_R2=raw_depth_GadR23_R2.Depth_data(Session_selection,:,:);
depth_GadR25_R1=raw_depth_GadR25_R1.Depth_data(Session_selection,:,:);
depth_GadR25_R2=raw_depth_GadR25_R2.Depth_data(Session_selection,:,:);
depth_GadR26_R1=raw_depth_GadR26_R1.Depth_data(Session_selection,:,:);
depth_GadR26_R2=raw_depth_GadR26_R2.Depth_data(Session_selection,:,:);


[Total_number_session, num_Stimulus_different, num_channel]=size(slope_GadR20_R1);
num_Stimulus_different16=size(slope_GadR20_R1,2);   



%%
%- Initialization of the data set


chid_GadR20 = [7,9,12,14];% 
chid_GadR22 = [7,9,13,15]; % R1 correct to ch13 and R2 measure ch7 again 40
chid_GadR23 = [8,11,13,16]; % 
chid_GadR25 = [7,9,12,14];
chid_GadR26 = [5,7,11,13];


raw_chid=[chid_GadR20; chid_GadR22; chid_GadR23; chid_GadR25; chid_GadR26];
Total_num_channel_analysis = size(raw_chid,2);

Slope_data = zeros(Total_number_session,num_Stimulus_different,num_channel,Total_num_animal,Total_num_PPI);


%%
Slope_data(:,:,:,1,1)=slope_GadR20_R1(:,:,:);
Slope_data(:,:,:,1,2)=slope_GadR20_R2(:,:,:);
Slope_data(:,:,:,2,1)=slope_GadR22_R1(:,:,:);
Slope_data(:,:,:,2,2)=slope_GadR22_R2(:,:,:);
Slope_data(:,:,:,3,1)=slope_GadR23_R1(:,:,:);
Slope_data(:,:,:,3,2)=slope_GadR23_R2(:,:,:);
Slope_data(:,:,:,4,1)=slope_GadR25_R1(:,:,:);
Slope_data(:,:,:,4,2)=slope_GadR25_R2(:,:,:);
Slope_data(:,:,:,5,1)=slope_GadR26_R1(:,:,:);
Slope_data(:,:,:,5,2)=slope_GadR26_R2(:,:,:);

Depth_data(:,:,:,1,1)=depth_GadR20_R1(:,:,:);
Depth_data(:,:,:,1,2)=depth_GadR20_R2(:,:,:);
Depth_data(:,:,:,2,1)=depth_GadR22_R1(:,:,:);
Depth_data(:,:,:,2,2)=depth_GadR22_R2(:,:,:);
Depth_data(:,:,:,3,1)=depth_GadR23_R1(:,:,:);
Depth_data(:,:,:,3,2)=depth_GadR23_R2(:,:,:);
Depth_data(:,:,:,4,1)=depth_GadR25_R1(:,:,:);
Depth_data(:,:,:,4,2)=depth_GadR25_R2(:,:,:);
Depth_data(:,:,:,5,1)=depth_GadR26_R1(:,:,:);
Depth_data(:,:,:,5,2)=depth_GadR26_R2(:,:,:);


% Initalization
chid=zeros(Total_num_animal,Total_num_channel_analysis,Total_number_session);

%% Chid size : (Total_num_animal,Effective_stimulus_experiment_number,Total_number_session)
for k=1:Total_number_session

    chid(:,:,k)=raw_chid;  
end





%% [3] Inter-animal Statistics
slope_ratio_animal= zeros(Total_number_session,Effective_intensity_max,Total_num_channel_analysis,Total_num_animal);
depth_ratio_animal= zeros(Total_number_session,Effective_intensity_max,Total_num_channel_analysis,Total_num_animal);
avg_slope_ratio_animal = zeros(Total_number_session,Effective_intensity_max,Total_num_channel_analysis);
std_slope_ratio_animal = zeros(Total_number_session,Effective_intensity_max,Total_num_channel_analysis);
avg_depth_ratio_animal = zeros(Total_number_session,Effective_intensity_max,Total_num_channel_analysis);
std_depth_ratio_animal = zeros(Total_number_session,Effective_intensity_max,Total_num_channel_analysis);


slope_animal= zeros(Total_number_session,Effective_intensity_max,Total_num_channel_analysis,Total_num_animal);
depth_animal= zeros(Total_number_session,Effective_intensity_max,Total_num_channel_analysis,Total_num_animal);
avg_slope_animal = zeros(Total_number_session,Effective_intensity_max,Total_num_channel_analysis);
std_slope_animal = zeros(Total_number_session,Effective_intensity_max,Total_num_channel_analysis);
avg_depth_animal = zeros(Total_number_session,Effective_intensity_max,Total_num_channel_analysis);
std_depth_animal = zeros(Total_number_session,Effective_intensity_max,Total_num_channel_analysis);

slope_animal_optic= zeros(Total_number_session,Effective_intensity_max,Total_num_channel_analysis,Total_num_animal);
depth_animal_optic= zeros(Total_number_session,Effective_intensity_max,Total_num_channel_analysis,Total_num_animal);
avg_slope_animal_optic = zeros(Total_number_session,Effective_intensity_max,Total_num_channel_analysis);
std_slope_animal_optic = zeros(Total_number_session,Effective_intensity_max,Total_num_channel_analysis);
avg_depth_animal_optic = zeros(Total_number_session,Effective_intensity_max,Total_num_channel_analysis);
std_depth_animal_optic = zeros(Total_number_session,Effective_intensity_max,Total_num_channel_analysis);



% Normalized data 2018/11/12
slope_animal_normalized=zeros(Total_number_session,Effective_intensity_max,Total_num_channel_analysis,Total_num_animal);
avg_slope_animal_normalized = zeros(Total_number_session,Effective_intensity_max,Total_num_channel_analysis);
std_slope_animal_normalized = zeros(Total_number_session,Effective_intensity_max,Total_num_channel_analysis);
avg_depth_animal_normalized = zeros(Total_number_session,Effective_intensity_max,Total_num_channel_analysis);
std_depth_animal_normalized = zeros(Total_number_session,Effective_intensity_max,Total_num_channel_analysis);

slope_animal_normalized_optic=zeros(Total_number_session,Effective_intensity_max,Total_num_channel_analysis,Total_num_animal);
avg_slope_animal_normalized_optic = zeros(Total_number_session,Effective_intensity_max,Total_num_channel_analysis);
std_slope_animal_normalized_optic = zeros(Total_number_session,Effective_intensity_max,Total_num_channel_analysis);
avg_depth_animal_normalized_optic = zeros(Total_number_session,Effective_intensity_max,Total_num_channel_analysis);
std_depth_animal_normalized_optic = zeros(Total_number_session,Effective_intensity_max,Total_num_channel_analysis);

normal_session_num = 1;


% STEP 1 Plot

for i=1:Total_number_session
    for j=1:Total_num_channel_analysis
        for l= 11:-1:1

             animal_index=[];  
             
             cnt=1;

             for k=1:Total_num_animal
%                 %%%%%%%%%%%%%%%%%%%%%%%%%% Context C PP2 non-data removal%%%%%%%%%%%%%%%%
%                 if (k==2)&&context_index(i)==3  % If PP2 and Context C(Conflict is selected)
% 
%                 else
%                 %%%%%%%%%%%%%%%%%%%%%%%%%% Context C PP2 non-data removal%%%%%%%%%%%%%%%%

                animal_index=[animal_index,k]; % Exclude PP2 with Context C : usually 5 or 4 in this case. 
                slope_animal(i,l,j,k)=Slope_data(i,l,chid(k,j,i),k,1);
                 depth_animal(i,l,j,k)=Depth_data(i,l,chid(k,j,i),k,1);
                 slope_animal_normalized(i,l,j,k)=slope_animal(i,l,j,k)/slope_animal(normal_session_num,11,j,k);
                 depth_animal_normalized(i,l,j,k)=depth_animal(i,l,j,k)/depth_animal(normal_session_num,11,j,k);
                
                 slope_animal_optic(i,l,j,k)=Slope_data(i,l,chid(k,j,i),k,2);
                 depth_animal_optic(i,l,j,k)=Depth_data(i,l,chid(k,j,i),k,2);
                 slope_animal_normalized_optic(i,l,j,k)=slope_animal_optic(i,l,j,k)/slope_animal(normal_session_num,11,j,k);
                 depth_animal_normalized_optic(i,l,j,k)=depth_animal_optic(i,l,j,k)/depth_animal(normal_session_num,11,j,k);
                             
                cnt=cnt+1;
                
%                 %%%%%%%%%%%%%%%%%%%%%%%%%% Context C PP2 non-data removal%%%%%%%%%%%%%%%%
%                 end
%                 %%%%%%%%%%%%%%%%%%%%%%%%%% Context C PP2 non-data removal%%%%%%%%%%%%%%%%
             end
                
                avg_slope_animal(i,l,j)=mean(slope_animal(i,l,j,1:cnt-1));
                std_slope_animal(i,l,j)=std(slope_animal(i,l,j,1:cnt-1))/sqrt(length(animal_index));
                avg_slope_animal_normalized(i,l,j)=mean(slope_animal_normalized(i,l,j,1:cnt-1));
                std_slope_animal_normalized(i,l,j)=std(slope_animal_normalized(i,l,j,1:cnt-1))/sqrt(length(animal_index));
    
                avg_slope_animal_optic(i,l,j)=mean(slope_animal_optic(i,l,j,1:cnt-1));
                std_slope_animal_optic(i,l,j)=std(slope_animal_optic(i,l,j,1:cnt-1))/sqrt(length(animal_index));
                avg_slope_animal_normalized_optic(i,l,j)=mean(slope_animal_normalized_optic(i,l,j,1:cnt-1));
                std_slope_animal_normalized_optic(i,l,j)=std(slope_animal_normalized_optic(i,l,j,1:cnt-1))/sqrt(length(animal_index));
    
                
                avg_depth_animal(i,l,j)=mean(depth_animal(i,l,j,1:cnt-1));
                std_depth_animal(i,l,j)=std(depth_animal(i,l,j,1:cnt-1))/sqrt(length(animal_index));
                avg_depth_animal_optic(i,l,j)=mean(depth_animal_optic(i,l,j,1:cnt-1));
                std_depth_animal_optic(i,l,j)=std(depth_animal_optic(i,l,j,1:cnt-1))/sqrt(length(animal_index));
  
                avg_depth_animal_normalized(i,l,j)=mean(depth_animal_normalized(i,l,j,1:cnt-1));
                std_depth_animal_normalized(i,l,j)=std(depth_animal_normalized(i,l,j,1:cnt-1))/sqrt(length(animal_index));
     
        end

    end
end





for i=1:Total_number_session
    for j=1:Total_num_channel_analysis
        
        for l=1:num_Stimulus_different16 %% for correcting normalization errors
            for k=1:Total_num_animal    
                
                slope_ratio_animal(i,l,j,k)=Slope_data(i,l,chid(k,j),k,2)/Slope_data(i,l,chid(k,j),k,1);
                depth_ratio_animal(i,l,j,k)=Depth_data(i,l,chid(k,j),k,2)/Depth_data(i,l,chid(k,j),k,1);
                
                if isnan(slope_ratio_animal(i,l,j,k))
                   slope_ratio_animal(i,l,j,k) = 0;
                elseif isnan(depth_ratio_animal(i,l,j,k))
                    
                end
            end
            
                avg_slope_ratio_animal(i,l,j)=mean(slope_ratio_animal(i,l,j,:));
                std_slope_ratio_animal(i,l,j)=std(slope_ratio_animal(i,l,j,:))/sqrt(Total_num_animal);
                avg_depth_ratio_animal(i,l,j)=mean(depth_ratio_animal(i,l,j,:));
                std_depth_ratio_animal(i,l,j)=std(depth_ratio_animal(i,l,j,:))/sqrt(Total_num_animal);
                
        end

    end
end






%% [3-1] Plot the inter-animal statistics

figure;
    for j=1:Total_num_channel_analysis
        
        for i=1:1:Total_number_session

        subplot(1,Total_num_channel_analysis,j)

        hold on
        test_error=errorbar(Intensity_list,avg_slope_animal(i,:,j),std_slope_animal(i,:,j));

        
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
        
        for i=1:1:Total_number_session

        subplot(1,Total_num_channel_analysis,j)

        hold on
        test_error=errorbar(Intensity_list,avg_slope_animal_optic(i,:,j),std_slope_animal_optic(i,:,j));

        
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
        
        for i=1:1:Total_number_session

        subplot(1,Total_num_channel_analysis,j)

        hold on
        test_error=errorbar(Intensity_list,avg_slope_animal_normalized(i,:,j),std_slope_animal_normalized(i,:,j));

        
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
        
        for i=1:1:Total_number_session

        subplot(1,Total_num_channel_analysis,j)

        hold on
        test_error=errorbar(Intensity_list,avg_slope_animal_normalized_optic(i,:,j),std_slope_animal_normalized_optic(i,:,j));

        
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
        
        for i=1:1:Total_number_session

        subplot(1,Total_num_channel_analysis,j)

        hold on
        test_error=errorbar(Intensity_list,avg_depth_animal(i,:,j),std_depth_animal(i,:,j));

        
        xlim([0, (Effective_intensity_max-1)*IntensityStep]);
        %ylim([-5 5]);
        axis square;
        hold off;    
        end
    end
    name_t='depth';
    text(15,22,name_t);
figure;
    for j=1:Total_num_channel_analysis
        
        for i=1:1:Total_number_session

        subplot(1,Total_num_channel_analysis,j)

        hold on
        test_error=errorbar(Intensity_list,avg_depth_animal_optic(i,:,j),std_depth_animal_optic(i,:,j));

        
        xlim([0, (Effective_intensity_max-1)*IntensityStep]);
        %ylim([-5 5]);
        axis square;
        hold off;    
        end
    end
    
%% [4] 1. Maximum slope 2. 50% of the maximum 3. Maximum slope variation.


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

max_slope_optic_auc           = zeros(Total_number_session,Total_num_channel_analysis,Total_num_animal); % 1 case
max_slope_optic_auc_avg       = zeros(Total_number_session,Total_num_channel_analysis); % 1 case avg
max_slope_optic_auc_std       = zeros(Total_number_session,Total_num_channel_analysis); % 1 case std


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
            max_slope_variation(i,j,k)=slope_ratio_animal(i,slope_50_int+1,j,k)-slope_ratio_animal(i,slope_50_int,j,k);
            
            avg_func_auc =@(x) interp1(Intensity_list,slope_animal(i,:,j,k),x,'linear');
            max_slope_auc(i,j,k)=integral(avg_func_auc,Intensity_list(1),Intensity_list(end));
            avg_func_auc_optic =@(x) interp1(Intensity_list,slope_animal_optic(i,:,j,k),x,'linear');
            max_slope_optic_auc(i,j,k)=integral(avg_func_auc_optic,Intensity_list(1),Intensity_list(end));
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

            max_slope_optic_auc_avg(i,j)=mean(max_slope_optic_auc(i,j,animal_index));
            max_slope_optic_auc_std(i,j)=std(max_slope_optic_auc(i,j,animal_index))/sqrt(num_animal_index);
        
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
upper_max_slope_optic_auc=reshape(max_slope_optic_auc(:,1,:),Total_number_session,Total_num_animal);
lower_max_slope_optic_auc=reshape(max_slope_optic_auc(:,4,:),Total_number_session,Total_num_animal);

