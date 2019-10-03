%get fEPSP (maximum slope of the excitatory sink of DS)
clear; close all; clc;

addpath('/Volumes/DATA/matlab/my_functions/');

drTYPE = 'DS_TYPE12/';
drART = '../MAT_ART/';
drSPEED = '../SPEED_TARGET_DIST/';
drCSDSINK = 'DS_CSD_ML_AVER/';
%drCSDSINK = 'DS_CSD_ML_AVER_200Hz/';

fnList = '../doc/fileList_v2.xlsx';
[~,~,RAW]=xlsread(fnList);
listChannels = RAW(2:end,4:12); %channels
listID = RAW(2:end,13);
listGen = RAW(2:end,14);
listExp = RAW(2:end,15);
listSinks = RAW(2:end,22); %alternative sinks


dataSetI = 1; %1=base, 2=exp, 3=2HR

if dataSetI == 1
    listFiles = RAW(2:end,1); %BASE
    listOK = RAW(2:end,23:26); %BASE
elseif dataSetI == 2
    listFiles = RAW(2:end,2); %EXP
    %listOK = RAW(2:end,27:30); %Exp
    listOK = RAW(2:end,23:26); %BASE
elseif dataSetI == 3
    listFiles = RAW(2:end,3); %2HR
    %listOK = RAW(2:end,31:34); %2HR
    listOK = RAW(2:end,23:26); %BASE
end


%BASE
%goodIDs = {{'PP5_2','PP6','PP7','PP12_2','PP13','PP15','PP18','PP19','PP20','PP23','PP25'},...
%           {'PP10','PP11','PP14','PP16','PP17','PP24','PP26','PP27'}};

%exclude PP19(noisy data), PP20 (no data)
%exclude PP11(no data ret), PP16 (no data ret), 
goodIDs = {{'PP5_2','PP6','PP7','PP12_2','PP13','PP15','PP18','PP23','PP25'},...
           {'PP10','PP11','PP14','PP16','PP17','PP24','PP26','PP27'}};

expLabels = {'PRETRAIN','TR1','TR2','TR3','RET','CON1','CON2','CON3'};
Nexp = length(expLabels);

eegFS = 2000;
winLenSec = 10*60;
winsLen = winLenSec*eegFS;
Nwin = floor(30*60 / winLenSec);

fr = 30; %frame rate
winLenMin = 1*fr; %length of window

stillRun = 1; %1=still,2=run,3=any

%avers = cell(2,2); %rows-WT,KO; columns-typeDS
%aversCSD = cell(2,4,Nexp); %rows-WT,KO; locations, Nexp
res = cell(2,2); %type sup/inf
resMouse = [];

resRaw = cell(2,2,Nexp);  %type sup/inf

for genI = 1%1:2
    
    animals = goodIDs{genI};
    
    for aI = 1:length(animals)
        
        animalID = animals{aI};
        disp(animalID)
        
        dExp = nan(2,2,Nexp); %type sup/inf exp
        
        for expI = 1:Nexp
            
            %exclude if bad file
            %k = strcmp(badFiles(:,1),animalID) & strcmp(badFiles(:,2),expLabels{expI});
            %if sum(k) > 0; continue; end;
            
            %find animal + exp combination
            k = strcmp(listID,animalID) & strcmp(listExp,expLabels{expI});
            if sum(k) == 0; continue; end;
            
            fnEEG = listFiles{k};
            if isnan(fnEEG); continue; end;
            
            channels = listChannels(k,:);
            channels = cell2mat(channels);
            chDG = channels(7);
            
            listOKexp = listOK(k,:);
            listOKexp = cell2mat(listOKexp);
            
            %if exist([drTYPE fnEEG],'file') == 0; continue; end;
            if exist([drSPEED fnEEG],'file') == 0; continue; end;
            if exist([drCSDSINK fnEEG],'file') == 0; continue; end;
            
            %load([drTYPE fnEEG]);
            %load([drART fnEEG]);
            %load([drSPEED fnEEG]);
            load([drCSDSINK fnEEG]);
            
            %{
            art = signalOK{chDG};
            
            %smooth speed (0.5s window = 15 frames)
            speedRoomTS = smth(speedRoomTS,15);
            
            if stillRun == 1
                ksp = speedRoomTS < 2;
            else
                ksp = speedRoomTS > 3;
            end
            %}
            
            %samples of type1 and type2
            %if length(samplesDS) ~= length(kType1); disp('size problem'); return; end;
            %samplesDS1 = samplesDS(kType1);
            %samplesDS2 = samplesDS(kType2);
            
            %DS1
            for typeI = 1:2
                
                for supInfI = 1:2 %sup,inf
                    
                    CSD = resultsCSDMLaver{supInfI,typeI};
                    if isempty(CSD); continue; end;
                    
                    %is file good
                    %sup1 sup2 inf1 inf2
                    ind = (supInfI-1)*2 + typeI;
                    xOK = listOKexp(ind);
                    
                    if isnan(xOK); continue; end;
                    if isempty(xOK); continue; end;
                    if xOK == 0; continue; end;
                    
                    m = CSD;
                    
                    %average
                    %m = nanmean(xCSD,1);
                    
                    resRaw{typeI,supInfI,expI} = cat(1,resRaw{typeI,supInfI,expI},m);
            
                    %find minima closest to middlepoint (sink)
                    %find preceding maxima
                    %find largest slope between them (minima in diff between two
                    %points)
          
                    kMaxes = find((m > [m(1) m(1:(end-1))]) & (m >= [m(2:end) m(end)]));
                    kMins = find((m < [m(1) m(1:(end-1))]) & (m <= [m(2:end) m(end)]));
                    N = length(m); %number of samples
                    [~,k] = min(abs(kMins-N/2));

                    kMin = kMins(k); %minima

                    kMax = kMaxes(kMaxes < kMin);
                    if isempty(kMax); continue; end;
                    kMax = kMax(end); %maxima
                    
                    d = diff(m); %differential
                    d = d(kMax:kMin);
                    d = min(d);
                    
                    %d = (m(kMax)-m(kMin)) / (kMax - kMin);
                    
                    dExp(typeI,supInfI,expI) = d; %type sup/inf exp
                end
            end
                
            %{
            %go through speed windows
            %break them to segments
            %find if they are artifact free
            %find islands of 1s
            ksp(1) = 0; ksp(end) = 0;
            kUP = find(ksp(1:end-1) == 0 & ksp(2:end) == 1); kUP = kUP + 1;
            kDOWN = find(ksp(1:end-1) == 1 & ksp(2:end) == 0);
          
            d = kDOWN - kUP; %length of overshoots
            kd = d > winLenMin; %1s in frames!
            kUP = kUP(kd);
            kDOWN = kDOWN(kd);
            
            r = [];
            
            %windows speed
            for xI = 1:length(kUP)
                
                stw = kUP(xI);
                edw = kDOWN(xI);
                
                %convert to samples
                stw = tsRoom(stw);
                edw = tsRoom(edw);
                
                %art
                if sum(art(stw:edw)==0) > 1; continue; end;
                
                T = (edw-stw+1)/eegFS;
                
                %get rates
                r1 = sum(samplesDS1 > stw & samplesDS1 < edw) / T;
                r2 = sum(samplesDS2 > stw & samplesDS2 < edw) / T;
                
                r = cat(2,r,[r1;r2]);
            end
            
            if isempty(r); continue; end;
            
            ratesExp(:,expI) = nanmean(r,2);
            %}
            

        end %exp
        
        %avers{genI} = cat(3,avers{genI},dExp);
        
        for typeI = 1:2
            for supInfI = 1:2
                x = dExp(typeI,supInfI,:);
                x = squeeze(x)';
                res{typeI,supInfI} = cat(1,res{typeI,supInfI},x); 
            end
        end
        resMouse = cat(1,resMouse,{animalID});
        
    end %animals
    
end %gen


figure('Position',[100 100 800 500]); hold on;

colors = [0.5 0.5 0.5; 1 0 0];

labels = {'SUP DS_L','SUP DS_M'; 'INF DS_L','INF DS_M'};

for typeI = 1:2
    for supInfI = 1:2
        
        subplot(2,2,(supInfI-1)*2+typeI); hold on;
        
        x = res{typeI,supInfI};
        x = x / 10000;

        m = nanmean(x,1);
        s = nanstd(x,[],1)/sqrt(size(x,1));
        
        errorbar(m,s,'Color',colors(1,:),'LineWidth',2);
        
        set(gca,'XTick',1:Nexp)
        set(gca,'XTickLabel',expLabels)

        ylabel('Max Sink');
        
        title(labels{supInfI,typeI});
        
    end
   
end


%pretrain vs RET RAW profile
figure('Position',[100 100 600 400]); hold on;

%expLabelsPlot = {'PRE','TR1','TR2','TR3','RET','CON1','CON2','CON3'};

expInds = [1,5]; %PRE-RET

for typeI = 1:2
    
    for supInfI = 1:2

        subplot(2,2,(supInfI-1)*2+typeI); hold on;
    
        for expI = 1:2
            
            expInd = expInds(expI);
        
            x = resRaw{typeI,supInfI,expInd};
            
            m = nanmean(x,1);
            s = nanstd(x,[],1)/sqrt(size(x,1));
            
            alphaVal = 0.3;
            Y = [m - s fliplr(m + s)]; %defines polygon
            X = [1:length(m) fliplr(1:length(m))]; %time from zero to end and back
            h = fill(X,Y,colors(expI,:),'HandleVisibility','off','LineStyle','none'); %fill polygon
            alpha(h, alphaVal); %set transparency

            %this removes it from legend entries
            set(get(get(h,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','off'); % Exclude line from legend

            plot (1:length(m), m,'Color',colors(expI,:),'LineWidth',2); %print mean	
        
        end
        title(labels{supInfI,typeI});
        
        axis([1 200 -15000 5000])

        set(gca,'XTick',[1,50,100,150,200])
        set(gca,'XTickLabel',{-50,-25,0,25,50})

        set(gca,'YTick',[-10000,-5000,0,5000])
        set(gca,'YTickLabel',{-1,-0.5,0,0.5})

        xlabel('Time [s]');
        ylabel('CSD a.u.');

    end
    
    
end

return;

%write stats data
labelsType = {'DS_L','DS_M'};
labelsLoc = {'Sup','Inf'};

fid = fopen('stats_maxSink_WT.dat', 'w');
fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', 'DS Type','Location','Mouse ID','Pre','Tr1','Tr2','Tr3','Ret','Con1','Con2','Con3');

for typeI = 1:2
    labType = labelsType{typeI};
    for supInfI = 1:2
        labLoc = labelsLoc{supInfI};
    
        x = res{typeI,supInfI};
        x = x/10000; %normalize

        %animals
        for i = 1:size(x,1);
            labID = resMouse{i};

            xi = x(i,:);

            fprintf(fid, '%s\t%s\t%s\t%6.6f\t%6.6f\t%6.6f\t%6.6f\t%6.6f\t%6.6f\t%6.6f\t%6.6f\n', ...
                labType,labLoc,labID,xi(1),xi(2),xi(3),xi(4),xi(5),xi(6),xi(7),xi(8));

        end
                
    end
end
fclose(fid);






