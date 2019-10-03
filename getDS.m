%get CSD of DS
%clear; close all; clc;
%indeces = 1..21
function getDS(indeces)


runLocal = 0;

if runLocal == 0 %linux
    addpath('/home/dd1348/matlab/my_functions/');
    addpath('/home/dd1348/matlab/OTC/');
    addpath('/home/dd1348/matlab/CSDplotter');
    addpath('/home/dd1348/matlab/CSDplotter/methods/');
    %addpath('/mnt/kant/dinod/matlab/my_functions/');
    %addpath('/mnt/kant/dinod/matlab/OTC/');
    %addpath('/mnt/kant/dinod/matlab/CSDplotter');
    %addpath('/mnt/kant/dinod/matlab/CSDplotter/methods/');
else %local
    addpath('/Volumes/DATA/matlab/my_functions/');
    addpath('/Volumes/DATA/matlab/OTC/');
    addpath('/Volumes/DATA/matlab/CSDplotter');
    addpath('/Volumes/DATA/matlab/CSDplotter/methods/');
end

%CSD
h = 50e-6; %electrode distance m
ex_cond = 0.3; %external conductivity S/m
top_cond = 0.3; %S/m
diam = 0.5e-3; %mm
gauss_sigma = 0.05e-3; %mm
filter_range = 5*gauss_sigma; % numeric filter must be finite in extent
%dt = 0.5; %sampling time ms (for 2k);
dt = 0.125; %sampling time in ms (for 8k);
scale_plot = 1; %focus
max_plot = 0;

CSDres = 200;
methodSpline = 1;

labelsFeat = {'sample','promZ','peakZ','peakRaw','promZbefore','promZafter',... %1..6
        'peakZbefore','peakZafter','minZbefore','minZafter',...%7..10
        'wDG','wBA','w50','tB','tA','tminB','tminA','acc','speed','max10_30before','max10_30after'}; %11..21

drMAT = '/scratch/dd1348/M32/MAT_EEG/';
drART = '/scratch/dd1348/M32/MAT_ART/';
drACC = '/scratch/dd1348/M32/MAT_ACC/';
drSP = '/scratch/dd1348/M32/SPEED_TARGET_DIST/';

drOUT = '/scratch/dd1348/M32/DS_detected_nolimit/';

fnList = 'fileList_v2.xlsx';
[~,~,RAW]=xlsread(fnList);
listBase = RAW(2:end,1); %base
listExp = RAW(2:end,2); %exp
list2HR = RAW(2:end,3); %2hr
listChannels = RAW(2:end,4:12); %channels
listID = RAW(2:end,13);
listGen = RAW(2:end,14);
listSinks = RAW(2:end,19); %alternative sinks
listDG = cell2mat(listChannels(:,7));

IDs = unique(listID);

%create filter for detection
eegFS = 2000;
[b, gd] = getFIRbandpass( 5,100,4,40,eegFS );

winLen = 100;

chAround = 5; %+/- channels around

winBA = [60 20]; %search window for maxima before/after in samples (-60..-40)

thMin = [];

%animals
for i = indeces
    ID = IDs{i};
    
    k = strcmp(listID,ID);
    filesBase = listBase(k);
    filesExp = listExp(k);
    files2HR = list2HR(k);
    channelsDG = listDG(k);
    
    %load all files to normalize them in the same way
    eeg = [];
    art = [];
    for fI = 1:length(filesBase)
        fnEEG = filesBase{fI};
        chDG = channelsDG(fI);
        
        if exist([drMAT fnEEG],'file') == 0; continue; end;
        
        load([drMAT fnEEG]);
        load([drART fnEEG]);
        
        eeg = cat(2,eeg,eegData(chDG,:));
        art = cat(2,art,signalOK{chDG});
    end
    
    art = logical(art);
    
    %get normalization values
    mx = nanmean(eeg(art));
    sx = nanstd(eeg(art));
        
    for typeI = 1:3
        if typeI == 1
            files = filesBase;
        elseif typeI == 2
            files = filesExp;
        elseif typeI == 3
            files = files2HR;
        end
            
        %process files independently
        for fI = 1:length(files)
            fnEEG = files{fI};
            chDG = channelsDG(fI);

            fnACC = [fnEEG(1:end-7) 'acc.mat'];
            if exist([drMAT fnEEG],'file') == 0; disp([fnEEG ' EEG does not exist']); continue; end;
            if exist([drACC fnACC],'file') == 0; disp([fnEEG ' ACC does not exist']); continue; end;
    
            if exist([drSP fnEEG],'file') == 0; 
                %disp([fnEEG ' Speed does not exist']); 
                speedExist = 0;
            else
                speedExist = 1;
            end
    
            if exist([drOUT fnEEG],'file') ~= 0; continue; end; %skip existing
    
            disp(fnEEG);

            load([drMAT fnEEG]);
            load([drACC fnACC]);
            load([drART fnEEG]);

            if speedExist == 1
                load([drSP fnEEG]);
            end

            eegData = double(eegData);

            Nchannels = size(eegData,1);

            %accelaration
            accX = (accData(1,:) - 1.6)/0.32;
            accY = (accData(2,:) - 1.6)/0.32;
            accZ = (accData(3,:) - 1.6)/0.32;
            acc = (accX.^2 + accY.^2 + accZ.^2).^(0.5);
            acc = acc-nanmean(acc); %should be close to 1

            %DG electrode
            eegDG = eegData(chDG,:);

            %2-100Hz
            chST = chDG - chAround;
            chED = chDG + chAround;
            if chED > Nchannels; chED = Nchannels; end;
            if chST < 1; chST = 1; end;
            channelsProfiles = chST:chED;
            chDGk = find(channelsProfiles == chDG);
            eegF = filter(b,1,eegData(chST:chED,:),[],2);
            eegF = cat(2,eegF(:,gd+1:end),zeros(size(eegF,1),gd));
            %DG channel
            eegFDG = eegF(chDGk,:);

            art = signalOK{chDG};

            %{
            %z-score normalize based on DG channel
            m = nanmean(eegFDG(art));
            s = nanstd(eegFDG(art));
            eegFDG = (eegFDG-m)/s;
            %}

            %global normalization values
            eegFDG = (eegFDG-mx)/sx;

            %DG - maxima and minima (extract width,prominence)
            [MaxesDG,kMaxesDG,wDG,pDG] = findpeaks(eegFDG);
            [~,kMinsDG] = findpeaks(1.01*max(eegFDG) - eegFDG);
            %MinsDG = eegFDG(kMinsDG);

            if speedExist == 1
                %check if tsRoom is increasing!
                dx = diff(tsRoom);
                if sum(dx == 0) > 0
                    disp('ZERO-DIFF TSROOM!');
                    speedExist = 0; %behave like there is no speed
                end
                if sum(dx < 0) > 0
                    disp('NEG-DIFF TSROOM!'); 
                    speedExist = 0; %behave like there is no speed
                end

                if speedExist == 1
                    %find Speed for all samples
                    edges = mean([tsRoom(2:end); tsRoom(1:end-1)]);
                    %d = edges(2)-edges(1);
                    d = 66; %33ms HARDCODED!
                    edges = cat(2,edges(1)-d, edges, edges(end)+d);
                    ind = discretize(kMaxesDG,edges);
                end
            end

            %for each maxima in DG
            feat = nan(length(kMaxesDG),23);

            for i = 10:length(kMaxesDG)-10 %skip first couple of minima for before/after min-max computations
                s = kMaxesDG(i); %sample
                pDGi = pDG(i); %prominence of peak in z-score signal
                wDGi = wDG(i)/eegFS*1000; %width of peak in z-score signal
                pDGB = pDG(i-1); %prominence of peak before in z-score signal
                pDGA = pDG(i+1); %prominence of peak after in z-score signal
                tDGB = (s - kMaxesDG(i-1))/eegFS*1000; %mseconds to maxima before
                tDGA = (kMaxesDG(i+1) - s)/eegFS*1000; %mseconds to maxima after
                mDG = MaxesDG(i); %max value in z-scored filtered signal
                acci = acc(s); %acceleration
                mDGB = MaxesDG(i-1); %preceding maxima in z-scored filtered signal
                mDGA = MaxesDG(i+1); %following maxima in z-scored filtered signal

                %indexes preceeding and following minima
                kB = kMinsDG(kMinsDG < s);
                kA = kMinsDG(kMinsDG > s);
                if isempty(kB) || isempty(kA); continue; end;
                kB = kB(end);
                kA = kA(1);
                %time to minima before
                tmB = (s - kB)/eegFS*1000; %mseconds to minima before
                tmA = (kA - s)/eegFS*1000; %mseconds to minima after

                %maxima in raw signal
                mr = max(eegDG(s-10:s+10));

                %peak to to minima before and after
                mB = eegFDG(kB);
                mA = eegFDG(kA);

                %speed
                spi = nan;
                if speedExist == 1 && ~isnan(ind(i))
                    spi = speedRoomTS(ind(i));
                end

                %width between minima
                wBA = (kA-kB)/eegFS*1000;

                %width 50%
                w50 = nan;
                k50B = nan;
                k50A = nan;
                if mB > mA %before higher
                    level = (mDG + mB)/2;
                else %after higher
                    level = (mDG + mA)/2;
                end
                x = eegFDG(kB:kA);
                %cut before maxima
                k = find(x(1:end-1) <= level & x(2:end) > level);
                if ~isempty(k)
                    k50B = k(1);
                end
                %cut after maxima
                k = find(x(1:end-1) >= level & x(2:end) < level);
                if ~isempty(k)
                    k50A = k(end);
                end
                if ~isnan(k50B) && ~isnan(k50A)
                    w50 = (k50A-k50B)/eegFS*1000;
                end

                %winBA = [60 20];
                %maxima before
                x = eegFDG(s-winBA(1):s-winBA(2));
                maxB = max(x);
                %maxima after
                x = eegFDG(s+winBA(2):s+winBA(1));
                maxA = max(x);

                %minima in raw before
                mrB = min(eegDG(kB-10:kB+10));

                %minima in raw after
                mrA = min(eegDG(kA-10:kA+10));

                feat(i,1) = s; %sample max sDG
                feat(i,2) = pDGi; %prominence of peak in z-score signal
                feat(i,3) = mDG;  %max value in z-scored filtered signal
                feat(i,4) = mr; %maxima in raw
                feat(i,5) = pDGB; %prominence of peak before in z-score signal
                feat(i,6) = pDGA; %prominence of peak after in z-score signal
                feat(i,7) = mDGB; %preceding maxima in z-scored filtered signal
                feat(i,8) = mDGA; %following maxima in z-scored filtered signal

                feat(i,9) = mB; %minima preceding peak in z-score signal
                feat(i,10) = mA; %minima following peak in z-score signal

                feat(i,11) = wDGi; %width of peak in z-score signal
                feat(i,12) = wBA; %width between minima
                feat(i,13) = w50; %width 50%

                feat(i,14) = tDGB; %mseconds to maxima before
                feat(i,15) = tDGA; %mseconds to maxima after
                feat(i,16) = tmB; %mseconds to minima before
                feat(i,17) = tmA; %mseconds to minima after

                feat(i,18) = acci; %acceleration
                feat(i,19) = spi; %speed

                feat(i,20) = maxB; %maxima 30-10ms before
                feat(i,21) = maxA; %maxima 10-30ms after

                feat(i,22) = mrB; %minima in raw before
                feat(i,23) = mrA; %minima in raw after
            end

            %remove nans
            knan = ~isnan(feat(:,1));
            feat = feat(knan,:);

            %apply thMin
            if ~isempty(thMin)
                k = feat(:,3) > thMin;
                %d = feat(:,3) - feat(:,9); %dist peak to minima before
                %k = d > 2;
                feat = feat(k,:);
            end

            %extract profiles around maxima
            s = feat(:,1); %samples;
            profilesDS = eegF(:,s);

            save([drOUT fnEEG],'feat','chDG','profilesDS','channelsProfiles','labelsFeat','mx','sx');
        end
    end
end

    


