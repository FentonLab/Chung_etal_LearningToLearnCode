%classify DS as significant from no_limit data
%classify DS into type1 and2
%do independently for each file
%clear; close all; clc;
%function getDStypesSink(fST,fED,dataSetI)
%fST = 1;
%fED = 50;
%dataSetI = 1;

drCSD = 'DS_CSD/';
drDS = 'DS_detected_nolimit/';
drOUT = 'DS_TYPE12/';
drART = '../MAT_ART/';
fnList = 'fileList_v2.xlsx';

[~,~,RAW]=xlsread(fnList);
listChannels = RAW(2:end,4:12); %channels
listID = RAW(2:end,13);
listGen = RAW(2:end,14);
listSinks = RAW(2:end,22); %alternative sinks

listBase = RAW(2:end,1); %file name
listExp = RAW(2:end,2); %file name
list2HR = RAW(2:end,3); %file name
     
for dataSetI = 1:3

    if dataSetI == 1
        files = listBase;
    elseif dataSetI == 2
        files = listExp;
    elseif dataSetI == 3
        files = list2HR;
    end
    
    for fI = 1:length(files)

        fnEEG = files{fI};
        disp(fnEEG);

        channels = listChannels(fI,:);
        channels = cell2mat(channels);
        chDG = channels(7);
        
        if isnan(fnEEG); continue; end;

        if exist([drOUT fnEEG],'file') ~= 0; continue; end;

        %if exist([drMAT fnEEG],'file') == 0; continue; end;
        if exist([drDS fnEEG],'file') == 0; continue; end;
        if exist([drCSD fnEEG],'file') == 0; continue; end;
        if  exist([drART fnEEG],'file') == 0; continue; end;

        wLenSearch = 4;
        if contains(fnEEG,'PP19') %double the minima search window for PP19
            wLenSearch = 8;
        end

        %get sinks
        sinks = listSinks{fI};
        if isnan(sinks); continue; end;
        if isempty(sinks); continue; end;
        s = strsplit(sinks,',');
        sinks = nan(1,length(s));
        for i = 1:length(s)
            sinks(i) = str2double(s{i});
        end

        load([drCSD fnEEG]);
        load([drDS fnEEG]);
        load([drART fnEEG]);

        %art
        art = signalOK{chDG};

        %labelsFeat = {'sample','promZ','peakZ','peakRaw','promZbefore','promZafter',... %1..6
        %    'peakZbefore','peakZafter','minZbefore','minZafter',...%7..10
        %    'wDG','wBA','w50','tB','tA','tminB','tminA','acc','speed','max10_30before','max10_30after'}; %11..21

        s = feat(:,1); %sample

        %exclude art
        kart = art(s)';
        feat = feat(kart,:);

        s = feat(:,1); %sample
        pwrZ = feat(:,3); %peakZ
        peakZbefore = feat(:,7); %peakZbefore
        peakZafter = feat(:,8); %peakZafter
        minZbefore = feat(:,9); %minZbefore
        minZafter = feat(:,10); %minZafter
        wDG = feat(:,11); %width

        diffFromMinBefore = pwrZ - minZbefore;
        diffFromMinAfter = pwrZ - minZafter;

        %get profiles
        kCH = find(channelsProfiles == chDG);
        %ratio to ch above
        ratioAbove = (profilesDS(kCH,:)./profilesDS(1,:))';
        ratioBelow = (profilesDS(kCH,:)./profilesDS(end,:))';

        %{
        DS features (updated Sept 12, 2017) 
        diffFromMinBefore > 2.5
        diffFromMinAfter > 2.5
        wDG 5-25
        ratioDecayUp > 1.5
        %}

        %k = diffFromMinBefore > 2 & diffFromMinAfter > 2 & wDG > 5 & wDG < 25;% & ratioAbove > 1.5;

        dB = feat(:,3) - feat(:,9); %distance to min before
        dA = feat(:,3) - feat(:,10); %distance to min before

        dB = log(dB); %log (gamma)
        dA = log(dA); %log (gamma)

        %mx = mean(dB);
        %sx = std(dB);
        dB = (dB-mxCSD)/sxCSD; %zscore
        dA = (dA-mxCSD)/sxCSD;

        %ratio of peak height vs height of maxima before
        fB = (feat(:,3) - feat(:,9)) ./ (feat(:,7) - feat(:,9));
        fA = (feat(:,3) - feat(:,10)) ./ (feat(:,8) - feat(:,10));

        kDS = dB > 0.75 & dA > 0.75 & wDG > 5 & wDG < 25;

        samplesDS = s(kDS);

        if length(samplesDS) ~= size(resultsCSD,2)
            disp('N DS ~= N CSD');
            return;
        end

        matches = zeros(4,size(resultsCSD,2));

        for i = 1:size(resultsCSD,2)

            x = resultsCSD(:,i)';

            %sources
            %kMax = find((x > [x(1) x(1:(end-1))]) & (x >= [x(2:end) x(end)])); 
            %kMax = kMax(kMax > 1 & kMax < 200);
            %xMax = x(kMax);

            %sinks
            kMin = find((x < [x(1) x(1:(end-1))]) & (x <= [x(2:end) x(end)]));
            kMin = kMin(kMin > 1 & kMin < 200);
            %xMin = x(kMin);
            
            %sinks (from top to bottom)
            for sI = 1:4
                
                s = sinks(sI);
                if isnan(s); continue; end;
                
                st = s-wLenSearch;
                ed = s+wLenSearch;
                if st < 1; st = 1; end;
                if ed > 200; ed = 200; end

                k = sum(kMin > st & kMin < ed);
                if k == 0; continue; end;

                matches(sI,i) = 1;
            end
        end

        %none
        %kNone = matches(1,:)==0 & matches(2,:)==0;
        %kBoth = matches(1,:)==1 & matches(2,:)==1;

        kType1sup = matches(1,:)==1 & matches(2,:)==0;
        kType2sup = matches(1,:)==0 & matches(2,:)==1;
        
        kType1inf = matches(4,:)==1 & matches(3,:)==0;
        kType2inf = matches(4,:)==0 & matches(3,:)==1;

        save([drOUT fnEEG],'kType1sup','kType1inf','kType2sup','kType2inf','kDS','samplesDS','matches');
    end 
end




