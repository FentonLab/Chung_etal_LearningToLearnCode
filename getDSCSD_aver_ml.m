%COMPUTE CSDS AND EXTRACT:
%1) ALL CSDs (time x space for each DS)
%2) Average DS
%3) DS in SINK level
function getDSCSD_aver_ml(fI)

runLocal = 0;

if runLocal == 0 %linux
    %addpath('/mnt/kant/dinod/matlab/my_functions/');
    %addpath('/mnt/kant/dinod/matlab/OTC/');
    addpath('/home/dd1348/matlab/my_functions/');
    addpath('/home/dd1348/matlab/OTC/');
    addpath('/home/dd1348/matlab/CSDplotter/');
    addpath('/home/dd1348/matlab/CSDplotter/methods/');
else %local
    addpath('/Volumes/DATA/matlab/my_functions/');
    addpath('/Volumes/DATA/matlab/OTC/');
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

drTYPE = '/scratch/dd1348/M32/DS_TYPE12/';
drMAT = '/scratch/dd1348/M32/MAT_EEG/';


drCSDAVER = '/scratch/dd1348/M32/DS_CSD_AVER/';
drCSDML = '/scratch/dd1348/M32/DS_CSD_ML/';
drCSDMLAVER = '/scratch/dd1348/M32/DS_CSD_ML_AVER/';

%{
drCSDAVER = '/scratch/dd1348/M32/DS_CSD_AVER_200Hz/';
drCSDML = '/scratch/dd1348/M32/DS_CSD_ML_200Hz/';
drCSDMLAVER = '/scratch/dd1348/M32/DS_CSD_ML_AVER_200Hz/';
%}

fnList = 'fileList_v2.xlsx';
[~,~,RAW]=xlsread(fnList);
listBase = RAW(2:end,1); %file name
listExp = RAW(2:end,2); %file name
list2Hr = RAW(2:end,3); %file name
listChannels = RAW(2:end,4:12); %channels
%listID = RAW(2:end,13);
%listGen = RAW(2:end,14);
listSinks = RAW(2:end,22);

%create filter for detection
eegFS = 2000;
[b, gd] = getFIRbandpass( 5,100,4,40,eegFS ); %100Hz
%[b, gd] = getFIRbandpass( 5,200,4,40,eegFS ); %200Hz

winLenHalf = 0.05*eegFS;

%disp(['found files ' num2str(length(listBase))]);
   
%base,exp,2hr
for dataSetI = 1:3
    
    if dataSetI == 1
        fnEEG = listBase{fI};
    elseif dataSetI == 2
        fnEEG = listExp{fI};
    else
        fnEEG = list2Hr{fI};
    end
    
    channels = listChannels(fI,:);
    channels = cell2mat(channels);
    %chDG = channels(7);
        
    %if ~(contains(fnEEG,'PRETRAIN_BASE') || contains(fnEEG,'PRE_BASE'))
    %    continue;
    %end
    
    chSel = [];
    labelsSel = [];
    zs = [];

    %if isnan(chFIS); return; end;
    %if isnan(chDGmi); chDGmi = 32; end;

    if isnan(fnEEG); continue; end;

    if exist([drTYPE fnEEG],'file') == 0; continue; end;
    if exist([drMAT fnEEG],'file') == 0; continue; end;
    %if exist([drOUT fnEEG],'file') ~= 0; continue; end;

    disp(fnEEG)

    load([drMAT fnEEG]);
    load([drTYPE fnEEG]);

    Nchannels = size(eegData,1);
    
    %get sinks
    sinks = listSinks{fI};
    if isnan(sinks); continue; end;
    if isempty(sinks); continue; end;
    s = strsplit(sinks,',');
    sinks = nan(1,length(s));
    for i = 1:length(s)
        sinks(i) = str2double(s{i});
    end

    if Nchannels == 32
        el_pos = 0.1:0.05:1.65; %our tetrode distances
        el_pos = el_pos*1e-3; %mm
    elseif Nchannels == 16
        el_pos = 0.1:0.05:0.85; %our tetrode distances
        el_pos = el_pos*1e-3; %mm
    end

    %filter
    eegF = filter(b,1,eegData,[],2);
    eegF = cat(2,eegF(:,gd+1:end),zeros(size(eegF,1),gd));

    eegData = eegF;

    resultsCSDaver = cell(2,2);
    resultsEEGaver = cell(2,2);
    resultsCSDML = cell(2,2); %sup/inf DS1,2
    resultsCSDMLaver = cell(2,2); %sup/inf DS1,2
    
    %superior/inferior sites
    for supinfI = 1:2
        
        %sup
        if supinfI == 1
            %samples of type1 and type2
            samplesDS1 = samplesDS(kType1sup);
            samplesDS2 = samplesDS(kType2sup);
        else %inf
            %samples of type1 and type2
            samplesDS1 = samplesDS(kType1inf);
            samplesDS2 = samplesDS(kType2inf);
        end

        for typeI = 1:2

            if typeI == 1
                sDS = samplesDS1;
            else
                sDS = samplesDS2;
            end
            
            %extract the sink
            %SUP1,SUP2,INF2,INF1
            if supinfI == 1 %SUP
               if typeI == 1 %DS1
                   sink = sinks(1);
               else %DS2
                   sink = sinks(2);
               end
            else %INF
               if typeI == 1 %DS1
                   sink = sinks(4);
               else %DS2
                   sink = sinks(3);
               end 
            end

            %{
            %select Nmax
            if length(sDS) > Nmax
                k = randsample(length(sDS),Nmax);
                sDS = sDS(k);
            end
            %}

            resultsCSD = zeros(CSDres,winLenHalf*2+1);
            resultsEEG = zeros(size(eegData,1),winLenHalf*2+1);
            resultsCSDMLx = nan(length(sDS),winLenHalf*2+1); %sink CSD responses
            counter = 0;

            %samples
            for sI = 1:length(sDS)

                s = sDS(sI);
                
                eeg = eegData(:,s-winLenHalf:s+winLenHalf);

                if sum(isnan(eeg(:))) > 0; continue; end;
                
                if methodSpline == 0
                    %compute CSD - standard method
                    %adding vaknin electrodes not to loose first and last channel
                    %duplicate first and last channel
                    eeg = cat(1,eeg(1,:),eeg,eeg(end,:));
                    eeg = eeg + eeg;

                    Nch = size(eeg,1);
                    CSD = -ex_cond*D1(Nch,h)*eeeegAvergX;

                    %filter iCSD
                    b0 = 0.54; %center
                    b1 = 0.23; %neighbor
                    [n1,n2]=size(CSD);            
                    CSD_add(1,:) = zeros(1,n2);   %add top and buttom row with zeros
                    CSD_add(n1+2,:)=zeros(1,n2);
                    CSD_add(2:n1+1,:)=CSD;        %CSD_add has n1+2 rows
                    CSD = S_general(n1+2,b0,b1)*CSD_add; % CSD has n1 row

                else
                    % compute spline iCSD:
                    Fcs = F_cubic_spline(el_pos,diam,ex_cond,top_cond);
                    [zs,CSD_cs] = make_cubic_splines(el_pos,eeg,Fcs);
                    if gauss_sigma~=0 %filter iCSD
                      [~,CSD_cs]=gaussian_filtering(zs,CSD_cs,gauss_sigma,filter_range);
                    end
                    CSD = CSD_cs;
                end

                if sum(isnan(CSD(:))) > 0; continue; end;

                resultsCSD = resultsCSD + CSD;
                resultsEEG = resultsEEG + eeg;
                counter = counter + 1;

                
                %save CSD at sink level
                if ~isnan(sink)
                    resultsCSDMLx(sI,:) = CSD(sink,:);
                end


                %if spline is used, find electrodes
                chs = 1:Nchannels;
                chsInd = zeros(1,length(chs));
                for chI = 1:Nchannels
                    [~,k] = min(abs(el_pos(chI) - zs));
                    chsInd(chI) = k;
                end

                %channel labels
                knan = ~isnan(channels);
                chSel = chsInd(channels(knan));
                labels = {'SP','SR','SLM','FIS','mDG','sDG','DG','iDG','mDG'};
                labelsSel = labels(knan);
            end

            resultsCSD = resultsCSD / counter;
            resultsEEG = resultsEEG / counter;

            resultsCSDaver{supinfI,typeI} = resultsCSD;
            resultsEEGaver{supinfI,typeI} = resultsEEG;
            resultsCSDML{supinfI,typeI} = resultsCSDMLx;
            resultsCSDMLaver{supinfI,typeI} = nanmean(resultsCSDMLx,1);
        end
    end

    %store average
    save([drCSDAVER fnEEG],'resultsCSDaver','resultsEEGaver','chSel','labelsSel','el_pos','zs');
   
    %store CSD ML
    save([drCSDML fnEEG],'resultsCSDML');
    
    %store CSD ML aver
    save([drCSDMLAVER fnEEG],'resultsCSDMLaver');
    
end

