%get CSD of DS
%clear; close all; clc;
function getDSCSD(indices)

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

drDS = '/scratch/dd1348/M32/DS_detected_nolimit/';

drOUT = '/scratch/dd1348/M32/DS_CSD/';

drMAT = '/scratch/dd1348/M32/MAT_EEG/';
drART = '/scratch/dd1348/M32/MAT_ART/';

fnList = 'fileList_v2.xlsx';
[~,~,RAW]=xlsread(fnList);
listBase = RAW(2:end,1); %base
listExp = RAW(2:end,2); %exp
list2HR = RAW(2:end,3); %2hr
listChannels = RAW(2:end,4:12); %channels
listID = RAW(2:end,13);
listGen = RAW(2:end,14);
%listSinks = RAW(2:end,22) %alternative sinks
listDG = cell2mat(listChannels(:,7));


%create filter for CSD filtering
eegFS = 2000;
[b, gd] = getFIRbandpass( 5,200,4,40,eegFS );

winLen = 10;

figure('Position',[100 100 200 800]);

IDs = unique(listID);

%animals
for i = indices
    ID = IDs{i};
    
    k = strcmp(listID,ID);
    filesBase = listBase(k);
    filesExp = listExp(k);
    files2HR = list2HR(k);
    channelsDG = listDG(k);
    chan = listChannels(k,:);
    
    %load all features, select threshold for DS
    featALL = [];
    for fI = 1:length(filesBase)
        fnEEG = filesBase{fI};
        if exist([drDS fnEEG],'file') == 0; continue; end;
        if exist([drART fnEEG],'file') == 0; continue; end;
        load([drDS fnEEG]);
        load([drART fnEEG]);
        
        chDG = channelsDG(fI);
        
        %art
        art = signalOK{chDG};

        %labelsFeat = {'sample','promZ','peakZ','peakRaw','promZbefore','promZafter',... %1..6
        %    'peakZbefore','peakZafter','minZbefore','minZafter',...%7..10
        %    'wDG','wBA','w50','tB','tA','tminB','tminA','acc','speed','max10_30before','max10_30after'}; %11..21

        s = feat(:,1); %sample
    
        %exclude art
        k = art(s);
        feat = feat(k,:);
        
        featALL = cat(1,featALL,feat);
    end
    
    %s = featALL(:,1); %sample
    wDG = featALL(:,11); %width

    dB = featALL(:,3) - featALL(:,9); %distance to min before
    dA = featALL(:,3) - featALL(:,10); %distance to min before

    dB = log(dB); %log (gamma)
    dA = log(dA); %log (gamma)

    %get normalization values
    mxCSD = mean(dB);
    sxCSD = std(dB);
    
    for typeI = 1:3
        if typeI == 1
            files = filesBase;
        elseif typeI == 2
            files = filesExp;
        elseif typeI == 3
            files = files2HR;
        end
    
        %process all files individually
        for fI = 1:length(files)
            fnEEG = files{fI};
            chDG = channelsDG(fI);

            channels = chan(fI,:);
            channels = cell2mat(channels);

            if isnan(fnEEG); continue; end;

            if exist([drOUT fnEEG],'file') ~= 0; continue; end;

            if exist([drMAT fnEEG],'file') == 0; continue; end;
            if exist([drDS fnEEG],'file') == 0; continue; end;
            if exist([drART fnEEG],'file') == 0; continue; end;

            disp(fnEEG)

            load([drMAT fnEEG]);
            load([drDS fnEEG]);
            load([drART fnEEG]);

            Nchannels = size(eegData,1);

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

            %art
            art = signalOK{chDG};

            %labelsFeat = {'sample','promZ','peakZ','peakRaw','promZbefore','promZafter',... %1..6
            %    'peakZbefore','peakZafter','minZbefore','minZafter',...%7..10
            %    'wDG','wBA','w50','tB','tA','tminB','tminA','acc','speed','max10_30before','max10_30after'}; %11..21

            s = feat(:,1); %sample

            %exclude art
            k = art(s);
            feat = feat(k,:);

            s = feat(:,1); %sample
            %pwrZ = feat(:,3); %peakZ
            %minZbefore = feat(:,9); %minZbefore
            %minZafter = feat(:,10); %minZafter
            wDG = feat(:,11); %width

            dB = feat(:,3) - feat(:,9); %distance to min before
            dA = feat(:,3) - feat(:,10); %distance to min before

            dB = log(dB); %log (gamma)
            dA = log(dA); %log (gamma)

            %processed for all files
            %mx = mean(dB);
            %sx = std(dB);
            dB = (dB-mxCSD)/sxCSD; %zscore
            dA = (dA-mxCSD)/sxCSD;

            k = dB > 0.75 & dA > 0.75 & wDG > 5 & wDG < 25;

            sk = s(k);

            resultsCSD = nan(CSDres,length(sk));

            %get CSDs
            for sI = 1:length(sk)

                s = sk(sI);

                eeg = eegData(:,s-winLen:s+winLen);

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

                CSDp = CSD(:,winLen+1);

                resultsCSD(:,sI) = CSDp;

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

            save([drOUT fnEEG],'resultsCSD','chSel','labelsSel','sk','mxCSD','sxCSD');

        end %files
    end
end



