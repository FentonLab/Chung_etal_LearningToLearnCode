%convert stim file recorded using Intan chip (32ch)
 clear; close all; %clc;
%TODO
%check continuity of eegTS

Nsamples8k = 1342; %8k sampling rate per stim
Nsamples10k = 1166; %10k

drOUT = '../MAT_STIM/';

if exist(drOUT,'dir') == 0; mkdir(drOUT); end;

fileOverwrite = 0; %1=overwrite existing files

%animals
drAnimals = dir('data*');
for aI = 1:length(drAnimals)
    
    drIN = [drAnimals(aI).name '/'];

    %dirs
    d = dir(drIN);
    isub = [d(:).isdir]; %# returns logical vector
    dirs = {d(isub).name}';
    dirs(ismember(dirs,{'.','..'})) = [];
    dirs =1; % dummy value to index : 20180508
     for dI = 1:length(dirs)
    
%         dr = dirs{dI};
        dr = []; % dummy value to not inex any sub folders : 20180508
        %files
        files = dir([drIN dr '/*stim.bin']);
        if length(files) < 1; continue; end;

        for fI = 1:length(files)

            fn = files(fI).name;
%             fnMAT = [fn(1:end-4) '.mat']; 
             fnMAT = [drIN(6:end-1) '.mat']; % Change to folder name:20180508
            if fileOverwrite == 0 %skip if file exists
               if exist([drOUT fnMAT],'file') ~= 0; continue; end;
            end

            disp(fn)

            szFile = files(fI).bytes;

            %read header
            szHeader = 0;
            stimInterTime = []; %older formats don't have it, define as empty so it can be checked at save
            fid = fopen([drIN dr '/' fn],'r');
            for i = 1:50
                tline = fgetl(fid);
                szHeader = szHeader + length(tline) + 2; %add new line characters (2)

                %{
                %HEADER_START
                SoftwareVersion;1.5
                SoftwareYear;2016
                DateTime;2/10/2014 11:08:56 PM
                SamplingRate;10000
                Channels;32
                BytesPerSample;66  %4B + 4B + 8*4B = 40B (time,io,data)
                RecordsPerSample;33 %1*time + 1*IO + 8 = 10
                Records;Time,IO,CH1,CH2,CH3,CH4,CH5,CH6,CH7,CH8...
                StimInterTime;15
                StimLength;125
                StimRepeats;5
                StimLevels;0,25,50,75,100,125,150,175,200,225,250,275,300,
                %HEADER_END
                %}

                if strfind(tline,'SamplingRate') %sampling rate
                    str = strsplit(tline,';');
                    stimFS = str2double(str{2});
                elseif strfind(tline,'Gain') %gain
                    str = strsplit(tline,';');
                    gain = str2double(str{2});
                elseif strfind(tline,'Channels') %channels
                    str = strsplit(tline,';');
                    Nchannels = str2double(str{2});
                elseif strfind(tline,'BytesPerSample') %bytes per sample
                    str = strsplit(tline,';');
                    bytesPerSample = str2double(str{2});
                elseif strfind(tline,'RecordsPerSample') %records per sample
                    str = strsplit(tline,';');
                    recordsPerSample = str2double(str{2});
                elseif strfind(tline,'StimInterTime') %records per sample
                    str = strsplit(tline,';');
                    stimInterTime = str2double(str{2});
                elseif strfind(tline,'StimLength') %records per sample
                    str = strsplit(tline,';');
                    stimLength = str2double(str{2});
                elseif strfind(tline,'StimRepeats') %records per sample
                    str = strsplit(tline,';');
                    stimRepeats = str2double(str{2});
                elseif strfind(tline,'StimLevels') %records per sample
                    str = strsplit(tline,';');
                    str = str{2}; %right from ;
                    str = strsplit(str,',');
                    stimLevels = [];
                    for sI = 1:length(str)
                        s = str{sI};
                        if isempty(s); continue; end;
                        stimLevels = cat(2,stimLevels,str2double(s));
                    end
                elseif strcmp(tline,'%HEADER_END') %last header line
                    break;
                end; 
            end
            frewind(fid);
            fread(fid,szHeader,'uchar'); %skip header

            %number of samples EEG
            szEEG = szFile - szHeader;
            NsamplesEEG = szEEG / bytesPerSample;
            if rem(NsamplesEEG,1) ~= 0
                disp('Non-integer number of samples !')
            end

            data = fread(fid,[recordsPerSample,NsamplesEEG],'uint16'); %read data
            fclose(fid);

            stimTS = data(1,:);
            stimData = data(2:end,:);

            %convert eeg data to voltage
            stimData = double(stimData);

            stimData = stimData - 32768;
            stimData = stimData/32768;
            stimData = stimData * 5e-3; %full scale Intan +/-5mV

            %number of recordings
            if rem(size(stimData,2),Nsamples8k) == 0
                Nsamples = Nsamples8k;
            elseif rem(size(stimData,2),Nsamples10k) == 0
                Nsamples = Nsamples10k;
            else
                disp('Nsamples non integer!');
                continue;
            end

            Nrec = size(stimData,2)/Nsamples;

            %reshape data into channels x time x stims
            stimData = reshape(stimData,Nchannels,Nsamples,Nrec);

            %{
            %replace missing samples with nan
            d = diff(eegTS);
            u = unique(d);
            if length(u) > 1 %more differences = lost data
                disp([fnEEG ' sync problem, fixing...']);
                %eegData + eegTS
                Nsamples = eegTS(end)-eegTS(1)+1;
                eegTSnew = 1:Nsamples;
                eegDataNew = nan(size(eegData,1),Nsamples);
                knan = eegTS - eegTS(1) + 1; %convert to index
                for chI = 1:size(eegData,1)
                    x = eegData(chI,:);
                    eegDataNew(chI,knan) = eegData(chI,:);
                end
                eegData = eegDataNew;
                eegTS = eegTSnew;

                %eegIO = finish
                ds = diff(eegSync);
                us = unique(ds); 
            end
            %}


            %save
            fnOut = [drOUT fnMAT];
            if isempty(stimInterTime) %was not present
                save(fnOut,'stimData','stimFS','stimTS');
            else
                save(fnOut,'stimData','stimFS','stimTS','stimInterTime','stimLength','stimLength','stimRepeats','stimLevels');
            end

        end

     end
end



