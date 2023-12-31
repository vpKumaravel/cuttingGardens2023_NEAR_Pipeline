% NEAR_steps.m -NEAR Tutorial Script
%
% See: https://github.com/vpKumaravel/vpKumaravel.github.io/wiki/Step-by-step-Tutorial-on-Newborns-EEG-Artifact-Removal-(NEAR)-pipeline
% Author: Velu Prabhakar Kumaravel (FBK; UNITN, Italy, 2018-2023)
% email: vpr.kumaravel@gmail.com
% Website: https://github.com/vpKumaravel/vpKumaravel.github.io/wiki/Step-by-step-Tutorial-on-Newborns-EEG-Artifact-Removal-(NEAR)-pipeline
% Feb 2022; 


%% Clear all variables and open EEGLAB

clc;
clear all;
close all;
eeglab;

%% Define variables
dname = 's5.set'; % dataset name with extension
dloc = '.\sample_data\NEAR_DATA'; % corresponding file location
%chanlocation_file = '.\sample_locs\GSN-HydroCel-129.sfp'; % channel location file location, if applicable
segt_file = 'NEAR_LookingTimes.xlsx';
segt_loc = '.\sample_data\NEAR_DATA';

%% Load data

EEG = pop_loadset('filename',dname,'filepath',[dloc filesep]); % import data
eeglab redraw

%% In case of ERP datasets (and already epoched), execute this

if(size(EEG.data, 3) > 1) % epoched data
   EEG = eeg_epoch2continuous(EEG); % make it continuous 
   isERP = 1; % set a flag
else
    isERP = 0;
end

%% Band-pass Filtering

EEG = clean_drifts(EEG,[0.15 0.3], []); % High pass filter with transition band [0.15 0.3] Hz
EEG = pop_eegfiltnew(EEG, [], 40, [], 0, [], 0); % Low pass filter @ 40 Hz
eeglab redraw

%% Segmentation based on visual attention

lookFile = importdata([segt_loc filesep segt_file]);
lookTimes = NEAR_getLookTimes(lookFile,'s5', 5000); % 's5' is the name of the sheet; 5 is the minimum looking time (in milliseconds) to consider for further analysis
EEG = pop_select(EEG, 'time', lookTimes);
origEEG = EEG; % making a copy for later use
eeglab redraw

%% Bad Channel Detection using Local Outlier Factor (LOF)

[EEG, flat_ch, lof_ch, periodo_ch, LOF_vec, thresh_lof_update] = NEAR_getBadChannels(EEG, 1, 5, 1, 2.5, 'seuclidean', 10, 0,[], [], [], [], 0);
badChans = sort(unique(union(flat_ch, lof_ch)));
EEG = pop_select(EEG, 'nochannel', badChans);
eeglab redraw

%% Bad Segments Rejection using Artifacts Subspace Reconstruction (ASR)

EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion','off','ChannelCriterion','off','LineNoiseCriterion','off', ...
    'Highpass','off','BurstCriterion',24,'WindowCriterion','off','BurstRejection','on','Distance','Euclidian');
eeglab redraw

%% Only for ERP

if(isERP) % set isERP = 1 if your input data is ERP
   EEG = pop_epoch( EEG, erp_event_markers, erp_epoch_duration, 'epochinfo', 'yes'); % markers (eg., {�Eyes Close�, �Eyes Open�}), duration (in ms, eg., [0 1000]) to be defined by user
   EEG = pop_rmbase( EEG, baseline_window ,[]); % baseline window (in ms, e.g., [0 200]) to be defined by user
end

%% Interpolate missing channels

EEG = pop_interp(EEG, origEEG.chanlocs, 'spherical'); % spherical interpolation
eeglab redraw

%% Average re-referencing (the user may edit according to the needs)

EEG = pop_reref( EEG, [],'refloc',struct('labels',{'Cz'},'Y',{0},...
    'X',{5.4492e-16},'Z',{8.8992},'sph_theta',{0},'sph_phi',{90},'sph_radius',{8.8992},...
    'theta',{0},'radius',{0},'type',{''},'ref',{''},'urchan',{132},'datachan',{0}));

eeglab redraw
%% Save Preprocessed data

if exist([dloc filesep 'NEAR_Processed'], 'dir') == 0
    mkdir([dloc filesep 'NEAR_Processed'])
end

% Save data
EEG = pop_saveset(EEG, 'filename',[erase(dname, '.set') '_NEAR_prep.set'],'filepath', [dloc filesep 'NEAR_Processed']);

%% Post-processing

% Author: Marco Buiatti, CIMeC (University of Trento, Italy), 2017.

type = 'DIN'; 

ev=1; % event counter
epoch=1;
while ev <= length(EEG.event)
    %     while strcmp(EEG.event(ev).type,'boundary') % if event(ev) is 'boundary', we test whether it is the beginning of the 'type' segment
    %         disp('boundary');
    if any(strcmp(EEG.event(ev).type,type)) % if event(ev+1) is 'type', ev=first time point
    %if EEG.event(ev).type == type % if event(ev+1) is 'type', ev=first time point
        % if ev is first event, there is no initial boundary, therefore the
        % epoch starts at time 0
        if ev==1 pointrange{epoch}(1)=0.004;
        else
            pointrange{epoch}(1)=EEG.event(ev-1).latency;
        end
        ev=ev+1;
        a=1;
        while a
            if ev == length(EEG.event) % if ev=last event
                pointrange{epoch}(2)=EEG.pnts; % ev=last time point of the dataset
                a=0;
            elseif ev == length(EEG.event)+1 % if ev=last event
                pointrange{epoch}(2)=EEG.pnts; % ev=last time point of the dataset
                a=0;
            elseif ~strcmp(EEG.event(ev).type,'boundary') % if event(ev) is not 'boundary', proceed to next ev
                ev=ev+1;
            else % if event(ev) is 'boundary'
                pointrange{epoch}(2)=EEG.event(ev).latency; % event(ev)=last time point
                a=0;
            end
        end
        % select data
        EEGsel = pop_select( EEG,'point',pointrange{epoch});
        dataNEAR{epoch}=EEGsel.data;
        epoch=epoch+1;
    end;
    ev=ev+1;
end

% For unprocessed (filtered, segmented) data - Repeat the computations

ev=1; % event counter
epoch=1;
while ev <= length(origEEG.event)
    %     while strcmp(origEEG.event(ev).type,'boundary') % if event(ev) is 'boundary', we test whether it is the beginning of the 'type' segment
    %         disp('boundary');
    if any(strcmp(origEEG.event(ev).type,type)) % if event(ev+1) is 'type', ev=first time point
    %if origEEG.event(ev).type == type % if event(ev+1) is 'type', ev=first time point
        % if ev is first event, there is no initial boundary, therefore the
        % epoch starts at time 0
        if ev==1 pointrange{epoch}(1)=0.004;
        else
            pointrange{epoch}(1)=origEEG.event(ev-1).latency;
        end
        ev=ev+1;
        a=1;
        while a
            if ev == length(origEEG.event) % if ev=last event
                pointrange{epoch}(2)=origEEG.pnts; % ev=last time point of the dataset
                a=0;
            elseif ev == length(origEEG.event)+1 % if ev=last event
                pointrange{epoch}(2)=origEEG.pnts; % ev=last time point of the dataset
                a=0;
            elseif ~strcmp(origEEG.event(ev).type,'boundary') % if event(ev) is not 'boundary', proceed to next ev
                ev=ev+1;
            else % if event(ev) is 'boundary'
                pointrange{epoch}(2)=origEEG.event(ev).latency; % event(ev)=last time point
                a=0;
            end
        end
        % select data
        EEGsel = pop_select( origEEG,'point',pointrange{epoch});
        dataRaw{epoch}=EEGsel.data;
        epoch=epoch+1;
    end
    ev=ev+1;
end


% Author: Marco Buiatti, CIMeC (University of Trento, Italy), 2016-.

orig_data = dataNEAR; %backup copy
wl = 2500;
srate = EEG.srate;
loc = [];

for ep=1:length(dataNEAR)
    loc(ep)=length(dataNEAR{ep}(1,:))/wl;
    disp(loc(ep));
end

dataNEAR(loc<1)=[];

if isempty(dataNEAR) %If there's not enough data, zero-pad
    count = 1; % Count of valid epochs, an epoch is valid with at least 625 samples
    for ep=1:length(orig_data)
        if(length(orig_data{ep}) >= wl/2) % It can be modified to 313 if required
            data_pad{count} = [orig_data{ep} zeros(size(orig_data{ep},1), wl - size(orig_data{ep},2))]; %Zero-padding
            count = count + 1;
        end
    end
    dataNEAR= data_pad;
    disp('Window longer than data! But we got it by zero-padding the data!');
end

% Taper option (default taper is square) 
option='s';
if  strcmp(option, 's')
   w=ones(wl,1)/wl;
elseif  strcmp(option, 'w')
    w=welch(wl);
elseif  strcmp(option, 'p')
   w=parzen(wl);
elseif  strcmp(option, 'h')
   w=hamming(wl);
end
w=w';

wss=wl*(sum(w.^2));	%window squared and summed

for el=1:size(dataNEAR{1},1)
    ap=0;
    for ep=1:length(dataNEAR)
        l=size(dataNEAR{ep}(el,:),2);
        if l-wl<floor(wl/4)  % use one window only if residual segment is shorter than wl/4
            l=wl;
            N=ceil(2*l/wl); % max number of consecutive HALF windows (non-overlapping if l is multiple of wl/2)
            dwl=0;
        else
            N=ceil(2*l/wl); % max number of consecutive HALF windows (non-overlapping if l is multiple of wl/2)
            dwl=floor((l-wl)/(N-2)); % step
        end
        n_loc(ep)=N-1; % number of consecutive full windows
        
        ap_loc=0;						%average power
        for i=1:n_loc(ep)
            wd=w.*dataNEAR{ep}(el,(i-1)*dwl+1:(i-1)*dwl+wl);
            ap_loc=ap_loc + abs(fft(wd)).^2;
        end
        ap=ap+ap_loc;
    end
    n=sum(n_loc);
    psNEAR(el,1)=ap(1)/(wss*n);
    psNEAR(el,2:wl/2)=(ap(2:wl/2)+ap(wl:-1:wl/2 +2))/(wss*n);
    psNEAR(el,wl/2 +1)=ap(wl/2 +1)/(wss*n);
    interval=0:1:(wl/2);
end
f=interval*srate/wl;

for ep=1:length(n_loc)
            disp(['Number of windows used for epoch ' num2str(ep) ' = ' num2str(n_loc(ep))]);
end

% Repeat this for Raw Data
orig_data = dataRaw; %backup copy
wl = 2500;
srate = EEG.srate;
loc = [];

for ep=1:length(dataRaw)
    loc(ep)=length(dataRaw{ep}(1,:))/wl;
    disp(loc(ep));
end

dataRaw(loc<1)=[];
data_pad = [];
if isempty(dataRaw) %If there's not enough data, zero-pad
    count = 1; % Count of valid epochs, an epoch is valid with at least 625 samples
    for ep=1:length(orig_data)
        if(length(orig_data{ep}) >= wl/2) % It can be modified to 313 if required
            data_pad{count} = [orig_data{ep} zeros(size(orig_data{ep},1), wl - size(orig_data{ep},2))]; %Zero-padding
            count = count + 1;
        end
    end
    dataRaw = data_pad;
    disp('Window longer than data! But we got it by zero-padding the data!');
end

% Taper option (default taper is square) 
option='s';
if  strcmp(option, 's')
   w=ones(wl,1)/wl;
elseif  strcmp(option, 'w')
    w=welch(wl);
elseif  strcmp(option, 'p')
   w=parzen(wl);
elseif  strcmp(option, 'h')
   w=hamming(wl);
end
w=w';

wss=wl*(sum(w.^2));	%window squared and summed

for el=1:size(dataRaw{1},1)
    ap=0;
    for ep=1:length(dataRaw)
        l=size(dataRaw{ep}(el,:),2);
        if l-wl<floor(wl/4)  % use one window only if residual segment is shorter than wl/4
            l=wl;
            N=ceil(2*l/wl); % max number of consecutive HALF windows (non-overlapping if l is multiple of wl/2)
            dwl=0;
        else
            N=ceil(2*l/wl); % max number of consecutive HALF windows (non-overlapping if l is multiple of wl/2)
            dwl=floor((l-wl)/(N-2)); % step
        end
        n_loc(ep)=N-1; % number of consecutive full windows
        
        ap_loc=0;						%average power
        for i=1:n_loc(ep)
            wd=w.*dataRaw{ep}(el,(i-1)*dwl+1:(i-1)*dwl+wl);
            ap_loc=ap_loc + abs(fft(wd)).^2;
        end
        ap=ap+ap_loc;
    end
    n=sum(n_loc);
    psRaw(el,1)=ap(1)/(wss*n);
    psRaw(el,2:wl/2)=(ap(2:wl/2)+ap(wl:-1:wl/2 +2))/(wss*n);
    psRaw(el,wl/2 +1)=ap(wl/2 +1)/(wss*n);
    interval=0:1:(wl/2);
end
f=interval*srate/wl;

for ep=1:length(n_loc)
            disp(['Number of windows used for epoch ' num2str(ep) ' = ' num2str(n_loc(ep))]);
end



ch = [62,70,71,72,75,76,77,78,83,84,90]; % channel indices from the ROI 
frange = 6:20; % frequency bins to plot
p = {psRaw};
figure; 
for i=1:length(p)
        plot(f(frange),mean(p{i}(ch,frange),1),'-*','linewidth',3);
    hold on;
end

axis tight
set(gca,'Fontsize',16);
title('Raw');
xlabel('Frequency (Hz)');
ylabel('Power Spectral Density \mu V^2');

p = {psNEAR};
figure; 
for i=1:length(p)
        plot(f(frange),mean(p{i}(ch,frange),1),'-*','linewidth',3);
    hold on;
end

axis tight
set(gca,'Fontsize',16);
ylim([0 24])
title('NEAR');
xlabel('Frequency (Hz)');
ylabel('Power Spectral Density \mu V^2');

