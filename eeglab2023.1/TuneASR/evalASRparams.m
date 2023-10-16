%%
% NEAR Pipeline Evaluation
%
% Function to evaluate all FTR values - returns ftr arrays and error log
%
% loc_in = input file location
% fname_in = file name
% k_in_array = list of K values for ASR (a.k.a., asr cut-off parameter)
% process_array = list of "Processing" values - on or off - length should be equal to k_in_array
% 'ON'  = 'ASR Removal'
% 'OFF' = 'ASR Correction' 
% rawEEG = EEG structure of a similar data but without removal of bad channels
%
% Velu Prabhakar Kumaravel, FBK/CIMeC (UNITN), Italy

function [metric, error_log] = evalASRparams(loc_in, fname_in, k_in_array, process_array)


EEG = pop_loadset('filename', fname_in,'filepath', loc_in);
metric = zeros(length(k_in_array), 1);

origEEG = load('chanlocs.mat'); % to interpolate removed channels
for p = 1:length(k_in_array)
    
    disp('Add other relevant preprocessing measures');
    
    EEG_current = pop_clean_rawdata(EEG, 'FlatlineCriterion','off', ...
        'ChannelCriterion','off','LineNoiseCriterion','off','Highpass',...
        'off','BurstCriterion',k_in_array(p),'WindowCriterion','off',...
        'BurstRejection',process_array{p},'Distance','Euclidian');
    
    EEG_current = pop_interp(EEG_current, origEEG.chanlocs, 'spherical'); % spherical interpolation
    EEG_current = pop_reref(EEG_current, [],'refloc',struct('labels',{'Cz'},'Y',{0},...
    'X',{5.4492e-16},'Z',{8.8992},'sph_theta',{0},'sph_phi',{90},'sph_radius',{8.8992},...
    'theta',{0},'radius',{0},'type',{''},'ref',{''},'urchan',{132},'datachan',{0}));

    try 
        
        disp('Use EEG1 to perform your post-processing here: epoching, computing psd measures etc.,');

        metric(p) = computeMetric(EEG_current); % Compute your measure using EEG_current and save it here.
        error_log{p} = 'Success';
        
        
    catch Error
        error_log{p} = Error.message;
        disp(Error);
        metric(p) = 0;
    end
    
end

function signal_to_noise_ratio = computeMetric(EEG)


    ch = [62,70,71,72,75,76,77,78,83,84,90]; % channel indices from the ROI
    signal_idx = 9; % 0.8 Hz
    noise_idx = [6:8, 10:12];
    
    ps = performPostProcessing(EEG);
    
    ftr = computeFTR(ps, ch, noise_idx, signal_idx);
    signal_to_noise_ratio = round(ftr *100)/100;



function [ftr] = computeFTR(ps, selch, noise_indices, signal_index)


    signal = mean(ps(selch,signal_index)); % mean across the channels for the given frequency index
    noise = mean(mean(ps(selch, noise_indices))); % mean across the channels for the mean of given frequency indices
    ftr = signal/noise;
    
function ps = performPostProcessing(EEG)

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
        evData{epoch}=EEGsel.data;
        epoch=epoch+1;
    end
    ev=ev+1;
end

orig_data = evData; %backup copy
wl = 2500;
srate = EEG.srate;
loc = [];

for ep=1:length(evData)
    loc(ep)=length(evData{ep}(1,:))/wl;
end

evData(loc<1)=[];
data_pad = [];
if isempty(evData) %If there's not enough data, zero-pad
    count = 1; % Count of valid epochs, an epoch is valid with at least 625 samples
    for ep=1:length(orig_data)
        if(length(orig_data{ep}) >= wl/2) % It can be modified to 313 if required
            data_pad{count} = [orig_data{ep} zeros(size(orig_data{ep},1), wl - size(orig_data{ep},2))]; %Zero-padding
            count = count + 1;
        end
    end
    if isempty(data_pad)
        disp('Not enough samples to employ Zero-padding. FTR metric will be set to Zero');
        return
    else
        disp('Window longer than data! Zero-padding technique is employed!');
        evData= data_pad;
    end
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

for el=1:size(evData{1},1)
    ap=0;
    for ep=1:length(evData)
        l=size(evData{ep}(el,:),2);
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
            wd=w.*evData{ep}(el,(i-1)*dwl+1:(i-1)*dwl+wl);
            ap_loc=ap_loc + abs(fft(wd)).^2;
        end
        ap=ap+ap_loc;
    end
    n=sum(n_loc);
    ps(el,1)=ap(1)/(wss*n);
    ps(el,2:wl/2)=(ap(2:wl/2)+ap(wl:-1:wl/2 +2))/(wss*n);
    ps(el,wl/2 +1)=ap(wl/2 +1)/(wss*n);
    interval=0:1:(wl/2);
end
f=interval*srate/wl;

for ep=1:length(n_loc)
    disp(['Number of windows used for epoch ' num2str(ep) ' = ' num2str(n_loc(ep))]);
end




