% Author: Marco Buiatti, CIMeC (University of Trento, Italy), 2017.


function NEAR_performPostProcessing(EEG, origEEG)

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
    end
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


end