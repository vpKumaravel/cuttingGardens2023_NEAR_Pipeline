% NEAR_pipeline.m -NEAR Tutorial Script (Wrapper for executing all functions)
%
%
% Author: Velu Prabhakar Kumaravel (FBK; UNITN, Italy, 2018-2023)
% email: vpr.kumaravel@gmail.com


%% Example:

% params.dname = 's5.set'; % dataset name with extension
% params.dloc = '.\sample_data\NEAR_DATA'; % corresponding file location
% params.chanlocation_file = '.\sample_locs\GSN-HydroCel-129.sfp'; % channel location file location, if applicable
% params.segt_file = 'NEAR_LookingTimes.xlsx';
% params.segt_loc = '.\sample_data\NEAR_DATA';
% params.hpf_tf = [0.15 0.30]; % in Hz
% params.lpf_fc = 40; % in Hz
% params.isFlat = 1;
% params.FlatSec = 5;
% params.isLOF = 1;
% params.LOFThreshold = 2.5;
% params.LOFMethod = 'seuclidean';
% params.isLOFAdapt = 10;
% params.isPeriodogram = 0;
% params.isPlot = 1;
% params.ASRCutoff = 24;
% params.ASRMode = 'on';
% params.erp_event_markers = [];
% params.erp_epoch_duration = [];
% params.baseline_window = [];
% params.interpolType = 'spherical';
% 
% outEEG = NEAR_pipeline(params);


function outEEG = NEAR_performPreprocessing(params)

%% Define variables
dname = params.dname; % dataset name with extension
dloc = params.dloc; % corresponding file location
%chanlocation_file = params.chanlocation_file; % channel location file location, if applicable
segt_file = params.segt_file;
segt_loc = params.segt_loc;

%% Load data

EEG = pop_loadset('filename',dname,'filepath',[dloc filesep]); % import data

%% In case of ERP datasets (and already epoched), execute this

if(size(EEG.data, 3) > 1) % epoched data
   EEG = eeg_epoch2continuous(EEG); % make it continuous 
   isERP = 1; % set a flag
else
    isERP = 0;
end

%% Band-pass Filtering

EEG = clean_drifts(EEG,params.hpf_tf, []); % High pass filter with transition band [0.15 0.3] Hz
EEG = pop_eegfiltnew(EEG, [], params.lpf_fc, [], 0, [], 0); % Low pass filter @ 40 Hz

%% Segmentation based on visual attention

if ~isempty(segt_loc)
    lookFile = importdata([segt_loc filesep segt_file]);
    dname_split = split(dname, '.');
    dname_only = dname_split{1};
    lookTimes = NEAR_getLookTimes(lookFile,dname_only, 5000); % 's5' is the name of the sheet; 5000ms is the minimum looking time to consider for further analysis
    EEG = pop_select(EEG, 'time', lookTimes);
end

%% Bad Channel Detection using Local Outlier Factor (LOF)
origEEG = EEG; % making a copy for later use

[EEG, flat_ch, lof_ch, ~, ~, ~] = ...
    NEAR_getBadChannels(EEG, params.isFlat, params.FlatSec, ...
    params.isLOF, params.LOFThreshold, params.LOFMethod, params.isLOFAdapt,...
    params.isPeriodogram,[], [], [], [], params.isPlot);
badChans = sort(unique(union(flat_ch, lof_ch)));
EEG = pop_select(EEG, 'nochannel', badChans);
eeglab redraw

%% Bad Segments Rejection using Artifacts Subspace Reconstruction (ASR)

EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion','off',...
                        'ChannelCriterion','off','LineNoiseCriterion','off','Highpass','off',...
                        'BurstCriterion',params.ASRCutoff,'WindowCriterion','off',...
                        'BurstRejection',params.ASRMode,'Distance','Euclidian');

%% Only for ERP

if(isERP) % set isERP = 1 if your input data is ERP
   EEG = pop_epoch(EEG, params.erp_event_markers, ...
                    params.erp_epoch_duration, 'epochinfo', 'yes'); % markers (eg., {‘Eyes Close’, ‘Eyes Open’}), duration (in ms, eg., [0 1000]) to be defined by user
   EEG = pop_rmbase( EEG, params.baseline_window ,[]); % baseline window (in ms, e.g., [0 200]) to be defined by user
end

%% Interpolate missing channels

EEG = pop_interp(EEG, origEEG.chanlocs, params.interpolType); % spherical interpolation

%% Average re-referencing (the user may edit according to the needs)

EEG = pop_reref(EEG, [],'refloc',struct('labels',{'Cz'},'Y',{0},...
    'X',{5.4492e-16},'Z',{8.8992},'sph_theta',{0},'sph_phi',{90},'sph_radius',{8.8992},...
    'theta',{0},'radius',{0},'type',{''},'ref',{''},'urchan',{132},'datachan',{0}));

%% Save Preprocessed data


if exist([dloc filesep 'NEAR_Processed'], 'dir') == 0
    mkdir([dloc filesep 'NEAR_Processed'])
end

% Save data
outEEG = pop_saveset(EEG, 'filename',[erase(dname, '.set') '_NEAR_prep.set'],'filepath', [dloc filesep 'NEAR_Processed']);

%% Post-processing (Valid for SSVEP dataset)

NEAR_performPostProcessing(outEEG, origEEG);

