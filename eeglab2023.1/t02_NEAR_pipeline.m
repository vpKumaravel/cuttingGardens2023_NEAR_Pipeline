params.dname = 's5.set'; % dataset name with extension
params.dloc = '.\sample_data\NEAR_DATA'; % corresponding file location
params.chanlocation_file = '.\sample_locs\GSN-HydroCel-129.sfp'; % channel location file location, if applicable
params.segt_file = 'NEAR_LookingTimes.xlsx';
params.segt_loc = '.\sample_data\NEAR_DATA';
params.hpf_tf = [0.15 0.30]; % in Hz
params.lpf_fc = 40; % in Hz
params.isFlat = 1; % if you want to remove flat channels
params.FlatSec = 5; % if so, what is the criteria in seconds
params.isLOF = 1; % if you want to employ LOF
params.LOFThreshold = 2.5; % LOF threshold
params.LOFMethod = 'seuclidean'; % Distance metric
params.isLOFAdapt = 10; % do you want to employ adaptive thresholding? if so, put % of channels you want to remove
params.isPeriodogram = 0; % set to 1 if you want to perform spectral analysis (optional)
params.isPlot = 1; % to plot intermediate results (if any)
params.ASRCutoff = 24; % ASR cut-off Parameter
params.ASRMode = 'on'; % ASR mode (rejection = 'on'; correction = 'off')
params.erp_event_markers = []; % ERP event markers
params.erp_epoch_duration = []; % ERP epoch duration
params.baseline_window = []; % ERP baseline correction
params.interpolType = 'spherical'; % Bad channel interpolation

outEEG = NEAR_performPreprocessing(params);