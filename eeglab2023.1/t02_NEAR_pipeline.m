params.dname = 's5.set'; % dataset name with extension
params.dloc = '.\sample_data\NEAR_DATA'; % corresponding file location
params.chanlocation_file = '.\sample_locs\GSN-HydroCel-129.sfp'; % channel location file location, if applicable
params.segt_file = 'NEAR_LookingTimes.xlsx';
params.segt_loc = '.\sample_data\NEAR_DATA';
params.hpf_tf = [0.15 0.30]; % in Hz
params.lpf_fc = 40; % in Hz
params.isFlat = 1;
params.FlatSec = 5;
params.isLOF = 1;
params.LOFThreshold = 2.5;
params.LOFMethod = 'seuclidean';
params.isLOFAdapt = 10;
params.isPeriodogram = 0;
params.isPlot = 1;
params.ASRCutoff = 24;
params.ASRMode = 'on';
params.erp_event_markers = [];
params.erp_epoch_duration = [];
params.baseline_window = [];
params.interpolType = 'spherical';

outEEG = NEAR_performPreprocessing(params);