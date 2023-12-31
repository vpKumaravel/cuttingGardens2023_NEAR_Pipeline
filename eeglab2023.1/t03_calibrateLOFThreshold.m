% 1) please ensure eeglab is accessible by MATLAB 
% 2) please ensure that NEAR_ChannelRejection-master folder is put inside
% the plugins folder of EEGLAB

clc; clear all; eeglab; 
%%
list_sub = {'s5_filt_segt' }; % add more files as comma separated values
ext = '.set';

list_labels = {'s5_labels'}; % add more files as comma separated values
filepath = '.\TuneLOF';
labelpath = filepath;


LOF_calibrate = []; % output struct file
counter = 1; % counter for output struct file

for eachfile = 1:numel(list_sub)


    %% Step 1) Import code - please modify this according to the filetype supported by EEGLAB
    EEG = pop_loadset('filename', [list_sub{eachfile} ext] ,'filepath', [filepath filesep]);
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    EEG = eeg_checkset( EEG );
    
    
    %% Step 2a) Add channel locations and other required preprocessing steps such as Filtering
    
    % Since the example data is already imported with channel locations,
    % skipping the step

    %% Step 3) Extract ground truth labels
    
    
    g_t = zeros(1,EEG.nbchan); % declaring ground truth vector
    disp(eachfile);
    T = readtable([labelpath filesep list_labels{eachfile} '.xlsx']); % read label file corresponding to the current EEG file

    groundTruth = T.labels'; % where x0, in this example, contains the list of 0s and 1s in the .csv file. (1 indicates bad channels)
    g_t(find(groundTruth)) = 1;

    %% Step 4) Perform NEAR Bad Channel Detection

    list_threshold = 1:0.5:5; % LOF optimal threshold grid-search
    
    for iThreshold = list_threshold % list of threshold values you want to explore. N.B. The lower limit should be >= 1.
        
        fprintf('\n Current Threshold is %f\n', iThreshold);
        p_t = zeros(1,EEG.nbchan); % declaring predicted labels vector
        
        [~, flat_ch, lof_ch, periodo_ch, ~] = NEAR_getBadChannels(EEG, 1, 5, 1, iThreshold, 'seuclidean',[], ...
            0, [], [], [], [], 0);
        NEARbadChans = sort(unique(union(union(flat_ch, lof_ch),periodo_ch)));
        p_t(NEARbadChans) = 1;
        
        % compute classification performance
        C = confusionmat(g_t,p_t);
        
        TN = C(1,1);
        FN = C(2,1);
        FP = C(1,2);
        TP = C(2,2);
        
        Sensitivity = (TP/(TP+FN));
        Specificity = (TN/(TN+FP));
        Accuracy = ((TN+TP)/(TN+TP+FN+FP));
        Precision = (TP/(TP+FP));
        Recall = Sensitivity;
        F1_scoreNEAR = (2*Precision*Recall)/(Precision+Recall);
        if(isnan(F1_scoreNEAR))
            F1_scoreNEAR = 0.0;
        end
        fprintf('F1 score is %f\n', F1_scoreNEAR);
        
        % Saving the metrics
        LOF_calibrate(counter).Subject = {list_sub{eachfile}};
        LOF_calibrate(counter).GroundTruth = find(g_t);
        LOF_calibrate(counter).LOF = iThreshold;
        LOF_calibrate(counter).NEAR  = NEARbadChans';
        LOF_calibrate(counter).Precision = Precision;
        LOF_calibrate(counter).Recall = Recall;
        LOF_calibrate(counter).F1_NEAR = F1_scoreNEAR;
        
        counter = counter + 1;
    end
    
end

    %% Step 5: Plot results

    vec = [LOF_calibrate.LOF];
    count = 1;
    mean_LOF = [];
    for iT = list_threshold
        idx = vec == iT;
        out = LOF_calibrate(idx);
        mean_LOF(count) = mean([out.F1_NEAR]);
        count = count +1;
        idx = [];
        out = [];
    end
    figure;plot(list_threshold,mean_LOF);
    xlabel('LOF Threshold','fontweight','bold','fontsize',24);
    ylabel('F1 Score','fontweight','bold','fontsize',24);
    set(get(gca, 'XAxis'), 'FontWeight', 'bold','fontsize',24);
    set(get(gca, 'YAxis'), 'FontWeight', 'bold','fontsize',24);

    lines = findobj(gcf,'Type','Line');
    for i = 1:numel(lines)
        lines(i).LineWidth = 2;
    end