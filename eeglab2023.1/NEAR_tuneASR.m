% Defining parameters
loc_in = 'C:\Users\velu.kumaravel\Desktop\Data Drive\Code\GIT\cuttingGardens2023_NEAR_Pipeline\eeglab2023.1\sample_data\NEAR_DATA';
fname_in = 's5.set';
k_in_array = [10,15,20];
process_array = repmat({'on'},length(k_in_array),1); % 'on' or 'off'

% Calling the wrapper
[metric, error_log] = evalASRparams(loc_in, fname_in, k_in_array, process_array);

% Plotting the results for interpretation
figure;
bar(k_in_array, metric);
% plot(k_in_array, metric); 
xlabel('ASR Parameter');
ylabel('Signal-to-Noise Ratio');
set(gca,'fontsize',28,'fontweight','bold');
set(gcf,'color','w');