% Defining parameters
loc_in = './TuneASR';
fname_in = 's5_filt_segt_badch.set'; % dataset with bad channels removed
k_in_array = [5,15,25]; % list of ASR parameters
process_array = repmat({'on'},length(k_in_array),1); % 'on' for Rejection or 'off' for Correction

% Calling the wrapper
[metric, error_log] = evalASRparams(loc_in, fname_in, k_in_array, process_array);

% Plotting the results for interpretation
figure;
% bar(k_in_array, metric);
plot(k_in_array, metric); 
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
lines(i).LineWidth = 2;
end
xlabel('ASR Parameter');
ylabel('Signal-to-Noise Ratio');
set(gca,'fontsize',28,'fontweight','bold');
set(gcf,'color','w');