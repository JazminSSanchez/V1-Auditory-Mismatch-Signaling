clear all
load('D:\Jazmin\MultichannelDataTanks\HIP\24_023\2\1-1-STIM-H\bst_KS4_MUA.mat')

% Parameters
unitOfInterest = 83;
binSize = 0.025; % 25 ms
timeWindow = [0 1]; % from stimulus onset to 1000 ms

% Step 1: Extract spike data for the unit
n = 1;
for i = 1:size(bst.Spikes.RasterSW,1)
    if bst.Spikes.unit(i) == unitOfInterest
        Spikes_Tr_ms(n,1) = bst.Spikes.TrialIdx(i);  % trial
        Spikes_Tr_ms(n,2) = bst.Spikes.RasterSW(i);  % spike time (s)
        n = n + 1;
    end
end

% Step 2: Identify STD and DEV trials
ODDRTrials = 21:420;
STDtrials = [];
DEVtrials = [];
m = 1; k = 1;
for j = ODDRTrials
    if contains(bst.Epocs.Values.type{j}, 'STD')
        STDtrials(m,1) = j; m = m + 1;
    elseif contains(bst.Epocs.Values.type{j}, 'DEV')
        DEVtrials(k,1) = j; k = k + 1;
    end
end

% Step 3: Define PSTH bins
edges = timeWindow(1):binSize:timeWindow(2); % bin edges from 0 to 1 s
binCenters = edges(1:end-1) + binSize/2;     % center of each bin

% Step 4: Collect spike times for each condition
spikeTimes_STD = Spikes_Tr_ms(ismember(Spikes_Tr_ms(:,1), STDtrials), 2);
spikeTimes_DEV = Spikes_Tr_ms(ismember(Spikes_Tr_ms(:,1), DEVtrials), 2);

% Step 5: Histogram counts normalized by trial count and bin width
[counts_STD, ~] = histcounts(spikeTimes_STD, edges);
counts_STD = counts_STD / length(STDtrials) / binSize; % spikes/s

[counts_DEV, ~] = histcounts(spikeTimes_DEV, edges);
counts_DEV = counts_DEV / length(DEVtrials) / binSize; % spikes/s

% Step 6: Plot PSTH
figure; hold on;
plot(binCenters*1000, counts_STD, 'b', 'LineWidth', 2); % convert x to ms
plot(binCenters*1000, counts_DEV, 'r', 'LineWidth', 2);
xlabel('Time (ms)');
ylabel('Firing rate (spikes/s)');
legend('STD', 'DEV');
title(['PSTH - Unit ', num2str(unitOfInterest), ' (0–1000 ms, 25 ms bins)']);
xlim([0 1000]);
set(gca, 'Box', 'off');

% Step 7: Save data and figure
outputDir = 'D:\Jazmin\MultichannelDataTanks\HIP\PSTH';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% Save figure
saveas(gcf, fullfile(outputDir, ['PSTH_Unit' num2str(unitOfInterest) '_0to1s_25ms.pdf']));
saveas(gcf, fullfile(outputDir, ['PSTH_Unit' num2str(unitOfInterest) '_0to1s_25ms.jpg']));

% Save data to CSV (in milliseconds)
PSTH_Data = table(binCenters'*1000, counts_STD', counts_DEV', ...
    'VariableNames', {'Time_ms', 'STD_FR', 'DEV_FR'});
writetable(PSTH_Data, fullfile(outputDir, ['PSTH_Unit' num2str(unitOfInterest) '_0to1s_25ms.csv']));
