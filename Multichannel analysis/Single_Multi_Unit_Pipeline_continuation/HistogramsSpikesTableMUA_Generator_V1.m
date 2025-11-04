clear all 
close all
BST = load('D:\Jazmin\Multichannel Data Tanks\Cortex\spiketable_KS4');
x1 = 250; %Start of window of analysis
x2 = 650; %End of window of analysis
%% 
for i = 1:size(BST.spiketable,1)
    % for k = 1:size(spiketable.oddData{1},1)
        [HistogramsSpikes] = SignificanceWindow_Spikes(BST.spiketable.oddData{i,1},x1,x2);
        
        HistogramsSpikesTable(i,:) = HistogramsSpikes;
        
end


%%%%%%%%%%%%%%%%%%%%%%%%%
savingPath = 'D:\Jazmin\Multichannel Data Tanks\Cortex';
%%%%%%%%%%%%%%%%%%%%%%%%%

cd (savingPath);
save('HistogramsSpikesTable.mat', 'HistogramsSpikesTable');
%% 
