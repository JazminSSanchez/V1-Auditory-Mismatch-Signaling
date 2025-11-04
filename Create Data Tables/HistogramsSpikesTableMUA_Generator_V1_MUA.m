clear all 
close all
BST = load('D:\Jazmin\MultichannelDataTanks\Cortex\spiketableKS4_MUA');
x1 = 200; %Start of window of analysis
x2 = 800; %End of window of analysis
%% 
for i = 1:size(BST.spiketable,1)
    % for k = 1:size(spiketable.oddData{1},1)
        [HistogramsSpikes] = SignificanceWindow_Spikes(BST.spiketable.oddData{i,1},x1,x2);
        
        HistogramsSpikesTableMUA(i,:) = HistogramsSpikes;
        
end


%%%%%%%%%%%%%%%%%%%%%%%%%
savingPath = 'D:\Jazmin\MultichannelDataTanks\Cortex';
%%%%%%%%%%%%%%%%%%%%%%%%%

cd (savingPath);
save('HistogramsSpikesTable_MUA.mat', 'HistogramsSpikesTableMUA');
%% 
