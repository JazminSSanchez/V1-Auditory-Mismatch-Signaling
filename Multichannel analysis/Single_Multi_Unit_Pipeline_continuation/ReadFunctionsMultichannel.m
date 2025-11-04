clear all 
close all
BST = load(['D:\Jazmin\MultichannelDataTanks\Cortex\spiketableKS4']);
% x1 = 200; %Start of window of analysis
% x2 = 800; %End of window of analysis
x1 = 89; %Start of window of analysis (200 ms)
x2 = 112; %End of window of analysis (800 ms)
%% 
for i = 2:size(BST.spiketable,1)
    % for k = 1:size(spiketable.oddData{1},1)
        [SignificanceWindowRow, NormalizedIndexes] = SignificanceWindow(BST.spiketable.oddData{i,1},x1,x2,25);
        
        SignificanceWindowTable(i,:) = SignificanceWindowRow;
        
        % [NormalizedIndexes] = IndexGeneratorDensityMC(BST.spiketable.oddData{i,1},x1,x2);

        % spiketable.oddData almacena todos los tiempos desde -1000 s hasta 1000 s
        % en los que hubo una espiga por bin
        
        NormalizedIndexesTable(i,:) = NormalizedIndexes;
    % end
end

SignificanceWindowArray = table2array(SignificanceWindowTable);
NormalizedIndexesArray = table2array(NormalizedIndexesTable);
%% 

n=1;
for i = 1:size(BST.spiketable,1)
    if sum(isnan(SignificanceWindowArray(i,:)), 2)>0
        RowsToDelete (n,1) = i;
        n = n + 1;
    end
end
%% 
for j = 1:size(BST.spiketable,1)
    if SignificanceWindowTable.SignificanceDEV1(j,1) == 1 && SignificanceWindowTable.SignificanceDEV2(j,1) == 1
        BST.spiketable.Auditory(j) = 1;
    else
        BST.spiketable.Auditory(j) = 0;
    end
end
%Periodic
for j = 1:size(BST.spiketable,1)
    if SignificanceWindowTable.SignificanceDEV1P(j,1) == 1 && SignificanceWindowTable.SignificanceDEV2P(j,1) == 1
        BST.spiketable.AuditoryP(j) = 1;
    else
        BST.spiketable.AuditoryP(j) = 0;
    end
end

%% 
[CSI CSI_P] = CSI_bootstrapping('D:\Jazmin\MultichannelDataTanks\Cortex\spiketableKS4.mat');
%
%% 
for i = 2:size(BST.spiketable,1)
    vals = CSI.Rango{i,1};  % extract [a b] from cell
    if vals(1) >= 0 && vals(2) >= 0 && BST.spiketable.Auditory(i) == 1
        BST.spiketable.Significance(i) = 1;
    % elseif CSI.Rango(i,1)<=0 && CSI.Rango(i,2)<=0
    %     BST.spiketable.Significance(j) = 1;
    else
        BST.spiketable.Significance(i) = 0;
    end
end

for i = 2:size(BST.spiketable,1)
    vals = CSI_P.Rango_P{i,1};  % extract [a b] from cell
    if vals(1) >= 0 && vals(2) >= 0 && BST.spiketable.AuditoryP(i) == 1
        BST.spiketable.SignificanceP(i) = 1;
    % elseif CSI_P.Rango(i,1)<=0 && CSI_P.Rango(i,2)<=0
    %     BST.spiketable.Significance(j) = 1;
    else
        BST.spiketable.SignificanceP(i) = 0;
    end
end
 %% 
% 
% SignificanceWindowArray(sum(isnan(SignificanceWindowArray), 2)>0, :) = [];
% NormalizedIndexesArray(sum(isnan(NormalizedIndexesArray), 2)>0, :) = [];
SignificanceWindowTableNoNan = SignificanceWindowTable;
NormalizedIndexesTableNoNan = NormalizedIndexesTable;
SignificanceWindowTableNoNan(RowsToDelete',:) = [];
NormalizedIndexesTableNoNan(RowsToDelete',:) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%
savingPath = 'D:\Jazmin\MultichannelDataTanks\Cortex\';
%%%%%%%%%%%%%%%%%%%%%%%%%
spiketable = BST.spiketable;
cd (savingPath);
save('spiketable_Aud_Sig.mat', 'spiketable');
save('NormalizedIndexesTable.mat', 'NormalizedIndexesTable');
save('SignificanceWindowTable.mat', 'SignificanceWindowTable');
%% 
