clear all 
close all
BST = load('D:\Jazmin\MultichannelDataTanks\Cortex\spiketableKS4_MUA_July');
% x1 = 200; %Start of window of analysis
% x2 = 800; %End of window of analysis
x1 = 89; %Start of window of analysis (200 ms)
x2 = 112; %End of window of analysis (800 ms)  
%% 
for i = 1:size(BST.spiketable,1)
    % for k = 1:size(spiketable.oddData{1},1)
        [SignificanceWindowRow, NormalizedIndexes] = SignificanceWindow(BST.spiketable.oddData{i,1},x1,x2,25);
        
        SignificanceWindowTableMUA(i,:) = SignificanceWindowRow;
        
        % [NormalizedIndexes] = IndexGeneratorDensityMC(BST.spiketable.oddData{i,1},x1,x2);

        % spiketable.oddData almacena todos los tiempos desde -1000 s hasta 1000 s
        % en los que hubo una espiga por bin
        
        NormalizedIndexesTableMUA(i,:) = NormalizedIndexes;
    % end
end

SignificanceWindowArrayMUA = table2array(SignificanceWindowTableMUA);
NormalizedIndexesArrayMUA = table2array(NormalizedIndexesTableMUA);
%% 

n=1;
for i = 1:size(BST.spiketable,1)
if sum(isnan(SignificanceWindowArrayMUA(i,:)), 2) > 0
    RowsToDeleteMUA (n,1) = i;
    n = n + 1;
end
end
%% 
for j = 1:size(BST.spiketable,1)
    if SignificanceWindowTableMUA.SignificanceDEV1(j,1) == 1 && SignificanceWindowTableMUA.SignificanceDEV2(j,1) == 1
        BST.spiketable.Auditory(j) = 1;
    else
        BST.spiketable.Auditory(j) = 0;
    end
end
%Periodic
for j = 1:size(BST.spiketable,1)
    if SignificanceWindowTableMUA.SignificanceDEV1P(j,1) == 1 && SignificanceWindowTableMUA.SignificanceDEV2P(j,1) == 1
        BST.spiketable.AuditoryP(j) = 1;
    else
        BST.spiketable.AuditoryP(j) = 0;
    end
end

%% 
[CSI CSI_P] = CSI_bootstrapping_MUA('D:\Jazmin\MultichannelDataTanks\Cortex\spiketableKS4_MUA');
%
%% 
for i = 1:size(BST.spiketable,1)
    vals = CSI{i,4};  % extract [a b] from cell
    if vals(1) >= 0 && vals(2) >= 0 && BST.spiketable.Auditory(i) == 1
        BST.spiketable.Significance(i) = 1;
    % elseif CSI.Rango(i,1)<=0 && CSI.Rango(i,2)<=0
    %     BST.spiketable.Significance(j) = 1;
    else
        BST.spiketable.Significance(i) = 0;
    end
end


for j = 1:size(BST.spiketable,1)
     vals = CSI_P{j,4};  % extract [a b] from cell
    if vals(1) >= 0 && vals(2) >= 0 && BST.spiketable.AuditoryP(j) == 1
        BST.spiketable.SignificanceP(j) = 1;
    % elseif CSI_P.Rango(i,1)<=0 && CSI_P.Rango(i,2)<=0
    %     BST.spiketable.Significance(j) = 1;
    else
        BST.spiketable.SignificanceP(j) = 0;
    end
end
%% 
% 
% SignificanceWindowArray(sum(isnan(SignificanceWindowArray), 2)>0, :) = [];
% NormalizedIndexesArray(sum(isnan(NormalizedIndexesArray), 2)>0, :) = [];

SignificanceWindowTableNoNanMUA = SignificanceWindowTableMUA;
NormalizedIndexesTableNoNanMUA = NormalizedIndexesTableMUA;
SignificanceWindowTableNoNanMUA(RowsToDeleteMUA',:) = [];
NormalizedIndexesTableNoNanMUA(RowsToDeleteMUA',:) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%
savingPath = 'D:\Jazmin\MultichannelDataTanks\Cortex\';
%%%%%%%%%%%%%%%%%%%%%%%%%
spiketable = BST.spiketable;
cd (savingPath);
save('spiketable_Aud_Sig_MUA_July.mat', 'spiketable');
save('NormalizedIndexesTable_MUA_July.mat', 'NormalizedIndexesTableMUA');
save('SignificanceWindowTable_MUA_July.mat', 'SignificanceWindowTableMUA');
%% 
