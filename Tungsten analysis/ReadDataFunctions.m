clear all

Sheet = 1; %Unit table with the information of interest 

[pathsSpikes tmp] = ReadPCAPaths('D:\PRJ_HIPPOCAMPUS_JSS\PCASorting\Cortex\spikes'); %Reads all the folders saved in the specified path

[UnitTable] = ReadUnitTable("A2",Sheet,'D:\PRJ_HIPPOCAMPUS_JSS\PCASorting\Cortex\UnitTable_V1_PCA_Junio2025.xlsx');
%Create a table named UnitTable with all the information of the unit table

[pathsSpikesPCA] = pathsPCA(tmp,pathsSpikes,UnitTable);
%Creates a table with all the directories of interest

[BST, SPK, SPK9] = ExtractMatFiles(UnitTable,pathsSpikesPCA);
%Read and save in a structure the complete spikedata, spikedata9 (when the deviant was presented after 9 standards) and the hole BST
x1 = 9; %200 Start of window of analysis
x2 = 32; %800 End of window of analysis
x1b = 37;
x2b = 40; 
%900 - 1000 ms
%%

for i = 1:size(UnitTable,1)
[NormalizedIndexes] = IndexGeneratorDensity(SPK(i).spikedata,x1,x2,x1b,x2b);
NormalizedIndexesTable(i,:) = NormalizedIndexes;
end

%%%%%%%%%%%%%%%%%%%%%%%%%
savingPath = 'D:\PRJ_HIPPOCAMPUS_JSS\PCASorting\Cortex';
%%%%%%%%%%%%%%%%%%%%%%%%%

cd ('D:\PRJ_HIPPOCAMPUS_JSS\PCASorting\Cortex');
save('NormalizedIndexesTableDG_Sig.mat', 'NormalizedIndexesTable');

%%
for i = 1:size(UnitTable,1)
[NormalizedIndexesC] = IndexGeneratorDensity(SPK(i).spikedata,x1,x2,x1b,x2b);
NormalizedIndexesTableC(i,:) = NormalizedIndexesC;
end

cd (savingPath);
save('NormalizedIndexesTableC.mat', ['NormalizedIndexesTableC']);

%%
for i = 1:size(UnitTable,1)
[NormalizedIndexes9] = IndexGeneratorDensity(SPK9(i).spikedata9,x1,x2,x1b,x2b);
NormalizedIndexesTable9(i,:) = NormalizedIndexes9;
end

cd (savingPath);
save('NormalizedIndexesTable9.mat', ['NormalizedIndexesTable9']);

%%
for i = 1:size(UnitTable,1)
[NormalizedIndexes9C] = IndexGeneratorDensity(SPK9(i).spikedata9,x1,x2,x1b,x2b);
NormalizedIndexesTable9C(i,:) = NormalizedIndexes9C;
end

cd (savingPath);
save('NormalizedIndexesTable9C.mat', ['NormalizedIndexesTable9C']);