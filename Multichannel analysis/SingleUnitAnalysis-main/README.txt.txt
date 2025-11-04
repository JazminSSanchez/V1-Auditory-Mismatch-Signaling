# SingleUnitAnalysis
Personal code, not public, and specific to CANELAB stimuli.

All code here works on data organised in the bst format, produced by the the [TDT_Single_Unit_Pipeline](https://github.com/adam-hockley/TDT_Single_Unit_Pipeline)

As previously, for the looping to work, data folders should be organised as: **AnimalTank/PositionNumber/RecordingBlocks**.

Modified to work for visual cortex data SOA 1000 ms by Jazmin S. Sanchez

Overview of code contained here:

## 1) FRA_Plot.m & FRA_Plot_MUA.m (MUA = multiunits)
Plots the FRA and spike shape for each neuron and save the figure.

## 2) AddTypesJazmin.m & AddTypesJazmin_MUA.m
Adds the types (e.g. ODD-DEV-ASC) to the bst stimuli table. Allows easy searching of spikes based on stimulus type in later analysis.
## Should be checked for accurate type labeling by all users!

## 3) CreateSpikeTable_Jazmin_Cortex.m & CreateSpikeTable_Jazmin_MUA_Cortex.m
Combines data from each bst and saves into one big table of spike data correspaonding to different stimuli used. Also produces PSTH plots for each neuron.

*** change line 111 depeding on your folder name to save the block name in the spiketable.

spiketableKS4.mat
spiketable_KS4.mat
spiketableKS4_MUA.mat

The restulting data table should contain all the data you need, and be small enough to quickly load, analyse and plot any way you want.
