% Create FRAs from sorted bst files
% Adam Hockley 11/11/2023
% Edited by Jazmin Sanchez 2025
clear
clc

Main_path = 'D:\Jazmin\MultichannelDataTanks\Cortex\';
tanks = {'24_023','24_024','24_030'};
Sorter = 'KS4';

SpikeWindow = [200 800]; % Only take spike from a set time window (ms)

for ta = 1:length(tanks)

    % Get electrode positions
    Positions = allfolders([Main_path tanks{ta}]);

    for pos = 1:length(Positions)

        tank_path = [Main_path tanks{ta} '\' Positions{pos}];

        % Use this bit of code to get all FRA blocks in the tank
        Blocks = allfolders(tank_path);
        Blocks = Blocks(contains(Blocks,'FRA'));

        for bl = 1:length(Blocks)

            block_path = [tank_path '\' Blocks{bl}];

            load([block_path '\bst_' Sorter '_MUA.mat'])

            Units = unique(bst.Spikes.unit);

            for un = 1:length(Units)

                %% FRA plotting
                figure

                % --- choose interpolation factor here so figure and Excel match ---
                interp = 5;

                [RateMatrix, Threshold, BF, SpntRte] = BST_FRA3(bst, SpikeWindow, Units(un), interp);

                beautify
                set(gcf, 'Color', [1,1,1], 'Position', [100 100 1000 400]);

                savpath = [tank_path '\' Blocks{bl} '\FRA_plots_MUA'];
                if ~exist(savpath,'dir'), mkdir(savpath); end

                % Save figure
                exportgraphics(gcf,[savpath '\' num2str(Units(un)) '_FRA_MUA.eps'])
                print(gcf,'-dpng','-r150',[savpath '\' num2str(Units(un)) '_FRA_MUA.png'])

                % ===== Excel export for FRA figure (same basename as the saved PNG) =====
                figPngPath = [savpath '\' num2str(Units(un)) '_FRA_MUA.png'];
                [~, figBase, ~] = fileparts(figPngPath);
                xlsPath = fullfile(savpath, [figBase '.xlsx']);

                % Recreate the exact axis vectors used in the plot
                Epoks  = bst.Epocs.Values;
                Freqs  = unique(Epoks.wfrq);        % all frequencies used
                Freqs  = Freqs(Freqs~=100);         % remove silent trials
                Levs   = unique(Epoks.us1_);        % all levels used
                LevsUD = flipud(Levs);              % plotting used flipped levels

                % Recreate the interpolated matrix (if any), to match the plotted image
                RateMatrixInterpOut = [];
                if interp > 1
                    RateMatrixInterpOut = imgaussfilt(RateMatrix, 1);
                    % Same call as inside BST_FRA3:
                    RateMatrixInterpOut = imresize(RateMatrixInterpOut, 'SCALE', [interp*2 interp/2]);
                end

                % Tick positions and labels exactly as used
                xt = 1:(6*(interp/2)):length(Freqs)*(interp/2);
                xtlbl_kHz = Freqs(1:6:end)/1000;

                yt = 2:(2*interp*2):length(Levs)*interp*2;
                ytlbl_dB = LevsUD(2:2:end);

                % Final color limits as used in plotting logic
                if interp > 1
                    maxVal = max(RateMatrixInterpOut, [], "all", "omitnan");
                else
                    maxVal = max(RateMatrix, [], "all", "omitnan");
                end
                if maxVal < 2
                    CLow = 0; CHigh = 2;
                else
                    CLow = 0; CHigh = maxVal;
                end

                % Write everything to Excel (separate sheets)
                % 1) Data matrices
                writematrix(RateMatrix, xlsPath, 'Sheet', 'RateMatrix', 'Range', 'A1');
                if ~isempty(RateMatrixInterpOut)
                    writematrix(RateMatrixInterpOut, xlsPath, 'Sheet', 'RateMatrixInterp', 'Range', 'A1');
                end

                % 2) Axes (Hz + kHz; Levels and flipped Levels)
                writetable(table(Freqs(:), Freqs(:)/1000, 'VariableNames', {'Freq_Hz','Freq_kHz'}), ...
                           xlsPath, 'Sheet', 'Freqs');
                writematrix(Levs(:),   xlsPath, 'Sheet', 'Levels_Raw',    'Range', 'A1');
                writematrix(LevsUD(:), xlsPath, 'Sheet', 'Levels_Flipped','Range', 'A1');

                % 3) Ticks/labels actually used in the plot
                writetable(table(xt(:), xtlbl_kHz(:), 'VariableNames', {'xticks','xticklabels_kHz'}), ...
                           xlsPath, 'Sheet', 'Xticks');
                writetable(table(yt(:), ytlbl_dB(:),  'VariableNames', {'yticks','yticklabels_dB'}), ...
                           xlsPath, 'Sheet', 'Yticks');

                % 4) Useful metadata to fully reproduce the panel
                writetable(table( ...
                    string(datetime('now')), interp, CLow, CHigh, Threshold, BF, SpntRte, ...
                    'VariableNames', {'SavedOn','interp','CLow','CHigh','Threshold_dB','BF_Hz','SpontRate_Hz'}), ...
                    xlsPath, 'Sheet', 'Meta');

                disp("Exported FRA vectors/matrices to: " + xlsPath);
                % =========================================================================

                close all

            end
        end
    end
end

% =======================================================================
% Trial select for BST
function [RateMatrix, Threshold, BF, SpntRte, max_idx, max_num] = BST_FRA3(bst, SpikeWindow, un, interp)

Epoks = bst.Epocs.Values;

% Get the list of all frequencies and levels used
Freqs = unique(Epoks.wfrq);
Levs = unique(Epoks.us1_);
LevsUD = flipud(Levs);

% Remove silent trials (wfrq = 100)
Freqs = Freqs(Freqs~=100);

% Create the empty rate matrix for the RF
clear RateMatrix
RateMatrix = zeros(length(Levs),length(Freqs));

% Find which rows contain the Freq and Lev required
for ii = 1:length(Freqs)
    for iii= 1:length(Levs)

        % Find the trials with the desired properties, then get the SW spike timings
        SelectedTrials = BST_TS3(bst,'wfrq',Freqs(ii),'us1_',Levs(iii));
        SpikeTimes = 1000 * BST_GS3(bst,SelectedTrials,un);

        % Just spikes from the stimulus, or all spikes in the sweep?
        if ~isempty(SpikeWindow)
            SpikeTimesL = (SpikeWindow(1) <= SpikeTimes) & (SpikeTimes <= SpikeWindow(2));
            SpikeTotal = sum(SpikeTimesL);
        else
            SpikeTotal = length(SpikeTimes);
        end

        % Convert spike total to spike rate, and add to the rate matrix in the correct position
        sweeplength = bst.Epocs.TSOn.bind(11) - bst.Epocs.TSOn.bind(10); % sweep length in seconds
        SpikeRate = (SpikeTotal / length(SelectedTrials)) / sweeplength;
        RateMatrix((length(Levs)-iii+1),ii) = SpikeRate;
    end
end

%% Get threshold and BF
% Finding the boundary of the RF...
toplim = mean(RateMatrix(end,:),'omitnan') + (2*std(RateMatrix(end,:),'omitnan')); % increase from baseline (0dB)
botlim = mean(RateMatrix(end,:),'omitnan') - (2*std(RateMatrix(end,:),'omitnan')); % decrease from baseline (0dB)

ThresholdRateTable = toplim; %#ok<NASGU>

siggrid = zeros(length(Levs),length(Freqs));
negsiggrid = zeros(length(Levs),length(Freqs));

% Fill the NEGATIVE significance grid
for m = 1:size(negsiggrid,1)
    for sg = 1:size(negsiggrid,2)
        if RateMatrix(m,sg) < botlim
            negsiggrid(m,sg) = 1;
        end
    end
end

% Fill the significance grid
for m = 1:size(siggrid,1)
    for sg = 1:size(siggrid,2)
        if RateMatrix(m,sg) > toplim
            siggrid(m,sg) = 1;
        end
    end
end

% Remove all not connected to others
siggrid = double(bwareaopen(siggrid, 3));
% heatmap(siggrid) % Check the stats match the RF

SpntRte = mean(RateMatrix(end,:));

% Get threshold and BF from the significance grid (Chooses higher freq if two equal max points)
[max_num,max_idx] = max(sum(siggrid));
max_idx = find(sum(siggrid)==max_num,1,'last');
BF = Freqs(max_idx);
if max_num > 0
    Threshold = LevsUD(max_num);
else
    Threshold = 100;
end

%% Plot
% Plot the RF
if interp > 1
    RateMatrixInterp = imgaussfilt(RateMatrix,1); % Gaussian 2d filter
    RateMatrixInterp = imresize(RateMatrixInterp, 'SCALE', [interp*2 interp/2]); % Increase resolution
    imagesc(RateMatrixInterp)
else
    imagesc(RateMatrix)
end
hold on
colormap(jet)

xticks(1:(6*(interp/2)):length(Freqs)*(interp/2))
xticklabels(Freqs(1:6:end)/1000)
yticks(2:(2*interp*2):length(Levs)*interp*2)
yticklabels(LevsUD(2:2:end))
set(gcf, 'Color', [1,1,1], 'Position', [200 200 700 500]); %
figtexts('','Level (dB SPL)','Frequency (kHz)');
set(findall(gcf,'-property','FontSize'),'FontSize',20)
set(gca,'TickDir','out');
h = colorbar;
h.Label.String = 'Spike rate (Hz)';

% optional, increase color limits of units with very little driving
colorlimits = clim;
if colorlimits(2) < 2
    set(gca,'CLim',[0 2]); % optional, set colorbars to same axis
else
    set(gca,'CLim',[0 colorlimits(2)]); % optional, set colorbars to same axis
end
% caxis([0 PeakRate]);

end

%% =======================================================================
function folders = allfolders(directory)
folders = dir(directory);
dirFlags = [folders.isdir] & ~strcmp({folders.name},'.') & ~strcmp({folders.name},'..');
folders = folders(dirFlags);
folders = {folders.name};
end

%%
function figtexts(tit,ylab,xlab)
% Just easier to have one line to set title, y and x
% Leave empty string for no text e.g. 
% figtexts('','MYYAXIS','')
title(gca,tit)
ylabel(gca,ylab)
xlabel(gca,xlab)
% box off
end
