% Create FRAs from sorted bst files
% Adam Hockley 11/11/2023
%


clear
clc

Main_path = 'D:\Jazmin\MultichannelDataTanks\Cortex\';
tanks = {'24_023','24_024','24_030'};
Sorter = 'KS4';

SpikeWindow = [250 650]; % Only take spike from a set time window (ms)

for ta = 1:length(tanks)

    % Get electrode positions
    Positions = allfolders([Main_path  tanks{ta} ]);
 
    for pos = 1:length(Positions)

        tank_path = [Main_path tanks{ta} '\' Positions{pos}];

        % Use this bit of code to get all FRA blocks in the tank
        Blocks = allfolders(tank_path);
        Blocks = Blocks(contains(Blocks,'FRA'));

        for bl = 1:length(Blocks)

            block_path = [tank_path '\' Blocks{bl}];

            load([block_path '\bst_' Sorter '.mat'])

            Units = unique(bst.Spikes.unit);

            for un = 1:length(Units)

                %% FRA plotting
                figure
                tiledlayout(1,2)
                nexttile
                [RateMatrix, Threshold, BF, SpntRte] = BST_FRA3(bst,SpikeWindow,Units(un),5);
                % nexttile
                % 
                % % And spike shape
                % timevector = (1:length(bst.SpikeShapes(1,:))) / 24.4140625;
                % plot(timevector, bst.SpikeShapes(un,:),'k')
                % figtexts('','','Time (ms)')
                % 
                % nspikes = length(find(bst.Spikes.unit == Units(un)));
                % sgtitle(['nspikes = ' num2str(nspikes)])

                beautify
                set(gcf, 'Color', [1,1,1], 'Position', [100 100 1000 400]);

                savpath = [tank_path '\' Blocks{bl} '\FRA_plots'];
                if ~exist(savpath,'dir'), mkdir(savpath); end

                exportgraphics(gcf,[savpath '\' num2str(Units(un)) '_FRA.eps'])% x = imresize(x,0.4811161);
                print(gcf,'-dpng','-r150',[savpath '\' num2str(Units(un)) '_FRA.png'])
                close all

            end
        end
    end
end

% Trial select for BST
function [RateMatrix, Threshold, BF, SpntRte, max_idx, max_num] = BST_FRA3(bst, SpikeWindow,un,interp)

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

        % Find the trials with the desired
        % properties, then get the SW spike
        % timings for each
        SelectedTrials = BST_TS3(bst,'wfrq',Freqs(ii),'us1_',Levs(iii));
        SpikeTimes = 1000 * BST_GS3(bst,SelectedTrials,un);
        % Just spikes from the stimulus, or all
        % spikes in the sweep?
        if ~isempty(SpikeWindow)
            SpikeTimesL = (SpikeWindow(1) <= SpikeTimes) & (SpikeTimes <= SpikeWindow(2));
            SpikeTotal = sum(SpikeTimesL);
        else
            SpikeTotal = length(SpikeTimes);
        end

        %Convert spike total to spike rate, and
        %add to the rate matrix in the correct
        %position
        sweeplength = bst.Epocs.TSOn.bind(11) - bst.Epocs.TSOn.bind(10); % Get the sweep length in seconds
        SpikeRate = (SpikeTotal / length(SelectedTrials)) / sweeplength;
        RateMatrix((length(Levs)-iii+1),ii) = SpikeRate;
    end
end

%% Get threshold and BF
% Finding the boundary of the RF...
toplim = mean(RateMatrix(end,:),'omitnan') + (2*std(RateMatrix(end,:),'omitnan')); %setting the limit for an increase from baseline (0dB)
botlim = mean(RateMatrix(end,:),'omitnan') - (2*std(RateMatrix(end,:),'omitnan')); %setting the limit for an decrease from baseline (0dB)

ThresholdRateTable = toplim;

siggrid = zeros(length(Levs),length(Freqs));
negsiggrid = zeros(length(Levs),length(Freqs));

% Fill the NEGATIVE significance grid
for m = 1:length(negsiggrid(:,1))
    for sg = 1:length(negsiggrid(1,:))
        if RateMatrix(m,sg) < botlim
            negsiggrid(m,sg) = 1;
        end
    end
end
%                 TypesStruct.TotalInhibition(ch,n) = sum(sum(negsiggrid));
%
% Fill the significance grid
for m = 1:length(siggrid(:,1))
    for sg = 1:length(siggrid(1,:))
        if RateMatrix(m,sg) > toplim
            siggrid(m,sg) = 1;
        end
    end
end

% Remove all not connected to others
siggrid = double(bwareaopen(siggrid, 3));
% heatmap(siggrid) % Check the stats match the RF

SpntRte = mean(RateMatrix(end,:));

% Get threshold and BF from the significance grid (Chooses the
% higher frequency point if there are 2 equal max points)
[max_num,max_idx] = max(sum(siggrid));
max_idx = find(sum(siggrid)==max_num,1,'last');
BF = Freqs(max_idx);
if max_num > 0
    Threshold = LevsUD(max_num);
else
    Threshold = 100;
end

%% Plot
%Plot the RF

if interp > 1
    RateMatrixInterp = imgaussfilt(RateMatrix,1); % Gaussian 2d filter
    %     RateMatrixInterp = imresize(RateMatrixInterp, interp); % Increase resolution by a factor of x;
    RateMatrixInterp = imresize(RateMatrixInterp, 'SCALE', [interp*2 interp/2]); % Increase resolution by a factor of x;
    %     RateMatrixInterp = imresize(RateMatrixInterp, [interp*2 interp/2]); % Increase resolution by a factor of x;
    imagesc(RateMatrixInterp)
    %     pc = pcolor(flipud(RateMatrixInterp)); % plot interpolated data
    %     pc.LineStyle = 'none';
else
    imagesc(RateMatrix)

end
%                 insertMarker(I,[1 2])
%                 set(gcf,'Visible', 'on');
hold on
scatter(max_idx*interp/2,max_num*interp*2,100,'go','filled')
colormap(inferno)

xticks(1:(6*(interp/2)):length(Freqs)*(interp/2))
xticklabels(Freqs(1:6:end)/1000)
%                 xtickangle(45)
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
%                 caxis([0 PeakRate]);

end

%%
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