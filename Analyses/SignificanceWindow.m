function [SignificanceWindowTable, NormalizedIndexes] = SignificanceWindow(spiketable,x1,x2,histwidth)
%   Detailed explanation goes here

    % histwidth = 25; %Width of the histogram bins
    % % 
    % % %Bin corresponding to the start of the window of analysis
    % % x1 = x1./histwidth + 83;
    % % x1 = round(x1);
    % % %Bin corresponding to the end of the window of analysis
    % % x2 = x2./histwidth + 108;
    % % x2 = round(x2);
    % %Bin corresponding to the end of the basal activity window    
    % x1b = 100/histwidth + 201;
    % x1b = round(x1b);
    %It stores the spikes corresponding to each of the specified labels in
    %the corresponding variable 

    dev2spikes = vertcat(spiketable.('ODD-ASC-DEV'){:});
    dev1spikes = vertcat(spiketable.('ODD-DES-DEV'){:});
    std1spikes = vertcat(spiketable.('ODD-ASC-STD'){:});
    std2spikes = vertcat(spiketable.('ODD-DES-STD'){:});
    dev2spikesP = vertcat(spiketable.('ODD-ASC-DEV-P'){:});
    dev1spikesP = vertcat(spiketable.('ODD-DES-DEV-P'){:});
    std1spikesP = vertcat(spiketable.('ODD-ASC-STD-P'){:});
    std2spikesP = vertcat(spiketable.('ODD-DES-STD-P'){:});
    ctr1caspikes = vertcat(spiketable.('CASC-ASC-F1'){:});
    ctr1cdspikes = vertcat(spiketable.('CASC-DES-F1'){:});
    ctr2caspikes = vertcat(spiketable.('CASC-ASC-F2'){:});
    ctr2cdspikes = vertcat(spiketable.('CASC-DES-F2'){:});
    ctr1MSspikes = vertcat(spiketable.('MS-F1'){:});
    ctr2MSspikes = vertcat(spiketable.('MS-F2'){:});
%% 
    %It generates the bin center values to use the fuction hist
    % histedges = 0:histwidth:1000;
    % histcenters = histedges(1:end-1)+histwidth/2;

    % DEVa = vertcat(oddData.('ODD-ASC-DEV'){:});
    % DEVb = vertcat(oddData.('ODD-DES-DEV'){:});
    % spikes = [DEVa; DEVb];

    histwidth = 25;
    histedges = -2000:histwidth:2000;
    histcenters = histedges(1:end-1)+histwidth/2;
    histcentersinterp = histcenters; % get new xaxis for interpolated data
    
    % dev2spikes = histcounts(dev2spikes,histedges);
    % dev2spikes = dev2spikes/((histwidth/1000)*80);

    tiledlayout(1,14); % Creates a tiled chart layout for displaying multiple plots
    nexttile

    spikeCountstd1 = hist(std1spikes,histcenters);
    % normalize and convert to density in spk/s
    histdatstd1 = 1000 * spikeCountstd1 / (histwidth*(size(std1spikes,1)));
    %Graph
    bar(histcenters,histdatstd1,'k','BarWidth',1);

%     histdatstd1 = histcounts(std1spikes,histedges);
%     histdatstd1 = (histdatstd1/40)/(histwidth/40); % normalise to Hz
%     bar(histcenters,histdatstd1,'k','BarWidth',1);

    hold on
    xticks(histedges(1:5:end));
    title('STD1');
    temp = get(gca,'ylim');
    ydat(1) = temp(2);

    nexttile
    
    spikeCountstd2 = hist(std2spikes,histcenters);
    % normalize and convert to density in spk/s
    histdatstd2 = 1000 * spikeCountstd2 / (histwidth*(size(std2spikes,1)));

    bar(histcenters,histdatstd2,'k','BarWidth',1);

%     histdatstd2 = histcounts(std2spikes,histedges);
%     histdatstd2 = (histdatstd2/40)/(histwidth/40); % normalise to Hz
%     bar(histcenters,histdatstd2,'k','BarWidth',1);

    hold on
    xticks(histedges(1:5:end));
    title('STD2');
    temp = get(gca,'ylim');
    ydat(2) = temp(2);

    nexttile

    spikeCountctr1 = hist(ctr1caspikes,histcenters);
    % normalize and convert to density in spk/s
    histdatctr1ca = 1000 * spikeCountctr1 / (histwidth*(size(ctr1caspikes,1)));

    bar(histcenters,histdatctr1ca,'k','BarWidth',1);
    
%     histdatctr1ca = histcounts(ctr1caspikes,histedges);
%     histdatctr1ca = (histdatctr1ca/40)/(histwidth/40); % normalise to Hz
%     bar(histcenters,histdatctr1ca,'k','BarWidth',1);

    hold on
    xticks(histedges(1:5:end));
    title('CTR1');
    temp = get(gca,'ylim');
    ydat(3) = temp(2);

    nexttile

    spikeCountctr2 = hist(ctr2caspikes,histcenters);
    % normalize and convert to density in spk/s
    histdatctr2ca = 1000 * spikeCountctr2 / (histwidth*(size(ctr2caspikes,1)));

    bar(histcenters,histdatctr2ca,'k','BarWidth',1);    

%     histdatctr2ca = histcounts(ctr2caspikes,histedges);
%     histdatctr2ca = (histdatctr2ca/40)/(histwidth/40); % normalise to Hz
%     bar(histcenters,histdatctr2ca,'k','BarWidth',1);

    hold on
    xticks(histedges(1:5:end));
    title('CTR2');
    temp = get(gca,'ylim');
    ydat(4) = temp(2);

     nexttile

    spikeCountctr1cd = hist(ctr1cdspikes,histcenters);
    % normalize and convert to density in spk/s
    histdatctr1cd = 1000 * spikeCountctr1cd / (histwidth*(size(ctr1cdspikes,1)));

    bar(histcenters,histdatctr1cd,'k','BarWidth',1);

%     histdatctr1cd = histcounts(ctr1cdspikes,histedges);
%     histdatctr1cd = (histdatctr1cd/40)/(histwidth/40); % normalise to Hz
%     bar(histcenters,histdatctr1cd,'k','BarWidth',1);

    hold on
    xticks(histedges(1:5:end));
    title('CTR1');
    temp = get(gca,'ylim');
    ydat(5) = temp(2);

    nexttile

    spikeCountctr2cd = hist(ctr2cdspikes,histcenters);
    % normalize and convert to density in spk/s
    histdatctr2cd = 1000 * spikeCountctr2cd / (histwidth*(size(ctr2cdspikes,1)));

    bar(histcenters,histdatctr2cd,'k','BarWidth',1);    

%     histdatctr2cd = histcounts(ctr2cdspikes,histedges);
%     histdatctr2cd = (histdatctr2cd/40)/(histwidth/40); % normalise to Hz
%     bar(histcenters,histdatctr2cd,'k','BarWidth',1);

    hold on
    xticks(histedges(1:5:end));
    title('CTR2');
    temp = get(gca,'ylim');
    ydat(6) = temp(2);

    nexttile

    spikeCountctr1MS = hist(ctr1MSspikes,histcenters);
    % normalize and convert to density in spk/s
    histdatctr1MS = 1000 * spikeCountctr1MS / (histwidth*(size(ctr1MSspikes,1)));

    bar(histcenters,histdatctr1MS,'k','BarWidth',1);    

%     histdatctr2cd = histcounts(ctr2cdspikes,histedges);
%     histdatctr2cd = (histdatctr2cd/40)/(histwidth/40); % normalise to Hz
%     bar(histcenters,histdatctr2cd,'k','BarWidth',1);

    hold on
    xticks(histedges(1:5:end));
    title('CTR1 MS');
    temp = get(gca,'ylim');
    ydat(7) = temp(2);
        
    nexttile

    spikeCountctr2MS = hist(ctr2MSspikes,histcenters);
    % normalize and convert to density in spk/s
    histdatctr2MS = 1000 * spikeCountctr2MS / (histwidth*(size(ctr2MSspikes,1)));

    bar(histcenters,histdatctr2MS,'k','BarWidth',1);    

%     histdatctr2cd = histcounts(ctr2cdspikes,histedges);
%     histdatctr2cd = (histdatctr2cd/40)/(histwidth/40); % normalise to Hz
%     bar(histcenters,histdatctr2cd,'k','BarWidth',1);

    hold on
    xticks(histedges(1:5:end));
    title('CTR2 MS');
    temp = get(gca,'ylim');
    ydat(8) = temp(2);


    nexttile

    spikeCountdev1 = hist(dev1spikes,histcenters);
    % normalize and convert to density in spk/s
    histdatdev1 = 1000 * spikeCountdev1 / (histwidth*(size(dev1spikes,1)));

    bar(histcenters,histdatdev1,'k','BarWidth',1);

%     histdatdev1 = histcounts(dev1spikes,histedges);
%     histdatdev1 = (histdatdev1/40)/(histwidth/40); % normalise to Hz
%     bar(histcenters,histdatdev1,'k','BarWidth',1);

    hold on
    xticks(histedges(1:5:end));
    title('DEV1');
    temp = get(gca,'ylim');
    ydat(9) = temp(2);

    nexttile

    spikeCountdev2 = hist(dev2spikes,histcenters);
    % normalize and convert to density in spk/s
    histdatdev2 = 1000 * spikeCountdev2 / (histwidth*(size(dev2spikes,1)));

    bar(histcenters,histdatdev2,'k','BarWidth',1);
    
    hold on
    xticks(histedges(1:5:end));
    title('DEV2');
    temp = get(gca,'ylim');
    ydat(10) = temp(2);

    nexttile

    spikeCountstd1P = hist(std1spikesP,histcenters);
    % normalize and convert to density in spk/s
    histdatstd1P = 1000 * spikeCountstd1P / (histwidth*(size(std1spikesP,1)));
    %Graph
    bar(histcenters,histdatstd1P,'k','BarWidth',1);

%     histdatstd1 = histcounts(std1spikes,histedges);
%     histdatstd1 = (histdatstd1/40)/(histwidth/40); % normalise to Hz
%     bar(histcenters,histdatstd1,'k','BarWidth',1);

    hold on
    xticks(histedges(1:5:end));
    title('STD1 P');
    temp = get(gca,'ylim');
    ydat(11) = temp(2);

    nexttile
    
    spikeCountstd2P = hist(std2spikesP,histcenters);
    % normalize and convert to density in spk/s
    histdatstd2P = 1000 * spikeCountstd2P / (histwidth*(size(std2spikesP,1)));

    bar(histcenters,histdatstd2P,'k','BarWidth',1);

%     histdatstd2 = histcounts(std2spikes,histedges);
%     histdatstd2 = (histdatstd2/40)/(histwidth/40); % normalise to Hz
%     bar(histcenters,histdatstd2,'k','BarWidth',1);

    hold on
    xticks(histedges(1:5:end));
    title('STD2 P');
    temp = get(gca,'ylim');
    ydat(12) = temp(2);

    nexttile
    
    spikeCountdev1P = hist(dev1spikesP,histcenters);
    % normalize and convert to density in spk/s
    histdatdev1P = 1000 * spikeCountdev1P / (histwidth*(size(dev1spikesP,1)));

    bar(histcenters,histdatdev1P,'k','BarWidth',1);

%     histdatdev2 = histcounts(dev2spikes,histedges);
%     histdatdev2 = (histdatdev2/40)/(histwidth/40); % normalise to Hz
%     bar(histcenters,histdatdev2,'k','BarWidth',1);

    hold on
    xticks(histedges(1:5:end));
    title('DEV1 P');
    temp = get(gca,'ylim');
    ydat(13) = temp(2);

    nexttile

    spikeCountdev2P = hist(dev2spikesP,histcenters);
    % normalize and convert to density in spk/s
    histdatdev2P = 1000 * spikeCountdev2P / (histwidth*(size(dev2spikesP,1)));

    bar(histcenters,histdatdev2P,'k','BarWidth',1);

    hold on
    xticks(histedges(1:5:end));
    title('DEV2P');
    temp = get(gca,'ylim');
    ydat(14) = temp(2);

    nexttile(1)
    ylim([0 max(ydat)])
    plot([0 75],[max(ydat) max(ydat)]*0.8,'r')
    nexttile(2)
    ylim([0 max(ydat)])
    plot([0 75],[max(ydat) max(ydat)]*0.8,'r')
    nexttile(3)
    ylim([0 max(ydat)])
    plot([0 75],[max(ydat) max(ydat)]*0.8,'r')
    nexttile(4)
    ylim([0 max(ydat)])
    plot([0 75],[max(ydat) max(ydat)]*0.8,'r')
    nexttile(5)
    ylim([0 max(ydat)])
    plot([0 75],[max(ydat) max(ydat)]*0.8,'r')
    nexttile(6)
    ylim([0 max(ydat)])
    plot([0 75],[max(ydat) max(ydat)]*0.8,'r')
    nexttile(7)
    ylim([0 max(ydat)])
    plot([0 75],[max(ydat) max(ydat)]*0.8,'r')
    nexttile(8)
    ylim([0 max(ydat)])
    plot([0 75],[max(ydat) max(ydat)]*0.8,'r')
    nexttile(9)
    ylim([0 max(ydat)])
    plot([0 75],[max(ydat) max(ydat)]*0.8,'r')
    nexttile(10)
    ylim([0 max(ydat)])
    plot([0 75],[max(ydat) max(ydat)]*0.8,'r')
    nexttile(11)
    ylim([0 max(ydat)])
    plot([0 75],[max(ydat) max(ydat)]*0.8,'r')
    nexttile(12)
    ylim([0 max(ydat)])
    plot([0 75],[max(ydat) max(ydat)]*0.8,'r')
    nexttile(13)
    ylim([0 max(ydat)])
    plot([0 75],[max(ydat) max(ydat)]*0.8,'r')
    nexttile(14)
    ylim([0 max(ydat)])
    plot([0 75],[max(ydat) max(ydat)]*0.8,'r')
    beautify
%% 
x = 85;
x3 = 89;
figure('visible','off')
% %BasalActivity reduction
    nd1 = size(histdatdev1,2);
    BasalActivityMeanyDev1 = mean(histdatdev1(x:x3),"omitnan");
    BasalActivitySTDDev1 = std(histdatdev1(x:x3),"omitnan");

    ns1 = size(histdatstd1,2);
    BasalActivityMeanySTD1 = mean(histdatstd1(x:x3),"omitnan");
    BasalActivitySTDSTD1 = std(histdatstd1(x:x3),"omitnan");

    nd2 = size(histdatdev2,2);
    BasalActivityMeanyDev2 = mean(histdatdev2(x:x3),"omitnan");
    BasalActivitySTDDev2 = std(histdatdev2(x:x3),"omitnan");

    ns2 = size(histdatstd2,2);
    BasalActivityMeanySTD2 = mean(histdatstd2(x:x3),"omitnan");
    BasalActivitySTDSTD2 = std(histdatstd2(x:x3),"omitnan");

    nd1P = size(histdatdev1P,2);
    BasalActivityMeanyDev1P = mean(histdatdev1P(x:x3),"omitnan");
    BasalActivitySTDDev1P = std(histdatdev1P(x:x3),"omitnan");

    ns1P = size(histdatstd1P,2);
    BasalActivityMeanySTD1P = mean(histdatstd1P(x:x3),"omitnan");
    BasalActivitySTDSTD1P = std(histdatstd1P(x:x3),"omitnan");

    nd2P = size(histdatdev2P,2);
    BasalActivityMeanyDev2P = mean(histdatdev2P(x:x3),"omitnan");
    BasalActivitySTDDev2P = std(histdatdev2P(x:x3),"omitnan");

    ns2P = size(histdatstd2P,2);
    BasalActivityMeanySTD2P = mean(histdatstd2P(x:x3),"omitnan");
    BasalActivitySTDSTD2P = std(histdatstd2P(x:x3),"omitnan");

    nca1 = size(histdatctr1ca,2);
    BasalActivityMeanyctr1cas = mean(histdatctr1ca(x:x3),"omitnan");
    BasalActivitySTDctr1cas = std(histdatctr1ca(x:x3),"omitnan");

    ncd1 = size(histdatctr1cd,2);
    BasalActivityMeanyctr1cd = mean(histdatctr1cd(x:x3),"omitnan");
    BasalActivitySTDctr1cd = std(histdatctr1cd(x:x3),"omitnan");

    nca2 = size(histdatctr2ca,2);
    BasalActivityMeanyctr2cas = mean(histdatctr2ca(x:x3),"omitnan");
    BasalActivitySTDctr2cas = std(histdatctr2ca(x:x3),"omitnan");

    ncd2 = size(histdatctr2cd,2);
    BasalActivityMeanyctr2cd = mean(histdatctr2cd(x:x3),"omitnan");
    BasalActivitySTDctr2cd = std(histdatctr2cd(x:x3),"omitnan");

    nMS1 = size(histdatctr1MS,2);
    BasalActivityMeanyctr1MS = mean(histdatctr1MS(x:x3),"omitnan");
    BasalActivitySTDctr1MS = std(histdatctr1MS(x:x3),"omitnan");

    nMS2 = size(histdatctr2MS,2);
    BasalActivityMeanyctr2MS = mean(histdatctr2MS(x:x3-1),"omitnan");
    BasalActivitySTDctr2MS = std(histdatctr2MS(x:x3-1),"omitnan");
% 
%     for i = 1:nd1
%         if histdatdev1(i) < BasalActivityMeanyDev1
%            histdatdev1(i) = 0;
%         else
%            histdatdev1(i) =  histdatdev1(i);
%         end
%     end
% 
%     for i = 1:nd2
%         if histdatdev2(i) < BasalActivityMeanyDev2
%            histdatdev2(i) = 0;
%         else
%            histdatdev2(i) = histdatdev2(i);
%         end
%     end
% 
%     for i = 1:ns1
%         if histdatstd1(i) < BasalActivityMeanySTD1
%            histdatstd1(i) = 0;
%         else
%            histdatstd1(i) = histdatstd1(i);
%         end
%     end
% 
%     for i = 1:ns2
%         if histdatstd2(i) < BasalActivityMeanySTD2
%            histdatstd2(i) = 0;
%         else
%             histdatstd2(i) = histdatstd2(i);
%         end
%     end
% 
%     for i = 1:nca1
%         if histdatctr1ca(i) < BasalActivityctr1cas
%            histdatctr1ca(i) = 0;
%         else
%            histdatctr1ca(i) =  histdatctr1ca(i);
%         end
%     end
% 
%     for i = 1:ncd1
%         if histdatctr1cd(i) < BasalActivityMeanyctr1cd
%            histdatctr1cd(i) = 0;
%         else
%            histdatctr1cd(i) = histdatctr1cd(i);
%         end
%     end
% 
%     for i = 1:nca2
%         if histdatctr2ca(i) < BasalActivityMeanyctr2cas
%            histdatctr2ca(i) = 0;
%         else
%            histdatctr2ca(i) = histdatctr2ca(i);
%         end
%     end
% 
%     for i = 1:ncd2
%         if histdatctr2cd(i) < BasalActivityMeanyctr2cd
%            histdatctr2cd(i) = 0;
%         else
%             histdatctr2cd(i) = histdatctr2cd(i);
%         end
%     end
    
    % if max(dev1spikes) > max(dev2spikes)
    %     ymax = max(dev1spikes) + 0.3;
    % else
    %     ymax = max(dev2spikes) + 0.3;
    % end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
% Average spikes in the window of analysis 
DEV1 = mean(histdatdev1(x1:x2),"omitnan");
STD1 = mean(histdatstd1(x1:x2),"omitnan");

DEV2 = mean(histdatdev2(x1:x2),"omitnan");
STD2 = mean(histdatstd2(x1:x2),"omitnan");

DEV1P = mean(histdatdev1P(x1:x2),"omitnan");
STD1P = mean(histdatstd1P(x1:x2),"omitnan");

DEV2P = mean(histdatdev2P(x1:x2),"omitnan");
STD2P = mean(histdatstd2P(x1:x2),"omitnan");

CTRCA1 = mean(histdatctr1ca(x1:x2),"omitnan");
CTRCA2 = mean(histdatctr2ca(x1:x2),"omitnan");

CTRCD1 = mean(histdatctr1cd(x1:x2),"omitnan");
CTRCD2 = mean(histdatctr2cd(x1:x2),"omitnan");

CTR1MS = mean(histdatctr1MS(x1:x2),"omitnan");
CTR2MS = mean(histdatctr2MS(x1:x2),"omitnan");
%% 
s=2;
if DEV1>(s*BasalActivitySTDDev1)
    SignificanceDEV1 = 1;
else
    SignificanceDEV1 = 0;
end

if DEV2>(s*BasalActivitySTDDev2)
    SignificanceDEV2 = 1;
else
    SignificanceDEV2 = 0;
end
if DEV1P>(s*BasalActivitySTDDev1P)
    SignificanceDEV1P = 1;
else
    SignificanceDEV1P = 0;
end

if DEV2P>(s*BasalActivitySTDDev2P)
    SignificanceDEV2P = 1;
else
    SignificanceDEV2P = 0;
end

if STD1>(s*BasalActivitySTDSTD1)
    SignificanceSTD1 = 1;
else
    SignificanceSTD1 = 0;
end

if STD2>(s*BasalActivitySTDSTD2)
    SignificanceSTD2 = 1;
else
    SignificanceSTD2 = 0;
end
if STD1P>(s*BasalActivitySTDSTD1P)
    SignificanceSTD1P = 1;
else
    SignificanceSTD1P = 0;
end

if STD2P>(s*BasalActivitySTDSTD2P)
    SignificanceSTD2P = 1;
else
    SignificanceSTD2P = 0;
end

if DEV1>(s*BasalActivitySTDDev1) && STD1>(s*BasalActivitySTDSTD1)
    SignificanceSI1 = 1;
else
    SignificanceSI1 = 0;
end

if DEV2>(s*BasalActivitySTDDev2) && STD2>(s*BasalActivitySTDSTD2)
    SignificanceSI2 = 1;
else
    SignificanceSI2 = 0;
end
if DEV1P>(s*BasalActivitySTDDev1P) && STD1P>(s*BasalActivitySTDSTD1P)
    SignificanceSI1P = 1;
else
    SignificanceSI1P = 0;
end

if DEV2P>(s*BasalActivitySTDDev2P) && STD2P>(s*BasalActivitySTDSTD2P)
    SignificanceSI2P = 1;
else
    SignificanceSI2P = 0;
end
    SignificanceWindowTable = table(SignificanceDEV1, SignificanceDEV2, SignificanceDEV1P, ...
        SignificanceDEV2P,SignificanceSTD1, SignificanceSTD2, SignificanceSTD1P, ...
        SignificanceSTD2P,DEV1, BasalActivityMeanyDev1, STD1, ...
        BasalActivityMeanySTD1, DEV2, BasalActivityMeanyDev2,   STD2, ...
        BasalActivityMeanySTD2, DEV1P, BasalActivityMeanyDev1P, STD1P, ...
        BasalActivityMeanySTD1P, DEV2P,   BasalActivityMeanyDev2P, STD2P, ...
        BasalActivityMeanySTD2P, CTRCA1, BasalActivityMeanyctr1cas, CTRCA2, ...
        BasalActivityMeanyctr2cas, CTRCD1, BasalActivityMeanyctr1cd, CTRCD2, ...
        BasalActivityMeanyctr2cd, CTR1MS,BasalActivityMeanyctr1MS, CTR2MS, ...
        BasalActivityMeanyctr2MS);
       %% Index Generator
    %%%%%%%%%%%%%%%%%%%%%% SI & CSI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SI1 = (DEV1-STD1)/(DEV1+STD1);
SI2 = (DEV2-STD2)/(DEV2+STD2);

CSI = ((DEV1+DEV2)-(STD1+STD2))/(DEV1+DEV2+STD1+STD2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalization
Nms1 = sqrt( DEV1^2 + STD1^2 + CTR1MS^2);
Nms2 = sqrt( DEV2^2 + STD2^2 + CTR1MS^2);

Nms1 = sqrt( DEV1^2 + STD1^2 + CTR1MS^2);
Nms2 = sqrt( DEV2^2 + STD2^2 + CTR2MS^2);

NormDEVms1 = DEV1/Nms1;
NormSTDms1 = STD1/Nms1;

NormCTRms1 = CTR1MS/Nms1;
NormCTRms1 = CTR1MS/Nms1;

NormDEVms2 = DEV2/Nms2;
NormSTDms2 = STD2/Nms2;

NormCTRms2 = CTR2MS/Nms2;
NormCTRms2 = CTR2MS/Nms2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalization
Nms1P = sqrt( DEV1P^2 + STD1P^2 + CTR1MS^2);
Nms2P = sqrt( DEV2P^2 + STD2P^2 + CTR2MS^2);

NormDEVPms1 = DEV1/Nms1P;
NormSTDPms1 = STD1/Nms1P;

NormCTRPms1 = CTR1MS/Nms1P;
NormCTRPms1 = CTR1MS/Nms1P;

NormDEVPms2 = DEV2/Nms2P;
NormSTDPms2 = STD2/Nms2P;

NormCTRPms2 = CTR2MS/Nms2P;
NormCTRPms2 = CTR2MS/Nms2P;

%%%%%%%%%%%%%%% CASCASC %%%%%%%%%%%%%%%%%

    iMM1ms = NormDEVms1 - NormSTDms1;
    iRS1ms = NormCTRms1 - NormSTDms1;
    iPE1ms = NormDEVms1 - NormCTRms1;


%%%%%%%%%%%%%%% CASCDESC %%%%%%%%%%%%%%%%%

    iMM2ms = NormDEVms2 - NormSTDms2;
    iRS2ms = NormCTRms2 - NormSTDms2;
    iPE2ms = NormDEVms2 - NormCTRms2;
    %%%%%%%%%%%%%%% CASCASC %%%%%%%%%%%%%%%%%

    iMM1msP = NormDEVPms1 - NormSTDPms1;
    iRS1msP = NormCTRPms1 - NormSTDPms1;
    iPE1msP = NormDEVPms1 - NormCTRPms1;


%%%%%%%%%%%%%%% CASCDESC %%%%%%%%%%%%%%%%%

    iMM2msP = NormDEVPms2 - NormSTDPms2;
    iRS2msP = NormCTRPms2 - NormSTDPms2;
    iPE2msP = NormDEVPms2 - NormCTRPms2;


NormalizedIndexes = table(CSI,SI1, iMM1ms, iRS1ms, iPE1ms,SI2, iMM2ms, iRS2ms, iPE2ms, NormDEVms1, NormSTDms1, NormCTRms1, NormDEVms2, NormSTDms2, NormCTRms2,  NormDEVPms1, NormSTDPms1, NormCTRPms1, NormDEVPms2, NormSTDPms2, NormCTRPms2, DEV1, STD1, CTRCA1, DEV2, STD2, CTRCD2, CTR1MS, CTR2MS);
close all
