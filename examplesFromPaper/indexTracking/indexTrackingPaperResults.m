
% clear variables
clear; clc;

% setup paths
run("../../starMOptSetup.m")


%% setup and run experiment

options                  = indexTrackingExperimentParameters();
options.initDateTrain    = '1-Jan-2023';
options.endDateTrain     = '31-May-2023';
options.initDateTest     = '1-Jun-2023';
options.endDateTest      = '31-Jul-2023';
options.maxIter          = 100;
options.alpha            = 1e-2;

% update filename
options.setFilename();

% get struct
options = obj2struct(options);


%%
% main 
saveFlag = 1;
indexTrackingRun(options, saveFlag);

% return;
%% 

% choose date of experiments
myDate = '10-Nov-2023';
% myDate = date;

% choose main directory
dirName = 'indexTrackingResults/';

% be sure to plot without docking
close all;
set(0,'DefaultFigureWindowStyle','normal')

% load results
load(sprintf([dirName,options.filename,'/',options.filename]))

%% 

[~,idx] = sort(results.sectors.names);


subDirname = [dirName,results.options.filename,'/','weightDistribution/'];

if ~exist(subDirname,'dir'), mkdir(subDirname); end

T = [];
for i = 1:length(results.sectorDistribution)
    T = cat(2,T,results.sectorDistribution{idx(i)});
end
T = array2table(T,'VariableNames',results.Mnames);
writetable(T,[subDirname,'sectorDistribution.csv']);

% save weight per stock
T = [];
for i = 1:length(results.X)
    T = cat(2,T,results.X{i}(:));
end

T = array2table(T,'VariableNames',results.Mnames);
writetable(T,[subDirname,'weightsPerStock.csv']);

% create pie charts
for i = 2:length(results.M)
    fig = figure(1); clf;

    pie(results.sectorDistribution{i}(idx))
    
    
    colormap jet;
    set(gca,'FontSize',18)
    exportgraphics(fig,[subDirname,'pie_',results.Mnames{i},'.png'],'BackgroundColor','none')
end

fig = figure(1); clf;
sectorNames = fields(results.sectors);
pie(1 / 11 * ones(1,11),sectorNames(1:11));

colormap jet;
set(gca,'FontSize',18)
exportgraphics(fig,[subDirname,'pie_legend.png'],'BackgroundColor','none')

fig = indexTrackingPlotWeights(results);
exportgraphics(fig,[subDirname,'weights_',results.Mnames{i},'.png'],'BackgroundColor','none')

%% create tracking plots

subDirname = [dirName,results.options.filename,'/','sectorTracking/'];
if ~exist(subDirname,'dir'), mkdir(subDirname); end

futureDate = '31-Jul-2023';
[A,~,~,sectors] = indexTrackingSetupData(results.options.initDateTrain, results.options.endDateTrain,'first');
[~,B]           = indexTrackingSetupData(results.options.initDateTrain, futureDate,'first');

sectorNames = fields(results.sectors);

for i = 1:11
    fig = indexTrackingPlotDataOneSectorTracking(results, results.options.initDateTrain, results.options.endDateTrain, 0, i, futureDate,A,B);
    exportgraphics(fig,[subDirname,'tracking_',sectorNames{i},'.png'],'BackgroundColor','none')

    T = indexTrackingWriteCSV(results, results.options.initDateTrain, results.options.endDateTrain, 0, i, futureDate,A,B);
    writetable(T,[subDirname,'tracking_',sectorNames{i},'.csv']);
end


%% create backtesting plots

subDirname = [dirName,results.options.filename,'/','backtesting/'];
if ~exist(subDirname,'dir'), mkdir(subDirname); end

fig = indexTrackingPlotBacktesting(results, results.options.initDateTrain, results.options.endDateTrain);
exportgraphics(fig,[subDirname,'train.png'],'BackgroundColor','none')

for d = {'31-Jul-2023','30-Sep-2023'}
    fig = indexTrackingPlotBacktesting(results, results.options.initDateTest, d{1});
    exportgraphics(fig,[subDirname,'test--',results.options.initDateTest,'--',d{1},'.png'],'BackgroundColor','none')
end






