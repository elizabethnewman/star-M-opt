

% clear variables
clear; clc;

% setup paths
run("../../starMOptSetup.m")


%%

% choose options


defaultOptions = romExperimentParameters();

options = {};
count = 1;

for MInit = {'C','I','Z'}
    
    tmpOptions              = defaultOptions;
    tmpOptions.k            = 2;
    tmpOptions.cN           = 50;
    tmpOptions.maxIter      = 1000;
    tmpOptions.MInitialize  = MInit{1};

    % update filename
    tmpOptions.setFilename();

    % get struct
    options{count} = obj2struct(tmpOptions);
    options{count}.filename = [date,'/',options{count}.filename];

    % update counter
    count = count + 1;
end

%%
% main loop
saveFlag = 1;
for i = 1:length(options)
    results = romRun(options{i}, saveFlag);
end

%% 

% choose date of experiments
myDate = date;

% choose main directory
dirName = 'romResults/';

% be sure to plot without docking
close all;
set(0,'DefaultFigureWindowStyle','normal')


%% create csv files
for MInit = {'C','I','Z'}
    load(sprintf([dirName,myDate,'/',myDate,'--k-02--init-%s'],MInit{1}));
    
    % convergence
    T = array2table(results.optInfo.values, 'VariableNames', results.optInfo.header);
    writetable(T,sprintf([dirName,myDate,'/convergence--init-%s.csv'],MInit{1}));

    
    % relative error per parameter
    c = linspace(results.options.cMin,results.options.cMax,results.options.cN);
    c = c(:);
    for i = 1:length(results.approxResults)
        errPerParam = results.approxResults{i}.errPerParam;
        errPerParam = errPerParam(:) ./ results.dataInfo.nrmAPerParam(:);
        c = cat(2,c,errPerParam);
    end
    T = array2table(c, 'VariableNames', cat(2,'c',results.Mnames));
    writetable(T,sprintf([dirName,myDate,'/error--init-%s.csv'],MInit{1}));

    for i = [3,4,6,7]
        frontalSlice   = [];
        SkHat = modeProduct(results.Sk{i},results.M{i});
        AHat  = modeProduct(results.A,results.M{i});

        for k = 1:size(results.A,3)
    
            sk = diag(SkHat(:,:,k));
            csk = cumsum(sk.^2) / fronorm(AHat(:,:,k))^2;
            frontalSlice = cat(2,frontalSlice,csk(:));
        end
        T = array2table([(1:numel(csk))',frontalSlice]);
        writetable(T,sprintf([dirName,myDate,'/sigma--M-%s--init-%s.csv'],...
            results.Mnames{i},MInit{1}));
    end
    
end

%% plot error per parameter

for MInit = {'C','I','Z'}
    load(sprintf([dirName,myDate,'/',myDate,'--k-02--init-%s'],MInit{1}));
    romPlotErrorPerParam(results)
    exportgraphics(gcf,[dirName,myDate,'/errPerParam_',MInit{1},'.png'],'BackgroundColor','none')
end


%% plot approixmation
load([dirName,myDate,'/',myDate,'--k-02--init-I'])
idx = 7;
romSaveApproximationPlot(results,idx)

%% visualize matrices

idx = [3,4,6];
MInit = {'C','I','Z'};
for i = 1:3
    load(sprintf([dirName,myDate,'/',myDate,'--k-02--init-%s'],MInit{i}));

    subDirName = [dirName,myDate,'/Mpics/'];
    if ~isfolder(subDirName), mkdir(subDirName); end
    romVisualizeMatrices(results,subDirName)
    [results.approxResults{idx(i)}.err,results.approxResults{end}.err] / results.dataInfo.nrmA
    exportgraphics(gcf,[subDirName,'matrix_',MInit{1},'.png'],'BackgroundColor','none')
end


%% visualize spectra
for MInit = {'C','I','Z'}
    load(sprintf([dirName,myDate,'/',myDate,'--k-02--init-%s'],MInit{1}));
    romVisualizeSpectrum(results)
    exportgraphics(gcf,[dirName,myDate,'/','spectrum_',MInit{1},'.png'],'BackgroundColor','none')
end


