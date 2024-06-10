
% clear variables
clear; clc;

% setup paths
run("../../starMOptSetup.m")


%% 

defaultOptions = digitsExperimentParameters();

options = {};
count   = 1;
kRange  = 1:5;

for k = kRange
    
    % choose parameters
    tmpOptions = defaultOptions;
    tmpOptions.k            = k;
    tmpOptions.maxIter      = 100;
    tmpOptions.nImg         = 20;
    tmpOptions.noiseLevel   = 0.0;
    tmpOptions.nTestRuns    = 20;
    tmpOptions.initializeM  = 'identity';
    tmpOptions.permuteFlag  = false;
    tmpOptions.transferDistribution = 'training';

    % update filename
    tmpOptions.setFilename();

    % update counter
    options{count} = obj2struct(tmpOptions);
    count = count + 1;
end

%%
% main loop
saveFlag = 1;
for i = 1:length(options)
    results = digitsRun(options{i}, saveFlag);
end

%% 

% choose date of experiments 
myDate = date;

% choose main directory
dirName = 'digitsResults/';

% be sure to plot without docking
close all;
set(0,'DefaultFigureWindowStyle','normal')


%% convergence results
fig = figure(1); clf;
for k = kRange
    load(sprintf([dirName,myDate,'--k-%0.2d--init-identity'],k));
    T = array2table(results.optInfo.values, 'VariableNames', results.optInfo.header);
    writetable(T,sprintf([dirName,myDate,'--k-%0.2d--init-identity.csv'],k));
    semilogy(results.optInfo.values(:,1),results.optInfo.values(:,3),'-o','DisplayName',['k = ', num2str(k)], 'LineWidth',3, 'MarkerSize', 10)
    hold on;
end
legend()
xlabel('iter')
ylabel('loss')
grid on
hold off;
set(gca,'FontSize',18)

exportgraphics(fig,[dirName,'convergence.png'],'BackgroundColor','none')


%% print relative error for LateX


T = [];
for k = kRange
    load(sprintf([dirName,myDate,'--k-%0.2d--init-identity'],k));

    tmp = cellfun(@(x) x.err, results.approxResults, 'UniformOutput',true);
    T = cat(1,T,[k,tmp]);

    fprintf('%d & %0.4f & %0.4f & %0.4f & %0.4f & %0.4f & %0.4f & %0.4f & %0.4f\\\\\n',[k,tmp])
end
headers = cat(2,'k',results.Mnames);
T = array2table(T, 'VariableNames', headers);
    
writetable(T,[dirName,myDate,'--relative_err--init-identity.csv']);

%% approximations

% create hyperplanes for select M
for k = kRange
    load(sprintf([dirName,myDate,'--k-%0.2d--init-identity'],k));
    digitsSaveApproximationPlot(results)
end

%% relative error

B = [];
nrmA = results.dataInfo.nrmA;
for k = kRange
    load(sprintf([dirName,myDate,'--k-%0.2d--init-identity'],k));

    for i = 1:length(results.approxResults)

        B(k,i) = results.approxResults{i}.err / nrmA;
    end
end

fig = heatmap(B);
fig.XDisplayLabels = results.Mnames;
exportgraphics(fig,[dirName,'relativeError.png'],'BackgroundColor','none')


%% transfer learning

T = [];
for k = kRange
    load(sprintf([dirName,myDate,'--k-%0.2d--init-identity'],k));
    T = cat(1,T,results.transferResults.meanErr);
end

T = array2table([(kRange)',T], 'VariableNames', cat(2,'k',results.transferResults.names));
writetable(T,[dirName,myDate,'--transfer--init-identity.csv']);

% fig = heatmap(T);
% fig.XDisplayLabels = results.transferResults.names;
% exportgraphics(fig,[dirName,'transfer.png'],'BackgroundColor','none')

%% features in transform domain

% create hyperplanes for select M
load(sprintf([dirName,myDate,'--k-%0.2d--init-identity'],1));

% ----------------------------------------------------------------------- %
rng(results.options.seed);

% setup data
[A,labels,nrmAPerClass,imgIdx] = digitsSetupData(results.options);

idx = [];
for i = unique(labels)
    idx1 = find(labels == i);
    idx  = cat(2,idx,idx1(1:3));
end

cMin = min(A(:));
cMax = max(A(:));
for i = [4,6,8]
    AHat = modeProduct(A,results.M{i});
    cMin  = min(cMin,min(AHat(:)));
    cMax  = max(cMax,max(AHat(:)));
end

for i = [4,6,8]
    subDirName = [dirName,'features/'];
    if ~exist(subDirName,'dir'), mkdir(subDirName); end

    AHat = modeProduct(A,results.M{i});

    tmpDir = [subDirName,results.Mnames{i},'/'];
    disp(tmpDir)
    
    if ~exist(tmpDir,'dir') || exist(tmpDir,'dir') == 7, mkdir(tmpDir); end

    count = 1;
    for j = idx
        fig = figure(1); clf;
        imagesc(AHat(:,:,j));
        axis('image')
        axis('off')
        colormap parula;
        % clim([cMin,cMax])
        % colormap hot;

        
        exportgraphics(fig,[tmpDir,'feature_',num2str(count),'.png'],'BackgroundColor','none')
        count = count + 1;
    end
end



