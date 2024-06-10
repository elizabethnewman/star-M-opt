
% clear variables
clear; clc;

% setup paths
run("../../starMOptSetup.m")

%% setup options
defaultOptions = lsExperimentParameters();

options = {};
count = 1;

for eta = [1e-1,1e-2,1e-3,1e-4,0]
    
    % choose parameters
    tmpOptions     = defaultOptions;
    tmpOptions.eta = eta;

    % update filename
    tmpOptions.setFilename();

    % update counter
    options{count} = obj2struct(tmpOptions);
    count = count + 1;
end

%% run experiments (can be done in parallel)

saveFlag = 1;
for i = 1:length(options)
    lsRun(options{i}, saveFlag);
end

%% plot data

% choose date of experiments
myDate = date;

% choose main directory
dirName = 'lsResults/';

% be sure to plot without docking
close all;
set(0,'DefaultFigureWindowStyle','normal')

%%
% save csv files
Q = @(theta) [cos(theta), sin(theta), -sin(theta), cos(theta)];

% all possible rotations
MTrue = reshape(Q(pi / 4),2,[]);

MCand = {MTrue, MTrue([2,1],:), -MTrue, -MTrue([2,1],:), diag([-1,1]) * MTrue, diag([1,-1]) * MTrue, diag([1,-1]) * MTrue([2,1],:), diag([-1,1]) * MTrue([2,1],:)};

for eta = [1e-1,1e-2,1e-3,1e-4,0]
    load(sprintf([dirName,myDate,'--eta-%0.0e/',myDate,'--eta-%0.0e'],eta,eta))

    T = array2table(results.optInfo.values, 'VariableNames', results.optInfo.header);
    writetable(T,sprintf([dirName,myDate,'--eta-%0.0e/',myDate,'--eta-%0.0e.csv'],eta,eta))

    tmp = results.optInfo.x';

    theta = atan(tmp(:,2) ./ tmp(:,1));

    % print(norm(tmp - QQ))
    tmp = cat(2,tmp,theta);
    

    for m = MCand
        % disp(m{1})
        d1 = vecnorm(tmp(:,1:4) - m{1}(:)',2,2);
        tmp = cat(2,tmp,d1);
    end
    tmp = cat(2,tmp,min(tmp(:,end-7:end),[],2));

    % add iteration counter
    tmp = cat(2,(0:size(tmp,1)-1)',tmp);

    T = array2table(tmp, 'VariableNames', {'iter','M11', 'M21', 'M12', 'M22', 'theta','D1','D2','D3','D4','D5','D6','D7','D8','Dmin'});
    writetable(T,sprintf([dirName,myDate,'--eta-%0.0e/',myDate,'--eta-%0.0e--matrix.csv'],eta,eta))

    disp(reshape(tmp(end,2:5),[],2))

end

%% plot data
load(sprintf([dirName,myDate,'--eta-%0.0e/',myDate,'--eta-%0.0e'],0,0))

[A,B] = lsSetupData(results.options);

fig = lsPlotHyperplanes(A,B,results.X{end},results.M{end},0,0);
saveas(fig,[dirName,myDate,'--transformDomain.png'])

fig = lsPlotHyperplanes(A,B,results.X{end},results.M{end},0,1);
saveas(fig,[dirName,myDate,'--spatialDomain.png'])

%% plot convergence

fig = figure(1); clf;
plotOptions = {'LineWidth',4,'MarkerSize',15};
markers = {'o', 'x', '^', '+', 's'};
count = 1;
for eta = [0,1e-4,1e-3,1e-2,1e-1]
    load(sprintf([dirName,myDate,'--eta-%0.0e/',myDate,'--eta-%0.0e'],eta,eta))

    semilogy(results.optInfo.values(:,1),results.optInfo.values(:,3), ['-',markers{count}],plotOptions{:}, 'DisplayName',sprintf('eta = %0.2e',eta))
    hold on;
    count = count + 1;
end

xlabel('iteration')
ylabel('$\overline{\Phi}(\mathbf{M})$','Interpreter','latex')
grid;
legend('Location','southwest')

hold off;

saveas(fig,['lsResults/',myDate,'--convergence.png'])

%% plot approximations

% create hyperplanes for select M
for eta = [0,1e-4,1e-3,1e-2,1e-1]
    localDirName = sprintf([dirName,myDate,'--eta-%0.0e/',myDate,'--eta-%0.0e'],eta,eta);
    if ~exist(localDirName,'dir'), mkdir(localDirName); end

    % load results
    load(localDirName)

    % setup data
    [A,B] = lsSetupData(results.options);

    for i = 1:length(results.M)
        fig = lsPlotHyperplanes(A,B,results.X{i},results.M{i},1,0);
        axis('off')
        subDirname = [localDirName,'/',results.Mnames{i}];
        if ~exist(subDirname,'dir'), mkdir(subDirname); end
        saveas(fig,[subDirname,'/','approx_transformDomain.png'])
    end

    fctn = objFctnTikhonovLeastSquaresVarPro(A,B,results.options.lambda,results.options.alpha);
    XTrue = fctn.solve(results.options.M);
    fig = lsPlotHyperplanes(A,B,XTrue,results.options.M,1,0);
    axis('off')

    subDirname = [localDirName,'/true'];
    if ~exist(subDirname,'dir'), mkdir(subDirname); end
    saveas(fig,[subDirname,'approx_transformDomain.png'])
end

%% plot training

dirName = sprintf([dirName,myDate,'--eta-%0.0e/',myDate,'--eta-%0.0e'],1e-2,1e-2);
subDirname = [dirName,'/',results.Mnames{end}];
if ~exist(subDirname,'dir'), mkdir(subDirname); end

% plot first 10 iterations
for i = 1:10
    fig = lsPlotTraining(results,i);
    axis('off')
    saveas(fig,sprintf([subDirname,'/','approx_transformDomain_iter_%d.png'],i))
end


