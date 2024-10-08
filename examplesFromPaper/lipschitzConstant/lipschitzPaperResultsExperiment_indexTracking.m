% clear variables
clear; clc;

% setup paths
run("../../starMOptSetup.m")


options          = lipschitzExperimentParameters();
options.probType = 'index-tracking';


%% run experiment

headers = {'L','G1G2','M1M2','meanG1G2', 'meanL','stdL','stdG1G2', 'meanM1M2', 'stdM1M2'};
results = [];

out = lipschitzRun(options);

L = out(:,1) ./ out(:,2);

tmp = [mean(L),std(L),mean(out(:,1)),std(out(:,1)),mean(out(:,2)),std(out(:,2))];
out = cat(2,L,out,ones(size(out,1),1) * tmp);

T = array2table(out,'VariableNames',headers);

if ~exist('./resultsLipschitz','dir'), mkdir('./resultsLipschitz'); end
writetable(T, 'resultsLipschitz/index_tracking.csv')



