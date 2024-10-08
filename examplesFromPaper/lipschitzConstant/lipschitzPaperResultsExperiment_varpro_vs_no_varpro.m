% clear variables
clear; clc;

% setup paths
run("../../starMOptSetup.m")



defaultOptions = lipschitzExperimentParameters();

options = {};
count   = 1;

for fctnType = {'ls-varpro', 'ls-no-varpro'}
    for condA = logspace(0,10,6)
        for myeps = logspace(0,-10,6)
            tmpOptions          = defaultOptions;
            tmpOptions.fctnType = fctnType{1};
            tmpOptions.condA    = condA;
            tmpOptions.probType = 'ill-conditioned';
            tmpOptions.MType    = 'close';
            tmpOptions.XType    = 'optimal';
            tmpOptions.perturbI = true;
            tmpOptions.eps      = myeps;
            options{count}      = obj2struct(tmpOptions);
            count               = count + 1;
        end
    end
end


%% run experiment

headers = {'eps', 'condA', 'meanL','stdL','meanG1G2', 'stdG1G2', 'meanM1M2', 'stdM1M2'};
results = [];
storeF  = {};
for i = 1:length(options)
    out = lipschitzRun(options{i});
    
    % Lipschitz
    L = out(:,1) ./ out(:,2);

    results = cat(1,results,[options{i}.eps, options{i}.condA, mean(L),std(L),...
        mean(out(:,1)), std(out(:,1)), ...
        mean(out(:,2)), std(out(:,2))]);
    storeF = cat(1,storeF,options{i}.fctnType);
end

T           = array2table(results,'VariableNames',headers);
T.fctnType  = storeF;

if ~exist('./resultsLipschitz','dir'), mkdir('./resultsLipschitz'); end
writetable(T, 'resultsLipschitz/vapro_vs_no_varpro.csv')

return;
%% 

for myCond = logspace(0,10,6)
    T = readtable('resultsLipschitz/vapro_vs_no_varpro.csv');
    T = T(T.condA == myCond,:);

    T.meanLPlus  = T.meanL + T.stdL;
    T.meanLMinus = T.meanL - T.stdL;

    T.meanG1G2Plus  = T.meanG1G2 + T.stdG1G2;
    T.meanG1G2Minus = T.meanG1G2 - T.stdG1G2;

    T.meanM1M2Plus  = T.meanM1M2 + T.stdM1M2;
    T.meanM1M2Minus = T.meanM1M2 - T.stdM1M2;

    
    for myF = {'ls-varpro','ls-no-varpro'}
        writetable(T(strcmp(T.fctnType,myF{1}),:), ...
            sprintf('resultsLipschitz/vapro_vs_no_varpro_%s_cond_%d.csv',myF{1},abs(int8(log10(myCond)))))
    end
end


