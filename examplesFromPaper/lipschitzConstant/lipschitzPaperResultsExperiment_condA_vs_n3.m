% clear variables
clear; clc;

% setup paths
run("../../starMOptSetup.m")



defaultOptions = lipschitzExperimentParameters();

options = {};
count   = 1;

for condA = logspace(0,10,6)
    for n3 = 2.^(1:8)
        for Mtype = {'random','close'}
            tmpOptions          = defaultOptions;
            tmpOptions.condA    = condA;
            tmpOptions.n3       = n3;
            tmpOptions.probType = 'ill-conditioned';
            tmpOptions.MType    = Mtype{1};
            tmpOptions.perturbI = true;
            options{count}      = obj2struct(tmpOptions);
            count               = count + 1;
        end
    end
end


%% run experiment

headers = {'n3', 'condA', 'meanL', 'stdL', 'meanG1G2', 'stdG1G2', 'meanM1M2', 'stdM1M2'};
results = [];
storeM  = {};
for i = 1:length(options)
    out = lipschitzRun(options{i});
    
    L = out(:,1) ./ out(:,2);

    results = cat(1,results,[options{i}.n3, options{i}.condA, mean(L),std(L), ...
        mean(out(:,1)), std(out(:,1)), ...
        mean(out(:,2)), std(out(:,2))]);
    storeM = cat(1,storeM,options{i}.MType);
end

T = array2table(results,'VariableNames',headers);
T.MType = storeM;

if ~exist('./resultsLipschitz','dir'), mkdir('./resultsLipschitz'); end
writetable(T, 'resultsLipschitz/condA_vs_n3.csv')

return;
%% 

for myCond = logspace(0,10,6)
    T = readtable('resultsLipschitz/condA_vs_n3.csv');
    T = T(T.condA == myCond,:);

    T.meanLPlus  = T.meanL + T.stdL;
    T.meanLMinus = T.meanL - T.stdL;

    T.meanG1G2Plus  = T.meanG1G2 + T.stdG1G2;
    T.meanG1G2Minus = T.meanG1G2 - T.stdG1G2;

    T.meanM1M2Plus  = T.meanM1M2 + T.stdM1M2;
    T.meanM1M2Minus = T.meanM1M2 - T.stdM1M2;

    
    for myM = {'random','close'}
        writetable(T(strcmp(T.MType,myM{1}),:), ...
            sprintf('resultsLipschitz/condA_vs_n3_M1M2%s_cond_%d.csv',myM{1},abs(int8(log10(myCond)))))
    end
end

