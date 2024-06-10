
% choose options

defaultOptions = lsExperimentParameters();

options = {};
count = 1;

% for eta = [1e-1,1e-2,1e-3,1e-4,0]
for eta = 1e-3
    
    % choose parameters
    tmpOptions = defaultOptions;
    tmpOptions.eta = eta;

    % update filename
    tmpOptions.setFilename();

    % update counter
    options{count} = obj2struct(tmpOptions);
    count = count + 1;
end

% main loop
saveFlag = 0;
for i = 1:length(options)
    results = lsRun(options{i}, saveFlag);
end

