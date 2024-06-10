function[results] = digitsRun(options, saveFlag)

if ~exist('options', 'var'), options = digitsExperimentParameters(); options = obj2struct(options); end
if ~exist('saveFlag', 'var') || isempty(saveFlag), saveFlag = false; end

if saveFlag
    if ~isfolder('digitsResults/'), mkdir('digitsResults/'); end
    diary(['digitsResults/',options.filename,'.txt'])
end



% ----------------------------------------------------------------------- %
rng(options.seed);

% setup data
[A,labels,nrmAPerClass,imgIdx] = digitsSetupData(options);

% setup objective function
fctn = objFctnLowRank(A,options.k,'avgFlag',0);

% names and solutions we store
names    = {'matrix1', 'matrix2', 'F', 'I', 'C', 'Z^T', 'Q', 'MOpt'};

% transformations we use ([] means no transformation yet)
[Z,~,~] = svd(modeUnfold(A), 'econ');
[Q,~]   = qr(randn(size(A,3)));
M       = {[], [], [], eye(size(A,3)), dctmtx(size(A,3)), Z', Q, []};

% ----------------------------------------------------------------------- %
% learn the best choice of M
opt = gradientDescent('verbose',1,'maxIter',options.maxIter,'manifold','Stiefel','logInterval',1);
opt.linesearch = armijoLinesearch('alpha',1e0,'manifold','Stiefel','gamma',1e-3,'maxIter',40);

% solve! TODO: should we try different initializations for M0?
rng(options.seedM);

if strcmp(options.initializeM,'identity')
    M0 = eye(size(A,3));
else
    [M0,~] = qr(randn(size(A,3)));
end
% 

startTime           = tic;    
[MSol,optInfo]      = opt.solve(fctn,M0);
endTime             = toc(startTime);
optInfo.totalTime   = endTime; 

% store learned transformation
M{end} = MSol;

% ----------------------------------------------------------------------- %
% results on training data

[Ak,Uk,Sk,Vk,approxResults] = digitsComputeResults(A,M,labels,options);

% ----------------------------------------------------------------------- %
% transfer learning statistics

transferResults = digitsTransfer(imgIdx,M,labels,names,options);


% ----------------------------------------------------------------------- %
% results

%
% store results
results.Ak                  = Ak;
results.Uk                  = Uk;
results.Sk                  = Sk;
results.Vk                  = Vk;
results.approxResults       = approxResults;
results.transferResults     = transferResults;
results.M                   = M;
results.Mnames              = names;
results.optInfo             = optInfo;
results.options             = options;
results.dataInfo            = struct('nrmA',fronorm(A),'nrmAPerClass',nrmAPerClass,...
    'sizeA',size(A),'numelA',numel(A),'storageUV',numel(Uk{end}) + numel(Vk{end}),'storageM',numel(M{end}));

% save results
if saveFlag
    save(['digitsResults/', options.filename], 'results')
end

end




