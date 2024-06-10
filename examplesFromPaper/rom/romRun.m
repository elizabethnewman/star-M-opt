function[results] = romRun(options, saveFlag)

if ~exist('options', 'var'), options = romExperimentParameters(); options = obj2struct(options); end
if ~exist('saveFlag', 'var') || isempty(saveFlag), saveFlag = false; end

if saveFlag
    if ~isfolder('romResults/'), mkdir('romResults/'); end
    if ~isfolder(['romResults/',options.date]), mkdir(['romResults/',options.date]); end
    diary(['romResults/',options.filename,'.txt'])
end

% ----------------------------------------------------------------------- %

% setup data
A = romSetupData(options);

% setup objective function
fctn = objFctnLowRank(A,options.k,'avgFlag',0);

% names and solutions we store
names   = {'matrix', 'F', 'I', 'C', 'Q', 'Z^T', 'MOpt'};

% transformations we use ([] means no transformation yet)
rng(options.seedM);
[Z,~,~] = svd(modeUnfold(A), 'econ');
[Q,~]   = qr(randn(size(A,3)));
M       = {[], [], eye(size(A,3)), dctmtx(size(A,3)), Q, Z', []};

% ----------------------------------------------------------------------- %
% learn the best choice of M
opt = gradientDescent('verbose',1,'maxIter',options.maxIter,'manifold','Stiefel','logInterval',1);
opt.linesearch = armijoLinesearch('alpha',1e0,'manifold','Stiefel','gamma',1e-3,'maxIter',40);

% solve! TODO: should we try different initializations for M0?
rng(options.seedM);

switch options.MInitialize
    case 'random'
        M0 = Q;
    case 'Z'
        M0 = Z';
    case 'C'
        M0 = dctmtx(size(A,3));
    case 'I'
        M0 = eye(size(A,3));
    otherwise
        error('Choose M intialization strategy')
end

startTime           = tic;    
[MSol,optInfo]      = opt.solve(fctn,M0);
endTime             = toc(startTime);
optInfo.totalTime   = endTime; 

% store learned transformation
M{end} = MSol;

% ----------------------------------------------------------------------- %
% results on training data

[Ak,Uk,Sk,Vk,approxResults] = romComputeResults(A,M,options.k,options);

% ----------------------------------------------------------------------- %
% results

% store results (including data A because time-consuming to generate)
results.A                   = A;
results.Ak                  = Ak;
results.Uk                  = Uk;
results.Sk                  = Sk;
results.Vk                  = Vk;
results.approxResults       = approxResults;
results.M                   = M;
results.Mnames              = names;
results.optInfo             = optInfo;
results.options             = options;
results.dataInfo            = struct('nrmA',fronorm(A),'nrmAPerParam',reshape(sqrt(sum(A.^2,[1,2])),1,[]),...
    'sizeA',size(A),'numelA',numel(A),'storageUV',numel(Uk{end}) + numel(Vk{end}),'storageM',numel(M{end}));


% save results
if saveFlag
    save(['romResults/', options.filename], 'results')
end

end




