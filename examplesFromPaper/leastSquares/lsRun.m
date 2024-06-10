function[results] = lsRun(options, saveFlag)

if ~exist('options', 'var'), options = lsExperimentParameters(); options = obj2struct(options); end
if ~exist('saveFlag', 'var') || isempty(saveFlag), saveFlag = false; end

if saveFlag
    dirName = ['lsResults/',options.filename,'/'];
    if ~isfolder(dirName), mkdir(dirName); end
    diary([dirName,options.filename,'.txt'])
end

% ----------------------------------------------------------------------- %
rng(options.seed);

% setup data
[A,B] = lsSetupData(options);

% setup objective function
fctn = objFctnTikhonovLeastSquaresVarPro(A,B,options.lambda,options.alpha);

% names and solutions we store
names    = {'matrix', 'F', 'I', 'C', 'Z^T', 'Q', 'MOpt'};

% transformations we use ([] means no transformation yet)
[Z,~,~] = svd(modeUnfold(A), 'econ');
[Q,~]   = qr(randn(size(A,3)));
M       = {[], [], eye(size(A,3)), dctmtx(size(A,3)), Z', Q, []};

% ----------------------------------------------------------------------- %
% learn the best choice of M
opt  = gradientDescent('verbose',1,'maxIter',options.maxIter,'manifold','Stiefel','logInterval',1,'storeIterates',1);
opt.linesearch = armijoLinesearch('alpha',1e2,'manifold','Stiefel','gamma',1e-3,'maxIter',40);

% solve! TODO: should we try different initializations for M0?
rng(options.seedM);
[M0,~] = qr(randn(size(A,3)));
% M0 = eye(size(A,3));

startTime           = tic;    
[MSol,optInfo]      = opt.solve(fctn,M0);
endTime             = toc(startTime);
optInfo.totalTime   = endTime; 

% store learned transformation
M{end} = MSol;

% ----------------------------------------------------------------------- %
% results on training data

[X,approxResults] = lsComputeResults(A,B,M,options);

% ----------------------------------------------------------------------- %
% compare with alternating

optAlt = alternatingDescent('manifoldX','Stiefel','manifoldY','Euclidean',...
    'linesearchX',copy(opt.linesearch),'linesearchY',copy(opt.linesearch),...
    'verbose',1,'maxIter',4 * opt.maxIter,'logInterval',1,'storeIterates',1);

fctnAlt = objFctnTikhonovLeastSquares(A,B,fctn.lambda,fctn.alpha);

% X0 = randn(size(X{end}));
X0 = solveForX(A,B,M0);
[MAlt,XAlt,infoAlt] = solve(optAlt,fctnAlt,M0,X0);

fprintf('M Alternating: ||A * X - B||_F / ||B||_F = %0.4e\n',fronorm(mprod(A,XAlt,MAlt) - B) / fronorm(B))

% ----------------------------------------------------------------------- %
% results

% store results
results.X                   = X;
results.approxResults       = approxResults;
results.M                   = M;
results.Mnames              = names;
results.optInfo             = optInfo;
results.options             = options;


results.alternating.optInfo = infoAlt;
results.alternating.M       = MAlt;
results.alternating.X       = XAlt;


% save results
if saveFlag
    save([dirName,options.filename], 'results')
end

end



function[X] = solveForX(A,B,M)

AHat  = modeProduct(A,M);
BHat  = modeProduct(B,M);
XHAt  = pagemldivide(AHat,BHat);
X     = modeProduct(XHAt,M');

end

