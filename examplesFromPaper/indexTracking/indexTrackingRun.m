function[results] = indexTrackingRun(options,saveFlag)

if ~exist('options', 'var'), options = indexTrackingExperimentParameters(); options = obj2struct(options); end
if ~exist('saveFlag', 'var') || isempty(saveFlag), saveFlag = false; end

if saveFlag
    dirName = ['indexTrackingResults/',options.filename,'/'];
    if ~isfolder(dirName), mkdir(dirName); end
    diary([dirName,options.filename,'.txt'])
end

% ----------------------------------------------------------------------- %
% setup problem

% get training data
[A,B,BMat,sectors,SP500] = indexTrackingSetupData(options.initDateTrain,options.endDateTrain,'first');

% setup objective function
fctn = objFctnConstrainedLeastSquaresVarPro(A,B,[],options.alpha);

% names and solutions we store
names               = {'matrix', 'I', 'C', 'Z^T', 'MOpt'};
X                   = cell(size(names));
evalInfo            = cell(size(names));
indexApprox         = cell(size(names));
trackingResults     = cell(size(names));
sectorDistribution  = cell(size(names));

% transformations we use ([] means no transformation yet)
[Z,~,~] = svd(modeUnfold(A), 'econ');
M       = {[], eye(size(A,3)), dctmtx(size(A,3)), Z', []};

% make sure initialization of X0 is the same for all methods
rng(options.seed);
X0 = rand(size(A,2),1,size(A,3));

%% learn the best choice of M
opt = gradientDescent('verbose',1,'maxIter',options.maxIter,'manifold','Stiefel','logInterval',1);
opt.linesearch = armijoLinesearch('alpha',1e0,'manifold','Stiefel','gamma',1e-3,'maxIter',40);

% solve! TODO: try different initializations for M0
rng(options.seedM);
[M0,~] = qr(randn(size(A,3)));

startTime       = tic;    
[MSol,optInfo]  = opt.solve(fctn,M0);
endTime         = toc(startTime);
optInfo.totalTime = endTime; 

% store learned transformation
M{end} = MSol;


%% Fit matrix data
fMatrix     = @(x) objFctnMatrix(x,A,BMat,fctn.alpha);
Aeq         = ones(1,numel(X0));
beq         = 1;
lb          = zeros(size(X0));
x           = fmincon(fMatrix,X0(:),[],[],Aeq,beq,lb,[],[],fctn.options);
X{1}        = x;
evalInfo{1} = fMatrix(x);
indexApprox{1}          = A(:,:) * x;
trackingResults{1}      = fronorm(indexApprox{1} - BMat);
sectorDistribution{1}   = reshape(sum(reshape(x,[],1,11),[1,2]),[],1);

% matrixResultsPerSector = struct('X',cell(1,11),'indexApprox',cell(1,11));

for i = 1:11
    fMatrix = @(x) objFctnMatrix(x,A(:,:,i),B(:,1,i),fctn.alpha);
    x0      = X0((i-1)*size(A,2)+1:i*size(A,2));
    Aeq     = ones(1,numel(x0));
    beq     = 1;
    lb      = zeros(size(x0));
    x       = fmincon(fMatrix,x0(:),[],[],Aeq,beq,lb,[],[],fctn.options);

    matrixResultsPerSector.X{i}                 = x;
    matrixResultsPerSector.indexApprox{i}       = A(:,:,i) * x;
    matrixResultsPerSector.trackingResults{i}   = fronorm(matrixResultsPerSector.indexApprox{i} - B(:,1,i));
end

%% Fit tensor data for known transformation
for i = 2:length(M)
    [~,tmpInfo]           = fctn.evaluate(M{i}, X0);
    X{i}                  = tmpInfo.X;
    evalInfo{i}           = fctn.evaluate(M{i}, X{i});
    indexApprox{i}        = mprod(A,X{i},M{i});
    trackingResults{i}    = fronorm(indexApprox{i} - B);
    sectorDistribution{i} = reshape(sum(X{i},[1,2]),[],1);
end


%% prediction using various indexTracking weights
[predictResults,trackingResultsTest,trackingResultsTestMatrixPerSector] = indexTrackingPrediction(X, M, matrixResultsPerSector, options.initDateTest, options.endDateTest);
matrixResultsPerSector.trackingResults = trackingResultsTestMatrixPerSector;

%% store results in one struct
results.X                   = X;
results.approx              = indexApprox;
results.sectors             = sectors;
results.SP500               = SP500;
results.sectorDistribution  = sectorDistribution;
results.M                   = M;
results.evalInfo            = evalInfo;
results.Mnames              = names;
results.optInfo             = optInfo;
results.optOptions          = fctn.options;
results.predictResults      = predictResults;
results.trackingResults     = trackingResultsTest;
results.options             = options;
results.sectors             = sectors;
results.matrixPerSector     = matrixResultsPerSector;

% save results
if saveFlag
    save([dirName, options.filename], 'results')
end

end 



function[f,df] = objFctnMatrix(x,A,B,alpha)

n   = 1 / size(A,1);
res = A(:,:) * x(:) - B(:);
f   = (0.5 * n) * fronorm(res)^2 + (0.5 * alpha) * fronorm(x)^2;
df  = (1 / n) * (A(:,:)' * res) + alpha * x(:);

end

