function[transferResults] = digitsTransfer(imgIdx,M,labels,names,options)

if ~exist('options', 'var'), options = digitsExperimentParameters(); options = obj2struct(options); end

rng(options.seedTransfer)

orthoFlag = {'orthoFlag',1};

transferResults.err = zeros(options.nTestRuns,length(M) + 3);

% allow for retraining with few iterations
opt = gradientDescent('verbose',1,'maxIter',options.innerMaxIter,'manifold','Stiefel','logInterval',1,'verbose',0);
opt.linesearch = armijoLinesearch('alpha',1e0,'manifold','Stiefel','gamma',1e-3,'maxIter',40);


% load data (rotated MNIST)
[XTrain,YTrain] = digitTrain4DArrayData;
classes = double(YTrain) - 1;

idx             = 1:size(XTrain,4);
idx(imgIdx)     = []; % remove training images
classes(imgIdx) = [];

for i = 1:options.nTestRuns
    % select a random batch
    switch options.transferDistribution
        case 'random'
            idxTest = randperm(length(idx),options.nImg * options.nClasses);
        case 'training'
            idxTest = [];
            for j = 0:options.nClasses-1
                idxLabels = find(classes == j);
                idxPerm   = randperm(length(idxLabels),options.nImg);
                idxTest   = cat(1,idxTest,idxLabels(idxPerm));
            end
    end

    % setup test data (normalized)
    ATest   = double(squeeze(XTrain(:,:,:,idx(idxTest))));
    ATest   = ATest / fronorm(ATest);
    ATest   = ATest + options.noiseLevel * randn(size(ATest));

    if options.permuteFlag
        ATest = permute(ATest,[1,3,2]);
    end

    % compute results with new Z^T
    [Z,~,~] = svd(modeUnfold(ATest),"econ");

    % retrain with few iterations
    fctn    = objFctnLowRank(ATest,options.k,'avgFlag',0);
    MSolI   = opt.solve(fctn,eye(size(M{end})));
    MSolM   = opt.solve(fctn,M{end});

    % compute results with test data
    [~,~,~,~,approxResults] = digitsComputeResults(ATest,cat(2,M,{Z',MSolI,MSolM}),labels,options);
    
    MRetrain   = {Z',MSolI,MSolM};
    errRetrain = zeros(1,length(MRetrain));
    for j = 1:length(MRetrain)
        [u,s,v]          = tSVDM(ATest,MRetrain{j},options.k);
        Ak               = mprod(u,mprod(s,tran(v),MRetrain{j},orthoFlag{:}),MRetrain{j},orthoFlag{:});
        errRetrain(j)    = fronorm(ATest - Ak);
    end

    % TODO: collect other approxResults outputs
    transferResults.idxTest{i} = idxTest;
    transferResults.labels{i}  = double(YTrain(idx(idxTest))) - 1;

    err = cellfun(@(x) x.err,approxResults);
    transferResults.err(i,:) = err;
    % transferResults.err(i,:) = [err,errRetrain];
end
transferResults.names   = cat(2,names,{'trainZ','trainI','trainMOpt'});
transferResults.meanErr = mean(transferResults.err,1);
transferResults.stdErr  = std(transferResults.err,1);


end

