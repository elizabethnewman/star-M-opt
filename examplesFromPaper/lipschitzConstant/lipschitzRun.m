function[out] = lipschitzRun(options)

% reproducibility
rng(options.seed);

% create data
switch options.probType
    case 'ill-conditioned'
        A = zeros(options.n1,options.m,options.n3);
        for i = 1:options.n3
            u = orth(randn(options.n1, options.m));
            v = orth(randn(options.m));
            A(:,:,i) = u * diag(logspace(0, -log10(options.condA), options.m)) * v';
        end
        A = options.nrmA * (A / fronorm(A));

    case 'index-tracking'
         options.fctnType = 'ls-constrained-varpro';
         [A,B]            = indexTrackingData();
         options.n1       = size(A,1);
         options.n2       = size(A,2);
         options.n3       = size(A,3);
         options.m        = size(B,2);

    otherwise
        A = randn(options.n1,options.n2,options.n3);
        A = options.nrmA * (A / fronorm(A));
end

% left-hand side
if ~strcmp(options.probType,'index-tracking') 
    B = randn(options.n1,options.m,options.n3);
    B = options.nrmB * (B / fronorm(B));
end

% choose objective function
switch options.fctnType
    case 'ls-varpro'
        fctn = objFctnTikhonovLeastSquaresVarPro(A,B);
    case 'ls-no-varpro'
        fctn = objFctnTikhonovLeastSquares(A,B);
    case 'ls-constrained-varpro'
        fctn = objFctnConstrainedLeastSquaresVarPro(A,B);
    case 'low-rank'
        fctn = objFctnLowRank(A, options.k);
    otherwise
        error('nyi')
end

switch options.XType
    case 'optimal'
        tmp = objFctnTikhonovLeastSquaresVarPro(A,B);
        X = tmp.solve(eye(options.n3));
    otherwise
        X = randn(size(A,2),size(B,2),options.n3);
end


%% main experiment

out = estimateLipschitzConstant(fctn,options,X);

end



%% helper functions

function[out] = estimateLipschitzConstant(fctn,options,X)

out = zeros(options.nRuns,2);
for i = 1:options.nRuns
    if options.perturbI
        % M1 = orth(eye(n3) + options.eps * randn(options.n3));
        Omega   = randn(options.n3);
        M1      = expm(options.eps * (Omega - Omega'));
    else
        M1 = orth(randn(options.n3));
    end
    
    switch options.MType
        case 'close'
            Omega = randn(options.n3);
            M2    = M1 * expm(options.eps  * (Omega - Omega'));
        otherwise
            M2 = orth(randn(options.n3));
    end

    if isa(fctn, 'objFctnTikhonovLeastSquares')
        [~,~,G1] = fctn.evaluate(M1,X);
        [~,~,G2] = fctn.evaluate(M2,X);
    else
        [~,~,G1] = fctn.evaluate(M1);
        [~,~,G2] = fctn.evaluate(M2);
    end

    out(i,1) = norm(G1 - G2,'fro');
    out(i,2) = norm(M1 - M2,'fro');

end

end

function[A,B] = indexTrackingData()
options                  = indexTrackingExperimentParameters();
options.initDateTrain    = '1-Jan-2023';
options.endDateTrain     = '31-May-2023';
options.initDateTest     = '1-Jun-2023';
options.endDateTest      = '31-Jul-2023';
options.maxIter          = 100;
options.alpha            = 1e-2;

% update filename
options.setFilename();

% get struct
options = obj2struct(options);

% load data
[A,B,~,~,~] = indexTrackingSetupData(options.initDateTrain,options.endDateTrain,'first');

end