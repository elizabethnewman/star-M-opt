


clear; clc;

% setup paths
run("../../starMOptSetup.m")

% setup directory
dirName    = 'convergenceResults/';

if ~exist(dirName,'dir'), mkdir(dirName); end
cd(dirName)

subDirName = [date,'/'];
if ~exist(subDirName,'dir'), mkdir(subDirName); end
cd(subDirName)

% return to main director
cd ../..


%% 

results.n = 1:4;
results.options = cell(1,length(results.n));

results.varpro.MSol             = cell(1,length(results.n));
results.varpro.M0               = cell(1,length(results.n));
results.varpro.optInfo          = cell(1,length(results.n));
results.varpro.approxResults    = cell(1,length(results.n));


results.alternating.MSol            = cell(1,length(results.n));
results.alternating.XSol            = cell(1,length(results.n));
results.alternating.optInfo         = cell(1,length(results.n));
results.alternating.approxResults   = cell(1,length(results.n));


for n = results.n 
    options = obj2struct(lsExperimentParameters('m',linspace(-1,1,2^n),'b',linspace(-1,1,2^n), 'maxIter',500 * n));
    results.options{n} = options;

    % setup data
    rng(options.seed)
    [A,B] = lsSetupData(options);
    
    % setup objective function
    fctn = objFctnTikhonovLeastSquaresVarPro(A,B,options.lambda,options.alpha);

    % ----------------------------------------------------------------------- %
    % learn the best choice of M using VarPro
    opt  = gradientDescent('verbose',0,'maxIter',options.maxIter,'manifold','Stiefel','logInterval',1,'storeIterates',1);
    opt.linesearch = armijoLinesearch('alpha',1e2,'manifold','Stiefel','gamma',1e-3,'maxIter',40);
    
    % solve!
    rng(options.seedM);
    [M0,~] = qr(randn(size(A,3)));
    % M0 = eye(2^n);
    
    % main iteration
    startTime           = tic;    
    [MSol,optInfo]      = opt.solve(fctn,M0);
    optInfo.totalTime   = toc(startTime); 

    XSol = fctn.solve(MSol);

    results.varpro.MSol{n}                 = MSol;
    results.varpro.M0{n}                   = M0;
    results.varpro.optInfo{n}              = optInfo;

    err = fronorm(mprod(A,XSol,MSol) - B);
    results.varpro.approxResults{n}.err         = err;
    results.varpro.approxResults{n}.relErr      = err / fronorm(B);
    results.varpro.approxResults{n}.totalTime   = optInfo.totalTime;

    fprintf('M VarPro: ||A * X - B||_F / ||B||_F = %0.4e\n',fronorm(mprod(A,XSol,MSol) - B) / fronorm(B))
    fprintf('Time VarPro: %0.4f\n', optInfo.totalTime);

    optInfo.header = {optInfo.header{1:2}, 'totalTime',optInfo.header{3:end}};
    optInfo.values = [optInfo.values(:,1:2), cumsum(optInfo.values(:,2)), optInfo.values(:,3:end)];

    T = array2table(optInfo.values, 'VariableNames',optInfo.header);
    writetable(T,sprintf([dirName,subDirName,'varpro_%0.3d.csv'],n))



    % ----------------------------------------------------------------------- %
    % learn the best choice of M using alternating descent
    optAlt = alternatingDescent('manifoldX','Stiefel','manifoldY','Euclidean',...
        'linesearchX',copy(opt.linesearch),'linesearchY',copy(opt.linesearch),...
        'verbose',0,'maxIter',opt.maxIter,'logInterval',1,'storeIterates',1);

    fctnAlt = objFctnTikhonovLeastSquares(A,B,fctn.lambda,fctn.alpha);
    
    % same initial guess
    X0 = modeProduct(pagemldivide(modeProduct(A,M0),modeProduct(B,M0)),M0');

    startTime           = tic;   
    [MAlt,XAlt,altInfo] = solve(optAlt,fctnAlt,M0,X0);
    altInfo.totalTime   = toc(startTime);

    results.alternating.MSol{n}                 = MAlt;
    results.alternating.XSol{n}                 = XAlt;
    results.alternating.optInfo{n}              = altInfo;

    err = fronorm(mprod(A,XAlt,MAlt) - B);
    results.alternating.approxResults{n}.err        = err;
    results.alternating.approxResults{n}.relErr     = err / fronorm(B);
    results.alternating.approxResults{n}.totalTime  = altInfo.totalTime;

    fprintf('M Alternating: ||A * X - B||_F / ||B||_F = %0.4e\n',fronorm(mprod(A,XAlt,MAlt) - B) / fronorm(B))
    fprintf('Time Alternating: %0.4f\n', altInfo.totalTime);

    altInfo.header{end-2} = 'alphaLS2';
    altInfo.header{end-1} = 'iterLS2';
    altInfo.header{end} = 'flagLS2';
    altInfo.header = {altInfo.header{1:2}, 'totalTime',altInfo.header{3:end}};
    altInfo.values = [altInfo.values(:,1:2), cumsum(altInfo.values(:,2)), altInfo.values(:,3:end)];
    T = array2table(altInfo.values, 'VariableNames',altInfo.header);
    writetable(T,sprintf([dirName,subDirName,'alter_%0.3d.csv'],n))

end

save([dirName,subDirName,'vp-vs-ad'], "results")






