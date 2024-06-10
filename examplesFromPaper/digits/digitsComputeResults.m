function[Ak,Uk,Sk,Vk,approxResults] = digitsComputeResults(A,M,labels,options)

if ~exist('options', 'var'), options = digitsExperimentParameters(); options = obj2struct(options); end

orthoFlag = {'orthoFlag',1};

% matrix
[u,s,v]               = svd(reshape(permute(A,[1,3,2]),[],size(A,2)),'econ');
[Uk{1},Sk{1},Vk{1}]   = deal(u(:,1:options.k), s(1:options.k,1:options.k), v(:,1:options.k));
Ak{1}                 = Uk{1} * Sk{1} * Vk{1}';
Ak{1}                 = permute(reshape(Ak{1},size(A,1),[],size(A,2)),[1,3,2]);
[err,errPerClass,res] = computeMetrics(A,Ak{1},labels,options);
approxResults{1}      = struct('err', err, 'res', res, 'errPerClass',errPerClass);

% matrix 2
[u,s,v]               = svd(reshape(A,[],size(A,3)),'econ');
[Uk{2},Sk{2},Vk{2}]   = deal(u(:,1:options.k), s(1:options.k,1:options.k), v(:,1:options.k));
Ak{2}                 = Uk{2} * Sk{2} * Vk{2}';
Ak{2}                 = reshape(Ak{2},size(A,1),[],size(A,3));
[err,errPerClass,res] = computeMetrics(A,Ak{2},labels,options);
approxResults{2}      = struct('err', err, 'res', res, 'errPerClass',errPerClass);


% t-product
[Uk{3},Sk{3},Vk{3}]   = tsvd(A,options.k);
Ak{3}                 = tprod(Uk{3},tprod(Sk{3},ttran(Vk{3})));
[err,errPerClass,res] = computeMetrics(A,Ak{3},labels,options);
approxResults{3}      = struct('err', err, 'res', res, 'errPerClass',errPerClass);
    
% remaining choices of M
for i = 4:length(M)
    [Uk{i},Sk{i},Vk{i}]   = tSVDM(A,M{i},options.k,orthoFlag{:});
    Ak{i}                 = mprod(Uk{i},mprod(Sk{i},tran(Vk{i}),M{i},orthoFlag{:}),M{i},orthoFlag{:});
    [err,errPerClass,res] = computeMetrics(A,Ak{i},labels,options);
    approxResults{i}      = struct('err', err, 'res', res, 'errPerClass',errPerClass);
end


end



function[err,errPerClass,res] = computeMetrics(A,Ak,labels,options)

if options.permuteFlag
    A  = permute(A,[1,3,2]);
    Ak = permute(Ak,[1,3,2]);
end

err = fronorm(A - Ak);

nClasses    = length(unique(labels));
errPerClass = zeros(1,nClasses);
for i = 0:nClasses-1
    idx = (labels == i);
    errPerClass(i + 1) = fronorm(A(:,:,idx) - Ak(:,:,idx));
end

res = reshape(sqrt(sum((A - Ak).^2,[1,2])),1,[]);

end


