function[Ak,U,S,V,approxResults] = romComputeResults(A,M,k,options)

orthoFlag = {'orthoFlag',1};

Ak  = cell(1,length(M));
U   = cell(1,length(M));
S   = cell(1,length(M));
V   = cell(1,length(M));

approxResults = cell(1,length(M));

% matrix
% r = min(k*n3,min(size(A(:,:),1), size(A(:,:),2)));

[U{1},S{1},V{1}] = svd(A(:,:),'econ');
Ak{1} = reshape(U{1}(:,1:k) * S{1}(1:k,1:k) * V{1}(:,1:k)',size(A));

% t-product
[U{2},S{2},V{2}] = tsvd(A);
Ak{2} = tprod(U{2}(:,1:k,:),tprod(S{2}(1:k,1:k,:), ttran(V{2}(:,1:k,:))));

% other choices of M
for i = 3:length(M)
    [U{i},S{i},V{i}] = tSVDM(A,M{i},[],orthoFlag{:});

    Ak{i} = mprod(U{i}(:,1:k,:),mprod(S{i}(1:k,1:k,:),tran(V{i}(:,1:k,:)),M{i},orthoFlag{:}),M{i},orthoFlag{:});
    
end


for i = 1:length(M)
    [err,errPerParam] = computeMetrics(A,Ak{i});
    approxResults{i}.err         = err;
    approxResults{i}.errPerParam = errPerParam;
end



end



function[err,errPerParam] = computeMetrics(A,Ak)

err         = fronorm(A - Ak);
errPerParam = reshape(sqrt(sum((A - Ak).^2, [1, 2])),1,[]);
end

