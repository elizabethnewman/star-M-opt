function[X,approxResults] = lsComputeResults(A,B,M,options)

AMat = reshape(permute(A, [2, 1, 3]), 2, [])';
BMat = reshape(permute(B, [2, 1, 3]), 1, [])';
XMat = AMat \ BMat;
X{1} = XMat;
approxResults{1} = struct('err',fronorm(AMat * XMat - BMat),...
    'relErr',fronorm(AMat * XMat - BMat) / fronorm(BMat),...
    'res',sqrt(sum((AMat * XMat - BMat).^2,2)));

AHat = ttransform(A);
BHat = ttransform(B);
XHat = pagemldivide(AHat,BHat);
X{2} = ttransform(XHat,1);
approxResults{2} = struct('err',fronorm(tprod(A,X{2}) - B),...
    'relErr',fronorm(tprod(A,X{2}) - B) / fronorm(B),...
    'res',reshape(sqrt(sum((tprod(A,X{2}) - B).^2,2)),[],1));


for i = 3:length(M)
    [X{i},err,relErr,res] = computeMetrics(A,B,M{i});
    approxResults{i} = struct('err',err, 'relErr',relErr,'res',res);
end

end


function[X,err,relErr,res] = computeMetrics(A,B,M)

AHat = modeProduct(A,M);
BHat = modeProduct(B,M);
XHat = pagemldivide(AHat,BHat);
X    = modeProduct(XHat,M'); % assume M is orthogonal

err     = fronorm(mprod(A,X,M) - B);
relErr  = err / fronorm(B);
res     = reshape(sqrt(sum((mprod(A,X,M) - B).^2,2)),[],1);

end
