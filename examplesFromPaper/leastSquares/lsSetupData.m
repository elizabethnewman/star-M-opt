function[A,B] = lsSetupData(options)

if ~exist('options', 'var'), options = lsExperimentParameters(); options = obj2struct(options); end

rng(options.seed)

% create tensors in transform domain
n3   = numel(options.m);        % number of frontal slices
AHat = ones(options.nPoints,2,n3);
BHat = zeros(options.nPoints,1,n3);

for i = 1:n3
    x           = randn(options.nPoints,1);
    AHat(:,2,i) = x;
    BHat(:,1,i) = options.m(i) * x + options.b(i);
end

% assume M is orthogonal, so M^{-1} = M^T
A = modeProduct(AHat, options.M');
B = modeProduct(BHat, options.M');
A = A + options.eta * randn(size(A));
B = B + options.eta * randn(size(B));

% A = modeProduct(AHat + (options.eta / fronorm(AHat)) * randn(size(AHat)), options.M');
% B = modeProduct(BHat + (options.eta / fronorm(BHat)) * randn(size(BHat)), options.M');


end

