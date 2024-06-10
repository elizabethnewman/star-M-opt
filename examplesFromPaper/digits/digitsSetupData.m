function[A,labels,nrmAPerClass,idx] = digitsSetupData(options)

if ~exist('options', 'var'), options = digitsExperimentParameters(); options = obj2struct(options); end

% get necessary parameters
seed        = options.seed;
nImg        = options.nImg;
nClasses    = options.nClasses;

% for reproducibility
rng(seed);

% load data (rotated MNIST)
[XTrain,YTrain] = digitTrain4DArrayData;

% get true labels
classes = double(YTrain) - 1;

% use only 300 images for illustrative purposes, 30 per class
idx = [];
for i = 0:nClasses-1
    idxLabels = find(classes == i);
    idxPerm   = randperm(length(idxLabels),nImg);
    idx       = cat(2,idx,idxLabels(idxPerm));
end
A = double(squeeze(XTrain(:,:,:,idx)));
labels = kron(0:nClasses-1,ones(1,nImg));

% info about A
nrmA = fronorm(A);

nrmAPerClass = zeros(1,nClasses);
for i = 0:nClasses-1
    nrmAPerClass(i + 1) = fronorm(A(:,:,labels == i));
end

% normalize (so error = relative error)
A = A / nrmA;

% shuffle one more time
% idxShuffle = randperm(size(A,3));
% A          = A(:,:,idxShuffle);
% labels     = labels(idxShuffle);

% add noise
A = A + options.noiseLevel * randn(size(A));

% permute
if options.permuteFlag
    A = permute(A,[1,3,2]);
end

end


