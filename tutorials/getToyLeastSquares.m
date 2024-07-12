function[A, B] = getToyLeastSquares(nPoints,n3,noiseLevel,M)
    defaultM = (nargin < 4);

    m = linspace(-1, 1, n3);
    b = linspace(-1, 1, n3);
    
    % generate points in the transform domain
    AHat = ones(nPoints, 2, n3);
    BHat = zeros(nPoints, 1, n3);
    for i = 1:n3
        xx = randn(nPoints, 1);
        AHat(:,2,i) = xx;
        BHat(:,1,i) = m(i) * xx + b(i) + noiseLevel * randn(nPoints, 1);
    end

    % if no transformation is specified by user, use a default orthogonal
    % matrix
    if defaultM
        M = 1 / sqrt(2) * [1, 1; 1, -1];
    end

    % move to spatial domain
    A = modeProduct(AHat, M, 'invFlag', 1) + noiseLevel * randn(size(AHat));
    B = modeProduct(BHat, M, 'invFlag', 1) + noiseLevel * randn(size(BHat));

end