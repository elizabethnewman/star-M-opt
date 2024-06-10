function[B,JacA,JacM] = modeProduct(A,M,varargin)
% Mode-k product
%   Apply a matrix M along the mode-k fibers of a tensor A
%       A x_k M
%
% Inputs:
%   A  : n1 x n2 x ... x n_{k-1} x n_{k} x n_{k+1} x ... x n_{d} array
%   M  : p x n_k matrix
%   k  : dimension along which to unfold (optional, default k = d)
%   varargin : additional flags
%       'transFlag'     : apply M' along mode-k (boolean, default = false)
%       'orthoFlag'     : apply M' along mode-k (boolean, default = false)
%       'invFlag'       : apply the inverse of M along mode k (boolean, default = false)
%
%   NOTE: 'orthoFlag' supersedes 'invFlag' (use cautiously with derivative checks)
%
% Outputs:
%   B    : A x_k M, n1 x n2 x ... x n_{k-1} x p x n_{k+1} x ... x n_{d} array
%   JacA : Jacobian of B w.r.t. A
%   JacM : Jacobian of B w.r.t. M
%

if nargin == 0, runMinimalExample; return; end


% parse inputs
for j = 1:2:length(varargin), eval([varargin{j},'= varargin{j+1};']); end

if ~exist('dim','var') || isempty(dim), dim = max(ndims(A),3); end

% get output size (will change if M is not square)
sB = size(A);

% initialize Jacobians
JacA = [];    doJacA = (nargout > 1); 
JacM = [];    doJacM = (nargout > 2);

% compute product
if (exist('transFlag','var') && transFlag == 1) || (exist('orthoFlag','var') && orthoFlag == 1)

    if ~isscalar(M) && size(M,1) ~= size(A,dim), error('incompatible sizes'); end

    % change the size of the k-th dimension of B
    if ~isscalar(M), sB(dim) = size(M,2); end

    % apply the transpose of M
    B = modeFold(M' * modeUnfold(A,dim),sB,dim);

    if doJacA
        JacA.A  = @(x) modeFold(M' * modeUnfold(x,dim),sB,dim);
        JacA.AT = @(y) modeFold(M * modeUnfold(y,dim),size(A),dim);
    end

    if doJacM
        JacM.A  = @(x) modeFold(x' * modeUnfold(A,dim),sB,dim);
        JacM.AT = @(y) modeUnfold(A,dim) * modeUnfold(y,dim)';
    end
  
elseif exist('invFlag','var') && invFlag == 1
    % apply the inverse of M

    if size(M,1) ~= size(M,2), error('M is not invertible'); end
    if ~isscalar(M) && size(M,2) ~= size(A,dim), error('incompatible sizes'); end

    B = modeFold(M \ modeUnfold(A,dim),sB,dim);

    if doJacA
        JacA.A  = @(x) modeFold(M \ modeUnfold(x,dim),sB,dim);
        JacA.AT = @(y) modeFold(M' \ modeUnfold(y,dim),size(A),dim);
    end

    if doJacM
        JacM.A  = @(x) modeFold((-M \ (x / M)) * modeUnfold(A,dim),sB,dim);
        JacM.AT = @(y) -M' \ ((modeUnfold(y,dim) * modeUnfold(A,dim)') / M');
    end

else
    if ~isscalar(M) && size(M,2) ~= size(A,dim), error('incompatible sizes'); end
    
    % change the size of the k-th dimension of B
    if ~isscalar(M), sB(dim) = size(M,1); end

    % apply the transformation
    B = modeFold(M * modeUnfold(A,dim),sB,dim);

    if doJacA
        JacA.A  = @(x) modeFold(M * modeUnfold(x,dim),sB,dim);
        JacA.AT = @(y) modeFold(M' * modeUnfold(y,dim),size(A),dim);
    end

    if doJacM
        JacM.A  = @(x) modeFold(x * modeUnfold(A,dim),sB,dim);
        JacM.AT = @(y) modeUnfold(y,dim) * modeUnfold(A,dim)';
    end

end

% store sizes in Jacobians
if doJacA
    JacA.inputSize  = size(A);
    JacA.outputSize = sB;
end

if doJacM
    JacM.inputSize  = size(M);
    JacM.outputSize = sB;
end



end


function[] = runMinimalExample()
A  = randn(5,4,3,2);
M1 = randn(5,5);
M2 = randn(4,1);
M3 = randn(7,3);
M4 = randn(2,2);

% ------------------------------------------------------------- %
% standard product
disp('all flags = 0')
[B,JacA,JacM] = modeProduct(A,M3,'dim',3);

fctn = @(a) modeProduct(a,M3,'dim',3);
CheckAdjoint(JacA);
CheckJacobian(fctn,JacA,[],0)
CheckGradient(fctn,JacA,[],0)

fctn = @(m) modeProduct(A,m,'dim',3);
CheckAdjoint(JacM);
CheckJacobian(fctn,JacM,[],0)
CheckGradient(fctn,JacM,[],0)

[B,JacA,JacM] = modeProduct(A,M4);

fctn = @(a) modeProduct(a,M4);
CheckAdjoint(JacA);
CheckJacobian(fctn,JacA,[],0)
CheckGradient(fctn,JacA,[],0)

fctn = @(m) modeProduct(A,m);
CheckAdjoint(JacM);
CheckJacobian(fctn,JacM,[],0)
CheckGradient(fctn,JacM,[],0)
disp('PASSED!')
disp(' ')


% % ------------------------------------------------------------- %
% % check transpose
disp('transFlag = 1')
[B,JacA,JacM] = modeProduct(A,M2,'dim',2,'transFlag',1);

fctn = @(a) modeProduct(a,M2,'dim',2,'transFlag',1);
CheckAdjoint(JacA);
CheckJacobian(fctn,JacA,[],0)
CheckGradient(fctn,JacA,[],0)

fctn = @(m) modeProduct(A,m,'dim',2,'transFlag',1);
CheckAdjoint(JacM);
CheckJacobian(fctn,JacM,[],0)
CheckGradient(fctn,JacM,[],0)
disp('PASSED!')
disp(' ')

% % ------------------------------------------------------------- %
disp('orthoFlag = 1')
[B,JacA,JacM] = modeProduct(A,M4,'orthoFlag',1);

fctn = @(a) modeProduct(a,M4,'orthoFlag',1);
CheckAdjoint(JacA);
CheckJacobian(fctn,JacA,[],0)
CheckGradient(fctn,JacA,[],0)

fctn = @(m) modeProduct(A,m,'orthoFlag',1);
CheckAdjoint(JacM);
CheckJacobian(fctn,JacM,[],0)
CheckGradient(fctn,JacM,[],0)
disp('PASSED!')
disp(' ')


% ------------------------------------------------------------- %
% check inverse
disp('invFlag = 1')
[B,JacA,JacM] = modeProduct(A,M1,'dim',1,'invFlag',1);


fctn = @(a) modeProduct(a,M1,'dim',1,'invFlag',1);
CheckAdjoint(JacA);
CheckJacobian(fctn,JacA,[],0)
CheckGradient(fctn,JacA,[],0)

fctn = @(m) modeProduct(A,m,'dim',1,'invFlag',1);
CheckAdjoint(JacM);

% need to give M1 as direction (special case because using implicit differentiation)
CheckJacobian(fctn,JacM,[],0,M1)
CheckGradient(fctn,JacM,[],0,M1)

disp('PASSED!')
disp(' ')

end
