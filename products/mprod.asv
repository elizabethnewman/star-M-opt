function[C,JacA,JacB,JacM] = mprod(A,B,M,varargin)
% starM product for third-order tensors 
%
% Inputs:
%   A  : n1 x p x n3 array
%   B  : p x n2 x n3 array
%   M  : n3 x n3 invertible matrix
%   varargin : additional flags for modeProduct
%       'transFlag'     : apply M' along mode-k (boolean, default = false)
%       'orthoFlag'     : apply M' along mode-k (boolean, default = false)
%       'invFlag'       : apply the inverse of M along mode k (boolean, default = false)
%
%   NOTE: 'orthoFlag' supersedes 'invFlag' (use cautiously with derivative checks)
%
% Outputs:
%   C    : A \starM B, n1 x n2 n3 array
%   JacA : Jacobian of C w.r.t. A
%   JacB : Jacobian of C w.r.t. B
%   JacM : Jacobian of C w.r.t. M
%

if nargin == 0, runMinimalExample; return; end

JacA = [];   doJacA = (nargout > 1); 
JacB = [];   doJacB = (nargout > 2); 
JacM = [];   doJacM = (nargout > 3); 

% ------------------------------------------- %
% move to transform domain (PART 1)
if doJacM
    [A,JacA1,JacM1A] = modeProduct(A,M);
elseif doJacA
    [A,JacA1] = modeProduct(A,M);
else
    A = modeProduct(A,M);
end

if doJacM 
    [B,JacB1,JacM1B] = modeProduct(B,M);
elseif doJacB
    [B,JacB1] = modeProduct(B,M);
else
    B = modeProduct(B,M);
end

% ------------------------------------------- %
% facewise multiply (PART 2)
if doJacM || doJacB
    [C,JacA2,JacB2] = facewise(A,B);
elseif doJacA
    [C,JacA2] = facewise(A,B);
else
    C = facewise(A,B);
end

% ------------------------------------------- %
% return to spatial domain (PART 3)
if doJacM
    [C,JacC3,JacM3] = modeProduct(C,M,'invFlag',1,varargin{:});
elseif doJacB || doJacA
    [C,JacC3] = modeProduct(C,M,'invFlag',1,varargin{:});
else
    C = modeProduct(C,M,'invFlag',1,varargin{:});
end

% ------------------------------------------- %
% put Jacobians together
if doJacA
    JacA.inputSize  = size(A);
    JacA.outputSize = size(C);
    JacA.A  = @(x) JacC3.A(JacA2.A(JacA1.A(x))); 
    JacA.AT = @(y) JacA1.AT(JacA2.AT(JacC3.AT(y)));
end

if doJacB
    JacB.inputSize  = size(B);
    JacB.outputSize = size(C);
    JacB.A  = @(x) JacC3.A(JacB2.A(JacB1.A(x))); 
    JacB.AT = @(y) JacB1.AT(JacB2.AT(JacC3.AT(y)));
end

if doJacM
    JacM.inputSize  = size(M);
    JacM.outputSize = size(C);
    JacM.A  = @(x) JacM3.A(x) + JacC3.A(JacA2.A(JacM1A.A(x))) + JacC3.A(JacB2.A(JacM1B.A(x)));
    JacM.AT = @(y) JacM3.AT(y) + JacM1A.AT(JacA2.AT(JacC3.AT(y))) + JacM1B.AT(JacB2.AT(JacC3.AT(y)));
end

end



function[] = runMinimalExample()
A = randn(3,4,5);
B = randn(4,7,5);

M = eye(5);
tol = 1e-14;
assert(fronorm(facewise(A,B) - mprod(A,B,M)) < tol, 'check implementation')


M = randn(5);

% check if M is inverible 
if cond(M) > 1e12, error('M close to singular, bad initialization, re-run test'); end

C = mprod(A,B,M);

disp('size(A) = ')
disp(size(A));

disp('size(B) = ')
disp(size(B));

disp('size(C) = ')
disp(size(C));


disp('test JacA')
[C,JacA] = mprod(A,B,M);
fctn = @(a) mprod(a,B,M);
CheckAdjoint(JacA);
CheckJacobian(fctn,JacA,[],0)
CheckGradient(fctn,JacA,[],0)
disp('PASSED!')
disp(' ')

disp('test JacB')
[C,JacA,JacB] = mprod(A,B,M);
fctn = @(b) mprod(A,b,M);
CheckAdjoint(JacB);
CheckJacobian(fctn,JacB,[],0)
CheckGradient(fctn,JacB,[],0)
disp('PASSED!')
disp(' ')

disp('test JacM')
[C,JacA,JacB,JacM] = mprod(A,B,M);
fctn = @(m) mprod(A,B,m);
CheckAdjoint(JacM);

% need to give M as direction (special case because using implicit differentiation)
CheckJacobian(fctn,JacM,[],0,M)
CheckGradient(fctn,JacM,[],0,M)
disp('PASSED!')
disp(' ')


end
