function[JacU,JacS,JacV] = svdJacobian(U,S,V)
% https://j-towns.github.io/papers/svd-derivative.pdf

% assume matrices are full and full-rank
% create functions (x is the direction in which we apply derivative)
% assume x is the size of A, m x n


if nargin == 0, runMinimalExample; return; end

% get sizes
m  = size(U,1);
n  = size(V,1);

overdeterminedFlag = (m > n);
if overdeterminedFlag
    mOld = m;
    m    = n;
    n    = mOld;

    UOld = U;
    U    = V;
    V    = UOld;
end


% Note: think about how to handle situations with equal singular values
s2     = diag(S).^2;
F      = -s2 + s2';
idx    = (F ~= 0);  % avoid small singular values
F(idx) = 1 ./ F(idx);

% get helper functions
dP = @(x) U' * x * V;
dC = @(x) F .* (dP(x) * S + S * dP(x)');
dD = @(x) F .* (S * dP(x) + dP(x)' * S);

% get jacobians
JacU.inputSize  = [m,n];
JacU.outputSize = size(U); 
JacU.A  = @(x) U * dC(x) + (eye(m) - U * U') * x * (V / S);
JacU.AT = @(y) (U * (F .* (U' * y - y' * U)) * S + (eye(m) - U * U') * (y / S)) * V';

JacS.inputSize  = [m,n];
JacS.outputSize = size(S); 
JacS.A  = @(x) eye(numel(s2)) .* dP(x);
JacS.AT = @(y) U * (eye(numel(s2)) .* y) * V';

JacV.inputSize  = [m,n];
JacV.outputSize = size(V); 
JacV.A  = @(x) V * dD(x) + (eye(n) - V * V') * x' * (U / S);
JacV.AT = @(y) U * (S * (F .* (V' * y - y' * V)) * V' + (S \ y') * (eye(n) - V * V'));


if overdeterminedFlag
    mOld = m;
    m    = n;
    n    = mOld;
    
    JacUOld = JacU;
    
    JacU.inputSize  = [m,n];
    JacU.outputSize = [m,n];
    JacU.A  = @(x) JacV.A(x');
    JacU.AT = @(x) JacV.AT(x)';
    
    JacS.inputSize  = [m,n];
    JacS.outputSize = [n,n];
    JacS.A  = @(x) JacS.A(x');
    JacS.AT = @(x) JacS.AT(x)';
    

    JacV.inputSize  = [m,n];
    JacV.outputSize = [n,n];
    JacV.A  = @(x) JacUOld.A(x');
    JacV.AT = @(x) JacUOld.AT(x)';
end


end


function runMinimalExample()

% ------------------------------------------ %
disp('UNDERDETERMINED')
A = randn(3,10);
[U,S,V] = svd(A,"econ");

[JacU,JacS,JacV] = svdJacobian(U,S,V);

disp('JacU')
% create function
function[u] = fctnU1(x)
    [u,~,~] = svd(x,"econ");
end
fctn = @(x) fctnU1(x);
CheckAdjoint(JacU);
CheckJacobian(fctn,JacU,[],0,A)
CheckGradient(fctn,JacU,[],0,A)

disp('JacS')
function[s] = fctnS1(x)
    [~,s,~] = svd(x,"econ");
end
fctn = @(x) fctnS1(x);
CheckAdjoint(JacS);
CheckJacobian(fctn,JacS,[],0,A)
CheckGradient(fctn,JacS,[],0,A)

disp('JacV')
function[v] = fctnV1(x)
    [~,~,v] = svd(x,"econ");
end
fctn = @(x) fctnV1(x);
CheckAdjoint(JacV);
CheckJacobian(fctn,JacV,[],0,A)
CheckGradient(fctn,JacV,[],0,A)

disp('PASSED!')


% ------------------------------------------ %
disp('OVERDETERMINED')
A = randn(20,5);
[U,S,V] = svd(A,"econ");

[JacU,JacS,JacV] = svdJacobian(U,S,V);

disp('JacU')
% create function
function[u] = fctnU3(x)
    [u,~,~] = svd(x,"econ");
end
fctn = @(x) fctnU3(x);
CheckAdjoint(JacU);
CheckJacobian(fctn,JacU,[],0,A)
CheckGradient(fctn,JacU,[],0,A)

disp('JacS')
function[s] = fctnS3(x)
    [~,s,~] = svd(x,"econ");
end
fctn = @(x) fctnS3(x);
CheckAdjoint(JacS);
CheckJacobian(fctn,JacS,[],0,A)
CheckGradient(fctn,JacS,[],0,A)

disp('JacV')
function[v] = fctnV3(x)
    [~,~,v] = svd(x,"econ");
end
fctn = @(x) fctnV3(x);
CheckAdjoint(JacV);
CheckJacobian(fctn,JacV,[],0,A)
CheckGradient(fctn,JacV,[],0,A)

disp('PASSED!')


% ------------------------------------------ %
disp('SQUARE')
A = randn(5,5);
[U,S,V] = svd(A,"econ");

[JacU,JacS,JacV] = svdJacobian(U,S,V);

disp('JacU')
% create function
function[u] = fctnU4(x)
    [u,~,~] = svd(x,"econ");
end
fctn = @(x) fctnU4(x);
CheckAdjoint(JacU);
CheckJacobian(fctn,JacU,[],0,A)
CheckGradient(fctn,JacU,[],0,A)

disp('JacS')
function[s] = fctnS4(x)
    [~,s,~] = svd(x,"econ");
end
fctn = @(x) fctnS4(x);
CheckAdjoint(JacS);
CheckJacobian(fctn,JacS,[],0,A)
CheckGradient(fctn,JacS,[],0,A)

disp('JacV')
function[v] = fctnV4(x)
    [~,~,v] = svd(x,"econ");
end
fctn = @(x) fctnV4(x);
CheckAdjoint(JacV);
CheckJacobian(fctn,JacV,[],0,A)
CheckGradient(fctn,JacV,[],0,A)

disp('PASSED!')




end







