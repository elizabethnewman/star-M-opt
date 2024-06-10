function[JacU,JacS,JacV,JacM] = tSVDMJacobian(U,S,V,M,transformFlag,JacUM1,JacUM2,Jac)
% https://j-towns.github.io/papers/svd-derivative.pdf

% assume matrices are full and full-rank
% create functions (x is the direction in which we apply derivative)
% assume x is the size of A, m x n


if nargin == 0, runMinimalExample; return; end


if ~exist('transformFlag','var') || isempty(transformFlag), transformFlag = 1; end

if transformFlag
    [U,JacUM1] = modeProduct(U,M);
    [S,JacSM1] = modeProduct(S,M);
    [V,JacVM1] = modeProduct(V,M);
end

[JacU2,JacS2,JacV2] = facewiseSVDJacobian(U,S,V);

% 

end


function runMinimalExample()

% ------------------------------------------ %
disp('UNDERDETERMINED')
A = randn(3,10);
[U,S,V] = svd(A,"econ");

[JacU,JacS,JacV] = facewiseSVDJacobian(U,S,V);

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
% 

% % ------------------------------------------ %
disp('OVERDETERMINED')
A = randn(20,5);
[U,S,V] = svd(A,"econ");

[JacU,JacS,JacV] = facewiseSVDJacobian(U,S,V);

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



end







