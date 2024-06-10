function[JacU,JacS,JacV] = facewiseSVDJacobian(U,S,V)
% https://j-towns.github.io/papers/svd-derivative.pdf

% assume matrices are full and full-rank
% create functions (x is the direction in which we apply derivative)
% assume x is the size of A, m x n


if nargin == 0, runMinimalExample; return; end

% get sizes
m  = size(U,1);
n  = size(V,1);
szA = size(U);

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
s      = facewiseDiag(S);
szS    = [1,m,szA(3:end)];
s2     = s.^2;
F      = -s2 + tran(s2);
idx    = (F ~= 0);  % avoid small singular values
F(idx) = 1 ./ F(idx);

% get helper functions
dP = @(x) facewise(tran(U),facewise(x,V));
dC = @(x) F .* (facewise(dP(x),S) + facewise(S,tran(dP(x))));
dD = @(x) F .* (facewise(S,dP(x)) + facewise(tran(dP(x)),S));


% get jacobians
UUT = facewise(U,tran(U));
JacU.inputSize  = [m,n,szA(3:end)];
JacU.outputSize = size(U); 
JacU.A  = @(x) facewise(U,dC(x)) ...
    + facewise(facewise(eye(m) - UUT,x),V ./ reshape(s,szS));

tmp1 = @(y) facewise(U,facewise(F .* (facewise(tran(U),y) - facewise(tran(y),U)),S));
tmp2 = @(y) facewise(eye(m) - UUT,y ./ reshape(s,szS));
JacU.AT = @(y) facewise(tmp1(y) + tmp2(y),tran(V));

JacS.inputSize  = [m,n,szA(3:end)];
JacS.outputSize = size(S); 
JacS.A  = @(x) eye(size(s2,1)) .* dP(x);
JacS.AT = @(y) facewise(facewise(U,(eye(m) .* y)),tran(V));

VVT = facewise(V,tran(V));
JacV.inputSize  = [m,n,szA(3:end)];
JacV.outputSize = size(V); 

JacV.A  = @(x) facewise(V,dD(x)) + facewise(eye(n) - VVT,facewise(tran(x),U ./ reshape(s,szS)));

tmp1 = @(y) facewise(facewise(S,(F .* (facewise(tran(V),y) - facewise(tran(y),V)))),tran(V));
tmp2 = @(y) facewise(tran(y ./ reshape(s,[1,m,szA(3:end)])),eye(n) - VVT);
JacV.AT = @(y) facewise(U,tmp1(y) + tmp2(y));


if overdeterminedFlag
    mOld = m;
    m    = n;
    n    = mOld;
    
    JacUOld = JacU;
    
    JacU.inputSize  = [m,n,szA(3:end)];
    JacU.outputSize = [m,n,szA(3:end)];
    JacU.A  = @(x) JacV.A(tran(x));
    JacU.AT = @(x) tran(JacV.AT(x));
    
    JacS.inputSize  = [m,n,szA(3:end)];
    JacS.outputSize = [n,n,szA(3:end)];
    JacS.A  = @(x) JacS.A(tran(x));
    JacS.AT = @(x) tran(JacS.AT(x));
    

    JacV.inputSize  = [m,n,szA(3:end)];
    JacV.outputSize = [n,n,szA(3:end)];
    JacV.A  = @(x) JacUOld.A(tran(x));
    JacV.AT = @(x) tran(JacUOld.AT(x));
end


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







