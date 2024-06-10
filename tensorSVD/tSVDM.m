function[U,S,V,JacU,JacS,JacV,JacMU,JacMS,JacMV] = tSVDM(A,M,k,varargin)
% Compute the t-SVDM using the starM-product for third-order tensors
%
% Inputs:
%   A: third-order tensor of size n1 x n2 x n3
%   M: cell array of function handles with {forward, inverse} or a matrix
%   k (optional): integer, truncation parameter 
%   varargin (optional): 
%
% Outputs:
%   U: n1 x k x n3 tensor containing basis
%   S: k x k x n3 tensor containing singular tubes
%   V: n2 x k x n3 tensor containing normalized coefficients
%

if nargin == 0, runMinimalExample; return; end

doJacUSV = (nargout > 3);   JacU  = []; JacS  = []; JacV  = [];   
doJacM   = (nargout > 6);   JacMU = []; JacMS = []; JacMV = [];          

for j = 1:2:length(varargin), eval([varargin{j},'= varargin{j+1};']); end

% set default
if ~exist('M','var') || isempty(M), M = 1; end  % identity
if ~exist('k','var'), k = min(size(A,1),size(A,2)); end

% apply transformation to tubes
if doJacM
    [A,JacA1,JacMA1] = modeProduct(A,M);
elseif doJacUSV
    [A,JacA1] = modeProduct(A,M);
else
    A = modeProduct(A,M);
end

if doJacUSV || doJacM
    [U,S,V] = facewiseSVD(A);
    [JacU2,JacS2,JacV2] = facewiseSVDJacobian(U,S,V);

    % TODO: don't truncate too early
    U = U(:,1:k,:);
    S = S(1:k,1:k,:);
    V = V(:,1:k,:);

else
    [U,S,V] = facewiseSVD(A,k);
end

% apply inverse transformation
if doJacM
    [U,JacU3,JacMU3] = modeProduct(U,M,'invFlag',1,varargin{:});
    [S,JacS3,JacMS3] = modeProduct(S,M,'invFlag',1,varargin{:});
    [V,JacV3,JacMV3] = modeProduct(V,M,'invFlag',1,varargin{:});
elseif doJacUSV
    [U,JacU3] = modeProduct(U,M,'invFlag',1,varargin{:});
    [S,JacS3] = modeProduct(S,M,'invFlag',1,varargin{:});
    [V,JacV3] = modeProduct(V,M,'invFlag',1,varargin{:});
else
    U = modeProduct(U,M,'invFlag',1,varargin{:});
    S = modeProduct(S,M,'invFlag',1,varargin{:});
    V = modeProduct(V,M,'invFlag',1,varargin{:});
end

if doJacUSV
    r = min(size(A,1),size(A,2));
    zeroColPad = @(y) cat(2,y,zeros(size(y,1),r-size(y,2),size(A,3)));
    zeroColCut = @(x) x(:,1:k,:);
    zeroPad = @(y) cat(1,zeroColPad(y),zeros(r-size(y,1),r,size(A,3)));
    zeroCut = @(x) x(1:k,1:k,:);

    JacU.inputSize  = size(A);
    JacU.outputSize = size(U);
    JacU.A          = @(x) JacU3.A(zeroColCut(JacU2.A(JacA1.A(x))));
    JacU.AT         = @(y) JacA1.AT(JacU2.AT(zeroColPad(JacU3.AT(y))));

    JacS.inputSize  = size(A);
    JacS.outputSize = size(S);
    JacS.A          = @(x) JacS3.A(zeroCut(JacS2.A(JacA1.A(x))));
    JacS.AT         = @(y) JacA1.AT(JacS2.AT(zeroPad(JacS3.AT(y))));

    JacV.inputSize  = size(A);
    JacV.outputSize = size(V);
    JacV.A          = @(x) JacV3.A(zeroColCut(JacV2.A(JacA1.A(x))));
    JacV.AT         = @(y) JacA1.AT(JacV2.AT(zeroColPad(JacV3.AT(y))));
end

if doJacM
    JacMU.inputSize  = size(M);
    JacMU.outputSize = size(U);
    JacMU.A          = @(x) JacMU3.A(x) + JacU3.A(zeroColCut(JacU2.A(JacMA1.A(x))));
    JacMU.AT         = @(y) JacMU3.AT(y) + JacMA1.AT(JacU2.AT(zeroColPad(JacU3.AT(y))));

    JacMS.inputSize  = size(M);
    JacMS.outputSize = size(S);
    JacMS.A          = @(x) JacMS3.A(x) + JacS3.A(zeroCut(JacS2.A(JacMA1.A(x))));
    JacMS.AT         = @(y) JacMS3.AT(y) + JacMA1.AT(JacS2.AT(zeroPad(JacS3.AT(y))));

    JacMV.inputSize  = size(M);
    JacMV.outputSize = size(V);
    JacMV.A          = @(x) JacMV3.A(x) + JacV3.A(zeroColCut(JacV2.A(JacMA1.A(x))));
    JacMV.AT         = @(y) JacMV3.AT(y) + JacMA1.AT(JacV2.AT(zeroColPad(JacV3.AT(y))));
end


end


function[] = runMinimalExample()
A = randn(3,4,5);
M = 1;
[U,S,V] = tSVDM(A,M);
[UF,SF,VF] = facewiseSVD(A);

fprintf('||U - UF||_F = %0.4e\n',fronorm(U - UF));
fprintf('||S - SF||_F = %0.4e\n',fronorm(S - SF));
fprintf('||V - VF||_F = %0.4e\n',fronorm(V - VF));

M = randn(5);

% check if M is inverible 
if cond(M) > 1e12, error('M close to singular, bad initialization, re-run test'); end


disp('FULL RANK')
[U,S,V] = tSVDM(A,M);

disp('size(A) = '); disp(size(A));
disp('size(U) = '); disp(size(U));
disp('size(S) = '); disp(size(S));
disp('size(V) = '); disp(size(V));


disp('ECONOMIC')
[U,S,V] = tSVDM(A,M,"econ");

disp('size(A) = '); disp(size(A));
disp('size(U) = '); disp(size(U));
disp('size(S) = '); disp(size(S));
disp('size(V) = '); disp(size(V));

disp('LOW RANK')
[U,S,V] = tSVDM(A,M,1);

disp('size(A) = '); disp(size(A));
disp('size(U) = '); disp(size(U));
disp('size(S) = '); disp(size(S));
disp('size(V) = '); disp(size(V));


disp('JACOBIAN')
A = randn(3,4,5);
M = randn(5);
k = 1;
[U,S,V,JacU,JacS,JacV,JacMU,JacMS,JacMV] = tSVDM(A,M,k);


disp('JacU')
function[u] = fctnU(x,M,k)
    [u,~,~] = tSVDM(x,M,k);
end
fctn = @(x) fctnU(x,M,k);
CheckAdjoint(JacU);
CheckJacobian(fctn,JacU,[],0,A)
CheckGradient(fctn,JacU,[],0,A)

disp('JacS')
function[s] = fctnS(x,M,k)
    [~,s,~] = tSVDM(x,M,k);
end
fctn = @(x) fctnS(x,M,k);
CheckAdjoint(JacS);
CheckJacobian(fctn,JacS,[],0,A)
CheckGradient(fctn,JacS,[],0,A)


disp('JacV')
function[v] = fctnV(x,M,k)
    [~,~,v] = tSVDM(x,M,k);
end
fctn = @(x) fctnV(x,M,k);
CheckAdjoint(JacV);
CheckJacobian(fctn,JacV,[],0,A)
CheckGradient(fctn,JacV,[],0,A)


A = randn(20,10,7);
M = randn(7);
k = 1;
[U,S,V,JacU,JacS,JacV,JacMU,JacMS,JacMV] = tSVDM(A,M,k);


disp('JacMU')
fctn = @(m) fctnU(A,m,k);
CheckAdjoint(JacMU);
CheckJacobian(fctn,JacMU,[],0,M)
CheckGradient(fctn,JacMU,[],0,M)

disp('JacMS')
fctn = @(m) fctnS(A,m,k);
CheckAdjoint(JacMS);
CheckJacobian(fctn,JacMS,[],0,M)
CheckGradient(fctn,JacMS,[],0,M)

disp('JacMV')
fctn = @(m) fctnV(A,m,k);
CheckAdjoint(JacMV);
CheckJacobian(fctn,JacMV,[],0,M)
CheckGradient(fctn,JacMV,[],0,M)




disp("PASSED!")

end


