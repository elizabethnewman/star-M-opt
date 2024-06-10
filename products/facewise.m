function[C,JacA,JacB] = facewise(A,B)
% facewise product
%   Multiply the frontal slices of tensors
%
% Inputs:
%   A  : n1 x p x n3 x ... x n_{d} array
%   B  : p x n2 x n3 x ... x n_{d} array
%
% Outputs:
%   C    : n1 x n2 x ... x n_{k-1} x p x n_{k+1} x ... x n_{d} array
%   JacA : Jacobian of C w.r.t. A
%   JacB : Jacobian of C w.r.t. B
%

if nargin == 0, runMinimalExample; return; end

JacA = [];   doJacA = (nargout > 1); 
JacB = [];   doJacB = (nargout > 2); 

% multiply
C = pagemtimes(A,B);

if doJacA
    JacA.inputSize  = size(A);
    JacA.outputSize = size(C);
    JacA.A  = @(x) pagemtimes(x,B); 
    JacA.AT = @(y) pagemtimes(y,tran(B));
end

if doJacB
    JacB.inputSize  = size(B);
    JacB.outputSize = size(C);
    JacB.A  = @(x) pagemtimes(A,x); 
    JacB.AT = @(y) pagemtimes(tran(A),y);
end

end



function[] = runMinimalExample()
    A = reshape(1:(5*4*3*2),[5,4,3,2]);
    B = reshape(1:(4*7*3*2),[4,7,3,2]);
    [C,JacA,JacB] = facewise(A,B);
    
    disp('size(A) = '); disp(size(A));
    disp('size(B) = '); disp(size(B));
    disp('size(C) = '); disp(size(C));
    
    
    disp('test JacA')
    fctn = @(a) facewise(a,B);
    CheckAdjoint(JacA);
    CheckJacobian(fctn,JacA,[],0)
    CheckGradient(fctn,JacA,[],0)
    disp('PASSED!')
    disp(' ')
    
    
    disp('test JacB')
    fctn = @(b) facewise(A,b);
    CheckAdjoint(JacB);
    CheckJacobian(fctn,JacB,[],0)
    CheckGradient(fctn,JacB,[],0)
    disp('PASSED!')
    disp(' ')

end
