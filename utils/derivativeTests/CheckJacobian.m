function CheckJacobian(fctn,Jac,seed,verbose,x,p)
% CheckJacobian
%   Check the Jacobian using Taylor approximation
%
% Inputs: 
%   fctn : function handle that takes in inputSize and returns outputSize
%   Jac  : Jacobian struct
%   seed       : random seed (optional)
%   verbose : boolean for printouts, optional
%
% Outputs: 
%

if exist('seed','var') && ~isempty(seed), rng(seed); end
if ~exist('verbose','var'), verbose = 0; end
if ~exist('x','var') || isempty(x), x = randn(Jac.inputSize); end

% evaluate
f0   = fctn(x);
nrm0 = fronorm(f0);

% directional derivative
if ~exist('p','var') || isempty(p), p = randn(size(x)); end
dfdp = Jac.A(p);

% main iteration
N   = 15;
err = zeros(N,3);

base = 2;
for k = 1:N
    h = base^(-k);
    
    ft = fctn(x + h * p);
    
    err0     = fronorm(f0 - ft) / nrm0;
    err1     = fronorm(f0 + h * dfdp - ft) / nrm0;
    err(k,:) = [h,err0,err1];

    [c0,d0] = convert2base(err0,base);
    [c1,d1] = convert2base(err1,base);
        
    if verbose
        fprintf('h=%0.2f x 2^(%0.2d)\tE0=%0.2f x 2^(%0.2d)\tE1=%0.2f x 2^(%0.2d)\n',1,-k,c0,d0,c1,d1);
    end
end

% first order
tol      = 0.1;
err1     = err(:,3);
order1   = err(1:end-1,3) ./ err(2:end,3);
gradFlag = or(sum(order1 > base^2 - tol) > (N / 3), mean(err1(err1 ~= Inf & ~isnan(err1))) < 1e-8);
assert(gradFlag,'Jacobian FAILED.')
if verbose, disp('Jacobian PASSED!'); end

end






