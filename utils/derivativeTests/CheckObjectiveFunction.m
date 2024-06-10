function CheckObjectiveFunction(fctn,x0,varargin)
% CheckObjectiveFunction
%   Check the gradient of objective function
%
% Inputs: 
%   fctn    : function handle that takes in inputSize and returns outputSize
%   Jac     : Jacobian struct
%   seed    : random seed (optional)
%   verbose : boolean for printouts, (optional)
%
% Outputs: 
%

for k = 1:2:length(varargin)
    eval([varargin{k}, ' = varargin{k + 1};'])
end

if ~exist('verbose','var'), verbose = 0; end
if exist('seed','var'), rng(seed); end

% evaluate function
[f0,~,df0] = fctn.evaluate(x0);
nrm0 = fronorm(f0);

% choose perturbation
p    = randn(size(x0));
dfdp = sum(df0 .* p,'all');

% main iteration
N   = 15;
err = zeros(N,3);

base = 2;
for k = 1:N
    h = base^(-k);

    ft = fctn.evaluate(x0 + h * p);
    
    err0     = fronorm(f0 - ft) / nrm0;
    err1     = fronorm(f0 + h * dfdp - ft) / nrm0;
    err(k,:) = [h,err0,err1];
        
    if verbose
        [c0,d0] = convert2base(err0,base);
        [c1,d1] = convert2base(err1,base);
        fprintf('h=%0.2f x 2^(%0.2d)\tE0=%0.2f x 2^(%0.2d)\tE1=%0.2f x 2^(%0.2d)\n',1,-k,c0,d0,c1,d1);
    end
end

% first order
tol      = 0.1;
err1     = err(:,3);
order1   = err(1:end-1,3) ./ err(2:end,3);
gradFlag = or(sum(order1 > base^2 - tol) > (N / 3), mean(err1(err1 ~= Inf & ~isnan(err1))) < 1e-8);

assert(gradFlag,'Gradient FAILED.')

if verbose, disp('Gradient PASSED!'); end

end


