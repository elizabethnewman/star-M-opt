function CheckGradient(fctn,Jac,seed,verbose,x,p)
% CheckGradient
%   Check the gradient (Jacobian transpose) in forward mode using Taylor approximation
%
% Inputs: 
%   fctn    : function handle that takes in inputSize and returns outputSize
%   Jac     : Jacobian struct
%   seed    : random seed (optional)
%   verbose : boolean for printouts, (optional)
%
% Outputs: 
%

if exist('seed','var') && ~isempty(seed), rng(seed); end
if ~exist('verbose','var'), verbose = 0; end
if ~exist('x','var') || isempty(x), x = randn(Jac.inputSize); end

% evaluate function
y = fctn(x);

% temporary test objective function
c    = randn(size(y));
f0   = 0.5 * fronorm(y - c).^2;
nrm0 = fronorm(f0);

% directional derivative
if ~exist('p','var') || isempty(p), p = randn(size(x)); end
df0  = Jac.AT(y - c);
dfdp = sum(df0 .* p,'all');

% main iteration
N   = 15;
err = zeros(N,3);

base = 2;
for k = 1:N
    h = base^(-k);
    
    yt = fctn(x + h * p);
    ft = 0.5 * fronorm(yt - c).^2;
    
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
assert(gradFlag,'Gradient FAILED.')
if verbose, disp('Gradient PASSED!'); end

end


