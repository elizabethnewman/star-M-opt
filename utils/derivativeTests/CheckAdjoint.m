function CheckAdjoint(Jac,seed,verbose)
% CheckAdjoint
%   Check the adjoint for the Jacobian
%       <Jac.A(v),w> = <v,Jac.AT(w)>
%
% Inputs: 
%   Jac  : Jacobian struct
%   seed : random seed (optional)
%
% Outputs: 
%

if exist('seed','var') && ~isempty(seed), rng(seed); end
if ~exist('verbose','var'), verbose = 0; end

% directions to apply Jacobian
v = randn(Jac.inputSize);
w = randn(Jac.outputSize);

% compute <Jac.A(v),w>
Av   = Jac.A(v);
w_Av = sum(w .* Av,'all');

% compute <v,Jac.AT(w)>
ATw = Jac.AT(w);
ATw_v = sum(ATw .* v,'all');


% compare results (should be equal, so error should be 0)
err = abs(w_Av - ATw_v);

tol = 1e-10;

assert(err < tol,'Adjoint FAILED.')
if verbose
    fprintf('ADJOINT PASSED! err = %0.4e\n',err)
end


end



