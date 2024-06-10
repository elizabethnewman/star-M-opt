function[tests] = svdJacobianUnitTests
    tests = functiontests(localfunctions);
end

function testUnderdetermined(testCase)
    A = randn([3,10]);
    [U,S,V] = svd(A,"econ");
    runAllTests(A,U,S,V);
end

function testSquare(testCase)
    A = randn(9);
    [U,S,V] = svd(A,"econ");
    runAllTests(A,U,S,V);
end

function testOverdetermined(testCase)
    A = randn([20,5]);
    [U,S,V] = svd(A,"econ");
    runAllTests(A,U,S,V);
end

% IN PROGRESS
% function testDiagonal(testCase)
%     A = diag(diag(randn(9)));
%     [U,S,V] = svd(A,"econ");
%     runAllTests(A,U,S,V);
% end

% function testLowRankUnderdetermined(testCase)
%     A = randn(4,20);
%     A = [A;A];
%     [U,S,V] = svd(A,"econ");
%     r = 4;
% 
%     runAllTests(A,U(:,1:r),S(1:r,1:r),V(:,1:r));
% end

function runAllTests(A,U,S,V)

    [JacU,JacS,JacV] = svdJacobian(U,S,V);

    fctnU = @(x) svdU(x,size(U,2));
    CheckAdjoint(JacU);
    CheckJacobian(fctnU,JacU,[],1,A)
    CheckGradient(fctnU,JacU,[],0,A)

    fctnS = @(x) svdS(x,size(S,1));
    CheckAdjoint(JacS);
    CheckJacobian(fctnS,JacS,[],0,A)
    CheckGradient(fctnS,JacS,[],0,A)
    
    fctnV = @(x) svdV(x,size(V,2));
    CheckAdjoint(JacV);
    CheckJacobian(fctnV,JacV,[],0,A)
    CheckGradient(fctnV,JacV,[],0,A)
end

function[u] = svdU(A,k)
    [u,~,~] = svd(A,"econ");
    if exist('k','var') && ~isempty(k)
        u = u(:,1:k);
    end
end

function[s] = svdS(A,k)
    [~,s,~] = svd(A,"econ");
    if exist('k','var') && ~isempty(k)
        s = s(1:k,1:k);
    end
end

function[v] = svdV(A,k)
    [~,~,v] = svd(A,"econ");
    if exist('k','var') && ~isempty(k)
        v = v(:,1:k);
    end
end

