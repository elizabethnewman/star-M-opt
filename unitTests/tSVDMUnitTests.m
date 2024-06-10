function[tests] = tSVDMUnitTests
    tests = functiontests(localfunctions);
end

function testTSVDMvsFacewise(testCase)
    A = randn(3,4,5);
    M = 1;
    [U,S,V] = tSVDM(A,M);
    [UF,SF,VF] = facewiseSVD(A);
    
    testCase.verifyEqual(U,UF);
    testCase.verifyEqual(S,SF)
    testCase.verifyEqual(V,VF)
end

function testTSVDMSize(testCase)
    A = randn(3,4,5);
    M = randn(5);
    k = 2;
    [U,S,V] = tSVDM(A,M,k);

    
    testCase.verifySize(U,[size(A,1),k,size(A,3)]);
    testCase.verifySize(S,[k,k,size(A,3)])
    testCase.verifySize(V,[size(V,1),k,size(A,3)])
end

function testJacU(testCase)
    A = randn(3,4,5);
    M = randn(5);
    [~,~,~,JacU,~,~] = tSVDM(A,M);
    
    
    fctn = @(x) tSVDMU(x,M);
    CheckAdjoint(JacU);
    CheckJacobian(fctn,JacU,[],0,A)
    CheckGradient(fctn,JacU,[],0,A)

end


function testJacS(testCase)
    A = randn(20,10,9);
    M = randn(9);
    [~,~,~,~,JacS,~] = tSVDM(A,M);
    
    
    fctn = @(x) tsvdMS(x,M);
    CheckAdjoint(JacS);
    CheckJacobian(fctn,JacS,[],0,A)
    CheckGradient(fctn,JacS,[],0,A)

end

function testJacV(testCase)
    A = randn(10,2,7);
    M = randn(7);
    [~,~,~,~,~,JacV] = tSVDM(A,M);
    
    
    fctn = @(x) tsvdMV(x,M);
    CheckAdjoint(JacV);
    CheckJacobian(fctn,JacV,[],0,A)
    CheckGradient(fctn,JacV,[],0,A)

end

function testJacMU(testCase)
    A = randn(3,4,5);
    M = randn(5);
    [~,~,~,~,~,~,JacMU,~,~] = tSVDM(A,M);
    
    
    fctn = @(m) tSVDMU(A,m);
    CheckAdjoint(JacMU);
    CheckJacobian(fctn,JacMU,[],0,M)
    CheckGradient(fctn,JacMU,[],0,M)

end


function testJacMS(testCase)
    A = randn(20,10,9);
    M = randn(9);
    [~,~,~,~,~,~,~,JacMS,~] = tSVDM(A,M);
    
    
    fctn = @(m) tsvdMS(A,m);
    CheckAdjoint(JacMS);
    CheckJacobian(fctn,JacMS,[],0,M)
    CheckGradient(fctn,JacMS,[],0,M)

end

function testJacMV(testCase)
    A = randn(10,2,7);
    M = randn(7);
    [~,~,~,~,~,~,~,~,JacMV] = tSVDM(A,M);
    
    
    fctn = @(m) tsvdMV(A,m);
    CheckAdjoint(JacMV);
    CheckJacobian(fctn,JacMV,[],0,M)
    CheckGradient(fctn,JacMV,[],0,M)

end


function[u] = tSVDMU(A,M)
    [u,~,~] = tSVDM(A,M);
end

function[s] = tsvdMS(A,M)
    [~,s,~] = tSVDM(A,M);
end

function[v] = tsvdMV(A,M)
    [~,~,v] = tSVDM(A,M);
end

