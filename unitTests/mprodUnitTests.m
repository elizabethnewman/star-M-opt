function[tests] = mprodUnitTests
    tests = functiontests(localfunctions);

function testMprodFacewise(testCase)
    A = randn(3,4,5);
    B = randn(4,7,5);
    
    M = eye(5);
    testCase.verifyEqual(facewise(A,B),mprod(A,B,M))

function testSize(testCase)
    A = randn(5,4,3);
    B = randn(4,7,3);
    M = randn(3,3);
    C = mprod(A,B,M);

    sC = [5,7,3];
    testCase.verifySize(C,sC); 

function testJacA(testCase)
    A = randn(3,4,5);
    B = randn(4,7,5);
    M = randn(5);

    [C,JacA] = mprod(A,B,M);
    fctn = @(a) mprod(a,B,M);
    CheckAdjoint(JacA);
    CheckJacobian(fctn,JacA,[],0)
    CheckGradient(fctn,JacA,[],0)


function testJacB(testCase)
    A = randn(3,4,5);
    B = randn(4,7,5);
    M = randn(5);

    [C,JacA,JacB] = mprod(A,B,M);
    fctn = @(b) mprod(A,b,M);
    CheckAdjoint(JacB);
    CheckJacobian(fctn,JacB,[],0)
    CheckGradient(fctn,JacB,[],0)


function testJacM(testCase)
    A = randn(3,4,5);
    B = randn(4,7,5);
    M = randn(5);

    [C,JacA,JacB,JacM] = mprod(A,B,M);
    fctn = @(m) mprod(A,B,m);
    CheckAdjoint(JacM);
    
    % need to give M as direction (special case because using implicit differentiation)
    CheckJacobian(fctn,JacM,[],0,M)
    CheckGradient(fctn,JacM,[],0,M)

