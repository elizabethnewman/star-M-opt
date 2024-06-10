function[tests] = facewiseUnitTests
    tests = functiontests(localfunctions);

function testFacewiseSize(testCase)
    A = randn(5,4,3,2);
    B = randn(4,7,3,2);
    C = facewise(A,B);

    sC = [5,7,3,2];
    testCase.verifySize(C,sC); 

function testFacewiseJacA(testCase)
    A = randn(5,4,3,2);
    B = randn(4,7,3,2);
    [C,JacA] = facewise(A,B);

    fctn = @(a) facewise(a,B);
    CheckAdjoint(JacA);
    CheckJacobian(fctn,JacA,[],0)
    CheckGradient(fctn,JacA,[],0)

function testFacewiseJacB(testCase)
    A = randn(5,4,3,2);
    B = randn(4,7,3,2);
    [C,JacA,JacB] = facewise(A,B);

    fctn = @(a) facewise(a,B);
    CheckAdjoint(JacA);
    CheckJacobian(fctn,JacA,[],0)
    CheckGradient(fctn,JacA,[],0)

    fctn = @(b) facewise(A,b);
    CheckAdjoint(JacB);
    CheckJacobian(fctn,JacB,[],0)
    CheckGradient(fctn,JacB,[],0)

