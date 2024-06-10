function[tests] = modeProductUnitTests
    tests = functiontests(localfunctions);

function testModeProductSize(testCase)
    A = randn(5,4,3,2,8);
    M1 = randn(7,5);
    B = modeProduct(A,M1,'dim',1);

    sB = [7,4,3,2,8];
    testCase.verifySize(B,sB); 

    M4 = randn(8);
    B = modeProduct(A,M4);
    sB = [5,4,3,2,8];
    testCase.verifySize(B,sB);

function testModeProductScalar(testCase)
    A  = randn(5,4,3,2,8);
    M3 = 1;
    B  = modeProduct(A,M3,'dim',3);
    sB = [5,4,3,2,8];
    testCase.verifySize(B,sB);
    testCase.verifyEqual(B,A);

function testModeProductNoFlags(testCase)

    A  = randn(5,4,3,2);
    M3 = randn(7,3);
    M4 = randn(2,2);
    
    % ------------------------------------------------------------- %
    % standard product
    [B,JacA,JacM] = modeProduct(A,M3,'dim',3);
    
    fctn = @(a) modeProduct(a,M3,'dim',3);
    CheckAdjoint(JacA);
    CheckJacobian(fctn,JacA,[],0)
    CheckGradient(fctn,JacA,[],0)
    
    fctn = @(m) modeProduct(A,m,'dim',3);
    CheckAdjoint(JacM);
    CheckJacobian(fctn,JacM,[],0)
    CheckGradient(fctn,JacM,[],0)
    
    [B,JacA,JacM] = modeProduct(A,M4);
    
    fctn = @(a) modeProduct(a,M4);
    CheckAdjoint(JacA);
    CheckJacobian(fctn,JacA,[],0)
    CheckGradient(fctn,JacA,[],0)
    
    fctn = @(m) modeProduct(A,m);
    CheckAdjoint(JacM);
    CheckJacobian(fctn,JacM,[],0)
    CheckGradient(fctn,JacM,[],0)

function testModeProductTransFlag(testCase)

    A  = randn(5,4,3,2);
    M1 = randn(5,5);
    M2 = randn(4,1);
    M3 = randn(7,3);
    M4 = randn(2,2);
    
    [B,JacA,JacM] = modeProduct(A,M2,'dim',2,'transFlag',1);

    fctn = @(a) modeProduct(a,M2,'dim',2,'transFlag',1);
    CheckAdjoint(JacA);
    CheckJacobian(fctn,JacA,[],0)
    CheckGradient(fctn,JacA,[],0)
    
    fctn = @(m) modeProduct(A,m,'dim',2,'transFlag',1);
    CheckAdjoint(JacM);
    CheckJacobian(fctn,JacM,[],0)
    CheckGradient(fctn,JacM,[],0)

function testModeProductOrthoFlag(testCase)

    A  = randn(5,4,3,2);
    M1 = randn(5,5);
    M2 = randn(4,1);
    M3 = randn(7,3);
    M4 = randn(2,2);
    
    [B,JacA,JacM] = modeProduct(A,M4,'orthoFlag',1);
    
    fctn = @(a) modeProduct(a,M4,'orthoFlag',1);
    CheckAdjoint(JacA);
    CheckJacobian(fctn,JacA,[],0)
    CheckGradient(fctn,JacA,[],0)
    
    fctn = @(m) modeProduct(A,m,'orthoFlag',1);
    CheckAdjoint(JacM);
    CheckJacobian(fctn,JacM,[],0)
    CheckGradient(fctn,JacM,[],0)

function testModeProductInvFlag(testCase)

    A  = randn(5,4,3,2);
    M1 = randn(5,5);
    M2 = randn(4,1);
    M3 = randn(7,3);
    M4 = randn(2,2);
    
    [B,JacA,JacM] = modeProduct(A,M1,'dim',1,'invFlag',1);
    
    
    fctn = @(a) modeProduct(a,M1,'dim',1,'invFlag',1);
    CheckAdjoint(JacA);
    CheckJacobian(fctn,JacA,[],0)
    CheckGradient(fctn,JacA,[],0)
    
    fctn = @(m) modeProduct(A,m,'dim',1,'invFlag',1);
    CheckAdjoint(JacM);
    
    % need to give M1 as direction (special case because using implicit differentiation)
    CheckJacobian(fctn,JacM,[],0,M1)
    CheckGradient(fctn,JacM,[],0,M1)
