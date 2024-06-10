function tests = utilsUnitTests
    tests = functiontests(localfunctions);

function testModeFoldUnFold(testCase)
    A  = randn(3,4,5,3);
    sA = size(A);

    testCase.assertEqual(modeFold(modeUnfold(A),sA),A)

    for k = 1:ndims(A)
        testCase.assertEqual(modeFold(modeUnfold(A,k),sA,k),A); 
    end

function testTran(testCase)
    A  = randn(3,4,5,3);

    testCase.assertEqual(tran(tran(A)),A); 

    testCase.assertEqual(tran(tran(A,1),1),A);


function testFacewiseEye(testCase)
    A = facewiseEye([3,5,6,2]);
    I = eye(size(A,1:2));

    tol = 1e-16;
    assert(fronorm(A - I) < tol); 
