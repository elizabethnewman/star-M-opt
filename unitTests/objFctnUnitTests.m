function[tests] = objFctnUnitTests
    tests = functiontests(localfunctions);
end

function testobjFctnAnonymous(testCase)
    f  = @(x) sum(x.^2, 'all');
    df = @(x) 2 * x;

    fctn = objFctnAnonymous(f, df);
    CheckObjectiveFunction(fctn, randn(3,4))
end

function testobjFctnTikhonovLeastSquaresM(testCase)
    for orthoFlag = [0, 1]
        for avgFlag = [0, 1]
            A = randn(randi([1,10]), randi([1,10]), randi([2,10]));
            B = randn(size(A,1), randi([1,10]), size(A,3));
            X = randn(size(A,2),size(B,2),size(A,3));
            
            M = randn(size(A,3));
            if orthoFlag
                [M,~] = qr(M);
            end

            lambda = rand(1);
            alpha  = rand(1);
            options = {'orthoFlag', orthoFlag, 'avgFlag', avgFlag};
            fctn = objFctnTikhonovLeastSquares(A,B,lambda,alpha,options{:});

            idxIn  = 1;
            idxOut = 3;
            fctnAnonM = objFctnAnonymous(fctn,[],idxIn,idxOut,X);
        
            CheckObjectiveFunction(fctnAnonM, M,'verbose',true)
        end 
    end
end

function testobjFctnTikhonovLeastSquaresX(testCase)

    for orthoFlag = [0, 1]
        for avgFlag = [0, 1]
            A = randn(randi([1,10]), randi([1,10]), randi([2,10]));
            B = randn(size(A,1), randi([1,10]), size(A,3));
            X = randn(size(A,2),size(B,2),size(A,3));
            
            M = randn(size(A,3));
            if orthoFlag
                [M,~] = qr(M);
            end

            lambda = rand(1);
            alpha  = rand(1);
            options = {'orthoFlag', orthoFlag, 'avgFlag', avgFlag};
            fctn = objFctnTikhonovLeastSquares(A,B,lambda,alpha,options{:});
                    
            idxIn  = 2;
            idxOut = 4;
            fctnAnonX = objFctnAnonymous(fctn,[],idxIn,idxOut,M);

            CheckObjectiveFunction(fctnAnonX, X)
        end 
    end

end


function testobjFctnTikhonovLeastSquaresVarPro(testCase)
    count = 1;
    for orthoFlag = [0,1]
        for avgFlag = [0,1]
            A = randn(randi([11,20]), randi([1,10]), randi([2,10]));
            B = randn(size(A,1), randi([1,10]), size(A,3));
            
            
            [M,~] = qr(randn(size(A,3)));
            

            lambda = rand(1);
            alpha  = rand(1);
            options = {'orthoFlag', orthoFlag, 'avgFlag', avgFlag};
            fctn = objFctnTikhonovLeastSquaresVarPro(A,B,lambda,alpha,options{:});
            disp(count)
            CheckObjectiveFunction(fctn,M,'verbose',true)
            count = count + 1;
        end 
    end
end


function testobjFctnLowRank(testCase)

    for k = 1:5
        A = randn(randi([10,20]), randi([10,20]), randi([2,10]));
        
        [M,~] = qr(randn(size(A,3)));
        
        options = {'orthoFlag', 1};
        fctn = objFctnLowRank(A,k,options{:});
        
        CheckObjectiveFunction(fctn,M)
    end

end



function testConstrainedLeastSquaresVarPro(testCase)

    A = rand(randi([1,5]),randi([6,10]),randi([2,10]));
    B = rand(size(A,1), randi([1,10]), size(A,3));
    [M,~] = qr(randn(size(A,3)));
    

    lambda = rand(1);
    alpha  = rand(1);
    options = {'orthoFlag', 1};
    fctn = objFctnConstrainedLeastSquaresVarPro(A,B,lambda,alpha,[],options{:});
    
    CheckObjectiveFunction(fctn,M,'verbose',0)

end

