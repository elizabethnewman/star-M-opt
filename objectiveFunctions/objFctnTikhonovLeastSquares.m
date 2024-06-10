classdef objFctnTikhonovLeastSquares < objFctn

    properties
        A
        B
        lambda = 0.0
        alpha  = 0.0
    end

    methods

        function[obj] = objFctnTikhonovLeastSquares(A,B,lambda,alpha,varargin)
            if size(A,1) ~= size(B,1) || size(A,3) ~= size(B,3)
                error('incompatible tensor sizes')
            end
    
            obj = obj@objFctn(varargin{:});
            obj.A = A;
            obj.B = B;

            if exist('lambda','var'), obj.lambda = lambda; end
            if exist('alpha','var'), obj.alpha = alpha; end

        end

        function[f,info,dfdM,dfdX] = evaluate(obj,M,X)
            info    = obj.initializeInfo();
            dfdM    = [];
            doGrad  = (nargout > 2);
            doGradX = (nargout > 3);

            if ~exist('X', 'var') && iscell(M)
                X = M{2};
                M = M{1};
            end
            
            if doGradX
                [AX,~,JacX,JacM] = mprod(obj.A,X,M,'orthoFlag',obj.orthoFlag);
            elseif doGrad
                [AX,~,~,JacM] = mprod(obj.A,X,M,'orthoFlag',obj.orthoFlag);
            else
                AX = mprod(obj.A,X,M,'orthoFlag',obj.orthoFlag);
            end
            
            R = AX - obj.B;
            
            n = 1.0;
            if obj.avgFlag, n = numel(R); end

            f = (0.5 / n) * fronorm(R).^2 + (0.5 * obj.alpha) * fronorm(X).^2 + (0.5 * obj.lambda) * fronorm(M).^2 ;

            if doGradX
                dfdX = JacX.AT((1 / n) * R) + obj.alpha * X;
            end
            
            if doGrad
                dfdM = JacM.AT((1 / n) * R) + obj.lambda * M;
            end

        end

    end
    


end


