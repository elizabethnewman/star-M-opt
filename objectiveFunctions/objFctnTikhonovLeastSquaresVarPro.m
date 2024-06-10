classdef objFctnTikhonovLeastSquaresVarPro < objFctn

    properties
        A
        B
        X
        lambda      = 0.0  % regularization on M
        alpha       = 0.0  % regularization on X
    end

    methods

        function[obj] = objFctnTikhonovLeastSquaresVarPro(A,B,lambda,alpha,varargin)
            if size(A,1) ~= size(B,1) || size(A,3) ~= size(B,3)
                error('incompatible tensor sizes')
            end
    
            obj = obj@objFctn(varargin{:});
            obj.A = A;
            obj.B = B;

            if exist('lambda','var'), obj.lambda = lambda; end
            if exist('alpha','var'), obj.alpha = alpha; end

            if ~obj.orthoFlag
                warning([class(obj), ' only supports orthogonal transforms - turning orthogonal flag on'])
                obj.orthoFlag = 1;
            end
        end

        function[f,info,dfdM,dfdX] = evaluate(obj,M)
            info    = obj.initializeInfo();
            dfdM    = [];
            doGradM = (nargout > 2);
            doGradX = (nargout > 3);

            options = {'orthoFlag',obj.orthoFlag};

            % solve for X
            obj.X = obj.solve(M);
            
            % compute residuals (gradient through X is zero!)
            if doGradX
                [AX,~,JacX,JacM] = mprod(obj.A,obj.X,M,options{:});
            elseif doGradM
                [AX,~,~,JacM] = mprod(obj.A,obj.X,M,options{:});
            else
                AX = mprod(obj.A,obj.X,M,options{:});
            end
            
            R = AX - obj.B;

            n = 1.0;
            if obj.avgFlag, n = numel(R); end
            
            f = (0.5 / n) * fronorm(R).^2 + (0.5 * obj.alpha) * fronorm(obj.X).^2 + (0.5 * obj.lambda) * fronorm(M).^2;

            if doGradX
                % should be 0
                dfdX = JacX.AT((1 / n) * R) + obj.alpha * obj.X;
            end
            
            if doGradM
                % skip the varpro part
                dfdM = JacM.AT((1 / n) * R) + obj.lambda * M;

                % % testing the varpro part for a sanity check
                % [H,~,~,JacHM]           = mprod(tran(obj.A),obj.A,M,options{:});
                % [rhs,JacAT,JacB,JacABM] = mprod(tran(obj.A),obj.B,M,options{:});
                % [lhs,JacH,JacX,JacHXM]  = mprod(H,obj.X,M,options{:});
                % 
                % 
                % tmp = JacABM.AT(rhs) - JacHM.AT(JacH.AT(obj.X));
                % 
                % % apply inverse
                % HHat    = modeProduct(H,M);
                % rhsHat  = modeProduct(rhs,M);
                % tHat    = pagemldivide(HHat,rhsHat);
                % t       = modeProduct(tHat,M');
            end

        end
        
        function[X] = solve(obj,M)

            % get sizes and normalizing constants
            [n1,n2,n3] = size(obj.A);

            n = 1.0;
            if obj.avgFlag, n = n1 * size(obj.B,2) * n3; end
            beta = 1 / sqrt(n);

            % move to transform domain
            AHat = modeProduct(obj.A,M);
            BHat = modeProduct(obj.B,M);
            
            % add regularization
            AHatL = cat(1,beta * AHat,obj.alpha * eye(n2) .* ones(1,1,n3));
            BHAtZ = cat(1,beta * BHat,zeros(n2,size(obj.B,2),n3));

            % solve
            XHat = pagemldivide(AHatL,BHAtZ);

            % return to spatial domain
            X    = modeProduct(XHat,M,'invFlag',1,'orthoFlag',obj.orthoFlag);
        end
    end
    


end


