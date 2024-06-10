classdef objFctnConstrainedLeastSquaresVarPro < objFctn

    properties
        A
        B
        lambda  = 0.0
        alpha   = 0.0
        options = optimoptions('fmincon','Display','none','Algorithm','interior-point', 'SpecifyObjectiveGradient', true, ...
            'MaxFunctionEvaluations', 5000, 'ConstraintTolerance', 1e-7, 'OptimalityTolerance', 1e-7);
    end

    methods

        function[obj] = objFctnConstrainedLeastSquaresVarPro(A,B,lambda,alpha,options,varargin)
            if size(A,1) ~= size(B,1) || size(A,3) ~= size(B,3)
                error('incompatible tensor sizes')
            end
    
            obj     = obj@objFctn(varargin{:});
            obj.A   = A;
            obj.B   = B;

            if exist('lambda','var') && ~isempty(lambda),   obj.lambda  = lambda;   end
            if exist('alpha','var') && ~isempty(alpha),    obj.alpha   = alpha;    end
            if exist('options','var') && ~isempty(options),  obj.options = options;  end

        end

        function[f,info,dfdM,dfdX] = evaluate(obj,M,X0,seed)
            info    = obj.initializeInfo();
            dfdM    = [];
            doGrad  = (nargout > 2);
            doGradX = (nargout > 3);

            sX = [size(obj.A,2),size(obj.B,2),size(obj.A,3)];
            
            if ~exist('X0', 'var') && iscell(M)
                X0 = M{2};
                M  = M{1};
            elseif ~exist('X0','var') || isempty(X0)
                if exist('seed', 'var') && ~isempty(seed), rng(seed); end
                X0 = rand(sX);
                X0 = X0 / sum(X0(:));
            end

            Aeq = ones(1,numel(X0));
            beq = 1;
            lb  = zeros(numel(X0),1);

            % solve for X (varpro step)
            fctn = @(x) obj.myFctn(x,M);
            x    = fmincon(fctn,X0(:),[],[],Aeq,beq,lb,[],[],obj.options);
            X    = reshape(x,sX);

            % store solution for X
            info.X = X;
            
            if doGradX
                [f,dfdX,dfdM] = obj.myFctn(X,M);
            elseif doGrad
                [f,~,dfdM] = obj.myFctn(X,M);
            else
                f = obj.myFctn(X,M);
            end

        end


        function[f,dfdX,dfdM] = myFctn(obj,X,M)
            X = reshape(X,size(obj.A,2),[],size(M,2));

            n = 1.0;
            if obj.avgFlag 
                % n = size(obj.A,1) * size(obj.A,3); 
                n = size(obj.A,1);
            end
            
            % TODO: only compute gradients if needed
            [AX,~,JacX,JacM] = mprod(obj.A,X,M,'orthoFlag',obj.orthoFlag);

            res = AX - obj.B;

            [r,drdX,drdM] = obj.myReg(X,M);
            
            f = (0.5 / n) * fronorm(res).^2 + r;

            if nargout > 1
                dfdX = JacX.AT((1 / n) * res) + drdX;
            end

            if nargout > 2
                dfdM = JacM.AT((1 / n) * res) + drdM;
            end

           
        end

        function[r,drdX,drdM] = myReg(obj,X,M)
            r = 0.5 * obj.alpha * fronorm(X).^2 + 0.5 * obj.lambda * fronorm(M).^2;

            if nargout > 1
                drdX = obj.alpha * X;
            end

            if nargout > 2
                drdM = obj.lambda * M;
            end


        end
    end
    


end


