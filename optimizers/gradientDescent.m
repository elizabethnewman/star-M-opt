classdef gradientDescent < optimizer

    properties
    end

    methods

        function[obj] = gradientDescent(varargin)
            obj = obj@optimizer(varargin{:})

            for k = 1:2:length(varargin)
                if isprop(obj,varargin{k})
                    obj.(varargin{k}) = varargin{k + 1}; 
                end
            end

            obj.linesearch.manifold = obj.manifold;

        end

        function[x,info] = solve(obj,fctn,x0)
            info = obj.initializeInfo(fctn);
            
            % initialize x and xOld
            [x,xOld] = deal(x0);
           
            % initial evaluation
            [f,infoFctn,df] = fctn.evaluate(x);
            
            nrmf = fronorm(df);
            nrm0 = nrmf;

            nrms  = fronorm(x * (0.5 * (x' * df - df' * x)));
            nrms0 = nrms;
            
            values = [0,0,f,nrm0,nrmf / nrm0, nrms0, nrms / nrms0,fronorm(x - xOld),obj.linesearch.alpha,0,0,infoFctn.values];
            if obj.verbose
                fprintf([repmat('%-12s',1,length(info.header)),'\n'], info.header{:});
                fprintf([info.frmt{:},'\n'],values); 
            end
            info.values = cat(1,info.values,values);
            if obj.storeIterates, info.x = cat(2,info.x,x(:)); end

            i = 1;
            while i <= obj.maxIter && (nrm0 > obj.absTol && nrmf / nrm0 > obj.relTol)
                startTime = tic;

                xOld = x;
                
                % take a step with line search
                switch obj.manifold
                    case 'Stiefel'
                        A = -x' * df;
                        s = 0.5 * (A - A');
                    otherwise
                        s = -df;
                end

                % update x with linesearch
                [x,infoLS] = obj.linesearch.step(fctn,x,s,f,df);
                endTime = toc(startTime);

                % re-evaluate
                f0 = f;
                [f,infoFctn,df]    = fctn.evaluate(x);
                nrmf        = fronorm(df);
                nrms  = fronorm(x * (0.5 * (x' * df - df' * x)));

                values      = [i,endTime,f,nrmf,nrmf / nrm0,nrms,nrms / nrms0,fronorm(x - xOld),infoLS.alpha,infoLS.lsIter,infoLS.flag,infoFctn.values];
                info.values = cat(1,info.values,values);
                if obj.storeIterates, info.x = cat(2,info.x,x(:)); end
               
                % if f0 < f
                %     disp('here')
                %     % set to previous xOld
                %     x = xOld;
                % end

                if obj.verbose && mod(i,obj.logInterval) == 0, fprintf([info.frmt{:},'\n'],values); end
                
            
                if infoLS.flag < 0
                    % disp('Line search break');
                    break;
                end
                
                i = i + 1;
            end


        end

        function test(~)
            rng(42);

            % create SPD matrix
            A = randn(100,2);
            A = A' * A + 1e-8 * eye(2);
            b = randn(2,1);

            xOpt = -A \ b;

            % create quadratic function            
            f    = @(x) 0.5 * x' * A * x + b' * x;
            df   = @(x) A * x + b;
            fctn = objFctnAnonymous(f,df);

            % solve
            x0 = randn(2,1);

            opt = gradientDescent('verbose',1,'maxIter',20);
            opt.linesearch = armijoLinesearch();

            [xSol,~] = opt.solve(fctn,x0);
              
            relErr = fronorm(xSol - xOpt) / fronorm(xOpt);
            fprintf('Relative Error = %0.4e\n',relErr);

            
            % create SPD matrix
            A = randn(2,100);
            [xOpt,~] = qr(randn(2));
            B = xOpt * A;
            
            % create quadratic function            
            f         = @(X) 0.5 * fronorm(X * A - B)^2;
            df        = @(X) (X * A - B) * A';
            fctnOrtho = objFctnAnonymous(f,df);

            % solve
            [x0,~] = qr(randn(2,2));

            opt = gradientDescent('verbose',1,'maxIter',20,'manifold','Stiefel');
            opt.linesearch = armijoLinesearch();

            [xSol,~] = opt.solve(fctnOrtho,x0);
              
            relErr = fronorm(xSol - xOpt) / fronorm(xOpt);
            fprintf('Relative Error = %0.4e\n',relErr);
            fprintf('Orthogonality = %0.4e\n',fronorm(xSol' * xSol - eye(size(xSol,2))));
        end

    end

end