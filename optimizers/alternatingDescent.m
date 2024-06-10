classdef alternatingDescent < optimizer

    properties
        linesearchX = constantLinesearch('alpha',1e-3)
        linesearchY = constantLinesearch('alpha',1e-3)
        manifoldX   = 'Euclidean'
        manifoldY   = 'Euclidean'
    end

    methods

        function[obj] = alternatingDescent(varargin)
            obj = obj@optimizer(varargin{:})

            for k = 1:2:length(varargin)
                if isprop(obj,varargin{k})
                    obj.(varargin{k}) = varargin{k + 1}; 
                end
            end

            obj.linesearchX.manifold = obj.manifoldX;
            obj.linesearchY.manifold = obj.manifoldY;

            obj.header = {'iter','time','Jc','|dJx|','|dJx|/|dJx0|','|x1-x0|','|dJy|','|dJy|/|dJy0|','|y1-y0|'};
            obj.frmt   = {'%-12.1f','%-12.2f','%-12.2e','%-12.2e','%-12.2e','%-12.2e','%-12.2e','%-12.2e','%-12.2e'};

            if ~isempty(obj.linesearch)
                obj.header = cat(2,obj.header,obj.linesearchX.header);
                obj.frmt   = cat(2,obj.frmt,obj.linesearchX.frmt);

                obj.header = cat(2,obj.header,obj.linesearchY.header);
                obj.frmt   = cat(2,obj.frmt,obj.linesearchY.frmt);
            end

        end

        function[x,y,info] = solve(obj,fctn,x0,y0)
            

            info = obj.initializeInfo();
            
            % initialize x and xOld
            [x,xOld] = deal(x0);
            [y,yOld] = deal(y0);
           
            % initial evaluation
            [f,~,dfX,dfY] = fctn.evaluate(x,y);

            [nrmfX,nrmfY] = deal(fronorm(dfX),fronorm(dfY));
            [nrm0X,nrm0Y] = deal(nrmfX,nrmfY);
            
            values = [0,0,f,nrm0X,nrmfX / nrm0X,fronorm(x - xOld),nrm0Y,nrmfY / nrm0Y,fronorm(y - yOld),obj.linesearchX.alpha,0,0,obj.linesearchY.alpha,0,0];
            if obj.verbose
                fprintf([repmat('%-12s',1,length(obj.header)),'\n'], obj.header{:});
                fprintf([obj.frmt{:},'\n'],values); 
            end
            info.values = cat(1,info.values,values);
            if obj.storeIterates
                info.x = cat(2,info.x,x(:)); 
                info.y = cat(2,info.y,y(:));
            end

            i = 1;
            while i <= obj.maxIter && ((nrm0X > obj.absTol && nrmfX / nrm0X > obj.relTol) || (nrm0Y > obj.absTol && nrmfY / nrm0Y > obj.relTol))
                % ------------------------------------------------------- %
                % update x
                startTime = tic;
                xOld = x;
                
                [f,~,dfX] = fctn.evaluate(x,y);

                % take a step with line search
                switch obj.manifoldX
                    case 'Stiefel'
                        A = -x' * dfX;
                        s = 0.5 * (A - A');
                    otherwise
                        s = -dfX;
                end

                % update x with linesearch
                fX          = objFctnAnonymous(@(x) fctn.evaluate(x,y), @(x) 0 * x);
                [x,infoLSX] = obj.linesearchX.step(fX,x,s,f,dfX);
                endTime = toc(startTime);

                % re-evaluate
                [f,~,dfX] = fctn.evaluate(x,y);
                nrmfX     = fronorm(dfX);
                values    = [i - 0.5,endTime,f,nrm0X,nrmfX / nrm0X,fronorm(x - xOld),nrm0Y,nrmfY / nrm0Y,fronorm(y - yOld),infoLSX.alpha,infoLSX.lsIter,infoLSX.flag,0,0,0];

                info.values = cat(1,info.values,values);
                
                if obj.storeIterates
                    info.x = cat(2,info.x,x(:)); 
                    info.y = cat(2,info.y,y(:));
                end

                if obj.verbose && mod(i,obj.logInterval) == 0, fprintf([obj.frmt{:},'\n'],values); end

                % ------------------------------------------------------- %
                % update y 
                startTime = tic;
                yOld = y;
                [f,~,~,dfY]    = fctn.evaluate(x,y);

                switch obj.manifoldY
                    case 'Stiefel'
                        A = -y' * dfY;
                        s = 0.5 * (A - A');
                    otherwise
                        s = -dfY;
                end

                % update y with linesearch
                fY          = objFctnAnonymous(@(y) fctn.evaluate(x,y), @(y) 0 * y);
                [y,infoLSY] = obj.linesearchY.step(fY,y,s,f,dfY);
                endTime = toc(startTime);

                % re-evaluate
                [f,~,~,dfY] = fctn.evaluate(x,y);
                nrmfY       = fronorm(dfY);
                values      = [i,endTime,f,nrm0X,nrmfX / nrm0X,fronorm(x - xOld),nrm0Y,nrmfY / nrm0Y,fronorm(y - yOld),0,0,0,infoLSY.alpha,infoLSY.lsIter,infoLSY.flag];

                info.values = cat(1,info.values,values);
                if obj.storeIterates
                    info.x = cat(2,info.x,x(:)); 
                    info.y = cat(2,info.y,y(:));
                end

               

                if obj.verbose && mod(i,obj.logInterval) == 0, fprintf([obj.frmt{:},'\n'],values); end
                
            
                if infoLSX.flag < 0 && infoLSY.flag < 0
                    disp('Line search break');
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