classdef armijoLinesearch < linesearch
    % Algorithm 4.2 in https://www.nicolasboumal.net/book/IntroOptimManifolds_Boumal_2023.pdf
    properties
        alpha       = 1e0
        beta        = 0.5
        gamma       = 1e-4
        c           = 1e-4
        maxIter     = 20
        alphaTol    = 1e-14
    end

    methods
        function[obj] = armijoLinesearch(varargin)
            obj = obj@linesearch(varargin{:});

            for k = 1:2:length(varargin)
                if isprop(obj,varargin{k})
                    obj.(varargin{k}) = varargin{k + 1}; 
                end
            end
        end

        function[x,info] = step(obj,fctn,x0,s,f0,g0)
            
            % check for descent direction (Euclidean case)
            tau = -s(:)' * s(:);
            % tau = g0(:)' * s(:);
            % 
            if tau >= 0
                disp('not a descent direction - switching to s = -g')
                s = -g0;
            end
            
            % initialize step size
            mu = obj.alpha;

            % main iteration
            iter = 1;
            while iter <= obj.maxIter

                % take step on manifold
                switch obj.manifold
                    case 'Stiefel'
                        x = x0 * expm(mu * s);
                    otherwise
                        x = x0 + mu * s;
                end
                
                ft = fctn.evaluate(x);

                if ft < f0 + mu * obj.gamma * tau
                    flag = 1;
                    break;
                end

                % update step size
                mu    = obj.beta * mu;
                iter  = iter + 1;
            end
            
            if iter > obj.maxIter
                flag = -1;
            end

            if mu < obj.alphaTol
                flag = -2;
            end

            info = struct('alpha',mu,'lsIter',iter,'flag',flag);
            
            % update step size for next time
            if iter == 1
                obj.alpha = 2 * mu;
            else
                obj.alpha = mu;
            end
        end
    end
end