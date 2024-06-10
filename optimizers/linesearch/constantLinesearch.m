classdef constantLinesearch < linesearch

    properties
        alpha = 1e-3
    end

    methods
        function[obj] = constantLinesearch(varargin)
            obj = obj@linesearch(varargin{:})
        end

        function[x,info] = step(obj,~,x0,s,varargin)

            switch obj.manifold
                case 'Stiefel'
                    x = x0 * expm(obj.alpha * s);
                otherwise
                    x = x0 + obj.alpha * s;
            end
            
            info = struct('alpha',obj.alpha,'lsIter',0,'flag',0);
        end
    end
end