classdef optimizer < handle

    properties
        maxIter         = 10
        absTol          = 1e-10
        relTol          = 1e-10
        linesearch      = constantLinesearch('alpha',1e-3)
        manifold        = 'Euclidean'
        header          = {'iter','time','Jc','|dJc|','|dJc|/|dJ0|','|gradJ|','|gradJ|/|gradJ0|','|x1-x0|'}
        frmt            = {'%-12d','%-12.2f','%-12.2e','%-12.2e','%-12.2e','%-12.2e','%-12.2e','%-12.2e'}
        verbose         = 0
        logInterval     = 1
        storeIterates   = 0
    end

    methods
        function[obj] = optimizer(varargin)
            for k = 1:2:length(varargin)
                if isprop(obj,varargin{k})
                    obj.(varargin{k}) = varargin{k + 1}; 
                end
            end

            % reset manifold
            if isempty(obj.manifold), obj.manifold = 'Euclidean'; end
            
            obj.linesearch.manifold = obj.manifold;
            if ~isempty(obj.linesearch)
                obj.header = cat(2,obj.header,obj.linesearch.header);
                obj.frmt   = cat(2,obj.frmt,obj.linesearch.frmt);
            end
        end

        function[info] = initializeInfo(obj,fctn)
            info        = struct('header',[],'frmt',[],'values',[]);
            info.header = obj.header;
            info.frmt   = obj.frmt;
            info.x      = [];
            info.y      = [];

            if exist('fctn','var') && ~isempty(fctn)
                info.header = cat(2,info.header,fctn.header);
                info.frmt   = cat(2,info.frmt,fctn.frmt);
            end

            info.values = zeros(0,length(info.header));
        end

        function[x,info] = solve(obj,fctn,x0)
            x    = [];
            info = struct();
            error('Not yet implemented')
        end
        

    end

end
