classdef linesearch < matlab.mixin.Copyable

    properties
        manifold    = 'Euclidean'
        absTol      = 1e-10
        relTol      = 1e-10           
        header      = {'alphaLS','iterLS','flagLS'}
        frmt        = {'%-12.2e','%-12d','%-12d'}
        verbose     = 0
        logInterval = 1
    end

    methods
        function[obj] = linesearch(varargin)
            for k = 1:2:length(varargin)
                if isprop(obj,varargin{k})
                    obj.(varargin{k}) = varargin{k + 1}; 
                end
            end

            % reset manifold
            if isempty(obj.manifold), obj.manifold = 'Euclidean'; end
        end

        function[info] = initializeInfo(obj)
            info        = struct('header',[],'frmt',[],'values',[]);
            info.header = obj.header;
            info.frmt   = obj.frmt;
            info.values = zeros(0,length(obj.header));
        end

        function[x,info] = step(obj,fctn,x0,s,varargin)
            x    = [];
            s    = [];
            info = struct();
            error('Not yet implemented')
        end

    end

end
