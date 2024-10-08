classdef lipschitzExperimentParameters < starMOptExperimentParameters

    properties
        seed    = 42

        probType = 'random';     % {'ill-conditioned', 'random', 'index-tracking'}
        fctnType = 'ls-varpro';  % {'ls-varpro', 'ls-no-varpro', 'ls-constrained-varpro', 'low-rank'}
        
        % number of runs
        nRuns    = 100;    
        
        MType    = 'random';    % {'random', 'close'}: 'close' means M2 = M1 * expm(eps * Omega)
        perturbI = false;       % M1 = expm(eps * Omega);
        eps      = 1e-8;
        
        % info about A
        n1       = 100;
        n2       = 7;
        n3       = 20;
        nrmA     = 1e1;
        condA    = 1e10;
        
        % info about B (if needed)
        m        = 10; % size(B,2)
        nrmB     = 1e0;
        
        % info about X (if needed)
        XType    = 'random';          % {'random', 'optimal'}
        
        % additional info for svd (if needed)
        k        = 4;
    end

    methods
        function[obj] = lipschitzExperimentParameters(varargin)
            obj = obj@starMOptExperimentParameters(varargin{:});
        end

        function setFilename(obj)
            obj.filename = sprintf([date, '--%s'],obj.probType);
        end

    end

end
