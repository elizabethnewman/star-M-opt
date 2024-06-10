classdef digitsExperimentParameters < starMOptExperimentParameters

    properties
        k               = 1
        seed            = 42
        seedM           = 1234
        seedTransfer    = 30
        maxIter         = 100
        innerMaxIter    = 10
        nImg            = 30
        nClasses        = 10
        nTestRuns       = 50
        noiseLevel      = 0.0
        initializeM     = 'identity'
        permuteFlag     = false
        transferDistribution = 'random'
    end

    methods
        function[obj] = digitsExperimentParameters(varargin)
            obj = obj@starMOptExperimentParameters(varargin{:});
        end

        function setFilename(obj)
            obj.filename = sprintf([date, '--k-%0.2d--init-%s'],obj.k,obj.initializeM);
        end

    end

end
