classdef lsExperimentParameters < starMOptExperimentParameters

    properties
        seed            = 42
        seedM           = 1234
        maxIter         = 50
        alpha           = 0.0
        lambda          = 0.0
        nPoints         = 100
        m               = [-1, 1]
        b               = [-1, 1]
        eta             = 0.0
        M               
    end

    methods
        function[obj] = lsExperimentParameters(varargin)
            obj = obj@starMOptExperimentParameters(varargin{:});

            obj.M = dctmtx(numel(obj.m));
        end

        function setFilename(obj)
            fname = sprintf([date,'--eta-%0.0e'],obj.eta);

            obj.filename = fname;

        end
    end

end

