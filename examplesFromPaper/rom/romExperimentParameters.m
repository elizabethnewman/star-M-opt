classdef romExperimentParameters < starMOptExperimentParameters

    properties
        k     = 2
        cN    = 50
        cMin  = 0.1
        cMax  = 5
        tN    = 31
        tMin  = 0
        tMax  = 5
        seed            = 42
        seedM           = 1234
        maxIter         = 100
        MInitialize     = 'I'
        date            
    end

    methods
        function[obj] = romExperimentParameters(varargin)
            obj = obj@starMOptExperimentParameters(varargin{:});

        end

        function setFilename(obj)
            obj.date     = date;
            obj.filename = sprintf([obj.date, '--k-%0.2d--init-%s'],obj.k, obj.MInitialize);
        end
    end

end
