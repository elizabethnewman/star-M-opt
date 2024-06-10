classdef indexTrackingExperimentParameters < starMOptExperimentParameters

    properties
        initDateTrain   = '13-Mar-2021'
        endDateTrain    = '13-Mar-2023'
        initDateTest    = '14-Mar-2023'
        endDateTest     = '14-Apr-2023'
        seed            = 42
        seedM           = 1234
        maxIter         = 20
        alpha           = 5e-3 % good results for normal defaults
        % alpha           = 3e-2
    end

    methods
        function[obj] = indexTrackingExperimentParameters(varargin)
            obj = obj@starMOptExperimentParameters(varargin{:});

        end

        function setFilename(obj)
            fname = sprintf([date, '--initDateTrain-%s--endDateTrain-%s--initDateTest-%s--endDateTest-%s'], ...
                obj.initDateTrain, obj.endDateTrain, obj.initDateTest, obj.endDateTest);

            obj.filename = fname;

        end
    end

end

