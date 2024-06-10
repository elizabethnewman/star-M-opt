classdef starMOptExperimentParameters < handle

    properties
        filename
    end


    methods
        function[obj] = starMOptExperimentParameters(varargin)

            for k = 1:2:length(varargin)
                obj.(varargin{k}) = varargin{k + 1};
            end
            % % overwrite defaults
            % fieldsInput = fields(options);
            % for k = 1:length(fieldsInput)
            %     obj.options.(fieldsInput{k}) = options.(fieldsInput{k});
            % end
            obj.setFilename()
            
        end

        function setFilename(obj)
            obj.filename = 'tmp';
        end

        function[options] = obj2struct(obj)
            tmp = fields(obj);
            options = struct();
            for i = 1:length(tmp)
                options.(tmp{i}) = obj.(tmp{i});
            end
        end
    end

end

