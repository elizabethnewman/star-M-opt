classdef objFctnAnonymous < objFctn
    properties
        f
        df
    end

    methods
        function[obj] = objFctnAnonymous(f,df,idxIn,idxOut,varargin)
            obj     = obj@objFctn(varargin{:});

            % df is optional and is a function handle
            % varargin are fixed inputs
            % idxIn is the index of the variable we will plug into our
            % function
            % idxOut is the idx of the output argument for gradient
            % information
            
            if isa(f,'objFctn')
                obj.f  = @(x) obj.getSpecificOutput(f,idxIn,1,x,varargin{:});
                obj.df = @(x) obj.getSpecificOutput(f,idxIn,idxOut,x,varargin{:});

            elseif isa(f, 'function_handle') && exist('df', 'var')
                obj.f   = f;
                obj.df  = df;
            else
                error('inputs must be objFctn or function_handles for f and df')
            end

        end

        function[f,info,dfdM] = evaluate(obj,M)
            info = obj.initializeInfo();
            dfdM = [];
            doGrad = (nargout > 2);
            
            f = obj.f(M);
            if doGrad
                dfdM = obj.df(M);
            end
        end

    end

    methods (Static)

        function[out] = getSpecificOutput(f,idxIn,idxOut,x,varargin)

            inputs = [varargin(1:idxIn-1),{x},varargin(idxIn:end)];
            switch idxOut
                case 2
                    [~,out] = f.evaluate(inputs{:});
                case 3
                    [~,~,out] = f.evaluate(inputs{:});
                case 4
                    [~,~,~,out] = f.evaluate(inputs{:});
                otherwise
                    out = f.evaluate(inputs{:});
            end
        end


    end

end