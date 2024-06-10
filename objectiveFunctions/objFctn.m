classdef objFctn < handle

    properties
        orthoFlag = 1
        avgFlag   = 1
        header
        frmt
    end

    methods
        function[obj] = objFctn(varargin)
            for k = 1:2:length(varargin)
                if isprop(obj,varargin{k})
                    obj.(varargin{k}) = varargin{k + 1};
                end
            end
        end

        function set.orthoFlag(obj, flag)
            obj.orthoFlag = flag;
        end

        function set.avgFlag(obj, flag)
            obj.avgFlag = flag;
        end

        function[info] = initializeInfo(obj)
            info        = struct('header',[],'frmt',[],'values',[]);
            info.header = obj.header;
            info.frmt   = obj.frmt;
            info.values = zeros(0,length(obj.header));
        end

        function[f,info,df] = evaluate(obj,x)
            f = [];
            info = struct();
            df = [];
            error('nyi');
        end

        function test(obj,x0,varargin)
            verbose    = 1;

            [f0,~,df0] = obj.evaluate(x0,varargin{:});

            nrm0 = fronorm(f0);

            p    = randn(size(x0));
            dfdp = sum(df0 .* p,'all');
            
            % main iteration
            N   = 15;
            err = zeros(N,3);
            
            base = 2;
            for k = 1:N
                h = base^(-k);

                ft = obj.evaluate(x0 + h * p,varargin{:});
                
                err0     = fronorm(f0 - ft) / nrm0;
                err1     = fronorm(f0 + h * dfdp - ft) / nrm0;
                err(k,:) = [h,err0,err1];
            
                [c0,d0] = convert2base(err0,base);
                [c1,d1] = convert2base(err1,base);
                    
                if verbose
                    fprintf('h=%0.2f x 2^(%0.2d)\tE0=%0.2f x 2^(%0.2d)\tE1=%0.2f x 2^(%0.2d)\n',1,-k,c0,d0,c1,d1);
                end
            end
        end
    end


end