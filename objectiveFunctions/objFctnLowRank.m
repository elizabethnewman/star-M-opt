classdef objFctnLowRank < objFctn

    properties
        A
        k
        nrmA
    end

    methods

        function[obj] = objFctnLowRank(A,k,varargin)
            obj = obj@objFctn(varargin{:});
            obj.A = A;
            obj.nrmA = fronorm(A);
            obj.k = k;
            
            if ~obj.orthoFlag
                warning([class(obj), ' only supports orthogonal transforms - turning orthogonal flag on'])
                obj.orthoFlag = 1;
            end
            obj.header = {'|A - Ak|','|A - Ak|/|A|'};
            obj.frmt   = {'%-15.2e','%-15.2e'};
        end

        function[f,info,dfdM] = evaluate(obj,M)
            info    = obj.initializeInfo();
            dfdM    = [];
            doGrad  = (nargout > 2);
            
            r           = min(size(obj.A,1),size(obj.A,2));
            zeroColPad  = @(y) cat(2,y,zeros(size(y,1),r-size(y,2),size(obj.A,3)));
            zeroColCut  = @(x) x(:,1:obj.k,:);
            zeroPad     = @(y) cat(1,zeroColPad(y),zeros(r-size(y,1),r,size(obj.A,3)));
            zeroCut     = @(x) x(1:obj.k,1:obj.k,:);

            if doGrad
                % [U,S,V,JacU,JacS,JacV,JacMU,JacMS,JacMV] = tSVDM(obj.A,M,obj.k);
                [AHat,~,JacMA1]     = modeProduct(obj.A,M);
                [U,S,V]             = facewiseSVD(AHat);
                [JacU2,JacS2,JacV2] = facewiseSVDJacobian(U,S,V);

                U = zeroColCut(U);
                S = zeroCut(S);
                V = zeroColCut(V);

                [B,JacU3,JacS3]     = facewise(U,S);
                [AkHat,JacB,JacV3]  = facewise(B,tran(V));
                [Ak,JacAk,JacMAk]   = modeProduct(AkHat,M,'invFlag',1,'orthoFlag',obj.orthoFlag);
            else
                AHat    = modeProduct(obj.A,M);
                [U,S,V] = facewiseSVD(AHat,obj.k);
                AkHat   = facewise(facewise(U,S),tran(V));
                Ak      = modeProduct(AkHat,M,'invFlag',1,'orthoFlag',obj.orthoFlag);
            end
            
            n = 1.0;
            if obj.avgFlag, n = size(obj.A,2); end

            R = obj.A - Ak;
            f = (0.5 / n) * fronorm(R).^2;
           
            % store information
            nrmR = fronorm(R);
            info.values = [nrmR, nrmR / obj.nrmA];

            if doGrad
                dR  = -(1 / n) * R;
                
                dfMinv  = JacMAk.AT(dR);
                dAkHat  = JacAk.AT(dR);
                dV3     = tran(JacV3.AT(dAkHat));
                dB3     = JacB.AT(dAkHat);
                dU3     = JacU3.AT(dB3);
                dS3     = JacS3.AT(dB3);


                % dU2 = facewise(facewise(dU3,V),S);
                dfU = JacU2.AT(zeroColPad(dU3));
               
                dfS = JacS2.AT(zeroPad(dS3));
                dfV = JacV2.AT(zeroColPad(dV3));

                dfAkHat = dfU + dfS + dfV;

                % dfAk = JacA1.AT(dfAkHat);
                dfM  = JacMA1.AT(dfAkHat);
                dfdM  = dfM + dfMinv;
                
%                 dU  = mprod(mprod(dR,V,M),tran(S),M);
%                 dfU = JacU.AT(dU);
% 
%                 dS  = mprod(mprod(tran(U),dR,M),V,M);
%                 dfS = JacS.AT(dS);
%                 
%                 dV  = mprod(mprod(tran(dR),U,M),S,M);
%                 dfV = JacV.AT(dV);
%
%                dfdM = dfU + dfS + dfV;
            end

        end

    end
    


end


