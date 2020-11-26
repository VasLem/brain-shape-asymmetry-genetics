
classdef AsymmetryComponentsAnalysis
    properties
        n
        nrV
        rep
        dim=3
    end
    methods(Access=private)
        function X =reshape(obj, SSvec)
            X = reshape(SSvec,obj.dim,length(SSvec)/obj.dim);
        end
        function X =normalizeInput(obj, Input)
            if (length(size(Input))==1)
                X = obj.reshape(input);
            else
                X = Input;
            end
            
        end
    end
    methods
        function obj=AsymmetryComponentsAnalysis(n,nrV, rep)
            obj.n = n;
            obj.nrV = nrV;
            obj.rep = rep;
        end
        function [F, TF]=flunctuatingMS(obj, SS_F)
            F = sum(SS_F);
            TF = sum(F);
            F = F./(obj.dim*(obj.n-1));
            TF = TF/((obj.n-1)*(obj.nrV-7));
        end
        
        function [ E, TE ]=errorMS(obj, SS_E)
            SS_E = obj.normalizeInput(SS_E);
            E = sum(SS_E);
            TE = sum(E);
            E = E./((obj.rep-1)*obj.n*2*obj.dim);
            TE = TE/((obj.rep-1)*obj.n*2*(obj.nrV-7));
        end
        
        function [ I, TI ]=individualMS(obj, SS_I)
            SS_I = obj.normalizeInput(SS_I);
            I = sum(SS_I);
            TI = sum(I);
            I = I./(obj.dim*(obj.n-1));
            TI = TI/((obj.nrV-7)*(obj.n-1));
        end
        
        function [D, TD ]=directionMS(obj, SS_D)
            SS_D = obj.normalizeInput(SS_D);
            D = sum(SS_D);
            TD = sum(D);
            D = (D./obj.dim)/obj.n;
            TD = (TD/(obj.nrV-7))/obj.n;
        end
        
        function [I, TI, F, TF, IF, TIF] = individualEffect(obj, SS_I, SS_F)
            % Getting Fluctuating MS as error term
            [F, TF]= obj.flunctuatingMS(SS_F);
            % Getting Individual MS
            [I, TI] = obj.individualMS(SS_I);
            % Getting F statistic
            IF = I./F;
            TIF = TI/TF;
        end
        
        
        function [D, TD, F, TF, DF, TDF] = directionEffect(obj,SS_D, SS_F)
            % Getting Fluctuating MS as error term
            [F, TF] = obj.flunctuatingMS(SS_F);
            % Getting Direction MS
            [D, TD] = obj.directionMS(SS_D);
            % Getting F statistic
            DF = D./F;
            TDF = TD/TF;
        end
        
        
        function [E, TE, F, TF, FF, TFF] = interactionEffect(obj,SS_E, SS_F)
            % getting error MS as error term
            [E, TE] = obj.errorMS(SS_E);
            % getting Fluctuating MS
            [F, TF] = obj.flunctuatingMS(SS_F);
            % getting F-statistic
            FF = F./E;
            TFF = TF/TE;
        end
        
        function [LM, Total]= apply(obj, SS)
            % Error
            [LM.E, Total.E]= obj.errorMS(SS(4,:));
            
            % Fluctuating
            [LM.F, Total.F] = obj.flunctuatingMS(SS(3,:));
            
            
            % getting F-statistic
            LM.FF = LM.F./LM.E;
            Total.FF =  Total.F/Total.E;
            
            % Individuals
            [LM.I, Total.I] = obj.individualMS(SS(2,:));
            
            % getting F-statistic
            LM.IF = LM.I./LM.F;
            Total.IF =  Total.I/Total.F;
            
            % Directional
            [LM.D, Total.D] = obj.directionMS(SS(1,:));
            
            % getting F-statistic
            LM.DF = LM.D./LM.F;
            Total.DF =  Total.D/Total.F;
            
            LM.DP = Ftest(LM.DF,obj.dim,obj.dim*(obj.n-1));
            LM.FP = Ftest(LM.FF,obj.dim*(obj.n-1),((obj.rep-1)*obj.n*2*obj.dim));
            LM.IP = Ftest(LM.IF,(obj.dim*(obj.n-1)),(obj.dim*(obj.n-1)));
            
            Total.DP = Ftest(Total.DF,obj.nrV-7,(obj.n-1)*(obj.nrV-7));
            Total.FP = Ftest(Total.FF,(obj.n-1)*(obj.nrV-7),(obj.rep-1)*obj.n*2*(obj.nrV-7));
            Total.IP = Ftest(Total.IF,((obj.n-1)*(obj.nrV-7)),(obj.n-1)*(obj.nrV-7));
        end
    end
end
