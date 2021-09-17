
classdef AsymmetryComponentsAnalysis
    properties
        n
        nrV
        rep
        dim=3
    end
    methods(Access=private)
        function X =normalizeInput(obj, Input)
            X = reshape(Input,obj.dim,length(Input)/obj.dim);
        end
    end
    methods
        
        function obj=AsymmetryComponentsAnalysis(n,nrV, rep)
            obj.n = n;
            obj.nrV = nrV;
            obj.rep = rep;
        end
        
        function dof = getFluctuatingDOF(obj, total)
            if ~total
                dof =((obj.n-1)*obj.dim);
            else
                dof = ((obj.n-1)*(obj.nrV-7));
            end
        end
        function dof = getErrorDOF(obj, total)
            if ~total
                dof = ((obj.rep-1)*obj.n*2*obj.dim);
            else
                dof = ((obj.rep-1)*obj.n*2*(obj.nrV-7));
            end
        end
        function dof = getIndividualDOF(obj, total)
            if ~total
                dof =((obj.n-1)*obj.dim);
            else
                dof = ((obj.n-1)*(obj.nrV-7));
            end
        end
        
        function dof = getDirectionDOF(obj, total)
            if ~total
                dof =(obj.dim);
            else
                dof = ((obj.nrV-7));
            end
        end
        
        
        function [F, TF]=fluctuatingMS(obj, SS_F)
            SS_F = obj.normalizeInput(SS_F);
            F = sum(SS_F);
            TF = sum(F);
            F = F./obj.getFluctuatingDOF(false);
            TF = TF/obj.getFluctuatingDOF(true);
        end
        
        function [ E, TE ]=errorMS(obj, SS_E)
            SS_E = obj.normalizeInput(SS_E);
            E = sum(SS_E);
            TE = sum(E);
            E = E./ obj.getErrorDOF(false);
            TE = TE/ obj.getErrorDOF(true);
        end
        
        function [ I, TI ]=individualMS(obj, SS_I)
            SS_I = obj.normalizeInput(SS_I);
            I = sum(SS_I);
            TI = sum(I);
            I = I./obj.getIndividualDOF(false);
            TI = TI/obj.getIndividualDOF(true);
        end
        
        function [D, TD ]=directionMS(obj, SS_D)
            SS_D = obj.normalizeInput(SS_D);
            D = sum(SS_D);
            TD = sum(D);
            D = D./obj.getDirectionDOF(false);
            TD = TD/obj.getDirectionDOF(true);
        end
        
        function [I, TI, F, TF, IF, TIF] = individualEffect(obj, SS_I, SS_F)
            % Getting Fluctuating MS as error term
            [F, TF]= obj.fluctuatingMS(SS_F);
            % Getting Individual MS
            [I, TI] = obj.individualMS(SS_I);
            % Getting F statistic
            IF = I./F;
            TIF = TI/TF;
        end
        
        function [D, TD, F, TF, DF, TDF] = directionEffect(obj,SS_D, SS_F)
            % Getting Fluctuating MS as error term
            [F, TF] = obj.fluctuatingMS(SS_F);
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
            [F, TF] = obj.fluctuatingMS(SS_F);
            % getting F-statistic
            FF = F./E;
            TFF = TF/TE;
            
        end
        
        function displayResultInScatter3D(Landmarks, result)
            shape = shape3D;
            shape.Vertices = Landmarks;
            shape.VertexSize = 20;
            shape.VertexValue = result;
            shape.S = [0 1 0];
            viewer(shape);
            drawnow;
        end
        
        function ret = shuffleColumnWise(obj, Set1, Set2)
            SSet1 = Set1;
            SSet2 = Set2;
            r = randi(2,obj.n,1);
            columnPermIndex = find(r==2);
            SSet1(:, columnPermIndex) = Set2(:, columnPermIndex);
            SSet2(:, columnPermIndex) = Set1(:, columnPermIndex);
            ret  = permute(cat(3, SSet1(:,:), SSet2(:, :)), [2,3,1]);
        end
        
        function ret = shuffleRowWise(obj, Set1, Set2)
            rowPermIndex = randperm(obj.n);
            SSet1 = Set1(:,:,rowPermIndex);
            ret =permute(cat(3, SSet1(:,:), Set2(:, :)), [2,3,1]);
        end
        
        function ret = shuffleResidual(obj, Set1, Set2)
            ret = permute(cat(3, Set1(:,:), Set2(:, :)), [2,3,1]);
            
            avgC = mean(ret, 1);
            avgR = mean(ret,2);
            avg = mean(avgC);
            ret = ret - avgC - avgR + avg;
            index = randperm(obj.n*obj.rep*2);
            ret = reshape(ret, obj.n*obj.rep * 2, obj.nrV);
            ret = reshape(ret(index,:),obj.n*obj.rep,2,obj.nrV);
        end

        
        function [LM, Total]= apply(obj, SS)
            % Error
            [LM.E, Total.E]= obj.errorMS(SS(4,:));
            [LM.F, Total.F] = obj.fluctuatingMS(SS(3,:));
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
            
            LM.DP = Ftest(LM.DF,obj.getDirectionDOF(false),obj.getFluctuatingDOF(false));
            LM.FP = Ftest(LM.FF,obj.getFluctuatingDOF(false),obj.getErrorDOF(false));
            LM.IP = Ftest(LM.IF,obj.getIndividualDOF(false),obj.getFluctuatingDOF(false));
            
            Total.DP = Ftest(Total.DF,obj.getDirectionDOF(true),obj.getFluctuatingDOF(true));
            Total.FP = Ftest(Total.FF,obj.getFluctuatingDOF(true),obj.getErrorDOF(true));
            Total.IP = Ftest(Total.IF,obj.getIndividualDOF(true),obj.getFluctuatingDOF(true));
        end
    end
end
