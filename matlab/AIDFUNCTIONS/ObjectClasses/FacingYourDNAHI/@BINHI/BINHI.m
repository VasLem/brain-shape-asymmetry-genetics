classdef BINHI < PPMHI
    properties
       VarName;
    end
    properties (Dependent = true)
       BV; 
       Labels;
       mLabel;
       MLabel;
       VarType;
       MinorityClass;
    end
    properties (Hidden = true)
        NanVal = nan;
    end
    properties(Hidden = true, Dependent = true)
        ID;
    end
    methods % CONSTRUCTOR
        function obj = BINHI(varargin)
            obj = obj@PPMHI(varargin{:});         
        end
    end
    methods % GENERAL SETTING/GETTING
        function set.BV(obj,in)
                 lab = unique(in(~isnan(in)));
                 if ~(length(lab)==2), error('TWO CLASSES ONLY');end
                 obj.XX = in;
        end
        function out = get.BV(obj)
                 out = obj.XX;
        end
        function out = get.Labels(obj)
                 out = unique(obj.BV(~isnan(obj.BV)))';
        end
        function out = get.VarName(obj) %#ok<*MANU>
           if isempty(obj.VarName), out = 'BIN'; return; end
           out = obj.VarName; 
        end
        function out = get.VarType(obj)
           out = 'BINARY VARIABLE';
        end
        function out = get.ID(obj)
            out = 333;
        end
        function out = get.mLabel(obj)
              tmp = obj.BV;
              lab = obj.Labels;
              one = sum(tmp==lab(1));
              two = sum(tmp==lab(2));
              switch one<=two
                  case true
                      out = lab(1);
                  case false
                      out = lab(2);
              end
        end
        function out = get.MLabel(obj)
              out = setdiff(obj.Labels,obj.mLabel);
        end
        function out = get.MinorityClass(obj)
                 out = 1;
        end
    end
    methods % CLASSIFIER
        function [G,COST,W,ind,Tind,Find] = prepTrainSVM(obj,X) %#ok<*INUSL>
            G = XX2GG(obj,X);
            ind = find(~(X==obj.NanVal));
            Find = find(G==-1*obj.MinorityClass);
            Tind = find(G==obj.MinorityClass);
            COST.ClassNames = [obj.MinorityClass -1*obj.MinorityClass];
            COST.ClassificationCosts = zeros(2,2);
            COST.ClassificationCosts(1,2) = 1-(length(Tind)/length(ind));
            COST.ClassificationCosts(2,1) = 1-(length(Find)/length(ind));
            W = nan*zeros(size(X));
            W(Tind) = 1-(length(Tind)/length(ind));
            W(Find) = 1-(length(Find)/length(ind));
       end
    end
    methods % BIOMETRICS
    end
    methods % INTERFACING
        function out = XX2XG(obj,in)
                 out = in;
                 out(in==obj.mLabel) = 1;% always code the minority class equal to plus one;
                 out(in==obj.MLabel) = -1;% always code the majority class equal to minus one;
        end
        function out = XG2GG(obj,in)
                 out = in;
        end
    end
    methods (Static = true)    
    end
end