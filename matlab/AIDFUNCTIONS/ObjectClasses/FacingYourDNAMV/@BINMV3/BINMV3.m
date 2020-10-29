classdef BINMV3 < CATMV3
    properties
       VarName;
    end
    properties (Dependent = true)
       VarType;
       SEX;
       MInd;
       FInd;
       MFBal;
    end
    properties (Hidden = true) 
    end
    properties(Hidden = true, Dependent = true)
        XREG;
        ID;
        PosClass;
    end
    methods % CONSTRUCTOR
        function obj = BINMV3(varargin)
            obj = obj@CATMV3(varargin{:});         
        end
    end
    methods % GENERAL SETTING/GETTING
        function out = get.VarName(obj) %#ok<*MANU>
           if isempty(obj.VarName), out = 'BIN'; return; end
           out = obj.VarName; 
        end
        function out = get.VarType(obj)
           out = 'BINARY VARIABLE';
        end
        function out = get.MInd(obj)
           if isempty(obj.X),out = [];return;end
           out = find(obj.X==-1);
        end
        function out = get.FInd(obj)
           if isempty(obj.X),out = [];return;end
           out = find(obj.X==1);
        end
        function out = get.MFBal(obj)
           if isempty(obj.X),out = [];return;end 
           out = length(obj.MInd)/length(obj.FInd); 
        end
        function out = get.SEX(obj)
           out = obj.X;
        end
        function set.SEX(obj,in)
           obj.X = in;
        end
        function out = get.XREG(obj)
            out = obj.X;
        end
        function out = get.ID(obj)
            out = 333;
        end
        function out = get.PosClass(obj)
            switch obj.MFBal<=1
                case true
                    out = -1;
                case false
                    out = 1;
                otherwise
                    out = [];
            end
        end
    end
    methods % CLASSIFIER
        function [X,COST,W,ind,Tind,Find] = prepTrainClassifier(obj,X) %#ok<*INUSL>
            ind = find(~isnan(X));
            Find = find(X==-1*obj.PosClass);
            Tind = find(X==obj.PosClass);
            COST.ClassNames = [obj.PosClass -1*obj.PosClass];
            COST.ClassificationCosts = zeros(2,2);
            COST.ClassificationCosts(1,2) = 1-(length(Tind)/length(ind));
            COST.ClassificationCosts(2,1) = 1-(length(Find)/length(ind));
            W = nan*zeros(size(X));
            W(Tind) = 1-(length(Tind)/length(ind));
            W(Find) = 1-(length(Find)/length(ind));
      end
      function [X,ind,Tind,Find,Tval] = prepTestClassifier(obj,X)
            ind = find(abs(X));
            Find = find(X==-1*obj.PosClass);
            Tind = find(X==obj.PosClass);
            Tval = obj.PosClass;
      end
    end
    methods % BIOMETRICS
        function out = getFreq(obj,in) %#ok<*INUSD>
           if isnan(in), out = nan; return; end
           %out = 0; return;
           switch in
              case -1
                  out = length(obj.MInd)/(length(obj.FInd)+length(obj.MInd));
              case 1
                  out = length(obj.FInd)/(length(obj.FInd)+length(obj.MInd));
              otherwise
                  out = nan;
           end
        end
    end
    methods % INTERFACING
    end
    methods (Static = true)    
    end
end