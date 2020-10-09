classdef SEXMOD < CATMOD
    properties
    end
    properties (Dependent = true)
       VarName;
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
    end
    methods % CONSTRUCTOR
        function obj = SEXMOD(varargin)
            obj = obj@CATMOD(varargin{:});         
        end
    end
    methods % GENERAL SETTING/GETTING
        function out = get.VarName(obj) %#ok<*MANU>
           out = 'SEX'; 
        end
        function out = get.VarType(obj)
           out = 'COVARIATE';
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
            out = 1;
        end
    end
    methods % BIOMETRICS
        function [TrX,TrRIP,TestX,TestRIP,posclass] = prepOppBiometricTest(obj,TestRIP,TestX)
            index = find(~isnan(TestX));
            TestRIP = TestRIP(index);TestX = TestX(index,:);
            TrX = obj.XMOD;
            TrRIP = obj.RIPMOD;
            posclass = obj.VarInfo.PosClass;
        end
    end
    methods % INTERFACING
        function out = getVarInfo(obj,DepVar)
            out.Range = 2;out.Mel = -1;
            out.MDepVar = mean(DepVar(obj.MInd,:),1);
            out.PDepVar = mean(DepVar(obj.FInd,:),1);
            switch obj.MFBal<=1
                case true
                    out.PosClass = -1;
                case false
                    out.PosClass = 1;
            end
        end
        function out = getFreq(obj,in) %#ok<*INUSD>
           out = 0; return;
           %switch in
           %    case -1
           %        out = length(obj.MInd)/(length(obj.FInd)+length(obj.MInd));
           %    case 1
           %        out = length(obj.FInd)/(length(obj.FInd)+length(obj.MInd));
           %    otherwise
           %        out = nan;
           %end
        end
        function [Distr,Freq,RIP] = X2RIP(obj,TrX,TrRIP,X)
            Distr = {};Freq = nan;RIP = nan;
            if isnan(X), return; end
            if isempty(obj.VarInfo), return; end
            Tind = find(TrX==X);% Inlier Distribution
            Find = find(TrX==-1*X);% Outlier Distribution
            [Distr.Tmu,Distr.Tsigma] = PPMMOD.getRIPStatistics(TrRIP(Tind)',obj.Kappa); %#ok<*FNDSB>
            RIP = Distr.Tmu;Distr.TpMax = normpdf(Distr.Tmu,Distr.Tmu,Distr.Tsigma);
            [Distr.Fmu,Distr.Fsigma] = PPMMOD.getRIPStatistics(TrRIP(Find)',obj.Kappa);
            Distr.FpMax = normpdf(Distr.Fmu,Distr.Fmu,Distr.Fsigma);
            Freq = getFreq(obj,X);
        end
        function postBuild(obj)
           return; 
        end
    end
    methods (Static = true)    
    end
end