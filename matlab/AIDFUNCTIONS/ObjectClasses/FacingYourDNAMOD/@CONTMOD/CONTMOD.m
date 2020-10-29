classdef CONTMOD < PPMMOD
    properties
        KN = 300;
        ON = 30;
    end
    properties (Dependent = true)
    end
    properties (Hidden = true) 
    end
    properties(Hidden = true, Dependent = true)
        XREG;
    end
    methods % CONSTRUCTOR
        function obj = CONTMOD(varargin)
            obj = obj@PPMMOD(varargin{:});         
        end
    end
    methods % GENERAL SETTING/GETTING
        function out = get.XREG(obj)
            out = obj.X;
        end
    end
    methods % BIOMETRICS
        function [TrX,TrRIP,TestX,TestRIP,posclass] = prepOppBiometricTest(obj,TestRIP,TestX)
            index = find(~isnan(TestX));
            TestRIP = TestRIP(index);
            TestX = TestX(index,:);
            TrX = obj.XMOD;
            TrRIP = obj.RIPMOD;
            posclass = 0;% redundant
        end
        function [EER,G,AUC,R] = oppBiometricTest(obj,TrX,TrRIP,TestX,TestRIP,posclass) %#ok<INUSD>
                nT = length(TestRIP);
                TMatches = nan*zeros(nT,1);
                FMatches = nan*zeros(nT,obj.ON);
                nrF = obj.ON;
                parfor t=1:1:nT
                    diff = abs(TestX-TestX(t));
                    [~,Find] = sort(diff);
                    Find = Find(end-nrF+1:end);
                    Distr = X2RIP(obj,TrX,TrRIP,TestX(t));
                    m = getMatch(obj,Distr,TestRIP);
                    TMatches(t) = m(t);
                    FMatches(t,:) = m(Find);
                end
                TScores = TMatches;FScores = FMatches;
                switch obj.Classify
                   case 'SOFT'
                       TScores = match2score(obj,TMatches);
                       FScores = match2score(obj,FMatches);
                   case 'HARD'
                       TScores = PPMMOD.class2score(match2class(obj,TMatches));
                       FScores = PPMMOD.class2score(match2class(obj,FMatches));
                end
                [EER,G,AUC] = PPMMOD.getEER(TScores,FScores(:));
                R = PPMMOD.getRANK(TScores,FScores(:));
        end
    end
    methods % INTERFACING
        function out = getVarInfo(obj,DepVar)
            [sortv,ind] = sort(obj.X(~isnan(obj.X)),'ascend');
           % lower 25%
            indl = round(0.25*length(ind));
            out.Mb = sortv(indl);
            Mindex = find(obj.X<=out.Mb);
            out.Mel = mean(obj.X(Mindex));
            out.Mrange = [min(obj.X(Mindex)) max(obj.X(Mindex))];
            out.MDepVar = mean(DepVar(Mindex,:));
           % upper 25%
            indu = round(0.75*length(ind));
            out.Pb = sortv(indu);
            Pindex = find(obj.X>=out.Pb);
            out.Pel = mean(obj.X(Pindex));
            out.Prange = [min(obj.X(Pindex)) max(obj.X(Pindex))];
            out.PDepVar = mean(DepVar(Pindex,:));
            out.Range = out.Pel-out.Mel;
        end
        function out = getFreq(obj,in)
           if isnan(in), out = nan; return; end
           out = 0;
           %[N,X] = hist(obj.XMOD);
           %out = interp1(X,N,in)/length(obj.XMOD);
        end
        function [Distr,Freq,RIP] = X2RIP(obj,TrX,TrRIP,X)
           Distr = {};Freq = nan;RIP = nan;
           if isnan(X), return; end
           if isempty(obj.VarInfo), return; end
           % eliminate nans 
            ind = find(~isnan(TrX));
            TrX = TrX(ind);
           % select closest samples
            W = abs(TrX-X);
            [~,ind2] = sort(W,'ascend');
           % TRUE DISTRIBUTION
            ind2T = ind2(1:obj.KN);indT = ind(ind2T);
            rip = TrRIP(indT)';
            WT = W(ind2T);WT = WT-min(WT);WT = WT/max(WT);WT = 1-WT';
            [Distr.Tmu,Distr.Tsigma] = PPMMOD.getWeightedRIPStatistics(rip,WT,obj.Kappa);
            Distr.TpMax = normpdf(Distr.Tmu,Distr.Tmu,Distr.Tsigma);
            RIP = Distr.Tmu;
           % FALSE DISTRIBUTION
            ind2F = ind2(end-obj.KN+1:end);indF = ind(ind2F);
            rip = TrRIP(indF)';
            WF = W(ind2F);WF = WF-min(WF);WF = WF/max(WF);WF = WF';
            [Distr.Fmu,Distr.Fsigma] = PPMMOD.getWeightedRIPStatistics(rip,WF,obj.Kappa);
            Distr.FpMax = normpdf(Distr.Fmu,Distr.Fmu,Distr.Fsigma);
            Freq = getFreq(obj,X);
        end
        function postBuild(obj) %#ok<MANU>
           return;
        end
    end
    methods (Static = true)    
    end
end