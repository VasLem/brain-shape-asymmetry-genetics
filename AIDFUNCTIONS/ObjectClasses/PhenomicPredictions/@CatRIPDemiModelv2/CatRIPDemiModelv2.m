classdef CatRIPDemiModelv2 < RIPDemiModel
    % the difference is in the pairwise classification instead of
    % straigtforward group classification.
    properties
        ScoreMethod = 'z-score';% z-score or classify / threshold or logistic regression
        ZBoost = 1;
        TrClassAcc;
        LOOClassAcc;
        TrInfoGain;
        LOOInfoGain;
    end
    properties (Dependent = true)
       ClassAcc;
       nrCat;
       CatLabel;
       CatDistr;
       LOOCatMu;
       LOOCatStd;
       TrCatMu;
       TrCatStd;
       nr2Class;
       maxScore;
       TrIr;
       LOOIr;
    end
    methods % Constructor
        function obj = CatRIPDemiModelv2(varargin)
            obj = obj@RIPDemiModel(varargin{:});
            obj.XType = 'Categorical';
        end
    end
    methods % Special Setting & Getting
        function out = get.nrCat(obj)
            out = length(obj.CatLabel);
        end
        function out = get.CatLabel(obj)
            out = unique(obj.TrX);
            out = out(~isnan(out))';
        end
        function out = get.LOOCatMu(obj)
            if isempty(obj.EstRIP); out = []; return; end
            out = zeros(1,obj.nrCat);
            for i=1:1:obj.nrCat
                out(i) =  nanmean(obj.EstRIP(obj.TrX==obj.CatLabel(i)));
            end
        end
        function out = get.LOOCatStd(obj)
             if isempty(obj.EstRIP); out = []; return; end
             out = zeros(1,obj.nrCat);
             for i=1:1:obj.nrCat
                out(i) =  nanstd(obj.EstRIP(obj.TrX==obj.CatLabel(i)));
             end
        end
        function out = get.TrCatMu(obj)
            if isempty(obj.TrRIP); out = []; return; end
            out = zeros(1,obj.nrCat);
            for i=1:1:obj.nrCat
                out(i) =  nanmean(obj.TrRIP(obj.TrX==obj.CatLabel(i)));
            end
        end
        function out = get.TrCatStd(obj)
             if isempty(obj.TrRIP); out = []; return; end
             out = zeros(1,obj.nrCat);
             for i=1:1:obj.nrCat
                out(i) =  nanstd(obj.TrRIP(obj.TrX==obj.CatLabel(i)));
             end
        end
        function out = get.nr2Class(obj)
           out = obj.nrCat-1; 
        end
        function out = get.maxScore(obj)
            out = obj.nr2Class;
        end
        function out = get.ClassAcc(obj)
            out = obj.LOOClassAcc;
        end
        function out = get.CatDistr(obj)
            out = zeros(1,obj.nrCat);
            index = ~isnan(obj.TrX);
            X = obj.TrX(index);
            for i=1:1:obj.nrCat
               out(i) =  (sum(X==obj.CatLabel(i))/length(X))*100;
            end
        end
        function out = get.TrIr(obj)
            if isempty(obj.TrInfoGain), out = [];return; end
            out = obj.TrInfoGain.Ir;        
        end
        function out = get.LOOIr(obj)
            if isempty(obj.LOOInfoGain), out = [];return; end
            out = obj.LOOInfoGain.Ir;     
        end
    end
    methods % General Interface Functions
        function out = getScore(obj,in,val)
           out = getRIP(obj,in);
           switch obj.ScoreMethod
               case 'z-score'
                   mu = obj.TrCatMu(obj.CatLabel==val);
                   st = obj.TrCatStd(obj.CatLabel==val);
                   out = abs(obj.ZBoost*mu-out)/st;
               case 'classify'
                   index = ~isnan(obj.TrX);
                   out = CatRIPDemiModelv2.classify2use(out',obj.TrRIP(index),obj.TrX(index));
                   out = 1-(out==val);
               case 'MisCost'
                   out1 = normpdf(out,obj.TrCatMu(1),obj.TrCatStd(1));
                   out2 = normpdf(out,obj.TrCatMu(2),obj.TrCatStd(2));
                   outsum = out1+out2;
                   out = ((val~=-1).*out1+(val~=1).*out2)./outsum;     
           end
        end
        function trainModel(obj)
           createMEstRIP(obj);
           getClassAcc(obj); 
        end
        function getClassAcc(obj)
            switch obj.ScoreMethod
                case 'z-score'
                    getAcc(obj);
                    obj.LOOClassAcc = obj.Acc;
                case 'classify'
                    index = ~isnan(obj.TrX);
                    out = CatRIPDemiModelv2.classify2use(obj.EstRIP(index)',obj.TrRIP(index),obj.TrX(index));
                    obj.LOOInfoGain = classifierInfoGain(out,obj.TrX(index));
                    out = 1-(out==obj.TrX(index)');
                    obj.LOOClassAcc = (1-sum(out)/length(out))*100;
                    out = CatRIPDemiModelv2.classify2use(obj.TrRIP(index),obj.TrRIP(index),obj.TrX(index));
                    obj.TrInfoGain = classifierInfoGain(out,obj.TrX(index));
                    out = 1-(out==obj.TrX(index)');
                    obj.TrClassAcc = (1-sum(out)/length(out))*100;
                    obj.Acc = obj.ClassAcc;
                case 'MisCost' % not needed
                    obj.LOOClassAcc = [100 100 100];
                    obj.TrClassAcc = [100 100 100];     
            end
        end
        function out = getOutSampleAcc(obj,Y,X)
            % observation and variables
            index = ~isnan(X);
            X = X(index);
            Y = Y(index,:);
            rip = getRIP(obj,Y');
            index = ~isnan(obj.TrX);
            err = CatRIPDemiModelv2.classify2use(rip',obj.TrRIP(index),obj.TrX(index));
            out.InfoGain = classifierInfoGain(err,X);
            err = 1-(err==X);
            err = (1-sum(err)/length(err))*100;
            out.err = err;
            out.mu = zeros(1,obj.nrCat);
            out.std = zeros(1,obj.nrCat);
            out.distr = zeros(1,obj.nrCat);
            for i=1:1:obj.nrCat
                out.mu(i) =  nanmean(rip(X==obj.CatLabel(i)));
                out.std(i) =  nanstd(rip(X==obj.CatLabel(i)));
                out.distr(i) =  (sum(X==obj.CatLabel(i))/length(X))*100;
            end
            out.XRIPCor = permCorr(X,rip,0);
            XG = cell(size(X));
            for i=1:1:obj.nrCat
                XG(X==obj.CatLabel(i)) = {num2str(i)};
            end
            [out.ANOVA.FA,out.ANOVA.pFA,out.ANOVA.ppermFA,out.ANOVA.T] = myPermAnova(XG,rip(:),10000);
            obj.OutSampleAcc = out;
            % is the model able to seperate groups
        end
        function deployModel(obj)
           trainModel(obj);
           %obj.TrY = [];
           %obj.TrRIP = [];
        end
    end
    methods (Static = true)
        function out = classify2use(TestFeatures,TrFeatures,TrClasses)
             out = classify(TestFeatures,TrFeatures,TrClasses,'quadratic')';
        end
    end
end