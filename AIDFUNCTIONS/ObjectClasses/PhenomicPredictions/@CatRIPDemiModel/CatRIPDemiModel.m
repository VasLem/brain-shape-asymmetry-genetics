classdef CatRIPDemiModel < RIPDemiModel
    properties
        ScoreMethod = 'z-score';% z-score or classify / threshold or logistic regression
        ClassifyMethod = 'classify';
        ZBoost = 1;
    end
    properties (Dependent = true)
       nrCat;
       CatLabel;
       CatMu;
       CatStd;
    end
    methods % Constructor
        function obj = CatRIPDemiModel(varargin)
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
            out = out(~isnan(out));
        end
        function out = get.CatMu(obj)
            if isempty(obj.EstRIP); out = []; return; end
            out = zeros(1,obj.nrCat);
            for i=1:1:obj.nrCat
                out(i) =  nanmean(obj.EstRIP(obj.TrX==obj.CatLabel(i)));
            end
        end
        function out = get.CatStd(obj)
             if isempty(obj.EstRIP); out = []; return; end
             out = zeros(1,obj.nrCat);
             for i=1:1:obj.nrCat
                out(i) =  nanstd(obj.EstRIP(obj.TrX==obj.CatLabel(i)));
             end
        end
    end
    methods % General Interface Functions
        function out = getScore(obj,in,val)
           out = getRIP(obj,in);
           switch obj.ScoreMethod
               case 'z-score'
                   mu = obj.CatMu(obj.CatLabel==val);
                   st = obj.CatStd(obj.CatLabel==val);
                   out = abs(obj.ZBoost*mu-out)/st;
               case 'classify'
                   switch obj.ClassifyMethod
                       case 'classify'
                           index = ~isnan(obj.TrX);
                           out = classify(out',obj.EstRIP(index)',obj.TrX(index),'quadratic')';
                           out = 1-(out==val);
                   end
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
                case 'classify'
                    errors = getScore(obj,obj.TrY',obj.TrX');
                    index = find(~isnan(errors));
                    obj.Acc = 100-((sum(errors(index))/length(index))*100);
            end
        end
        function deployModel(obj)
           trainModel(obj);
           obj.TrY = [];
           obj.TrRIP = [];
        end
    end
end