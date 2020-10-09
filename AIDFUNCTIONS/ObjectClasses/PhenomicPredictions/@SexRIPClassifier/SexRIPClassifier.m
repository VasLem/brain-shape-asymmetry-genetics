classdef SexRIPClassifier < superClass
    % the difference is in the pairwise classification instead of
    % straigtforward group classification.
    properties
        TrX;
        TrY;
        TrC;
        M;
        Method = 'ROC';% ROC, Naief Bayes, ...
        ROCT;
    end
    properties (Dependent = true)
       nrCat;
       nrC;
       CatLabel;
       CatDistr;
    end
    methods % Constructor
        function obj = SexRIPClassifier(varargin)
            obj = obj@superClass(varargin{:});
        end
    end
    methods % Special Setting & Getting
        function out = get.nrCat(obj)
            out = length(obj.CatLabel);
        end
        function out = get.nrC(obj)
            if isempty(obj.TrC), out = 0; return; end
            out = size(obj.TrC,2);
        end
        function out = get.CatLabel(obj)
            out = unique(obj.TrX);
            out = out(~isnan(out))';
        end
        function out = get.CatDistr(obj)
            out = zeros(1,obj.nrCat);
            index = ~isnan(obj.TrX);
            X = obj.TrX(index);
            for i=1:1:obj.nrCat
               out(i) =  (sum(X==obj.CatLabel(i))/length(X))*100;
            end
        end
    end
    methods % General Interface Functions
        function trainModel(obj)
           createMEstRIP(obj); 
        end
        function deployModel(obj)
           trainModel(obj);
           %obj.TrY = [];
           %obj.TrRIP = [];
        end
        function createMEstRIP(obj)          
            if obj.nrC==0
               RIP = obj.TrRIP;
            else
               RIP = [obj.TrC obj.TrRIP];
               disp('With TrC');
            end
            [RIP,Y] = eliminateNAN(RIP,obj.TrY);
            [~,~,~,~,obj.TrM] = plsregress(RIP,Y,size(RIP,2));       
            TrY = obj.TrY;
            EstRIP = zeros(1,obj.nrTr); %#ok<*PROP>
            Ref = obj.YAvg;
            std = obj.YEigStd;
            nrTr = obj.nrTr;
            %tic;
            parfor i=1:nrTr
               build = setdiff((1:nrTr),i); 
               [forRIP,Y] = eliminateNAN(RIP(build,:),TrY(build,:));
               [~,~,~,~,M] = plsregress(forRIP,Y,size(forRIP,2));
               out = RIPDemiModel.getRIPStatic(Ref,std,TrY(i,:)',M);
               EstRIP(i) = out;
            end
            obj.EstRIP = EstRIP;
            
            % check if a sign flip occured
            if sign(obj.XTrCor)==sign(obj.XLOOCor), obj.Flipped = false; return; end
            disp('Flipping');
            obj.EstRIP = -1*obj.EstRIP;
            obj.Flipped = true;
            %toc;
            %obj.EstRIP = getRIP(obj,obj.TrY');
            %obj.EstRIP = obj.TrRIP';
        end
    end
    methods (Static = true)
        function out = classify2use(TestFeatures,TrFeatures,TrClasses)
             out = classify(TestFeatures,TrFeatures,TrClasses,'quadratic')';
        end
    end
end