classdef SNPRIPDemiModel < RIPDemiModel
    % the difference is in the pairwise classification instead of
    % straigtforward group classification.
    properties
        ScoreMethod = 'z-score';% z-score or classify / threshold or logistic regression
        ClassifyMethod = 'ClassifyPW'; % ClassifyAll ClassifyPW
        Classifications = cell(1,3);
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
        function obj = SNPRIPDemiModel(varargin)
            obj = obj@RIPDemiModel(varargin{:});
            obj.XType = 'Categorical';
            obj.Classifications{1} = 'AABB';
            obj.Classifications{2} = 'AAAB';
            obj.Classifications{3} = 'BBAB';
        end
    end
    methods % Special Setting & Getting
        function out = get.nrCat(obj)
            out = length(obj.CatLabel);
        end
        function out = get.CatLabel(obj)
            out = unique(obj.TrX)';
            out = out(~isnan(out));
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
            %out(out<=0) = 0;
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
            out = zeros(1,obj.nrCat);
            for i=1:1:obj.nrCat
                out(i) = obj.TrInfoGain{i}.Ir;
            end        
        end
        function out = get.LOOIr(obj)
            if isempty(obj.LOOInfoGain), out = [];return; end
            out = zeros(1,obj.nrCat);
            for i=1:1:obj.nrCat
                out(i) = obj.LOOInfoGain{i}.Ir;
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
                       case 'ClassifyAll'
                           index = ~isnan(obj.TrX);
                           out = classify2use(out',obj.TrRIP(index)',obj.TrX(index));
                           out = 1-(out==val);
                       case 'ClassifyPW' % Naive Bayes classifier
                           index = find(obj.TrX~=0);
                           out1 = SNPRIPDemiModel.classify2use(out',obj.TrRIP(index),obj.TrX(index))';
                           out1 = (1-(out1==val))*obj.ClassAcc(1);
                           index = find(obj.TrX~=1);
                           out2 = SNPRIPDemiModel.classify2use(out',obj.TrRIP(index),obj.TrX(index))';
                           out2 = (1-(out2==val))*obj.ClassAcc(2);
                           index = find(obj.TrX~=-1);
                           out3 = SNPRIPDemiModel.classify2use(out',obj.TrRIP(index),obj.TrX(index))';
                           out3 = (1-(out3==val))*obj.ClassAcc(3);
                           out = ((val~=0).*out1+(val~=1).*out2+(val~=-1).*out3)/2;
                       case 'MisCost'  % Misclassification Cost
                            out1 = normpdf(out,obj.TrCatMu(1),obj.TrCatStd(1));
                            out2 = normpdf(out,obj.TrCatMu(2),obj.TrCatStd(2));
                            out3 = normpdf(out,obj.TrCatMu(3),obj.TrCatStd(3));
                            outsum = out1+out2+out3;
                            out = ((val~=-1).*out1+(val~=0).*out2+(val~=1).*out3)./outsum;
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
                    obj.ClassAcc = obj.Acc;
                case 'classify'
                    switch obj.ClassifyMethod
                        case 'ClassifyAll'
                            error('NOT IMPLEMENTED: PROBABLY OBSOLETE');
                            %errors = getScore(obj,obj.TrY',obj.TrX');
                            %index = find(~isnan(errors));
                            %obj.TrClassAcc = 100-((sum(errors(index))/length(index))*100);
                            %obj.TrAcc = obj.ClassAcc;
                            % TO BE REIMPLEMENTED, BUT PROBABLY OBSOLETE
                        case 'ClassifyPW'
                            obj.LOOInfoGain = cell(1,obj.nrCat);
                            obj.TrInfoGain = cell(1,obj.nrCat);
                            % -1 versus 1
                            index = find(obj.TrX~=0);
                            val = obj.TrX(index)';
                            out = SNPRIPDemiModel.classify2use(obj.EstRIP(index)',obj.TrRIP(index),val);
                            obj.LOOInfoGain{1} = classifierInfoGain(out,val');
                            scores = 1-(val==out);
                            obj.LOOClassAcc(1) = (1-sum(scores)/length(scores));
                            out = SNPRIPDemiModel.classify2use(obj.TrRIP(index),obj.TrRIP(index),val);
                            obj.TrInfoGain{1} = classifierInfoGain(out,val');
                            scores = 1-(val==out);
                            obj.TrClassAcc(1) = (1-sum(scores)/length(scores));
                            % -1 versus 0
                            index = find(obj.TrX~=1);
                            val = obj.TrX(index)';
                            out = SNPRIPDemiModel.classify2use(obj.EstRIP(index)',obj.TrRIP(index),val);
                            obj.LOOInfoGain{2} = classifierInfoGain(out,val');
                            scores = 1-(val==out);
                            obj.LOOClassAcc(2) = (1-sum(scores)/length(scores));
                            out = SNPRIPDemiModel.classify2use(obj.TrRIP(index),obj.TrRIP(index),val);
                            obj.TrInfoGain{2} = classifierInfoGain(out,val');
                            scores = 1-(val==out);
                            obj.TrClassAcc(2) = (1-sum(scores)/length(scores));
                            % 1 versus 0
                            index = find(obj.TrX~=-1);
                            val = obj.TrX(index)';
                            out = SNPRIPDemiModel.classify2use(obj.EstRIP(index)',obj.TrRIP(index),val);
                            obj.LOOInfoGain{3} = classifierInfoGain(out,val');
                            scores = 1-(val==out);
                            obj.LOOClassAcc(3) = (1-sum(scores)/length(scores));
                            out = SNPRIPDemiModel.classify2use(obj.TrRIP(index),obj.TrRIP(index),val);
                            obj.TrInfoGain{3} = classifierInfoGain(out,val');
                            scores = 1-(val==out);
                            obj.TrClassAcc(3) = (1-sum(scores)/length(scores));
                            obj.LOOClassAcc = (obj.LOOClassAcc-0.5)*2;% look at the fraction that is better than chance
                            obj.LOOClassAcc(obj.LOOClassAcc<0) = 0;
                            obj.TrClassAcc = (obj.TrClassAcc-0.5)*2;% look at the fraction that is better than chance
                            obj.TrClassAcc(obj.TrClassAcc<0) = 0;
                        case 'MisCost' % not needed
                            obj.LOOClassAcc = [100 100 100];
                            obj.TrClassAcc = [100 100 100];
                    end
            end
        end
        function test = getOutSampleAcc(obj,Y,X)
            test.acc = zeros(1,3);
            % observation and variables
            index = ~isnan(X);
            X = X(index);
            Y = Y(index,:);
            rip = getRIP(obj,Y');
            % Distribution statistics
            test.mu = zeros(1,obj.nrCat);
            test.std = zeros(1,obj.nrCat);
            test.distr = zeros(1,obj.nrCat);
            for i=1:1:obj.nrCat
                test.mu(i) =  nanmean(rip(X==obj.CatLabel(i)));
                test.std(i) =  nanstd(rip(X==obj.CatLabel(i)));
                test.distr(i) =  (sum(X==obj.CatLabel(i))/length(X))*100;
            end
            test.XRIPCor = permCorr(X,rip,0);
            % can we seperate any of the groups
            XG = cell(size(X));
            for i=1:1:obj.nrCat
                XG(X==obj.CatLabel(i)) = {num2str(i)};
            end
            [test.ANOVA.FA,test.ANOVA.pFA,test.ANOVA.ppermFA,test.ANOVA.T] = myPermAnova(XG,rip(:),10000);
            % -1 versus 1
            indexin = find(obj.TrX~=0);
            indexout = find(X~=0);
            out = SNPRIPDemiModel.classify2use(rip(indexout)',obj.TrRIP(indexin),obj.TrX(indexin))';
            test.InfoGain(1) = classifierInfoGain(out,X(indexout));
            tmp = (out==X(indexout));
            test.acc(1) = sum(tmp)/length(tmp);
            % -1 versus 0
            indexin = find(obj.TrX~=1);
            indexout = find(X~=1);
            out = SNPRIPDemiModel.classify2use(rip(indexout)',obj.TrRIP(indexin),obj.TrX(indexin))';
            test.InfoGain(2) = classifierInfoGain(out,X(indexout));
            tmp = (out==X(indexout));
            test.acc(2) = sum(tmp)/length(tmp);
            % 1 versus 0
            indexin = find(obj.TrX~=-1);
            indexout = find(X~=-1);
            out = SNPRIPDemiModel.classify2use(rip(indexout)',obj.TrRIP(indexin),obj.TrX(indexin))';
            test.InfoGain(3) = classifierInfoGain(out,X(indexout));
            tmp = (out==X(indexout));
            test.acc(3) = sum(tmp)/length(tmp);
            test.acc = (test.acc-0.5)*2;
            test.acc(test.acc<0) = 0;
            test.Ir = zeros(1,3);
            for i=1:1:obj.nrCat
               test.Ir(i) = test.InfoGain(i).Ir;
            end
            obj.OutSampleAcc = test;
        end
        function deployModel(obj)
           trainModel(obj);
           %obj.TrY = [];
           %obj.TrRIP = [];
        end
        function CatBoxPlot(obj)
           figure;boxplot(obj.EstRIP,obj.TrX);        
        end
    end
    methods (Static = true)
        function out = classify2use(TestFeatures,TrFeatures,TrClasses)
             %out = classify(TestFeatures,TrFeatures,TrClasses,'quadratic','empirical')';
             out = classify(TestFeatures,TrFeatures,TrClasses,'quadratic')';
        end
    end
end