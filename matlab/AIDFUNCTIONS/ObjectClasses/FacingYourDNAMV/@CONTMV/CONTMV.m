classdef CONTMV < PPMMV
    properties
        RegressorType = 'SVMRBF';
        X;
    end
    properties (Dependent = true)
        Predictor;
    end
    properties (Hidden = true)
        Regressor;
        FST=0.1;
        REGFST = 0.1;
        Sigma = 0.1;
        Kappa = 2;
        ON = 100;
        IBTH = 0.5;
        Bins = 10;
    end
    properties(Hidden = true, Dependent = true)
        XREG;
        treeTempl;
        predictortype;
    end
    methods % CONSTRUCTOR
        function obj = CONTMV(varargin)
            obj = obj@PPMMV(varargin{:});         
        end
    end
    methods % GENERAL SETTING/GETTING
        function out = get.XREG(obj)
            out = obj.X;
        end
        function out = get.Predictor(obj)
           out = obj.Regressor; 
        end
        function out = get.treeTempl(obj)
            if obj.nFeatures==0, out = []; return; end
            out = templateTree('MaxNumSplits',obj.nrSplits,'Type','Regression');
        end
        function out = get.predictortype(obj) %#ok<MANU>
           out = 'regressor';
        end
    end
    methods % PREDICTOR    
        function trainPredictor(obj,X,REGF,PCF,ADDF)
           % PREPARATION
             [X,ind,W] = prepPredictorTraining(obj,X);
             %W = ones(size(X));
           % FEATURE SELECTION
             obj.REGFInd = CONTMV.featureSelection(REGF,X,ind,obj.nLC,obj.REGFST); 
             obj.PCFInd = CONTMV.featureSelection(PCF,X,ind,obj.maxPCF,obj.FST);
             obj.ADDFInd = CONTMV.featureSelection(ADDF,X,ind,obj.maxADDF,obj.FST);
             if obj.nFeatures==0, return; end
           % ACCUMULATE FEATURES
             FS = accumulateFeatures(obj,REGF,PCF,ADDF);
           % TRAINING
            switch obj.RegressorType
                case 'SVM'
                    obj.Regressor = fitrsvm(FS(ind,:),X(ind),'Standardize',true,'Weights',W(ind));
                case 'SVMRBF'
                    obj.Regressor = fitrsvm(FS(ind,:),X(ind),'Standardize',true,'KernelFunction','RBF','KernelScale','auto','Weights',W(ind));
                case 'LSBOOST'
                    obj.Regressor = fitensemble(FS(ind,:),X(ind),'LSBoost',obj.nrTrees,obj.treeTempl,'Weights',W(ind));
                case 'BAG'
                    obj.Regressor = fitensemble(FS(ind,:),X(ind),'bag',obj.nrTrees,obj.treeTempl,'Type','Regression','Weights',W(ind));
                otherwise 
                    error('WRONG TYPE OF REGRESSOR');
            end
            learnSigma(obj,X,ind,REGF,PCF,ADDF);
        end
        function out = testPredictor(obj,X,REGF,PCF,ADDF,display)
            % PREPARATION
              [X,ind] = prepPredictorTesting(obj,X);
            % PREDICTION
              PX = predict(obj,REGF,PCF,ADDF);
            % EVALUATION
              tmp = corrcoef(X(ind),PX(ind));
              out.CORR = tmp(1,2);
              [~,HM] = Pred2Match(obj,X(ind),PX(ind));
              out.ACC = 1-sum(HM)/length(HM);
            % OPPOSITE AXIS TESTING
              X = X(ind);PX = PX(ind);
              nT = length(ind); 
              TMatches = nan*zeros(nT,1);
              FMatches = nan*zeros(nT,obj.ON);
              TClass = nan*zeros(nT,1);
              FClass = nan*zeros(nT,obj.ON);
              nrF = obj.ON;
              parfor t=1:1:nT
                 diff = abs(X-X(t));
                 [~,Find] = sort(diff);
                 Find = Find(end-nrF+1:end);
                 [MS,MH] = Pred2Match(obj,X(t),PX);
                 TMatches(t) = MS(t);
                 FMatches(t,:) = MS(Find);
                 TClass(t) = MH(t);
                 FClass(t,:) = MH(Find);
              end
              [out.EER,out.G,out.AUC,out.XX,out.YY,out.TH,out.Y] = PPMMV.getEER(TMatches,FMatches(:));% CHECK ORDER
              [out.HPREC,out.HREC,out.HG,out.HTNF,out.HX,out.HY,out.HACC] = PPMMV.getHardClass(1-TClass,1-FClass(:));          
            % VISUALIZATION
              if display
                 figure;hold on;grid on;
                 plot(X,PX,'b.');
                 title(['CORR: ' num2str(out.CORR)]);
                 figure;hold on;grid on;plot(0:0.1:1,0:0.1:1,'k-');
                 plot(out.XX,out.YY,'b-','LineWidth',1.5);
                 plot(out.HX,out.HY,'r.','MarkerSize',20);
                 title(['EER: ' num2str(out.EER) ' HG: ' num2str(out.HG)]);
              end
        end
        function [out,emt] = predict(obj,REGF,PCF,ADDF)
            % ACCUMULATE FEATURES
              FS = accumulateFeatures(obj,REGF,PCF,ADDF);
            % PREDICT  
              out = predict(obj.Predictor,FS);
              emt = [];
        end
        function learnSigma(obj,X,ind,REGF,PCF,ADDF)
              PX = predict(obj,REGF,PCF,ADDF);
              obj.Sigma = std(X(ind)-PX(ind));
              %obj.Sigma = obj.StdX;
        end
    end
    methods % BIOMETRICS
        function [SM,HM] = Pred2Match(obj,X,PX,varargin)
            if isscalar(X)
               if isnan(X), SM = nan*zeros(size(PX));HM = SM; return; end 
               X = repmat(X,size(PX));
            end
            IB = CONTMV.getInlierBelief(X-PX,obj.Sigma,obj.Kappa);
            HM = IB<obj.IBTH; % 0 = correct, 1 is not 
            SM = 1-IB;% The lower the better
            %SM = abs(X-PX)/obj.Sigma;
        end
        function out = prepBiometricTest(obj,X)
           out = X; 
       end
    end
    methods % INTERFACING
        function [X,ind,W] = prepPredictorTraining(obj,X)
            ind = find(~isnan(X));
            W = getW(obj,X);
        end
        function [X,ind] = prepPredictorTesting(obj,X)
            ind = find(~isnan(X));
        end
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
        function out = getFreq(obj,in) %#ok<*INUSL>
           out = 0;
           %out = nan*zeros(size(in));
           %ind =  find(~isnan(in));
           %if isempty(in), return; end
           %[N,x] = hist(obj.XMOD,min(min(obj.XMOD),min(in)):(max(obj.XMOD)-min(obj.XMOD))/obj.Bins:max(max(obj.XMOD),max(in)));
           %out(ind) = interp1(x,N,in(ind))/length(obj.XMOD);
        end
        function out = getW(obj,in)
           out = nan*zeros(size(in));
           ind =  find(~isnan(in));
           if isempty(in), return; end
           [N,x] = hist(in(ind),obj.Bins);
           out(ind) = 1-interp1(x,N,in(ind))/length(ind); 
        end
        function postFeatureLearning(obj) %#ok<MANU>
           return;
        end
    end
    methods (Static = true)
        function [out,COR] = featureSelection(FS,X,ind,maxF,FST)
             if isempty(FS), out = []; return; end
             nFS = size(FS,2);
             COR = nan*zeros(1,nFS);
             parfor f=1:nFS
                 tmp = corrcoef(X(ind),FS(ind,f)); %#ok<*PFBNS>
                 COR(f) = tmp(1,2);
             end
             COR = abs(COR);
             [sortedCOR,sortind] = sort(COR,'descend');
             featind = find(sortedCOR>=FST);
             if length(featind)>maxF,featind = featind(1:maxF);end
             out = sortind(featind);
             %figure;hist(COR);
        end
        function out = getInlierBelief(res,sigma,kappa)
                     IP = normpdf(res,0,sigma);
                     L = (1/(sqrt(2*pi)*sigma))*exp(-0.5*kappa^2);
                     out = IP./(IP+L);
        end
        function out = getInlierBeliefv2(res,sigma,kappa)
                     IP = normpdf(res,0,sigma);
                     out = IP./normpdf(0,0,sigma);
                     %L = (1/(sqrt(2*pi)*sigma))*exp(-0.5*kappa^2);
                     %out = IP./(IP+L);
        end
    end
end