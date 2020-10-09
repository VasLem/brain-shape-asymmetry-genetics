classdef CATMV3TEST < PPMMV3TEST
    % GENERAL PROPERTIES
    properties
        X;
    end
    properties (Dependent = true)
    end
    properties (Hidden = true)
        
    end
    properties(Hidden = true, Dependent = true)
    end
    % CLASSIFIER
    properties
        Classifier = [];
        Standardize = true;
        Optimize = false;
    end
    properties (Dependent = true)
    end
    properties (Hidden = true)
        FSIndex = [];
        nTraining = 0;
        CV;
        OffsetPXC;
        ScalePXC;
        AvgFS;
        StdFS;
    end
    properties(Hidden = true, Dependent = true)
    end
    % ABSTRACT PROPERTIES
    properties (Abstract = true)
        PosClass; 
    end
    methods % CONSTRUCTOR
        function obj = CATMV3TEST(varargin)
            obj = obj@PPMMV3TEST(varargin{:});         
        end
    end
    methods % GENERAL SETTING/GETTING
        function obj = set.X(obj,in)
           % Conversion for memory saving 
           %ind = isnan(in); in = int8(in);in(ind) = -128;
           if strcmp(obj.Type,'SNPMV'), obj.X = single(in);return; end
           if strcmp(obj.Type,'SEXMV'), obj.X = single(in);return; end
           obj.X = in;
        end
    end
    methods % CLASSIFIER ANALYSIS
        function trainClassifier(obj,X,FS,varargin)
            % READING INPUT
              Input = find(strcmpi(varargin, 'K'));
              if isempty(Input), K = 4; else K = varargin{Input+1};end
              Input = find(strcmpi(varargin, 'optimize'));
              if isempty(Input), obj.Optimize = false; else obj.Optimize = varargin{Input+1};end
              Input = find(strcmpi(varargin, 'standardize'));
              if isempty(Input), obj.Standardize = true; else obj.Standardize = varargin{Input+1};end
            % PREPARATION
              [X,COST,W,ind,~,~] = prepTrainClassifier(obj,X);
              obj.nTraining = length(ind);
            % SVM LEARNING
              if ~isobject(K), rng(1);K = cvpartition(length(ind),'KFold',K);end
              out = CATMV3TEST.trainSVMRBF(FS(ind,:),X(ind),COST,W(ind),K,obj.Optimize,obj.Standardize,obj.PosClass);
              obj.Classifier = out.Classifier;
              obj.CV.Loss = out.cv;
              obj.CV.AUC = out.AUC;
              obj.CV.pAUC = out.pAUC;
              obj.OffsetPXC = out.OffsetPXC;
              obj.ScalePXC = out.ScalePXC;
              obj.AvgFS = out.AvgFS;
              obj.StdFS = out.StdFS;
        end
        function out = testClassifier(obj,X,FS,varargin)
            % READING INPUT
              Input = find(strcmpi(varargin, 'display'));
              if isempty(Input), display = false; else display = varargin{Input+1};end
              if isempty(obj.Classifier), warning('NO CLASSIFIER'); return;end
            % PREPARATION
              [~,~,Tind,Find,Tval] = prepTestClassifier(obj,X);
            % PREDICTION
              [PX,~,PXCN] = predict(obj,FS,length(X));   
              switch Tval
                  case -1
                      PXCN = PXCN(:,2);
                  case 1
                      PXCN = PXCN(:,1);
              end
            % EVALUATION
              [out.EER,out.G,out.AUC,out.XX,out.YY,out.TH,out.Y,out.pAUC] = CATMV3TEST.getEER(PXCN(Tind),PXCN(Find));
              res = PX==Tval;
              [out.HPREC,out.HREC,out.HG,out.HTNF,out.HX,out.HY,out.HACC] = CATMV3TEST.getHardClass(res(Tind),res(Find));
            % VISUALIZATION
              if display
                  figure;hold on;grid on;grid minor;plot(0:0.1:1,0:0.1:1,'k-');
                  plot(0:0.1:1,1:-0.1:0,'k--');
                  plot(out.XX,out.YY,'b-','LineWidth',1.5);
                  plot(out.HX,out.HY,'r.','MarkerSize',20);
                  xlabel('FPR');ylabel('TPR');
                  title(['G: ' num2str(out.HG) ' AUC: ' num2str(out.AUC) ' pAUC: ' num2str(out.pAUC)]);
                  drawnow;
              end
        end
        function [PX,PXC,PXCN] = predict(obj,FS,nS)
            % ACCUMULATE FEATURES
              FS = getFS(obj,FS);
            % STANDARDIZE FEATURES
              if obj.Standardize
                 FS = FS-repmat(obj.AvgFS,nS,1); 
                 FS = FS./repmat(obj.StdFS,nS,1); 
              end
            % PREDICT  
              [PX,PXC] = predict(obj.Classifier,FS);
              PXCN = CATMV3TEST.statNormPXC(PXC,obj.OffsetPXC,obj.ScalePXC);
        end
        function out = getFS(obj,FS)
            if ~iscell(FS), out = FS; return; end
            out = getCellData(obj.HI,FS,obj.FSIndex); 
        end
        function getFSIndex(obj,pCrit,val)
           if nargin<3, val = obj.pT;end
           if nargin<2, pCrit = 1;end
           out = [];
           for ind=1:1:obj.nLC
               % ind=1;
               if val(ind)>pCrit, continue;end
               [lev,cl] = getLevClust(obj,ind);
%                list = getAllChildren(obj,lev,cl);
%                test = sum(val(list)<val(ind));
               list = getChildren(obj,lev,cl);
               if ~isempty(list)
                  list = [getClind(obj,lev+1,list(1)) getClind(obj,lev+1,list(2))];
                  test = sum(val(list)<val(ind));
               else
                  test = 0;
               end
               if test>0, continue;end
               out = [out ind getAllChildren(obj,lev,cl)]; %#ok<AGROW>
               %Children = getChildren(obj,lev,cl);
               %if isempty(Children),continue;end
               %out = [out getClind(obj,lev+1,Children(1)) getClind(obj,lev+1,Children(2))];
           end
           out = unique(out);
           obj.FSIndex = out;
        end
    end
    methods % BIOMETRICS
    end
    methods % EVOMORPH
    end
    methods % INTERFACING
    end
    methods (Static = true)
        function out = trainSVMRBF(FS,X,COST,W,K,opt,Stand,Tval)
            %FS = FS(ind,:);X = X(ind);W = W(ind);Tval = obj.PosClass;
            if Stand
               out.Standardize = true;
               [FS,out.AvgFS,out.StdFS] = PPMMV3TEST.standardize(FS);
            else
               out.Standarize = false;
               out.AvgFS = [];
               out.StdFS = [];
            end
            if opt
               COST.ClassificationCosts(1,2) = 1;
               COST.ClassificationCosts(2,1) = 0.5;
            end
            %COST.ClassificationCosts(1,2) = 1;
            %COST.ClassificationCosts(2,1) = 1;
            C = 1;
            tmp = fitcsvm(FS,X,'Standardize',false,'KernelFunction','rbf','BoxConstraint',C,...
                          'KernelScale','auto','Cost',COST,'Weights',W,'Prior','empirical');
            sigma = tmp.KernelParameters.Scale;          
            if ~isobject(K), K = cvpartition(length(X),'KFold',K); end          
            switch opt
                case false
                    out.Classifier = fitcsvm(FS,X,'Standardize',false,'KernelFunction','rbf','BoxConstraint',C,...
                                             'KernelScale',sigma,'Cost',COST,'Weights',W,'Prior','empirical');                                       
                case true
                    init.scale = log(sigma);
                    init.c = log(1);
                    [out.cv,sigma,C] = CATMV3TEST.optSVMRBF(FS,X,COST,W,K,init);
                    out.Classifier = fitcsvm(FS,X,'Standardize',false,'KernelFunction','rbf','BoxConstraint',C,...
                                             'KernelScale',sigma,'Cost',COST,'Weights',W,'Prior','empirical');                     
            end
            % CROSS VALIDATION LOSS
             out.cv = kfoldLoss(fitcsvm(FS,X,'CVPartition',K,...
                               'KernelFunction','rbf','BoxConstraint',C,...
                               'Cost',COST,'Weights',W,'Prior','empirical',...
                               'KernelScale',sigma));
            % CROSS VALIDATION EER
             [out.AUC,out.pAUC] = CATMV3TEST.kfoldEER(FS,X,K,COST,W,C,sigma,Tval,false);
            % PREDICTION OFFSET AND SCALE LEARNING  
             [~,PXC] = predict(out.Classifier,FS);
             out.OffsetPXC  = min(PXC);
             out.ScalePXC = max(PXC)-min(PXC);
        end
        function [fval,sigma,C] = optSVMRBF(FS,X,COST,W,K,init)
            if nargin<5, K = 4; end
            if ~isobject(K), K = cvpartition(length(X),'KFold',K); end
            minfn = @(z)kfoldLoss(fitcsvm(FS,X,'CVPartition',K,'Standardize',false,...
                       'KernelFunction','rbf','BoxConstraint',exp(z(2)),...
                       'Cost',COST,'Weights',W,...
                       'KernelScale',exp(z(1)),'Prior','empirical'));
            options = psoptimset('UseParallel', true, 'CompletePoll', 'on', 'Vectorized', 'off','TolMesh',5e-5,'Display','off');
            [searchmin, fval] = patternsearch(minfn,[init.scale, init.c],[],[],[],[],-5*ones(2,1),5*ones(2,1),options);
            z = exp(searchmin);
            sigma = z(1);
            C = z(2);
        end
        function [auc,pauc] = kfoldEER(fs,x,K,COST,w,C,sigma,Tval,display)
                 if nargin < 9, display = false; end
                 auc = zeros(1,K.NumTestSets);
                 pauc = zeros(1,K.NumTestSets);
                 if display
                   figure;hold on;grid on;plot(0:0.1:1,0:0.1:1,'k-');
                   plot(0:0.1:1,1:-0.1:0,'k--');
                 end
                 for i=1:K.NumTestSets
                     TrInd = K.training(i);
                     TestInd = K.test(i);
                     % BUILD CLASSIFIER
                     Classifier = fitcsvm(fs(TrInd,:),x(TrInd),'Standardize',false,'KernelFunction','rbf','BoxConstraint',C,...
                                          'KernelScale',sigma,'Cost',COST,'Weights',w(TrInd),'Prior','empirical');
                     % CLASSIFY                 
                     [~,PXCN] = predict(Classifier,fs(TestInd,:));
                     Find = find(x(TestInd)==-1*Tval);
                     Tind = find(x(TestInd)==Tval);  
                     switch Tval
                          case -1
                              PXCN = PXCN(:,2);
                          case 1
                              PXCN = PXCN(:,1);
                     end
                     [~,~,auc(i),X,Y,~,~,pauc(i)] = CATMV3TEST.getEER(PXCN(Tind),PXCN(Find));% CHECK ORDER
                     if display
                         plot(X,Y,'b-','LineWidth',1.5);
                     end
                 end
        end
        function [EER,G,AUC,x,y,TH,Y,pAUC] = getEER(tmatches,fmatches)
                g = [ones(1,length(tmatches)), -1*ones(1,length(fmatches))];
                [sorted,order] = sort([tmatches;fmatches],'ascend');
                g = g(order);
                true_neg = g == -1;nn = sum(true_neg);fpf = scale(cumsum(true_neg),nn);dx = diff(fpf);
                true_pos = g == 1;na = sum(true_pos);tpf = scale(cumsum(true_pos),na);dy = diff(tpf);
                y = tpf(1:end-1)+dy./2;
                x = fpf(1:end-1)+dx./2;
                [Y,indy] = max(y-x);
                yn = 1-y;d = abs(x-yn);
                [~,ind] = min(d);
                EER = ((x(ind)+yn(ind))/2);
                AUC = sum(dx.*y);
                FN = tpf(ind+1)*na;TN = fpf(ind+1)*nn;
                TP = na-FN;FP = nn-TN;
                G = 1-sqrt((TP/(TP+FN))*(TN/(TN+FP)));
                T(1,1) = TP;
                T(1,2) = FP;
                T(2,1) = FN;
                T(2,2) = TN;
                TH(1) = 1-sorted(ind);
                TH(2) = 1-sorted(indy);
                se = PPMMV2.standardError(AUC,na,nn);
                pAUC = normpdf((AUC-0.5)/se,0,1);
        end
        function [PREC,REC,G,TNF,XH,YH,ACC] = getHardClass(Tclass,Fclass)
                Pc = length(Tclass);Nc = length(Fclass);
                TP = length(find(Tclass==1));
                TN = length(find(Fclass==0));TNF = TN/Nc;
                FN = Pc-TP;FP = Nc-TN;
                PREC = TP/(TP+FP);
                REC = TP/(TP+FN);
                G = sqrt((TP/(TP+FN))*(TN/(TN+FP)));
                YH = TP/Pc;
                XH = 1-TN/Nc;
                ACC = (TP+TN)/(Pc+Nc);
        end
        function out = statNormPXC(PXC,OffsetPXC,ScalePXC)
                n = size(PXC,1);
                out =  (PXC-repmat(OffsetPXC,n,1))./repmat(ScalePXC,n,1);
                out(out<0) = 0;
                out(out>1) = 1; 
        end
    end
end



% CVO = cvpartition(species,'KFold',10);
%       err = zeros(CVO.NumTestSets,1);
%       for i = 1:CVO.NumTestSets
%           trIdx = CVO.training(i);
%           teIdx = CVO.test(i);
%           ytest = classify(meas(teIdx,:),meas(trIdx,:),species(trIdx,:));
%           err(i) = sum(~strcmp(ytest,species(teIdx)));
%       end
%       cvErr = sum(err)/sum(CVO.TestSize);
    
%     function out = crossvalidate(obj,X,FS,K)
%            if nargin<4, K = 10;end
%            n = length(X);
%            F = crossvalind('Kfold',n,K);
%            numfs = 10:10:obj.maxFS;
%            nfs = length(numfs);
%            EER = nan*zeros(K,nfs);
%            HG = nan*zeros(K,nfs);
%            parfor k=1:K
%                Tind = find(F==k);
%                Trind = setdiff(1:n,Tind);
%                forobj = clone(obj);
%                forX = X(Trind); %#ok<PFBNS>
%                forFS = FS(Trind,:); %#ok<PFBNS>
%                forTX = X(Tind);
%                forTFS = FS(Tind,:);
%                % PREPARATION
%                forBKX = forX;
%                [forX,COST,W,ind,~,~] = prepPredictorTraining(forobj,forX);
%                % FEATURE SELECTION
%                featureSelection(forobj,forX(ind),forBKX(ind),forFS(ind,:));
%                fsind = forobj.FSInd;
%                tmpEER = nan*zeros(1,nfs);
%                tmpHG = nan*zeros(1,nfs);
%                for i=1:1:nfs
%                    % PREP FS INDEX
%                      forobj.FSInd = fsind(1:numfs(i)); %#ok<PFBNS>
%                    % ACCUMULATE FEATURES
%                      sFS = selectFeatures(forobj,forFS); 
%                    % TRAINING
%                      fitClassifier(forobj,sFS(ind,:),forX(ind),COST,W(ind));
%                    % LEARNING HOW TO NORMALIZE  
%                      [~,PXC] = predict(forobj,forFS);
%                      learnPXCNorm(forobj,PXC);
%                    % TESTING  
%                      test = testPredictor(forobj,forTX,forTFS,false);
%                      tmpEER(i) = test.EER;
%                      tmpHG(i)  = test.HG;
%                end
%                EER(k,:) = tmpEER;
%                HG(k,:) = tmpHG;
%            end
%            out.EER = EER;
%            out.HG = HG;
%            [minEER,ind] = min(mean(EER,1));
%            out.nFSEER = numfs(ind);
%            out.minEER = minEER;
%            out.maxEER = max(mean(EER,1));
%            out.minEERHG = mean(HG(:,ind));
%            tmp = EER(:,ind);
%            out.FailedEER = sum(tmp>=0.5);      
%            [maxHG,ind] = max(mean(HG,1));
%            out.nFSHG = numfs(ind);
%            out.maxHG = maxHG;
%            out.minHG = min(mean(HG,1));
%            out.maxHGEER = mean(EER(:,ind));
%            tmp = HG(:,ind);
%            out.FailedHG = sum(tmp<=0.5);
%         end


% function fitClassifier(obj,FS,X,COST,W)
%             rng(1);
%             switch obj.ClassifierType
%                 case 'SVM'
%                     obj.Classifier = fitcsvm(FS,X,'Standardize',true,'Cost',COST,'Weights',W);
%                     %obj.Classifier = fitPosterior(obj.Classifier);
%                 case 'SVMRBF'
%                     switch obj.OptSVMRBF
%                         case false  
%                             %COST.ClassificationCosts(1,2) = 1;
%                             %COST.ClassificationCosts(2,1) = 1;
%                             tmp = fitcsvm(FS,X,'Standardize',false,'KernelFunction','rbf','BoxConstraint',1,...
%                                                      'KernelScale','auto','Cost',COST,'Weights',W,'Prior','empirical');
%                             obj.Classifier = fitcsvm(FS,X,'Standardize',false,'KernelFunction','rbf','BoxConstraint',1,...
%                                                      'KernelScale',tmp.KernelParameters.Scale,'Cost',COST,'Weights',W,'Prior','empirical');
%                            cv = crossval(obj.Classifier,'Kfold',5);
%                            obj.CVloss = kfoldLoss(cv);
%                         case true
%                             COST.ClassificationCosts(1,2) = 1;
%                             COST.ClassificationCosts(2,1) = 1;
%                             %COST.ClassificationCosts
%                             %[obj.CVloss,sigma,C,nFS] = CATMV3TEST.optSVMRBF(FS,X,COST,W,5);
%                             %obj.FSInd = obj.FSInd(1:nFS);
%                             %obj.Classifier = fitcsvm(FS(:,1:nFS),X,'KernelFunction','rbf',...
%                             %                         'KernelScale',sigma,'BoxConstraint',C,...
%                             %                         'Cost',COST,'Weights',W);
%                             [obj.CVloss,sigma,C] = CATMV3TEST.optSVMRBF1(FS,X,COST,W,5);
%                             obj.Classifier = fitcsvm(FS,X,'KernelFunction','rbf','Standardize',false,...
%                                                      'KernelScale',sigma,'BoxConstraint',C,...
%                                                      'Cost',COST,'Weights',W);
%                             
%                     end
%                     %cv = crossval(obj.Classifier,'Kfold',5);
%                     %obj.CVloss = kfoldLoss(cv);
%                 case 'ROBUSTBOOST'
%                     obj.Classifier = fitensemble(FS,X,'RobustBoost',obj.nrTrees,obj.treeTempl,'Cost',COST,'Weights',W);
%                 case 'RUSBOOST'
%                     obj.Classifier = fitensemble(FS,X,'RUSBoost',obj.nrTrees,obj.treeTempl,'Cost',COST,'Weights',W);
%                  otherwise
%                     error('WRONG TYPE OF CLASSIFIER');
%             end
%         end