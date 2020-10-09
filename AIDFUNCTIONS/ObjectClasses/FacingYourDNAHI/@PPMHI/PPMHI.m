classdef PPMHI <superClassLight
% PROPERTIES
    properties % GENERAL INTERFACING
    end
    properties (Dependent = true)
        nLev;
        BAL;
        FREQ;
        mFREQ;
        MFREQ;
        CALL;
    end
    properties (Hidden = true)
       HI;
       XX;
    end
    properties (Hidden = true, Dependent = true)
       nModules;
       nS;
       nSNotNan;
       XG;
       GG;
    end
    properties % FORWARD ASSOCIATION TESTING
       CCA_F;
       CCA_pF;
       CCA_CHI;
       CCA_pCHI;
       CCA_Wilks;
       CCA_R;
       CCA_B;
       CCA_STATS;
       CorrectY = true;
    end
    properties (Dependent = true)
    end
    properties (Hidden = true)
    end
    properties (Hidden = true, Dependent = true)
    end
    properties % FORWARD REPLICATION TESTING
    end
    properties (Dependent = true)
    end
    properties (Hidden = true)
    end
    properties(Hidden = true, Dependent = true)  
    end
    properties % CLASSIFIER
    end
    properties (Dependent = true)
        MaxScoreMod;
    end
    properties (Hidden = true)
        ModuleIndex = [];
        CLASS = {};
        gg;
        ngg;
        AUCS
        pAUCS;
        score;
        SVMOptions = [];
        SVMClassifier = [];
        FITCOptions = [];
        FITCClassifier = [];
        FITCSTATS = [];
        
        FITCCASOptions = [];
        FITCCAS = [];
        FITCCASSTATS = [];
        FITCCASFUSER = [];
        
        CONVNET = [];
        
%         SVMOptimize = false;
%         SVMStandardize = true;
%         SVMModules = [];
%         SVMAdasyn = false;
%         SVMFeatureSelection = false;
%         SVMFeatureIndex = [];
    end
    properties (Hidden = true, Dependent = true)
        nCLASS;
    end
    properties (Abstract = true)
        NanVal;
    end
% OBJECT METHODS    
    methods % CONSTRUCTOR
       function obj = PPMHI(varargin)
            obj = obj@superClassLight(varargin{:});
            obj.HI = HierarchicalInterface;
        end
    end
    methods % GENERAL GETTING
       function out = get.HI(obj)
           out = obj.HI;
           if~superClassLight.isH(out), out = []; end
       end
       function out = get.nLev(obj)
           if isempty(obj.HI), out = []; return; end
           out = obj.HI.nL;
       end
       function out = get.nModules(obj)
           if isempty(obj.HI), out = []; return; end
           out = obj.HI.nLC;
       end
       function out = get.nS(obj)
           out = length(obj.XX); 
       end
       function out = get.nSNotNan(obj)
           out = length(obj.XX)-sum(obj.XX==obj.NanVal); 
       end
       function out = get.BAL(obj)
           if isempty(obj.XX), out = []; return; end
           out = sum(obj.GG==obj.MinorityClass)/sum(obj.GG==-1*obj.MinorityClass);
       end
       function out = get.FREQ(obj)
           Tind = find(obj.GG==obj.MinorityClass);
           Find = find(obj.GG==-1*obj.MinorityClass);
           tot = length(Tind) + length(Find);
           out = [length(Tind)./tot length(Find)./tot];
       end
       function out = get.mFREQ(obj)
           out = min(obj.FREQ);
       end
       function out = get.MFREQ(obj)
           out = max(obj.FREQ);
       end
       function out = get.CALL(obj)
           if isempty(obj.XX), out = []; return; end
           out = (obj.nS-sum(obj.XX==obj.NanVal))./obj.nS;
       end
       function out = get.GG(obj)
                 out = XG2GG(obj,obj.XG); 
        end
       function out = get.XG(obj)
                 out = XX2XG(obj,obj.XX);
        end
       function out = get.nCLASS(obj)
                 out = length(obj.ModuleIndex);
        end
       function out = get.MaxScoreMod(obj)
                if isempty(obj.score), out = []; return; end
                [~,ind] = max(obj.score);
                out = ind(1);
       end
    end
    methods % GENERAL SETTING
       function obj = set.nLev(obj,in)
           if isempty(obj.HI), obj.HI=HierarchicalInterface;end
           obj.HI.L = single(in);
       end
    end
    methods % FORWARD ASSOCIATION ANALYSIS
        function out = CCA(obj,X,Y,C,varargin)
                 % reading options
                   options = PPMHI.readVarargin(varargin{:});
                 % Initializing
                   F = nan*zeros(1,obj.nModules);pF = nan*zeros(1,obj.nModules);
                   CHI = nan*zeros(1,obj.nModules);pCHI = nan*zeros(1,obj.nModules);
                   R = nan*zeros(1,obj.nModules);Wilks = nan*zeros(1,obj.nModules);
                   STATS = cell(1,obj.nModules); 
                   B = cell(1,obj.nModules);
                 % Converting X to G  
                   G = XX2GG(obj,X);
                   for i=1:1:obj.nModules
                       [F(i),pF(i),B{i},R(i),CHI(i),pCHI(i),Wilks(i),STATS{i}] =  PPMHI.parCCA(G,getCellData(obj.HI,Y,i),C,options.correcty);
                   end
                   out.F = F;out.pF = pF;out.CHI = CHI;out.pCHI = pCHI;out.Wilks = Wilks;out.R = R;out.B = B;out.STATS = STATS;
                   obj.CCA_F = F;obj.CCA_pF = pF;obj.CCA_CHI = CHI;obj.CCA_pCHI = pCHI;obj.CCA_Wilks = Wilks;obj.CCA_R = R;
                   obj.CCA_B = B;obj.CCA_STATS = STATS;
                   if options.display
                      plotHierValues3Dv2(-log10(out.pF),'crit',-log10(5*10^-8),'range',[0 max([-log10(5*10^-8) min(-log10(out.pF) -log10(1*10^-15))])],'title',num2str(min(out.pF))); 
                   end
        end
    end
    methods % FORWARD REPLICATION ANALYSIS
    end
    methods % BACKWARDS PREDICTION ANALYSIS
        function out = PLSR(obj,X,Y,C,varargin)
            % reading options
                   options = PPMHI.readVarargin(varargin{:});
                 % Initializing
                   COR = nan*zeros(1,obj.nModules);
                   pCOR = nan*zeros(1,obj.nModules);
                   PCTVAR = nan*zeros(1,obj.nModules);
                   F = nan*zeros(1,obj.nModules);
                   pF = nan*zeros(1,obj.nModules);
                 % Converting X to G
                   G = XX2GG(obj,X);
                   parfor i=1:1:obj.nModules
                       YFOR = [];
                       switch options.Hierarchical %#ok<*PFBNS>
                           case true
                               modindex = [i getAllChildren(obj.HI,i)];
                               YFOR = getCellData(obj.HI,Y,modindex);
                           case false
                               YFOR = getCellData(obj.HI,Y,i);
                       end
                       [COR(i),pCOR(i),PCTVAR(i),F(i),pF(i)] = PPMHI.predPLSR(G,YFOR,C,options);
                   end
                   out.COR = COR;out.pCOR = pCOR;out.PCTVAR = PCTVAR;out.F = F;out.pF = pF;
                   if options.display
                      figure;
                      ax = subplot(1,2,1);hold on;
                      plotHierValues3Dv2(-log10(out.pCOR),'crit',-log10(5*10^-8),'range',[0 max([-log10(5*10^-8) min(-log10(out.pCOR) -log10(1*10^-15))])],'title',num2str(min(out.pCOR)),'peer',ax); 
                      ax = subplot(1,2,2);hold on;
                      plotHierValues3Dv2(out.COR,'crit',max(out.COR)-0.01,'range',[-1 1],'title',num2str(max(out.COR)),'peer',ax);
                      
                      figure;
                      ax = subplot(1,2,1);hold on;
                      %plotHierValues3Dv2(-log10(out.pF),'crit',-log10(5*10^-8),'range',[0 max([-log10(5*10^-8) min(-log10(out.pF) -log10(1*10^-15))])],'title',num2str(min(out.pF)),'peer',ax); 
                      plotHierValues3Dv2(-log10(out.pF),'crit',-log10(5*10^-3),'range',[0 max([-log10(5*10^-3) min(-log10(out.pF) -log10(1*10^-5))])],'title',num2str(min(out.pF)),'peer',ax); 
                      ax = subplot(1,2,2);hold on;
                      plotHierValues3Dv2(out.F,'crit',max(out.F)-0.01,'range',[0 3],'title',num2str(max(out.F)),'peer',ax);
                   end
        end
        function out = FITC(obj,X,Y,C,varargin)
            % reading options
                   options = PPMHI.readVarargin(varargin{:});
                 % Initializing
                   CHI2 = nan*zeros(1,obj.nModules);
                   pCHI2 = nan*zeros(1,obj.nModules);
                   F = nan*zeros(1,obj.nModules);
                   pF = nan*zeros(1,obj.nModules);
                 % Converting X to G
                   G = XX2GG(obj,X);
                   if options.Weighted
                      W = GG2WW(obj,G);
                   else
                      W = ones(size(G));
                   end
                   parfor i=1:1:obj.nModules
                       YFOR = [];
                       switch options.Hierarchical %#ok<*PFBNS>
                           case true
                               modindex = [i getAllChildren(obj.HI,i)];
                               YFOR = getCellData(obj.HI,Y,modindex);
                           case false
                               YFOR = getCellData(obj.HI,Y,i);
                       end
                       [CHI2(i),pCHI2(i),F(i),pF(i)] = PPMHI.predFITC(G,YFOR,C,options,W);
                   end
                   out.CHI2 = CHI2;out.pCHI2 = pCHI2;out.F = F;out.pF = pF;
                   if options.display
                      figure;
                      ax = subplot(1,2,1);hold on;
                      plotHierValues3Dv2(-log10(out.pCHI2),'crit',-log10(5*10^-8),'title',num2str(min(out.pCHI2)),'peer',ax); 
                      ax = subplot(1,2,2);hold on;
                      plotHierValues3Dv2(out.CHI2,'title',num2str(max(out.CHI2)),'peer',ax);
                      
                      figure;
                      ax = subplot(1,2,1);hold on;
                      %plotHierValues3Dv2(-log10(out.pF),'crit',-log10(5*10^-8),'range',[0 max([-log10(5*10^-8) min(-log10(out.pF) -log10(1*10^-15))])],'title',num2str(min(out.pF)),'peer',ax); 
                      plotHierValues3Dv2(-log10(out.pF),'crit',-log10(5*10^-8),'title',num2str(min(out.pF)),'peer',ax); 
                      ax = subplot(1,2,2);hold on;
                      plotHierValues3Dv2(out.F,'crit',max(out.F)-0.01,'title',num2str(max(out.F)),'peer',ax);
                   end
        end
    end
    methods % CLASSIFIERS
        function trainFITC(obj,X,Y,index,varargin)
                 options = PPMHI.readVarargin(varargin{:});
                 obj.FITCOptions = PPMHI.readFITCVarargin(varargin{:});
                 obj.FITCOptions.Modules = index;
                 G = XX2GG(obj,X);
                 useindex = PPMHI.notNAN(G);
                 G = G(useindex);
                 Y = reduceCellData(obj.HI,Y,useindex);
                 if options.Weighted
                    W = GG2WW(obj,G);
                 else
                    W = ones(size(G));
                 end
                 FS = getCellData(obj.HI,Y,obj.FITCOptions.Modules); 
                 if obj.FITCOptions.Standardize,[FS,obj.FITCOptions.AvgFS,obj.FITCOptions.StdFS] = PPMHI.standardize(FS);end        
                 [STATS.CHI2,STATS.pCHI2,STATS.F,STATS.pF,STATS.f,STATS.pf,STATS.LAB,STATS.posterior,obj.FITCClassifier] = PPMHI.predFITC(G,FS,[],options,W); 
                 tmpG = obj.GG;uind = PPMHI.notNAN(tmpG);
                 gg = unique(tmpG(uind));ngg = length(gg);
                 AUC = zeros(1,ngg);
                 pAUC = zeros(1,ngg);  
                 if options.display
                    figure;hold on;grid on;plot(0:0.1:1,0:0.1:1,'k-');
                    plot(0:0.1:1,1:-0.1:0,'k--');
                 end 
                 for i=1:ngg
                     %i=1;
                     Find = find(~(G==gg(i)));
                     Tind = find(G==gg(i));
                     PXCN = -1*STATS.posterior(:,i);
                     [~,~,AUC(i),X,Y,~,~,pAUC(i)] = PPMHI.getEER(PXCN(Tind,1),PXCN(Find,1));
                     if options.display,plot(X,Y,'LineWidth',1.5);end
                 end
                 STATS.AUC = AUC;
                 STATS.pAUC = pAUC;
                 obj.FITCSTATS = STATS;   
        end
        function STATS = testFITC(obj,X,Y,varargin)
                 options = PPMHI.readVarargin(varargin{:});
                 G = XX2GG(obj,X);
                 useindex = PPMHI.notNAN(G);
                 G = G(useindex);
                 Y = reduceCellData(obj.HI,Y,useindex);
                 FS = getCellData(obj.HI,Y,obj.FITCOptions.Modules);

                 [STATS.LAB,STATS.posterior] = predict(obj.FITCClassifier,FS);
                 [STATS.CHI2,STATS.pCHI2,STATS.F,STATS.pF,STATS.f,STATS.pf,STATS.TABLE] = PPMHI.getCHIF(G,STATS.LAB,STATS.posterior);
                 
                 gg = unique(G);ngg = length(gg);
                 AUC = zeros(1,ngg);
                 pAUC = zeros(1,ngg);  
                 if options.display
                    figure;hold on;grid on;plot(0:0.1:1,0:0.1:1,'k-');
                    plot(0:0.1:1,1:-0.1:0,'k--');
                 end 
                 for i=1:ngg
                     %i=1;
                     Find = find(~(G==gg(i)));
                     Tind = find(G==gg(i));
                     PXCN = -1*STATS.posterior(:,i);
                     [~,~,AUC(i),X,Y,~,~,pAUC(i)] = PPMHI.getEER(PXCN(Tind,1),PXCN(Find,1));
                     plot(X,Y,'LineWidth',1.5);
                 end
                 STATS.AUC = AUC;
                 STATS.pAUC = pAUC;
        end 
        function out = getFITCMatches(obj,X,Y,varargin)
                 %options = PPMHI.readVarargin(varargin{:});
                 G = XX2GG(obj,X);
                 n = length(G);
                 FS = getCellData(obj.HI,Y,obj.FITCOptions.Modules);
                 %FS = FS(:,obj.SVMOptions.FeatureIndex);
                 % CLASSIFY
                 [PXC,PXCN] = predict(obj.FITCClassifier,FS);
                 out = nan*zeros(1,n,n);
                 for i=1:1:n
                    index = find(obj.FITCClassifier.ClassNames==G(i));
                    if ~isempty(index), out(1,:,i) = PXCN(:,index); end
                 end
        end
        function trainCLASS(obj,X,Y,index,varargin)
                 if isempty(index), index = 1:obj.nModules;end
                 options = PPMHI.readVarargin(varargin{:});
                 obj.ModuleIndex = index;
                 CLASS = cell(1,obj.nCLASS);
                 G = XX2GG(obj,X);
                 useindex = PPMHI.notNAN(G);
                 G = G(useindex);
                 Y = reduceCellData(obj.HI,Y,useindex);
                 obj.gg = unique(G);obj.ngg = length(obj.gg);
                 if options.Weighted
                    W = GG2WW(obj,G);
                 else
                    W = ones(size(G));
                 end
                 options.correcty = false;
                 paucs = nan*zeros(obj.ngg,obj.nCLASS);
                 aucs = nan*zeros(obj.ngg,obj.nCLASS);
                 parfor i=1:obj.nCLASS
                    FS = getCellData(obj.HI,Y,index(i)); 
                    [CLASS{i}.STATS.CHI2,CLASS{i}.STATS.pCHI2,CLASS{i}.STATS.F,...
                     CLASS{i}.STATS.pF,CLASS{i}.STATS.f,CLASS{i}.STATS.pf,CLASS{i}.STATS.LAB,CLASS{i}.STATS.posterior,...
                     CLASS{i}.MDL] = PPMHI.predFITC(G,FS,[],options,W); %#ok<*PROPLC>
                     AUC = zeros(1,obj.ngg);pAUC = zeros(1,obj.ngg);
                     xx = cell(1,obj.ngg);
                     yy = cell(1,obj.ngg);
                     for j=1:obj.ngg
                         %i=1;
                         Find = find(~(G==obj.gg(j)));
                         Tind = find(G==obj.gg(j));
                         PXCN = -1*CLASS{i}.STATS.posterior(:,j);
                         [~,~,AUC(j),xx{j},yy{j},~,~,pAUC(j)] = PPMHI.getEER(PXCN(Tind,1),PXCN(Find,1));
                     end
                     CLASS{i}.AUC = AUC;
                     CLASS{i}.pAUC = pAUC;
                     CLASS{i}.xx = xx;
                     CLASS{i}.yy = yy;
                     paucs(:,i) = pAUC(:);
                     aucs(:,i) = AUC(:);
                 end
                 obj.pAUCS = paucs;
                 obj.AUCS = aucs;
                 
                 score = nan*zeros(1,obj.nCLASS);
                 for i=1:1:obj.nCLASS
                     score(i) = -log10(pfast(obj.pAUCS(:,i)));
                 end
                 obj.score = score;      
                 if options.display
                    figure;hold on;grid on;plot(0:0.1:1,0:0.1:1,'k-');
                    plot(0:0.1:1,1:-0.1:0,'k--');
                    cmap = colormap('lines');
                    str = {'-' ':' '-.'};
                    for i=1:obj.nCLASS
                        for j=1:obj.ngg
                            plot(CLASS{i}.xx{j},CLASS{i}.yy{j},'LineWidth',1.5,'LineStyle',str{j},'Color',cmap(i,:));
                        end
                    end
                 end
                 obj.CLASS = CLASS;
        end
        function out = testCLASS(obj,X,Y,varargin)
                 options = PPMHI.readVarargin(varargin{:});
                 G = XX2GG(obj,X);
                 useindex = PPMHI.notNAN(G);
                 G = G(useindex);
                 Y = reduceCellData(obj.HI,Y,useindex);
                 %posterior = zeros(length(G),obj.ngg,obj.nCLASS);
                 TEST = cell(1,obj.nCLASS);
                 for i=1:obj.nCLASS
                     FS = getCellData(obj.HI,Y,obj.ModuleIndex(i)); 
                     [~,posterior] = predict(obj.CLASS{i}.MDL,FS);
                     AUC = zeros(1,obj.ngg);pAUC = zeros(1,obj.ngg);
                     xx = cell(1,obj.ngg);
                     yy = cell(1,obj.ngg);
                     for j=1:obj.ngg
                         %i=1;
                         Find = find(~(G==obj.gg(j)));
                         Tind = find(G==obj.gg(j));
                         PXCN = -1*posterior(:,j);
                         [~,~,AUC(j),xx{j},yy{j},~,~,pAUC(j)] = PPMHI.getEER(PXCN(Tind,1),PXCN(Find,1));
                     end
                     TEST{i}.AUC = AUC;
                     TEST{i}.pAUC = pAUC;
                     TEST{i}.xx = xx;
                     TEST{i}.yy = yy;
                     TEST{i}.Score = pfast(pAUC);
                 end
                 out.TEST = TEST;
                 score = zeros(1,obj.nCLASS);
                 for i=1:obj.nCLASS
                    score(i) = -log10(out.TEST{i}.Score);
                 end
                 out.score = score;
                 [out.MaxScore,ind] = max(score);
                 if options.display
                    figure;hold on;grid on;plot(0:0.1:1,0:0.1:1,'k-');
                    plot(0:0.1:1,1:-0.1:0,'k--');
                    cmap = colormap('lines');
                    str = {'-' ':' '-.'};
                    for i=1:obj.nCLASS
                        for j=1:obj.ngg
                            if i==ind(1)
                                plot(TEST{i}.xx{j},TEST{i}.yy{j},'LineWidth',2,'LineStyle',str{j},'Color',cmap(i,:));
                            else
                                plot(TEST{i}.xx{j},TEST{i}.yy{j},'LineWidth',1,'LineStyle',str{j},'Color',cmap(i,:));
                            end
                        end
                    end
                 end       
        end
        function out = getCLASSMatches(obj,X,Y,what,varargin)
                 %options = PPMHI.readVarargin(varargin{:});
                 G = XX2GG(obj,X);
                 n = length(G);
                 switch lower(what)
                     case 'all'
                         out = nan*zeros(obj.nCLASS,n,n);
                         parfor i=1:1:obj.nCLASS
                             FS = getCellData(obj.HI,Y,obj.ModuleIndex(i)); 
                             [~,posterior] = predict(obj.CLASS{i}.MDL,FS);
                             forout = nan*zeros(n,n);
                             for j=1:1:n
                                index = find(obj.CLASS{i}.MDL.ClassNames==G(j));
                                if ~isempty(index), forout(:,j) = posterior(:,index); end
                             end
                             out(i,:,:) = forout;
                         end
                     case 'best'
                         [~,ind] = max(obj.score);
                         %disp(num2str(ind));
                         out = nan*zeros(1,n,n);
                         FS = getCellData(obj.HI,Y,obj.ModuleIndex(ind)); 
                         [~,posterior] = predict(obj.CLASS{ind}.MDL,FS);
                         forout = nan*zeros(n,n);
                         for j=1:1:n
                             index = find(obj.CLASS{ind}.MDL.ClassNames==G(j));
                             if ~isempty(index), forout(:,j) = posterior(:,index); end
                         end
                         out(1,:,:) = forout;                  
                 end
        end
        
        function trainDL(obj,X,Y,varargin)
                 options = PPMHI.readVarargin(varargin{:});
                 G = XX2GG(obj,X);
                 useindex = PPMHI.notNAN(G);
                 G = G(useindex);
                 Y = Y(:,:,:,useindex);
                 if options.Weighted
                    W = GG2WW(obj,G);
                 else
                    W = ones(size(G));
                 end
                 RES = size(Y,1);
                 gg = unique(obj.GG);ngg = length(gg);
                 
                 layers = [imageInputLayer([RES RES 3]);
                 convolution2dLayer(15,20);
                 reluLayer();
                 maxPooling2dLayer(2,'Stride',2);
                 dropoutLayer(0.5);
                 fullyConnectedLayer(10);
                 reluLayer();
                 fullyConnectedLayer(ngg);
                 softmaxLayer();
                 classificationLayer()];
                 
                 dloptions = trainingOptions('sgdm','MaxEpochs',30,...
	                                       'InitialLearnRate',0.001,'MiniBatchSize',30,'ExecutionEnvironment','gpu');
            
                 obj.CONVNET = trainNetwork(Y,categorical(G),layers,dloptions);
                 
                 [XTest,ERR] = classify(obj.CONVNET,Y);
                                  
                 AUC = zeros(1,ngg);
                 pAUC = zeros(1,ngg);  
                 if options.display
                    figure;hold on;grid on;plot(0:0.1:1,0:0.1:1,'k-');
                    plot(0:0.1:1,1:-0.1:0,'k--');
                 end 
                 for i=1:ngg
                     %i=1;
                     Find = find(~(G==gg(i)));
                     Tind = find(G==gg(i));
                     PXCN = -1*ERR(:,i);
                     [~,~,AUC(i),X,Y,~,~,pAUC(i)] = PPMHI.getEER(PXCN(Tind,1),PXCN(Find,1));
                     if options.display, plot(X,Y,'LineWidth',1.5);end
                 end
        end
        function out = testDL(obj,X,Y,varargin)
                 options = PPMHI.readVarargin(varargin{:});
                 G = XX2GG(obj,X);
                 useindex = PPMHI.notNAN(G);
                 G = G(useindex);
                 Y = Y(:,:,:,useindex);
                 
                 gg = unique(obj.GG);ngg = length(gg);
                 
                 [GTest,ERR] = classify(obj.CONVNET,Y);
                 
                 out.Accuracy = sum(GTest == categorical(G))/numel(G);   
                  
                 AUC = zeros(1,ngg);
                 pAUC = zeros(1,ngg);  
                 if options.display
                    figure;hold on;grid on;plot(0:0.1:1,0:0.1:1,'k-');
                    plot(0:0.1:1,1:-0.1:0,'k--');
                 end 
                 for i=1:ngg
                     %i=1;
                     Find = find(~(G==gg(i)));
                     Tind = find(G==gg(i));
                     PXCN = -1*ERR(:,i);
                     [~,~,AUC(i),X,Y,~,~,pAUC(i)] = PPMHI.getEER(PXCN(Tind,1),PXCN(Find,1));
                     if options.display, plot(X,Y,'LineWidth',1.5);end
                 end
                 out.AUC = AUC;
                 out.pAUC = pAUC;
        end
        function out = getDLMatches(obj,X,Y)
                 %options = PPMHI.readVarargin(varargin{:});
                 G = XX2GG(obj,X);
                 n = length(G);
                 out = nan*zeros(1,n,n);
                 [~,posterior] = classify(obj.CONVNET,Y);
                 forout = nan*zeros(n,n);
                 gg = obj.CONVNET.Layers(end).ClassNames;
                 ngg = length(gg);
                 classnames = zeros(1,ngg);
                 for i=1:1:ngg
                    classnames(i) = str2double(gg{i});
                 end
                 for j=1:1:n
                     index = find(classnames==G(j));
                     if ~isempty(index), forout(:,j) = posterior(:,index); end                  
                 end
                 out(1,:,:) = forout;
        end
    end
    methods % INTERFACIN
       function reset(obj)
           obj.pF = [];
           obj.F = [];
           obj.R = []; %#ok<*MCHV3>
       end
       function out = XX2GG(obj,in)
                out = XG2GG(obj,XX2XG(obj,in));
       end
       function out = GG2WW(obj,in)
                out = zeros(size(in));
                index = PPMHI.notNAN(in);
                gg = unique(in(index));
                ngg = length(gg);
                %nN = length(index);
                for n=1:1:ngg
                    ind = find(in==gg(n));
                    %out(ind) = (nN-length(ind))./nN;
                    out(ind) = 1/length(ind);
                end
                out = out/ngg;
       end
    end
    methods (Abstract = true)
        out = XX2XG(obj,in);
        out = XG2GG(obj,in);
    end
% STATIC METHODS    
    methods (Static = true) % ANALYSIS
        function [F,pF,B,R,CHI,pCHI,Wilks,STATS] = parCCA(X,Y,C,correctY)
               if nargin<4, correctY = true;end
               [X,Y] = PPMHI.prepareXYC(X,Y,C,correctY);
               [~,B,R,~,~,STATS] = canoncorr(double(X),double(Y));
               F = STATS.F(end);
               pF = STATS.pF(end);
               CHI = STATS.chisq(end);
               pCHI = STATS.pChisq(end);
               Wilks = STATS.Wilks(end);
        end
        function [out,M] = getResiduals(X,Y)
                [~,~,~,~,M] = plsregress(X,Y,min(size(X,2),size(Y,2)));
                Y_est = [ones(size(X,1),1) X]*M;
                out = Y-Y_est;
        end
        function [C,pC,PCTVAR,F,pF] = predPLSR(X,Y,C,options)
                 G = X;
                 [X,Y,~,index] = PPMHI.prepareXYC(X,Y,C,options.correcty);
                 G = G(index);
                 [~,~,~,~,M,PCTVAR] = plsregress(Y,X,1);
                 PCTVAR = PCTVAR(2);
                 if options.LXO==0
                     X_est = [ones(size(Y,1),1) Y]*M;
                 else
                    n = length(X);
                    nF = floor(n/options.LXO);
                    F = crossvalind('Kfold',n,nF);
                    X_est = nan*zeros(size(X));
                    for f=1:1:max(F)
                        TeInd = find(F==f);
                        TrInd = setdiff(1:n,TeInd);
                        [~,~,~,~,M] = plsregress(Y(TrInd,:),X(TrInd),1);
                        X_est(TeInd) = [ones(size(Y(TeInd,:),1),1) Y(TeInd,:)]*M;
                    end
                 end
                 [C,pC] = corr(X,X_est);
                 [pF,anovatab] = anova1(X_est,G,'off');
                 F = anovatab{2,5};
                 %figure;boxplot(X_est,G);
        end
        function [Mdl] = getFITC(Y,G,W,type,NN)
                 if nargin<5, NN = 10; end
                 switch type
                         case 'LDA'
                             Mdl = fitcdiscr(Y,G,'Weights',W,'DiscrimType','linear');
                         case 'QDA'
                             Mdl = fitcdiscr(Y,G,'Weights',W,'DiscrimType','pseudoquadratic');
                         case 'NB'
                             Mdl = fitcnb(Y,G,'Weights',W);
                         case 'KNN'
                             Mdl = fitcknn(Y,G,'Weights',W,'NumNeighbors',NN);
                         case 'SVM'
                             Mdl = fitcsvm(Y,G,'Weights',W,'KernelFunction','rbf','KernelScale','auto');
                 end
        end
        function [CHI2,pCHI2,F,pF,f,pf,LAB,posterior,MDL] = predFITC(X,Y,C,options,W)
                 %if nargin<6, W = ones(size(X)); end
                 G = X;
                 [X,Y,~,index] = PPMHI.prepareXYC(X,Y,C,options.correcty);
                 G = G(index);
                 gg = unique(G);ngg=length(gg);
                 MDL = PPMHI.getFITC(Y,G,W(index),options.FitcType,options.NumNeighbors);
                 if options.LXO==0
                     [LAB,posterior] = predict(MDL,Y);
                 else
                    n = length(X);
                    nF = floor(n/options.LXO);
                    rng(10);F = crossvalind('Kfold',n,nF);
                    LAB = nan*zeros(size(X));
                    posterior = nan*zeros(length(index),ngg);
                    for f=1:1:max(F)
                        %disp(num2str(f));
                        TeInd = find(F==f);
                        TrInd = setdiff(1:n,TeInd);
                        Mdl = PPMHI.getFITC(Y(TrInd,:),G(TrInd),W(TrInd),options.FitcType,options.NumNeighbors);
                        [LAB(TeInd),posterior(TeInd,:)] = predict(Mdl,Y(TeInd,:));
                    end
                 end
                 [CHI2,pCHI2,F,pF,f,pf] = PPMHI.getCHIF(G,LAB,posterior);
        end
        function [CHI2,pCHI2,F,pF,f,pf,TABLE] = getCHIF(G,LAB,posterior)
                 gg = unique(G);ngg=length(gg);
                 [TABLE,CHI2,pCHI2] = crosstab(G,LAB);
                 f = nan*zeros(1,ngg);
                 pf = nan*zeros(1,ngg);
                 for n=1:1:ngg
                     [pf(n),anovatab] = anova1(posterior(:,n),G,'off');
                     f(n) = anovatab{2,5};
                 end
                 [pF,ind] = min(pf);
                 F = f(ind);
        end
        function [X,Y,C,index] = prepareXYC(X,Y,C,correctY)
            if ~isempty(C)
                  index = intersect(PPMHI.notNAN(C),PPMHI.notNAN(X));
                  if correctY
                     Y = PPMMV3.getResiduals(C(index,:),Y(index,:));
                  else
                     Y = Y(index,:);
                  end
                  X = PPMMV3.getResiduals(C(index,:),X(index,:));
               else
                 index = PPMHI.notNAN(X);   
                 Y = Y(index,:);
                 X = X(index,:);
            end
        end
    end
    methods (Static = true) % CLASSIFIER
        function [out,inavg,instd] = standardize(in)
                [nO,~] = size(in); 
                inavg = mean(in,1);
                instd = std(in,[],1);
                out = in-repmat(inavg,nO,1);
                out = out./repmat(instd,nO,1);
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
                if AUC>=0.5
                   pAUC = normpdf((AUC-0.5)/se,0,1);
                else
                   pAUC = 1;
                end
        end
    end
    methods (Static = true) % INTERFACING
       function out = notNAN(in)
           index = (1:size(in,1));
           [i,~] = find(isnan(in));
           out = setdiff(index,unique(i));
       end
       function out = getL1CL1DepVar(in)
           if ~iscell(in), out = in;return;end
           out = in{1}.DepVar{1};
       end 
       function out = convertUInt(in)
           if isempty(in), out = []; return; end
           M = max(in(:));
           if M<=intmax('uint8'),out = uint8(in);return;end
           if M<=intmax('uint16'), out = uint16(in);return;end
           if M<=intmax('uint32'), out = uint32(in);return;end
           out = uint64(in);
       end
       function out = CellList2Structure(in)
                out = [];
                if isempty(in), return; end
                nrCells = length(in);
                succes = 0;
                for i=1:1:nrCells
                   if isempty(in{i}),continue;end 
                   names = fieldnames(in{i});
                   for n=1:1:length(names)
                       eval(['out.' names{n} ' = [];']);
                   end
                   succes = 1;
                   break;
                end
                if ~succes, return;end
                for i=1:1:nrCells
                   if isempty(in{i})
                      for n=1:1:length(names)
                        eval(['out.' names{n} ' = [out.' names{n} ' nan];']);
                      end
                   else
                      for n=1:1:length(names)
                        eval(['out.' names{n} ' = [out.' names{n} ' in{i}.' names{n} '];']);
                      end
                   end
                end
       end
       function options = readVarargin(varargin)
                 Input = find(strcmpi(varargin, 'options'));
                 if ~isempty(Input), options = varargin{Input+1}; return;end    
                 Input = find(strcmpi(varargin, 'display'));
                 if isempty(Input),options.display = false;else, options.display = varargin{Input+1};end
                 Input = find(strcmpi(varargin, 'k'));
                 if isempty(Input),options.K = 3;else, options.K = varargin{Input+1};end
                 Input = find(strcmpi(varargin, 'LXO'));
                 if isempty(Input),options.LXO = 10;else, options.LXO = varargin{Input+1};end
                 Input = find(strcmpi(varargin, 'hierarchical'));
                 if isempty(Input),options.Hierarchical = false;else, options.Hierarchical = varargin{Input+1};end
                 Input = find(strcmpi(varargin, 'fitctype'));
                 if isempty(Input),options.FitcType = false;else, options.FitcType = varargin{Input+1};end
                 Input = find(strcmpi(varargin, 'correcty'));
                 if isempty(Input),options.correcty = false;else, options.correcty = varargin{Input+1};end
                 Input = find(strcmpi(varargin, 'NumNeighbors'));
                 if isempty(Input),options.NumNeighbors = 10;else, options.NumNeighbors = varargin{Input+1};end
                 Input = find(strcmpi(varargin, 'Weighted'));
                 if isempty(Input),options.Weighted = false;else, options.Weighted = varargin{Input+1};end
       end
       function options = readSVMVarargin(varargin)
                 Input = find(strcmpi(varargin, 'Options'));
                 if ~isempty(Input), options = varargin{Input+1}; return;end
                 Input = find(strcmpi(varargin, 'K'));
                 if isempty(Input),options.K = 3;else, options.K = varargin{Input+1};end
                 Input = find(strcmpi(varargin, 'Optimize'));
                 if isempty(Input), options.Optimize = false; else, options.Optimize = varargin{Input+1};end
                 Input = find(strcmpi(varargin, 'Standardize'));
                 if isempty(Input), options.Standardize = true; else, options.Standardize = varargin{Input+1};end
                 Input = find(strcmpi(varargin, 'Modules'));
                 if isempty(Input), options.Modules = []; else, options.Modules = varargin{Input+1};end
                 Input = find(strcmpi(varargin, 'Adasyn'));
                 if isempty(Input), options.Adasyn = []; else, options.Adasyn = varargin{Input+1};end
                 Input = find(strcmpi(varargin, 'FeatureSelection'));
                 if isempty(Input), options.FeatureSelection = []; else, options.FeatureSelection = varargin{Input+1};end
                 Input = find(strcmpi(varargin, 'FeatureIndex'));
                 if isempty(Input), options.FeatureIndex = []; else, options.FeatureIndex = varargin{Input+1};end
                 Input = find(strcmpi(varargin, 'Prior'));
                 if isempty(Input), options.Prior = 'empirical'; else, options.Prior = varargin{Input+1};end
                 Input = find(strcmpi(varargin, 'C'));
                 if isempty(Input), options.C = 1; else, options.C = varargin{Input+1};end
       end
       function options = readFITCVarargin(varargin)
                 Input = find(strcmpi(varargin, 'Options'));
                 if ~isempty(Input), options = varargin{Input+1}; return;end
                 Input = find(strcmpi(varargin, 'K'));
                 if isempty(Input),options.K = 3;else, options.K = varargin{Input+1};end
                 Input = find(strcmpi(varargin, 'Optimize'));
                 if isempty(Input), options.Optimize = false; else, options.Optimize = varargin{Input+1};end
                 Input = find(strcmpi(varargin, 'Standardize'));
                 if isempty(Input), options.Standardize = false; else, options.Standardize = varargin{Input+1};end
                 Input = find(strcmpi(varargin, 'Modules'));
                 if isempty(Input), options.Modules = []; else, options.Modules = varargin{Input+1};end
                 Input = find(strcmpi(varargin, 'Adasyn'));
                 if isempty(Input), options.Adasyn = []; else, options.Adasyn = varargin{Input+1};end
                 Input = find(strcmpi(varargin, 'FeatureSelection'));
                 if isempty(Input), options.FeatureSelection = []; else, options.FeatureSelection = varargin{Input+1};end
                 Input = find(strcmpi(varargin, 'FeatureIndex'));
                 if isempty(Input), options.FeatureIndex = []; else, options.FeatureIndex = varargin{Input+1};end
                 Input = find(strcmpi(varargin, 'Prior'));
                 if isempty(Input), options.Prior = 'empirical'; else, options.Prior = varargin{Input+1};end
                 Input = find(strcmpi(varargin, 'C'));
                 if isempty(Input), options.C = 1; else, options.C = varargin{Input+1};end
                 Input = find(strcmpi(varargin, 'fitctype'));
                 if isempty(Input),options.FitcType = false;else, options.FitcType = varargin{Input+1};end
                 Input = find(strcmpi(varargin, 'Weighted'));
                 if isempty(Input),options.Weighted = false;else, options.Weighted = varargin{Input+1};end
                 Input = find(strcmpi(varargin, 'fuser'));
                 if isempty(Input),options.fuser = 'SUM';else, options.fuser = varargin{Input+1};end
       end
    end
    methods (Static = true) % IMAGING
        function v = setupViewer(scan)
            if nargin<1
               v = viewer3DObj;
            else
               v = viewer(scan);
            end
            v.BackgroundColor = [1 1 1];v.AxesVisible = false;v.AxesWallColor = [1 1 1];
            v.AxesXColor = [0 0 0];v.AxesYColor = [0 0 0];v.AxesZColor = [0 0 0];
            v.SceneLightVisible = true;
            v.SceneLightLinked = true;
        end
        
    end
end


% function out = getCLASSMatch(obj,X,Y,varargin)
%                  %options = PPMHI.readVarargin(varargin{:});
%                  G = XX2GG(obj,X);
%                  n = length(G);
%                  [M,ind] = max(obj.score);
%                  out = nan*zeros(1,n,n);
%                  FS = getCellData(obj.HI,Y,obj.ModuleIndex(ind)); 
%                  [~,posterior] = predict(obj.CLASS{ind}.MDL,FS);
%                  forout = nan*zeros(n,n);
%                  for j=1:1:n
%                      index = find(obj.CLASS{ind}.MDL.ClassNames==G(j));
%                      if ~isempty(index), forout(:,j) = posterior(:,index); end
%                  end
%                  out(1,:,:) = forout;
%         end
% 


% function trainSVM(obj,X,Y,varargin)
%                 % Reading options
%                  options = PPMHI.readVarargin(varargin{:});
%                  obj.SVMOptions = PPMHI.readSVMVarargin(varargin{:});
%                  
% %                  options = PPMHI.readVarargin('display',true);
% %                  obj.SVMOptions = PPMHI.readSVMVarargin();          
%                  
%                 % Preparing X for classification
%                  [G,COST,W,ind,Tind,Find] = prepTrainSVM(obj,X);
%                 % Setting up Cross Validation Partition    
%                  if (~isobject(obj.SVMOptions.K)&&~(obj.SVMOptions.K==0)),rng(1);obj.SVMOptions.K = cvpartition(length(ind),'KFold',obj.SVMOptions.K);end
%                 % Feature Accumulation and selection
%                  if isempty(obj.SVMOptions.Modules), obj.SVMOptions.Modules = 1:obj.nModules; end
%                  FS = getCellData(obj.HI,Y,obj.SVMOptions.Modules);
%                  if obj.SVMOptions.FeatureSelection
%                      % do something
%                      nMod = length(obj.SVMOptions.Modules);
%                  else
%                     obj.SVMOptions.FeatureIndex = 1:size(FS,2); 
%                  end
%                  FS = FS(ind,obj.SVMOptions.FeatureIndex);
%                  G = G(ind);W = W(ind);
%                 % Feature Standardization
%                  if obj.SVMOptions.Standardize,[FS,obj.SVMOptions.AvgFS,obj.SVMOptions.StdFS] = PPMHI.standardize(FS);end
%                 % SVM TRAINING
%                  if ~obj.SVMOptions.Optimize
%                      obj.SVMClassifier = fitcsvm(FS,G,'Standardize',false,'KernelFunction','rbf',...
%                                                  'BoxConstraint',obj.SVMOptions.C,'KernelScale','auto','Cost',COST,...
%                                                  'Weights',W,'Prior',obj.SVMOptions.Prior);
%                      obj.SVMOptions.Scale = obj.SVMClassifier.KernelParameters.Scale;
%                      if isobject(obj.SVMOptions.K),obj.SVMOptions.CV = cvSVM(obj,G,FS,obj.SVMOptions.K,COST,W,obj.SVMOptions.C,obj.SVMOptions.Scale,options.display);end 
%                  else
%                      % check out bayesopt and other optimization options!
%                  end
%         end
%         function out = cvSVM(obj,G,FS,K,COST,W,C,Scale,display)
%                  if ~isobject(K), K = cvpartition(length(G),'KFold',K); end
%                  auc = zeros(1,K.NumTestSets);
%                  pauc = zeros(1,K.NumTestSets);
%                  if display
%                    figure;hold on;grid on;plot(0:0.1:1,0:0.1:1,'k-');
%                    plot(0:0.1:1,1:-0.1:0,'k--');
%                  end
%                  for i=1:K.NumTestSets
%                      TrInd = K.training(i);
%                      TestInd = K.test(i);
%                      % BUILD CLASSIFIER
%                      Classifier = fitcsvm(FS(TrInd,:),G(TrInd),'Standardize',false,'KernelFunction','rbf','BoxConstraint',C,...
%                                           'KernelScale',Scale,'Cost',COST,'Weights',W(TrInd),'Prior',obj.SVMOptions.Prior);
%                      % CLASSIFY                 
%                      [~,PXCN] = predict(Classifier,FS(TestInd,:));
%                      Find = find(G(TestInd)==-1*obj.MinorityClass);
%                      Tind = find(G(TestInd)==obj.MinorityClass);  
%                      switch obj.MinorityClass
%                           case -1
%                               PXCN = PXCN(:,2);
%                           case 1
%                               PXCN = PXCN(:,1);
%                      end
%                      [~,~,auc(i),X,Y,~,~,pauc(i)] = PPMHI.getEER(PXCN(Tind,1),PXCN(Find,1));% CHECK ORDER
%                      if display
%                          plot(X,Y,'b-','LineWidth',1.5);
%                      end
%                  end
%                  out.auc = auc;
%                  out.pauc = pauc;
%         end
%         function out = testSVM(obj,X,Y,varargin)
%                  options = PPMHI.readVarargin(varargin{:});
%                  G = XX2GG(obj,X);
%                  FS = getCellData(obj.HI,Y,obj.SVMOptions.Modules);
%                  FS = FS(:,obj.SVMOptions.FeatureIndex);
%                  % CLASSIFY
%                  [~,PXCN] = predict(obj.SVMClassifier,FS);
%                  Find = find(G==-1*obj.MinorityClass);
%                  Tind = find(G==obj.MinorityClass);
%                  switch obj.MinorityClass
%                       case -1
%                          PXCN = PXCN(:,2);
%                       case 1
%                          PXCN = PXCN(:,1);
%                  end
%                  [~,~,out.auc,X,Y,~,~,out.pauc] = PPMHI.getEER(PXCN(Tind,1),PXCN(Find,1));%#ok<*FNDSB> % CHECK ORDER
%                  if options.display
%                     figure;hold on;grid on;plot(0:0.1:1,0:0.1:1,'k-');
%                     plot(0:0.1:1,1:-0.1:0,'k--'); 
%                     plot(X,Y,'b-','LineWidth',1.5);
%                  end
%         end
%         function out = getSVMMatches(obj,X,Y,varargin)
%                  %options = PPMHI.readVarargin(varargin{:});
%                  G = XX2GG(obj,X);
%                  n = length(G);
%                  FS = getCellData(obj.HI,Y,obj.SVMOptions.Modules);
%                  FS = FS(:,obj.SVMOptions.FeatureIndex);
%                  % CLASSIFY
%                  [PXC,PXCN] = predict(obj.SVMClassifier,FS);
%                  out = nan*zeros(n,n);
%                  for i=1:1:n
%                     switch G(i)
%                         case 1
%                             %out(:,i) = PXCN(:,1);
%                             out(:,i) = PXCN(:,2);
%                         case -1
%                             out(:,i) = PXCN(:,1);
%                             %out(:,i) = PXCN(:,2);
%                         otherwise
%                             continue;
%                     end
%                  end
%         end