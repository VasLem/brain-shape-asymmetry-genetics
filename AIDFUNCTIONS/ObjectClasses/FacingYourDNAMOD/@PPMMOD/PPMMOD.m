classdef PPMMOD <superClassLight
    properties
       X = [];
       Level = 1;
       Cluster = 1;
    end
    properties (Dependent = true)
       nS;
       AvgX;
       StdX;
       PercNan;
    end
    properties % ANALYSIS
       R2 = [];
       pR2 = [];
       A = [];
       pA = [];
       pAR = [];
       EER = [];
    end
    properties % MODELING
       HTest = true;
       HTestp = 0.0001;
       RIPNorm = true;
    end
    properties % MATCHING
       Match = 'QD';
       Classify = 'SOFT';
       ScoreWeights = 'P';
       Kappa = 6;
       TH = 0.5;
    end
    properties (Hidden = true)
       VarInfo;
       M;
       H;
       HM;
       MRIP;
       NormMRIP;
       HMRIP;
       NormHMRIP;
    end
    properties(Hidden = true, Dependent = true)  
       XMOD; % X used in model (excluding nan values in XREG)
       XMODInd; % Index of non nan XREG values
       MMOD;
       RIPMOD;
       DIM;
    end
    properties (Abstract = true)
       XREG; % X USED FOR REGRESSION/ SUBCLASS DEPENDENT
    end
    methods % CONSTRUCTOR
        function obj = PPMMOD(varargin)
            obj = obj@superClassLight(varargin{:});         
        end
    end
    methods % GENERAL SETTING/GETTING
        function out = get.nS(obj)
           out = length(obj.X); 
        end
        function out = get.XMODInd(obj)
           out = find(~isnan(obj.XREG)); 
        end
        function out = get.XMOD(obj)
           if isempty(obj.X), out = []; return; end
           out = obj.XREG(obj.XMODInd);
        end
        function out = get.RIPMOD(obj)
          if strcmp(obj.Type,'SNPMOD')% TEMPORARY HACK!!! TO BE REMOVED
              if obj.HTest, 
                 if obj.RIPNorm, out = obj.NormHMRIP; return; end
                 out = obj.HMRIP; return;
              end
              if obj.RIPNorm, out = obj.NormMRIP; return; end
              out = obj.MRIP;
          else
              if obj.HTest, 
                 if obj.RIPNorm, out = obj.NormHMRIP(obj.XMODInd); return; end
                 out = obj.HMRIP(obj.XMODInd); return;
              end
              if obj.RIPNorm, out = obj.NormMRIP(obj.XMODInd); return; end
              out = obj.MRIP(obj.XMODInd);
          end
        end
        function out = get.MMOD(obj)
           if obj.HTest, out = obj.HM; return; end
           out = obj.M;
        end
        function out = get.PercNan(obj)
           if isempty(obj.X), out = []; return; end
           out = round((sum(isnan(obj.X))/obj.nS)*100);
        end
        function out = get.DIM(obj)
           if isempty(obj.M), out = []; return; end
           out = length(obj.M);
        end
    end
    methods % TESTING
        function out = runR2Test(obj,DepVar,COV,t)
           [out.R2,out.pR2] = PPMMOD.testPartialPLSR(obj.XREG,getDepVar(obj,DepVar),COV,t);
           obj.R2 = out.R2;obj.pR2 = out.pR2;
        end
        function out = runAngleTest(obj,DepVar,COV,maxM,AngleTable,t)
           [out.A,out.M,out.STAT] = PPMMOD.AngleTest([COV obj.XREG],getDepVar(obj,DepVar),AngleTable,t,maxM,0);
           obj.A = out.STAT.AvgA;obj.pA = out.STAT.pA(1);obj.pAR = out.STAT.pAR;
        end
        function out = runBiometricTest(obj,TestDepVar,TestX,t)
            if ~isempty(TestDepVar)
                TestDepVar = getDepVar(obj,TestDepVar);
                TestRIP = getRIP(obj,TestDepVar);
                [TrX,TrRIP,TestX,TestRIP,posclass] = prepOppBiometricTest(obj,TestRIP,TestX);
            else
                [TrX,TrRIP,TestX,TestRIP,posclass] = prepOppBiometricTest(obj,[],[]);
            end   
            [out.EER,out.G,out.AUC,out.x,out.y,out.R,out.PREC,out.REC,out.GH,out.TNF,out.TH,out.Y,out.XH,out.YH] = oppBiometricTest(obj,TrX,TrRIP,TestX,TestRIP,posclass);
            %[out.EER,out.G,out.AUC,out.R,out.PREC,out.REC,out.GH,out.TNF] = oppBiometricTest(obj,TrX,TrRIP,TrX,TrRIP,posclass);
            if t==0, return; end
            N = length(TestX);
            EER = nan*ones(1,t);
            for i=1:1:t
                ind = randsample(N,N,true);
                EER(i) = oppBiometricTest(obj,TrX,TrRIP,TestX(ind),TestRIP(ind),posclass);
                if mod(i,100)==0
                   index = find(~isnan(EER)); 
                   tmppEER = (sum(EER(index)>=0.5)+1)/(length(index)+1); 
                   acc = 10/length(index);
                   if tmppEER>acc, break;end
                end
            end
            index = find(~isnan(EER));
            EER = EER(index); %#ok<*PROP>
            EER = sort(EER);
            out.AvgEER = nanmean(EER);
            out.UpperEER = EER(round(0.975*length(index)));
            out.LowerEER = EER(round(0.025*length(index)+1));
            out.pEER = (sum(EER>=0.5)+1)/(length(index)+1);
        end
    end
    methods % BUILDING
        function buildFromAngleTest(obj,DepVar,ATest)
           DepVar = getDepVar(obj,DepVar);
          % VARIABLE INFORMATION
           obj.VarInfo = getVarInfo(obj,DepVar);
          % M CONSTRUCTION 
           fullM = [squeeze(ATest.M(:,1,:)) squeeze(ATest.M(:,2,:))];
           obj.M = nanmedian(fullM,2)';
           obj.H = PPMMOD.wilcoxonM(fullM',obj.HTestp);
           obj.HM = obj.M.*obj.H;
          % RIP COMPUTATIONS
           obj.MRIP = PPMMOD.updateRIP(DepVar,obj.M)';
           obj.NormMRIP = PPMMOD.normalizeRIP(obj.MRIP,obj.M,obj.VarInfo);
           obj.HMRIP = PPMMOD.updateRIP(DepVar,obj.HM)';
           obj.NormHMRIP = PPMMOD.normalizeRIP(obj.HMRIP,obj.HM,obj.VarInfo);
           postBuild(obj);
        end
    end
    methods % BIOMETRICS
        function out = getMatch(obj,distr,rip)
          out = PPMMOD.getStaticMatch(distr,rip,obj.Match);  
        end
        function out = match2score(obj,matches)
          switch obj.Match
              case 'QD'
                  out = -log(matches);
              case 'INLIER'
                  out = 1-matches;
          end
        end
        function out = match2class(obj,matches)
          switch obj.Match
              case 'QD'
                  T = 1;
              case 'INLIER'
                  T = obj.TH;
          end
          out = PPMMOD.matchthresholding(matches,T);
        end
        function [QDTMatches,QDFMatches,INTMatches,INFMatches,wF,wA,wP,wPR] = biometricMatchesWeights(obj,TestDepVar,TestX,AngleTable)
           % INITIALIZE 
            TestDepVar = getDepVar(obj,TestDepVar);
            [nrT,dim] = size(TestDepVar);
           % ANGLE BASED WEIGHTS 
            wA = max(obj.A,0);
            wP = PPMMOD.pAWeight(obj.pA,0.05);
            wPR = PPMMOD.pRAWeight(obj.A,dim,AngleTable);
           % MATCHES and F WEIGHTS 
            RIP = getRIP(obj,TestDepVar);
            QDTMatches = nan*zeros(1,nrT);
            QDFMatches = nan*zeros(1,nrT,nrT-1);
            INTMatches = nan*zeros(1,nrT);
            INFMatches = nan*zeros(1,nrT,nrT-1);
            wF = nan*zeros(1,nrT);
            parfor t=1:nrT
                % t=2;
                tInd = setdiff(1:nrT,t);
                [Distr,freq] = X2RIP(obj,obj.XMOD,obj.RIPMOD,TestX(t));
                wF(t) = 1-freq;
                if isempty(Distr), continue; end
                m = PPMMOD.getStaticMatch(Distr,RIP,'QD');
                QDTMatches(:,t) = m(t);
                QDFMatches(:,t,:) = m(tInd);
                m = PPMMOD.getStaticMatch(Distr,RIP,'INLIER');
                INTMatches(:,t) = m(t);
                INFMatches(:,t,:) = m(tInd);
            end
        end
    end
    methods % INTERFACING
       function obj = reduceSamples(obj,index)
            if nargout==1, obj = clone(obj); end
            obj.X = obj.X(index);
            if~isempty(obj.MRIP), obj.MRIP = obj.MRIP(index); end
            if~isempty(obj.HMRIP), obj.HMRIP = obj.HMRIP(index); end
       end
       function out = getDepVar(obj,in)
           if ~iscell(in), out = in;return;end
           if ~obj.Level==0,out = in{obj.Level}.DepVar{obj.Cluster};return;end
           out = [];
           for l=1:1:length(in)
              for i=1:1:length(in{l}.DepVar)
                  out = [out, in{l}.DepVar{i}]; %#ok<AGROW>
              end
           end
       end
       function out = getSymSpace(obj,in)
           if ~iscell(in), out = in;return;end
           if ~obj.Level==0,out = in{obj.Level}.SymSpace{obj.Cluster};return;end
           out = [];% cannot retrieve a single SymSpace for level 0;
       end
       function out = getRIP(obj,DepVar)
           out = PPMMOD.updateRIP(DepVar,obj.MMOD);
           if obj.RIPNorm, out = PPMMOD.normalizeRIP(out,obj.MMOD,obj.VarInfo);end
           out = out';
       end
    end
    methods (Abstract = true)
        out = getVarInfo(obj,DepVar);
        [TrX,TrRIP,TestX,TestRIP,posclass] = prepOppBiometricTest(obj,TestRIP,TestX);
        [Distr,Freq,RIP] = X2RIP(obj,TrX,TrRIP,X);
        postBuild(obj);
    end
    methods (Static = true)
       function [R2,pR2] = testPartialPLSR(X,Y,C,t)
           index = intersect(PPMMOD.notNAN(C),PPMMOD.notNAN(X));
           if ~isempty(C)
            E = PPMMOD.getResiduals(C(index,:),Y(index,:));
            X = PPMMOD.getResiduals(C(index,:),X(index,:));
           else
            E = Y(index,:);
            X = X(index,:);
           end
           [~,~,~,~,~,var] = plsregress(X,E,1);R2 = var(2);
           if t==0, pR2 = nan; return; end
           R2Count = false(t,1);
           % incremental permutation
           for p=1:1:t
               pind = randperm(length(index));
               [~,~,~,~,~,forvar] = plsregress(X,E(pind,:),1);
               R2Count(p) = forvar(2)>=R2;
               if mod(p,100)==0
                  pR2 = (sum(R2Count(1:p))+1)/(p+1); 
                  acc = 10/p;
                  if pR2>acc, return;end
               end
           end
           pR2 = (sum(R2Count)+1)/(t+1);
       end 
       function [A,M,STAT] = AngleTest(IndVar,DepVar,AngleTable,t,maxM,Aval)
               if nargin<6,Aval=0;end
               if nargin<5,maxM=t;end
               A = nan*zeros(t,1);
               [nS,ShapeDim] = size(DepVar);
               M = nan*zeros(ShapeDim,2,maxM);
               el = 0;
               for i=1:t
                   el = el+1;
                   F = crossvalind('Kfold',nS,2);
                   forM = nan*zeros(ShapeDim,2);
                   for k=1:1:2
                       ind = find(F==k);
                       forM(:,k) = PPMMOD.getRegression(IndVar(ind,:),DepVar(ind,:));
                   end
                   A(i) = angle(forM(:,1),forM(:,2));
                   M(:,:,el) = forM;
                   if mod(i,maxM)==0, el = 0; end
                   if mod(i,100)==0
                      tmppA = (length(find(A(1:i)<=0))+1)/(i+1); 
                      acc = 10/i;
                      if tmppA>acc, break;end
                   end
               end
               STAT.AvgA = mean(A(1:i));
               STAT.pA = zeros(length(Aval),1);
               for j=1:1:length(Aval)
                   STAT.pA(j) = (length(find(A(1:i)<=Aval(j)))+1)/(i+1);  
               end
               STAT.pAR = PPMMOD.lookUppRA(STAT.AvgA,ShapeDim,AngleTable);
               STAT.MedA = median(A(1:i));
               STAT.StdA = std(A(1:i));
               STAT.MadA = mad(A(1:i));
               sortA = sort(A(1:i));
               STAT.UpperA = sortA(round(0.975*i));
               STAT.LowerA = sortA(round(0.025*i+1));
               if i<maxM,M = M(:,:,1:i);end
       end
       function M = getRegression(A,B)
           [A,B] = eliminateNAN(A,B);
           [~,~,~,~,M] = plsregress(A,B,size(A,2));
           M = M(end,:);
       end
       function out = getResiduals(X,Y)
                [~,~,~,~,~,~,~,stats] = plsregress(X,Y,min(size(X,2),size(Y,2)));
                out = stats.Yresiduals;
       end
       function out = lookUppRA(A,Dim,Table)
           out = interp1(Table.Angles,Table.pAngles(Dim,:),A);
       end
       function out = pRAWeight(A,Dim,Table)
                min = -log(BASEMOD.lookUppRA(0,Dim,Table));
                out = -log(BASEMOD.lookUppRA(A,Dim,Table));
                out = max(out-min,0);
       end
       function out = pAWeight(pA,minpA)
                min = -log(minpA);
                out = -log(pA);
                out = max(out-min,0);
       end
       function out = notNAN(in)
           index = (1:size(in,1));
           [i,~] = find(isnan(in));
           out = setdiff(index,unique(i));
       end
       function out = getL1CL1DepVar(in)
           if ~iscell(in), out = in;return;end
           out = in{1}.DepVar{1};
       end
       function H = wilcoxonM(M,pT)
            nrD = size(M,2);
            P = zeros(1,nrD);
            for j=1:1:nrD  
               P(j) = signrank(M(:,j));% wilcoxon test for median  
            end
            H = P<=pT;
       end
       function [out] = normalizeRIP(rip,M,var)
             Mrip = dot(var.MDepVar',M'/norm(M'));
             Prip = dot(var.PDepVar',M'/norm(M'));
             out = ((rip-Mrip)/(Prip-Mrip))*var.Range+var.Mel;
       end
       function [out] = updateRIP(in,M)
            out = dot(in',repmat(M'/norm(M'),1,size(in,1)));
       end
       function [EER,G,AUC,x,y,TH,Y] = getEER(tmatches,fmatches)
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
                %disp(num2str(T));
                %disp(num2str(1-sorted(ind)));
                TH(1) = 1-sorted(ind);
                TH(2) = 1-sorted(indy);
                %disp(['EER TH: ' num2str(TH(1))]);
                %disp(['Y TH: ' num2str(TH(2))]);
       end
       function out = getRANK(tmatches,fmatches)
                R = sum(repmat(tmatches,1,length(fmatches))>repmat(fmatches',length(tmatches),1),2);
                R = R./(length(fmatches)+1);
                RM = mean(R);
                R1 = sum(R<=0.01)/length(R);
                R10 = sum(R<=0.10)/length(R);
                R20 = sum(R<=0.20)/length(R);
                out = [R1 R10 R20 RM];
       end
       function [AVG,STD,MED,MAD] = getRIPStatistics(rip,kappa)
                  MED = median(rip);MAD = mad(rip);
                  if nargin<2,AVG = mean(rip);STD = std(rip);return;end 
                  nrO = length(rip);W = ones(1,nrO);
                  for i=1:1:10
                     AVG = sum(W.*rip)/sum(W); 
                     STD = sqrt(sum(W.*((rip-AVG).^2))/sum(W));
                     IP = normpdf(rip,AVG,STD);
                     L = (1/(sqrt(2*pi)*STD))*exp(-0.5*kappa^2);
                     W = IP./(IP+L);
                  end
       end
       function [AVG,STD] = getWeightedRIPStatistics(rip,W2,kappa)
                if nargin<3
                   AVG = sum(W2.*rip)/sum(W2);
                   STD = sqrt(sum(W2.*((rip-AVG).^2))/sum(W2));
                   return;
                end
                nrO = length(rip);
                W1 = ones(1,nrO);  
                for i=1:1:10
                    W = W1.*W2;
                    AVG = sum(W.*rip)/sum(W); 
                    STD = sqrt(sum(W.*((rip-AVG).^2))/sum(W));
                    IP = normpdf(rip,AVG,STD);
                    L = (1/(sqrt(2*pi)*STD))*exp(-0.5*kappa^2);
                    W1 = IP./(IP+L);
                end
       end
       function out = matchthresholding(matches,T)
                out = matches>=T; 
       end
       function out = class2score(classes)
                out = 1-classes; 
       end
       function out = getStaticMatch(distr,rip,type)
                if isempty(distr), out = nan*ones(size(rip));return;end
                pT = normpdf(rip,distr.Tmu,distr.Tsigma)';
                pF = normpdf(rip,distr.Fmu,distr.Fsigma)';
                switch type
                    case 'QD'
                        out = pT./pF;
                    case 'INLIER'
                        out = pT./(pT+pF);
                end  
       end
       function [PREC,REC,G,TNF,XH,YH] = getHardClass(Tclass,Fclass)
                Pc = length(Tclass);Nc = length(Fclass);
                TP = length(find(Tclass==1));
                TN = length(find(Fclass==0));TNF = TN/Nc;
                FN = Pc-TP;FP = Nc-TN;
                PREC = TP/(TP+FP);
                REC = TP/(TP+FN);
                G = sqrt((TP/(TP+FN))*(TN/(TN+FP)));
                YH = TP/Pc;
                XH = 1-TN/Nc;
       end
    end
end