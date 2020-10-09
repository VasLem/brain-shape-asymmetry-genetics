function [obj,out] = NextGenBRIM(obj,BootInd,CondInd,PartInd)
% [obj,NIndVar,Iter] = NextGenBRIM(obj,BRIMindex,Condindex,Partindex)
% BootInd, index of variables to bootstrap (BRIM)
% CondInd, index of conditioning variables (not to BRIM)
% PartInd, index of the variable based upon which to partition the data,
% when empty default kfold paritioting will take place
if nargin < 4,PartInd = [];end
if ~checkNrObservations(obj), error('different amount of observations between X and Y'); end
if obj.MaxIterations == 0
   Bootstrap = false; 
else
   Bootstrap = true;  
end
% Building Independent and Dependent Variables
    n = obj.nrO;Ind = (1:n);
    nr2Condition = length(CondInd);
    var.nr2Boot = length(BootInd);
    obj.RIPX(:,CondInd) = obj.X(:,CondInd);% copying the ones not to BRIM
    IndVar = obj.X(:,[CondInd BootInd]);
    var.Booting = (nr2Condition+1:nr2Condition+var.nr2Boot);
% Partitioning Variable used to define the Folds
    if isempty(PartInd)
        G = ones(n,1);
    else
        G = obj.X(:,PartInd);
    end
% Dependent Shape Variables
    switch obj.ShapeDistance
        case 'Mahalanobis'
            DepVar = obj.Y./repmat(obj.Model.EigStd',n,1);
        case 'Euclidean'
            DepVar = obj.Y;
    end
    out.IndVar = IndVar;
    out.DepVar = DepVar;
% examining variables to Boot and create reference faces, the lattter is
% needed to rescale RIP variables onto common scales
    for i=1:1:var.nr2Boot
       v = IndVar(:,var.Booting(i));
       el = unique(v(~isnan(v)));
       if length(el)<6%Categorical parameter, currently allowing up to 6 categories  
          var.Info{i}.Type = 'Categorical';
          var.Info{i}.El = el;var.Info{i}.Mel = min(el);var.Info{i}.Pel = max(el);
          var.Info{i}.Range = var.Info{i}.Pel-var.Info{i}.Mel;
          % getting categorical reference faces
          var.Info{i}.Mindex = find(v==var.Info{i}.Mel);
          var.Info{i}.MDepVar = mean(DepVar(var.Info{i}.Mindex,:)); %#ok<*FNDSB>
          var.Info{i}.Pindex = find(v==var.Info{i}.Pel);
          var.Info{i}.PDepVar = mean(DepVar(var.Info{i}.Pindex,:));
       else% variable is continous
          var.Info{i}.Type = 'Continous';
          index = find(~isnan(v));
          [sortv,ind] = sort(v(index),'ascend');
          % lower 25%
          indl = round(0.25*length(ind));
          var.Info{i}.Mb = sortv(indl);
          var.Info{i}.Mindex = find(v<=var.Info{i}.Mb);
          var.Info{i}.Mel = mean(v(var.Info{i}.Mindex));
          var.Info{i}.MDepVar = mean(DepVar(var.Info{i}.Mindex,:));
          % upper 25%
          indu = round(0.75*length(ind));
          var.Info{i}.Pb = sortv(indu);
          var.Info{i}.Pindex = find(v>=var.Info{i}.Pb);
          var.Info{i}.Pel = mean(v(var.Info{i}.Pindex));
          var.Info{i}.PDepVar = mean(DepVar(var.Info{i}.Pindex,:));
          var.Info{i}.Range = var.Info{i}.Pel-var.Info{i}.Mel;
       end     
    end
    out.Var = var;
% Outer partitioning of data
    Fouter = DOB_SCV(obj.OuterFold,DepVar,G); 
    FoldResults = cell(1,obj.OuterFold);
    disp('Starting OuterFolds...');
    %parfor_progress(obj.OuterFold);
    parfor fo=1:obj.OuterFold
        FoTestInd = find(Fouter==fo);
        FoldResults{fo}.TestInd = FoTestInd;
        FoTrInd = setdiff(Ind,FoTestInd);
        FoldResults{fo}.TrInd = FoTrInd;
        FoTrIndVar = IndVar(FoTrInd,:); %#ok<*PFBNS>
        FoTrDepVar = DepVar(FoTrInd,:);
        FoTestIndVar = IndVar(FoTestInd,:);
        FoTestDepVar = DepVar(FoTestInd,:);
        FoldResults{fo}.DepVar = FoTestDepVar;
        FoTrOrigIndVar = FoTrIndVar;
        FoTrNewIndVar = FoTrIndVar;
        BootProgress = zeros(var.nr2Boot,obj.MaxIterations,2);
        ContBoot = true;counter = 0;
        while ContBoot && Bootstrap>0
           % seperate data into inner folds
           Finner = DOB_SCV(obj.InnerFold,FoTrDepVar,G(FoTrInd));
           counter = counter + 1;
           for fi=1:obj.InnerFold
                FiTestInd = find(Finner==fi);
                FiTrInd = setdiff(1:length(FoTrInd),FiTestInd);
                [M,H] = getBootSampleRegression(FoTrIndVar(FiTrInd,:),FoTrDepVar(FiTrInd,:),var,obj.SamplingRuns);
                if obj.Htest, M = H.*M;end% setting close to zero partial regression coefficients to zero
                rip = updateRIP(FoTrDepVar(FiTestInd,:),M)';
                if obj.RIPNormalize, rip = normalizeRIP(rip,M,var);end
                FoTrNewIndVar(FiTestInd,var.Booting) = rip;
           end
           % monitoring progress of Booting
           ToStop = zeros(1,var.nr2Boot);
           for b=1:1:var.nr2Boot
               V1 = FoTrIndVar(:,var.Booting(b));
               V2 = FoTrNewIndVar(:,var.Booting(b));
               index =  ~isnan(V1);
               tmp = corrcoef(V1(index),V2(index));
               BootProgress(b,counter,2) = tmp(1,2);
               if tmp(1,2)>=obj.StopCorr, ToStop(b) = 1;end
               V1 = FoTrOrigIndVar(:,var.Booting(b));
               index =  ~isnan(V1);
               tmp = corrcoef(V1(index),V2(index));
               BootProgress(b,counter,1) = tmp(1,2);
           end
           if (sum(ToStop)/var.nr2Boot==1), ContBoot = false; end
           if counter >= obj.MaxIterations, ContBoot = false; end
           FoTrIndVar = FoTrNewIndVar;
        end
        if counter>0, FoldResults{fo}.BootProgress = BootProgress(:,1:counter,:);end
        FoldResults{fo}.BootIter = counter;
        % getting the final M from outer training data
        [M,H] = getBootSampleRegression(FoTrIndVar,FoTrDepVar,var,obj.SamplingRuns);
        if obj.Htest, M = H.*M;end
        FoldResults{fo}.M = M;FoldResults{fo}.H = H;
        % estimating rip values outer test data
        rip = updateRIP(FoTestDepVar,M)';
        if obj.RIPNormalize, rip = normalizeRIP(rip,M,var);end
        FoldResults{fo}.RIP = rip;
        FoldResults{fo}.IndVar = FoTestIndVar(:,var.Booting);
        % perform within outer fold statistics
        [FoldResults{fo}.stat,FoldResults{fo}.statType,FoldResults{fo}.statP] = performStat(FoldResults{fo}.IndVar,rip,var,1000);
        %parfor_progress;
    end
    %parfor_progress(0);
    disp('CrossFold Testing...');
    % crossfold results
    Stat = zeros(length(var.Booting),obj.OuterFold);
    StatP = zeros(length(var.Booting),obj.OuterFold);
    StatType = cell(length(var.Booting),obj.OuterFold);
    BootIter = zeros(1,obj.OuterFold);
    M = zeros(var.nr2Boot,size(DepVar,2),obj.OuterFold);
    H = zeros(var.nr2Boot,size(DepVar,2),obj.OuterFold);
    CM = zeros(obj.OuterFold,obj.OuterFold,var.nr2Boot);
    CStat = zeros(obj.OuterFold,obj.OuterFold,var.nr2Boot);  
    CRIP = zeros(n,var.nr2Boot);
    for i=1:obj.OuterFold% does not fit in parfor loop because of indexing
        CRIP(FoldResults{i}.TestInd,:) = FoldResults{i}.RIP;      
    end
    %parfor_progress(obj.OuterFold);
    parfor i=1:obj.OuterFold
        CMfor = zeros(obj.OuterFold,length(var.Booting));
        CStatfor = zeros(obj.OuterFold,length(var.Booting));
        Stat(:,i) = FoldResults{i}.stat;
        StatType(:,i) = FoldResults{i}.statType;
        StatP(:,i) = FoldResults{i}.statP;
        M(:,:,i) = FoldResults{i}.M;
        H(:,:,i) = FoldResults{i}.H;
        BootIter(i) = FoldResults{i}.BootIter;
        M1 = FoldResults{i}.M;
        for j=1:1:obj.OuterFold
            M2 = FoldResults{j}.M;
            % correlation between regression paths
            for k=1:1:length(var.Booting)
                CMfor(j,k) = vectorCorr(M1(k,:)',M2(k,:)');
            end
            % cross fold statistics
            rip = updateRIP(FoldResults{j}.DepVar,M1)';
            if obj.RIPNormalize, rip = normalizeRIP(rip,M1,var);end
            stat = performStat(FoldResults{j}.IndVar,rip,var,1);
            CStatfor (j,:)= stat;
        end 
        CM(i,:,:) = CMfor;
        CStat(i,:,:) = CStatfor;
        %parfor_progress;
    end
    %parfor_progress(0);
    out.Stat = Stat;out.StatType = StatType;out.StatP = StatP;
    out.StatP(out.StatP==0) = 0.001;
%     out.FisherP = ones(1,var.nr2Boot);
%     out.FisherChi = ones(1,var.nr2Boot);
    out.pfast = ones(1,var.nr2Boot);
    for i=1:1:var.nr2Boot
        %[out.FisherChi(i),out.FisherP(i)] = myFisher(out.StatP(i,:));
        out.pfast(i) = pfast(out.StatP(i,:));
    end
    out.CM = CM;out.CStat = CStat;out.M = M;out.H = H;out.BootIter = BootIter;
    % performing overall statistics
    [stat,statType,statP] = performStat(IndVar(:,var.Booting),CRIP,var,1000);
    out.CRIP = CRIP;
    out.TotStatType = statType;out.TotStat = stat;out.TotStatP = statP;
    % extracting final M
    out.MedM = median(M,3);%out.MedH = median(H,3);
    rip = updateRIP(DepVar,out.MedM)';
    if obj.RIPNormalize, rip = normalizeRIP(rip,out.MedM,var);end
    out.MedRIP = rip;
    if obj.OuterFold>=8% when there are not enough folds this test cannot be performed
        out.H = wilcoxonM(M,var);
    else
        %disp('Not enough folds skipping H testing');
        out.H = ones(size(out.MedM));
    end
    out.HMedM = out.H.*out.MedM;
    rip = updateRIP(DepVar,out.HMedM)';
    if obj.RIPNormalize, rip = normalizeRIP(rip,out.HMedM,var);end
    out.HMedRIP = rip;
    if obj.Show
       nrPC = size(DepVar,2);
       for i=1:1:var.nr2Boot
           Mfor = squeeze(M(i,:,:));
           f = figure;set(f,'Name',obj.XNames{BootInd(i)});hold on;
           boxplot(Mfor');
           plot(1:nrPC,out.HMedM(i,:),'mo');
           plot(1:nrPC,zeros(1,nrPC),'k-'); 
       end  
    end
    drawnow;
    if obj.CreateMorphs
        disp('Creating Morphs');
       % morphs for the overall path
       [out.HMedM_FminMorphs,out.HMedM_RegMorphs] = extractMorphs(out.HMedM,out.DepVar,out.Var,out.IndVar,obj.BF,obj.Model);
       % morphs from the different folds
       out.Fold_FminMorphs = cell(obj.OuterFold,2,var.nr2Boot);
       out.Fold_RegMorphs = cell(obj.OuterFold,2,var.nr2Boot);
       for i=1:1:obj.OuterFold
           [out.Fold_FminMorphs(i,:,:),out.Fold_RegMorphs(i,:,:)] = extractMorphs(out.M(:,:,i),out.DepVar,out.Var,out.IndVar,obj.BF,obj.Model);
       end
    end
    % linking results to BRIMShapeModel
    obj.RIPX(:,BootInd) = rip;
    disp('Done');
end
%% HELP FUNCTIONS
function [M,R2] = getRegression(A,B,Booting)
         if nargin < 3, Booting = (1:size(A,2));end
         [A,B] = eliminateNAN(A,B);
         [~,~,~,~,M] = plsregress(A,B,size(A,2));
         M = M(1+Booting,:);
         if nargout < 2, return; end
         R2 = zeros(size(M));
         for i=1:1:length(Booting)
             % building the reduced model
             index = setdiff(1:size(A,2),Booting(i));
             AR = A(:,index);
             [~,~,~,~,~,~,~,statsR] = plsregress(AR,B,size(AR,2));
             E = statsR.Yresiduals;
             [~,~,~,~,~,~,~,statsR] = plsregress(AR,A(:,Booting(i)));
             Afor = statsR.Yresiduals;
             [~,~,~,~,~,~,~,stats] = plsregress(Afor,E,1);
             [n,~] = size(E);
             P = E-stats.Yresiduals;
             avg = repmat(mean(E),n,1);
             SST = sum((E-avg).^2);
             SSR = sum((P-avg).^2);
             R2(i,:) = (SSR./SST);
         end
end
function [out] = updateRIP(in,M)
            % in this implementation I take the reference as the origin
            n2 = size(in,1); % determine input size
            in = in';
            n1 = size(M,1);
            out = nan*zeros(n1,n2);% allocate memory
            for j=1:n1
                coeff2 = M(j,:)';
                coeff2 = coeff2/norm(coeff2);
                out(j,:) = dot(in,repmat(coeff2,1,n2));
            end
end
function [stat,statType,statP] = performStat(vIn,vripIn,var,t)
        stat = zeros(1,length(var.Booting));
        statType = cell(1,length(var.Booting));
        statP = zeros(1,length(var.Booting));
        for s=1:1:length(var.Booting)
           v = vIn(:,s);
           vrip = vripIn(:,s);
           index = find(~isnan(v));
           switch var.Info{s}.Type
               case 'Categorical'
                   statType{s} = 'F';
                   XG = cell(size(v(index)));
                   for l=1:length(var.Info{s}.El)
                       XG((v(index)==var.Info{s}.El(l))) = {num2str(l)};
                   end
                   [stat(s),~,statP(s)] = myPermAnova(XG,vrip(index),t);
               case 'Continous'
                   statType{s} = 'Cor';
                   [stat(s),statP(s)] = permCorr(v(index),vrip(index),t);
           end
        end  
end
function [out] = normalizeRIP(rip,M,var)
     out = rip;
     for b= 1:1:var.nr2Boot % testing reference faces and rescale accordingly
        Mrip = updateRIP(var.Info{b}.MDepVar,M)';
        Prip = updateRIP(var.Info{b}.PDepVar,M)';
        tmp = rip(:,b);
        range = Prip(b)-Mrip(b);
        tmp = (tmp-Mrip(b))/range;           
        out(:,b) = tmp*var.Info{b}.Range+var.Info{b}.Mel;  
     end
end
function [M,H] = getBootSampleRegression(IndVar,DepVar,var,runs)
        nrS = size(IndVar,1);
        M = zeros(var.nr2Boot,size(DepVar,2),runs);
        for s=1:1:runs
            if runs==1% no random sampling
               SampleFi = 1:nrS; 
            else % random sampling with replacement
               SampleFi = randsample(nrS,nrS,true);
            end
            M(:,:,s) = getRegression(IndVar(SampleFi,:),DepVar(SampleFi,:),var.Booting);
        end
        if runs==1,M = squeeze(M);H = ones(size(M));return;end
        H = wilcoxonM(M,var);
        M = median(M,3);       
end
function H = wilcoxonM(M,var)
    H = zeros(var.nr2Boot,size(M,2));
    for i=1:1:var.nr2Boot
       Mfor = squeeze(M(i,:,:)); 
       for j=1:1:size(Mfor,1)
          [~,H(i,j)] = signrank(Mfor(j,:));    
       end
    end
end
function [FminMorphs,RegMorphs] = extractMorphs(M,Pop,var,IndVar,BF,shape)
          FminMorphs = cell(2,var.nr2Boot);
          RegMorphs = cell(2,var.nr2Boot);
          RIP = updateRIP(Pop,M)';
          IndVar(:,var.Booting) = RIP;
          Mreg = getRegression(IndVar,Pop,(1:size(IndVar,2)));
          avgIndVar = nanmean(IndVar);
          avgY = nanmean(Pop);
          for i=1:1:var.nr2Boot
              avgRIP = nanmean(RIP(:,i));
              stdRIP = nanstd(RIP(:,i));
              Mval = avgRIP-BF(i)*stdRIP;
              Pval = avgRIP+BF(i)*stdRIP;
              % fmin version
              C = getFminMorph(M(i,:),Mval);
              FminMorphs{1,i} = getScan(shape,C);
              C = getFminMorph(M(i,:),Pval);
              FminMorphs{2,i} = getScan(shape,C);
              % regression version
              X = avgIndVar;
              X(var.Booting(i)) = Mval;
              deltaX = X-avgIndVar;
              dY = deltaX*Mreg;
              C = avgY+dY;
              RegMorphs{1,i} = getScan(shape,C);
              X(var.Booting(i)) = Pval;
              deltaX = X-avgIndVar;
              dY = deltaX*Mreg;
              C = avgY+dY;
              RegMorphs{2,i} = getScan(shape,C);
          end
end
function [out,optX] = getFminMorph(M,val)      
         optX = fminsearch(@(X) MorphError(X,M,val),val);
         out = optX*M;
end
function out = MorphError(X,M,val)        
         C = X*M;
         test = updateRIP(C,M);
         out = abs(test-val);
end