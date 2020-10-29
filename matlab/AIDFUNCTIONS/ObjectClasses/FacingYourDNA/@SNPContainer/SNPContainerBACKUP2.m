classdef SNPContainer < superClass
   properties
      TrackID = 'X';
      SNPs = [];
   end
   properties (Hidden = true)
      ParCont = [];
   end
   properties (Dependent = true)
      %ShapeSpace;
      %RedShapeSpace;
      %Cov;
      %CovNames;
      nrSNP;
   end
   % Stage 0 QUALITY CHECKING
   properties (Hidden = true, Dependent = true)
      ST0MinSample;
   end
   % Stage 1 LOCATION BASE TESTS 
   properties (Hidden = true, Dependent = true)
      ST1PT;
      ST1MinHits;
   end
   % Stage 2 SNPBRIM BASE TESTS
   properties (Hidden = true, Dependent = true)
      ST2PFFoldT;
      ST2PGFoldT;
      ST2PFT;
      ST2PGT;
      ST2NrFoldT;
      ST2MinHits;
   end
   % Stage 3 SNPBRIM COVARIATES TEST
   properties (Hidden = true, Dependent = true)
      ST3PFFoldT;
      ST3PGFoldT;
      ST3PFT;
      ST3PGT;
      ST3NrFoldT;
   end
   properties (Hidden = true, Dependent = true) % SNPBRIM
      RegRuns; % Number of Regression runs (SNPBRIM)
      RegSampling;% Bootstrap/None (SNPBRIM)
      RegBalance; % balance the data yes or no (SNPBRIM)
      RegBalanceTH; % desired balance factor (SNPBRIM)
      RegBalanceMethod; % UpSample/DownSample/ADASYN (SNPBRIM)
      OuterFold; % number of Outer Folds (SNPBRIM)
      InnerFold; % number of Inner Folds (SNPBRIM)
      UseRedShape; % Use reduced space or covariates (SNPBRIM)
      MinFoldSampleSize; % Min samples per Outer Fold and Inner Fold (SNPBRIM)
      MaxIterations; % Maximum number of BRIM iterations (SNPBRIM)
      Htest; % Perform Wilcoxon test om partial regression coefficients (SNPBRIM)
      RIPNormalize; % normalize RIP values based on group averages (SNPBRIM)
      StopCorr; % Stopping correlation between subsequent iterations (SNPBRIM)
      RIPStatPerm; % Number of permutations in testing significance of RIP values (SNPBRIM)
      GroupDef;% Definition of groups, GG = according to Test performed, GT = according to orginal genotypes (SNPBRIM)
   end
   % SNPContainer ANALYSIS
   properties
       R2Total;
       MCorr;
       RIPCorr;
       SimilarT = 0.4;
       nrCL = 10;
       CL = [];
       CLM = [];
   end
   properties (Dependent = true);
       AvgMCorr;
       nrSimilar;
   end
   methods % Constructor
        function obj = SNPContainer(varargin)
            obj = obj@superClass(varargin{:});
            if isempty(obj.ParCont), obj.ParCont = ParContainer;end
        end
   end
   methods % General GETTING/SETTING
       function out = get.ParCont(obj)
           out = obj.ParCont;
           if ~superClass.isH(out), out = []; end
       end
       function out = get.nrSNP(obj)
           out = length(obj.SNPs);
       end  
   end
   methods % STAGE 0 GETTING/SETTING
       function out = get.ST0MinSample(obj)
           out = obj.ParCont.ST0MinSample;
       end
       function obj = set.ST0MinSample(obj,in)
           obj.ParCont.ST0MinSample = in; 
       end
   end
   methods % STAGE 1 GETTING/SETTING
       function out = get.ST1PT(obj)
           out = obj.ParCont.ST1PT;
       end
       function obj = set.ST1PT(obj,in)
           obj.ParCont.ST1PT = in; 
       end
       function out = get.ST1MinHits(obj)
           out = obj.ParCont.ST1MinHits;
       end
       function obj = set.ST1MinHits(obj,in)
           obj.ParCont.ST1MinHits = in; 
       end
   end
   methods % STAGE 2 GETTING/SETTING
       function out = get.ST2PFFoldT(obj)
           out = obj.ParCont.ST2PFFoldT;
       end
       function obj = set.ST2PFFoldT(obj,in)
           obj.ParCont.ST2PFFoldT = in; 
       end
       function out = get.ST2PGFoldT(obj)
           out = obj.ParCont.ST2PGFoldT;
       end
       function obj = set.ST2PGFoldT(obj,in)
           obj.ParCont.ST2PGFoldT = in; 
       end
       function out = get.ST2PFT(obj)
           out = obj.ParCont.ST2PFT;
       end
       function obj = set.ST2PFT(obj,in)
           obj.ParCont.ST2PFT = in; 
       end
       function out = get.ST2PGT(obj)
           out = obj.ParCont.ST2PGT;
       end
       function obj = set.ST2PGT(obj,in)
           obj.ParCont.ST2PGT = in; 
       end
       function out = get.ST2NrFoldT(obj)
           out = obj.ParCont.ST2NrFoldT;
       end
       function obj = set.ST2NrFoldT(obj,in)
           obj.ParCont.ST2NrFoldT = in; 
       end
       function out = get.ST2MinHits(obj)
           out = obj.ParCont.ST2MinHits;
       end
       function obj = set.ST2MinHits(obj,in)
           obj.ParCont.ST2MinHits = in; 
       end
   end
   methods % STAGE 3 GETTING/SETTING
       function out = get.ST3PFFoldT(obj)
           out = obj.ParCont.ST3PFFoldT;
       end
       function obj = set.ST3PFFoldT(obj,in)
           obj.ParCont.ST3PFFoldT = in; 
       end
       function out = get.ST3PGFoldT(obj)
           out = obj.ParCont.ST3PGFoldT;
       end
       function obj = set.ST3PGFoldT(obj,in)
           obj.ParCont.ST3PGFoldT = in; 
       end
       function out = get.ST3PFT(obj)
           out = obj.ParCont.ST3PFT;
       end
       function obj = set.ST3PFT(obj,in)
           obj.ParCont.ST3PFT = in; 
       end
       function out = get.ST3PGT(obj)
           out = obj.ParCont.ST3PGT;
       end
       function obj = set.ST3PGT(obj,in)
           obj.ParCont.ST3PGT = in; 
       end
       function out = get.ST3NrFoldT(obj)
           out = obj.ParCont.ST3NrFoldT;
       end
       function obj = set.ST3NrFoldT(obj,in)
           obj.ParCont.ST3NrFoldT = in; 
       end
   end
   methods % SNPBRIM GETTING/SETTING
       function out = get.RegRuns(obj)
           out = obj.ParCont.RegRuns;
       end
       function obj = set.RegRuns(obj,in)
           obj.ParCont.RegRuns = in; 
       end
       function out = get.RegSampling(obj)
           out = obj.ParCont.RegSampling;
       end
       function obj = set.RegSampling(obj,in)
           obj.ParCont.RegSampling = in; 
       end
       function out = get.RegBalance(obj)
           out = obj.ParCont.RegBalance;
       end
       function obj = set.RegBalance(obj,in)
           obj.ParCont.RegBalance = in; 
       end
       function out = get.RegBalanceTH(obj)
           out = obj.ParCont.RegBalanceTH;
       end
       function obj = set.RegBalanceTH(obj,in)
           obj.ParCont.RegBalanceTH = in; 
       end
       function out = get.RegBalanceMethod(obj)
           out = obj.ParCont.RegBalanceMethod;
       end
       function obj = set.RegBalanceMethod(obj,in)
           obj.ParCont.RegBalanceMethod = in; 
       end
       function out = get.OuterFold(obj)
           out = obj.ParCont.OuterFold;
       end
       function obj = set.OuterFold(obj,in)
           obj.ParCont.OuterFold = in; 
       end
       function out = get.InnerFold(obj)
           out = obj.ParCont.InnerFold;
       end
       function obj = set.InnerFold(obj,in)
           obj.ParCont.InnerFold = in; 
       end
       function out = get.UseRedShape(obj)
           out = obj.ParCont.UseRedShape;
       end
       function obj = set.UseRedShape(obj,in)
           obj.ParCont.UseRedShape = in; 
       end
       function out = get.MinFoldSampleSize(obj)
           out = obj.ParCont.MinFoldSampleSize;
       end
       function obj = set.MinFoldSampleSize(obj,in)
           obj.ParCont.MinFoldSampleSize = in; 
       end
       function out = get.MaxIterations(obj)
           out = obj.ParCont.MaxIterations;
       end
       function obj = set.MaxIterations(obj,in)
           obj.ParCont.MaxIterations = in; 
       end
       function out = get.Htest(obj)
           out = obj.ParCont.Htest;
       end
       function obj = set.Htest(obj,in)
           obj.ParCont.Htest = in; 
       end
       function out = get.RIPNormalize(obj)
           out = obj.ParCont.RIPNormalize;
       end
       function obj = set.RIPNormalize(obj,in)
           obj.ParCont.RIPNormalize = in; 
       end
       function out = get.StopCorr(obj)
           out = obj.ParCont.StopCorr;
       end
       function obj = set.StopCorr(obj,in)
           obj.ParCont.StopCorr = in; 
       end
       function out = get.RIPStatPerm(obj)
           out = obj.ParCont.RIPStatPerm;
       end
       function obj = set.RIPStatPerm(obj,in)
           obj.ParCont.RIPStatPerm = in; 
       end
       function out = get.GroupDef(obj)
           out = obj.ParCont.GroupDef;
       end
       function obj = set.GroupDef(obj,in)
           obj.ParCont.GroupDef = in; 
       end      
   end
   methods % ANALYSIS GETTING/SETTING
       function out = get.AvgMCorr(obj)
          if isempty(obj.MCorr), out = []; return; end
          out = triu(obj.MCorr);
          out = mean(abs(out(:)));
       end
       function out = get.nrSimilar(obj)
           if isempty(obj.MCorr), out = []; return; end
           out = obj.MCorr>=obj.SimilarT;
           out = sum(out,2)-1;
           out = (out/(obj.nrSNP-1))*100;
       end
   end
   methods % Interface functions
       function out = getRS(obj,index)
           if obj.nrSNP==0, out = []; return;end
           if nargin < 2, index = 1:obj.nrSNP; end
           out = cell(1,length(index));
           for i=1:1:length(index)
              if isempty(obj.SNPs{index(i)}), continue; end 
              if ~strcmp(obj.SNPs{index(i)}.Type,'SNP'), continue; end
              out{i} = obj.SNPs{index(i)}.RS;
           end
       end
       function out = getValues(obj,prop)
           if obj.nrSNP==0, out = []; return;end
           counter = 0;
           for i=1:1:obj.nrSNP % look for first valid SNP
               if ~isempty(obj.SNPs{i})&&strcmp(obj.SNPs{i}.Type,'SNP')
                   tmp = eval(['obj.SNPs{i}.' prop]);
                   break;
               end
               counter = counter + 1;
           end
           if counter==obj.nrSNP, out = []; return; end% there were no valid SNPs
           out = nan*zeros(size(tmp,1),size(tmp,2),obj.nrSNP);
           for j=counter+1:1:obj.nrSNP
              if isempty(obj.SNPs{j}), continue; end 
              if ~strcmp(obj.SNPs{j}.Type,'SNP'), continue; end
              out(:,:,j) = eval(['obj.SNPs{j}.' prop]);
           end
           out = squeeze(out);
       end
       function parseSNPs(obj,Names,GT)
          nrSNP = length(Names);
          obj.SNPs = cell(1,nrSNP);
          disp('PARSING SNPS');
          %f = waitbar(0,'Parsing SNPs');drawnow;
          for i=1:1:nrSNP
            snp = SNP;snp.RS = Names{i};snp.GT = GT(:,i);
            snp.ParCont = obj.ParCont;
            obj.SNPs{i} = snp;
            %waitbar(i/nrSNP,f);drawnow;
          end
          %delete(f);drawnow;
       end
       function runStage(obj,stage,BaseCont,t,varargin)
           nrSNP = obj.nrSNP;
           if nrSNP==0, return; end
           obj.RIPStatPerm = t;
           [path,ID] = setupParForProgress(nrSNP);
           parcont = clone(obj.ParCont);
           unlinkContainer(obj);snps = obj.SNPs;
           done = cell(1,nrSNP);
           % preparing minimum data transfer
           switch stage
               case 1
                   inputdata.DepVar = BaseCont.RedDepVar;
               case 2
                   inputdata.DepVar = BaseCont.DepVar;
                   inputdata.Cov = [BaseCont.CovRIP, BaseCont.GBRIP];
               case 3
                   inputdata.BaseCont = BaseCont;
                   inputdata.UseRedShape = varargin{1};
           end
           parfor i=1:nrSNP
               switch stage
                   case 1
                       done{i} = runStage1(snps{i},t,inputdata.DepVar,parcont); %#ok<*PFBNS>
                   case 2
                       done{i} = runStage2(snps{i},inputdata.DepVar,inputdata.Cov,parcont);
                   case 3
                       done{i} = runStage3(snps{i},inputdata.BaseCont,inputdata.UseRedShape,parcont);
               end
               parfor_progress;
           end
           closeParForProgress(path,ID);
           obj.SNPs = done;linkContainer(obj);
       end
       function prepareStage3Splitting(obj,pcuttof)
           testindex = [1 2 3 4 5 6 7];
           new = cell(1,length(testindex)*obj.nrSNP);
           counter = 0;
           %f = waitbar(0,'SPLITTING SNPs');drawnow;
           for i=1:1:obj.nrSNP
               forsnp = obj.SNPs{i};
               % each SNP is copied X times, X = number of tests that it survived in stage 2
               for j=1:1:length(testindex)
                  if forsnp.ST2P(testindex(j))<=pcuttof% Does it survive the test
                     counter = counter + 1;
                     tmp = clone(forsnp);
                     tmp.ST3TestInd = testindex(j);
                     tmp.ST3F = tmp.ST2F(testindex(j));
                     tmp.ST3P = tmp.ST2P(testindex(j));
                     tmp.ST3FFold = tmp.ST2FFold(:,testindex(j));
                     tmp.ST3PFold = tmp.ST2PFold(:,testindex(j));
                     new{counter} = tmp;
                  end
               end
               %waitbar(i/obj.nrSNP,f);drawnow;
           end
           %delete(f);drawnow;
           disp('DELETING');
           deleteSNPs(obj);
           disp('ADDING');
           obj.SNPs = new(1:counter);linkContainer(obj);
       end
       function obj = reduceSNPs(obj,index)
          unlinkContainer(obj);
          if nargout==1, cobj = obj; obj = clone(cobj);linkContainer(cobj);end
          tmpSNPs = obj.SNPs(index);
          obj.SNPs = tmpSNPs;
          linkContainer(obj);
       end
       function addSNPs(obj,new)
           obj.SNPs = [obj.SNPs new];
           linkContainer(obj);
       end
       function obj = mergeContainers(obj,addcont)
                [ind12,~] = vlookup(getRS(obj),getRS(addcont));
                addindex = setdiff(1:addcont.nrSNP,ind12);
                unlinkContainer(obj);unlinkContainer(addcont);
                if nargout==1, cobj = obj; obj = clone(cobj);linkContainer(cobj);end
                tmpsnps = cell(1,obj.nrSNP+length(addindex));
                tmpsnps(1:obj.nrSNP) = obj.SNPs;
                f = waitbar(0,'Merging');drawnow;
                for i=1:length(addindex)
                    tmpsnps{obj.nrSNP+i} = clone(addcont.SNPs{addindex(i)});
                    waitbar(i/obj.nrSNP,f);drawnow;
                end
                obj.SNPs = tmpsnps;
                linkContainer(obj);linkContainer(addcont);
                delete(f);drawnow;
       end
       function linkContainer(obj) 
          nrSNP = obj.nrSNP;
          for i=1:1:nrSNP
            if isempty(obj.SNPs{i}), continue; end
            if ~strcmp(obj.SNPs{i}.Type,'SNP'), continue; end  
            obj.SNPs{i}.ParCont = obj.ParCont;
          end
       end
       function unlinkContainer(obj) 
          nrSNP = obj.nrSNP;
          for i=1:1:nrSNP
            if isempty(obj.SNPs{i}), continue; end
            if ~strcmp(obj.SNPs{i}.Type,'SNP'), continue; end  
            obj.SNPs{i}.ParCont = [];
          end
       end
       function deleteSNPs(obj)
          %f = waitbar(0,'Deleting SNPs');drawnow;  
          nrSNP = obj.nrSNP;
          for i=1:1:nrSNP
            if isempty(obj.SNPs{i}), continue; end
            if ~strcmp(obj.SNPs{i}.Type,'SNP'), continue; end  
            obj.SNPs{i}.ParCont = [];
            delete(obj.SNPs{i});
            %waitbar(i/nrSNP,f);drawnow;
          end
          obj.SNPs = {};
          %delete(f);drawnow;
       end
       function splitAllSNPsForST3(obj)
           testindex = [1 2 3 5 6 7];
           new = cell(1,length(testindex)*obj.nrSNP);
           counter = 0;
           f = waitbar(0,'SPLITTING SNPs');drawnow;
           for i=1:1:obj.nrSNP
               for j=1:1:length(testindex) % each SNP is copied X times, X = number of tests
                  counter = counter + 1; 
                  tmp = clone(obj.SNPs{i});
                  tmp.ST3TestInd = testindex(j);
                  new{counter} = tmp;
               end
               waitbar(i/obj.nrSNP,f);drawnow;
           end
           delete(f);drawnow;
           deleteSNPs(obj);obj.SNPs = new;linkContainer(obj);
       end
       function delete(obj)
                deleteSNPs(obj);
                obj.SNPs = {};
                delete@superClass(obj);
       end
       function [IB,pTIB,pT] = matchFaces(obj,faces,genotypes,BF,ind)    
                nrSNP = obj.nrSNP;nrFaces = size(faces,1);
                if nargin<5, ind = (1:nrSNP); end
                snps = obj.SNPs;
                IB = nan*zeros(length(ind),nrFaces);
                pTIB = nan*zeros(length(ind),nrFaces);
                pT = nan*zeros(length(ind),nrFaces);
                for i=1:1:length(ind)
                   [IB(i,:),pTIB(i,:),pT(i,:)] = matchFaces(snps{i},faces,genotypes(i),BF);  
                end
       end
       function [RIP,Distr] = GT2RIP(obj,genotypes)
                RIP = nan*zeros(1,obj.nrSNP);Distr = cell(1,obj.nrSNP);
                for i=1:1:obj.nrSNP
                    [RIP(i),Distr(i)] = GT2RIP(obj.SNPs{i},genotypes(i)); 
                end
       end
       function out = getGTFreq(obj,genotypes)
                out = nan*zeros(1,obj.nrSNP);
                for s=1:1:obj.nrSNP
                  out(s) = GTTypicality(obj.SNPs{s},genotypes(s));
                end
       end
       function obj = reduceSamples(obj,index)
                if nargout == 1, obj = clone(obj); end
                for i=1:1:obj.nrSNP
                   reduceSamples(obj.SNPs{i},index);
                end
                linkContainer(obj);
       end
   end
   methods % SNPCONTAINER ANALYSIS
       function getR2Total(obj)
           for i=1:1:obj.nrSNP
              if i==1, obj.R2Total = obj.SNPs{i}.ST4MOD.R2; continue; end
              obj.R2Total = obj.R2Total+obj.SNPs{i}.ST4MOD.R2;
           end
           obj.R2Total = obj.R2Total./obj.nrSNP;
       end
       function illustrateR2Total(obj,refscan)
          for i = 1:1:3
              scan = clone(refscan);
              scan.Value = obj.R2Total(i,:);
              scan.ColorMode = 'Indexed';
              v = viewer(scan);
              colormap(v.RenderAxes,'jet');
              set(v.RenderAxes,'clim',[0 0.03]);
              v.BackgroundColor = [1 1 1];
          end
       end
       function getMCorr(obj,dim,nTR)
           [M,RIP] = getMAndRIP(obj,dim,nTR);
           obj.MCorr = abs(-1*(squareform(pdist(M,'cosine'))-1));
           obj.RIPCorr = abs(-1*(squareform(pdist(M,'correlation'))-1));
           %obj.MCorr = nan*zeros(obj.nrSNP,obj.nrSNP);
           %for i=1:1:obj.nrSNP
           %    for j=1:1:obj.nrSNP
                   %obj.MCorr(i,j) = vectorCorr(obj.SNPs{i}.ST4MOD.M',obj.SNPs{j}.ST4MOD.M');
           %        obj.MCorr(i,j) = angle(obj.SNPs{i}.ST4MOD.M',obj.SNPs{j}.ST4MOD.M');
           %    end
           %end
           %obj.MCorr = abs((obj.MCorr+obj.MCorr')/2);
       end
       function illustrateMCorr(obj)
           obj.AvgMCorr
           figure;imagesc(abs(obj.MCorr));set(gca,'clim',[0.3 1]);colormap(gca,'jet');
           figure;hist(obj.nrSimilar);
       end
       function clusterSNP(obj)
           obj.CL = mySpecClustering(obj.MCorr,obj.nrCL);
       end
       function out = clusterMembers(obj)
           out = zeros(1,obj.nrCL);
           for i=1:1:obj.nrCL
               out(i) = length(find(obj.CL==i));
           end
           obj.CLM = out;
       end
       function out = extractRelativeCW(obj)
           out = ones(obj.nrSNP,1);
           for i=1:1:obj.nrCL
              out(obj.CL==i) = 1/obj.CLM(i); 
           end
       end
       function out = extractRelativeCorrW(obj,cl)
           out = obj.MCorr>=cl;
           out = sum(out,2);
           out = 1./out;
       end
       function [M,RIP] = getMAndRIP(obj,dim,nTR)
           M = nan*zeros(obj.nrSNP,dim);
           RIP = nan*zeros(obj.nrSNP,nTR);
           for i=1:1:obj.nrSNP
               M(i,:) = obj.SNPs{i}.ST4MOD.M;
               RIP(i,:) = obj.SNPs{i}.ST4MOD.CompRIP;
           end 
       end
   end
   methods (Static = true)
       function obj = loadobj(struc)
            disp('Loading SNP Container');
            obj = struc2obj(eval(struc.Type),struc);
            % relinking SNP models
            linkContainer(obj);
       end
   end
end

%        function runStage(obj,stage,t)
%            if obj.nrSNP==0, return; end
%            if nargin<3, t = 1000; end
%            nrSNP = obj.nrSNP;
%            f = waitbar(0,'Processing SNPs');drawnow;
%            snps = obj.SNPs;parcont = obj.ParCont;
%            done = cell(1,nrSNP);
%            for i=1:nrSNP
%                switch stage
%                    case 1
%                        done{i} = runStage1(snps{i},t,parcont); %#ok<*PFBNS>
%                    case 2
%                        done{i} = runStage2(snps{i},parcont);
%                    case 3
%                        done{i} = runStage3(snps{i},parcont);
%                    case 4
%                        done{i} = runStage4(snps{i},parcont);
%                end
%                waitbar(i/nrSNP,f);drawnow;
%            end
%            delete(f);drawnow;
%            obj.SNPs = done;linkContainer(obj);
%        end




%        function out = get.ShapeSpace(obj)
%            out = obj.ParCont.ShapeSpace;
%        end
%        function obj = set.ShapeSpace(obj,in)
%            obj.ParCont.ShapeSpace = in; 
%        end
%        function out = get.RedShapeSpace(obj)
%            out = obj.ParCont.RedShapeSpace;
%        end
%        function obj = set.RedShapeSpace(obj,in)
%            obj.ParCont.RedShapeSpace = in; 
%        end

%        function out = get.Cov(obj)
%            out = obj.ParCont.Cov;
%        end
%        function obj = set.Cov(obj,in)
%            obj.ParCont.Cov = in; 
%        end
%        function out = get.CovNames(obj)
%            out = obj.ParCont.CovNames;
%        end
%        function obj = set.CovNames(obj,in)
%            obj.ParCont.CovNames = in; 
%        end

%        function createReducedSpace(obj)
%                 createReducedSpace(obj.ParCont);
%        end


