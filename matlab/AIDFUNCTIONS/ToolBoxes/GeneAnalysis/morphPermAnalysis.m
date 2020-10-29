function out = morphPermAnalysis(A,shape,RefScan,BF,CondValues,t)
           RefScan.SingleColor = [0.8 0.8 0.8];
         % multiple regression to obtain multiple R-Squared
            [~,~,~,~,~,pctvar] = plsregress(A,shape,size(A,2));
            R2Tot = sum(pctvar(2,:));
            STATS.R2Tot = R2Tot;
         % Conditioning regression to obtain reference face according to CondValues
            AC = A(:,1:end-1);
            [~,~,~,~,MC,pctvarC,~,statsC] = plsregress(AC,shape,size(AC,2));
            E = statsC.Yresiduals;
            mC = mean(AC);
            mShape = mean(shape);
            if isempty(CondValues)
               CondValues = mC; 
            end
            DX = mC-CondValues;
            DY = DX*MC(2:end,:);
            NY = mShape-DY;
            scanC = clone(RefScan);
            delete(scanC.PoseLM);
            scanC.Vertices = reshape(NY,3,length(NY)/3);
            R2Red = sum(pctvarC(2,:));
            R2Add = R2Tot-R2Red;
            STATS.R2Add = R2Add;
         % reduced model construction   
            AR = A(:,end);
            [~,~,~,~,~,~,~,statsR] = plsregress(AC,AR,1);
            AR = statsR.Yresiduals;
            [~,~,~,~,MR,pctvar,~,stats] = plsregress(AR,E,1);
            R2Part = pctvar(2);
            STATS.R2Part = R2Part;
            [LocalE,LocalS] = getLocalStatistic(MR,stats.Yresiduals,E,'shape');
            STATS.LocalE = LocalE;
            STATS.LocalS = LocalS;
          % creating morphs in reduced model space   
            Y = mean(E);
            m = mean(AR);mv = std(AR);
            range(1) = BF*(m-3*mv);
            range(2) = BF*(m+3*mv);
            counter = 0;
            for k=1:1:2
                V = range(k);
                counter = counter+1;
                DX = m-V;
                DY = DX*MR(2:end,:);
                NY = Y-DY;
                scan = clone(RefScan);
                delete(scan.PoseLM);
                scan.Vertices = scanC.Vertices + reshape(NY,3,length(NY)/3);
                eval(['morph' num2str(counter) ' = clone(scan);']);
                if ~isempty(RefScan.PoseLM)
                    eval(['transferPoseLM(morph' num2str(counter) ',RefScan);']);
                end
            end
            out.Morph1 = morph1;
            out.Morph2 = morph2;
          % compare morphs
            STATS.comparison = compareMorphs(morph1,morph2);
            if t<1, out.STATS = STATS; return; end
          % Initialize perm counters
           % effect size and partial effect
            LocalSCount = false(scanC.nrV,t);
            R2PartCount = false(t,1);
           % Distances along the normal
            H1PosNormalDistCount = false(scanC.nrV,t);
            H1NegNormalDistCount = false(scanC.nrV,t);
            H2NormalDistCount = false(scanC.nrV,t);
           % Local Area ratio
            H1PosAreaCount = false(scanC.nrV,t);
            H1NegAreaCount = false(scanC.nrV,t);
            H2AreaCount = false(scanC.nrV,t);
           % Local Curvature ratio
            H1PosCurvRatioCount = false(scanC.nrV,t);
            H1NegCurvRatioCount = false(scanC.nrV,t);
            H2CurvRatioCount = false(scanC.nrV,t);
           % Local Curvature Diff
            H1PosCurvDiffCount = false(scanC.nrV,t);
            H1NegCurvDiffCount = false(scanC.nrV,t);
            H2CurvDiffCount = false(scanC.nrV,t);    
            n = size(A,1);
            tic;
            parfor i=1:t % permutation loop
                    ind = randperm(n);
                    Eperm = E(ind,:);  %#ok<PFBNS> Permute the errors in the reduced model
                    [~,~,~,~,Mperm,pctvarperm,~,statsperm] = plsregress(AR,Eperm,1);
                    R2perm = pctvarperm(2);
                    R2PartCount(i) = (R2perm>=R2Part);
                    [~,LocalSperm] = getLocalStatistic(Mperm,statsperm.Yresiduals,Eperm,'shape');
                    LocalSCount(:,i) = (LocalSperm >= LocalS);
                    % create morphs under permutation
                        V = range(1);
                        DX = m-V;
                        DY = DX*Mperm(2:end,:);
                        NY = Y-DY;
                        morphperm1 = clone(RefScan);
                        delete(morphperm1.PoseLM);
                        morphperm1.Vertices = scanC.Vertices + reshape(NY,3,length(NY)/3);
                        if ~isempty(RefScan.PoseLM)
                           transferPoseLM(morphperm1,RefScan);
                        end
                        V = range(2);
                        DX = m-V;
                        DY = DX*Mperm(2:end,:);
                        NY = Y-DY;
                        morphperm2 = clone(RefScan);
                        delete(morphperm2.PoseLM);
                        morphperm2.Vertices = scanC.Vertices + reshape(NY,3,length(NY)/3);
                        if ~isempty(RefScan.PoseLM)
                           transferPoseLM(morphperm2,RefScan);
                        end
                    % compare morphs under permutation
                      comparisonperm = compareMorphs(morphperm1,morphperm2);
                    % Distances along the normal
                      H1PosNormalDistCount(:,i) = (comparisonperm.NormalDistances >= STATS.comparison.NormalDistances);
                      H1NegNormalDistCount(:,i) = (comparisonperm.NormalDistances <= STATS.comparison.NormalDistances);
                      H2NormalDistCount(:,i) = (abs(comparisonperm.NormalDistances) >= abs(STATS.comparison.NormalDistances));
                    % Local Area ratio
                      H1PosAreaCount(:,i) = (comparisonperm.AreaRatios >= STATS.comparison.AreaRatios);
                      H1NegAreaCount(:,i) = (comparisonperm.AreaRatios <= STATS.comparison.AreaRatios);
                      H2AreaCount(:,i) = (abs(comparisonperm.AreaRatios) >= abs(STATS.comparison.AreaRatios));
                    % Curvature ratio
                      H1PosCurvRatioCount(:,i) = (comparisonperm.CurvatureRatios >= STATS.comparison.CurvatureRatios);
                      H1NegCurvRatioCount(:,i) = (comparisonperm.CurvatureRatios <= STATS.comparison.CurvatureRatios);
                      H2CurvRatioCount(:,i) = (abs(comparisonperm.CurvatureRatios) >= abs(STATS.comparison.CurvatureRatios));
                    % Curvature Diff
                      H1PosCurvDiffCount(:,i) = (comparisonperm.signedCurvatureDiff >= STATS.comparison.signedCurvatureDiff);
                      H1NegCurvDiffCount(:,i) = (comparisonperm.signedCurvatureDiff <= STATS.comparison.signedCurvatureDiff);
                      H2CurvDiffCount(:,i) = (abs(comparisonperm.signedCurvatureDiff) >= abs(STATS.comparison.signedCurvatureDiff));
            end % end permutation loop
            toc;
           % getting p values
            STATS.pR2Part = sum(R2PartCount)/t;
            STATS.pLocalS = sum(LocalSCount,2)/t;
            % Normal distances
            STATS.comparison.pH1PosNormalDistances = sum(H1PosNormalDistCount,2)/t;
            STATS.comparison.pH1NegNormalDistances = sum(H1NegNormalDistCount,2)/t;
            STATS.comparison.pH2NormalDistances = sum(H2NormalDistCount,2)/t;
            % Area Ratios
            STATS.comparison.pH1PosAreaRatios = sum(H1PosAreaCount,2)/t;
            STATS.comparison.pH1NegAreaRatios = sum(H1NegAreaCount,2)/t;
            STATS.comparison.pH2AreaRatios = sum(H2AreaCount,2)/t;
            % Curvature Ratios
            STATS.comparison.pH1PosCurvRatios = sum(H1PosCurvRatioCount,2)/t;
            STATS.comparison.pH1NegCurvRatios = sum(H1NegCurvRatioCount,2)/t;
            STATS.comparison.pH2CurvRatios = sum(H2CurvRatioCount,2)/t;
            % Curvature differences
            STATS.comparison.pH1PosCurvDiff = sum(H1PosCurvDiffCount,2)/t;
            STATS.comparison.pH1NegCurvDiff = sum(H1NegCurvDiffCount,2)/t;
            STATS.comparison.pH2CurvDiff = sum(H2CurvDiffCount,2)/t;
            out.STATS = STATS;
end

% function [LocalE,LocalS] = getLocalStatistic(M,R,B,type)
%          [n,nB] = size(B); 
%          P = B-R;
%          A = repmat(mean(B),n,1);
%          SST = sum((B-A).^2);
%          SSR = sum((P-A).^2);
%          LocalE = zeros(4,nB/3);
%          LocalE(1:3,:) = reshape(M(2,:),3,nB/3);
%          LocalE(4,:) = sqrt(sum(LocalE(1:3,:).^2));
%          SSR = reshape(SSR,3,nB/3);
%          SSR = sum(SSR);
%          SST = reshape(SST,3,nB/3);
%          SST = sum(SST);
%          LocalS = SSR./SST;
% end

%          % MorphCreation
%            %out.morphs = createMorphs(A,shape,RefScan,scale,CondValues);
%            [out.morphs,out.STATS] = createMorphsReducedModel(A,shape,RefScan,scale,CondValues);
%          % Compute image MorphAnalysis  
%            out.comparison = compareMorphs(out.morphs);
%          % compute effect and effect-size


%   ind = randperm(n);
%                     Eperm = E(ind,:);  %#ok<PFBNS> Permute the errors in the reduced model
%                     [~,~,~,~,Mperm,pctvarperm,~,statsperm] = plsregress(AR,Eperm,1);
%                     R2perm = pctvarperm(2);
%                     R2PartCount(i) = (R2perm>=R2Part);
%                     [~,LocalSperm] = getLocalStatistic(Mperm,statsperm.Yresiduals,Eperm,'shape');
%                     LocalSCount(:,i) = (LocalSperm >= LocalS);
%                     % create morphs under permutation
%                       permcounter = 0;
%                       for k=1:1:2
%                         V = range(k);
%                         permcounter = permcounter+1;
%                         DX = m-V;
%                         DY = DX*Mperm(2:end,:);
%                         NY = Y-DY;
%                         scan = clone(RefScan);
%                         delete(scan.PoseLM);
%                         scan.Vertices = scanC.Vertices + reshape(NY,3,length(NY)/3);
%                         eval(['morphperm' num2str(permcounter) ' = clone(scan);']);
%                         if ~isempty(RefScan.PoseLM)
%                             eval(['transferPoseLM(morphperm' num2str(permcounter) ',RefScan);']);
%                         end
%                       end
%                     % compare morphs under permutation
%                       comparisonperm = compareMorphs(morphperm1,morphperm2);
%                     % Distances along the normal
%                       H1PosNormalDistCount(:,i) = (comparisonperm.NormalDistances >= STATS.comparison.NormalDistances);
%                       H1NegNormalDistCount(:,i) = (comparisonperm.NormalDistances <= STATS.comparison.NormalDistances);
%                       H2NormalDistCount(:,i) = (abs(comparisonperm.NormalDistances) >= abs(STATS.comparison.NormalDistances));
