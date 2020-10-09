function [out] = sampledPlsRegress(A,B,t,p,type)         
         n = size(A,1);
         nA = size(A,2);
         nB = size(B,2);
         ncomp = min(n,nA);
         if naring <5 type = 'Default'; end
         if nargin < 4, p = (1:nA); end
         if nargin < 3, t = 0; end
         [~,~,~,~,M,pctvar,~,stats] = plsregress(A,B,ncomp);
         Differences = zeros(3,nB/3,nA);
         SS = zeros(nA,nB/3);
         MS = zeros(nA,nB/3);
         Distances = zeros(nA,nB/3);
         PD = zeros(nA,1);
         RMSE = zeros(nA,1);
         for i=1:1:nA
             Differences(:,:,i) = reshape(M(1+i,:),3,nB/3);
             SS(i,:) = sum(Differences(:,:,i).^2);
             Distances(i,:) = sqrt(SS(i,:));
             MS(i,:) = SS(i,:)/3;
             PD(i) = sum(SS(i,:));
             RMSE(i) = sqrt((mean(Distances(i,:).^2)));
         end       
         out.M = M;
         out.pctvar = pctvar;
         out.Tpctvar = sum(out.pctvar,2);
         out.Differences = Differences;
         out.SS = SS;
         out.MS = MS;
         out.PD = PD;
         out.RMSE = RMSE;
         out.Stats = stats;
         out.Distances = Distances;
         if t==0, return; end
         f = statusbar('Permuting');
         drawnow;
         for k=1:1:length(p);
             DistancesCount = false(nB/3,t);
             Dist = Distances(p(k),:);
             PDCount = false(t,1);
             pd = PD(p(k));
             % creating reduced model
             index = setdiff((1:nA),p(k));
             AR = A(:,index);
             [~,~,~,~,~,pctvar,~,statsR] = plsregress(AR,B,ncomp-1);
             E = statsR.Yresiduals;
             Afor = A(:,p(k));
             parfor i=1:t % permutation loop
                    ind = randperm(n);
                    Bfor = E(ind,:);  %#ok<PFBNS>
                    [~,~,~,~,Mfor,~,~,~] = plsregress(Afor,Bfor,1);
                    Differencesfor = reshape(Mfor(2,:),3,nB/3);
                    Distancesfor = sqrt(sum(Differencesfor.^2));
                    PDfor = sum(sum(Differencesfor.^2));
                    DistancesCount(:,i) = (Distancesfor >= Dist);
                    PDCount(i) = PDfor >= pd;
              end % end permutation loop
              PermTest(p(k)).Tpctvar = sum(pctvar,2);
              PermTest(p(k)).pPD = (sum(PDCount)/t); %#ok<*AGROW>
              PermTest(p(k)).pDistances = sum(DistancesCount,2)/t;
              PermTest(p(k)).p1Significant = PermTest(p(k)).pDistances <= 0.1;
              PermTest(p(k)).p1Perc = (sum(PermTest(p(k)).p1Significant,1)/(nB/3))*100;
              PermTest(p(k)).p05Significant = PermTest(p(k)).pDistances <= 0.05;
              PermTest(p(k)).p05Perc = (sum(PermTest(p(k)).p05Significant,1)/(nB/3))*100;
              PermTest(p(k)).p001Significant = PermTest(p(k)).pDistances <= 0.001;
              PermTest(p(k)).p001Perc = (sum(PermTest(p(k)).p001Significant,1)/(nB/3))*100;
              statusbar(k/length(p),f);
              drawnow;
         end 
         delete(f);
         out.PermTest = PermTest;
end



% function [out] = ProcrustesPlsRegressv3(A,B,t,p)         
%          n = size(A,1);
%          nA = size(A,2);
%          nB = size(B,2);
%          ncomp = min(n,nA);
%          if nargin < 4, p = (1:nA); end
%          if nargin < 3, t = 0; end
%          [~,~,~,~,M,pctvar,~,stats] = plsregress(A,B,ncomp);
%          Differences = zeros(3,nB/3,nA);
%          SS = zeros(nA,nB/3);
%          Distances = zeros(nA,nB/3);
%          PD = zeros(nA,1);
%          RMSE = zeros(nA,1);
%          for i=1:1:nA
%              Differences(:,:,i) = reshape(M(1+i,:),3,nB/3);
%              SS(i,:) = sum(Differences(:,:,i).^2);
%              Distances(i,:) = sqrt(SS(i,:));
%              %Distances(i,:) = (SS(i,:))/3;
%              PD(i) = sum(SS(i,:));
%              RMSE(i) = sqrt((mean(Distances(i,:).^2)));
%          end       
%          DistancesCount = false(nB/3,nA,t);
%          PDCount = false(t,nA);
%          % creating reduced model
%          index = setdiff((1:p),p);
%          AR = A(:,index);
%          [~,~,~,~,~,~,~,statsR] = plsregress(AR,B,ncomp-1);
%          E = statsR.Yresiduals;
%          Afor = A(:,p);
%          if ~(t==0);
%              tic
%              parfor i=1:t % permutation loop
%                  ind = randperm(n);
%                  Bfor = E(ind,:); 
%                  [~,~,~,~,Mfor,~,~,~] = plsregress(Afor,Bfor,1);
%                  Distancesfor = zeros(size(Distances));
%                  PDfor = zeros(size(PD));
%                  k=p;
%                     Differencesfor = reshape(Mfor(2,:),3,nB/3);
%                     Distancesfor(k,:) = sqrt(sum(Differencesfor.^2));
%                     %Distancesfor(k,:) = (sum(Differencesfor.^2))/3;
%                     PDfor(k) = sum(sum(Differencesfor.^2));
%                  DistancesCount(:,:,i) = (Distancesfor >= Distances)';
%                  PDCount(i,:) = PDfor >= PD;
%              end % end permutation loop
%              toc
%          end 
%          out.M = M;
%          out.pctvar = pctvar;
%          out.Differences = Differences;
%          out.SS = SS;
%          out.PD = PD;
%          out.RMSE = RMSE;
%          out.Stats = stats;
%          out.pPD = (sum(PDCount,1)/t);
%          out.Distances = Distances;
%          out.pDistances = sum(DistancesCount,3)/t;
%          out.p1Significant = out.pDistances <= 0.1;
%          out.p1Perc = (sum(out.p1Significant,1)/(nB/3))*100;
%          out.p05Significant = out.pDistances <= 0.05;
%          out.p05Perc = (sum(out.p05Significant,1)/(nB/3))*100;
%          out.p001Significant = out.pDistances <= 0.001;
%          out.p001Perc = (sum(out.p001Significant,1)/(nB/3))*100;         
% end
% 
