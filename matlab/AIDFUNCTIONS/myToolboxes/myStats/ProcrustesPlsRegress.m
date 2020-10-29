function [out] = ProcrustesPlsRegress(A,B,t,p,pu)         
         n = size(A,1);
         nA = size(A,2);
         nB = size(B,2);
         ncomp = min(n,nA);
         if nargin < 5, pu = zeros(1,nA); end% all unpaired data input
         if nargin < 4, p = (1:nA); end
         if nargin < 3, t = 0; end
         %if t< 1000, t = 0; end
         [~,~,~,~,M,pctvar,~,stats] = plsregress(A,B,ncomp);
         error = sum(stats.Yresiduals.^2);
         error = reshape(error,3,nB/3);
         error = sum(error)/(3*(n-1));
         Terror = sum(error);
         Differences = zeros(3,nB/3,nA);
         SS = zeros(nA,nB/3);
         Distances = zeros(nA,nB/3);
         PD = zeros(nA,1);
         RMSE = zeros(nA,1);
         for i=1:1:nA
             Differences(:,:,i) = reshape(M(1+i,:),3,nB/3);
             SS(i,:) = sum(Differences(:,:,i).^2);
             Distances(i,:) = sqrt(SS(i,:));
             PD(i) = sum(SS(i,:));
             RMSE(i) = sqrt((mean(Distances(i,:).^2)));
         end       
         Tpctvar = sum(pctvar(2,:));
         pctvarcount = false(2,ncomp,t);
         Tcount = false(1,t);
         DistancesCount = false(nB/3,nA,t);
         ErrorCount = false(nB/3,t);
         PDCount = false(t,nA);
         TerrorCount = false(t,1);
         if ~(t==0);
             tic
             parfor i=1:t % permutation loop               
                 Afor = A;
                 for k=1:1:length(p)
                     switch pu(k)
                         case 0 % unpaired randomization
                             ind=randperm(n);
                             Afor(:,p(k)) = Afor(ind,p(k));
                         case 1 % paired randomization
                             r = randi(2,n/2,1);
                             A1for = zeros(1,length(A(:,p(k)))/2);
                             A2for = zeros(1,length(A(:,p(k)))/2);
                             for l=1:1:n/2
                                switch r(l)
                                    case 1
                                        A1for(l) = A(l,p(k));
                                        A2for(l) = A(l+n/2,p(k));
                                    case 2
                                        A1for(l) = A(l+n/2,p(k));
                                        A2for(l) = A(l,p(k));
                                end
                             end
                            Afor(:,p(k)) = [A1for A2for]';
                     end
                 end
                 [~,~,~,~,Mfor,pctvarfor,~,stats] = plsregress(Afor,B,ncomp);
                 errorfor = sum(stats.Yresiduals.^2);
                 errorfor = reshape(errorfor,3,nB/3);
                 errorfor = sum(errorfor)/(3*(n-1));
                 Terrorfor = sum(errorfor);
                 pctvarcount(:,:,i) = pctvarfor>=pctvar;
                 Tfor = sum(pctvarfor(2,:));
                 Tcount(i) = Tfor>=Tpctvar;
                 Distancesfor = zeros(size(Distances));
                 PDfor = zeros(size(PD));
                 for k=1:1:nA
                    Differencesfor = reshape(Mfor(1+k,:),3,nB/3);
                    Distancesfor(k,:) = sqrt(sum(Differencesfor.^2));
                    PDfor(k) = sum(sum(Differencesfor.^2));
                 end
                 DistancesCount(:,:,i) = (Distancesfor >= Distances)';
                 ErrorCount(:,i) = (errorfor <= error)';
                 PDCount(i,:) = PDfor >= PD;
                 TerrorCount(i) = Terrorfor <= Terror;
             end % end permutation loop
             toc
         end
         pctvarP = (sum(pctvarcount,3)/t);
         TpctvarP = (sum(Tcount)/t);    
         out.M = M;
         out.pctvar = pctvar*100;
         out.pctvarP = pctvarP;
         out.Tpctvar = Tpctvar*100;
         out.TpctvarP = TpctvarP;
         out.Differences = Differences;
         out.SS = SS;
         out.PD = PD;
         out.RMSE = RMSE;
         out.pPD = (sum(PDCount,1)/t);
         out.Error = error;
         out.pError = (sum(ErrorCount,2)/t)';
         out.Terror = Terror;
         out.pTerror = sum(TerrorCount)/t;
         out.Distances = Distances;
         out.pDistances = sum(DistancesCount,3)/t;
         out.p1Significant = out.pDistances <= 0.1;
         out.p1Perc = (sum(out.p1Significant,1)/(nB/3))*100;
         out.p05Significant = out.pDistances <= 0.05;
         out.p05Perc = (sum(out.p05Significant,1)/(nB/3))*100;
         out.p001Significant = out.pDistances <= 0.001;
         out.p001Perc = (sum(out.p001Significant,1)/(nB/3))*100;         
end

% function [out] = ProcrustesPlsRegress(A,B,t,p,pu)         
%          n = size(A,1);
%          nA = size(A,2);
%          nB = size(B,2);
%          ncomp = min(n,nA);
%          if nargin < 5, pu = zeros(1,nA); end% all unpaired data input
%          if nargin < 4, p = (1:nA); end
%          if nargin < 3, t = 0; end
%          %if t< 1000, t = 0; end
%          [~,~,~,~,M,pctvar,~,stats] = plsregress(A,B,ncomp);
%          error = sum(stats.Yresiduals.^2);
%          error = reshape(error,3,nB/3);
%          error = sum(error)/(3*(n-1));
%          Differences = zeros(3,nB/3,nA);
%          SS = zeros(nA,nB/3);
%          Distances = zeros(nA,nB/3);
%          PD = zeros(nA,1);
%          RMSE = zeros(nA,1);
%          for i=1:1:nA
%              Differences(:,:,i) = reshape(M(1+i,:),3,nB/3);
%              SS(i,:) = sum(Differences(:,:,i).^2);
%              Distances(i,:) = sqrt(SS(i,:));
%              PD(i) = sum(SS(i,:));
%              RMSE(i) = sqrt((mean(Distances(i,:).^2)));
%          end       
%          Tpctvar = sum(pctvar(2,:));
%          pctvarcount = false(2,ncomp,t);
%          Tcount = false(1,t);
%          DistancesCount = false(nB/3,nA,t);
%          PDCount = false(t,nA);
%          if ~(t==0);
%              tic
%              parfor i=1:t % permutation loop               
%                  Afor = A;
%                  for k=1:1:length(p)
%                      switch pu(k)
%                          case 0 % unpaired randomization
%                              ind=randperm(n);
%                              Afor(:,p(k)) = Afor(ind,p(k));
%                          case 1 % paired randomization
%                              r = randi(2,n/2,1);
%                              A1for = zeros(1,length(A(:,p(k)))/2);
%                              A2for = zeros(1,length(A(:,p(k)))/2);
%                              for l=1:1:n/2
%                                 switch r(l)
%                                     case 1
%                                         A1for(l) = A(l,p(k));
%                                         A2for(l) = A(l+n/2,p(k));
%                                     case 2
%                                         A1for(l) = A(l+n/2,p(k));
%                                         A2for(l) = A(l,p(k));
%                                 end
%                              end
%                             Afor(:,p(k)) = [A1for A2for]';
%                      end
%                  end
%                  [~,~,~,~,Mfor,pctvarfor,~,stats] = plsregress(Afor,B,ncomp);
%                  pctvarcount(:,:,i) = pctvarfor>=pctvar;
%                  Tfor = sum(pctvarfor(2,:));
%                  Tcount(i) = Tfor>=Tpctvar;
%                  Distancesfor = zeros(size(Distances));
%                  PDfor = zeros(size(PD));
%                  for k=1:1:nA
%                     Differencesfor = reshape(Mfor(1+k,:),3,nB/3);
%                     Distancesfor(k,:) = sqrt(sum(Differencesfor.^2));
%                     PDfor(k) = sum(sum(Differencesfor.^2));
%                  end
%                  DistancesCount(:,:,i) = (Distancesfor >= Distances)';
%                  PDCount(i,:) = PDfor >= PD;
%              end % end permutation loop
%              toc
%          end
%          pctvarP = (sum(pctvarcount,3)/t);
%          TpctvarP = (sum(Tcount)/t);    
%          out.M = M;
%          out.pctvar = pctvar*100;
%          out.pctvarP = pctvarP;
%          out.Tpctvar = Tpctvar*100;
%          out.TpctvarP = TpctvarP;
%          out.Differences = Differences;
%          out.SS = SS;
%          out.PD = PD;
%          out.RMSE = RMSE;
%          out.pPD = (sum(PDCount,1)/t);
%          out.Error = error;
%          out.Distances = Distances;
%          out.pDistances = sum(DistancesCount,3)/t;
%          out.p1Significant = out.pDistances <= 0.1;
%          out.p1Perc = (sum(out.p1Significant,1)/(nB/3))*100;
%          out.p05Significant = out.pDistances <= 0.05;
%          out.p05Perc = (sum(out.p05Significant,1)/(nB/3))*100;
%          out.p001Significant = out.pDistances <= 0.001;
%          out.p001Perc = (sum(out.p001Significant,1)/(nB/3))*100;         
% end