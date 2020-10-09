function out = NonParametricManova2(X1,X2,NL,NLsets,t)
         if nargin < 5, t = 0; end
         if nargin < 4, NLsets = 2; end
         if nargin < 3, NL = 0; end
         [nX1,nV] = size(X1);
         % Noise injection
         X1new = zeros(NLsets*nX1,nV);
         X2new = zeros(NLsets*nX1,nV);
         ind = (1:NLsets*nX1);
         for i=1:1:NLsets
             tmp = find(mod(ind,NLsets)==i-1);
             if i==1
                X1new(tmp,:) = X1;
                X2new(tmp,:) = X2;
             else
                X1new(tmp,:) = X1+ NL*randn(size(X1));
                X2new(tmp,:) = X2+ NL*randn(size(X2));
             end
         end
         % Building the Distance Matrix
         D = squareform(pdist([X1new;X2new],'euclidean'));
         tic;
         %Following the nomenclature from M.J. Anderson, A new method for
         %non-parametric multivariate analysis of variance
         a = 2;
         b = nX1;
         n = NLsets;
         % 
         % Creating Distance Matrix;
         parfor i=1:nrV
            Set1 = X1(:,i,:); 
            Set1 = reshape(Set1,n,NLsets)'; 
            Set2 = X2(:,i,:); 
            Set2 = reshape(Set2,n,NLsets)';
            X = [Set1(:) Set2(:)];
            [~,TABLE,STATS] = anova2(X,NLsets,'off');
            SS(:,i) =  cell2mat(TABLE(2:5,2));
            Means(:,i) = STATS.colmeans';
         end
         toc;
         SS_D = reshape(SS(1,:),3,length(SS(1,:))/3);
         SS_I = reshape(SS(2,:),3,length(SS(2,:))/3);
         SS_F = reshape(SS(3,:),3,length(SS(3,:))/3);
         SS_E = reshape(SS(4,:),3,length(SS(4,:))/3);
         
         MeanX1 = reshape(Means(1,:),3,length(Means(1,:))/3);
         MeanX2 = reshape(Means(2,:),3,length(Means(2,:))/3);
         MeanDiff = MeanX2 - MeanX1;
         % Error
         LM.E = sum(SS_E);
         Total.E = sum(LM.E);
         LM.E = LM.E./((NLsets-1)*n*2*3);
         Total.E = Total.E/((NLsets-1)*n*2*(nrV-7));
         % FLuctuating
         LM.F = sum(SS_F);
         Total.F = sum(LM.F);
         LM.F = LM.F./(3*(n-1));
         Total.F = Total.F/((n-1)*(nrV-7));
         LM.FF = LM.F./LM.E;
         LM.FP = Ftest(LM.FF,3*(n-1),((NLsets-1)*n*2*3));
         Total.FF =  Total.F/Total.E;
         Total.FP = Ftest(Total.FF,(n-1)*(nrV-7),(NLsets-1)*n*2*(nrV-7));
         % Individuals
         LM.I = sum(SS_I);
         Total.I = sum(LM.I);
         LM.I = LM.I./(3*(n-1));
         Total.I = Total.I/((n-1)*(nrV-7));
         LM.IF = LM.I./LM.F;
         LM.IP = Ftest(LM.IF,(3*(n-1)),(3*(n-1)));
         Total.IF =  Total.I/Total.F;
         Total.IP = Ftest(Total.IF,((n-1)*(nrV-7)),(n-1)*(nrV-7));         
         % Directional
         LM.D = sum(SS_D);
         Total.D = sum(LM.D);
         LM.D = LM.D./3;
         Total.D = Total.D/(nrV-7);
         LM.DF = LM.D./LM.F;
         LM.DP = Ftest(LM.DF,3,3*(n-1));
         Total.DF =  Total.D/Total.F;
         Total.DP = Ftest(Total.DF,nrV-7,(n-1)*(nrV-7));
         out.Total = Total;
         out.LM = LM;
         out.MeanX1 = MeanX1;
         out.MeanX2 = MeanX2;
         out.MeanDiff = MeanDiff;         
         % Disp, permuting test
         if t==0, return; end
         FFCount = false(t,nrV/3);
         TFFCount = false(t,1);
         IFCount = false(t,nrV/3);
         TIFCount = false(t,1);
         DFCount = false(t,nrV/3);
         TDFCount = false(t,1);
         f = statusbar('Permuting');
         drawnow;
         for k=1:1:t
             SSF = zeros(4,nrV);
             SSI = zeros(4,nrV);
             SSD = zeros(4,nrV);
             parfor i=1:nrV
                Set1 = X1(:,i,:); 
                Set2 = X2(:,i,:);
                % Colom-wise shuffeling of cells for Directional effect
                    Set1Copy = Set1;
                    Set2Copy = Set2;
                    r = randi(2,n,1);
                    index = find(r==2);
                    Set1Copy(index,:,:) = Set2(index,:,:);
                    Set2Copy(index,:,:) = Set1(index,:,:);
                    Set1Copy = reshape(Set1Copy,n,NLsets)'; 
                    Set2Copy = reshape(Set2Copy,n,NLsets)';
                    X = [Set1Copy(:) Set2Copy(:)];
                    [~,TABLE] = anova2(X,NLsets,'off');
                    SSD(:,i) =  cell2mat(TABLE(2:5,2));
                % Row-wise shuffeling for Individual effect
                    Set1Copy = reshape(Set1,n,NLsets)'; 
                    Set2Copy = reshape(Set2,n,NLsets)';
                    index = randperm(n);
                    Set1Copy = Set1Copy(:,index);
                    X = [Set1Copy(:) Set2Copy(:)];
                    [~,TABLE] = anova2(X,NLsets,'off');
                    SSI(:,i) =  cell2mat(TABLE(2:5,2));
                % Residual shuffeling for Interaction effect
                    X = [Set1Copy(:) Set2Copy(:)];
                    avgC = mean(X,1);
                    avgR = mean(X,2);
                    avg = mean(X(:));
                    X = X - repmat(avgC,n*NLsets,1) - repmat(avgR,1,2) + repmat(avg,n*NLsets,2);
                    index = randperm(n*NLsets*2);
                    X = reshape(X(index),n*NLsets,2);
                    [~,TABLE] = anova2(X,NLsets,'off');
                    SSF(:,i) =  cell2mat(TABLE(2:5,2));
             end
             % analyzing Direction effect
                 SS = SSD;
                 SS_D = reshape(SS(1,:),3,length(SS(1,:))/3);
                 SS_F = reshape(SS(3,:),3,length(SS(3,:))/3);
                 % Getting Fluctuating MS as error term
                 F = sum(SS_F); TF = sum(F);
                 F = F./(3*(n-1)); TF = TF/((n-1)*(nrV-7));
                 % Getting Direction MS
                 D = sum(SS_D); TD = sum(D);
                 D = D./3;TD = TD/(nrV-7);
                 % Getting F statistic
                 DF = D./F;
                 TDF = TD/TF;
                 DFCount(k,:) = (DF>=LM.DF);
                 TDFCount(k) = TDF>=Total.DF;
             % analyzing Individual effect / compared to the total error
             % now and not the interaction
                 SS = SSI;
                 SS_I = reshape(SS(2,:),3,length(SS(2,:))/3);
                 SS_F = reshape(SS(3,:),3,length(SS(3,:))/3);
                 SS_E = reshape(SS(4,:),3,length(SS(4,:))/3);
                 E = sum(SS_E);
                 TE = sum(E);
                 E = E./((NLsets-1)*n*2*3);
                 TE = TE/((NLsets-1)*n*2*(nrV-7));
                 % Getting Fluctuating MS as error term
                 F = sum(SS_F); TF = sum(F);
                 F = F./(3*(n-1)); TF = TF/((n-1)*(nrV-7));
                 % Getting Direction MS
                 I = sum(SS_I); TI = sum(I);
                 I = I./(3*(n-1));TI = TI/((nrV-7)*(n-1));
                 % Getting F statistic
                 IF = I./F;
                 TIF = TI/TF;
                 IFCount(k,:) = (IF>=LM.IF);
                 TIFCount(k) = TIF>=Total.IF;
            % analyzing interaction effect
                 SS = SSF;
                 SS_F = reshape(SS(3,:),3,length(SS(3,:))/3);
                 SS_E = reshape(SS(4,:),3,length(SS(4,:))/3);
                 % getting error MS as error term
                 E = sum(SS_E);
                 TE = sum(E);
                 E = E./((NLsets-1)*n*2*3);
                 TE = TE/((NLsets-1)*n*2*(nrV-7));
                 % getting Fluctuating MS
                 F = sum(SS_F); TF = sum(F);
                 F = F./(3*(n-1)); TF = TF/((n-1)*(nrV-7));
                 % getting F-statistic
                 FF = F./E;
                 TFF = TF/TE;
                 FFCount(k,:) = (FF>=LM.FF);
                 TFFCount(k) = TFF>=Total.FF;
                 statusbar(k/t,f);
                 drawnow;
         end
         delete(f);
         out.LM.permFF = sum(FFCount,1)/t;
         out.Total.permFF = sum(TFFCount)/t;
         out.LM.permDF = sum(DFCount,1)/t;
         out.Total.permDF = sum(TDFCount)/t;
         out.LM.permIF = sum(IFCount,1)/t;
         out.Total.permIF = sum(TIFCount)/t;      
end

% function [MSA,MSB,MSAB,MSR] = extractSS(D,a,b,n)
%          N = a*b*n;
%          T = 0;A = 0;B = 0;R = 0;
%          for i=1:1:(N-1)
%              for j=(i+1):1:N
%                  d = D(i,j)^2;
%                  T = T+d;
%                  if inA(i,N,a)==inA(j,N,a)
%                     A = A+d;
%                  end
%                  if inB(i,N,a,b)==inB(i,N,a,b)
%                     B = B+d;
%                     if inA(i,N,a)==inA(i,N,a)
%                        R = R+d;
%                     end
%                  end
%              end
%          end
%          SST = T/N;
%          SSA = A/(n*b);
%          SSB = B/(n*a);
%          SSR = R/n;
%          SSAB = SST-SSA-SSB-SSR;
%          MSA = SSA/(a-1);
%          MSB = SSB/(b-1);
%          MSAB = SSAB/((a-1)*(b-1));
%          MSR = SSR/(N-a*b);     
% end
% function out = inA(i,N,a)
%          nA = N/a;
%          out = ceil(i/nA);
% end
% function out = inB(i,N,a,b)
%          nA = N/a;
%          nB = (N/a)/b;
%          i = i-((inA(i,N,a)-1)*nA);
%          out = ceil(i/nB);
% end
