function out = ProcrustesAnova2WayAsymmetryMEM(X1,X2,t)
         if nargin < 3, t = 0; end
         [n,nrV,rep] = size(X1);
         SS = zeros(4,nrV);
         Means = zeros(2,nrV);
         factor = 10000;
         tic;
         [path,ID] = setupParForProgress(nrV);
         parfor i=1:nrV
            Set1 = squeeze(single(X1(:,i,:))/factor)'; 
            %Set1 = reshape(Set1,n,rep)'; 
            Set2 = squeeze(single(X2(:,i,:))/factor)'; 
            %Set2 = reshape(Set2,n,rep)';
            X = [Set1(:) Set2(:)];
            [~,TABLE,STATS] = anova2(X,rep,'off');
            ss = zeros(4,1);
            for j=1:4
                ss(j) = TABLE{j+1,2};
            end
            SS(:,i) =  ss(:);
            Means(:,i) = STATS.colmeans';
            parfor_progress;
         end
         closeParForProgress(path,ID)
         toc;
%         pSet1 = reshape(X1, [], 1);
%         pSet2= reshape(X2, [], 1);
%         pX = [pSet1(:) pSet2(:)];
%          [~, PROP_TABLE, PROP_STATS] = anova2(pX, t);
        
         % Return stats to 3d setup
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
         LM.E = LM.E./((rep-1)*n*2*3);
         Total.E = Total.E/((rep-1)*n*2*(nrV-7));
         % FLuctuating
         LM.F = sum(SS_F);
         Total.F = sum(LM.F);
         LM.F = LM.F./(3*(n-1));
         Total.F = Total.F/((n-1)*(nrV-7));
         LM.FF = LM.F./LM.E;
         LM.FP = Ftest(LM.FF,3*(n-1),((rep-1)*n*2*3));
         Total.FF =  Total.F/Total.E;
         Total.FP = Ftest(Total.FF,(n-1)*(nrV-7),(rep-1)*n*2*(nrV-7));
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
         %%%% To be verified in R
         LM.D = sum(SS_D);
         Total.D = sum(LM.D);
         LM.D = (LM.D./3)/n;
         Total.D = (Total.D/(nrV-7))/n;
         LM.DF = LM.D./LM.F;
         LM.DP = Ftest(LM.DF,3,3*(n-1));
         Total.DF =  Total.D/Total.F;
         Total.DP = Ftest(Total.DF,nrV-7,(n-1)*(nrV-7));
         %%%%
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
         %f = statusbar('Permuting');
         [path,ID] = setupParForProgress(t);
         parfor k=1:1:t
             %k=1;
             %disp(num2str(k));
             SSF = zeros(4,nrV);
             SSI = zeros(4,nrV);
             SSD = zeros(4,nrV);
             for i=1:nrV
                %i=1;
                Set1 = single(X1(:,i,:))/factor;  
                Set2 = single(X2(:,i,:))/factor;
                % Colom-wise shuffeling of cells for Directional effect
                    Set1Copy = Set1;
                    Set2Copy = Set2;
                    r = randi(2,n,1);
                    index = find(r==2);
                    Set1Copy(index,:,:) = Set2(index,:,:);
                    Set2Copy(index,:,:) = Set1(index,:,:);
                    Set1Copy = reshape(Set1Copy,n,rep)'; 
                    Set2Copy = reshape(Set2Copy,n,rep)';
                    X = [Set1Copy(:) Set2Copy(:)];
                    [~,TABLE] = anova2(X,rep,'off');
                    ss = zeros(4,1);
                    for j=1:4
                        ss(j) = TABLE{j+1,2};
                    end
                    SSD(:,i) =  ss;
                % Row-wise shuffeling for Individual effect
                    Set1Copy = reshape(Set1,n,rep)'; 
                    Set2Copy = reshape(Set2,n,rep)';
                    index = randperm(n);
                    Set1Copy = Set1Copy(:,index);
                    X = [Set1Copy(:) Set2Copy(:)];
                    [~,TABLE] = anova2(X,rep,'off');
                    ss = zeros(4,1);
                    for j=1:4
                        ss(j) = TABLE{j+1,2};
                    end
                    SSI(:,i) =  ss(:);
                % Residual shuffeling for Interaction effect
                    X = [Set1Copy(:) Set2Copy(:)];
                    avgC = mean(X,1);
                    avgR = mean(X,2);
                    avg = mean(X(:));
                    X = X - repmat(avgC,n*rep,1) - repmat(avgR,1,2) + repmat(avg,n*rep,2);
                    index = randperm(n*rep*2);
                    X = reshape(X(index),n*rep,2);
                    [~,TABLE] = anova2(X,rep,'off');
                    ss = zeros(4,1);
                    for j=1:4
                        ss(j) = TABLE{j+1,2};
                    end
                    SSF(:,i) =  ss(:);
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
                 D = (D./3)/n;TD = (TD/(nrV-7))/n;
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
                 E = E./((rep-1)*n*2*3);
                 TE = TE/((rep-1)*n*2*(nrV-7));
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
                 E = E./((rep-1)*n*2*3);
                 TE = TE/((rep-1)*n*2*(nrV-7));
                 % getting Fluctuating MS
                 F = sum(SS_F); TF = sum(F);
                 F = F./(3*(n-1)); TF = TF/((n-1)*(nrV-7));
                 % getting F-statistic
                 FF = F./E;
                 TFF = TF/TE;
                 FFCount(k,:) = (FF>=LM.FF);
                 TFFCount(k) = TFF>=Total.FF;
                 parfor_progress;
         end
         closeParForProgress(path,ID);
         out.LM.permFF = (sum(FFCount,1)+1)/(t+1);
         out.Total.permFF = (sum(TFFCount)+1)/(t+1);
         out.LM.permDF = (sum(DFCount,1)+1)/(t+1);
         out.Total.permDF = (sum(TDFCount)+1)/(t+1);
         out.LM.permIF = (sum(IFCount,1)+1)/(t+1);
         out.Total.permIF = (sum(TIFCount)+1)/(t+1);      
end

