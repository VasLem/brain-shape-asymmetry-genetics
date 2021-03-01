function out = ProcrustesAnova2WayAsymmetryMEM(X1,X2,t)
if nargin < 3, t = 0; end
[n,nrV,rep] = size(X1);
SSs = zeros(4,nrV);
Means = zeros(2,nrV);


factor = 10000;
tic;
[path,ID] = setupParForProgress(nrV);
parfor i=1:nrV
    Set1 = squeeze(single(X1(:,i,:))/factor)';
    Set2 = squeeze(single(X2(:,i,:))/factor)';
    X = [Set1(:) Set2(:)];
    [~,TABLE,STATS] = anova2(X,rep,'off');
    ss = zeros(4,1);
    for j=1:4
        ss(j) = TABLE{j+1,2};
    end
    SSs(:,i) =  ss(:);
    Means(:,i) = STATS.colmeans';
    parfor_progress;
end
closeParForProgress(path,ID)
toc;

MeanX1 = reshape(Means(1,:),3,length(Means(1,:))/3);
MeanX2 = reshape(Means(2,:),3,length(Means(2,:))/3);
MeanDiff = MeanX2 - MeanX1;
asm = AsymmetryComponentsAnalysis(n,nrV, rep);

[LM, Total] = asm.apply(SSs);

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
% for k=1:1:t
    %k=1;
    %disp(num2str(k));
    SSF = zeros(4,nrV);
    SSI = zeros(4,nrV);
    SSD = zeros(4,nrV);
    for i=1:nrV
        %i=1;
        Set1 = squeeze(single(X1(:,i,:)))/factor;
        Set2 = squeeze(single(X2(:,i,:)))/factor;
        % Column-wise shuffling of cells for Directional effect
        Set1Copy = Set1';
        Set2Copy = Set2';
        r = randi(2,n,1);
        index = find(r==2);
        Set1Copy(:, index) = Set2Copy(:, index);
        Set2Copy(:, index) = Set1Copy(:, index);
        X = [Set1Copy(:) Set2Copy(:)];
        [~,TABLE] = anova2(X,rep,'off');
        ss = zeros(4,1);
        for j=1:4
            ss(j) = TABLE{j+1,2};
        end
        SSD(:,i) =  ss;
        % Row-wise shuffling for Individual effect
        Set1Copy = Set1';
        Set2Copy = Set2';
        index = randperm(n);
        Set1Copy = Set1Copy(:,index);
        X = [Set1Copy(:) Set2Copy(:)];
        [~,TABLE] = anova2(X,rep,'off');
        ss = zeros(4,1);
        for j=1:4
            ss(j) = TABLE{j+1,2};
        end
        SSI(:,i) =  ss(:);
        % Residual shuffling for Interaction effect
        Set1Copy = Set1';
        Set2Copy = Set2';
        X = [Set1Copy(:) Set2Copy(:)];
        avgC = mean(X,1);
        avgR = mean(X,2);
        avg = mean(X(:));
        X = X - avgC - avgR + avg;
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
    SS_F= SS(3,:);
    SS_D = SS(1,:);
    asm = AsymmetryComponentsAnalysis(n,nrV, rep);
    [~, ~, ~, ~, DF, TDF] = asm.directionEffect(SS_D, SS_F);
    DFCount(k,:) = (DF>=LM.DF);
    TDFCount(k) = TDF>=Total.DF;
    
    
    % analyzing Individual effect / compared to the total error
    % now and not the interaction
    SS = SSI;
    SS_I = SS(2,:);
    SS_F = SS(3,:);
    [~, ~, ~, ~, IF, TIF] = asm.individualEffect(SS_I, SS_F);
    IFCount(k,:) = (IF>=LM.IF);
    TIFCount(k) = TIF>=Total.IF;
    
    % analyzing interaction effect
    SS = SSF;
    SS_E = SS(4,:);
    SS_F = SS(3,:);
    [~, ~, ~, ~, FF, TFF] = asm.interactionEffect(SS_E, SS_F);
    FFCount(k,:) = (FF>=LM.FF);
    TFFCount(k) = TFF>=Total.FF;
    parfor_progress;
end
closeParForProgress(  path,ID);
out.LM.permFF = (sum(FFCount,1)+1)/(t+1);
out.Total.permFF = (sum(TFFCount)+1)/(t+1);
out.LM.permDF = (sum(DFCount,1)+1)/(t+1);
out.Total.permDF = (sum(TDFCount)+1)/(t+1);
out.LM.permIF = (sum(IFCount,1)+1)/(t+1);
out.Total.permIF = (sum(TIFCount)+1)/(t+1);


end

