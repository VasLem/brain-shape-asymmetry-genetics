function out = ProcrustesAnova2WayAsymmetryMEM(X1,X2,nIter,factor,nSplits, showProgress)
%Factor is the value required to divide X1 and X2 to get the original
%value, it will be assumed to be 10000 (for shake of backwards
%compatibility) if not provided
%
%nSplits is supplied to split the computation into parts in case it is
%memory intensive
arguments
X1
X2
nIter=0
factor=10000;
nSplits=1;
showProgress=true;
end
[n,nrV,rep] = size(X1);
SSs = zeros(4,nrV);
Means = zeros(2,nrV);
pX1 = permute(X1, [3,2,1]);
pX2 = permute(X2, [3,2,1]);

if showProgress
    tic;
    ppb = ParforProgressbar(nrV);
end
parfor i=1:nrV
    Set1 = squeeze(single(pX1(:,i,:))/factor);
    Set2 = squeeze(single(pX2(:,i,:))/factor);    
    X = [Set1(:) Set2(:)];
    [~,TABLE,STATS] = anova2(X,rep,'off');
    ss = zeros(4,1);
    for j=1:4 
        ss(j) = TABLE{j+1,2};
    end
    SSs(:,i) =  ss(:);
    Means(:,i) = STATS.colmeans';
    if showProgress, ppb.increment(); end
end
if showProgress
    delete(ppb)
    toc;
end
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
if nIter==0, return; end
FFCount = false(nIter,nrV/3);
TFFCount = false(nIter,1);
IFCount = false(nIter,nrV/3);
TIFCount = false(nIter,1);
DFCount = false(nIter,nrV/3);
TDFCount = false(nIter,1);
%f = statusbar('Permuting');
pX1 = permute(pX1, [2,1,3]);
pX2 = permute(pX2, [2,1,3]);
% pX1 = permute(X1, [3,2,1]);
% pX2 = permute(X2, [3,2,1]);
% pX1 = pX1(:,:,1:n);
% pX2 = pX2(:,:,1:n);

splitSize = ceil(nrV / nSplits);
if showProgress
tic;
ppb = ParforProgressbar(nIter);    
end
parfor k=1:nIter
    asm = AsymmetryComponentsAnalysis(n,nrV, rep);
    SSI = zeros(4, nrV);
    SSD = zeros(4, nrV);
    SSF = zeros(4, nrV);
    for i=1:nSplits
        assignedInds =  ((i-1) * splitSize + 1): (min(i*splitSize,nrV));
        
        Set1 = single(pX1(assignedInds,:, :)) / factor;
        Set2 = single(pX2(assignedInds,:, :)) / factor;
        X = asm.shuffleColumnWise(Set1, Set2);
        SSD(:, assignedInds) = computeAnova2SS(X,rep);
        X = asm.shuffleRowWise(Set1, Set2);
        SSI(:, assignedInds)  = computeAnova2SS(X,rep);
        X = asm.shuffleResidual(Set1, Set2);
        SSF(:, assignedInds) = computeAnova2SS(X,rep);
    end
    
    asm = AsymmetryComponentsAnalysis(n,nrV, rep);
    % analyzing Direction effect
    SS = SSD;
    SS_D = SS(1,:);
    SS_F= SS(3,:);
    [~, ~, ~, ~, DF, TDF] = asm.directionEffect(SS_D, SS_F);
    DFCount(k,:) = DF>=LM.DF;
    TDFCount(k) = TDF>=Total.DF;
    
    % analyzing Individual effect / compared to the total error
    % now and not the interaction
    SS = SSI;
    SS_I = SS(2,:);
    SS_F= SS(3,:);
    [~, ~, ~, ~, IF, TIF] = asm.individualEffect(SS_I, SS_F);
    IFCount(k,:) = IF>=LM.IF;
    TIFCount(k) = TIF>=Total.IF;
    
    % analyzing interaction effect
    SS = SSF;
    SS_E = SS(4,:);
    SS_F= SS(3,:);
    [~, ~, ~, ~, FF, TFF] = asm.interactionEffect(SS_E, SS_F);
    FFCount(k,:) = FF>=LM.FF;
    TFFCount(k) = TFF>=Total.FF;
    if showProgress; ppb.increment(); end
end
if showProgress
    delete(ppb);
    toc;
end
out.LM.permFF = (sum(FFCount,1)+1)/(nIter+1);
out.Total.permFF = (sum(TFFCount)+1)/(nIter+1);
out.LM.permDF = (sum(DFCount,1)+1)/(nIter+1);
out.Total.permDF = (sum(TDFCount)+1)/(nIter+1);
out.LM.permIF = (sum(IFCount,1)+1)/(nIter+1);
out.Total.permIF = (sum(TIFCount)+1)/(nIter+1);
end