function [out] = ProcrustesMVDispersion(X1,X2,t)
% This function tests the equality of group dispersion
% X is a DataMatrix 
if nargin < 3, t = 100; end
tic;
nX1 = size(X1,1);
nX2 = size(X2,1);
AvgX1 = mean(X1);
AvgX2 = mean(X2);
% Take the within group Distances
PDX1 = sum((X1-repmat(AvgX1,nX1,1)).^2,2);
PDX2 = sum((X2-repmat(AvgX2,nX2,1)).^2,2);
out.ProcrDispX1 = mean(PDX1);
out.ProcrDispX2 = mean(PDX2);
X = [PDX1(:);PDX2(:)];
G = [ones(nX1,1);zeros(nX2,1)];
[P,TAB] = ANOVA1(X,G,'off');
F = TAB{2,5};
out.F = F;
out.P = P;
% Fstatistic
FCount = false(1,t);
nT = nX1+nX2;
disp('Permuting');
parfor i=1:t
    ind = randperm(nT);
    [~,TABfor] = anova1(X(ind),G,'off');
    FCount(i) = TABfor{2,5} >= F;
end
disp('Done');
toc;
out.pF = (sum(FCount)/t);
end