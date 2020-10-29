function [out] = ShapeModelMVDispersion(X1,X2,AM,percvar,type,t)
% This function tests the equality of group dispersion
% X is a DataMatrix
if nargin < 3, t = 100; end
tic;
nX1 = size(X1,1);
nX2 = size(X2,1);
% Building Shape model
disp('Building Shape Model');
SM = shapePCA;
SM.RefScan = AM;
getAverage(SM,[X1',X2']);
getModel(SM,[X1',X2']);
stripPercVar(SM,percvar);
CX1 = SM.Tcoeff(1:nX1,:);
CX2 = SM.Tcoeff(nX1+1:end,:);
if strcmpi(type,'mahalanobis');
    CX1 = CX1./(repmat(SM.EigStd',nX1,1));
    CX2 = CX2./(repmat(SM.EigStd',nX2,1));
end
out.ShapeModel = SM;
out.CX1 = CX1;
out.CX2 = CX2;
AvgX1 = mean(CX1);
AvgX2 = mean(CX2);
% Take the within group Distances
PDX1 = sqrt(sum((CX1-repmat(AvgX1,nX1,1)).^2,2));
PDX2 = sqrt(sum((CX2-repmat(AvgX2,nX2,1)).^2,2));
out.DispX1 = mean(PDX1);
out.DispX2 = mean(PDX2);
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