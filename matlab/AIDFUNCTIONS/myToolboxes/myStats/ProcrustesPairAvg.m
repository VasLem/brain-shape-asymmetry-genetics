function [out] = ProcrustesPairAvg(X1,X2,t)
% This function computed the Directional Procrustes Differences based on
% the averaged of balanced (paired) input arguments X1 & X2;
if nargin < 3, t = 100; end
nX1 = size(X1,1);
nrV = size(X1,2)/3;
AvgX1 = mean(X1);
AvgX2 = mean(X2);
Differences = reshape(AvgX1-AvgX2,3,nrV);
SS = sum(Differences.^2);
Distances = sqrt(SS);
PD = sum(SS);
RMSE = sqrt((mean(Distances.^2)));
DistancesCount = false(nrV,t);
PDCount = false(1,t);
tic
disp('Permuting');
parfor i=1:t
    r = randi(2,nX1,1);
    X1for = zeros(size(X1));
    X2for = zeros(size(X2));
    for k=1:1:nX1
        switch r(k)
            case 1
                X1for(k,:) = X1(k,:);
                X2for(k,:) = X2(k,:);
            case 2
                X1for(k,:) = X2(k,:);
                X2for(k,:) = X1(k,:);
        end
    end
    AvgX1for = mean(X1for);
    AvgX2for = mean(X2for);
    Differencesfor = reshape(AvgX1for-AvgX2for,3,nrV);
    SSfor = sum(Differencesfor.^2);
    Distancesfor = sqrt(SSfor);
    PDfor = sum(SSfor);
    DistancesCount(:,i) = Distancesfor >= Distances;
    PDCount(i) = PDfor >= PD;
end
disp('Done');
toc;
out.AvgX1 = AvgX1;
out.AvgX2 = AvgX2;
out.Differences = Differences;
out.Distances = Distances;
out.SS = SS;
out.PD = PD;
out.RMSE = RMSE;
out.pPD = (sum(PDCount)/t);
out.pDistances = sum(DistancesCount,2)/t;
out.p1Significant = out.pDistances <= 0.1;
out.p1Perc = (sum(out.p1Significant)/nrV)*100;
out.p05Significant = out.pDistances <= 0.05;
out.p05Perc = (sum(out.p05Significant)/nrV)*100;
out.p001Significant = out.pDistances <= 0.001;
out.p001Perc = (sum(out.p001Significant)/nrV)*100;
% Getting the magnitude of the Fluctuating differences
% Getting the paired differences
% out.FL.Difference = (X1-X2)- repmat(out.DI.Differences(:)',nX1,1);
end