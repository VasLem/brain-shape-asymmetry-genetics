function out = evaluateShapeInModel(obj,shape)
         c = getCoeff(obj,shape);
         recshape = getScan(obj,c);
         out.ProjDistances = vDistances(recshape,shape);
         out.ProjRMSE = sqrt(mean(out.ProjDistances.^2));
         [out.pF,out.pFR] = getPlausibilty(obj,c);
         out.MD = getDistance(obj,[],c,'Mahalanobis');
         out.ED = getDistance(obj,[],c,'Euclidean');
end

% function [out,percout] = getPlausibilty(obj,coeff)
% % out = getPlausibilty(obj,coeff)
% % Determine the plausibility of an entity in PCA space
% % INPUT
% % obj = PCA space object
% % coeff = PCA coefficients of the entity
% % OUTPUT
% % out = plausibility
% % percout = plausibility expressed in percentage according to maximum
% % plausibility
% % created by Peter Claes
%      if size(coeff,1)==1, coeff = coeff';end
%      out = (1/sqrt(2*pi))*exp(-0.5*(sum((coeff./obj.EigStd).^2)));
%      maxp = (1/sqrt(2*pi))*exp(-0.5*(sum((zeros(size(coeff))./obj.EigStd).^2)));
%      percout = (out/maxp)*100;
% end


%     eval([testset{i} '.Scores = nan*zeros(n,k);']);
%     eval([testset{i} '.ScoresRank = nan*zeros(n,1);']);
%     eval([testset{i} '.RMSE = nan*zeros(n,k);']);
%     eval([testset{i} '.LocalFit = nan*zeros(n,Model.Average.nrV,k);']);
%     eval([testset{i} '.MD = nan*zeros(n,k);']);
%     eval([testset{i} '.MDRank = nan*zeros(n,1);']);
%     eval([testset{i} '.ED = nan*zeros(n,k);']);
%     eval([testset{i} '.EDRank = nan*zeros(n,1);']);
%     eval([testset{i} '.nrEV = nan*zeros(n,1);']);
%     eval([testset{i} '.Diversity = nan*zeros(n,1);']);
%     eval([testset{i} '.BF.LocalError = nan*zeros(n,Model.Average.nrV,k);']);
%     eval([testset{i} '.BF.Angle = nan*zeros(n,k);']);
%     eval([testset{i} '.BF.Distinctiveness = nan*zeros(n,1);']);
%     eval([testset{i} '.WAF.LocalError = nan*zeros(n,Model.Average.nrV,k);']);
%     eval([testset{i} '.WAF.Angle = nan*zeros(n,k);']);
%     eval([testset{i} '.WAF.Distinctiveness = nan*zeros(n,1);']);