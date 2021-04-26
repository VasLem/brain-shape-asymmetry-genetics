function out = computeAmmiModel(Shapes)
inputShapes = Shapes;
n = size(inputShapes,3)/2;
LH = inputShapes(:,:,1:n);
RH =  inputShapes(:,:,n+1:2*n);
avgShapes = (LH + RH)/2;
sideDif =LH  - RH ;
%%
D = sqrt(sum(mean(sideDif,3).^2,2));
%%
avgShape = mean(avgShapes,3);
indDif = avgShapes - avgShape ;

I = mean(sqrt(sum(indDif.^2,2)),3);
%%
flucMat = LH - RH + avgShape;
flucMat= squeeze(sqrt(sum((flucMat).^2,2)));
[latent, ~, latentHigh]  = pa_test(flucMat', 100,0.05);
numComponents = sum(latent > latentHigh');
 [coeff, score, ~] =  pca(flucMat','Centered',false);
rec = coeff(:, 1:numComponents) * score(:, 1:numComponents)' ;
 F = mean(rec,2);
out.LM.I = I';
out.LM.D = D';
out.LM.F = F';
out.Raw.F = rec;
end

function [latent, latentLow, latentHigh] = pa_test(x, nShuffle, alpha, varargin)
% Parallel Analysis (PA) to for determining the number of components to
% retain from PCA. component is retained if the associated eigenvalue is bigger than the 95th of the distribution of eigenvalues derived from the random data.
% Syntax:
% ======
% pa_test(x, nShuffle, alpha, princomp_parameters[ ])
%       x - the data matrix (nXp where n is the number of observation and p is dimension of each observation)
%       nShuffle - number of shuffles. optional, default = 100
%       alpha - significance level. optional, default 0.05
%       princomp_parameters - parameters to pass to the princomp function (see help princomp). optional, default ={true,'Centered',false}
% Background:
% ==========
% From Wikipedia: http://en.wikipedia.org/wiki/Factor_analysis
% Horn's Parallel Analysis (PA):
% A Monte-Carlo based simulation method that compares the observed eigenvalues with those obtained from uncorrelated normal variables.
% A factor or component is retained if the associated eigenvalue is bigger than the 95th of the distribution of eigenvalues derived from the random data.
% PA is one of the most recommendable rules for determining the number of components to retain, but only few programs include this option.
% References:
% * Ledesma, R.D.; Valero-Mora, P. (2007). "Determining the Number of Factors to Retain in EFA: An easy-to-use computer program for carrying out Parallel Analysis". Practical Assessment Research & Evaluation 12 (2): 1�11.
if isempty(varargin)
        varargin = {'Centered',true};
end
if ~exist('nShuffle', 'var')
        nShuffle = 100;
end
if ~exist('alpha', 'var')
        alpha = 0.05;
end
[~,~,latent, ~] = pca(x, varargin{:});
xShuffle = x;
latentShuffle = zeros(length(latent), nShuffle);
for iShuffle = 1:nShuffle
        for dim = 1:size(x,2)
                xShuffle(:, dim) = x(randperm(size(x,1)), dim);
        end
        [~, ~ ,latentShuffle(:,iShuffle)] = pca(xShuffle, varargin{:});
end
latentHigh = quantile(latentShuffle', 1-alpha);
latentLow = quantile(latentShuffle', alpha);
if nargout == 0
        plot(latent,'bo-');
        hold on;
        plot(latentHigh,'ro-');
        plot(latentLow,'go-');
        hold off;
        legend('data',sprintf('shuffled %d%%', 100*(1-alpha)),sprintf('shuffled %d%%', 100*alpha));
end
end

