function [lowci,upperci,corrs,corrobs]=permSpearman(x,y,t)
% Permutation test for the Spearman correlation between two variables;
% We will fix y and permute x;
% Then calculate the correlation between x & y;
% Calcuate the m correlations;
warning off; %#ok<*WNOFF>
if nargin < 3, t = 1000; end
if t< 1000, t = 1000; end
n=length(x);
corrobs=Spearman(x,y);
corrsamp=zeros(1,t);
parfor i=1:t,
    ind=randperm(n);
    xfor=x(ind);
    corrsamp(i)=Spearman(xfor,y);
end;
warning on; %#ok<*WNON>
% Get a 95% confidence intervals;
corrs = corrsamp;
[sortcorr,~]=sort(corrsamp);
lowci=sortcorr(0.975*t);
upperci=sortcorr(0.025*t+1);
end
