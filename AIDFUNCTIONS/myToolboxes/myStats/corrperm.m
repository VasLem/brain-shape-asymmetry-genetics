function [lowci,upperci,corrs,corrobs]=corrperm(x,y,m)
% Permutation test for the correlation between two variables;
% We will fix y and permute x;
% Then calculate the correlation between x & y;
% Calcuate the m correlations;
n=length(x);
corrobs=corrcoef(x,y);
corrsamp=zeros(1,m);

parfor i=1:m,
    ind=randperm(n);
    xfor=x(ind);
    r=corrcoef(xfor,y);
    corrsamp(i)=r(2,1);
end;

% Get a 95% confidence intervals;
corrs = corrsamp;
[sortcorr,ignore]=sort(corrsamp);
lowci=sortcorr(0.975*m);
upperci=sortcorr(0.025*m+1);
end
