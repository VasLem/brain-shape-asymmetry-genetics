function [c,pc]=permCorr(x,y,t)
% Permutation test for the correlation between two variables;
% We will fix y and permute x;
% Then calculate the correlation between x & y;
% Calcuate the m correlations;
if nargin < 3, t = 1000; end
if t< 1000, t = 1000; end
n=length(x);
corrobs=corrcoef(x,y);
c = corrobs(2,1);
corrCount=zeros(1,t);
if t==0, pc = nan; return; end
parfor i=1:t,
    ind=randperm(n);
    xfor=x(ind);
    r=corrcoef(xfor,y);
    corrsamp(i)=r(2,1);
    switch sign(c)
        case 1
            corrCount(i) = corrsamp(i)>=c;
        case -1
            corrCount(i) = corrsamp(i)<=c;
    end
end;
% parfor i=1:t,
%     ind=randperm(n);
%     xfor=x(ind);
%     r=corrcoef(xfor,y);
%     corrsamp(i)=r(2,1);
%     corrCount(i) = abs(corrsamp(i))>=abs(c)
% end;
pc = sum(corrCount)/t;
% % Get a 95% confidence intervals;
% corrs = corrsamp;
% [sortcorr]=sort(corrsamp);
% lowci=sortcorr(0.975*t);
% upperci=sortcorr(0.025*t+1);
end



% function [lowci,upperci,corrs,corrobs]=permCorr(x,y,t)
% % Permutation test for the correlation between two variables;
% % We will fix y and permute x;
% % Then calculate the correlation between x & y;
% % Calcuate the m correlations;
% if nargin < 3, t = 1000; end
% if t< 1000, t = 1000; end
% n=length(x);
% corrobs=corrcoef(x,y);
% corrsamp=zeros(1,t);
% 
% parfor i=1:t,
%     ind=randperm(n);
%     xfor=x(ind);
%     r=corrcoef(xfor,y);
%     corrsamp(i)=r(2,1);
% end;
% 
% % Get a 95% confidence intervals;
% corrs = corrsamp;
% [sortcorr]=sort(corrsamp);
% lowci=sortcorr(0.975*t);
% upperci=sortcorr(0.025*t+1);
% end
