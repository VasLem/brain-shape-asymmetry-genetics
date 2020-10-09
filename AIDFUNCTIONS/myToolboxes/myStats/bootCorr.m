function [avg,upper,lower]=bootCorr(x,y,t)
% Bootstrap for the correlation between two variables;
% We will fix y and permute x;
% Then calculate the correlation between x & y;
% Calcuate the m correlations;
        if nargin < 3, t = 1000; end
        if t< 1000, t = 1000; end
        n=length(x);
        corrsamp = zeros(1,t);
        parfor i=1:t,
            ind = randsample((1:n),n,true);
            xfor=x(ind); %#ok<*PFBNS>
            yfor = y(ind);
            r=corrcoef(xfor,yfor);
            corrsamp(i)=r(2,1);
        end;
        corrsamp = sort(corrsamp);
        upper = corrsamp(0.975*t);
        lower= corrsamp(0.025*t+1);
        avg = mean(corrsamp);
end