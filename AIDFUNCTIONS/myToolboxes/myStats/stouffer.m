function pcomb = stouffer(p)
% Stouffer et al's (1949) unweighted method for combination of 
% independent p-values via z's 
       k = length(p);
       if k==1, pcomb = p; return;end
       pcomb = erfc(sum(sqrt(2) * erfcinv(2*p))/sqrt(2*k))/2;  
end

%     if length(p)==0
%         error('pfast was passed an empty array of p-values')
%         pcomb=1;
%     else
%         pcomb = (1-erf(sum(sqrt(2) * erfinv(1-2*p))/sqrt(2*length(p))))/2;
%     end
%        pcomb = erfc(sum(sqrt(2) * erfcinv(2*p))/sqrt(2*length(p)))/2;
%        THIS WORKS AS WELL


%erfcinv(x) = erf(1-x);

%        Z = sum(sqrt(2) * erfcinv(2*p))/sqrt(2*length(p));
%        pcomb = 1-erf(Z)/2
%        
%        pcomb = erfc(Z)/2
%        
%        z = norminv(1-p);
%        Z = sum(z)/sqrt(length(p));
%        pcomb = normcdf(-Z);