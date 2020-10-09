function [pcomb,B,SE,Z] = MetaUVREG(p,beta,se,dfe,type)
% Stouffer et al's (1949) unweighted method for combination of 
% independent p-values via z's 
       k = length(p);%number of studies
       if k==1, pcomb = p; return;end% only one study
       if nargin<5, type = 'inverse variance';end
       B = [];SE = [];Z = [];
       switch lower(type)
           case 'inverse variance'
               w = 1./(se.^2);
               SE = sqrt(1/sum(w));
               B = sum(beta.*w)/sum(w);
               Z = B/SE;
               pcomb = 2*normcdf(-1*abs(Z),0,1);
           case 'sample size'
               z = norminv(1-(p./2),0,1).*sign(beta);
               w = sqrt(dfe);
               Z = sum(z.*w)/sqrt(sum(w.^2));
               pcomb = 2*normcdf(-1*abs(Z),0,1);
           case 'stouffer'
               pcomb = erfc(sum(sqrt(2) * erfcinv(2*p))/sqrt(2*k))/2;
           case 'fisher'
               pcomb = pfast(p);
       end   
end

function p = pfast(p)
% Fisher's (1925) method for combination of independent p-values
% Code adapted from Bailey and Gribskov (1998)
    product=prod(p);
    n=length(p);
    if n<=0
        error('pfast was passed an empty array of p-values')
    elseif n==1
        p = product;
        return
    elseif product == 0
        p = 0;
        return
    else
        x = -log(product);
        t=product;
        p=product;
        for i = 1:n-1
            t = t * x / i;
            p = p + t;
        end
    end  
end