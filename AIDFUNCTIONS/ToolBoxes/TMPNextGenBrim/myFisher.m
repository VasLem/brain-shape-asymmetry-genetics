function [chi,p] = myFisher(p)
         Lp = log(p);
         chi = -2*sum(Lp);
         p = chi2pdf(chi,2*length(p));
end