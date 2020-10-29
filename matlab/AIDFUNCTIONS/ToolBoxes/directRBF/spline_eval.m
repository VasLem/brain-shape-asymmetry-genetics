function s=spline_eval(cent,coeff,u)
%
% Function which directly evaluates the RBF spline
%
% s(u)=linear poly +\sum_i=1^n coeff(i)*phi(u-cent(i,:) 
%
% at the point u where the linear polynomial is given  
% in terms of its monomial cofficients. 
%
% Syntax     s=eval_direct(cent,coeff,u)
%
% Inputs
%    cent   n by dim array of coordinates for the centres
%    coeff  n+dim +1 vector of coefficients with the linear
%           polynomial part last.
%    u      row vector    point at which to evaluate
%
% Output
%    s    Value of the RBF at position u
%
% Code is written for clarity rather than Matlab efficiency
[n dim] = size(cent);
s=0;
for j=1:n
   s = s + coeff(j)*phi(  norm(u - cent(j,:) )   );
end
% Now the linear polynomial bit
s=s+coeff(n+1);               % The constant
for i=1:dim
   s=s+coeff(i+n+1)*u(i);     % The various components of u
end

