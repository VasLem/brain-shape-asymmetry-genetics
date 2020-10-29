function [Gam,Lam,lnlik,lnlika,csq]=FG_algorithm(S,n,Gam_guess);
% FG algorithm to fit CPC model
% See Flury & Gautschi (1986), SIAM J. St. Comp., 7, 169-184.
%
% Gam: p x p matrix of eigen-vectors
% Lam{i} = diagonal matrix of eigen-values for group i 
% lnlika: maximized log likelihood under saturated model
% lnlik: maximized log likelihood under CPC model
% csq: chi squared test statistic for lack of fit
%
g=length(n);
p=length(S{1});
if nargin==3
  Gam=Gam_guess;
else
  Gam=eye(p);
end;
del = 1;
Ip=eye(p);
while del > 1.e-5
  del=0;
  for i1=1:p
    ei1=Ip(:,i1);
    for i2=i1+1:p
      ei2=Ip(:,i2);
      Gami=Gam(:,[i1 i2]');
      Ei=[ei1 ei2];
      T=[];
      for i=1:g
	T{i}=Gami'*S{i}*Gami;
      end;
      Q=eye(2);
      del2=1;
      sold=0;
      while del2> 1.e-5
	R=0;
	for i=1:g
	  d1=Q(:,1)'*T{i}*Q(:,1);
	  d2=Q(:,2)'*T{i}*Q(:,2);
	  R=R+n(i)*(d1-d2)*T{i}/(d1*d2);
	end;
        if abs(R(1,2))< 100*eps
	  c=1;
	  s=0;
	  del2=0;
	else
	  d=R(2,2)-R(1,1);
	  t=(d-sqrt(d^2+4*R(1,2)^2))/(2*R(1,2));
	  if abs(t)> 1; t=-1/t; end;
	  c=1/sqrt(1+t^2);
	  s=c*t;
	  Q=[c -s; s c];
	  del2=abs(s-sold);
	  sold=s;
	end;
      end;
      del=del+abs(Q(1,2));
      J=eye(p) + Ei*(Q-eye(2))*Ei';
      Gam=Gam*J;
    end;
  end;
end
lnlik=-sum(n)*p/2;
lnlika=lnlik;
for j=1:g
  Lam{j}=diag(Gam'*S{j}*Gam);
  lnlik=lnlik-n(j)*sum(log(Lam{j}))/2;
  lnlika=lnlika-n(j)*log(det(S{j}))/2;
end
csq=2*(lnlika-lnlik);
