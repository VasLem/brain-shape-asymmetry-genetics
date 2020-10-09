function [theta,U,V] = mySubspaceAngles(QF,QG)
%SUBSPACEA angles between subspaces
%  subspacea(F,G,A)
%  Finds all min(size(orth(F),2),size(orth(G),2)) principal angles
%  between two subspaces spanned by the columns of matrices F and G 
  threshold=sqrt(2)/2; % Define threshold for determining when an angle is small
  q = min(size(QF,2),size(QG,2));
  %q = min(q,ncomp);  
  QF = QF(:,1:q);
  QG = QG(:,1:q);
  [Ys,s,Zs] = svd(QF'*QG,0);
  if size(s,1)==1
    % make sure s is column for output
    s=s(1);
  end
  s = min(diag(s),1);
  theta = max(acos(s),0);
  U = QF*Ys;
  V = QG*Zs;
  indexsmall = s > threshold;
  if max(indexsmall) % Check for small angles and recompute only small   
    RF = U(:,indexsmall); 
    RG = V(:,indexsmall); 
    %[Yx,x,Zx] = svd(RG-RF*(RF'*RG),0);
    [~,x,Zx] = svd(RG-QF*(QF'*RG),0); % Provides more accurate results
    if size(x,1)==1
      % make sure x is column for output
      x=x(1);
    end
    Tmp = fliplr(RG*Zx);
    V(:,indexsmall) = Tmp(:,indexsmall);
    U(:,indexsmall) = RF*(RF'*V(:,indexsmall))*...   
    diag(1./s(indexsmall)); 
    x = diag(x);               
    thetasmall=flipud(max(asin(min(x,1)),0));
    theta(indexsmall) = thetasmall(indexsmall);
  end     
end