function [Z,lambda,VAR,CVAR,Q]=pcoa(Dist)
% Gower's Principal Coordinate Analysis.
% Adapted from Leslie Marcus's PRCRND2.m program
% in Reyment RA, Joreskog KG {Marcus LF append}(1993)
% Applied factor analysis in the natural sciences.
% Cambridge University Press, New York.
% format: [Z,lambda,VAR,CVAR,Q]=PCOA(Dist);
% where, Dist is any dissimilarity matrix.
% Z=principal coordinates
% Similarity matrices should be converted to
% dissimilarity matrices (e.g., Dist=ones(N,N)-NNESS) prior
% to calling PCOA(Dist);
% output:
% lambda=eigenvalues
% Z=Station coordinates for positive eigenvalues.
% (n.b., for non-metric indices only the first few columns
% of Z may be useful - check the lambda vector for negative
% eigenvalues!)
% VAR, CVAR=Variance and cumulative variance
% Q=standardized association matrix.
% Refs: Gower, J. C. 1966. Biometrika 53:325.
% L Marcus's appendix in Reyment & Joreskog. 1993. Applied
% Factor Analysis, Cambridge U. Press.,
% Adapted by E. Gallagher Environmental Sciences Program
% Eugene.Gallagher@umb.edu, revised: 6/2/95
[N,N]=size(Dist);
Dist=Dist.^2; % Must convert Euclidean distances to Dist.^2
               % prior to transformation in order to have PCO distances
               % match Euclidean distances, see discussion of equ. (4)
               % in Gower (1966)
cmean=mean(Dist); % column means of squared distance matrix
cm=cmean(ones(N,1),:); % creates order N matrix, each row=cmean;
% create matrix of grand means: gm
gm=mean(cmean);gm=gm(ones(N,1),:);gm=gm';gm=gm(ones(N,1),:);
Q=-.5*(Dist-cm-cm'+gm);
[V,S]=eig(Q); % eigenanalysis of Q.
% these three statements sort eigenvalues in ascending order and
% sort eigenvectors accordingly:
[lambda,k]=sort(diag(S));
lambda=flipud(lambda); % rearrange to descending order
V=V(:,k); % sort eigenvectors in descending order
V=fliplr(V); % eigenvectors corresponding to sorted lambda
VAR=abs(lambda)/sum(abs(lambda))*100; % modulus
CVAR=cumsum(VAR);
posl=find(lambda>0);
Z=V(:,posl)*diag(sqrt(lambda(posl)));
end

