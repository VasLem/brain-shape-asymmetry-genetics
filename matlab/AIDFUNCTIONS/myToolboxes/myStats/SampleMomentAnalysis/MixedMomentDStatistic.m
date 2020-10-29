function [out,fig] = MixedMomentDStatistic(X1,X2,t,mcp,corr)
% function to test mixed moment differences (difference of correlation structure, subspaces) between shape samples X1 and X2
% X1 = N by K matrix, N number of objects in sample 1, K the amount of
% shape landmarks times 3; Therefore Shape landmark data is to be
% vectorized: [Xa Ya Za Xb Yb Zb ... Xz Yz Zz]; for landmarks a to z.
% X2 = M by K matrix, M number of objects in sample 2, K the amount of
% shape landmarks times 3; Therefore Shape landmark data is to be
% vectorized: [Xa Ya Za Xb Yb Zb ... Xz Yz Zz]; for landmarks a to z.
% t = the amount of permutations, best used in combination with matlab parpool and parallel computing.
% if t==0, then permutation testing is ommited.
% mcp, maximum number of principal components to dedscribe subspaces,
% default is minimum(N,M);
% corr = correction of variable standard deviations to 1, default is false.
% 
% out a structure containing the statistical results
%
% Created by Peter Claes, peter.claes@kuleuven.be
% This routine was written for and used in  
% Claes, P., et al. (2012). "Sexual Dimorphism in Multiple Aspects of 3D Facial Symmetry & Asymmetry defined by Spatially-dense Geometric Morphometrics." Journal of Anatomy 221(2): 97-114.
% Claes, P., et al. (2015). "An investigation of matching symmetry in the human pinnae with possible implications for 3D ear recognition and sound localization." Journal of Anatomy 226(1): 60-72.
% Please cite these works when using this function, thank you
         if nargin < 5, corr = false; end
         if nargin < 4, mcp = +inf; end
         if nargin < 3, t = 0; end
         nX1 = size(X1,1);
         nX2 = size(X2,1);
         mcp = min([nX1 nX2 mcp]);
         % Standarizing the Data
           % Centering
           disp('centering');
           X1 = X1-repmat(mean(X1),nX1,1);
           X2 = X2-repmat(mean(X2),nX2,1);
           if corr% scaling
              disp('scaling');
              X1 = X1./repmat(std(X1),nX1,1);
              X2 = X2./repmat(std(X2),nX2,1);
           end
         % Create Subspaces
           [SSpaceX1,EigValX1] = createSubspaceSVD(X1);
           [SSpaceX2,EigValX2] = createSubspaceSVD(X2);
           SSpaceX1 = SSpaceX1(:,1:mcp);EigValX1 = EigValX1(1:mcp);
           SSpaceX2 = SSpaceX2(:,1:mcp);EigValX2 = EigValX2(1:mcp);
         % Getting Distance Statistic
           [DStat,A] = getSubspaceDistances(SSpaceX1,SSpaceX2,'projection');
           EigValDStat = cumsum(abs(EigValX1-EigValX2));
         % generating Effect output  
           out.SSpaceX1 = SSpaceX1;
           out.EigValX1 = EigValX1;
           out.SSpaceX2 = SSpaceX2;
           out.EigValX2 = EigValX2;
           out.DStat = DStat;
           out.A = A;
           out.EigValDStat = EigValDStat;
         % Permutation test  
           if t<=0, return; end
           StatCount = false(length(DStat),t);
           EigValStatCount = false(length(EigValDStat),t);
           DStatpermC = zeros(t,length(DStat));
           EigValDStatpermC = zeros(t,length(EigValDStat));
           nT = nX1+nX2;
           X = [X1; X2];
           disp('Permuting');
           tic;
           parfor i=1:t
                  ind = randperm(nT);
                  X1perm = X(ind(1:nX1),:); %#ok<*PFBNS>
                  X2perm = X(ind(nX1+1:end),:);
                  [SSpaceX1perm,EigValX1perm] = createSubspaceSVD(X1perm);
                  [SSpaceX2perm,EigValX2perm] = createSubspaceSVD(X2perm);
                  SSpaceX1perm = SSpaceX1perm(:,1:mcp);EigValX1perm = EigValX1perm(1:mcp);
                  SSpaceX2perm = SSpaceX2perm(:,1:mcp);EigValX2perm = EigValX2perm(1:mcp);
                  DStatperm = getSubspaceDistances(SSpaceX1perm,SSpaceX2perm,'projection');
                  EigValDStatperm = cumsum(abs(EigValX1perm-EigValX2perm));
                  StatCount(:,i) = (DStatperm>=DStat)';
                  EigValStatCount(:,i) = EigValDStatperm>=EigValDStat;
                  DStatpermC(i,:) = DStatperm;
                  EigValDStatpermC(i,:) = EigValDStatperm;
           end
           toc;
           fig = figure;
           set(gca,'Xlim',[1 mcp]);
           hold on;
           for i=1:1:t
               plot(DStatpermC(i,:),'r-.');
               drawnow;   
           end
           plot(DStat,'b','linewidth',3);
           out.pperm = ((sum(StatCount,2)+1)/(t+1))';
           out.EigValpperm = ((sum(EigValStatCount,2)+1)/(t+1))';
           out.DStatpermC = DStatpermC;
           out.EigValDStatpermC = EigValDStatpermC;
end

function [QF,EF] = createSubspaceSVD(F)
         if size(F,1) < size(F,2)
            [QF,~,EF] = my_princomp(F,'econ');
         else
            [QF,~,EF] = my_princomp(F);
         end
end

function [coeff, score, latent, tsquare] = my_princomp(x,econFlag)
%PRINCOMP Principal Components Analysis.
%   COEFF = PRINCOMP(X) performs principal components analysis on the N-by-P
%   data matrix X, and returns the principal component coefficients, also
%   known as loadings.  Rows of X correspond to observations, columns to
%   variables.  COEFF is a P-by-P matrix, each column containing coefficients
%   for one principal component.  The columns are in order of decreasing
%   component variance.
%
%   [COEFF, SCORE] = PRINCOMP(X) returns the principal component scores,
%   i.e., the representation of X in the principal component space.  Rows
%   of SCORE correspond to observations, columns to components.
%
%   [COEFF, SCORE, LATENT] = PRINCOMP(X) returns the principal component
%   variances, i.e., the eigenvalues of the covariance matrix of X, in
%   LATENT.
%
%   [COEFF, SCORE, LATENT, TSQUARED] = PRINCOMP(X) returns Hotelling's
%   T-squared statistic for each observation in X.
%
%   When N <= P, SCORE(:,N:P) and LATENT(N:P) are necessarily zero, and the
%   columns of COEFF(:,N:P) define directions that are orthogonal to X.
%
%   [...] = PRINCOMP(X,'econ') returns only the elements of LATENT that are
%   not necessarily zero, i.e., when N <= P, only the first N-1, and the
%   corresponding columns of COEFF and SCORE.  This can be significantly
%   faster when P >> N.
%
%   See also FACTORAN, PCACOV, PCARES.

%   References:
%     [1] Jackson, J.E., A User's Guide to Principal Components,
%         Wiley, 1988.
%     [2] Jolliffe, I.T. Principal Component Analysis, 2nd ed.,
%         Springer, 2002.
%     [3] Krzanowski, W.J., Principles of Multivariate Analysis,
%         Oxford University Press, 1988.
%     [4] Seber, G.A.F., Multivariate Observations, Wiley, 1984.

%   Copyright 1993-2004 The MathWorks, Inc.
%   $Revision: 2.9.2.6 $  $Date: 2004/07/28 04:39:01 $

    [n,p] = size(x);
    r = min(n-1,p); % max possible rank of X0
    %disp('ok');
    % When X has more variables than observations, the default behavior is to
    % return all the pc's, even those that have zero variance.  When econFlag
    % is 'econ', those will not be returned.
    if nargin < 2, econFlag = 0; end

    % The principal component coefficients are the eigenvectors of
    % S = X0'*X0./(n-1), but computed using SVD.
    [U,sigma,coeff] = svd(x,econFlag); % put in 1/sqrt(n-1) later

    if nargout < 2
        % When econFlag is 'econ', only (n-1) components should be returned.
        % See comment below.
        if (n <= p) && isequal(econFlag, 'econ')
            coeff(:,n) = [];
        end

    else
        % Project X0 onto the principal component axes to get the scores.
        if n == 1 % sigma might have only 1 row
            sigma = sigma(1);
        else
            sigma = diag(sigma);
        end
        score = U .* repmat(sigma',n,1); % == x0*coeff
        sigma = sigma./ sqrt(n-1);
        % When X has at least as many variables as observations, eigenvalues
        % n:p of S are exactly zero.
        if n <= p
            % When econFlag is 'econ', nothing corresponding to the zero
            % eigenvalues should be returned.  svd(,'econ') won't have
            % returned anything corresponding to components (n+1):p, so we
            % just have to cut off the n-th component.
            if isequal(econFlag, 'econ')
                sigma(n) = [];
                coeff(:,n) = [];
                score(:,n) = [];

            % Otherwise, set those eigenvalues and the corresponding scores to
            % exactly zero.  svd(,0) won't have returned columns of U
            % corresponding to components (n+1):p, need to fill those out.
            else
                sigma(n:p,1) = 0; % make sure this extends as a column
                score(:,n:p) = 0;
            end
        end
        % The variances of the pc's are the eigenvalues of S = X0'*X0./(n-1).
        latent = sigma.^2;
        % Hotelling's T-squared statistic is the sum of squares of the
        % standardized scores, i.e., Mahalanobis distances.  When X appears to
        % have column rank < r, ignore components that are orthogonal to the
        % data.
        if nargout == 4
            q = sum(sigma > max(n,p).*eps(sigma(1)));
            if q < r
                warning('stats:princomp:colRankDefX', ...
                        ['Columns of X are linearly dependent to within machine precision.\n' ...
                         'Using only the first %d components to compute TSQUARED.'],q);
            end
            tsquare = (n-1) .* sum(U(:,1:q).^2,2); % == sum(score*diag(1./sigma),2)
        end
    end
end

function [d, A, U, V] = getSubspaceDistances(QF,QG,type)
         if nargin < 3, type = 'Projection Frobenius'; end
         k1 = size(QF,2);
         k2 = size(QG,2);
         [A, U, V] = mySubspaceAngles(QF,QG);
         nA = length(A);
         d = zeros(1,nA);
         for i=1:1:nA
             k= i;
             switch lower(type)
                 case 'projection frobenius'
                     d(i) = sqrt(k1+k2-2*sum(cos(A(1:i)).^2));
                 case 'projection'
                     d(i) = sqrt(k-sum(cos(A(1:i)).^2));
                 case 'procrustes'
                     d(i) = 2*sqrt(sum(sin(A(1:i)/2).^2));
                 case 'mitteroecker'
                     d(i) = sqrt(sum(log(cos(A(1:i))).^2));
                 case 'average'
                     d(i) = k/sum(cos(A(1:i)));
                 case 'krzanowski'
                     d(i) = k-sum(cos(A(1:i)).^2);
             end
         end
end

