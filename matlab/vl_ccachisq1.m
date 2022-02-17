function [chisq, pChisq, dfs, Q2, T22, rankY] = vl_ccachisq1(Xs,Y, Q2, T22, rankY)
% An expansion of canoncorr to reuse/return Y related information and not
% recompute it if not needed. Also, it does not compute  U = (X - mean(X))*A and
%      V = (Y - mean(Y))*B.
%CANONCORR Canonical correlation analysis.
%   [A,B] = CANONCORR(X,Y) computes the sample canonical coefficients for
%   the N-by-P1 and N-by-P2 data matrices X and Y.  X and Y must have the
%   same number of observations (rows) but can have different numbers of
%   variables (cols).  A and B are P1-by-D and P2-by-D matrices, where D =
%   min(rank(X),rank(Y)).  The jth columns of A and B contain the canonical
%   coefficients, i.e. the linear combination of variables making up the
%   jth canonical variable for X and Y, respectively.  Columns of A and B
%   are scaled to make COV(U) and COV(V) (see below) the identity matrix.
%   If X or Y are less than full rank, CANONCORR gives a warning and
%   returns zeros in the rows of A or B corresponding to dependent columns
%   of X or Y.
%
%   [A,B,R] = CANONCORR(X,Y) returns the 1-by-D vector R containing the
%   sample canonical correlations.  The jth element of R is the correlation
%   between the jth columns of U and V (see below).
%
%
%   [A,B,R, STATS] = CANONCORR(X,Y) returns a structure containing
%   information relating to the sequence of D null hypotheses H0_K, that
%   the (K+1)st through Dth correlations are all zero, for K = 0:(D-1).
%   STATS contains seven fields, each a 1-by-D vector with elements
%   corresponding to values of K:
%
%      Wilks:    Wilks' lambda (likelihood ratio) statistic
%      chisq:    Bartlett's approximate chi-squared statistic for H0_K,
%                with Lawley's modification
%      pChisq:   the right-tail significance level for CHISQ
%      F:        Rao's approximate F statistic for H0_K
%      pF:       the right-tail significance level for F
%      df1:      the degrees of freedom for the chi-squared statistic,
%                also the numerator degrees of freedom for the F statistic
%      df2:      the denominator degrees of freedom for the F statistic
%
%   Example:
%
%      load carbig;
%      X = [Displacement Horsepower Weight Acceleration MPG];
%      nans = sum(isnan(X),2) > 0;
%      [A B r U V] = canoncorr(X(~nans,1:3), X(~nans,4:5));
%
%      plot(U(:,1),V(:,1),'.');
%      xlabel('0.0025*Disp + 0.020*HP - 0.000025*Wgt');
%      ylabel('-0.17*Accel + -0.092*MPG')
%
%   See also PCA, MANOVA1.

%   The STATS output contains two additional fields that remain for
%   backward compatibility but are not intended to be used any longer.

%   References:
%     [1] Krzanowski, W.J., Principles of Multivariate Analysis,
%         Oxford University Press, Oxford, 1988.
%     [2] Seber, G.A.F., Multivariate Observations, Wiley, New York, 1984.

%   Copyright 1993-2014 The MathWorks, Inc.
if nargin < 3
    Q2 =nan;
    T22 =nan;
    perm2  =nan;
    rankY =nan;
end
if nargin < 2
    error(message('stats:canoncorr:TooFewInputs'));
end

[n,p1,~] = size(Xs);
if size(Y,1) ~= n
    error(message('stats:canoncorr:InputSizeMismatch'));
elseif n == 1
    error(message('stats:canoncorr:NotEnoughData'));
end
p2 = size(Y,2);
Xs = double(Xs);
% Center the variables
Xs = Xs - mean(Xs,1);

if size(Xs,2)==1
    Q1 = Xs ./ sqrt(sum(Xs.^2));
    rankXs = ones(size(Xs, 3),1);
else
    % Factor the inputs, and find a full rank set of columns if necessary
    rankXs = zeros(size(Xs, 3),1);
    Q1 = zeros([size(Xs,1),  size(Xs,2), size(Xs,3)]);
    for i=1:size(Xs, 3)
        [Q, T11] = qr(Xs(:,:,i),0);
        rankXs(i) = sum(abs(diag(T11)) > eps(abs(T11(1)))*max(n,p1));
        Q = Q(:, 1:rankXs(i));
        Q1(1:size(Q,1),1:size(Q,2),i) = Q;
    end
end

if any(rankXs == 0)
    error(message('stats:canoncorr:BadData', 'X'));
elseif any(rankXs < p1)
    warning(message('stats:canoncorr:NotFullRank', 'X'));
end


if nargin < 3
    Y = double(Y);
    Y = Y - mean(Y,1);
    [Q2,T22] = qr(Y,0);
    rankY = sum(abs(diag(T22)) > eps(abs(T22(1)))*max(n,p2));
    if rankY == 0
        error(message('stats:canoncorr:BadData', 'Y'));
    elseif rankY < p2
        warning(message('stats:canoncorr:NotFullRank', 'Y'));
        Q2 = Q2(:,1:rankY); T22 = T22(1:rankY,1:rankY);
    end
end

ds = min(rankXs,rankY);
Ds = pagesvd(pagemtimes(permute(Q1, [2,1,3]), reshape(Q2, [size(Q2), 1])),'econ');
rs = min(max(Ds, 0), 1); % remove roundoff errs
parfor i=1:size(Ds,3)
    r = rs(:,:,i)';
    d = ds(i);
    r = r(1:d);
    rankX = rankXs(i);
    k = 0:(d-1);
    d1k = (rankX-k);
    d2k = (rankY-k);
    nondegen = r < 1;
    logLambda = -Inf( 1, d);
    logLambda(nondegen) = cumsum(log(1 - r(nondegen).^2), 'reverse');
    df1 = d1k .* d2k;
    chi = (-(n - k - .5*(rankX+rankY+3) + cumsum([0 1./r(1:(d-1))].^2)) .* logLambda)';
    ret = chi2pval(chi', df1);
    pChisq(i) = ret(1);
    chisq(i) = chi(1);
    dfs(i) = df1;
end
end


function p = chi2pval(x,v)
%FPVAL Chi-square distribution p-value function.
%   P = CHI2PVAL(X,V) returns the upper tail of the chi-square cumulative
%   distribution function with V degrees of freedom at the values in X.  If X
%   is the observed value of a chi-square test statistic, then P is its
%   p-value.
%
%   The size of P is the common size of the input arguments.  A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   See also CHI2CDF, CHI2INV.

%   References:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.4.

%   Copyright 2009 The MathWorks, Inc.


if nargin < 2
    error(message('stats:chi2pval:TooFewInputs'));
end

[errorcode,x,v] = distchck(2,x,v);

if errorcode > 0
    error(message('stats:chi2pval:InputSizeMismatch'));
end

% Return NaN for out of range parameters.
v(v <= 0) = NaN;
x(x < 0) = 0;

p = gammainc(x/2,v/2,'upper');
end
