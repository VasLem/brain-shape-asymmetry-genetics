function D = mydistfun(X0,X)
% calculation of distance
%X0 is a 1-by-p vector containing a single point (observation) of the input data matrix X.
%X is an n-by-p matrix of points.
%D is an n-by-1 vector of distances, and D(k) is the distance between observations X0 and X(k,:).
idx = ismember(X0,X,'rows');
D = X(idx,:)';
end