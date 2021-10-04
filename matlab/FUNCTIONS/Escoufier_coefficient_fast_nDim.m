function RV = Escoufier_coefficient_fast_nDim(covMatrix)
% The RV coefficient is a scalar measure for the strength of association between 2 sets of variables
% input:
% indices1, indices2: LM-indices [1 7150], each index refers to three
% coordinates (x,y,z per each LM). OR (t); (x,y,z,t)
% covMatrix: covariance matrix 21450 x 21450 of the 7150 LMs on x,y,z
% coordinates
ndim = size(covMatrix,1)/2;
numRV = trace(covMatrix(1:ndim,ndim+1:end)*covMatrix(ndim+1:end,1:ndim)); % ok

var1sum = sum(sum(covMatrix(1:ndim,1:ndim).^2));
var2sum = sum(sum(covMatrix(ndim+1:end,ndim+1:end).^2));
denRV = sqrt(var1sum*var2sum);

RV = numRV/denRV;

end

% test for only thickness(t) ndim =1
% RV1 = Escoufier_coefficient_fast_nDim(covMatrix1); 
% sqrt(RV1) % = abs(off-diag value of corrcov) always positive values
% corrcov(covMatrix1)
