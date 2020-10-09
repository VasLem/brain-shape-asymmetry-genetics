function RV = Escoufier_coefficient_fast(covMatrix)
% The RV coefficient is a scalar measure for the strength of association between 2 sets of variables
% input:
% indices1, indices2: LM-indices [1 7150], each index refers to three
% coordinates (x,y,z per each LM).
% covMatrix: covariance matrix 21450 x 21450 of the 7150 LMs on x,y,z
% coordinates

numRV = trace(covMatrix(1:3,3+1:end)*covMatrix(3+1:end,1:3)); % ok

var1sum = sum(sum(covMatrix(1:3,1:3).^2));
var2sum = sum(sum(covMatrix(3+1:end,3+1:end).^2));
denRV = sqrt(var1sum*var2sum);

RV = numRV/denRV;

end