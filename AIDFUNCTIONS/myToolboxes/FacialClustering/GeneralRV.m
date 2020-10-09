function RV = GeneralRV(A,B)
% The RV coefficient is a scalar measure for the strength of association between 2 sets of variables
% input:
% indices1, indices2: LM-indices [1 7150], each index refers to three
% coordinates (x,y,z per each LM).
% covMatrix: covariance matrix 21450 x 21450 of the 7150 LMs on x,y,z
% coordinates

[~,nA] = size(A);
%[~,nB] = size(B);


covMatrix = cov([A,B]);


numRV = trace(covMatrix(1:nA,nA+1:end)*covMatrix(nA+1:end,1:nA)); % ok

var1sum = sum(sum(covMatrix(1:nA,1:nA).^2));
var2sum = sum(sum(covMatrix(nA+1:end,nA+1:end).^2));
denRV = sqrt(var1sum*var2sum);

RV = numRV/denRV;

end